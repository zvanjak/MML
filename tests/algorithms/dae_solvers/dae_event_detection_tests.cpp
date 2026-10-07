///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///  DAE zero-crossing events, terminal actions, and consistent reinitialization      ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/DAESolvers.h>
#endif

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::DAEEventDetectionTests {

	class TimedEventDAE : public IODESystemDAEWithEvents {
		Real _eventTime;
		EventAction _action;

	public:
		TimedEventDAE(Real eventTime, EventAction action)
			: _eventTime(eventTime), _action(action) {}

		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }

		void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& dxdt) const override {
			dxdt[0] = REAL(1.0);
		}

		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y,
			Vector<Real>& constraints) const override {
			constraints[0] = x[0] + y[0] - REAL(1.0);
		}

		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = REAL(0.0);
		}
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = REAL(0.0);
		}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = REAL(1.0);
		}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = REAL(1.0);
		}

		int getNumEvents() const override { return 1; }
		Real eventFunction(int, Real t, const Vector<Real>&, const Vector<Real>&) const override {
			return t - _eventTime;
		}
		EventDirection getEventDirection(int) const override { return EventDirection::Increasing; }
		EventAction getEventAction(int) const override { return _action; }

		void handleEvent(int, Real, Vector<Real>& x, Vector<Real>& y) const override {
			x[0] = REAL(0.25);
			y[0] = REAL(42.0); // Deliberately inconsistent; the solver must reinitialize y.
		}
	};

	DAESolverConfig EventTestConfig() {
		DAESolverConfig config;
		config.step_size = REAL(0.2);
		config.newton_tol = TOL(1e-12, 1e-5);
		config.constraint_tol = TOL(1e-10, 1e-5);
		return config;
	}

	TEST_CASE("DAE events locate terminal root and stop", "[dae][events]") {
		TimedEventDAE system(REAL(0.35), EventAction::Stop);
		Vector<Real> x0{REAL(0.0)}, y0{REAL(1.0)};
		auto result = SolveDAEBackwardEulerWithEvents(system, REAL(0.0), x0, y0,
			REAL(1.0), EventTestConfig());

		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE(result.terminated_by_event);
		REQUIRE(result.events.size() == 1);
		REQUIRE_THAT(result.events[0].time, WithinAbs(REAL(0.35), TOL(1e-9, 1e-4)));
		REQUIRE(result.events[0].direction == EventDirection::Increasing);
		REQUIRE_THAT(result.final_time, WithinAbs(REAL(0.35), TOL(1e-9, 1e-4)));
	}

	TEST_CASE("DAE non-terminal event continues integration", "[dae][events]") {
		TimedEventDAE system(REAL(0.35), EventAction::Continue);
		Vector<Real> x0{REAL(0.0)}, y0{REAL(1.0)};
		auto result = SolveDAEBackwardEulerWithEvents(system, REAL(0.0), x0, y0,
			REAL(1.0), EventTestConfig());

		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE_FALSE(result.terminated_by_event);
		REQUIRE(result.events.size() == 1);
		REQUIRE_THAT(result.final_time, WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.final_differential_state[0], WithinAbs(REAL(1.0), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.final_algebraic_state[0], WithinAbs(REAL(0.0), TOL(1e-9, 1e-4)));
	}

	TEST_CASE("DAE restart event recomputes consistent algebraic state", "[dae][events]") {
		TimedEventDAE system(REAL(0.5), EventAction::Restart);
		Vector<Real> x0{REAL(0.0)}, y0{REAL(1.0)};
		auto result = SolveDAEBackwardEulerWithEvents(system, REAL(0.0), x0, y0,
			REAL(1.0), EventTestConfig());

		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE(result.events.size() == 1);
		REQUIRE_THAT(result.events[0].differential_state[0], WithinAbs(REAL(0.5), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.final_differential_state[0], WithinAbs(REAL(0.75), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.final_algebraic_state[0], WithinAbs(REAL(0.25), TOL(1e-9, 1e-4)));
		REQUIRE(result.integration.final_constraint_norm <= EventTestConfig().constraint_tol);
	}

} // namespace MML::Tests::Algorithms::DAEEventDetectionTests
