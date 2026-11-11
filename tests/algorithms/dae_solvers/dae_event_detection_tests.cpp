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
			dxdt[0] = static_cast<MML::Real>(1.0);
		}

		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y,
			Vector<Real>& constraints) const override {
			constraints[0] = x[0] + y[0] - static_cast<MML::Real>(1.0);
		}

		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = static_cast<MML::Real>(0.0);
		}
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = static_cast<MML::Real>(0.0);
		}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = static_cast<MML::Real>(1.0);
		}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override {
			value(0, 0) = static_cast<MML::Real>(1.0);
		}

		int getNumEvents() const override { return 1; }
		Real eventFunction(int, Real t, const Vector<Real>&, const Vector<Real>&) const override {
			return t - _eventTime;
		}
		EventDirection getEventDirection(int) const override { return EventDirection::Increasing; }
		EventAction getEventAction(int) const override { return _action; }

		void handleEvent(int, Real, Vector<Real>& x, Vector<Real>& y) const override {
			x[0] = static_cast<MML::Real>(0.25);
			y[0] = static_cast<MML::Real>(42.0); // Deliberately inconsistent; the solver must reinitialize y.
		}
	};

	DAESolverConfig EventTestConfig() {
		DAESolverConfig config;
		config.step_size = static_cast<MML::Real>(0.2);
		config.newton_tol = TOL(1e-12, 1e-5);
		config.constraint_tol = TOL(1e-10, 1e-5);
		return config;
	}

	TEST_CASE("DAE events locate terminal root and stop", "[dae][events]") {
		TimedEventDAE system(static_cast<MML::Real>(0.35), EventAction::Stop);
		Vector<Real> x0{static_cast<MML::Real>(0.0)}, y0{static_cast<MML::Real>(1.0)};
		auto result = SolveDAEBackwardEulerWithEvents(system, static_cast<MML::Real>(0.0), x0, y0,
			static_cast<MML::Real>(1.0), EventTestConfig());

		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE(result.terminated_by_event);
		REQUIRE(result.events.size() == 1);
		REQUIRE_THAT(result.events[0].time, WithinAbs(static_cast<MML::Real>(0.35), TOL(1e-9, 1e-4)));
		REQUIRE(result.events[0].direction == EventDirection::Increasing);
		REQUIRE_THAT(result.final_time, WithinAbs(static_cast<MML::Real>(0.35), TOL(1e-9, 1e-4)));
	}

	TEST_CASE("DAE non-terminal event continues integration", "[dae][events]") {
		TimedEventDAE system(static_cast<MML::Real>(0.35), EventAction::Continue);
		Vector<Real> x0{static_cast<MML::Real>(0.0)}, y0{static_cast<MML::Real>(1.0)};
		auto result = SolveDAEBackwardEulerWithEvents(system, static_cast<MML::Real>(0.0), x0, y0,
			static_cast<MML::Real>(1.0), EventTestConfig());

		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE_FALSE(result.terminated_by_event);
		REQUIRE(result.events.size() == 1);
		REQUIRE_THAT(result.final_time, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.final_differential_state[0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.final_algebraic_state[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-9, 1e-4)));
		REQUIRE(result.integration.lu_factorizations > 0);
		REQUIRE(result.integration.linear_solves >= result.integration.lu_factorizations);
		REQUIRE(result.integration.jacobian_evaluations == result.integration.lu_factorizations);
	}

	TEST_CASE("DAE restart event recomputes consistent algebraic state", "[dae][events]") {
		TimedEventDAE system(static_cast<MML::Real>(0.5), EventAction::Restart);
		Vector<Real> x0{static_cast<MML::Real>(0.0)}, y0{static_cast<MML::Real>(1.0)};
		auto result = SolveDAEBackwardEulerWithEvents(system, static_cast<MML::Real>(0.0), x0, y0,
			static_cast<MML::Real>(1.0), EventTestConfig());

		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE(result.events.size() == 1);
		REQUIRE_THAT(result.events[0].differential_state[0], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.final_differential_state[0], WithinAbs(static_cast<MML::Real>(0.75), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.final_algebraic_state[0], WithinAbs(static_cast<MML::Real>(0.25), TOL(1e-9, 1e-4)));
		REQUIRE(result.integration.final_constraint_norm <= EventTestConfig().constraint_tol);
	}

	TEST_CASE("DAE event root budget exhaustion is reported instead of an incoherent event",
		"[dae][events][MinimalMathLibrary-hrdm]") {
		TimedEventDAE system(static_cast<MML::Real>(0.3), EventAction::Continue);
		Vector<Real> x0{static_cast<MML::Real>(0.0)}, y0{static_cast<MML::Real>(1.0)};
		DAESolverConfig config = EventTestConfig();
		config.step_size = static_cast<MML::Real>(1.0);
		DAEEventConfig eventConfig;
		eventConfig.max_root_iterations = 1;

		auto result = SolveDAEBackwardEulerWithEvents(system, static_cast<MML::Real>(0.0), x0, y0,
			static_cast<MML::Real>(1.0), config, eventConfig);

		INFO(result.integration.error_message);
		REQUIRE(result.integration.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(result.events.empty());
		REQUIRE_FALSE(result.terminated_by_event);
		REQUIRE(result.final_time == static_cast<MML::Real>(0.0));
		REQUIRE(result.integration.error_message.find("max_root_iterations") != std::string::npos);

		SECTION("Sufficient budget yields a coherent time/state/residual triple") {
			eventConfig.max_root_iterations = 80;
			auto ok = SolveDAEBackwardEulerWithEvents(system, static_cast<MML::Real>(0.0), x0, y0,
				static_cast<MML::Real>(1.0), config, eventConfig);
			REQUIRE(ok.integration.status == AlgorithmStatus::Success);
			REQUIRE(ok.events.size() == 1);
			REQUIRE_THAT(ok.events[0].time, WithinAbs(static_cast<MML::Real>(0.3), TOL(1e-9, 1e-4)));
			REQUIRE_THAT(ok.events[0].differential_state[0], WithinAbs(ok.events[0].time, TOL(1e-9, 1e-4)));
			REQUIRE_THAT(ok.events[0].event_value, WithinAbs(static_cast<MML::Real>(0.0), eventConfig.event_tolerance));
		}
	}

	TEST_CASE("DAE tiny-magnitude event functions are not lost to sign-product underflow",
		"[dae][events][MinimalMathLibrary-hrdm]") {
		if constexpr (!std::is_same_v<Real, float>) {
			class ScaledEventDAE : public TimedEventDAE {
			public:
				using TimedEventDAE::TimedEventDAE;
				Real eventFunction(int i, Real t, const Vector<Real>& x, const Vector<Real>& y) const override {
					return static_cast<MML::Real>(1e-200) * TimedEventDAE::eventFunction(i, t, x, y);
				}
			};
			ScaledEventDAE system(static_cast<MML::Real>(0.5), EventAction::Stop);
			Vector<Real> x0{static_cast<MML::Real>(0.0)}, y0{static_cast<MML::Real>(1.0)};
			DAESolverConfig config = EventTestConfig();
			config.step_size = static_cast<MML::Real>(1.0);
			DAEEventConfig eventConfig;
			eventConfig.event_tolerance = static_cast<MML::Real>(1e-210);

			auto result = SolveDAEBackwardEulerWithEvents(system, static_cast<MML::Real>(0.0), x0, y0,
				static_cast<MML::Real>(1.0), config, eventConfig);
			INFO(result.integration.error_message);
			REQUIRE(result.integration.status == AlgorithmStatus::Success);
			REQUIRE(result.terminated_by_event);
			REQUIRE(result.events.size() == 1);
			REQUIRE_THAT(result.events[0].time, WithinAbs(static_cast<MML::Real>(0.5), static_cast<MML::Real>(1e-9)));
		}
	}

} // namespace MML::Tests::Algorithms::DAEEventDetectionTests
