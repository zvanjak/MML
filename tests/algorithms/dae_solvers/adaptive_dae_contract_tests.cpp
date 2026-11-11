#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/DAESolvers/DAEAdaptive.h>
#include <mml/algorithms/DAESolvers/DAERODAS4.h>
#endif

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace MML;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace {
	class AdaptiveDecayDAE : public IODESystemDAEWithJacobian {
	public:
		Real constraintScale = 1;
		mutable int evaluations = 0;
		bool forbidEvaluation = false;
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void Count() const {
			++evaluations;
			if (forbidEvaluation) throw std::logic_error("Unexpected DAE evaluation");
		}
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
			Count();
			f[0] = -x[0];
		}
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			Count();
			g[0] = y[0] - constraintScale * x[0];
		}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -constraintScale; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
	};

	class AdaptiveQuadraticDAE : public AdaptiveDecayDAE {
	public:
		void diffEqs(Real t, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = t * t; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
	};

	class AdaptiveConstraintDAE : public AdaptiveDecayDAE {
	public:
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
		void algConstraints(Real t, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& g) const override {
			g[0] = y[0] * y[0] - (1 + t);
		}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>& y, Matrix<Real>& j) const override { j(0, 0) = 2 * y[0]; }
	};

	DAESolverConfig FixedAdaptiveConfig(Real step) {
		DAESolverConfig config;
		config.step_size = config.min_step_size = config.max_step_size = step;
		config.abs_tolerance = Precision::ZeroTolerance;
		config.rel_tolerance = 0;
		config.newton_tol = Precision::ZeroTolerance;
		return config;
	}

	void CheckAdaptiveHistory(const DAESolverResult& result) {
		CHECK(result.total_steps == result.accepted_steps);
		CHECK(result.solution.size() == result.accepted_steps + 1);
		CHECK(result.solution.getNumStepsOK() == result.accepted_steps);
		CHECK(result.solution.getNumStepsBad() == result.rejected_steps);
	}

	void CheckInitialOnly(const DAESolverResult& result, Real x = 1, Real y = 1) {
		CheckAdaptiveHistory(result);
		CHECK(result.accepted_steps == 0);
		REQUIRE(result.solution.size() == 1);
		CHECK(result.solution.getTValue(0) == 0);
		CHECK(result.solution.getXValuesAtEnd()[0] == x);
		CHECK(result.solution.getYValuesAtEnd()[0] == y);
		CHECK(result.final_constraint_norm == 0);
		CHECK_FALSE(result.error_message.empty());
	}
}

TEST_CASE("Adaptive DAE never forces an inaccurate minimum step into the trajectory", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	const auto config = FixedAdaptiveConfig(TOL(0.1, 0.5));
	for (Real direction : {Real(1), Real(-1)}) {
		CAPTURE(direction);
		auto result = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, direction * config.step_size, config);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(result.failure_reason == DAEFailureReason::StepSizeUnderflow);
		CHECK(result.rejected_steps == 1);
		CHECK(result.last_error_norm > 1);
		CHECK_THAT(result.error_message, ContainsSubstring("minimum step"));
		CheckInitialOnly(result);
		auto control = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{1}, direction * config.step_size, config);
		CHECK(control.status == AlgorithmStatus::NumericalInstability);
		CHECK(control.failure_reason == DAEFailureReason::StepSizeUnderflow);
		CHECK(control.rejected_steps == 1);
		CHECK(control.last_error_norm > 1);
		CheckInitialOnly(control);
	}
}

TEST_CASE("Adaptive DAE clipped terminal steps still require an acceptable error", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	for (Real direction : {Real(1), Real(-1)}) {
		auto config = FixedAdaptiveConfig(Real(0.125));
		const Real end = direction / 32;
		auto failed = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, end, config);
		CHECK(failed.status == AlgorithmStatus::NumericalInstability);
		CHECK(failed.failure_reason == DAEFailureReason::StepSizeUnderflow);
		CHECK(failed.rejected_steps == 1);
		CheckInitialOnly(failed);
		config.abs_tolerance = 1;
		auto accepted = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, end, config);
		CHECK(accepted.status == AlgorithmStatus::Success);
		CHECK(accepted.accepted_steps == 1);
		CHECK(accepted.rejected_steps == 0);
		CheckAdaptiveHistory(accepted);
		CHECK(accepted.solution.getTValue(1) == end);
		const Real expected = 1 / ((1 + end / 2) * (1 + end / 2));
		CHECK_THAT(accepted.solution.getXValuesAtEnd()[0], WithinAbs(expected, TOL(1e-12, 1e-5)));
	}
}

TEST_CASE("Adaptive DAE minimum-step failure retains earlier accepted states and counters", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveQuadraticDAE system;
	const Real step = Real(0.125);
	auto config = FixedAdaptiveConfig(step);
	config.abs_tolerance = step * step * step / 2;
	for (Real direction : {Real(1), Real(-1)}) {
		auto result = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{0}, Vector<Real>{0}, direction * 3 * step, config);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(result.failure_reason == DAEFailureReason::StepSizeUnderflow);
		CHECK(result.accepted_steps == 1);
		CHECK(result.rejected_steps == 1);
		CheckAdaptiveHistory(result);
		REQUIRE(result.solution.size() == 2);
		CHECK(result.solution.getTValue(1) == direction * step);
		// Two BE half steps integrate t^2 to (h/2)^3 + (h/2)*h^2 = 5*h^3/8.
		const Real expected = direction * 5 * step * step * step / 8;
		CHECK_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(expected, TOL(1e-12, 1e-5)));
		CHECK_THAT(result.solution.getYValuesAtEnd()[0], WithinAbs(expected, TOL(1e-12, 1e-5)));
		CHECK(result.last_error_norm > 1);
		CHECK(result.max_error_norm >= result.last_error_norm);
		CHECK(result.final_constraint_norm == 0);
	}
}

TEST_CASE("Adaptive DAE successful retries satisfy the estimator and the analytic decay control", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	auto config = FixedAdaptiveConfig(Real(0.1));
	config.min_step_size = Precision::DAEMinStepSize;
	config.abs_tolerance = Real(1e-5);
	for (Real direction : {Real(1), Real(-1)}) {
		auto result = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, direction, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		CHECK(result.rejected_steps > 0);
		CHECK(result.accepted_steps > 1);
		CHECK(result.max_error_norm > 1);
		CHECK(result.last_error_norm <= 1);
		CheckAdaptiveHistory(result);
		CHECK(result.solution.getTValue(result.solution.size() - 1) == direction);
		CHECK_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(std::exp(-direction), TOL(1e-2, 5e-2)));
		for (int i = 1; i < result.solution.size(); ++i) {
			const Real h = result.solution.getTValue(i) - result.solution.getTValue(i - 1);
			const Real start = result.solution.getXValue(i - 1, 0);
			const Real full = start / (1 + h);
			const Real halves = start / ((1 + h / 2) * (1 + h / 2));
			CHECK(direction * h > 0);
			CHECK_THAT(result.solution.getXValue(i, 0), WithinAbs(halves, TOL(1e-12, 1e-5)));
			CHECK(std::abs(halves - full) <= config.abs_tolerance * (1 + TOL(1e-7, 1e-3)));
			CHECK_THAT(result.solution.getYValue(i, 0), WithinAbs(halves, TOL(1e-12, 1e-5)));
		}
	}
}

TEST_CASE("Adaptive DAE counts accepted and rejected attempts against the budget", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	auto config = FixedAdaptiveConfig(Real(0.125));
	config.min_step_size = Precision::DAEMinStepSize;
	config.max_steps = 1;
	auto rejected = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	CHECK(rejected.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(rejected.failure_reason == DAEFailureReason::MaxStepsExceeded);
	CHECK(rejected.rejected_steps == 1);
	CheckInitialOnly(rejected);
	config.abs_tolerance = 1;
	config.max_steps = 2;
	auto accepted = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	CHECK(accepted.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(accepted.failure_reason == DAEFailureReason::MaxStepsExceeded);
	CHECK(accepted.accepted_steps == 2);
	CHECK(accepted.rejected_steps == 0);
	CheckAdaptiveHistory(accepted);
	CHECK(accepted.solution.getTValue(2) == Real(0.25));
}

TEST_CASE("Adaptive DAE rejects nonfinite error scales instead of treating them as zero error", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	auto config = FixedAdaptiveConfig(Real(0.125));
	config.abs_tolerance = config.rel_tolerance = std::numeric_limits<Real>::max();
	auto result = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{2}, Vector<Real>{2}, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.failure_reason == DAEFailureReason::NonFiniteState);
	CHECK(result.rejected_steps == 1);
	CheckInitialOnly(result, 2, 2);
}

TEST_CASE("Adaptive DAE algebraic error controls acceptance independently of differential error", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	auto config = FixedAdaptiveConfig(Real(0.125));
	config.abs_tolerance = 1;
	auto accepted = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, config.step_size, config);
	CHECK(accepted.status == AlgorithmStatus::Success);
	system.constraintScale = 1000;
	// Reach the truncation-error gate rather than demanding a sub-ULP algebraic residual.
	config.newton_tol = std::max(config.newton_tol, 4 * Precision::Eps * system.constraintScale);
	auto rejected = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1000}, config.step_size, config);
	CHECK(rejected.status == AlgorithmStatus::NumericalInstability);
	CHECK(rejected.failure_reason == DAEFailureReason::StepSizeUnderflow);
	CHECK(rejected.last_error_norm > 1);
	CHECK(rejected.rejected_steps == 1);
	CheckInitialOnly(rejected, 1, 1000);
}

TEST_CASE("Adaptive DAE rejected constraint candidates do not replace accepted-state diagnostics", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveConstraintDAE system;
	auto config = FixedAdaptiveConfig(Real(0.125));
	config.newton_tol = Real(1e-4);
	config.constraint_tol = Precision::ZeroTolerance;
	config.abs_tolerance = 1;
	auto result = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{0}, Vector<Real>{1}, config.step_size, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.failure_reason == DAEFailureReason::ConstraintViolation);
	CHECK(result.rejected_steps == 1);
	CheckInitialOnly(result, 0, 1);
	CHECK(result.max_constraint_violation > config.constraint_tol);
}

TEST_CASE("Adaptive DAE validates every numerical control before evaluating the model", "[dae][adaptive][adaptive-dae-contract]") {
	AdaptiveDecayDAE system;
	system.forbidEvaluation = true;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real infinity = std::numeric_limits<Real>::infinity();
	for (Real DAESolverConfig::* field : {&DAESolverConfig::step_size, &DAESolverConfig::min_step_size,
		&DAESolverConfig::max_step_size, &DAESolverConfig::abs_tolerance,
		&DAESolverConfig::newton_tol, &DAESolverConfig::constraint_tol}) {
		for (Real value : {Real(0), Real(-1), nan, infinity, -infinity}) {
			DAESolverConfig config;
			config.*field = value;
			CAPTURE(value);
			CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
		}
	}
	for (Real value : {Real(-1), nan, infinity, -infinity}) {
		DAESolverConfig config;
		config.rel_tolerance = value;
		CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
	}
	for (Real value : {Real(0), Real(-1), Real(1), Real(2), nan, infinity}) {
		DAESolverConfig config;
		config.step_safety = value;
		CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
	}
	for (int DAESolverConfig::* field : {&DAESolverConfig::max_steps, &DAESolverConfig::max_newton_iter}) {
		for (int value : {0, -1}) {
			DAESolverConfig config;
			config.*field = value;
			CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
		}
	}
	DAESolverConfig reversedBounds;
	reversedBounds.max_step_size = reversedBounds.min_step_size / 2;
	CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, reversedBounds), ArgumentError);
	CHECK(system.evaluations == 0);
}

TEST_CASE("Adaptive DAE rejects invalid intervals and dimensions before model evaluation", "[dae][adaptive][adaptive-dae-contract]") {
	class DimensionDAE : public AdaptiveDecayDAE {
	public:
		int differential = 1, algebraic = 1;
		int getDiffDim() const override { return differential; }
		int getAlgDim() const override { return algebraic; }
	};
	DimensionDAE system;
	system.forbidEvaluation = true;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real infinity = std::numeric_limits<Real>::infinity();
	const Real largest = std::numeric_limits<Real>::max();
	for (const auto& times : std::vector<std::pair<Real, Real>>{
		{Real(0), Real(0)}, {nan, Real(1)}, {Real(0), nan}, {infinity, Real(1)}, {Real(0), infinity}, {-largest, largest}}) {
		CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, times.first, Vector<Real>{1}, Vector<Real>{1}, times.second), ArgumentError);
	}
	for (const auto& dims : std::vector<std::array<int, 4>>{
		{0, 1, 0, 1}, {-1, 1, 0, 1}, {1, -1, 1, 0}, {1, 1, 0, 1},
		{1, 1, 1, 0}, {1, 1, 2, 1}, {1, 1, 1, 2}}) {
		system.differential = dims[0];
		system.algebraic = dims[1];
		Vector<Real> x(dims[2]), y(dims[3]);
		CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, x, y, 1), ArgumentError);
	}
	CHECK(system.evaluations == 0);
}

TEST_CASE("Adaptive DAE Newton failures retain partial history and initial failures remain explicit", "[dae][adaptive][adaptive-dae-contract]") {
	class LateFailureDAE : public AdaptiveDecayDAE {
	public:
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			if (t > Real(0.125)) f[0] = std::numeric_limits<Real>::quiet_NaN();
			else AdaptiveDecayDAE::diffEqs(t, x, y, f);
		}
	} late;
	auto config = FixedAdaptiveConfig(Real(0.125));
	config.abs_tolerance = 1;
	auto partial = SolveDAEBackwardEulerAdaptive(late, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	CHECK(partial.status == AlgorithmStatus::NumericalInstability);
	CHECK(partial.failure_reason == DAEFailureReason::NonFiniteState);
	CHECK(partial.accepted_steps == 1);
	CHECK(partial.rejected_steps == 1);
	CheckAdaptiveHistory(partial);
	REQUIRE(partial.solution.size() == 2);
	CHECK(partial.solution.getTValue(1) == Real(0.125));
	CHECK(partial.final_constraint_norm == 0);
	AdaptiveDecayDAE system;
	config.max_newton_iter = 1;
	AdaptiveConstraintDAE nonlinear;
	auto exhausted = SolveDAEBackwardEulerAdaptive(nonlinear, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	CHECK(exhausted.failure_reason == DAEFailureReason::NewtonFailure);
	CHECK(exhausted.rejected_steps == 1);
	CheckInitialOnly(exhausted);
	auto inconsistent = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{0}, 1);
	CHECK(inconsistent.status == AlgorithmStatus::InvalidInput);
	CHECK(inconsistent.failure_reason == DAEFailureReason::InconsistentInitialConditions);
	CHECK(inconsistent.rejected_steps == 0);
	CheckAdaptiveHistory(inconsistent);
	CHECK(inconsistent.final_constraint_norm == 1);
	auto nonfinite = SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{std::numeric_limits<Real>::quiet_NaN()}, Vector<Real>{1}, 1);
	CHECK(nonfinite.status == AlgorithmStatus::NumericalInstability);
	CHECK(nonfinite.failure_reason == DAEFailureReason::NonFiniteState);
	CheckAdaptiveHistory(nonfinite);
	system.forbidEvaluation = true;
	CHECK_THROWS_AS(SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1), std::logic_error);
}

TEST_CASE("Adaptive DAE supports zero algebraic dimensions and reports stalled clocks", "[dae][adaptive][adaptive-dae-contract]") {
	class ODEOnlyDAE : public AdaptiveDecayDAE {
	public:
		int getAlgDim() const override { return 0; }
		void algConstraints(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&) const override {}
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
	} ode;
	auto config = FixedAdaptiveConfig(Real(0.125));
	config.abs_tolerance = 1;
	for (Real initial : {config.min_step_size / 2, config.max_step_size * 2}) {
		config.step_size = initial;
		auto result = SolveDAEBackwardEulerAdaptive(ode, 0, Vector<Real>{1}, Vector<Real>{}, config.max_step_size, config);
		CHECK(result.status == AlgorithmStatus::Success);
		CHECK(result.accepted_steps == 1);
		CHECK(result.solution.getAlgDim() == 0);
		CheckAdaptiveHistory(result);
	}
	AdaptiveDecayDAE system;
	const Real clock = 1 / std::numeric_limits<Real>::epsilon();
	auto stalled = SolveDAEBackwardEulerAdaptive(system, clock, Vector<Real>{1}, Vector<Real>{1}, clock - 4, config);
	CHECK(stalled.status == AlgorithmStatus::Stalled);
	CHECK(stalled.failure_reason == DAEFailureReason::StepSizeUnderflow);
	CHECK(stalled.accepted_steps == 0);
	CHECK(stalled.rejected_steps == 0);
	CheckAdaptiveHistory(stalled);
	CHECK(stalled.solution.getTValue(0) == clock);
}
