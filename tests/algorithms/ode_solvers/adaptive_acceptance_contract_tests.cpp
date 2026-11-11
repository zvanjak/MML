#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
	class AcceptanceODE : public IODESystem {
	public:
		Real value = 1;
		mutable std::vector<Real> attempts;
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>&, Vector<Real>& rhs) const override { rhs[0] = value; }
	};

	class StationaryZeroODE : public IODESystem {
	public:
		mutable int calls = 0;
		mutable Real minimumTime = std::numeric_limits<Real>::infinity();
		mutable Real maximumTime = -std::numeric_limits<Real>::infinity();
		int getDim() const override { return 1; }
		void derivs(Real time, const Vector<Real>&, Vector<Real>& rhs) const override {
			++calls;
			minimumTime = std::min(minimumTime, time);
			maximumTime = std::max(maximumTime, time);
			rhs[0] = time;
		}
	};

	class RecordingStepper {
		const AcceptanceODE& _system;
		bool _reject = true;
	public:
		explicit RecordingStepper(const IODESystem& system) : _system(dynamic_cast<const AcceptanceODE&>(system)) {}
		void resetFSAL() { _reject = true; }
		StepResult doStep(Real, Vector<Real>& x, Vector<Real>& derivative, Real h, Real) {
			_system.attempts.push_back(h);
			StepResult result;
			result.hDone = h;
			result.hNext = _reject ? h / 2 : h * 8;
			result.accepted = !_reject;
			if (result.accepted) {
				x[0] += h;
				derivative[0] = 1;
			}
			_reject = false;
			return result;
		}
		Vector<Real> interpolate(Real time) const { return Vector<Real>{time}; }
	};

	class StageFailureODE : public AcceptanceODE {
	public:
		mutable int calls = 0;
		mutable std::vector<Real> evaluationTimes;
		int failAt = 0;
		Real invalid = std::numeric_limits<Real>::quiet_NaN();
		mutable bool invalidInput = false;
		void derivs(Real time, const Vector<Real>& state, Vector<Real>& rhs) const override {
			evaluationTimes.push_back(time);
			invalidInput = invalidInput || !std::isfinite(time) || !std::isfinite(state[0]);
			rhs[0] = ++calls == failAt ? invalid : Real(1);
		}
	};

	class InvalidResultODE : public AcceptanceODE {
	public:
		int fault = 0;
	};

	class InvalidResultStepper {
		const InvalidResultODE& _system;
	public:
		explicit InvalidResultStepper(const IODESystem& system) : _system(dynamic_cast<const InvalidResultODE&>(system)) {}
		void resetFSAL() {}
		StepResult doStep(Real, Vector<Real>& x, Vector<Real>& derivative, Real h, Real) {
			const Real nan = std::numeric_limits<Real>::quiet_NaN();
			StepResult result;
			result.accepted = true;
			result.hDone = h;
			result.hNext = h;
			x[0] += h;
			switch (_system.fault) {
			case 0: result.errMax = nan; break;
			case 1: result.hDone = nan; break;
			case 2: result.hNext = nan; break;
			case 3: x[0] = nan; break;
			case 4: derivative[0] = nan; break;
			case 5: result.hDone = 2 * h; break;
			case 7: result.errMax = 2; break;
			case 8: x.resize(0); break;
			}
			return result;
		}
		Vector<Real> interpolate(Real time) const {
			return Vector<Real>{_system.fault == 6 ? std::numeric_limits<Real>::quiet_NaN() : time};
		}
	};

	template<class Test>
	void ForEachStepper(Test test) {
		test.template operator()<DormandPrince5_Stepper>("DP5");
		test.template operator()<CashKarp_Stepper>("Cash-Karp");
		test.template operator()<DormandPrince8_Stepper>("DOP853");
		test.template operator()<BulirschStoer_Stepper>("Bulirsch-Stoer");
		test.template operator()<BulirschStoerRational_Stepper>("Bulirsch-Stoer alternate sequence");
	}
}

TEST_CASE("Adaptive acceptance rejects a NaN RHS rather than certifying a trajectory", "[ode-acceptance-baseline][ode-acceptance]") {
	AcceptanceODE system;
	system.value = std::numeric_limits<Real>::quiet_NaN();
	ODEAdaptiveConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	config.integrator.initial_step_size = Real(0.125);
	config.integrator.output_interval = Real(0.125);
	const auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{0}, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(result.IsSuccess());
	CHECK_FALSE(result.error_message.empty());
}

TEST_CASE("Adaptive RK error scales retain history and resolve higher-order endpoint growth", "[ode-acceptance]") {
	for (Real step : {Real(-0.125), Real(0.125)}) {
		const Real history = 2 + std::abs(3 * step) + Precision::DivisionSafetyThreshold;
		CHECK(ODEStepperDetail::RKErrorScale(2, 3, step, 1) == history);
		CHECK(ODEStepperDetail::RKErrorScale(2, 3, step, -4) == 4);
		CHECK(ODEStepperDetail::RKErrorScale(0, 0, step, step * step / 2) == step * step / 2);
	}
	for (Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
		CHECK_THROWS_AS(ODEStepperDetail::RKErrorScale(0, 0, Real(0.125), invalid), ODESolverError);
	CHECK_THROWS_AS(ODEStepperDetail::RKErrorScale(std::numeric_limits<Real>::max(),
		std::numeric_limits<Real>::max(), 1, 0), ODESolverError);
}

TEST_CASE("Adaptive steppers accept quadratic emergence from a stationary zero", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real step : {Real(-0.1), Real(0.1)}) {
			CAPTURE(step);
			StationaryZeroODE system;
			Stepper stepper(system);
			Vector<Real> state{0}, derivative{0};
			const auto result = stepper.doStep(0, state, derivative, step, Precision::ODEDefaultTolerance);
			REQUIRE(result.accepted);
			CHECK(result.hDone == step);
			CHECK(result.errMax <= 1);
			CHECK_THAT(state[0], WithinRel(step * step / 2, Precision::EqualityTolerance));
			CHECK(derivative[0] == step);
			CHECK(system.calls > 0);
			CHECK(result.funcEvals == system.calls);
			CHECK(system.minimumTime >= std::min(Real(0), step));
			CHECK(system.maximumTime <= std::max(Real(0), step));
		}
	});
}

TEST_CASE("Adaptive trajectories preserve quadratic clocks from a stationary zero in both directions", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real end : {Real(-1), Real(1)}) {
			CAPTURE(end);
			StationaryZeroODE system;
			ODEAdaptiveConfig config;
			config.integrator.initial_step_size = Real(0.1);
			config.integrator.output_interval = Real(0.125);
			const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, end, config);
			REQUIRE(result.IsSuccess());
			REQUIRE(result.solution.has_value());
			const auto& solution = *result.solution;
			REQUIRE(solution.size() == 9);
			CHECK(system.calls > 0);
			CHECK(result.function_evaluations == system.calls);
			CHECK(system.minimumTime >= std::min(Real(0), end));
			CHECK(system.maximumTime <= std::max(Real(0), end));
			for (int point = 0; point < solution.size(); ++point) {
				const Real time = end * point / 8;
				CHECK(solution.getTValue(point) == time);
				CHECK_THAT(solution.getXValue(point, 0), WithinAbs(time * time / 2, Precision::EqualityTolerance));
			}
		}
	});
}

TEST_CASE("Direct DP5 rejects nonfinite stage derivatives without changing state", "[ode-acceptance-baseline][ode-acceptance]") {
	AcceptanceODE system;
	system.value = std::numeric_limits<Real>::quiet_NaN();
	DormandPrince5_Stepper stepper(system);
	Vector<Real> state{0}, derivative{1};
	CHECK_THROWS_AS(stepper.doStep(0, state, derivative, Real(0.125), Precision::ODEDefaultTolerance), ODESolverError);
	CHECK(state[0] == Real(0));
	CHECK(derivative[0] == Real(1));
}

TEST_CASE("Adaptive maximum step bounds initial attempts retries and growth in either direction", "[ode-acceptance-baseline][ode-acceptance]") {
	for (const Real direction : {Real(-1), Real(1)}) {
		AcceptanceODE system;
		ODEAdaptiveIntegrator<RecordingStepper> integrator(system);
		ODEIntegratorConfig config;
		config.initial_step_size = Real(0.5);
		config.max_step_size = Real(0.125);
		config.output_interval = Real(0.125);
		const auto solution = integrator.integrate(Vector<Real>{0}, 0, direction, config);
		REQUIRE_FALSE(system.attempts.empty());
		for (const Real step : system.attempts) {
			CHECK(std::abs(step) <= config.max_step_size);
			CHECK(step * direction > 0);
		}
		CHECK(integrator.getStatistics().maxStepSize <= config.max_step_size);
		CHECK(integrator.getStatistics().rejectedSteps == 1);
		CHECK_THAT(solution.getXValuesAtEnd()[0], WithinAbs(direction, Precision::EqualityTolerance));
	}
}

TEST_CASE("Every adaptive stepper rejects each nonfinite RHS stage before changing caller state", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		StageFailureODE control;
		Stepper controlStepper(control);
		Vector<Real> state{0}, derivative{1};
		REQUIRE(controlStepper.doStep(0, state, derivative, Real(0.125), Precision::ODEDefaultTolerance).accepted);
		REQUIRE(control.calls > 0);
		for (int stage = 1; stage <= control.calls; ++stage) {
			for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
				std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
				INFO("callback stage " << stage);
				StageFailureODE system;
				system.failAt = stage;
				system.invalid = invalid;
				Stepper stepper(system);
				Vector<Real> x{0}, dx{1};
				CHECK_THROWS_AS(stepper.doStep(0, x, dx, Real(0.125), Precision::ODEDefaultTolerance), ODESolverError);
				CHECK(system.calls == stage);
				CHECK_FALSE(system.invalidInput);
				CHECK(x[0] == Real(0));
				CHECK(dx[0] == Real(1));
			}
		}
	});
}

TEST_CASE("Adaptive finite safety is mandatory under both Detailed check modes and exception policies", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (bool check : {false, true}) {
			for (Real initialStep : {Real(0), Real(0.125)}) {
				for (int failAt : {1, 3}) {
					StageFailureODE system;
					system.failAt = failAt;
					ODEAdaptiveConfig config;
					config.check_finite = check;
					config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
					config.integrator.initial_step_size = initialStep;
					config.integrator.output_interval = Real(0.125);
					auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, 1, config);
					CHECK(result.status == AlgorithmStatus::NumericalInstability);
					CHECK_FALSE(result.IsSuccess());
					CHECK_FALSE(result.error_message.empty());
					CHECK_FALSE(system.invalidInput);
					system.calls = 0;
					config.exception_policy = EvaluationExceptionPolicy::Propagate;
					CHECK_THROWS_AS(ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, 1, config), ODESolverError);
				}
			}
		}
	});
}

TEST_CASE("Adaptive integration validates initial values and bounds before invoking callbacks", "[ode-acceptance]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		StageFailureODE system;
		ODEAdaptiveConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{invalid}, 0, 1, config);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(system.calls == 0);
		config.integrator.max_step_size = invalid;
		if (invalid != std::numeric_limits<Real>::infinity()) {
			result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{0}, 0, 1, config);
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			CHECK(system.calls == 0);
		}
	}
	for (Real maximum : {Real(0), Real(-1), Precision::ODEMinStepSize / 2}) {
		StageFailureODE system;
		ODEAdaptiveConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.integrator.max_step_size = maximum;
		CHECK(ODEAdaptiveIntegrateDetailed(system, Vector<Real>{0}, 0, 1, config).status ==
			AlgorithmStatus::NumericalInstability);
		CHECK(system.calls == 0);
	}
	StageFailureODE system;
	DormandPrince5Integrator integrator(system);
	CHECK_THROWS_AS(integrator.integrate(Vector<Real>{0, 1}, 0, 1, Real(0.125)), ODESolverError);
	CHECK_THROWS_AS(integrator.integrateAt(Vector<Real>{0, 1}, Vector<Real>{0, 1}), ODESolverError);
	CHECK(system.calls == 0);
}

TEST_CASE("Direct adaptive steps reject invalid inputs and nonfinite arithmetic", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		AcceptanceODE system;
		for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity()}) {
			Stepper stepper(system);
			Vector<Real> state{0}, derivative{1};
			CHECK_THROWS_AS(stepper.doStep(0, state, derivative, Real(0.125), invalid), ODESolverError);
		}
		for (Real invalid : {Real(0), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity()}) {
			Stepper stepper(system);
			Vector<Real> state{0}, derivative{1};
			CHECK_THROWS_AS(stepper.doStep(0, state, derivative, invalid, Precision::ODEDefaultTolerance), ODESolverError);
		}
		Stepper stepper(system);
		Vector<Real> state{0}, derivative{std::numeric_limits<Real>::quiet_NaN()};
		CHECK_THROWS_AS(stepper.doStep(0, state, derivative, Real(0.125), Precision::ODEDefaultTolerance), ODESolverError);
		system.value = std::numeric_limits<Real>::max();
		derivative[0] = system.value;
		CHECK_THROWS_AS(stepper.doStep(0, state, derivative, Real(8), Precision::ODEDefaultTolerance), ODESolverError);
		CHECK(state[0] == Real(0));
		CHECK(derivative[0] == system.value);
	});
}

TEST_CASE("Built-in adaptive methods enforce the configured cap with explicit and automatic startup", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			for (Real initialStep : {Real(0), Real(0.5)}) {
				AcceptanceODE system;
				ODEAdaptiveConfig config;
				config.integrator.initial_step_size = initialStep;
				config.integrator.max_step_size = Real(0.125);
				config.integrator.output_interval = Real(0.125);
				const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, direction, config);
				REQUIRE(result.IsSuccess());
				REQUIRE(result.solution.has_value());
				CHECK(result.max_step_size <= config.integrator.max_step_size);
				CHECK(result.min_step_size > 0);
				for (int index = 0; index < result.solution->size(); ++index) {
					CHECK_THAT(result.solution->getXValue(index, 0),
						WithinAbs(result.solution->getTValue(index), Precision::EqualityTolerance));
				}
			}
		}
	});
}

TEST_CASE("Adaptive legacy paths and stepper reuse retain finite failure closure", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		StageFailureODE system;
		Stepper stepper(system);
		Vector<Real> state{0}, derivative{1};
		REQUIRE(stepper.doStep(0, state, derivative, Real(0.125), Precision::ODEDefaultTolerance).accepted);
		const Real saved = state[0];
		system.failAt = system.calls + 1;
		CHECK_THROWS_AS(stepper.doStep(Real(0.125), state, derivative, Real(0.125), Precision::ODEDefaultTolerance), ODESolverError);
		CHECK(state[0] == saved);
		CHECK(derivative[0] == Real(1));
		system.failAt = 0;
		REQUIRE(stepper.doStep(Real(0.125), state, derivative, Real(0.125), Precision::ODEDefaultTolerance).accepted);
		CHECK_THAT(state[0], WithinAbs(Real(0.25), Precision::EqualityTolerance));
		system.calls = 0;
		system.failAt = 2;
		ODEAdaptiveIntegrator<Stepper> integrator(system);
		CHECK_THROWS_AS(integrator.integrate(Vector<Real>{0}, 0, 1, Real(0.125), Precision::ODEDefaultTolerance, Real(0.125)), ODESolverError);
		system.calls = 0;
		CHECK_THROWS_AS(integrator.integrateAt(Vector<Real>{0}, Vector<Real>{0, Real(0.5), 1},
			Precision::ODEDefaultTolerance, Real(0.125)), ODESolverError);
	});
}

TEST_CASE("Adaptive error normalization and control cannot turn nonfinite values into zero error", "[ode-acceptance]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(ODEStepperDetail::ErrorRatio(invalid, Real(1)), ODESolverError);
		CHECK_THROWS_AS(ODEStepperDetail::ErrorRatio(Real(1), invalid), ODESolverError);
		StepSizeController controller;
		CHECK_THROWS_AS(controller.acceptedStep(Real(0.125), invalid), ODESolverError);
		CHECK_THROWS_AS(controller.rejectedStep(Real(0.125), invalid), ODESolverError);
	}
	const Real underflow = std::numeric_limits<Real>::min();
	CHECK_THROWS_AS(ODEStepperDetail::ErrorRatio(std::numeric_limits<Real>::max(), underflow), ODESolverError);
	CHECK_THROWS_AS(ODEStepperDetail::ErrorRatio(Real(1), Real(0)), ODESolverError);
	CHECK_THROWS_AS(ODEStepperDetail::ErrorRatio(Real(1), Real(-1)), ODESolverError);
	StepSizeController controller;
	CHECK_THROWS_AS(controller.acceptedStep(std::numeric_limits<Real>::max() / 2, Real(0)), ODESolverError);
}

TEST_CASE("Adaptive automatic startup probes honor a cap below the uncapped Euler estimate", "[ode-acceptance]") {
	ForEachStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			StageFailureODE system;
			ODEIntegratorConfig config;
			config.max_step_size = Real(1) / 1024;
			config.output_interval = Real(1) / 64;
			ODEAdaptiveIntegrator<Stepper> integrator(system);
			const auto solution = integrator.integrate(Vector<Real>{1}, 0, direction / 16, config);
			REQUIRE(system.evaluationTimes.size() >= 3);
			CHECK(system.evaluationTimes[0] == Real(0));
			CHECK(system.evaluationTimes[1] == Real(0));
			CHECK(std::abs(system.evaluationTimes[2]) <= config.max_step_size);
			CHECK(direction * system.evaluationTimes[2] > 0);
			CHECK(integrator.getStatistics().maxStepSize <= config.max_step_size);
			CHECK_THAT(solution.getXValuesAtEnd()[0], WithinAbs(1 + direction / 16, Precision::EqualityTolerance));
		}
	});
}

TEST_CASE("Adaptive driver checks custom stepper acceptance and dense payloads before storing them", "[ode-acceptance]") {
	for (int fault = 0; fault <= 8; ++fault) {
		INFO("fault " << fault);
		InvalidResultODE system;
		system.fault = fault;
		ODEAdaptiveConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.integrator.initial_step_size = Real(0.5);
		config.integrator.output_interval = Real(0.125);
		const auto result = ODEAdaptiveIntegrateDetailed<InvalidResultStepper>(system, Vector<Real>{0}, 0, 1, config);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK_FALSE(result.IsSuccess());
		CHECK_FALSE(result.error_message.empty());
	}
}
