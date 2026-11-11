#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class OutputCallbackError : public std::runtime_error {
	public:
		OutputCallbackError() : std::runtime_error("output callback failure") {}
	};

	class OutputODE : public IODESystem {
	public:
		Real slope = 1;
		int faultAttempt = 0;
		int fault = 0;
		bool rejectFirst = false;
		mutable std::vector<Real> acceptedTimes;
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>&, Vector<Real>& rhs) const override { rhs[0] = slope; }
	};

	class OutputStepper {
		const OutputODE& _system;
		int _attempt = 0;
		Real _time = 0, _state = 0;
	public:
		explicit OutputStepper(const IODESystem& system) : _system(dynamic_cast<const OutputODE&>(system)) {
			if (_system.fault == 5) throw OutputCallbackError();
		}
		void resetFSAL() { _attempt = 0; }
		StepResult doStep(Real time, Vector<Real>& state, Vector<Real>&, Real h, Real) {
			++_attempt;
			_time = time;
			_state = state[0];
			StepResult result;
			result.accepted = !(_system.rejectFirst && _attempt == 1);
			result.hDone = h;
			result.hNext = h;
			result.funcEvals = 3;
			if (_attempt == _system.faultAttempt) {
				state[0] = std::numeric_limits<Real>::quiet_NaN();
				if (_system.fault == 0) throw OutputCallbackError();
				if (_system.fault == 1) throw ResourceLimitError("injected storage limit");
				if (_system.fault == 2) {
					state[0] = _state;
					result.hDone = std::copysign(std::numeric_limits<Real>::denorm_min(), h);
				}
			} else if (result.accepted) {
				state[0] += h * _system.slope;
				_system.acceptedTimes.push_back(time + h);
			} else state[0] = std::numeric_limits<Real>::quiet_NaN();
			return result;
		}
		Vector<Real> interpolate(Real time) const {
			if (_system.fault == 4) throw OutputCallbackError();
			return Vector<Real>{_state + (time - _time) * _system.slope};
		}
	};

	class CallbackOutputODE : public OutputODE {
	public:
		Real direction = 1;
		int failAtCall = 0;
		mutable int calls = 0;
		void derivs(Real time, const Vector<Real>&, Vector<Real>& rhs) const override {
			++calls;
			const bool fail = failAtCall > 0 ? calls == failAtCall : direction * time > Real(0.25);
			if (fail && fault == 0) throw OutputCallbackError();
			if (fail && fault == 1) throw ResourceLimitError("callback resource failure");
			rhs[0] = fail ? std::numeric_limits<Real>::quiet_NaN() : Real(1);
		}
	};

	template<class Test>
	void ForEachOutputStepper(Test test) {
		test.template operator()<DormandPrince5_Stepper>("DP5");
		test.template operator()<CashKarp_Stepper>("Cash-Karp");
		test.template operator()<DormandPrince8_Stepper>("DOP853");
		test.template operator()<BulirschStoer_Stepper>("Bulirsch-Stoer");
		test.template operator()<BulirschStoerRational_Stepper>("Bulirsch-Stoer alternate sequence");
	}

	void CheckOutput(const ODESystemSolution& solution, Real direction, Real slope = 1) {
		REQUIRE(solution.size() > 0);
		for (int i = 0; i < solution.size(); ++i) {
			CHECK(std::isfinite(solution.getTValue(i)));
			CHECK_THAT(solution.getXValue(i, 0),
				WithinAbs(solution.getTValue(i) * slope, Precision::EqualityTolerance));
			if (i > 0) CHECK(direction * (solution.getTValue(i) - solution.getTValue(i - 1)) > 0);
		}
	}
}

TEST_CASE("Adaptive output rejects an unrepresentable fixed clock and retains the initial point", "[ode-output-baseline][ode-output]") {
	ForEachOutputStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			OutputODE system;
			system.slope = 0;
			ODEAdaptiveConfig config;
			config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
			config.integrator.initial_step_size = 2;
			config.integrator.output_interval = Real(0.1);
			const Real start = std::ldexp(Real(1), std::numeric_limits<Real>::digits);
			const Real end = start + direction * 2;
			REQUIRE(std::isfinite(end));
			REQUIRE(direction * (end - start) > 0);
			REQUIRE(start + direction * config.integrator.output_interval == start);
			auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{3}, start, end, config);
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			REQUIRE(result.solution.has_value());
			REQUIRE(result.solution->size() == 1);
			CHECK(result.solution->getTValue(0) == start);
			CHECK(result.solution->getXValue(0, 0) == Real(3));
			CHECK(result.accepted_steps == 0);
			config.integrator.output_interval = 2;
			result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{3}, start, end, config);
			REQUIRE(result.IsSuccess());
			REQUIRE(result.solution->size() == 2);
			CHECK(result.solution->getTValue(1) == end);
		}
	});
}

TEST_CASE("Zero output interval stores accepted integration steps and excludes rejected trials", "[ode-output-baseline][ode-output]") {
	for (Real direction : {Real(-1), Real(1)}) {
		OutputODE system;
		system.rejectFirst = true;
		ODEAdaptiveIntegrator<OutputStepper> integrator(system);
		ODEIntegratorConfig config;
		config.initial_step_size = Real(0.25);
		const auto solution = integrator.integrate(Vector<Real>{0}, 0, direction, config);
		REQUIRE(solution.size() == static_cast<int>(system.acceptedTimes.size()) + 1);
		for (int i = 1; i < solution.size(); ++i)
			CHECK(solution.getTValue(i) == system.acceptedTimes[i - 1]);
		CHECK(solution.getNumStepsOK() == 4);
		CHECK(solution.getNumStepsBad() == 1);
		CheckOutput(solution, direction);
	}
}

TEST_CASE("Explicit tiny output times keep their exact clocks and corresponding states", "[ode-output-baseline][ode-output]") {
	ForEachOutputStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			OutputODE system;
			system.slope = Real(1e18);
			const Vector<Real> times{0, direction * Real(1e-18), direction * Real(2e-18), direction * Real(1e-17)};
			ODEAdaptiveIntegrator<Stepper> integrator(system);
			const auto solution = integrator.integrateAt(Vector<Real>{0}, times,
				Precision::ODEDefaultTolerance, Real(1e-18), 0);
			REQUIRE(solution.size() == times.size());
			for (int i = 0; i < times.size(); ++i) {
				CHECK(solution.getTValue(i) == times[i]);
				const Real expected = times[i] * system.slope;
				CHECK_THAT(solution.getXValue(i, 0), WithinAbs(expected,
					Precision::EqualityTolerance * std::max(Real(1), std::abs(expected))));
			}
		}
	});
}

TEST_CASE("Converted adaptive failures retain the last accepted off-grid endpoint and completed counters", "[ode-output-baseline][ode-output]") {
	for (Real direction : {Real(-1), Real(1)}) {
		for (int fault : {0, 1, 3}) {
			OutputODE system;
			system.faultAttempt = 2;
			system.fault = fault;
			ODEAdaptiveConfig config;
			config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
			config.integrator.initial_step_size = Real(0.25);
			config.integrator.output_interval = Real(0.5);
			const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, direction, config);
			CHECK(result.status == (fault == 0 ? AlgorithmStatus::AlgorithmSpecificFailure :
				fault == 1 ? AlgorithmStatus::ResourceLimitExceeded : AlgorithmStatus::NumericalInstability));
			REQUIRE(result.solution.has_value());
			REQUIRE(result.solution->size() == 2);
			CHECK(result.solution->getTValue(1) == direction * Real(0.25));
			CHECK(result.accepted_steps == 1);
			CHECK(result.rejected_steps == 0);
			CHECK(result.total_func_evals == (fault == 3 ? 7 : 4));
			CHECK(result.solution->getNumStepsOK() == 1);
			CheckOutput(*result.solution, direction);
		}
	}
}

TEST_CASE("Adaptive dense output appends nondividing endpoints exactly once", "[ode-output-baseline][ode-output]") {
	ForEachOutputStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			OutputODE system;
			ODEAdaptiveIntegrator<Stepper> integrator(system);
			const auto solution = integrator.integrate(Vector<Real>{0}, 0, direction, Real(0.3),
				Precision::ODEDefaultTolerance, Real(0.25));
			REQUIRE(solution.size() == 5);
			CHECK(solution.getTValue(4) == direction);
			CheckOutput(solution, direction);
		}
	});
}

TEST_CASE("All adaptive methods retain accepted endpoints on real callback failures", "[ode-output]") {
	ForEachOutputStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			for (int fault : {0, 1, 3}) {
				CallbackOutputODE system;
				system.direction = direction;
				system.fault = fault;
				ODEAdaptiveConfig config;
				config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
				config.integrator.initial_step_size = Real(0.25);
				config.integrator.max_step_size = Real(0.25);
				config.integrator.output_interval = Real(0.5);
				const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, direction, config);
				CHECK(result.status == (fault == 0 ? AlgorithmStatus::AlgorithmSpecificFailure :
					fault == 1 ? AlgorithmStatus::ResourceLimitExceeded : AlgorithmStatus::NumericalInstability));
				CHECK_FALSE(result.error_message.empty());
				REQUIRE(result.solution.has_value());
				REQUIRE(result.solution->size() == 2);
				CHECK(result.solution->getTValue(1) == direction * Real(0.25));
				CHECK(result.accepted_steps == 1);
				CHECK(result.solution->getNumStepsOK() == 1);
				CHECK(result.total_func_evals > 0);
				CHECK(result.function_evaluations == result.total_func_evals);
				CHECK(result.function_evaluations == system.calls);
				CHECK(result.min_step_size == Real(0.25));
				CHECK(result.max_step_size == Real(0.25));
				CheckOutput(*result.solution, direction);
			}
		}
	});
}

TEST_CASE("Initial RHS and automatic estimator failures retain only valid initial data", "[ode-output]") {
	for (int call : {1, 2, 3}) {
		for (int fault : {0, 1, 3}) {
			CallbackOutputODE system;
			system.failAtCall = call;
			system.fault = fault;
			ODEAdaptiveConfig config;
			config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
			const auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{0}, 0, 1, config);
			CHECK_FALSE(result.IsSuccess());
			REQUIRE(result.solution.has_value());
			REQUIRE(result.solution->size() == 1);
			CHECK(result.solution->getTValue(0) == Real(0));
			CHECK(result.solution->getXValue(0, 0) == Real(0));
			CHECK(result.accepted_steps == 0);
			CHECK(result.rejected_steps == 0);
			CHECK(result.min_step_size == Real(0));
			CHECK(result.function_evaluations == system.calls);
		}
	}
}

TEST_CASE("Adaptive exception policies preserve callback types and legacy resource contracts", "[ode-output]") {
	OutputODE system;
	system.faultAttempt = 2;
	ODEAdaptiveConfig config;
	config.integrator.initial_step_size = Real(0.25);
	config.integrator.output_interval = Real(0.5);
	CHECK_THROWS_AS(ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config), OutputCallbackError);
	ODEAdaptiveIntegrator<OutputStepper> integrator(system);
	CHECK_THROWS_AS(integrator.integrate(Vector<Real>{0}, 0, 1, config.integrator), OutputCallbackError);
	system.fault = 3;
	CHECK_THROWS_AS(ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config), ODESolverError);
	system.fault = 1;
	const auto partial = integrator.integrate(Vector<Real>{0}, 0, 1, config.integrator);
	REQUIRE(partial.size() == 2);
	CHECK(partial.getTValue(1) == Real(0.25));
	CHECK(integrator.getStatistics().status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK_THROWS_AS(integrator.integrate(Vector<Real>{0}, 0, 1, Real(0.5),
		Precision::ODEDefaultTolerance, Real(0.25)), ResourceLimitError);
	const auto detailedResource = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config);
	CHECK(detailedResource.status == AlgorithmStatus::ResourceLimitExceeded);
	REQUIRE(detailedResource.solution.has_value());
	CHECK(detailedResource.solution->getTValue(1) == Real(0.25));
}

TEST_CASE("Adaptive attempt budgets retain off-grid endpoints and rejection counters", "[ode-output]") {
	for (Real direction : {Real(-1), Real(1)}) {
		OutputODE system;
		system.rejectFirst = true;
		ODEAdaptiveConfig config;
		config.integrator.initial_step_size = Real(0.25);
		config.integrator.output_interval = Real(0.5);
		config.integrator.max_steps = 2;
		const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, direction, config);
		CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(result.accepted_steps == 1);
		CHECK(result.rejected_steps == 1);
		CHECK(result.total_func_evals == 7);
		REQUIRE(result.solution.has_value());
		REQUIRE(result.solution->size() == 2);
		CHECK(result.solution->getNumStepsOK() == 1);
		CHECK(result.solution->getNumStepsBad() == 1);
		CHECK(result.solution->getTValue(1) == direction * Real(0.25));
		CheckOutput(*result.solution, direction);
	}
}

TEST_CASE("Dense interpolation failures retain valid samples and the accepted endpoint", "[ode-output]") {
	OutputODE system;
	system.fault = 4;
	ODEAdaptiveConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	config.integrator.initial_step_size = Real(0.25);
	config.integrator.output_interval = Real(0.125);
	const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
	REQUIRE(result.solution.has_value());
	REQUIRE(result.solution->size() == 2);
	CHECK(result.solution->getTValue(1) == Real(0.25));
	CHECK(result.accepted_steps == 1);
	CheckOutput(*result.solution, 1);
}

TEST_CASE("Fixed output clocks stay origin anchored and include exact divisible endpoints", "[ode-output]") {
	for (Real direction : {Real(-1), Real(1)}) {
		OutputODE system;
		ODEAdaptiveIntegrator<OutputStepper> integrator(system);
		const Real spacing = Real(0.1);
		const auto solution = integrator.integrate(Vector<Real>{0}, 0, direction * 10, spacing,
			Precision::ODEDefaultTolerance, Real(0.25));
		REQUIRE(solution.size() == 101);
		for (int i = 1; i < 100; ++i)
			CHECK(solution.getTValue(i) == std::fma(direction * spacing, Real(i), Real(0)));
		CHECK(solution.getTValue(100) == direction * 10);
		CheckOutput(solution, direction);
	}
}

TEST_CASE("Adaptive output validates spacing and checks actual accepted clock progress", "[ode-output]") {
	OutputODE system;
	ODEAdaptiveConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		config.integrator.output_interval = invalid;
		const auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{0}, 0, 1, config);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK_FALSE(result.solution.has_value());
		CHECK(result.accepted_steps == 0);
	}
	system.fault = 2;
	system.faultAttempt = 1;
	config.integrator.output_interval = 0;
	config.integrator.initial_step_size = Real(0.25);
	const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{1}, 1, 2, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.accepted_steps == 0);
	REQUIRE(result.solution.has_value());
	REQUIRE(result.solution->size() == 1);
	CHECK(result.solution->getTValue(0) == Real(1));
}

TEST_CASE("All built-in adaptive methods save exactly one point per accepted step at zero interval", "[ode-output]") {
	ForEachOutputStepper([]<class Stepper>(const char* name) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			OutputODE system;
			ODEAdaptiveConfig config;
			config.integrator.initial_step_size = Real(0.25);
			config.integrator.max_step_size = Real(0.25);
			const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, direction, config);
			REQUIRE(result.IsSuccess());
			REQUIRE(result.solution.has_value());
			CHECK(result.solution->size() == result.accepted_steps + 1);
			CHECK(result.solution->getNumStepsOK() == result.accepted_steps);
			CHECK(result.solution->getNumStepsBad() == result.rejected_steps);
			CHECK(result.solution->getTValue(result.solution->size() - 1) == direction);
			CheckOutput(*result.solution, direction);
		}
	});
}

TEST_CASE("A fixed grid that later loses representability fails without duplicate samples", "[ode-output]") {
	OutputODE system;
	system.slope = 0;
	ODEAdaptiveConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	config.integrator.initial_step_size = 4;
	config.integrator.output_interval = Real(0.5);
	const Real boundary = std::ldexp(Real(1), std::numeric_limits<Real>::digits - 1);
	const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(
		system, Vector<Real>{3}, boundary - 2, boundary + 2, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.accepted_steps == 1);
	REQUIRE(result.solution.has_value());
	REQUIRE(result.solution->size() == 6);
	for (int i = 1; i < result.solution->size(); ++i) {
		CHECK(result.solution->getTValue(i) > result.solution->getTValue(i - 1));
		CHECK(result.solution->getXValue(i, 0) == Real(3));
	}
	CHECK(result.solution->getTValue(5) == boundary + 2);
}

TEST_CASE("Next-step minimum failures preserve the completed endpoint and prior dense samples", "[ode-output]") {
	OutputODE system;
	ODEAdaptiveConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	config.integrator.initial_step_size = Real(0.25);
	config.integrator.min_step_size = Real(0.5);
	config.integrator.output_interval = Real(0.125);
	const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.accepted_steps == 1);
	REQUIRE(result.solution.has_value());
	REQUIRE(result.solution->size() == 3);
	CHECK(result.solution->getTValue(1) == Real(0.125));
	CHECK(result.solution->getTValue(2) == Real(0.25));
	CheckOutput(*result.solution, 1);

	ODEAdaptiveIntegrator<OutputStepper> integrator(system);
	const auto zero = integrator.integrate(Vector<Real>{0}, 0, 1, Real(0),
		Precision::ODEDefaultTolerance, Real(0.25));
	CHECK(zero.size() == 5);
	CHECK(zero.getNumStepsOK() == 4);
	CheckOutput(zero, 1);
}

TEST_CASE("Detailed exception policy also covers stepper construction", "[ode-output]") {
	OutputODE system;
	system.fault = 5;
	ODEAdaptiveConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto result = ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK_FALSE(result.solution.has_value());
	CHECK(result.accepted_steps == 0);
	CHECK_FALSE(result.error_message.empty());
	config.exception_policy = EvaluationExceptionPolicy::Propagate;
	CHECK_THROWS_AS(ODEAdaptiveIntegrateDetailed<OutputStepper>(system, Vector<Real>{0}, 0, 1, config), OutputCallbackError);
}
