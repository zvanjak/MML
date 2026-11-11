#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#include <mml/algorithms/ODESolvers/ODESolverEventDetection.h>
#endif

#include <limits>
#include <cmath>
#include <stdexcept>
#include <type_traits>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class WorkCallbackError : public std::runtime_error {
	public:
		WorkCallbackError() : std::runtime_error("RHS work contract callback failure") {}
	};

	class WorkODE : public IODESystem {
	public:
		mutable int calls = 0;
		int failAt = 0;
		int fault = 0;
		bool growth = false;
		Real rate = 1;
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>& state, Vector<Real>& rhs) const override {
			if (++calls == failAt) {
				if (fault == 0) throw WorkCallbackError();
				if (fault == 1) throw ResourceLimitError("RHS work contract resource failure");
				rhs[0] = std::numeric_limits<Real>::quiet_NaN();
			} else rhs[0] = growth ? rate * state[0] : Real(1);
		}
	};

	class InvalidWorkODE : public WorkODE {
	public:
		int declaredCalls = 0;
		int declaredAttempts = 1;
	};

	class InvalidWorkStepper {
		const InvalidWorkODE& _system;
	public:
		explicit InvalidWorkStepper(const IODESystem& system) : _system(dynamic_cast<const InvalidWorkODE&>(system)) {}
		void resetFSAL() {}
		StepResult doStep(Real, Vector<Real>&, Vector<Real>&, Real h, Real) {
			StepResult result;
			result.accepted = true;
			result.hDone = result.hNext = h;
			result.funcEvals = _system.declaredCalls;
			result.attemptedSteps = _system.declaredAttempts;
			return result;
		}
		Vector<Real> interpolate(Real time) const { return Vector<Real>{time}; }
	};

	class WorkEventODE : public IODESystemWithEvents {
	public:
		mutable int calls = 0;
		mutable int failAt = 0;
		mutable bool restarted = false;
		bool restart = false;
		bool failAfterRestart = false;
		bool growth = false;
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>& state, Vector<Real>& rhs) const override {
			if (++calls == failAt) throw WorkCallbackError();
			rhs[0] = growth ? 10 * state[0] : Real(1);
		}
		int getNumEvents() const override { return restart ? 1 : 0; }
		Real eventFunction(int, Real time, const Vector<Real>&) const override {
			return restarted ? Real(1) : time - Real(0.5);
		}
		EventAction getEventAction(int) const override { return EventAction::Restart; }
		void handleEvent(int, Real, Vector<Real>&) const override {
			restarted = true;
			if (failAfterRestart) failAt = calls + 2;
		}
	};

	template<class Test>
	void ForEachWorkStepper(Test test) {
		test.template operator()<DormandPrince5_Stepper>("DP5", 6);
		test.template operator()<CashKarp_Stepper>("Cash-Karp", 6);
		test.template operator()<DormandPrince8_Stepper>("DOP853", 15);
		test.template operator()<BulirschStoer_Stepper>("Bulirsch-Stoer", 7);
		test.template operator()<BulirschStoerRational_Stepper>("Bulirsch-Stoer alternate", 7);
	}

	template<class Operation>
	void CheckWorkFailure(Operation operation, int fault) {
		if (fault == 0) CHECK_THROWS_AS(operation(), WorkCallbackError);
		else if (fault == 1) CHECK_THROWS_AS(operation(), ResourceLimitError);
		else CHECK_THROWS_AS(operation(), ODESolverError);
	}
}

TEST_CASE("Adaptive direct steps count executed stages and dense setup but not supplied derivatives", "[ode-work]") {
	ForEachWorkStepper([]<class Stepper>(const char* name, int expectedCalls) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			WorkODE system;
			Stepper stepper(system);
			Vector<Real> state{0}, derivative{1};
			const Real h = direction / 8;
			for (int i = 0; i < 2; ++i) {
				const int before = system.calls;
				const auto result = stepper.doStep(i * h, state, derivative, h, Precision::ODEDefaultTolerance);
				REQUIRE(result.accepted);
				CHECK(result.hDone == h);
				CHECK(result.funcEvals == expectedCalls);
				CHECK(result.funcEvals == system.calls - before);
				CHECK(stepper.getLastFunctionEvaluations() == result.funcEvals);
				CHECK(result.attemptedSteps == 1);
				CHECK(result.rejectedSteps() == 0);
				CHECK(stepper.getLastAttemptedSteps() == 1);
				CHECK(stepper.getLastRejectedSteps() == 0);
				CHECK_THAT(state[0], WithinAbs((i + 1) * h, Precision::EqualityTolerance));
				const int beforeInterpolation = system.calls;
				CHECK_THAT(stepper.interpolate((i + Real(0.5)) * h)[0],
					WithinAbs((i + Real(0.5)) * h, Precision::EqualityTolerance));
				CHECK(system.calls == beforeInterpolation);
			}
			CHECK_THROWS_AS(stepper.doStep(2 * h, state, derivative, 0, Precision::ODEDefaultTolerance), ODESolverError);
			CHECK(stepper.getLastFunctionEvaluations() == 0);
			CHECK(stepper.getLastAttemptedSteps() == 0);
			CHECK(stepper.getLastRejectedSteps() == 0);
			CHECK(system.calls == 2 * expectedCalls);
		}
	});
}

TEST_CASE("Adaptive solve surfaces count startup estimation and every executed step exactly once", "[ode-work]") {
	ForEachWorkStepper([]<class Stepper>(const char* name, int stepCalls) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			for (bool automatic : {false, true}) {
				for (int surface : {0, 1, 2, 3}) {
					CAPTURE(direction, automatic, surface);
					WorkODE system;
					ODEAdaptiveIntegrator<Stepper> integrator(system);
					const Real h0 = automatic ? Real(0) : Real(0.125);
					ODEIntegratorConfig config;
					config.initial_step_size = h0;
					config.output_interval = Real(0.125);
					int evaluations = 0, accepted = 0, rejected = 0;
					if (surface == 3) {
						ODEAdaptiveConfig detailed;
						detailed.integrator = config;
						const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, direction, detailed);
						REQUIRE(result.IsSuccess());
						REQUIRE(result.solution.has_value());
						CHECK(result.solution->getTValue(result.solution->size() - 1) == direction);
						CHECK(result.function_evaluations == result.total_func_evals);
						evaluations = result.total_func_evals;
						accepted = result.accepted_steps;
						rejected = result.rejected_steps;
					} else {
						const auto solution = surface == 0
							? integrator.integrate(Vector<Real>{0}, 0, direction, Real(0.125), config.tolerance, h0)
							: surface == 1 ? integrator.integrate(Vector<Real>{0}, 0, direction, config)
							: integrator.integrateAt(Vector<Real>{0}, Vector<Real>{0, direction / 4, direction / 2, direction},
								config.tolerance, h0);
						CHECK(solution.getTValue(solution.size() - 1) == direction);
						evaluations = integrator.getStatistics().totalFuncEvals;
						accepted = integrator.getStatistics().acceptedSteps;
						rejected = integrator.getStatistics().rejectedSteps;
					}
					REQUIRE(accepted > 0);
					CHECK(rejected == 0);
					CHECK(evaluations == system.calls);
					CHECK(evaluations == 1 + (automatic ? 2 : 0) + accepted * stepCalls);
				}
			}
		}
	});
}

TEST_CASE("Adaptive Detailed counts every throwing or nonfinite RHS stage including dense setup", "[ode-work]") {
	ForEachWorkStepper([]<class Stepper>(const char* name, int stepCalls) {
		INFO(name);
		for (Real direction : {Real(-1), Real(1)}) {
			for (int call = 1; call <= 1 + stepCalls; ++call) {
				for (int fault : {0, 1, 2}) {
					CAPTURE(direction, call, fault);
					WorkODE system;
					system.failAt = call;
					system.fault = fault;
					ODEAdaptiveConfig config;
					config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
					config.integrator.initial_step_size = Real(0.125);
					const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, direction / 4, config);
					CHECK(result.status == (fault == 0 ? AlgorithmStatus::AlgorithmSpecificFailure :
						fault == 1 ? AlgorithmStatus::ResourceLimitExceeded : AlgorithmStatus::NumericalInstability));
					CHECK_FALSE(result.error_message.empty());
					CHECK(result.function_evaluations == call);
					CHECK(result.total_func_evals == system.calls);
					CHECK(system.calls == call);
					CHECK(result.accepted_steps == 0);
					CHECK(result.attempted_steps == (call == 1 ? 0 : 1));
					CHECK(result.rejected_steps == 0);
					REQUIRE(result.solution.has_value());
					REQUIRE(result.solution->size() == 1);
					CHECK(result.solution->getXValue(0, 0) == 0);
				}
			}
		}
	});
}

TEST_CASE("Adaptive startup failures include initial and automatic-estimator callback work", "[ode-work]") {
	ForEachWorkStepper([]<class Stepper>(const char* name, int) {
		INFO(name);
		for (int call : {1, 2, 3}) {
			for (int fault : {0, 1, 2}) {
				WorkODE system;
				system.failAt = call;
				system.fault = fault;
				ODEAdaptiveConfig config;
				config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
				config.check_finite = false;
				const auto result = ODEAdaptiveIntegrateDetailed<Stepper>(system, Vector<Real>{0}, 0, 1, config);
				CHECK_FALSE(result.IsSuccess());
				CHECK(result.function_evaluations == call);
				CHECK(result.total_func_evals == call);
				CHECK(system.calls == call);
				CHECK(result.accepted_steps == 0);
				REQUIRE(result.solution.has_value());
				CHECK(result.solution->size() == 1);
			}
		}
	});
}

TEST_CASE("Adaptive legacy failures retain actual work from accepted and incomplete attempts", "[ode-work]") {
	ForEachWorkStepper([]<class Stepper>(const char* name, int stepCalls) {
		INFO(name);
		for (int surface : {0, 1, 2}) {
			for (int fault : {0, 1, 2}) {
				WorkODE system;
				system.failAt = stepCalls + 2;
				system.fault = fault;
				ODEAdaptiveIntegrator<Stepper> integrator(system);
				ODEIntegratorConfig config;
				config.initial_step_size = Real(0.125);
				config.max_step_size = Real(0.125);
				auto integrate = [&]() {
					if (surface == 0)
						return integrator.integrate(Vector<Real>{0}, 0, 1, Real(0.125),
							Precision::ODEDefaultTolerance, Real(0.125));
					if (surface == 1) return integrator.integrate(Vector<Real>{0}, 0, 1, config);
					return integrator.integrateAt(Vector<Real>{0}, Vector<Real>{0, Real(0.125), Real(0.25)},
						Precision::ODEDefaultTolerance, Real(0.125));
				};
				if (surface == 1 && fault == 1) {
					const auto partial = integrate();
					REQUIRE(partial.size() == 2);
					CHECK(partial.getTValue(1) == Real(0.125));
					CHECK(integrator.getStatistics().status == AlgorithmStatus::ResourceLimitExceeded);
				} else CheckWorkFailure(integrate, fault);
				CHECK(system.calls == system.failAt);
				CHECK(integrator.getStatistics().totalFuncEvals == system.calls);
				CHECK(integrator.getStatistics().acceptedSteps == 1);
				CHECK(integrator.getStatistics().attemptedSteps == 2);
				CHECK(integrator.getStatistics().rejectedSteps == 0);
			}
		}
	});
}

TEMPLATE_TEST_CASE("Adaptive RK retry work excludes supplied derivatives on every attempt", "[ode-work]",
	DormandPrince5_Stepper, CashKarp_Stepper, DormandPrince8_Stepper) {
	for (Real direction : {Real(-1), Real(1)}) {
		WorkODE system;
		system.growth = true;
		TestType stepper(system);
		Vector<Real> state{1}, derivative{1};
		const Real h = direction * (std::is_same_v<Real, float> ? 16 : 4);
		const auto result = stepper.doStep(0, state, derivative, h, Precision::ODEDefaultTolerance);
		REQUIRE(result.accepted);
		CHECK(std::abs(result.hDone) < std::abs(h));
		CHECK(result.funcEvals == system.calls);
		CHECK(stepper.getLastFunctionEvaluations() == system.calls);
		constexpr int acceptedOnly = std::is_same_v<TestType, CashKarp_Stepper> ? 1 :
			std::is_same_v<TestType, DormandPrince8_Stepper> ? 3 : 0;
		constexpr int attemptCalls = std::is_same_v<TestType, CashKarp_Stepper> ? 5 :
			std::is_same_v<TestType, DormandPrince8_Stepper> ? 12 : 6;
		CHECK(result.funcEvals > attemptCalls + acceptedOnly);
		CHECK((result.funcEvals - acceptedOnly) % attemptCalls == 0);
		CHECK(result.attemptedSteps == (system.calls - acceptedOnly) / attemptCalls);
		CHECK(result.rejectedSteps() == result.attemptedSteps - 1);
		CHECK(stepper.getLastAttemptedSteps() == result.attemptedSteps);
		CHECK(stepper.getLastRejectedSteps() == result.rejectedSteps());
	}
}

TEST_CASE("DP5 internal rejected trials remain visible on every adaptive solve surface", "[ode-work][ode-attempts]") {
	for (int surface : {0, 1, 2, 3}) {
		CAPTURE(surface);
		WorkODE system;
		system.growth = true;
		system.rate = 10;
		DormandPrince5Integrator integrator(system);
		ODEIntegratorConfig config;
		config.initial_step_size = Real(0.1);
		config.tolerance = MML::Testing::Tol(Real(1e-12), Real(1e-5));
		int accepted = 0, rejected = 0, evaluations = 0, storedRejected = 0;
		if (surface == 3) {
			ODEAdaptiveConfig detailed;
			detailed.integrator = config;
			const auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{1}, 0, 1, detailed);
			REQUIRE(result.IsSuccess());
			REQUIRE(result.solution.has_value());
			accepted = result.accepted_steps;
			rejected = result.rejected_steps;
			evaluations = result.total_func_evals;
			storedRejected = result.solution->getNumStepsBad();
			CHECK(result.function_evaluations == system.calls);
			CHECK(result.attempted_steps == accepted + rejected);
		} else {
			const auto solution = surface == 0
				? integrator.integrate(Vector<Real>{1}, 0, 1, 0, config.tolerance, config.initial_step_size)
				: surface == 1 ? integrator.integrate(Vector<Real>{1}, 0, 1, config)
				: integrator.integrateAt(Vector<Real>{1}, Vector<Real>{0, 1}, config.tolerance, config.initial_step_size);
			accepted = integrator.getStatistics().acceptedSteps;
			rejected = integrator.getStatistics().rejectedSteps;
			evaluations = integrator.getStatistics().totalFuncEvals;
			storedRejected = solution.getNumStepsBad();
			CHECK(integrator.getStatistics().attemptedSteps == accepted + rejected);
			CHECK(integrator.getStatistics().totalSteps() == accepted + rejected);
			CHECK(solution.getTotalNumSteps() == accepted + rejected);
		}
		CHECK(evaluations == system.calls);
		REQUIRE(accepted > 0);
		REQUIRE(system.calls > 1 + 6 * accepted);
		CHECK(rejected == (system.calls - 1) / 6 - accepted);
		CHECK(rejected > 0);
		CHECK(storedRejected == rejected);
	}
}

TEMPLATE_TEST_CASE("Adaptive RK max_steps bounds internal attempts rather than accepted doStep returns",
	"[ode-work][ode-attempts]", DormandPrince5_Stepper, CashKarp_Stepper, DormandPrince8_Stepper) {
	for (Real direction : {Real(-1), Real(1)}) {
		WorkODE system;
		system.growth = true;
		ODEAdaptiveConfig config;
		config.integrator.initial_step_size = std::is_same_v<Real, float> ? 16 : 4;
		config.integrator.max_steps = 1;
		const auto result = ODEAdaptiveIntegrateDetailed<TestType>(system, Vector<Real>{1}, 0,
			2 * direction * config.integrator.initial_step_size, config);
		CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(result.accepted_steps == 0);
		CHECK(result.rejected_steps == 1);
		CHECK(result.attempted_steps == 1);
		constexpr int callsPerAttempt = std::is_same_v<TestType, CashKarp_Stepper> ? 5 :
			std::is_same_v<TestType, DormandPrince8_Stepper> ? 12 : 6;
		CHECK(system.calls == 1 + callsPerAttempt);
		CHECK(result.total_func_evals == system.calls);
		REQUIRE(result.solution.has_value());
		CHECK(result.solution->size() == 1);
		CHECK(result.solution->getNumStepsBad() == 1);
	}
}

TEST_CASE("Interrupted RK retries retain completed rejections and the incomplete attempt separately", "[ode-work][ode-attempts]") {
	for (int fault : {0, 1, 2}) {
		WorkODE system;
		system.growth = true;
		system.failAt = 8;
		system.fault = fault;
		ODEAdaptiveConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.integrator.initial_step_size = std::is_same_v<Real, float> ? 16 : 4;
		const auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{1}, 0,
			2 * config.integrator.initial_step_size, config);
		CHECK(result.status == (fault == 0 ? AlgorithmStatus::AlgorithmSpecificFailure :
			fault == 1 ? AlgorithmStatus::ResourceLimitExceeded : AlgorithmStatus::NumericalInstability));
		CHECK(result.attempted_steps == 2);
		CHECK(result.accepted_steps == 0);
		CHECK(result.rejected_steps == 1);
		CHECK(result.function_evaluations == 8);
		CHECK(result.total_func_evals == system.calls);
		REQUIRE(result.solution.has_value());
		CHECK(result.solution->size() == 1);
		CHECK(result.solution->getNumStepsBad() == 1);
	}
	WorkODE system;
	system.growth = true;
	system.failAt = 7;
	DormandPrince5_Stepper stepper(system);
	Vector<Real> state{1}, derivative{1};
	const Real h = std::is_same_v<Real, float> ? 16 : 4;
	CHECK_THROWS_AS(stepper.doStep(0, state, derivative, h, Precision::ODEDefaultTolerance), WorkCallbackError);
	CHECK(stepper.getLastAttemptedSteps() == 2);
	CHECK(stepper.getLastRejectedSteps() == 1);
	CHECK(stepper.getLastFunctionEvaluations() == 7);
	CHECK(state[0] == 1);
	CHECK_THROWS_AS(stepper.doStep(0, state, derivative, 0, Precision::ODEDefaultTolerance), ODESolverError);
	CHECK(stepper.getLastAttemptedSteps() == 0);
	CHECK(stepper.getLastRejectedSteps() == 0);
	CHECK(system.calls == 7);
}

TEST_CASE("Event integration counts startup, internal retries and restart estimation without resetting history", "[ode-work][ode-attempts]") {
	for (bool automatic : {false, true}) {
		for (bool restart : {false, true}) {
			WorkEventODE system;
			system.restart = restart;
			system.growth = !restart;
			DormandPrince5EventIntegrator integrator(system);
			const auto result = integrator.integrateWithEvents(system, Vector<Real>{restart ? Real(0) : Real(1)},
				0, 1, Real(0.125), MML::Testing::Tol(Real(1e-12), Real(1e-5)),
				Precision::EventTolerance, automatic ? Real(0) : Real(0.1));
			const auto& stats = integrator.getStatistics();
			CHECK(result.status == AlgorithmStatus::Success);
			CHECK(result.finalTime == 1);
			CHECK(stats.totalFuncEvals == system.calls);
			CHECK(stats.attemptedSteps == stats.acceptedSteps + stats.rejectedSteps);
			CHECK(result.solution.getNumStepsOK() == stats.acceptedSteps);
			CHECK(result.solution.getNumStepsBad() == stats.rejectedSteps);
			REQUIRE(stats.acceptedSteps > 0);
			if (restart) {
				CHECK(system.restarted);
				CHECK(result.events.size() == 1);
				CHECK(stats.acceptedSteps > 1);
			} else if (!automatic) CHECK(stats.rejectedSteps > 0);
		}
	}
}

TEST_CASE("Event integration retains actual callback work before stage and restart-probe failures", "[ode-work][ode-attempts]") {
	for (int fault : {0, 1, 2}) {
		WorkEventODE system;
		system.restart = fault == 2;
		system.failAfterRestart = fault == 2;
		system.failAt = fault == 0 ? 1 : fault == 1 ? 3 : 0;
		DormandPrince5EventIntegrator integrator(system);
		CHECK_THROWS_AS(integrator.integrateWithEvents(system, Vector<Real>{0}, 0, 1, Real(0.125),
			Precision::ODEDefaultTolerance, Precision::EventTolerance, Real(0.125)), WorkCallbackError);
		const auto& stats = integrator.getStatistics();
		CHECK(stats.totalFuncEvals == system.calls);
		CHECK(stats.rejectedSteps == 0);
		if (fault == 0) CHECK(stats.attemptedSteps == 0);
		else if (fault == 1) CHECK(stats.attemptedSteps == 1);
		else {
			CHECK(system.restarted);
			CHECK(stats.acceptedSteps > 1);
			CHECK(stats.attemptedSteps == stats.acceptedSteps);
		}
	}
}

TEST_CASE("Custom step attempt metadata is validated and defaults preserve single-trial compatibility", "[ode-work][ode-attempts]") {
	for (int attempts : {-1, 0, std::numeric_limits<int>::max(), 2}) {
		InvalidWorkODE system;
		system.declaredAttempts = attempts;
		ODEAdaptiveConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.integrator.initial_step_size = Real(0.125);
		config.integrator.max_steps = 2;
		const auto result = ODEAdaptiveIntegrateDetailed<InvalidWorkStepper>(system, Vector<Real>{0}, 0, 1, config);
		CHECK(result.total_func_evals == 1);
		CHECK(system.calls == 1);
		REQUIRE(result.solution.has_value());
		if (attempts == 2) {
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK(result.attempted_steps == 2);
			CHECK(result.accepted_steps == 1);
			CHECK(result.rejected_steps == 1);
			CHECK(result.solution->getNumStepsBad() == 1);
		} else {
			CHECK(result.status == (attempts < 1 ? AlgorithmStatus::NumericalInstability : AlgorithmStatus::ResourceLimitExceeded));
			CHECK(result.attempted_steps == 0);
			CHECK(result.accepted_steps == 0);
			CHECK(result.rejected_steps == 0);
		}
	}
	StepResult legacy;
	CHECK(legacy.attemptedSteps == 1);
	CHECK(legacy.rejectedSteps() == 1);
	legacy.accepted = true;
	CHECK(legacy.rejectedSteps() == 0);
}

TEST_CASE("Adaptive work resets per operation and invalid input executes no callbacks", "[ode-work]") {
	ForEachWorkStepper([]<class Stepper>(const char* name, int) {
		INFO(name);
		WorkODE system;
		ODEAdaptiveIntegrator<Stepper> integrator(system);
		integrator.integrate(Vector<Real>{0}, 0, 1, Real(0.125), Precision::ODEDefaultTolerance, Real(0.125));
		CHECK(integrator.getStatistics().totalFuncEvals == system.calls);
		const int beforeEstimator = system.calls;
		CHECK(integrator.estimateInitialStep(0, Vector<Real>{0}, 1, Precision::ODEDefaultTolerance) > 0);
		CHECK(system.calls - beforeEstimator == 2);
		CHECK(integrator.getStatistics().totalFuncEvals == 2);
		const int beforeInvalid = system.calls;
		CHECK_THROWS_AS(integrator.integrateAt(Vector<Real>{0}, Vector<Real>{0, 0}), ODESolverError);
		CHECK(integrator.getStatistics().totalFuncEvals == 0);
		CHECK(system.calls == beforeInvalid);
		integrator.integrateAt(Vector<Real>{0}, Vector<Real>{0, Real(0.5), 1},
			Precision::ODEDefaultTolerance, Real(0.125));
		CHECK(integrator.getStatistics().totalFuncEvals == system.calls - beforeInvalid);
	});
}

TEST_CASE("Adaptive counters reject exhausted budgets and malformed custom work without overflow", "[ode-work]") {
	WorkODE system;
	int evaluations = 2;
	ODEStepperDetail::CountingODESystem counted(system, evaluations, 3);
	Vector<Real> state{0}, derivative(1);
	counted.derivs(0, state, derivative);
	CHECK(evaluations == 3);
	CHECK(system.calls == 1);
	CHECK_THROWS_AS(counted.derivs(0, state, derivative), ResourceLimitError);
	CHECK(evaluations == 3);
	CHECK(system.calls == 1);
	evaluations = std::numeric_limits<int>::max();
	ODEStepperDetail::CountingODESystem exhausted(system, evaluations);
	CHECK_THROWS_AS(exhausted.derivs(0, state, derivative), ResourceLimitError);
	CHECK(system.calls == 1);
	for (int declaredCalls : {-1, std::numeric_limits<int>::max()}) {
		InvalidWorkODE invalid;
		invalid.declaredCalls = declaredCalls;
		ODEAdaptiveConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.integrator.initial_step_size = Real(0.125);
		const auto result = ODEAdaptiveIntegrateDetailed<InvalidWorkStepper>(invalid, state, 0, 1, config);
		CHECK(result.status == (declaredCalls < 0 ? AlgorithmStatus::NumericalInstability : AlgorithmStatus::ResourceLimitExceeded));
		CHECK(result.function_evaluations == 1);
		CHECK(result.total_func_evals == invalid.calls);
		CHECK(result.accepted_steps == 0);
		REQUIRE(result.solution.has_value());
		CHECK(result.solution->size() == 1);
	}
}
