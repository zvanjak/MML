#include "../../TestPrecision.h"
#include <test_beds/stiff_solver_benchmarks.h>
#include <catch2/catch_test_macros.hpp>

using namespace MML;
using namespace MML::TestBeds;

TEST_CASE("Stiff benchmark references - fixed-step refinement on classical problems",
	"[StiffSolvers][stiff-benchmark]")
{
	for (const auto& problem : StiffBenchmarkProblems())
		for (auto method : {StiffBenchmarkMethod::BackwardEuler, StiffBenchmarkMethod::BDF2}) {
			CAPTURE(problem.name, StiffBenchmarkMethodName(method));
			StiffSolverConfig config;
			config.max_steps = 50000;
			config.step_size = problem.smokeEndTime / Real(100);
			config.newton_tol = Precision::ZeroTolerance;
			auto coarse = RunStiffBenchmark(problem, method, config, true);
			config.step_size /= 4;
			auto fine = RunStiffBenchmark(problem, method, config, true);
			INFO(coarse.message);
			INFO(fine.message);
			REQUIRE(coarse.status == AlgorithmStatus::Success);
			REQUIRE(fine.status == AlgorithmStatus::Success);
			REQUIRE(fine.acceptedSteps >= coarse.acceptedSteps * 3);
			REQUIRE(fine.rhsCalls > coarse.rhsCalls);
			REQUIRE(fine.jacobianCalls > 0);
			REQUIRE(fine.maxAbsoluteError < coarse.maxAbsoluteError);
			REQUIRE(fine.maxAbsoluteError < Precision::CoarseTolerance * Real(10));
			if (method == StiffBenchmarkMethod::BDF2)
				REQUIRE(fine.maxAbsoluteError < Precision::CoarseTolerance);
			REQUIRE(fine.finalState.size() == problem.initial.size());
		}
}

TEST_CASE("Stiff benchmark references - ROS2 accuracy and measured work",
	"[StiffSolvers][stiff-benchmark][ROS2]")
{
	for (const auto& problem : StiffBenchmarkProblems()) {
		CAPTURE(problem.name);
		StiffSolverConfig config;
		config.max_steps = 50000;
		config.step_size = problem.smokeEndTime / Real(100);
		config.tolerance = Precision::CoarseTolerance;
		auto coarse = RunStiffBenchmark(problem, StiffBenchmarkMethod::ROS2, config, true);
		config.tolerance /= Real(16);
		config.step_size /= 4;
		auto fine = RunStiffBenchmark(problem, StiffBenchmarkMethod::ROS2, config, true);
		INFO(coarse.message);
		INFO(fine.message);
		REQUIRE(coarse.status == AlgorithmStatus::Success);
		REQUIRE(fine.status == AlgorithmStatus::Success);
		REQUIRE(fine.maxAbsoluteError < coarse.maxAbsoluteError);
		REQUIRE(fine.maxAbsoluteError < Precision::CoarseTolerance);
		REQUIRE(fine.rhsCalls > coarse.rhsCalls);
		REQUIRE(fine.jacobianCalls > 0);
		REQUIRE(fine.acceptedSteps > 0);
		REQUIRE(fine.rhsCalls == 3 * fine.jacobianCalls);
		REQUIRE(fine.jacobianCalls == static_cast<std::uint64_t>(fine.acceptedSteps + fine.rejectedSteps));
		REQUIRE(fine.luFactorizations == fine.jacobianCalls);
		REQUIRE(fine.linearSolves == 2 * fine.luFactorizations);
	}
}

TEST_CASE("Stiff benchmark - exhausted budgets are explicit, not success-shaped endpoints",
	"[StiffSolvers][stiff-benchmark]")
{
	auto problem = StiffBenchmarkProblems().front();
	StiffSolverConfig config;
	config.step_size = problem.smokeEndTime / Real(100);
	config.max_steps = 1;
	auto measurement = RunStiffBenchmark(problem, StiffBenchmarkMethod::BackwardEuler, config, true);
	REQUIRE(measurement.status == AlgorithmStatus::MaxIterationsExceeded);
	REQUIRE_FALSE(measurement.message.empty());
	REQUIRE(std::isnan(measurement.maxAbsoluteError));
	REQUIRE(std::isnan(measurement.maxScaledError));
	REQUIRE(measurement.rhsCalls > 0);
	REQUIRE(measurement.acceptedSteps == 1);
	REQUIRE(measurement.reachedTime == config.step_size);
}
