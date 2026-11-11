#include <catch2/catch_all.hpp>
#include <mml/core/Integration/Integration1D.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <mml/core/Integration/IntegrationParallel.h>
#include <mml/core/Integration/MonteCarloParallel.h>
#include "../../TestPrecision.h"
#include <atomic>
#include <barrier>
#include <limits>
#include <set>
#include <mutex>
#include <sstream>

namespace MML::Tests::ParallelIntegration {
	TEST_CASE("Fixed GK estimates do not claim untested convergence", "[integration][parallel][GK-status]")
	{
		auto f = [](Real x) { return x * x; };
		const auto fixed = Integration::IntegrateGK21Core(f, 0, 1);
		REQUIRE(fixed.status == AlgorithmStatus::Completed);
		REQUIRE_FALSE(fixed.converged);
		REQUIRE(fixed.error_estimated);
		REQUIRE(fixed.IsUsable());
		REQUIRE_NOTHROW(fixed.ValueOrThrow("GK21"));
		std::ostringstream printed;
		printed << fixed;
		REQUIRE(printed.str().find("not converged") == std::string::npos);
		const auto adaptive = Integration::IntegrateGKAdaptiveCore(f, 0, 1);
		REQUIRE(adaptive.converged);
		REQUIRE(adaptive.status == AlgorithmStatus::Success);
		for (const auto& result : {Integration::IntegrateGK15Core(f, 0, 1),
			Integration::IntegrateGK31Core(f, 0, 1)}) {
			REQUIRE(result.status == AlgorithmStatus::Completed);
			REQUIRE_FALSE(result.converged);
			REQUIRE(result.IsUsable());
		}
		const auto limited = Integration::IntegrateGKAdaptiveCore(
			[](Real x) { return std::abs(x - Real(0.3)); }, 0, 1,
			std::numeric_limits<Real>::min(), std::numeric_limits<Real>::min(), 0);
		REQUIRE(limited.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE_FALSE(limited.converged);
		REQUIRE_FALSE(limited.IsUsable());
	}

	TEST_CASE("Monte Carlo zero mean cannot bypass statistical uncertainty", "[montecarlo][parallel][MC-status]")
	{
		class Alternating : public IScalarFunction<1> {
			mutable int calls = 0;
		public:
			Real operator()(const VectorN<Real, 1>&) const override { return ++calls % 2 ? Real(-1) : Real(1); }
		} function;
		MonteCarloIntegrator<1> integrator(42);
		MonteCarloConfig config;
		config.num_samples = 2;
		const auto result = integrator.integrate(function, config);
		REQUIRE(result.value == Real(0));
		REQUIRE(result.error_estimate > 0);
		REQUIRE_FALSE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE_FALSE(result.IsUsable());
		config.absolute_error = 2;
		const auto absolute = integrator.integrate(function, config);
		REQUIRE(absolute.converged);
		REQUIRE(absolute.status == AlgorithmStatus::Success);
	}

	TEST_CASE("Parallel quadrature has worker-independent ordered results", "[integration][parallel]")
	{
		ParallelIntegrationConfig config;
		config.partitions = 12;
		config.method = GAUSS10KRONROD21;
		config.thread_count = 1;
		auto f = [](Real x) { return std::exp(x); };
		const auto sequential = IntegrateParallelDetailed(f, 0, 1, config);
		config.thread_count = 4;
		const auto parallel = IntegrateParallelDetailed(f, 0, 1, config);
		REQUIRE(parallel.value == sequential.value);
		REQUIRE(parallel.error_estimate == sequential.error_estimate);
		REQUIRE(parallel.status == AlgorithmStatus::Completed);
		REQUIRE_FALSE(parallel.converged);
		REQUIRE(parallel.IsUsable());
		REQUIRE(parallel.function_evaluations == 12 * 21);
		REQUIRE(parallel.partitions_used == 12);
		REQUIRE(parallel.threads_used == 4);
		REQUIRE(IsNearlyEqual(parallel.value, std::exp(Real(1)) - 1, TOL(1e-12, 1e-5)));
		const auto reversed = IntegrateParallelDetailed(f, 1, 0, config);
		REQUIRE(IsNearlyEqual(reversed.value, -parallel.value, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Parallel integration preserves fixed and checked numerical semantics", "[integration][parallel]")
	{
		ParallelIntegrationConfig config;
		config.partitions = 8;
		config.thread_count = 3;
		config.method = GAUSS10;
		const auto fixed = IntegrateParallelDetailed([](Real x) { return x * x; }, 0, 1, config);
		REQUIRE(fixed.status == AlgorithmStatus::Completed);
		REQUIRE_FALSE(fixed.error_estimated);
		REQUIRE(std::isnan(fixed.error_estimate));
		REQUIRE(fixed.function_evaluations == 80);
		config.method = SIMPSON;
		const auto smooth = IntegrateParallelDetailed([](Real x) { return x * x; }, 0, 1, config);
		REQUIRE(smooth.converged);
		REQUIRE(smooth.status == AlgorithmStatus::Success);
		REQUIRE(IsNearlyEqual(smooth.value, Real(1) / 3, TOL(1e-12, 1e-5)));
		const auto cancellation = IntegrateParallelDetailed([](Real x) { return x; }, -1, 1, config);
		REQUIRE(IsWithinAbs(cancellation.value, 0, TOL(1e-12, 1e-5)));
		const auto oscillation = IntegrateParallelDetailed([](Real x) { return std::sin(40 * x); }, 0, Constants::PI, config);
		REQUIRE(oscillation.status == AlgorithmStatus::Success);
		REQUIRE(IsWithinAbs(oscillation.value, 0, TOL(1e-10, 1e-4)));
		const auto saved = AlgorithmContext::Get();
		AlgorithmContext::Get().simpsonIntegrationMaxSteps = 1;
		const auto failed = IntegrateParallelDetailed([](Real x) { return std::exp(x); }, 0, 1, config);
		AlgorithmContext::Get() = saved;
		REQUIRE(failed.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE_FALSE(failed.converged);
		REQUIRE_FALSE(failed.IsUsable());
		REQUIRE(failed.error_message.find("partition 0") != std::string::npos);
	}

	TEST_CASE("Parallel quadrature accepts interfaces and avoids pools for trivial work", "[integration][parallel]")
	{
		class Square : public IRealFunction {
		public: Real operator()(Real x) const override { return x * x; }
		} function;
		const IRealFunction& interface = function;
		ParallelIntegrationConfig config;
		config.thread_count = 0;
		config.partitions = 1;
		const auto one = IntegrateParallelDetailed(interface, 0, 1, config);
		REQUIRE(one.threads_used == 1);
		config.partitions = 2;
		config.min_parallel_tasks = 3;
		const auto small = IntegrateParallelDetailed(interface, 0, 1, config);
		REQUIRE(small.threads_used == 1);
		const auto zero = IntegrateParallelDetailed([](Real) -> Real { throw std::runtime_error("must not run"); }, 2, 2, config);
		REQUIRE(zero.value == Real(0));
		REQUIRE(zero.status == AlgorithmStatus::Completed);
		REQUIRE(zero.function_evaluations == 0);
		REQUIRE(zero.threads_used == 0);
	}

	TEST_CASE("Parallel batches really run concurrently and restore numerical context", "[parallel][threadpool]")
	{
		ParallelExecutionConfig config;
		config.thread_count = 2;
		std::barrier gate(2);
		std::mutex mutex;
		std::set<std::thread::id> threads;
		const auto original = AlgorithmContext::Get();
		const int expectedBudget = original.simpsonIntegrationMaxSteps;
		std::atomic<int> observed{0};
		ParallelDetail::Execute(2, config, [&](size_t) {
			if (AlgorithmContext::Get().simpsonIntegrationMaxSteps == expectedBudget) ++observed;
			gate.arrive_and_wait();
			std::lock_guard<std::mutex> lock(mutex);
			threads.insert(std::this_thread::get_id());
			AlgorithmContext::Get().simpsonIntegrationMaxSteps = 1;
		});
		REQUIRE(threads.size() == 2);
		REQUIRE(observed == 2);
		REQUIRE(AlgorithmContext::Get().simpsonIntegrationMaxSteps == expectedBudget);
		config.thread_count = 1;
		ParallelDetail::Execute(3, config, [&](size_t) { AlgorithmContext::Get().simpsonIntegrationMaxSteps = 1; });
		REQUIRE(AlgorithmContext::Get().simpsonIntegrationMaxSteps == expectedBudget);
	}

	TEST_CASE("Parallel task failures drain the batch and rethrow by index", "[parallel][threadpool]")
	{
		ParallelExecutionConfig config;
		config.thread_count = 3;
		std::atomic<int> finished{0};
		REQUIRE_THROWS_WITH(ParallelDetail::Execute(9, config, [&](size_t index) {
			++finished;
			if (index == 0) throw std::runtime_error("first partition");
			if (index == 1) throw std::runtime_error("second partition");
		}), "first partition");
		REQUIRE(finished == 9);
		ParallelIntegrationConfig integration;
		integration.thread_count = 2;
		integration.partitions = 4;
		integration.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto result = IntegrateParallelDetailed([](Real) -> Real { throw std::runtime_error("callback failed"); }, 0, 1, integration);
		REQUIRE(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
		REQUIRE(result.error_message == "callback failed");
	}

	TEST_CASE("Parallel integration rejects invalid and unrepresentable work", "[integration][parallel]")
	{
		ParallelIntegrationConfig config;
		config.partitions = std::numeric_limits<size_t>::max();
		REQUIRE_THROWS_AS(IntegrateParallelDetailed([](Real x) { return x; }, 0, 1, config), ArgumentError);
		config.partitions = 4;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.tolerance = 0;
		REQUIRE(IntegrateParallelDetailed([](Real x) { return x; }, 0, 1, config).status == AlgorithmStatus::InvalidInput);
		config.tolerance = Precision::IntegrationTolerance;
		const auto nonfinite = IntegrateParallelDetailed([](Real) { return std::numeric_limits<Real>::infinity(); }, 0, 1, config);
		REQUIRE(nonfinite.status == AlgorithmStatus::NumericalInstability);
		config.tolerance = Precision::IntegrationTolerance;
		const auto collapsed = IntegrateParallelDetailed([](Real x) { return x; },
			Real(1), std::nextafter(Real(1), Real(2)), config);
		REQUIRE(collapsed.status == AlgorithmStatus::InvalidInput);
	}

	TEST_CASE("Parallel Monte Carlo is reproducible across workers and leaves the caller RNG untouched", "[montecarlo][parallel]")
	{
		ParallelMonteCarloConfig config;
		config.num_samples = 20000;
		config.block_size = 257;
		config.seed = 123;
		config.thread_count = 1;
		const VectorN<Real, 2> lower{0, 0}, upper{1, 1};
		auto function = [](const VectorN<Real, 2>& x) { return x[0] * x[1]; };
		auto savedEngine = Random::Engine();
		const auto reference = IntegrateMonteCarloParallelDetailed(function, lower, upper, config);
		config.thread_count = 4;
		const auto parallel = IntegrateMonteCarloParallelDetailed(function, lower, upper, config);
		const auto repeated = IntegrateMonteCarloParallelDetailed(function, lower, upper, config);
		REQUIRE(parallel.value == reference.value);
		REQUIRE(parallel.variance == reference.variance);
		REQUIRE(parallel.error_estimate == reference.error_estimate);
		REQUIRE(repeated.value == parallel.value);
		REQUIRE(parallel.samples_used == 20000);
		REQUIRE(parallel.function_evaluations == 20000);
		REQUIRE(parallel.seed_used == 123);
		REQUIRE(parallel.threads_used == 4);
		REQUIRE(parallel.blocks_used == reference.blocks_used);
		REQUIRE(Random::Engine()() == savedEngine());
		REQUIRE(std::abs(parallel.value - Real(1) / 4) <= 6 * parallel.error_estimate);
	}

	TEST_CASE("Parallel Welford merge preserves between-block variation", "[montecarlo][parallel]")
	{
		Detail::WelfordAccumulator a, b, all;
		for (Real value : {Real(1000000), Real(1000001), Real(1000002), Real(1000003)}) {
			all.add(value);
			if (value < Real(1000002)) a.add(value); else b.add(value);
		}
		a.merge(b);
		REQUIRE(a.count == all.count);
		REQUIRE(a.mean == all.mean);
		REQUIRE(IsNearlyEqual(a.sampleVariance(), all.sampleVariance(), Precision::EqualityTolerance));
	}

	TEST_CASE("Parallel antithetic statistics count pairs but retain the exact evaluation budget", "[montecarlo][parallel]")
	{
		ParallelMonteCarloConfig config;
		config.use_antithetic = true;
		config.num_samples = 10000;
		config.seed = 42;
		config.block_size = 100;
		config.thread_count = 3;
		const auto result = IntegrateMonteCarloParallelDetailed([](const VectorN<Real, 1>& x) { return x[0]; },
			VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, config);
		REQUIRE(result.value == Real(0.5));
		REQUIRE(result.error_estimate == 0);
		REQUIRE(result.variance == 0);
		REQUIRE(result.samples_used == 10000);
		REQUIRE(result.function_evaluations == 10000);
		REQUIRE(result.status == AlgorithmStatus::Success);
		config.num_samples = 5;
		REQUIRE_THROWS_AS(IntegrateMonteCarloParallelDetailed([](const VectorN<Real, 1>& x) { return x[0]; },
			VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, config), ArgumentError);
	}

	TEST_CASE("Parallel Monte Carlo handles zero domains small work and expected failures", "[montecarlo][parallel]")
	{
		ParallelMonteCarloConfig config;
		config.seed = 7;
		config.num_samples = 10;
		config.block_size = 2;
		config.thread_count = 3;
		const VectorN<Real, 1> lower{0}, upper{1};
		auto constant = [](const VectorN<Real, 1>&) { return Real(2); };
		const auto small = IntegrateMonteCarloParallelDetailed(constant, lower, upper, config);
		REQUIRE(small.threads_used == 1);
		REQUIRE(small.blocks_used == 5);
		const auto zero = IntegrateMonteCarloParallelDetailed(constant, lower, lower, config);
		REQUIRE(zero.value == 0);
		REQUIRE(zero.status == AlgorithmStatus::Completed);
		REQUIRE(zero.samples_used == 0);
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto failed = IntegrateMonteCarloParallelDetailed([](const VectorN<Real, 1>&) -> Real { throw std::runtime_error("sample failure"); }, lower, upper, config);
		REQUIRE(failed.status == AlgorithmStatus::AlgorithmSpecificFailure);
		REQUIRE(failed.error_message == "sample failure");
		config.block_size = 0;
		REQUIRE(IntegrateMonteCarloParallelDetailed(constant, lower, upper, config).status == AlgorithmStatus::InvalidInput);
		config.block_size = 2;
		const auto nonfinite = IntegrateMonteCarloParallelDetailed(
			[](const VectorN<Real, 1>&) { return std::numeric_limits<Real>::quiet_NaN(); }, lower, upper, config);
		REQUIRE(nonfinite.status == AlgorithmStatus::NumericalInstability);
		REQUIRE_THROWS_AS(IntegrateMonteCarloParallelDetailed(constant, upper, lower), ArgumentError);
		config.target_error = std::numeric_limits<Real>::min();
		const auto missed = IntegrateMonteCarloParallelDetailed([](const VectorN<Real, 1>& x) { return x[0]; },
			lower, upper, config);
		REQUIRE(missed.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE_FALSE(missed.converged);
		REQUIRE(missed.error_estimate > 0);
		REQUIRE_FALSE(missed.error_message.empty());
	}

	TEST_CASE("Serial Monte Carlo rejects false uncertainty and nonfinite samples", "[montecarlo][MC-status]")
	{
		class Constant : public IScalarFunction<1> {
		public: Real operator()(const VectorN<Real, 1>&) const override { return Real(1); }
		} function;
		MonteCarloIntegrator<1> integrator(42);
		MonteCarloConfig config;
		config.num_samples = 1;
		REQUIRE_THROWS_AS(integrator.integrate(function, config), ArgumentError);
		config.use_antithetic = true;
		for (size_t n : {size_t(2), size_t(3), size_t(5)}) {
			config.num_samples = n;
			REQUIRE_THROWS_AS(integrator.integrate(function, config), ArgumentError);
		}
	}
}
