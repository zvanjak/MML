#include <MML.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Single header uses the base pool for deterministic parallel integrals", "[single_header][parallel]")
{
	using namespace MML;
	ParallelIntegrationConfig config;
	config.thread_count = 2;
	config.partitions = 4;
	config.method = GAUSS10KRONROD21;
	const auto result = IntegrateParallelDetailed([](Real x) { return x * x; }, 0, 1, config);
	REQUIRE(result.status == AlgorithmStatus::Completed);
	REQUIRE_FALSE(result.converged);
	REQUIRE(result.IsUsable());
	REQUIRE(result.function_evaluations == 84);
	REQUIRE(IsNearlyEqual(result.value, Real(1) / 3, Precision::EqualityTolerance));
	ParallelMonteCarloConfig mc;
	mc.num_samples = 10000;
	mc.seed = 42;
	mc.thread_count = 2;
	const auto estimate = IntegrateMonteCarloParallelDetailed([](const VectorN<Real, 1>&) { return Real(2); },
		VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, mc);
	REQUIRE(estimate.value == Real(2));
	REQUIRE(estimate.status == AlgorithmStatus::Success);
	REQUIRE(estimate.samples_used == 10000);
}
