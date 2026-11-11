#include <mml/core/Integration.h>
#include <iostream>
#include <stdexcept>

void Docs_Demo_ParallelIntegration()
{
	using namespace MML;
	ParallelIntegrationConfig quadrature;
	quadrature.thread_count = 4;
	quadrature.partitions = 32;
	quadrature.method = GAUSS10KRONROD21;
	const auto integral = IntegrateParallelDetailed([](Real x) { return x * x; }, 0, 1, quadrature);
	if (integral.status != AlgorithmStatus::Completed ||
		!IsNearlyEqual(integral.value, Real(1) / 3, Precision::EqualityTolerance) ||
		integral.function_evaluations != 32 * 21)
		throw std::runtime_error("Composite parallel quadrature demo failed");
	std::cout << "Parallel quadrature: value=" << integral.value << ", partitions="
	          << integral.partitions_used << ", workers=" << integral.threads_used
	          << ", evaluations=" << integral.function_evaluations
	          << ", status=" << to_string(integral.status) << '\n';

	ParallelMonteCarloConfig sampling;
	sampling.seed = 123;
	sampling.num_samples = 20000;
	sampling.thread_count = 1;
	sampling.use_antithetic = true;
	auto function = [](const VectorN<Real, 2>& x) { return x[0] + x[1]; };
	const VectorN<Real, 2> lower{0, 0}, upper{1, 1};
	const auto serialBlocks = IntegrateMonteCarloParallelDetailed(function, lower, upper, sampling);
	sampling.thread_count = 4;
	const auto parallel = IntegrateMonteCarloParallelDetailed(function, lower, upper, sampling);
	if (parallel.status != AlgorithmStatus::Success || !IsNearlyEqual(parallel.value, Real(1), Precision::EqualityTolerance) ||
		parallel.value != serialBlocks.value || parallel.variance != serialBlocks.variance ||
		parallel.samples_used != sampling.num_samples)
		throw std::runtime_error("Deterministic parallel Monte Carlo demo failed");
	std::cout << "Parallel Monte Carlo: value=" << parallel.value << ", standard error="
	          << parallel.error_estimate << ", evaluations=" << parallel.function_evaluations
	          << ", blocks=" << parallel.blocks_used << ", workers=" << parallel.threads_used
	          << ", seed=" << parallel.seed_used << ", status=" << to_string(parallel.status) << '\n';
}
