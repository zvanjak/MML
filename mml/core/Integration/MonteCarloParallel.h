#ifndef MML_MONTE_CARLO_PARALLEL_H
#define MML_MONTE_CARLO_PARALLEL_H

#include <mml/base/ParallelExecution.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <cstdint>
#include <random>

namespace MML {
	struct ParallelMonteCarloConfig : MonteCarloConfig, ParallelExecutionConfig {
		size_t block_size = ParallelDefaults::MonteCarloBlockSize; ///< Independent samples per indexed RNG block; pairs for antithetic sampling.
		size_t min_parallel_samples = ParallelDefaults::MonteCarloMinSamples; ///< Smaller evaluation budgets execute the same blocks inline.
		EvaluationExceptionPolicy exception_policy = EvaluationExceptionPolicy::Propagate;
	};

	struct ParallelMonteCarloResult : MonteCarloResult {
		size_t blocks_used = 0;
		size_t threads_used = 0;
		unsigned int seed_used = 0;
	};

	template<int N, class Function>
	ParallelMonteCarloResult IntegrateMonteCarloParallelDetailed(Function&& function,
		const VectorN<Real, N>& lower, const VectorN<Real, N>& upper,
		const ParallelMonteCarloConfig& config = {}) {
		EvaluationConfigBase policy;
		policy.exception_policy = config.exception_policy;
		return IntegrationDetail::ExecuteIntegrationDetailed<ParallelMonteCarloResult>(
			"MonteCarloParallel", policy, [&](ParallelMonteCarloResult& result) {
				Detail::ValidateMonteCarloConfig(config);
				if (config.block_size == 0) throw ArgumentError("Parallel Monte Carlo block_size must be positive");
				const Real volume = Detail::ValidateMonteCarloBounds(lower, upper, "MonteCarloParallel");
				if (volume == 0) {
					result.status = AlgorithmStatus::Completed;
					result.error_estimate = 0;
					return;
				}
				const size_t samples = config.num_samples / (config.use_antithetic ? 2 : 1);
				const size_t blocks = 1 + (samples - 1) / config.block_size;
				std::vector<Detail::WelfordAccumulator> statistics(blocks);
				result.seed_used = config.seed ? config.seed : std::random_device{}();
				ParallelExecutionConfig execution = config;
				if (config.num_samples < config.min_parallel_samples) execution.thread_count = 1;
				result.threads_used = ParallelDetail::Execute(blocks, execution, [&](size_t block) {
					const std::uint64_t index = block;
					std::seed_seq seed{result.seed_used, static_cast<unsigned int>(index),
						static_cast<unsigned int>(index >> 32)};
					std::mt19937 engine(seed);
					std::uniform_real_distribution<Real> uniform(0, 1);
					const size_t start = block * config.block_size;
					const size_t length = std::min(config.block_size, samples - start);
					VectorN<Real, N> point, reflected;
					for (size_t i = 0; i < length; ++i) {
						for (int d = 0; d < N; ++d) {
							const Real u = uniform(engine);
							point[d] = std::lerp(lower[d], upper[d], u);
							reflected[d] = std::lerp(lower[d], upper[d], Real(1) - u);
						}
						const Real value = std::invoke(function, point);
						if (!std::isfinite(value)) throw NumericalMethodError("Parallel Monte Carlo: non-finite integrand");
						if (config.use_antithetic) {
							const Real opposite = std::invoke(function, reflected);
							if (!std::isfinite(opposite)) throw NumericalMethodError("Parallel Monte Carlo: non-finite antithetic value");
							statistics[block].add(value / 2 + opposite / 2);
						} else statistics[block].add(value);
					}
				});
				Detail::WelfordAccumulator combined;
				for (const auto& block : statistics) combined.merge(block);
				static_cast<MonteCarloResult&>(result) = Detail::FinishMonteCarlo(combined, volume, config.num_samples, config);
				result.algorithm_name = "MonteCarloParallel";
				result.blocks_used = blocks;
			});
	}

	template<int N, class Function>
	Real IntegrateMonteCarloParallel(Function&& function,
		const VectorN<Real, N>& lower, const VectorN<Real, N>& upper,
		const ParallelMonteCarloConfig& config = {}) {
		return IntegrateMonteCarloParallelDetailed(std::forward<Function>(function), lower, upper, config)
			.ValueOrThrow("IntegrateMonteCarloParallel");
	}
}
#endif
