#ifndef MML_INTEGRATION_PARALLEL_H
#define MML_INTEGRATION_PARALLEL_H

#include <mml/base/ParallelExecution.h>
#include <mml/core/Integration/Integration1D.h>
#include <limits>
#include <numeric>

namespace MML {
	struct ParallelIntegrationConfig : IntegrationConfig, ParallelExecutionConfig {
		IntegrationMethod method = SIMPSON;
		Real tolerance = Precision::IntegrationTolerance; ///< Partition-local tolerance, not a global relative bound.
		size_t partitions = 0; ///< Zero selects four indexed partitions per requested worker.
	};

	struct ParallelIntegrationResult : IntegrationResult {
		size_t partitions_used = 0;
		size_t threads_used = 0;
	};

	template<RealFunctionCallable Function>
	ParallelIntegrationResult IntegrateParallelDetailed(Function&& function, Real a, Real b,
		const ParallelIntegrationConfig& config = {}) {
		return IntegrationDetail::ExecuteIntegrationDetailed<ParallelIntegrationResult>(
			"IntegrateParallel", config, [&](ParallelIntegrationResult& result) {
				if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(config.tolerance) || config.tolerance <= 0)
					throw ArgumentError("Parallel integration requires finite endpoints and positive finite tolerance");
				if (config.method < TRAP || config.method > GAUSS10KRONROD21)
					throw ArgumentError("Parallel integration: invalid quadrature method");
				const auto& context = AlgorithmContext::Get();
				if ((config.method == TRAP && context.trapezoidIntegrationMaxSteps <= 0) ||
					(config.method == SIMPSON && context.simpsonIntegrationMaxSteps <= 0) ||
					(config.method == ROMBERG && context.rombergIntegrationMaxSteps <= 0))
					throw ArgumentError("Parallel integration requires a positive refinement budget");
				if (a == b) {
					static_cast<IntegrationResult&>(result) = IntegrationResult::FixedRule(0);
					result.algorithm_name = "IntegrateParallel";
					return;
				}
				const size_t maxCount = static_cast<size_t>(std::numeric_limits<int>::max());
				const size_t workers = ParallelDetail::WorkerCount(config.thread_count, maxCount);
				const size_t factor = ParallelDefaults::PartitionsPerWorker;
				const size_t partitions = config.partitions ? config.partitions : std::min(workers, maxCount / factor) * factor;
				if (partitions > maxCount)
					throw ArgumentError("Parallel integration partition count exceeds supported counter range");
				std::vector<IntegrationResult> results(partitions);
				result.threads_used = ParallelDetail::Execute(partitions, config, [&](size_t index) {
					const Real left = std::lerp(a, b, Real(index) / Real(partitions));
					const Real right = index + 1 == partitions ? b : std::lerp(a, b, Real(index + 1) / Real(partitions));
					if (left == right) throw ArgumentError("Parallel integration partitions cannot advance floating-point coordinates");
					int evaluations = 0;
					auto counted = [&](Real x) {
						if (evaluations == std::numeric_limits<int>::max())
							throw ResourceLimitError("Parallel integration function counter exhausted");
						++evaluations;
						const Real value = std::invoke(function, x);
						if (!std::isfinite(value)) throw NumericalMethodError("Parallel integration: non-finite integrand");
						return value;
					};
					results[index] = IntegrateCore(counted, left, right, config.method, config.tolerance);
					results[index].function_evaluations = evaluations;
				});
				Real correction = 0, error = 0;
				bool allConverged = true, allEstimated = true;
				result.status = AlgorithmStatus::Success;
				std::int64_t evaluations = 0, iterations = 0;
				for (size_t index = 0; index < partitions; ++index) {
					const auto& local = results[index];
					// Neumaier compensated reduction in fixed partition order.
					const Real next = result.value + local.value;
					correction += std::abs(result.value) >= std::abs(local.value)
						? (result.value - next) + local.value : (local.value - next) + result.value;
					result.value = next;
					allConverged = allConverged && local.converged;
					allEstimated = allEstimated && local.error_estimated;
					if (local.error_estimated) error += std::abs(local.error_estimate);
					evaluations += local.function_evaluations;
					iterations += local.iterations;
					if (!local.IsUsable() && MML::IsUsable(result.status)) {
						result.status = local.status;
						result.error_message = "Parallel integration partition " + std::to_string(index) +
							": " + (local.error_message.empty() ? "quadrature did not converge" : local.error_message);
					}
				}
				result.value += correction;
				if (!std::isfinite(result.value) || (allEstimated && !std::isfinite(error)))
					throw NumericalMethodError("Parallel integration reduction overflow");
				if (evaluations > std::numeric_limits<int>::max() || iterations > std::numeric_limits<int>::max())
					throw ResourceLimitError("Parallel integration aggregate counter exceeds supported range");
				result.function_evaluations = static_cast<int>(evaluations);
				result.iterations = static_cast<int>(iterations);
				result.partitions_used = partitions;
				result.error_estimated = allEstimated;
				result.error_estimate = allEstimated ? error : NoEstimate;
				result.converged = allConverged && MML::IsUsable(result.status);
				if (MML::IsUsable(result.status))
					result.status = result.converged ? AlgorithmStatus::Success : AlgorithmStatus::Completed;
			});
	}

	template<RealFunctionCallable Function>
	Real IntegrateParallel(Function&& function, Real a, Real b,
		const ParallelIntegrationConfig& config = {}) {
		return IntegrateParallelDetailed(std::forward<Function>(function), a, b, config).ValueOrThrow("IntegrateParallel");
	}
}
#endif
