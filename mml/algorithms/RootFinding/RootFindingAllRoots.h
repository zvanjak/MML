///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingAllRoots.h                                               ///
///  Description: Isolate and refine all detectable real roots on an interval         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOT_FINDING_ALL_ROOTS_H
#define MML_ROOT_FINDING_ALL_ROOTS_H

#include <mml/algorithms/RootFinding/RootFindingBase.h>
#include <mml/algorithms/RootFinding/RootFindingMethods.h>
#include <mml/algorithms/RootFinding/RootIsolation.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace MML::RootFinding {

	struct FindAllRealRootsConfig {
		RootIsolationConfig isolation;
		RootFindingConfig refinement;
		Real merge_tolerance = 0.0;
		bool keep_rejected_candidates = true;
	};

	struct RealRootResult {
		RootCandidate candidate;
		RootFindingResult refinement;
	};

	struct FindAllRealRootsResult : public EvaluationResultBase {
		std::vector<Real> roots;
		std::vector<RealRootResult> accepted;
		std::vector<RealRootResult> rejected;
		int candidates_found = 0;
		int candidates_accepted = 0;
		int roots_accepted = 0;
		int candidates_rejected = 0;
		bool complete = false;
	};

	namespace Detail {
		class CountingRealFunction : public IRealFunction {
		public:
			explicit CountingRealFunction(const IRealFunction& function) : _function(function) {}

			Real operator()(Real x) const override {
				++_evaluations;
				return _function(x);
			}

			[[nodiscard]] int Evaluations() const noexcept { return _evaluations; }

		private:
			const IRealFunction& _function;
			mutable int _evaluations = 0;
		};

		inline RootFindingResult RefineExactCandidate(const IRealFunction& function,
			const RootCandidate& candidate, const RootFindingConfig& config) {
			RootFindingResult result;
			result.root = candidate.estimate;
			result.function_value = function(candidate.estimate);
			result.function_evaluations = 1;
			result.algorithm_name = "ExactSample";
			result.achieved_tolerance = std::abs(result.function_value);
			if (!std::isfinite(result.function_value)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Non-finite exact-sample residual";
			} else if (std::abs(result.function_value) <= EffectiveFTolerance(config)) {
				result.converged = true;
				result.status = AlgorithmStatus::Success;
			} else {
				result.status = AlgorithmStatus::ToleranceUnachievable;
				result.error_message = "Exact sample did not satisfy final residual tolerance";
			}
			return result;
		}

		inline RootFindingResult RefineTangentRoot(const IRealFunction& function,
			const RootCandidate& candidate, const FindAllRealRootsConfig& config) {
			RootFindingResult result;
			RootFindingResultFinalizer finalizer(result, "TangentResidualMinimization");
			const int iterations = std::max(1, config.isolation.tangent_refinement_iterations);
			RootIsolationConfig localConfig = config.isolation;
			localConfig.tangent_refinement_iterations = iterations;
			RootCandidate refined = RefineTangentCandidate(function,
				candidate.interval.lower, candidate.interval.upper,
				candidate.sample_index, localConfig);
			result.root = refined.estimate;
			result.function_value = function(refined.estimate);
			result.function_evaluations = iterations + 5;
			result.iterations_used = iterations;
			result.x_error = refined.interval.Width();
			result.achieved_tolerance = std::abs(result.function_value);
			if (!std::isfinite(result.function_value)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Non-finite tangent candidate residual";
				return result;
			}
			if (std::abs(result.function_value) <= EffectiveFTolerance(config.refinement)) {
				result.converged = true;
				return result;
			}
			result.status = AlgorithmStatus::ToleranceUnachievable;
			result.error_message = "Tangent candidate did not satisfy final residual tolerance";
			return result;
		}
	} // namespace Detail

	inline bool IsValidConfig(const FindAllRealRootsConfig& config) {
		return IsValidConfig(config.isolation) && IsValidConfig(config.refinement)
			&& std::isfinite(config.merge_tolerance) && config.merge_tolerance >= 0.0;
	}

	inline FindAllRealRootsResult FindAllRealRootsInInterval(const IRealFunction& function,
		Real lower, Real upper, const FindAllRealRootsConfig& config = {}) {
		AlgorithmTimer timer;
		FindAllRealRootsResult result;
		result.algorithm_name = "FindAllRealRootsInInterval";
		if (!IsValidConfig(config) || !std::isfinite(lower) || !std::isfinite(upper) || lower >= upper) {
			result.status = AlgorithmStatus::InvalidInput;
			result.error_message = "FindAllRealRootsInInterval: invalid configuration or interval";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		Detail::CountingRealFunction counted(function);
		std::vector<RootCandidate> candidates;
		try {
			candidates = FindRootIntervals(counted, lower, upper, config.isolation);
		} catch (const std::exception& error) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.error_message = error.what();
			result.function_evaluations = counted.Evaluations();
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		result.candidates_found = static_cast<int>(candidates.size());

		int failedRefinements = 0;
		for (const RootCandidate& candidate : candidates) {
			RootFindingResult refinement;
			try {
				if (candidate.type == RootCandidateType::SignChange) {
					refinement = FindRootBrent(counted, candidate.interval.lower,
						candidate.interval.upper, config.refinement);
				} else if (candidate.type == RootCandidateType::ExactSample) {
					refinement = Detail::RefineExactCandidate(counted, candidate, config.refinement);
				} else {
					refinement = Detail::RefineTangentRoot(counted, candidate, config);
				}
			} catch (const std::exception& error) {
				refinement.root = candidate.estimate;
				refinement.function_value = std::numeric_limits<Real>::quiet_NaN();
				refinement.status = AlgorithmStatus::NumericalInstability;
				refinement.algorithm_name = "CandidateRefinement";
				refinement.error_message = error.what();
			}

			RealRootResult rootResult{candidate, refinement};
			if (refinement.IsSuccess()) result.accepted.push_back(rootResult);
			else {
				++failedRefinements;
				if (config.keep_rejected_candidates) result.rejected.push_back(rootResult);
			}
		}

		const Real mergeTolerance = config.merge_tolerance > 0.0
			? config.merge_tolerance
			: EffectiveXTolerance(config.refinement, std::max(std::abs(lower), std::abs(upper)));
		std::sort(result.accepted.begin(), result.accepted.end(),
			[](const RealRootResult& left, const RealRootResult& right) {
				return left.refinement.root < right.refinement.root;
			});
		for (const RealRootResult& accepted : result.accepted) {
			if (result.roots.empty() || std::abs(accepted.refinement.root - result.roots.back()) > mergeTolerance) {
				result.roots.push_back(accepted.refinement.root);
			} else if (std::abs(accepted.refinement.function_value)
				< std::abs(counted(result.roots.back()))) {
				result.roots.back() = accepted.refinement.root;
			}
		}

		result.roots_accepted = static_cast<int>(result.roots.size());
		result.candidates_accepted = static_cast<int>(result.accepted.size());
		result.candidates_rejected = failedRefinements;
		result.complete = result.candidates_rejected == 0;
		result.status = result.complete ? AlgorithmStatus::Success : AlgorithmStatus::ToleranceUnachievable;
		if (!result.complete)
			result.error_message = "One or more isolated candidates failed final refinement";
		result.function_evaluations = counted.Evaluations();
		result.elapsed_time_ms = timer.elapsed_ms();
		return result;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	inline FindAllRealRootsResult FindAllRealRootsInInterval(Function&& function,
		Real lower, Real upper, const FindAllRealRootsConfig& config = {}) {
		MML::Detail::RealFunctionCallableAdapter<Function> adapter(function);
		return FindAllRealRootsInInterval(adapter, lower, upper, config);
	}

} // namespace MML::RootFinding

#endif // MML_ROOT_FINDING_ALL_ROOTS_H
