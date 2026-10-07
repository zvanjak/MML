///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootIsolation.h                                                     ///
///  Description: Structured scalar root isolation over finite intervals             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOT_ISOLATION_H
#define MML_ROOT_ISOLATION_H

#include <mml/MMLBase.h>
#include <mml/base/Function.h>
#include <mml/interfaces/IFunction.h>
#include <mml/algorithms/RootFinding/RootFindingBase.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace MML::RootFinding {

	enum class RootCandidateType {
		SignChange,
		ExactSample,
		TangentCandidate
	};

	struct RootInterval {
		Real lower = 0.0;
		Real upper = 0.0;
		Real f_lower = 0.0;
		Real f_upper = 0.0;

		[[nodiscard]] Real Width() const noexcept { return upper - lower; }
		[[nodiscard]] bool IsPoint() const noexcept { return lower == upper; }
	};

	struct RootCandidate {
		RootCandidateType type = RootCandidateType::SignChange;
		RootInterval interval;
		Real estimate = 0.0;
		Real residual = 0.0;
		int sample_index = -1;
	};

	struct RootIsolationConfig {
		int num_intervals = 100;
		Real zero_tolerance = 0.0;
		Real merge_tolerance = 0.0;
		Real tangent_tolerance = std::sqrt(Constants::Eps);
		int tangent_refinement_iterations = 48;
		bool detect_tangent_roots = true;
		bool check_finite = true;
	};

	inline bool IsValidConfig(const RootIsolationConfig& config) {
		return config.num_intervals > 0
			&& std::isfinite(config.zero_tolerance) && config.zero_tolerance >= 0.0
			&& std::isfinite(config.merge_tolerance) && config.merge_tolerance >= 0.0
			&& std::isfinite(config.tangent_tolerance) && config.tangent_tolerance >= 0.0
			&& config.tangent_refinement_iterations > 0;
	}

	namespace Detail {
		inline RootCandidate RefineTangentCandidate(const IRealFunction& function,
			Real lower, Real upper, int sampleIndex, const RootIsolationConfig& config) {
			constexpr Real goldenSection = Real(0.6180339887498948482);
			Real left = lower;
			Real right = upper;
			Real xLeft = right - goldenSection * (right - left);
			Real xRight = left + goldenSection * (right - left);
			Real fLeft = function(xLeft);
			Real fRight = function(xRight);
			if (config.check_finite && (!std::isfinite(fLeft) || !std::isfinite(fRight)))
				throw RootFindingError("FindRootIntervals: non-finite tangent refinement value");

			for (int iteration = 0; iteration < config.tangent_refinement_iterations; ++iteration) {
				if (std::abs(fLeft) < std::abs(fRight)) {
					right = xRight;
					xRight = xLeft;
					fRight = fLeft;
					xLeft = right - goldenSection * (right - left);
					fLeft = function(xLeft);
					if (config.check_finite && !std::isfinite(fLeft))
						throw RootFindingError("FindRootIntervals: non-finite tangent refinement value");
				} else {
					left = xLeft;
					xLeft = xRight;
					fLeft = fRight;
					xRight = left + goldenSection * (right - left);
					fRight = function(xRight);
					if (config.check_finite && !std::isfinite(fRight))
						throw RootFindingError("FindRootIntervals: non-finite tangent refinement value");
				}
			}

			const bool useLeft = std::abs(fLeft) <= std::abs(fRight);
			const Real estimate = useLeft ? xLeft : xRight;
			const Real value = useLeft ? fLeft : fRight;
			const Real fLower = function(lower);
			const Real fUpper = function(upper);
			if (config.check_finite && (!std::isfinite(fLower) || !std::isfinite(fUpper)))
				throw RootFindingError("FindRootIntervals: non-finite tangent interval endpoint");
			return {RootCandidateType::TangentCandidate,
				{lower, upper, fLower, fUpper},
				estimate, std::abs(value), sampleIndex};
		}

		inline std::vector<RootCandidate> FindRootIntervalsImpl(const IRealFunction& function,
			const IRealFunction* derivative, Real lower, Real upper,
			const RootIsolationConfig& config) {
		if (!IsValidConfig(config) || !std::isfinite(lower) || !std::isfinite(upper) || lower >= upper)
			throw RootFindingError("FindRootIntervals: invalid configuration or interval");

		const Real step = (upper - lower) / config.num_intervals;
		const Real mergeTolerance = config.merge_tolerance > 0.0
			? config.merge_tolerance
			: Real(8.0) * Constants::Eps * std::max({Real(1.0), std::abs(lower), std::abs(upper)});
		std::vector<Real> points(config.num_intervals + 1);
		std::vector<Real> values(config.num_intervals + 1);
		std::vector<Real> derivativeValues;
		if (derivative != nullptr) derivativeValues.resize(config.num_intervals + 1);
		for (int sampleIndex = 0; sampleIndex <= config.num_intervals; ++sampleIndex) {
			points[sampleIndex] = sampleIndex == config.num_intervals
				? upper : lower + sampleIndex * step;
			values[sampleIndex] = function(points[sampleIndex]);
			if (config.check_finite && !std::isfinite(values[sampleIndex]))
				throw RootFindingError("FindRootIntervals: non-finite function value");
			if (derivative != nullptr) {
				derivativeValues[sampleIndex] = (*derivative)(points[sampleIndex]);
				if (config.check_finite && !std::isfinite(derivativeValues[sampleIndex]))
					throw RootFindingError("FindRootIntervals: non-finite derivative value");
			}
		}

		std::vector<RootCandidate> rawCandidates;
		for (int sampleIndex = 0; sampleIndex <= config.num_intervals; ++sampleIndex) {
			const Real xCurrent = points[sampleIndex];
			const Real fCurrent = values[sampleIndex];

			const bool currentIsRoot = std::abs(fCurrent) <= config.zero_tolerance;
			if (currentIsRoot) {
				rawCandidates.push_back({RootCandidateType::ExactSample,
					{xCurrent, xCurrent, fCurrent, fCurrent}, xCurrent, std::abs(fCurrent), sampleIndex});
			} else if (sampleIndex > 0
				&& std::abs(values[sampleIndex - 1]) > config.zero_tolerance
				&& std::signbit(values[sampleIndex - 1]) != std::signbit(fCurrent)) {
				rawCandidates.push_back({RootCandidateType::SignChange,
					{points[sampleIndex - 1], xCurrent, values[sampleIndex - 1], fCurrent},
					Real(0.5) * (points[sampleIndex - 1] + xCurrent),
					std::min(std::abs(values[sampleIndex - 1]), std::abs(fCurrent)), sampleIndex - 1});
			}
		}

		if (config.detect_tangent_roots) {
			for (int sampleIndex = 1; sampleIndex < config.num_intervals; ++sampleIndex) {
				const Real absPrevious = std::abs(values[sampleIndex - 1]);
				const Real absCurrent = std::abs(values[sampleIndex]);
				const Real absNext = std::abs(values[sampleIndex + 1]);
				if (absCurrent > absPrevious || absCurrent > absNext
					|| absCurrent <= config.zero_tolerance) continue;
				const bool leftSignChange = std::signbit(values[sampleIndex - 1])
					!= std::signbit(values[sampleIndex]);
				const bool rightSignChange = std::signbit(values[sampleIndex])
					!= std::signbit(values[sampleIndex + 1]);
				if (leftSignChange || rightSignChange) continue;

				if (derivative != nullptr) {
					const bool derivativeTurns = std::signbit(derivativeValues[sampleIndex - 1])
						!= std::signbit(derivativeValues[sampleIndex + 1]);
					const bool derivativeIsSmall = std::abs(derivativeValues[sampleIndex])
						<= std::sqrt(Constants::Eps);
					if (!derivativeTurns && !derivativeIsSmall) continue;
				}

				RootCandidate candidate = RefineTangentCandidate(function,
					points[sampleIndex - 1], points[sampleIndex + 1], sampleIndex, config);
				if (candidate.residual <= config.tangent_tolerance)
					rawCandidates.push_back(candidate);
			}
		}

		std::sort(rawCandidates.begin(), rawCandidates.end(),
			[](const RootCandidate& left, const RootCandidate& right) {
				return left.estimate < right.estimate;
			});
		std::vector<RootCandidate> candidates;
		for (const RootCandidate& candidate : rawCandidates) {
			if (candidates.empty()
				|| std::abs(candidate.estimate - candidates.back().estimate) > mergeTolerance) {
				candidates.push_back(candidate);
				continue;
			}
			RootCandidate& previous = candidates.back();
			if (candidate.residual < previous.residual
				|| (candidate.type == RootCandidateType::ExactSample
					&& previous.type != RootCandidateType::ExactSample)) previous = candidate;
		}
		return candidates;
		}
	} // namespace Detail

	inline std::vector<RootCandidate> FindRootIntervals(const IRealFunction& function,
		Real lower, Real upper, const RootIsolationConfig& config = {}) {
		return Detail::FindRootIntervalsImpl(function, nullptr, lower, upper, config);
	}

	inline std::vector<RootCandidate> FindRootIntervals(const IRealFunction& function,
		const IRealFunction& derivative, Real lower, Real upper,
		const RootIsolationConfig& config = {}) {
		return Detail::FindRootIntervalsImpl(function, &derivative, lower, upper, config);
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	inline std::vector<RootCandidate> FindRootIntervals(Function&& function,
		Real lower, Real upper, const RootIsolationConfig& config = {}) {
		MML::Detail::RealFunctionCallableAdapter<Function> adapter(function);
		return FindRootIntervals(adapter, lower, upper, config);
	}

	template<RealFunctionCallable Function, RealFunctionCallable Derivative>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>
			&& !std::derived_from<std::remove_cvref_t<Derivative>, IRealFunction>)
	inline std::vector<RootCandidate> FindRootIntervals(Function&& function,
		Derivative&& derivative, Real lower, Real upper,
		const RootIsolationConfig& config = {}) {
		MML::Detail::RealFunctionCallableAdapter<Function> functionAdapter(function);
		MML::Detail::RealFunctionCallableAdapter<Derivative> derivativeAdapter(derivative);
		return FindRootIntervals(functionAdapter, derivativeAdapter, lower, upper, config);
	}

} // namespace MML::RootFinding

#endif // MML_ROOT_ISOLATION_H
