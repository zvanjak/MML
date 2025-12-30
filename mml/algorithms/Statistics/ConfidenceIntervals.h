///////////////////////////////////////////////////////////////////////////////////////////
// ConfidenceIntervals.h
// 
// Confidence interval functions for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Requires: Statistics.h (for basic stats functions and distributions)
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_CONFIDENCE_INTERVALS_H
#define MML_CONFIDENCE_INTERVALS_H

#include "../Statistics.h"
#include "CoreDistributions.h"

namespace MML
{
	namespace Statistics
	{
		/**
		 * @brief Structure to hold confidence interval results
		 */
		struct ConfidenceInterval
		{
			Real estimate;           // Point estimate
			Real lowerBound;         // Lower confidence limit
			Real upperBound;         // Upper confidence limit
			Real marginOfError;      // Half-width of interval
			Real confidenceLevel;    // Confidence level (e.g., 0.95 for 95%)
			std::string parameter;   // What we're estimating (e.g., "Mean", "Proportion")

			ConfidenceInterval(Real est, Real lower, Real upper, Real margin, Real conf, const std::string& param)
				: estimate(est), lowerBound(lower), upperBound(upper), 
				  marginOfError(margin), confidenceLevel(conf), parameter(param) {}
		};

		/**
		 * @brief Confidence Interval for Population Mean (single sample)
		 * 
		 * Uses t-distribution for unknown population variance.
		 * CI = x̄ ± t(α/2, df) × (s/√n)
		 * 
		 * @param sample Sample data
		 * @param confidenceLevel Confidence level (default: 0.95 for 95% CI)
		 * @return ConfidenceInterval with bounds
		 */
		inline ConfidenceInterval ConfidenceIntervalMean(
			const Vector<Real>& sample,
			Real confidenceLevel = 0.95
		)
		{
			if (sample.size() < 2)
				throw StatisticsError("ConfidenceIntervalMean requires at least 2 samples");

			if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0)
				throw StatisticsError("Confidence level must be in (0, 1)");

			Real sampleMean = Mean(sample);
			Real sampleStd = StdDev(sample);
			int n = static_cast<int>(sample.size());
			int df = n - 1;

			// Standard error
			Real standardError = sampleStd / std::sqrt(static_cast<Real>(n));

			// Critical value from t-distribution
			Real alpha = 1.0 - confidenceLevel;
			TDistribution tDist(df);
			Real tCritical = tDist.criticalValue(alpha);

			// Margin of error
			Real marginOfError = tCritical * standardError;

			return ConfidenceInterval(
				sampleMean,
				sampleMean - marginOfError,
				sampleMean + marginOfError,
				marginOfError,
				confidenceLevel,
				"Mean"
			);
		}

		/**
		 * @brief Confidence Interval for Difference of Two Means (independent samples)
		 * 
		 * Uses pooled variance assuming equal population variances.
		 * CI = (x̄₁ - x̄₂) ± t(α/2, df) × SE
		 * 
		 * @param sample1 First sample
		 * @param sample2 Second sample
		 * @param confidenceLevel Confidence level (default: 0.95)
		 * @return ConfidenceInterval for μ₁ - μ₂
		 */
		inline ConfidenceInterval ConfidenceIntervalMeanDifference(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real confidenceLevel = 0.95
		)
		{
			if (sample1.size() < 2 || sample2.size() < 2)
				throw StatisticsError("ConfidenceIntervalMeanDifference requires at least 2 samples per group");

			if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0)
				throw StatisticsError("Confidence level must be in (0, 1)");

			Real mean1 = Mean(sample1);
			Real mean2 = Mean(sample2);
			Real var1 = Variance(sample1);
			Real var2 = Variance(sample2);
			int n1 = static_cast<int>(sample1.size());
			int n2 = static_cast<int>(sample2.size());

			// Pooled variance
			Real pooledVar = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2);
			Real pooledStd = std::sqrt(pooledVar);

			// Standard error
			Real standardError = pooledStd * std::sqrt(1.0 / n1 + 1.0 / n2);

			// Degrees of freedom
			int df = n1 + n2 - 2;

			// Critical value
			Real alpha = 1.0 - confidenceLevel;
			TDistribution tDist(df);
			Real tCritical = tDist.criticalValue(alpha);

			// Margin of error
			Real marginOfError = tCritical * standardError;

			// Difference
			Real difference = mean1 - mean2;

			return ConfidenceInterval(
				difference,
				difference - marginOfError,
				difference + marginOfError,
				marginOfError,
				confidenceLevel,
				"Mean Difference"
			);
		}

		/**
		 * @brief Confidence Interval for Population Proportion
		 * 
		 * Uses normal approximation (valid for large samples).
		 * CI = p̂ ± z(α/2) × √(p̂(1-p̂)/n)
		 * 
		 * @param successes Number of successes
		 * @param trials Total number of trials
		 * @param confidenceLevel Confidence level (default: 0.95)
		 * @return ConfidenceInterval for proportion
		 */
		inline ConfidenceInterval ConfidenceIntervalProportion(
			int successes,
			int trials,
			Real confidenceLevel = 0.95
		)
		{
			if (trials <= 0)
				throw StatisticsError("Number of trials must be positive");

			if (successes < 0 || successes > trials)
				throw StatisticsError("Successes must be between 0 and trials");

			if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0)
				throw StatisticsError("Confidence level must be in (0, 1)");

			// Sample proportion
			Real pHat = static_cast<Real>(successes) / trials;

			// Standard error
			Real standardError = std::sqrt(pHat * (1.0 - pHat) / trials);

			// Critical value from standard normal
			Real alpha = 1.0 - confidenceLevel;
			NormalDistribution norm(0.0, 1.0);
			Real zCritical = norm.inverseCdf(1.0 - alpha / 2.0);

			// Margin of error
			Real marginOfError = zCritical * standardError;

			// Bounds (clamped to [0, 1])
			Real lower = std::max(0.0, pHat - marginOfError);
			Real upper = std::min(1.0, pHat + marginOfError);

			return ConfidenceInterval(
				pHat,
				lower,
				upper,
				marginOfError,
				confidenceLevel,
				"Proportion"
			);
		}

		/**
		 * @brief Confidence Interval for Difference of Two Proportions
		 * 
		 * Uses normal approximation.
		 * CI = (p̂₁ - p̂₂) ± z(α/2) × √(p̂₁(1-p̂₁)/n₁ + p̂₂(1-p̂₂)/n₂)
		 * 
		 * @param successes1 Successes in first sample
		 * @param trials1 Trials in first sample
		 * @param successes2 Successes in second sample
		 * @param trials2 Trials in second sample
		 * @param confidenceLevel Confidence level (default: 0.95)
		 * @return ConfidenceInterval for p₁ - p₂
		 */
		inline ConfidenceInterval ConfidenceIntervalProportionDifference(
			int successes1, int trials1,
			int successes2, int trials2,
			Real confidenceLevel = 0.95
		)
		{
			if (trials1 <= 0 || trials2 <= 0)
				throw StatisticsError("Number of trials must be positive");

			if (successes1 < 0 || successes1 > trials1 || successes2 < 0 || successes2 > trials2)
				throw StatisticsError("Successes must be between 0 and trials");

			if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0)
				throw StatisticsError("Confidence level must be in (0, 1)");

			// Sample proportions
			Real p1 = static_cast<Real>(successes1) / trials1;
			Real p2 = static_cast<Real>(successes2) / trials2;

			// Standard error
			Real standardError = std::sqrt(
				p1 * (1.0 - p1) / trials1 + 
				p2 * (1.0 - p2) / trials2
			);

			// Critical value
			Real alpha = 1.0 - confidenceLevel;
			NormalDistribution norm(0.0, 1.0);
			Real zCritical = norm.inverseCdf(1.0 - alpha / 2.0);

			// Margin of error
			Real marginOfError = zCritical * standardError;

			// Difference
			Real difference = p1 - p2;

			return ConfidenceInterval(
				difference,
				difference - marginOfError,
				difference + marginOfError,
				marginOfError,
				confidenceLevel,
				"Proportion Difference"
			);
		}

		/**
		 * @brief Confidence Interval for Mean of Paired Differences
		 * 
		 * For paired data (before-after, matched pairs).
		 * CI = d̄ ± t(α/2, df) × (s_d/√n)
		 * 
		 * @param before First measurements
		 * @param after Second measurements
		 * @param confidenceLevel Confidence level (default: 0.95)
		 * @return ConfidenceInterval for mean difference
		 */
		inline ConfidenceInterval ConfidenceIntervalPairedDifference(
			const Vector<Real>& before,
			const Vector<Real>& after,
			Real confidenceLevel = 0.95
		)
		{
			if (before.size() != after.size())
				throw StatisticsError("ConfidenceIntervalPairedDifference requires equal sample sizes");

			if (before.size() < 2)
				throw StatisticsError("ConfidenceIntervalPairedDifference requires at least 2 pairs");

			if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0)
				throw StatisticsError("Confidence level must be in (0, 1)");

			// Compute differences
			Vector<Real> differences(before.size());
			for (size_t i = 0; i < before.size(); i++)
				differences[i] = after[i] - before[i];

			// Use standard mean CI on differences
			return ConfidenceIntervalMean(differences, confidenceLevel);
		}

	} // namespace Statistics
}  // namespace MML

#endif // MML_CONFIDENCE_INTERVALS_H
