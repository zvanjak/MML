///////////////////////////////////////////////////////////////////////////////////////////
// HypothesisTesting.h
// 
// Statistical hypothesis testing functions for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Requires: Statistics.h (for basic stats functions and distributions)
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_HYPOTHESIS_TESTING_H
#define MML_HYPOTHESIS_TESTING_H

#include "../Statistics.h"
#include "CoreDistributions.h"

namespace MML
{
	namespace Statistics
	{
		/**
		 * @brief Result structure for hypothesis tests
		 * 
		 * Contains all relevant information from a hypothesis test:
		 * - Test statistic and its distribution
		 * - P-value and decision at given confidence level
		 * - Degrees of freedom if applicable
		 */
		struct HypothesisTestResult
		{
			Real testStatistic;      ///< Computed test statistic (t, z, chi-square, F, etc.)
			Real pValue;             ///< P-value for the test
			Real criticalValue;      ///< Critical value at the given significance level
			bool rejectNull;         ///< Whether to reject the null hypothesis
			Real confidenceLevel;    ///< Confidence level used (e.g., 0.95 for 95%)
			int degreesOfFreedom;    ///< Degrees of freedom (if applicable, -1 otherwise)
			std::string testName;    ///< Descriptive name of the test

			/**
			 * @brief Construct a hypothesis test result
			 */
			HypothesisTestResult(
				Real statistic = 0.0,
				Real p = 0.0,
				Real critical = 0.0,
				bool reject = false,
				Real confidence = 0.95,
				int df = -1,
				const std::string& name = "Hypothesis Test"
			) : testStatistic(statistic),
			    pValue(p),
			    criticalValue(critical),
			    rejectNull(reject),
			    confidenceLevel(confidence),
			    degreesOfFreedom(df),
			    testName(name)
			{}
		};

		/**
		 * @brief One-sample t-test
		 * 
		 * Tests the null hypothesis H0: μ = μ0 against H1: μ ≠ μ0
		 * where μ is the population mean.
		 * 
		 * Test statistic: t = (x̄ - μ0) / (s / √n)
		 * where x̄ is sample mean, s is sample standard deviation, n is sample size
		 * 
		 * @param sample Vector of sample values
		 * @param mu0 Hypothesized population mean
		 * @param alpha Significance level (default: 0.05 for 95% confidence)
		 * @return HypothesisTestResult with test details
		 */
		inline HypothesisTestResult OneSampleTTest(
			const Vector<Real>& sample,
			Real mu0,
			Real alpha = 0.05
		)
		{
			if (sample.size() < 2)
				throw StatisticsError("OneSampleTTest requires at least 2 samples");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Compute sample statistics
			Real sampleMean = Mean(sample);
			Real sampleStd = StdDev(sample);
			int n = static_cast<int>(sample.size());
			int df = n - 1;

			// Handle zero variance edge case
			if (sampleStd == 0.0) {
				// All values are identical
				Real diff = std::abs(sampleMean - mu0);
				if (diff < std::numeric_limits<Real>::epsilon()) {
					// Perfect match with mu0 - cannot reject null
					return HypothesisTestResult(0.0, 1.0, 0.0, false, 1.0 - alpha, df, "One-Sample t-Test");
				}
				else {
					// All values identical but different from mu0 - strong evidence against null
					// Return large finite t-statistic instead of infinity
					Real largeT = (sampleMean > mu0) ? 1000.0 : -1000.0;
					return HypothesisTestResult(largeT, 0.0, 0.0, true, 1.0 - alpha, df, "One-Sample t-Test");
				}
			}

			// Compute t-statistic
			Real standardError = sampleStd / std::sqrt(static_cast<Real>(n));
			Real tStatistic = (sampleMean - mu0) / standardError;

			// Create t-distribution with appropriate degrees of freedom
			TDistribution tDist(df);

			// Compute two-tailed p-value
			Real pValue = tDist.twoTailedPValue(tStatistic);

			// Critical value for two-tailed test
			Real criticalValue = tDist.criticalValue(alpha);

			// Decision: reject if |t| > critical value (or equivalently, p < alpha)
			bool rejectNull = std::abs(tStatistic) > criticalValue;

			return HypothesisTestResult(
				tStatistic,
				pValue,
				criticalValue,
				rejectNull,
				1.0 - alpha,
				df,
				"One-Sample t-Test"
			);
		}

		/**
		 * @brief Two-sample t-test with equal variances (pooled)
		 * 
		 * Tests H0: μ1 = μ2 against H1: μ1 ≠ μ2
		 * Assumes equal population variances (use Welch's t-test if unequal).
		 * 
		 * Test statistic: t = (x̄1 - x̄2) / (sp · √(1/n1 + 1/n2))
		 * where sp² = ((n1-1)s1² + (n2-1)s2²) / (n1 + n2 - 2) is pooled variance
		 * 
		 * @param sample1 First sample
		 * @param sample2 Second sample
		 * @param alpha Significance level (default: 0.05)
		 * @return HypothesisTestResult with test details
		 */
		inline HypothesisTestResult TwoSampleTTest(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real alpha = 0.05
		)
		{
			if (sample1.size() < 2 || sample2.size() < 2)
				throw StatisticsError("TwoSampleTTest requires at least 2 samples in each group");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Compute sample statistics
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

			// t-statistic
			Real tStatistic = (mean1 - mean2) / standardError;

			// Degrees of freedom
			int df = n1 + n2 - 2;

			// Create t-distribution
			TDistribution tDist(df);

			// Two-tailed p-value
			Real pValue = tDist.twoTailedPValue(tStatistic);

			// Critical value
			Real criticalValue = tDist.criticalValue(alpha);

			// Decision
			bool rejectNull = std::abs(tStatistic) > criticalValue;

			return HypothesisTestResult(
				tStatistic,
				pValue,
				criticalValue,
				rejectNull,
				1.0 - alpha,
				df,
				"Two-Sample t-Test (Pooled)"
			);
		}

		/**
		 * @brief Welch's t-test (unequal variances)
		 * 
		 * Tests H0: μ1 = μ2 against H1: μ1 ≠ μ2
		 * Does NOT assume equal variances (more robust than pooled t-test).
		 * 
		 * Test statistic: t = (x̄1 - x̄2) / √(s1²/n1 + s2²/n2)
		 * Degrees of freedom use Welch-Satterthwaite equation
		 * 
		 * @param sample1 First sample
		 * @param sample2 Second sample
		 * @param alpha Significance level (default: 0.05)
		 * @return HypothesisTestResult with test details
		 */
		inline HypothesisTestResult WelchTTest(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real alpha = 0.05
		)
		{
			if (sample1.size() < 2 || sample2.size() < 2)
				throw StatisticsError("WelchTTest requires at least 2 samples in each group");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Compute sample statistics
			Real mean1 = Mean(sample1);
			Real mean2 = Mean(sample2);
			Real var1 = Variance(sample1);
			Real var2 = Variance(sample2);
			int n1 = static_cast<int>(sample1.size());
			int n2 = static_cast<int>(sample2.size());

			// Standard error (Welch's formula)
			Real se1 = var1 / n1;
			Real se2 = var2 / n2;
			Real standardError = std::sqrt(se1 + se2);

			// t-statistic
			Real tStatistic = (mean1 - mean2) / standardError;

			// Welch-Satterthwaite degrees of freedom
			Real numerator = (se1 + se2) * (se1 + se2);
			Real denominator = (se1 * se1) / (n1 - 1) + (se2 * se2) / (n2 - 1);
			int df = static_cast<int>(std::floor(numerator / denominator));

			// Ensure df is at least 1
			if (df < 1) df = 1;

			// Create t-distribution
			TDistribution tDist(df);

			// Two-tailed p-value
			Real pValue = tDist.twoTailedPValue(tStatistic);

			// Critical value
			Real criticalValue = tDist.criticalValue(alpha);

			// Decision
			bool rejectNull = std::abs(tStatistic) > criticalValue;

			return HypothesisTestResult(
				tStatistic,
				pValue,
				criticalValue,
				rejectNull,
				1.0 - alpha,
				df,
				"Welch's t-Test (Unequal Variances)"
			);
		}

		/**
		 * @brief Paired t-test
		 * 
		 * Tests H0: μd = 0 against H1: μd ≠ 0
		 * where μd is the mean difference between paired observations.
		 * 
		 * Equivalent to one-sample t-test on differences.
		 * 
		 * @param before First measurement (e.g., before treatment)
		 * @param after Second measurement (e.g., after treatment)
		 * @param alpha Significance level (default: 0.05)
		 * @return HypothesisTestResult with test details
		 */
		inline HypothesisTestResult PairedTTest(
			const Vector<Real>& before,
			const Vector<Real>& after,
			Real alpha = 0.05
		)
		{
			if (before.size() != after.size())
				throw StatisticsError("PairedTTest requires equal sample sizes");

			if (before.size() < 2)
				throw StatisticsError("PairedTTest requires at least 2 pairs");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Compute differences
			Vector<Real> differences(before.size());
			for (size_t i = 0; i < before.size(); i++)
				differences[i] = after[i] - before[i];

			// Perform one-sample t-test on differences (testing mean = 0)
			auto result = OneSampleTTest(differences, 0.0, alpha);

			// Update test name
			result.testName = "Paired t-Test";

			return result;
		}

		/**
		 * @brief Chi-Square Goodness-of-Fit Test
		 * 
		 * Tests whether observed frequencies match expected frequencies.
		 * H0: Data follows the expected distribution
		 * H1: Data does not follow the expected distribution
		 * 
		 * Test statistic: χ² = Σ((O_i - E_i)² / E_i)
		 * 
		 * @param observed Vector of observed frequencies
		 * @param expected Vector of expected frequencies (must sum to same as observed)
		 * @param alpha Significance level (default: 0.05)
		 * @return HypothesisTestResult with test details
		 * @throws StatisticsError if sizes don't match, frequencies invalid, or expected has zeros
		 */
		inline HypothesisTestResult ChiSquareGoodnessOfFit(
			const Vector<Real>& observed,
			const Vector<Real>& expected,
			Real alpha = 0.05
		)
		{
			int n = observed.size();
			if (n < 2)
				throw StatisticsError("ChiSquareGoodnessOfFit requires at least 2 categories");

			if (expected.size() != n)
				throw StatisticsError("Observed and expected must have the same size");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Validate frequencies
			Real obsSum = 0.0, expSum = 0.0;
			for (int i = 0; i < n; i++) {
				if (observed[i] < 0.0)
					throw StatisticsError("Observed frequencies must be non-negative");
				if (expected[i] <= 0.0)
					throw StatisticsError("Expected frequencies must be positive");
				obsSum += observed[i];
				expSum += expected[i];
			}

			// Check that sums match (within tolerance)
			if (std::abs(obsSum - expSum) > 1e-6 * std::max(obsSum, expSum))
				throw StatisticsError("Observed and expected frequencies must sum to the same value");

			// Compute chi-square statistic
			Real chiSquare = 0.0;
			for (int i = 0; i < n; i++) {
				Real diff = observed[i] - expected[i];
				chiSquare += (diff * diff) / expected[i];
			}

			// Degrees of freedom = number of categories - 1
			int df = n - 1;

			// Create chi-square distribution
			ChiSquareDistribution chiDist(df);

			// Compute p-value (right-tail)
			Real pValue = chiDist.rightTailPValue(chiSquare);

			// Critical value
			Real criticalValue = chiDist.criticalValue(alpha);

			// Decision: reject if χ² > critical value
			bool rejectNull = chiSquare > criticalValue;

			return HypothesisTestResult(
				chiSquare,
				pValue,
				criticalValue,
				rejectNull,
				1.0 - alpha,
				df,
				"Chi-Square Goodness-of-Fit"
			);
		}

		/**
		 * @brief Chi-Square Test of Independence
		 * 
		 * Tests whether two categorical variables are independent.
		 * H0: Variables are independent
		 * H1: Variables are dependent (associated)
		 * 
		 * @param contingencyTable Matrix of observed frequencies (rows × columns)
		 * @param alpha Significance level (default: 0.05)
		 * @return HypothesisTestResult with test details
		 * @throws StatisticsError if table too small or has invalid frequencies
		 */
		inline HypothesisTestResult ChiSquareIndependence(
			const Matrix<Real>& contingencyTable,
			Real alpha = 0.05
		)
		{
			int rows = contingencyTable.RowNum();
			int cols = contingencyTable.ColNum();

			if (rows < 2 || cols < 2)
				throw StatisticsError("ChiSquareIndependence requires at least 2x2 table");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Compute row and column totals
			Vector<Real> rowTotals(rows);
			Vector<Real> colTotals(cols);
			Real grandTotal = 0.0;

			for (int i = 0; i < rows; i++) {
				rowTotals[i] = 0.0;
				for (int j = 0; j < cols; j++) {
					if (contingencyTable(i, j) < 0.0)
						throw StatisticsError("Frequencies must be non-negative");
					rowTotals[i] += contingencyTable(i, j);
					grandTotal += contingencyTable(i, j);
				}
			}

			for (int j = 0; j < cols; j++) {
				colTotals[j] = 0.0;
				for (int i = 0; i < rows; i++) {
					colTotals[j] += contingencyTable(i, j);
				}
			}

			if (grandTotal <= 0.0)
				throw StatisticsError("Table must have positive total frequency");

			// Compute chi-square statistic
			Real chiSquare = 0.0;
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					// Expected frequency under independence
					Real expected = (rowTotals[i] * colTotals[j]) / grandTotal;
					if (expected < 1.0) {
						// Warning: expected frequency < 1 may give unreliable results
						// But we'll continue (user should be aware)
					}
					Real observed = contingencyTable(i, j);
					Real diff = observed - expected;
					chiSquare += (diff * diff) / expected;
				}
			}

			// Degrees of freedom = (rows - 1) * (cols - 1)
			int df = (rows - 1) * (cols - 1);

			// Create chi-square distribution
			ChiSquareDistribution chiDist(df);

			// Compute p-value (right-tail)
			Real pValue = chiDist.rightTailPValue(chiSquare);

			// Critical value
			Real criticalValue = chiDist.criticalValue(alpha);

			// Decision: reject if χ² > critical value
			bool rejectNull = chiSquare > criticalValue;

			return HypothesisTestResult(
				chiSquare,
				pValue,
				criticalValue,
				rejectNull,
				1.0 - alpha,
				df,
				"Chi-Square Independence Test"
			);
		}

		/**
		 * @brief One-Way ANOVA (Analysis of Variance)
		 * 
		 * Tests whether the means of three or more groups are equal.
		 * H0: μ₁ = μ₂ = μ₃ = ... = μₖ (all group means equal)
		 * H1: At least one mean differs
		 * 
		 * F-statistic = (Between-group variance) / (Within-group variance)
		 *             = (SSB/df_between) / (SSW/df_within)
		 * 
		 * @param groups Vector of groups, where each group is a vector of observations
		 * @param alpha Significance level (default: 0.05)
		 * @return HypothesisTestResult with F-statistic and p-value
		 * @throws StatisticsError if fewer than 2 groups, any group too small, or empty groups
		 */
		inline HypothesisTestResult OneWayANOVA(
			const std::vector<Vector<Real>>& groups,
			Real alpha = 0.05
		)
		{
			int k = static_cast<int>(groups.size());  // number of groups
			if (k < 2)
				throw StatisticsError("OneWayANOVA requires at least 2 groups");

			if (alpha <= 0.0 || alpha >= 1.0)
				throw StatisticsError("Significance level alpha must be in (0, 1)");

			// Validate groups and compute total sample size
			int N = 0;  // total number of observations
			for (int i = 0; i < k; i++) {
				if (groups[i].size() < 2)
					throw StatisticsError("OneWayANOVA requires at least 2 observations per group");
				N += static_cast<int>(groups[i].size());
			}

			// Compute group means and grand mean
			Vector<Real> groupMeans(k);
			Vector<int> groupSizes(k);
			Real grandSum = 0.0;

			for (int i = 0; i < k; i++) {
				groupMeans[i] = Mean(groups[i]);
				groupSizes[i] = static_cast<int>(groups[i].size());
				grandSum += groupMeans[i] * groupSizes[i];
			}

			Real grandMean = grandSum / N;

			// Compute Sum of Squares Between groups (SSB)
			Real SSB = 0.0;
			for (int i = 0; i < k; i++) {
				Real diff = groupMeans[i] - grandMean;
				SSB += groupSizes[i] * diff * diff;
			}

			// Compute Sum of Squares Within groups (SSW)
			Real SSW = 0.0;
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < groupSizes[i]; j++) {
					Real diff = groups[i][j] - groupMeans[i];
					SSW += diff * diff;
				}
			}

			// Degrees of freedom
			int df_between = k - 1;           // between groups
			int df_within = N - k;            // within groups
			int df_total = N - 1;             // total

			// Mean squares
			Real MSB = SSB / df_between;      // Mean square between
			Real MSW = SSW / df_within;       // Mean square within

			// F-statistic
			Real F = MSB / MSW;

			// Create F-distribution
			FDistribution fDist(df_between, df_within);

			// Compute p-value (right-tail)
			Real pValue = fDist.rightTailPValue(F);

			// Critical value
			Real criticalValue = fDist.criticalValue(alpha);

			// Decision: reject if F > critical value
			bool rejectNull = F > criticalValue;

			return HypothesisTestResult(
				F,
				pValue,
				criticalValue,
				rejectNull,
				1.0 - alpha,
				df_between,  // Store between-groups df in result
				"One-Way ANOVA"
			);
		}

	} // namespace Statistics
}  // namespace MML

#endif // MML_HYPOTHESIS_TESTING_H
