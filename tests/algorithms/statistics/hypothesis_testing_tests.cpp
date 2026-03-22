#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics.h"
#include "algorithms/Statistics/Distributions.h"
#include "algorithms/Statistics/StatisticsHypothesis.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::HypothesisTests
{
	/*********************************************************************/
	/*****          Hypothesis Testing - Chi-Square Tests           *****/
	/*********************************************************************/

	TEST_CASE("Statistics::ChiSquareDistribution_Properties", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		Statistics::ChiSquareDistribution chi5(5);
		
		// Mean = df
		REQUIRE_THAT(chi5.mean(), WithinAbs(5.0, 0.001));
		
		// Variance = 2*df
		REQUIRE_THAT(chi5.variance(), WithinAbs(10.0, 0.001));
		
		// Stddev = sqrt(2*df)
		REQUIRE_THAT(chi5.stddev(), WithinAbs(std::sqrt(10.0), 0.001));
		
		// PDF at mode (df-2 for df > 2)
		Real mode = 3.0;  // df - 2 = 5 - 2 = 3
		REQUIRE(chi5.pdf(mode) > 0.0);
		
		// CDF properties
		REQUIRE_THAT(chi5.cdf(0.0), WithinAbs(0.0, 0.001));
		REQUIRE(chi5.cdf(100.0) > 0.99);  // Very high x approaches 1
		
		// Right tail
		REQUIRE_THAT(chi5.rightTailPValue(5.0), WithinAbs(1.0 - chi5.cdf(5.0), 0.001));
	}

	TEST_CASE("Statistics::ChiSquareGoodnessOfFit_UniformDie", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// Fair die: 60 rolls, expect 10 per face
		Vector<Real> observed({12.0, 8.0, 11.0, 9.0, 10.0, 10.0});
		Vector<Real> expected({10.0, 10.0, 10.0, 10.0, 10.0, 10.0});
		
		auto result = Statistics::ChiSquareGoodnessOfFit(observed, expected, 0.05);
		
		// Should not reject (data consistent with fair die)
		REQUIRE_FALSE(result.rejectNull);
		REQUIRE(result.pValue > 0.05);
		
		REQUIRE(result.degreesOfFreedom == 5);  // 6 categories - 1
		REQUIRE(result.testName == "Chi-Square Goodness-of-Fit");
	}

	TEST_CASE("Statistics::ChiSquareGoodnessOfFit_BiasedDie", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// Biased die: heavily favors 6
		Vector<Real> observed({5.0, 5.0, 5.0, 5.0, 5.0, 35.0});
		Vector<Real> expected({10.0, 10.0, 10.0, 10.0, 10.0, 10.0});
		
		auto result = Statistics::ChiSquareGoodnessOfFit(observed, expected, 0.05);
		
		// Should reject (clearly biased)
		REQUIRE(result.rejectNull);
		REQUIRE(result.pValue < 0.001);
		REQUIRE(result.testStatistic > 40.0);  // Large chi-square
	}

	TEST_CASE("Statistics::ChiSquareGoodnessOfFit_CustomDistribution", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// Test against non-uniform expected distribution
		Vector<Real> observed({30.0, 40.0, 20.0, 10.0});
		Vector<Real> expected({25.0, 50.0, 15.0, 10.0});  // Different expected proportions
		
		auto result = Statistics::ChiSquareGoodnessOfFit(observed, expected, 0.05);
		
		// Check result structure
		REQUIRE(result.testStatistic >= 0.0);
		REQUIRE(result.pValue >= 0.0);
		REQUIRE(result.pValue <= 1.0);
		REQUIRE(result.degreesOfFreedom == 3);
	}

	TEST_CASE("Statistics::ChiSquareIndependence_Independent", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// 2x2 contingency table: Gender vs Preference (independent)
		// Proportions are similar across rows
		Matrix<Real> table(2, 2);
		table(0, 0) = 30;  // Male, Likes A
		table(0, 1) = 20;  // Male, Likes B
		table(1, 0) = 30;  // Female, Likes A
		table(1, 1) = 20;  // Female, Likes B
		
		auto result = Statistics::ChiSquareIndependence(table, 0.05);
		
		// Perfect independence: should not reject
		REQUIRE_FALSE(result.rejectNull);
		REQUIRE_THAT(result.testStatistic, WithinAbs(0.0, 0.01));
		REQUIRE(result.pValue > 0.90);  // Very high p-value
		
		REQUIRE(result.degreesOfFreedom == 1);  // (2-1)*(2-1) = 1
		REQUIRE(result.testName == "Chi-Square Independence Test");
	}

	TEST_CASE("Statistics::ChiSquareIndependence_Dependent", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// 2x2 contingency table: strong association
		Matrix<Real> table(2, 2);
		table(0, 0) = 40;  // Group A, Outcome Yes
		table(0, 1) = 10;  // Group A, Outcome No
		table(1, 0) = 10;  // Group B, Outcome Yes
		table(1, 1) = 40;  // Group B, Outcome No
		
		auto result = Statistics::ChiSquareIndependence(table, 0.05);
		
		// Strong dependence: should reject
		REQUIRE(result.rejectNull);
		REQUIRE(result.pValue < 0.001);
		REQUIRE(result.testStatistic > 18.0);  // Large chi-square
	}

	TEST_CASE("Statistics::ChiSquareIndependence_3x3Table", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// 3x3 table: Education level vs Income bracket
		Matrix<Real> table(3, 3);
		table(0, 0) = 20; table(0, 1) = 10; table(0, 2) = 5;
		table(1, 0) = 15; table(1, 1) = 25; table(1, 2) = 10;
		table(2, 0) = 5;  table(2, 1) = 15; table(2, 2) = 25;
		
		auto result = Statistics::ChiSquareIndependence(table, 0.05);
		
		// Should show some association
		REQUIRE(result.degreesOfFreedom == 4);  // (3-1)*(3-1) = 4
		REQUIRE(result.testStatistic >= 0.0);
		REQUIRE(result.pValue >= 0.0);
		REQUIRE(result.pValue <= 1.0);
	}

	TEST_CASE("Statistics::ChiSquare_EdgeCases", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		// Minimum size (2 categories)
		Vector<Real> obs2({10.0, 20.0});
		Vector<Real> exp2({15.0, 15.0});
		auto result_min = Statistics::ChiSquareGoodnessOfFit(obs2, exp2, 0.05);
		REQUIRE(result_min.degreesOfFreedom == 1);
		
		// Mismatched sizes should throw
		Vector<Real> obs3({10.0, 20.0, 30.0});
		REQUIRE_THROWS_AS(
			Statistics::ChiSquareGoodnessOfFit(obs3, exp2, 0.05),
			StatisticsError
		);
		
		// Negative observed should throw
		Vector<Real> obs_neg({-5.0, 15.0});
		REQUIRE_THROWS_AS(
			Statistics::ChiSquareGoodnessOfFit(obs_neg, exp2, 0.05),
			StatisticsError
		);
		
		// Zero/negative expected should throw
		Vector<Real> exp_zero({0.0, 30.0});
		REQUIRE_THROWS_AS(
			Statistics::ChiSquareGoodnessOfFit(obs2, exp_zero, 0.05),
			StatisticsError
		);
		
		// Contingency table too small
		Matrix<Real> table_1x2(1, 2);
		table_1x2(0, 0) = 10; table_1x2(0, 1) = 20;
		REQUIRE_THROWS_AS(
			Statistics::ChiSquareIndependence(table_1x2, 0.05),
			StatisticsError
		);
	}

	TEST_CASE("Statistics::ChiSquare_AlphaLevels", "[statistics][hypothesis][chi-square]")
	{
		TEST_PRECISION_INFO();
		
		Vector<Real> observed({15.0, 20.0, 25.0, 40.0});
		Vector<Real> expected({25.0, 25.0, 25.0, 25.0});
		
		auto result_05 = Statistics::ChiSquareGoodnessOfFit(observed, expected, 0.05);
		auto result_01 = Statistics::ChiSquareGoodnessOfFit(observed, expected, 0.01);
		
		// Same test statistic
		REQUIRE_THAT(result_05.testStatistic, WithinAbs(result_01.testStatistic, 0.001));
		
		// Different critical values (0.01 is more stringent)
		REQUIRE(result_01.criticalValue > result_05.criticalValue);
		
		// Same p-value
		REQUIRE_THAT(result_05.pValue, WithinAbs(result_01.pValue, 0.001));
	}

	/*********************************************************************/
	/*****          Hypothesis Testing - ANOVA                       *****/
	/*********************************************************************/

	TEST_CASE("Statistics::FDistribution_Properties", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		Statistics::FDistribution f_5_10(5, 10);
		
		// Mean = df2/(df2-2) for df2 > 2
		REQUIRE_THAT(f_5_10.mean(), WithinAbs(10.0/8.0, 0.001));
		
		// Variance defined for df2 > 4
		REQUIRE(f_5_10.variance() > 0.0);
		
		// PDF is positive for x > 0
		REQUIRE(f_5_10.pdf(1.0) > 0.0);
		REQUIRE(f_5_10.pdf(2.0) > 0.0);
		
		// CDF properties
		REQUIRE_THAT(f_5_10.cdf(0.0), WithinAbs(0.0, 0.001));
		REQUIRE(f_5_10.cdf(100.0) > 0.99);
		
		// Right tail = 1 - CDF
		Real f_val = 2.5;
		REQUIRE_THAT(f_5_10.rightTailPValue(f_val), WithinAbs(1.0 - f_5_10.cdf(f_val), 0.01));
	}

	TEST_CASE("Statistics::OneWayANOVA_NoEffect", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		// Three groups with same mean (around 10), similar variance
		std::vector<Vector<Real>> groups;
		groups.push_back(Vector<Real>({9.0, 10.0, 11.0, 10.5, 9.5}));
		groups.push_back(Vector<Real>({9.5, 10.5, 10.0, 9.0, 11.0}));
		groups.push_back(Vector<Real>({10.0, 9.5, 10.5, 11.0, 9.0}));
		
		auto result = Statistics::OneWayANOVA(groups, 0.05);
		
		// Should not reject (groups have similar means)
		REQUIRE_FALSE(result.rejectNull);
		REQUIRE(result.pValue > 0.05);
		REQUIRE(result.testStatistic >= 0.0);  // F is always non-negative
		
		// Degrees of freedom: between = k-1 = 2, within = N-k = 15-3 = 12
		REQUIRE(result.degreesOfFreedom == 2);
		REQUIRE(result.testName == "One-Way ANOVA");
	}

	TEST_CASE("Statistics::OneWayANOVA_SignificantEffect", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		// Three groups with clearly different means
		std::vector<Vector<Real>> groups;
		groups.push_back(Vector<Real>({10.0, 11.0, 9.0, 10.5, 10.5}));    // mean â‰ˆ 10
		groups.push_back(Vector<Real>({20.0, 21.0, 19.0, 20.5, 20.5}));   // mean â‰ˆ 20
		groups.push_back(Vector<Real>({30.0, 31.0, 29.0, 30.5, 30.5}));   // mean â‰ˆ 30
		
		auto result = Statistics::OneWayANOVA(groups, 0.05);
		
		// Should strongly reject (groups very different)
		REQUIRE(result.rejectNull);
		REQUIRE(result.pValue < 0.001);
		REQUIRE(result.testStatistic > 100.0);  // Very large F
	}

	TEST_CASE("Statistics::OneWayANOVA_TwoGroups", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		// ANOVA with 2 groups is equivalent to two-sample t-test
		std::vector<Vector<Real>> groups;
		groups.push_back(Vector<Real>({10.0, 12.0, 14.0, 16.0, 18.0}));
		groups.push_back(Vector<Real>({20.0, 22.0, 24.0, 26.0, 28.0}));
		
		auto result = Statistics::OneWayANOVA(groups, 0.05);
		
		// Should reject (clearly different means)
		REQUIRE(result.rejectNull);
		REQUIRE(result.pValue < 0.01);
		
		// df between = 2-1 = 1
		REQUIRE(result.degreesOfFreedom == 1);
		
		// For 2 groups, F = tÂ² (approximately)
		auto tResult = Statistics::TwoSampleTTest(groups[0], groups[1], 0.05);
		Real tSquared = tResult.testStatistic * tResult.testStatistic;
		REQUIRE_THAT(result.testStatistic, WithinAbs(tSquared, 0.1));
	}

	TEST_CASE("Statistics::OneWayANOVA_UnequalSampleSizes", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		// Groups with different sample sizes
		std::vector<Vector<Real>> groups;
		groups.push_back(Vector<Real>({10.0, 11.0}));                    // n=2
		groups.push_back(Vector<Real>({15.0, 16.0, 17.0, 18.0}));        // n=4
		groups.push_back(Vector<Real>({20.0, 21.0, 22.0, 23.0, 24.0}));  // n=5
		
		auto result = Statistics::OneWayANOVA(groups, 0.05);
		
		// Should work with unequal sizes
		REQUIRE(result.testStatistic >= 0.0);
		REQUIRE(result.pValue >= 0.0);
		REQUIRE(result.pValue <= 1.0);
		
		// df between = 3-1 = 2, within = 11-3 = 8
		REQUIRE(result.degreesOfFreedom == 2);
	}

	TEST_CASE("Statistics::OneWayANOVA_FourGroups", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		// Four groups with moderate differences
		std::vector<Vector<Real>> groups;
		groups.push_back(Vector<Real>({10.0, 11.0, 12.0}));
		groups.push_back(Vector<Real>({13.0, 14.0, 15.0}));
		groups.push_back(Vector<Real>({16.0, 17.0, 18.0}));
		groups.push_back(Vector<Real>({19.0, 20.0, 21.0}));
		
		auto result = Statistics::OneWayANOVA(groups, 0.05);
		
		// df between = 4-1 = 3
		REQUIRE(result.degreesOfFreedom == 3);
		REQUIRE(result.testStatistic >= 0.0);
	}

	TEST_CASE("Statistics::OneWayANOVA_EdgeCases", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		// Single group should throw
		std::vector<Vector<Real>> single_group;
		single_group.push_back(Vector<Real>({10.0, 11.0, 12.0}));
		REQUIRE_THROWS_AS(
			Statistics::OneWayANOVA(single_group, 0.05),
			StatisticsError
		);
		
		// Group with only 1 observation should throw
		std::vector<Vector<Real>> tiny_group;
		tiny_group.push_back(Vector<Real>({10.0, 11.0}));
		tiny_group.push_back(Vector<Real>({15.0}));  // Only 1
		REQUIRE_THROWS_AS(
			Statistics::OneWayANOVA(tiny_group, 0.05),
			StatisticsError
		);
		
		// Invalid alpha
		std::vector<Vector<Real>> valid_groups;
		valid_groups.push_back(Vector<Real>({10.0, 11.0}));
		valid_groups.push_back(Vector<Real>({15.0, 16.0}));
		REQUIRE_THROWS_AS(
			Statistics::OneWayANOVA(valid_groups, 0.0),
			StatisticsError
		);
		REQUIRE_THROWS_AS(
			Statistics::OneWayANOVA(valid_groups, 1.0),
			StatisticsError
		);
	}

	TEST_CASE("Statistics::OneWayANOVA_AlphaLevels", "[statistics][hypothesis][anova]")
	{
		TEST_PRECISION_INFO();
		
		std::vector<Vector<Real>> groups;
		groups.push_back(Vector<Real>({10.0, 11.0, 12.0, 13.0}));
		groups.push_back(Vector<Real>({14.0, 15.0, 16.0, 17.0}));
		groups.push_back(Vector<Real>({18.0, 19.0, 20.0, 21.0}));
		
		auto result_05 = Statistics::OneWayANOVA(groups, 0.05);
		auto result_01 = Statistics::OneWayANOVA(groups, 0.01);
		
		// Same F-statistic
		REQUIRE_THAT(result_05.testStatistic, WithinAbs(result_01.testStatistic, 0.001));
		
		// Same p-value
		REQUIRE_THAT(result_05.pValue, WithinAbs(result_01.pValue, 0.001));
		
		// Different critical values (0.01 more stringent)
		REQUIRE(result_01.criticalValue > result_05.criticalValue);
	}

	// ================================================================================
	// HYPOTHESIS TESTING: T-TESTS
	// ================================================================================

	TEST_CASE("Statistics::OneSampleTTest_BasicUsage", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Sample data: [10, 12, 15, 18, 20]
		// Mean = 15, std â‰ˆ 4.18, n = 5
		// Test H0: Î¼ = 15 (should NOT reject)
		Vector<Real> sample({10.0, 12.0, 15.0, 18.0, 20.0});
		
		auto result = Statistics::OneSampleTTest(sample, 15.0, 0.05);
		
		// Testing against the true mean should give t â‰ˆ 0
		REQUIRE_THAT(result.testStatistic, WithinAbs(0.0, 0.01));
		
		// P-value should be close to 1 (no evidence against H0)
		REQUIRE(result.pValue > 0.9);
		
		// Should NOT reject null hypothesis
		REQUIRE_FALSE(result.rejectNull);
		
		// Check metadata
		REQUIRE(result.confidenceLevel == 0.95);
		REQUIRE(result.degreesOfFreedom == 4);
		REQUIRE(result.testName == "One-Sample t-Test");
	}

	TEST_CASE("Statistics::OneSampleTTest_RejectNull", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Sample clearly from population with mean > 15
		// [20, 22, 25, 28, 30] mean = 25
		Vector<Real> sample({20.0, 22.0, 25.0, 28.0, 30.0});
		
		// Test H0: Î¼ = 15 (should reject)
		auto result = Statistics::OneSampleTTest(sample, 15.0, 0.05);
		
		// Large positive t-statistic (sample mean >> hypothesized mean)
		REQUIRE(result.testStatistic > 3.0);
		
		// Very small p-value
		REQUIRE(result.pValue < 0.05);
		
		// Should reject null hypothesis
		REQUIRE(result.rejectNull);
	}

	TEST_CASE("Statistics::OneSampleTTest_KnownValues", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Example from statistics textbook
		// Sample: [5, 7, 9, 11, 13] mean = 9, std = sqrt(10) â‰ˆ 3.162, n = 5
		// Test H0: Î¼ = 10
		Vector<Real> sample({5.0, 7.0, 9.0, 11.0, 13.0});
		
		Real sampleMean = Statistics::Mean(sample);
		Real sampleStd = Statistics::StdDev(sample);
		int n = 5;
		
		// Manual calculation: t = (9 - 10) / (3.162 / sqrt(5)) = -1 / 1.414 â‰ˆ -0.707
		Real expectedT = (sampleMean - 10.0) / (sampleStd / std::sqrt(static_cast<Real>(n)));
		
		auto result = Statistics::OneSampleTTest(sample, 10.0, 0.05);
		
		REQUIRE_THAT(result.testStatistic, WithinAbs(expectedT, 0.01));
		
		// df = 4, critical value at Î±=0.05 is Â±2.776
		REQUIRE_THAT(result.criticalValue, WithinAbs(2.776, 0.01));
		
		// |t| < critical, so don't reject
		REQUIRE_FALSE(result.rejectNull);
	}

	TEST_CASE("Statistics::TwoSampleTTest_EqualMeans", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Two samples from same distribution
		Vector<Real> sample1({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> sample2({11.0, 13.0, 15.0, 17.0, 19.0});
		
		// Both have mean â‰ˆ 14-15, should not reject H0: Î¼1 = Î¼2
		auto result = Statistics::TwoSampleTTest(sample1, sample2, 0.05);
		
		// t-statistic should be small
		REQUIRE(std::abs(result.testStatistic) < 1.0);
		
		// Should NOT reject null
		REQUIRE_FALSE(result.rejectNull);
		
		// P-value should be large
		REQUIRE(result.pValue > 0.1);
		
		// df = n1 + n2 - 2 = 8
		REQUIRE(result.degreesOfFreedom == 8);
	}

	TEST_CASE("Statistics::TwoSampleTTest_DifferentMeans", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Clearly different populations
		Vector<Real> sample1({10.0, 12.0, 14.0, 16.0, 18.0});  // mean â‰ˆ 14
		Vector<Real> sample2({20.0, 22.0, 24.0, 26.0, 28.0});  // mean â‰ˆ 24
		
		// Should reject H0: Î¼1 = Î¼2
		auto result = Statistics::TwoSampleTTest(sample1, sample2, 0.05);
		
		// Large negative t-statistic (sample1 < sample2)
		REQUIRE(result.testStatistic <= -5.0);
		
		// Very small p-value
		REQUIRE(result.pValue < 0.01);
		
		// Should reject null
		REQUIRE(result.rejectNull);
		
		REQUIRE(result.testName == "Two-Sample t-Test (Pooled)");
	}

	TEST_CASE("Statistics::TwoSampleTTest_UnequalSizes", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Different sample sizes
		Vector<Real> sample1({10.0, 15.0, 20.0});           // n1 = 3
		Vector<Real> sample2({12.0, 14.0, 16.0, 18.0});    // n2 = 4
		
		// Should handle unequal sizes correctly
		auto result = Statistics::TwoSampleTTest(sample1, sample2, 0.05);
		
		// df = 3 + 4 - 2 = 5
		REQUIRE(result.degreesOfFreedom == 5);
		
		// Computation should complete without error
		REQUIRE(std::isfinite(result.testStatistic));
		REQUIRE(std::isfinite(result.pValue));
	}

	TEST_CASE("Statistics::WelchTTest_BasicUsage", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Two samples with potentially different variances
		Vector<Real> sample1({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> sample2({15.0, 20.0, 25.0, 30.0, 35.0});  // Higher variance
		
		auto result = Statistics::WelchTTest(sample1, sample2, 0.05);
		
		// Welch's test should handle unequal variances
		REQUIRE(std::isfinite(result.testStatistic));
		REQUIRE(result.pValue >= 0.0);
		REQUIRE(result.pValue <= 1.0);
		
		// df should be computed via Welch-Satterthwaite
		REQUIRE(result.degreesOfFreedom >= 1);
		REQUIRE(result.degreesOfFreedom < 8);  // Less than n1 + n2 - 2
		
		REQUIRE(result.testName == "Welch's t-Test (Unequal Variances)");
	}

	TEST_CASE("Statistics::WelchTTest_VsPooled", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// When variances are equal, Welch and pooled should give similar results
		Vector<Real> sample1({10.0, 11.0, 12.0, 13.0, 14.0});
		Vector<Real> sample2({15.0, 16.0, 17.0, 18.0, 19.0});
		
		auto pooled = Statistics::TwoSampleTTest(sample1, sample2, 0.05);
		auto welch = Statistics::WelchTTest(sample1, sample2, 0.05);
		
		// Test statistics should be similar (but not identical)
		REQUIRE_THAT(pooled.testStatistic, WithinAbs(welch.testStatistic, 0.5));
		
		// Both should reach same conclusion
		REQUIRE(pooled.rejectNull == welch.rejectNull);
	}

	TEST_CASE("Statistics::PairedTTest_NoEffect", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Before and after measurements with no real change
		Vector<Real> before({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> after({11.0, 13.0, 15.0, 17.0, 19.0});  // +1 consistently
		
		// Mean difference = 1.0, but we need to check if it's significant
		auto result = Statistics::PairedTTest(before, after, 0.05);
		
		// Perfect pairing should detect the difference
		// But let's test a case with no difference first
		Vector<Real> after_same = before;
		auto result_same = Statistics::PairedTTest(before, after_same, 0.05);
		
		// No difference: t â‰ˆ 0
		REQUIRE_THAT(result_same.testStatistic, WithinAbs(0.0, 0.01));
		REQUIRE_FALSE(result_same.rejectNull);
		
		REQUIRE(result_same.testName == "Paired t-Test");
	}

	TEST_CASE("Statistics::PairedTTest_SignificantChange", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Clear improvement after treatment
		Vector<Real> before({100.0, 105.0, 110.0, 115.0, 120.0});
		Vector<Real> after({110.0, 115.0, 120.0, 125.0, 130.0});  // Consistent +10
		
		auto result = Statistics::PairedTTest(before, after, 0.05);
		
		// Should detect the significant change
		REQUIRE(result.testStatistic > 0.0);  // Positive change
		REQUIRE(result.pValue < 0.05);
		REQUIRE(result.rejectNull);
		
		// df = n - 1 = 4
		REQUIRE(result.degreesOfFreedom == 4);
	}

	TEST_CASE("Statistics::PairedTTest_VsTwoSample", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Paired test is more powerful than two-sample when data is paired
		Vector<Real> before({10.0, 50.0, 30.0, 70.0, 20.0});
		Vector<Real> after({15.0, 55.0, 35.0, 75.0, 25.0});  // Each increased by 5
		
		auto paired = Statistics::PairedTTest(before, after, 0.05);
		auto twoSample = Statistics::TwoSampleTTest(before, after, 0.05);
		
		// Paired test should detect the consistent +5 difference
		// Two-sample test might not due to high variance
		
		// Both should have finite results
		REQUIRE(std::isfinite(paired.testStatistic));
		REQUIRE(std::isfinite(twoSample.testStatistic));
		
		// Paired test typically more sensitive (higher |t|)
		// (though not guaranteed in all cases)
	}

	TEST_CASE("Statistics::TTest_EdgeCases", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Minimum sample size
		Vector<Real> tiny_sample({10.0, 20.0});
		
		auto result = Statistics::OneSampleTTest(tiny_sample, 15.0, 0.05);
		REQUIRE(result.degreesOfFreedom == 1);
		REQUIRE(std::isfinite(result.testStatistic));
		
		// Should throw with insufficient data
		Vector<Real> too_small({10.0});
		REQUIRE_THROWS_AS(
			Statistics::OneSampleTTest(too_small, 15.0, 0.05),
			StatisticsError
		);
		
		// Invalid alpha
		REQUIRE_THROWS_AS(
			Statistics::OneSampleTTest(tiny_sample, 15.0, 0.0),
			StatisticsError
		);
		
		REQUIRE_THROWS_AS(
			Statistics::OneSampleTTest(tiny_sample, 15.0, 1.0),
			StatisticsError
		);
	}

	TEST_CASE("Statistics::PairedTTest_EdgeCases", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		Vector<Real> before({10.0, 20.0, 30.0});
		Vector<Real> after({15.0, 25.0});  // Different size
		
		// Should throw with mismatched sizes
		REQUIRE_THROWS_AS(
			Statistics::PairedTTest(before, after, 0.05),
			StatisticsError
		);
		
		// Too small
		Vector<Real> tiny_before({10.0});
		Vector<Real> tiny_after({15.0});
		
		REQUIRE_THROWS_AS(
			Statistics::PairedTTest(tiny_before, tiny_after, 0.05),
			StatisticsError
		);
	}

	TEST_CASE("Statistics::TTest_AlphaLevels", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		Vector<Real> sample({10.0, 12.0, 14.0, 16.0, 18.0});
		
		// Different alpha levels should give different critical values
		auto result_05 = Statistics::OneSampleTTest(sample, 20.0, 0.05);
		auto result_01 = Statistics::OneSampleTTest(sample, 20.0, 0.01);
		
		// More stringent alpha (0.01) requires larger critical value
		REQUIRE(result_01.criticalValue > result_05.criticalValue);
		
		// Confidence levels should match
		REQUIRE(result_05.confidenceLevel == 0.95);
		REQUIRE(result_01.confidenceLevel == 0.99);
		
		// Same test statistic regardless of alpha
		REQUIRE_THAT(result_05.testStatistic, WithinAbs(result_01.testStatistic, 0.001));
	}

	/******************************************************************************/
	/*****             Detailed API Tests - Hypothesis Testing                *****/
	/******************************************************************************/

	TEST_CASE("Statistics::OneSampleTTestDetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> sample({10.0, 12.0, 14.0, 16.0, 18.0});
		auto simple = Statistics::OneSampleTTest(sample, 15.0, 0.05);
		auto detailed = Statistics::OneSampleTTestDetailed(sample, 15.0, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "OneSampleTTest");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
		REQUIRE_THAT(detailed.criticalValue, WithinAbs(simple.criticalValue, 1e-15));
		REQUIRE(detailed.rejectNull == simple.rejectNull);
		REQUIRE(detailed.degreesOfFreedom == simple.degreesOfFreedom);
		REQUIRE(detailed.elapsed_time_ms >= 0.0);
	}

	TEST_CASE("Statistics::TwoSampleTTestDetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> s1({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> s2({11.0, 13.0, 15.0, 17.0, 19.0});
		auto simple = Statistics::TwoSampleTTest(s1, s2, 0.05);
		auto detailed = Statistics::TwoSampleTTestDetailed(s1, s2, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "TwoSampleTTest");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
		REQUIRE(detailed.rejectNull == simple.rejectNull);
	}

	TEST_CASE("Statistics::WelchTTestDetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> s1({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> s2({15.0, 17.0, 19.0});
		auto simple = Statistics::WelchTTest(s1, s2, 0.05);
		auto detailed = Statistics::WelchTTestDetailed(s1, s2, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "WelchTTest");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
	}

	TEST_CASE("Statistics::PairedTTestDetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> before({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> after({11.0, 13.0, 15.0, 17.0, 19.0});
		auto simple = Statistics::PairedTTest(before, after, 0.05);
		auto detailed = Statistics::PairedTTestDetailed(before, after, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "PairedTTest");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
	}

	TEST_CASE("Statistics::ChiSquareGoodnessOfFitDetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> observed({20.0, 15.0, 25.0, 18.0, 22.0, 20.0});
		Vector<Real> expected({20.0, 20.0, 20.0, 20.0, 20.0, 20.0});
		auto simple = Statistics::ChiSquareGoodnessOfFit(observed, expected, 0.05);
		auto detailed = Statistics::ChiSquareGoodnessOfFitDetailed(observed, expected, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ChiSquareGoodnessOfFit");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
	}

	TEST_CASE("Statistics::ChiSquareIndependenceDetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		Matrix<Real> table(2, 2);
		table(0, 0) = 50; table(0, 1) = 30;
		table(1, 0) = 20; table(1, 1) = 100;
		auto simple = Statistics::ChiSquareIndependence(table, 0.05);
		auto detailed = Statistics::ChiSquareIndependenceDetailed(table, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ChiSquareIndependence");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
	}

	TEST_CASE("Statistics::OneWayANOVADetailed - values match simple API", "[statistics][hypothesis][Detailed]")
	{
		std::vector<Vector<Real>> groups = {
			Vector<Real>({10.0, 12.0, 14.0}),
			Vector<Real>({15.0, 17.0, 19.0}),
			Vector<Real>({20.0, 22.0, 24.0})
		};
		auto simple = Statistics::OneWayANOVA(groups, 0.05);
		auto detailed = Statistics::OneWayANOVADetailed(groups, 0.05);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "OneWayANOVA");
		REQUIRE_THAT(detailed.testStatistic, WithinAbs(simple.testStatistic, 1e-15));
		REQUIRE_THAT(detailed.pValue, WithinAbs(simple.pValue, 1e-15));
	}

	TEST_CASE("Statistics::OneSampleTTestDetailed - ConvertToStatus on invalid input", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> empty;
		Statistics::StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = Statistics::OneSampleTTestDetailed(empty, 0.0, 0.05, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(result.error_message.empty());
	}

	TEST_CASE("Statistics::OneSampleTTestDetailed - Propagate policy throws", "[statistics][hypothesis][Detailed]")
	{
		Vector<Real> empty;
		Statistics::StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::Propagate;

		REQUIRE_THROWS(Statistics::OneSampleTTestDetailed(empty, 0.0, 0.05, config));
	}

	TEST_CASE("Statistics::OneWayANOVADetailed - ConvertToStatus on invalid input", "[statistics][hypothesis][Detailed]")
	{
		std::vector<Vector<Real>> oneGroup = { Vector<Real>({1.0, 2.0, 3.0}) };
		Statistics::StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = Statistics::OneWayANOVADetailed(oneGroup, 0.05, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
	}

}