#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics.h"
#include "algorithms/Statistics/Distributions.h"
#include "algorithms/Statistics/StatisticsConfidence.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Statistics::ConfidenceIntervalsTests
{

	// ================================================================================
	// PHASE 5: CONFIDENCE INTERVALS
	// ================================================================================

	TEST_CASE("Statistics - Confidence Interval for Mean", "[statistics][confidence][ci_mean]")
	{
		using namespace MML::Statistics;

		// Basic 95% CI for mean
		Vector<Real> sample({10.0, 12.0, 14.0, 16.0, 18.0});
		ConfidenceInterval ci = ConfidenceIntervalMean(sample, 0.95);

		REQUIRE_THAT(ci.estimate, WithinAbs(14.0, 0.001));  // Mean
		REQUIRE(ci.parameter == "Mean");
		REQUIRE_THAT(ci.confidenceLevel, WithinAbs(0.95, 0.001));
		REQUIRE(ci.lowerBound < ci.estimate);
		REQUIRE(ci.upperBound > ci.estimate);
		REQUIRE_THAT(ci.marginOfError, WithinAbs(ci.estimate - ci.lowerBound, 1e-6));
		REQUIRE_THAT(ci.marginOfError, WithinAbs(ci.upperBound - ci.estimate, 1e-6));

		// 99% CI should be wider than 95%
		ConfidenceInterval ci99 = ConfidenceIntervalMean(sample, 0.99);
		REQUIRE(ci99.marginOfError > ci.marginOfError);
		REQUIRE(ci99.lowerBound < ci.lowerBound);
		REQUIRE(ci99.upperBound > ci.upperBound);

		// Larger sample should give narrower CI
		Vector<Real> largeSample({10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0});
		ConfidenceInterval ciLarge = ConfidenceIntervalMean(largeSample, 0.95);
		REQUIRE(ciLarge.marginOfError < ci.marginOfError);
	}

	TEST_CASE("Statistics - Confidence Interval for Mean Difference", "[statistics][confidence][ci_mean_diff]")
	{
		using namespace MML::Statistics;

		// Two samples with different means
		Vector<Real> sample1({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> sample2({15.0, 17.0, 19.0, 21.0, 23.0});

		ConfidenceInterval ci = ConfidenceIntervalMeanDifference(sample1, sample2, 0.95);

		REQUIRE(ci.parameter == "Mean Difference");
		REQUIRE_THAT(ci.confidenceLevel, WithinAbs(0.95, 0.001));
		REQUIRE_THAT(ci.estimate, WithinAbs(-5.0, 0.001));  // mean1 - mean2 = -5
		REQUIRE(ci.lowerBound < 0.0);  // Difference should be negative
		REQUIRE(ci.upperBound < 0.0);  // Both bounds negative since sample2 > sample1

		// Identical samples should give CI containing zero
		Vector<Real> identical1({10.0, 12.0, 14.0});
		Vector<Real> identical2({10.1, 11.9, 14.0});  // Very similar
		ConfidenceInterval ciZero = ConfidenceIntervalMeanDifference(identical1, identical2, 0.95);
		REQUIRE(ciZero.lowerBound < 0.1);
		REQUIRE(ciZero.upperBound > -0.1);
	}

	TEST_CASE("Statistics - Confidence Interval for Proportion", "[statistics][confidence][ci_proportion]")
	{
		using namespace MML::Statistics;

		// 50% success rate
		ConfidenceInterval ci = ConfidenceIntervalProportion(50, 100, 0.95);

		REQUIRE(ci.parameter == "Proportion");
		REQUIRE_THAT(ci.estimate, WithinAbs(0.5, 0.001));
		REQUIRE_THAT(ci.confidenceLevel, WithinAbs(0.95, 0.001));
		REQUIRE(ci.lowerBound > 0.4);
		REQUIRE(ci.upperBound < 0.6);
		REQUIRE(ci.lowerBound >= 0.0);  // Clamped to [0,1]
		REQUIRE(ci.upperBound <= 1.0);

		// High success rate (90%)
		ConfidenceInterval ciHigh = ConfidenceIntervalProportion(90, 100, 0.95);
		REQUIRE_THAT(ciHigh.estimate, WithinAbs(0.9, 0.001));
		REQUIRE(ciHigh.lowerBound > 0.8);
		REQUIRE(ciHigh.upperBound <= 1.0);  // Should be clamped at 1.0

		// Low success rate (10%)
		ConfidenceInterval ciLow = ConfidenceIntervalProportion(10, 100, 0.95);
		REQUIRE_THAT(ciLow.estimate, WithinAbs(0.1, 0.001));
		REQUIRE(ciLow.lowerBound >= 0.0);  // Should be clamped at 0.0
		REQUIRE(ciLow.upperBound < 0.2);

		// Larger sample gives narrower CI
		ConfidenceInterval ciLarge = ConfidenceIntervalProportion(500, 1000, 0.95);
		REQUIRE_THAT(ciLarge.estimate, WithinAbs(0.5, 0.001));
		REQUIRE(ciLarge.marginOfError < ci.marginOfError);
	}

	TEST_CASE("Statistics - Confidence Interval for Proportion Difference", "[statistics][confidence][ci_prop_diff]")
	{
		using namespace MML::Statistics;

		// Two proportions
		ConfidenceInterval ci = ConfidenceIntervalProportionDifference(
			60, 100,  // p1 = 0.6
			40, 100,  // p2 = 0.4
			0.95
		);

		REQUIRE(ci.parameter == "Proportion Difference");
		REQUIRE_THAT(ci.estimate, WithinAbs(0.2, 0.001));  // 0.6 - 0.4
		REQUIRE(ci.lowerBound > 0.0);  // Positive difference
		REQUIRE(ci.upperBound < 0.4);
		REQUIRE_THAT(ci.confidenceLevel, WithinAbs(0.95, 0.001));

		// Similar proportions should include zero
		ConfidenceInterval ciSimilar = ConfidenceIntervalProportionDifference(
			50, 100,  // p1 = 0.5
			48, 100,  // p2 = 0.48
			0.95
		);
		REQUIRE(ciSimilar.lowerBound < 0.0);
		REQUIRE(ciSimilar.upperBound > 0.0);  // CI contains zero
	}

	TEST_CASE("Statistics - Confidence Interval for Paired Differences", "[statistics][confidence][ci_paired]")
	{
		using namespace MML::Statistics;

		// Before-after measurements
		Vector<Real> before({100.0, 105.0, 110.0, 115.0, 120.0});
		Vector<Real> after ({110.0, 115.0, 125.0, 130.0, 135.0});

		ConfidenceInterval ci = ConfidenceIntervalPairedDifference(before, after, 0.95);

		REQUIRE(ci.parameter == "Mean");  // Uses mean CI internally
		REQUIRE(ci.estimate > 0.0);  // Positive difference (after > before)
		REQUIRE(ci.lowerBound > 0.0);  // All differences positive
		REQUIRE(ci.upperBound > ci.estimate);
		REQUIRE_THAT(ci.confidenceLevel, WithinAbs(0.95, 0.001));

		// Check mean difference
		Real meanDiff = Mean(after) - Mean(before);
		REQUIRE_THAT(ci.estimate, WithinAbs(meanDiff, 0.1));

		// No change should give CI containing zero
		Vector<Real> noChange1({10.0, 12.0, 14.0});
		Vector<Real> noChange2({10.1, 11.9, 14.0});
		ConfidenceInterval ciNoChange = ConfidenceIntervalPairedDifference(noChange1, noChange2, 0.95);
		REQUIRE(ciNoChange.lowerBound < 0.1);
		REQUIRE(ciNoChange.upperBound > -0.1);
	}

	TEST_CASE("Statistics - Confidence Interval Edge Cases", "[statistics][confidence][edge_cases]")
	{
		using namespace MML::Statistics;

		// Small sample (n=2)
		Vector<Real> small({10.0, 20.0});
		ConfidenceInterval ciSmall = ConfidenceIntervalMean(small, 0.95);
		REQUIRE_THAT(ciSmall.estimate, WithinAbs(15.0, 0.001));
		REQUIRE(ciSmall.marginOfError > 5.0);  // Very wide for small sample

		// Too small sample should throw
		Vector<Real> tooSmall({10.0});
		REQUIRE_THROWS_AS(
			ConfidenceIntervalMean(tooSmall, 0.95),
			StatisticsError
		);

		// Invalid confidence level
		Vector<Real> sample({10.0, 12.0, 14.0});
		REQUIRE_THROWS_AS(
			ConfidenceIntervalMean(sample, 0.0),  // At boundary
			StatisticsError
		);
		REQUIRE_THROWS_AS(
			ConfidenceIntervalMean(sample, 1.0),  // At boundary
			StatisticsError
		);
		REQUIRE_THROWS_AS(
			ConfidenceIntervalMean(sample, -0.5),  // Negative
			StatisticsError
		);
		REQUIRE_THROWS_AS(
			ConfidenceIntervalMean(sample, 1.5),  // > 1
			StatisticsError
		);

		// Proportion edge cases
		REQUIRE_THROWS_AS(
			ConfidenceIntervalProportion(-1, 100, 0.95),  // Negative successes
			StatisticsError
		);
		REQUIRE_THROWS_AS(
			ConfidenceIntervalProportion(101, 100, 0.95),  // successes > trials
			StatisticsError
		);
		REQUIRE_THROWS_AS(
			ConfidenceIntervalProportion(50, 0, 0.95),  // Zero trials
			StatisticsError
		);

		// Paired differences - unequal sizes
		Vector<Real> v1({10.0, 12.0});
		Vector<Real> v2({15.0, 17.0, 19.0});
		REQUIRE_THROWS_AS(
			ConfidenceIntervalPairedDifference(v1, v2, 0.95),
			StatisticsError
		);
	}

	/******************************************************************************/
	/*****             Detailed API Tests - Confidence Intervals              *****/
	/******************************************************************************/

	TEST_CASE("Statistics - ConfidenceIntervalMeanDetailed values match simple API", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;
		Vector<Real> sample({10.0, 12.0, 14.0, 16.0, 18.0, 20.0});

		auto simple = ConfidenceIntervalMean(sample, 0.95);
		auto detailed = ConfidenceIntervalMeanDetailed(sample, 0.95);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ConfidenceIntervalMean");
		REQUIRE_THAT(detailed.estimate, WithinAbs(simple.estimate, 1e-12));
		REQUIRE_THAT(detailed.lowerBound, WithinAbs(simple.lowerBound, 1e-12));
		REQUIRE_THAT(detailed.upperBound, WithinAbs(simple.upperBound, 1e-12));
		REQUIRE_THAT(detailed.marginOfError, WithinAbs(simple.marginOfError, 1e-12));
		REQUIRE_THAT(detailed.confidenceLevel, WithinAbs(simple.confidenceLevel, 1e-12));
		REQUIRE(detailed.elapsed_time_ms >= 0.0);
	}

	TEST_CASE("Statistics - ConfidenceIntervalMeanDifferenceDetailed values match simple API", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;
		Vector<Real> s1({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> s2({11.0, 13.0, 15.0, 17.0, 19.0});

		auto simple = ConfidenceIntervalMeanDifference(s1, s2, 0.95);
		auto detailed = ConfidenceIntervalMeanDifferenceDetailed(s1, s2, 0.95);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ConfidenceIntervalMeanDifference");
		REQUIRE_THAT(detailed.estimate, WithinAbs(simple.estimate, 1e-12));
		REQUIRE_THAT(detailed.lowerBound, WithinAbs(simple.lowerBound, 1e-12));
		REQUIRE_THAT(detailed.upperBound, WithinAbs(simple.upperBound, 1e-12));
	}

	TEST_CASE("Statistics - ConfidenceIntervalProportionDetailed values match simple API", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;

		auto simple = ConfidenceIntervalProportion(60, 100, 0.95);
		auto detailed = ConfidenceIntervalProportionDetailed(60, 100, 0.95);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ConfidenceIntervalProportion");
		REQUIRE_THAT(detailed.estimate, WithinAbs(simple.estimate, 1e-12));
		REQUIRE_THAT(detailed.lowerBound, WithinAbs(simple.lowerBound, 1e-12));
		REQUIRE_THAT(detailed.upperBound, WithinAbs(simple.upperBound, 1e-12));
	}

	TEST_CASE("Statistics - ConfidenceIntervalProportionDifferenceDetailed values match simple API", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;

		auto simple = ConfidenceIntervalProportionDifference(60, 100, 40, 100, 0.95);
		auto detailed = ConfidenceIntervalProportionDifferenceDetailed(60, 100, 40, 100, 0.95);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ConfidenceIntervalProportionDifference");
		REQUIRE_THAT(detailed.estimate, WithinAbs(simple.estimate, 1e-12));
		REQUIRE_THAT(detailed.lowerBound, WithinAbs(simple.lowerBound, 1e-12));
		REQUIRE_THAT(detailed.upperBound, WithinAbs(simple.upperBound, 1e-12));
	}

	TEST_CASE("Statistics - ConfidenceIntervalPairedDifferenceDetailed values match simple API", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;
		Vector<Real> before({10.0, 12.0, 14.0, 16.0, 18.0});
		Vector<Real> after({11.0, 14.0, 15.0, 18.0, 20.0});

		auto simple = ConfidenceIntervalPairedDifference(before, after, 0.95);
		auto detailed = ConfidenceIntervalPairedDifferenceDetailed(before, after, 0.95);

		REQUIRE(detailed.IsSuccess());
		REQUIRE(detailed.algorithm_name == "ConfidenceIntervalPairedDifference");
		REQUIRE_THAT(detailed.estimate, WithinAbs(simple.estimate, 1e-12));
		REQUIRE_THAT(detailed.lowerBound, WithinAbs(simple.lowerBound, 1e-12));
		REQUIRE_THAT(detailed.upperBound, WithinAbs(simple.upperBound, 1e-12));
	}

	TEST_CASE("Statistics - ConfidenceIntervalMeanDetailed ConvertToStatus on invalid input", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;
		Vector<Real> tooSmall({10.0});
		StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = ConfidenceIntervalMeanDetailed(tooSmall, 0.95, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(result.error_message.empty());
	}

	TEST_CASE("Statistics - ConfidenceIntervalMeanDetailed Propagate policy throws", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;
		Vector<Real> tooSmall({10.0});
		StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::Propagate;

		REQUIRE_THROWS(ConfidenceIntervalMeanDetailed(tooSmall, 0.95, config));
	}

	TEST_CASE("Statistics - ConfidenceIntervalProportionDetailed ConvertToStatus on invalid input", "[statistics][confidence][Detailed]")
	{
		using namespace MML::Statistics;
		StatisticsConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		auto result = ConfidenceIntervalProportionDetailed(-1, 100, 0.95, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
	}
}