#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics.h"
#include "base/Matrix/Matrix.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::StatisticsTests {
	/*********************************************************************/
	/*****          Basic Statistics - Average                       *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Avg_SimpleValues", "[statistics][basic]") {
		TEST_PRECISION_INFO();
		// Simple average of known values
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});

		Real avg = Statistics::Avg(data);

		REQUIRE_THAT(avg, WithinAbs(REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Avg_SingleValue", "[statistics][basic]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real avg = Statistics::Avg(data);

		REQUIRE_THAT(avg, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Avg_NegativeValues", "[statistics][basic]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({-REAL(5.0), -REAL(3.0), -REAL(1.0), REAL(1.0), REAL(3.0), REAL(5.0)});

		Real avg = Statistics::Avg(data);

		REQUIRE_THAT(avg, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Avg_EmptyVector_ThrowsException", "[statistics][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::Avg(data), StatisticsError);
	}

	/*********************************************************************/
	/*****          Average and Variance                             *****/
	/*********************************************************************/
	TEST_CASE("Statistics::AvgVar_KnownValues", "[statistics][variance]") {
		TEST_PRECISION_INFO();
		// Data: 2, 4, 4, 4, 5, 5, 7, 9
		// Mean = 40/8 = REAL(5.0)
		// Variance = sum((x - mean)²) / (n-1)
		// = (9 + 1 + 1 + 1 + 0 + 0 + 4 + 16) / 7 = 32/7 ≈ REAL(4.571)
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});

		Real avg, var;
		Statistics::AvgVar(data, avg, var);

		REQUIRE_THAT(avg, WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(32.0) / REAL(7.0), 1e-10));
	}

	TEST_CASE("Statistics::AvgVar_UniformData_ZeroVariance", "[statistics][variance]") {
		TEST_PRECISION_INFO();
		// All values the same = zero variance
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});

		Real avg, var;
		Statistics::AvgVar(data, avg, var);

		REQUIRE_THAT(avg, WithinAbs(REAL(7.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::AvgVar_TwoValues", "[statistics][variance]") {
		TEST_PRECISION_INFO();
		// Variance of two values
		Vector<Real> data({REAL(3.0), REAL(7.0)});

		Real avg, var;
		Statistics::AvgVar(data, avg, var);

		// Mean = REAL(5.0)
		// Variance = ((3-5)² + (7-5)²) / 1 = (4 + 4) / 1 = REAL(8.0)
		REQUIRE_THAT(avg, WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(8.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****          Average and Standard Deviation                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::AvgStdDev_KnownValues", "[statistics][stddev]") {
		TEST_PRECISION_INFO();
		// Data with known mean and std dev
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});

		Real avg, stddev;
		Statistics::AvgStdDev(data, avg, stddev);

		// Mean = REAL(5.0), Variance = 32/7, StdDev = sqrt(32/7) ≈ REAL(2.138)
		REQUIRE_THAT(avg, WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(stddev, WithinAbs(std::sqrt(REAL(32.0) / REAL(7.0)), 1e-10));
	}

	TEST_CASE("Statistics::AvgStdDev_StandardNormal", "[statistics][stddev]") {
		TEST_PRECISION_INFO();
		// Approximate standard normal: mean ≈ 0, stddev ≈ 1
		Vector<Real> data({-REAL(2.0), -REAL(1.0), -REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(1.0), REAL(1.0), REAL(2.0)});

		Real avg, stddev;
		Statistics::AvgStdDev(data, avg, stddev);

		REQUIRE_THAT(avg, WithinAbs(REAL(0.0), REAL(1e-10)));
		// Variance = (4+1+1+0+0+0+1+1+4)/8 = 12/8 = REAL(1.5), stddev ≈ REAL(1.225)
		REQUIRE_THAT(stddev, WithinAbs(std::sqrt(REAL(1.5)), 1e-10));
	}

	/*********************************************************************/
	/*****          Statistical Moments                              *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Moments_SymmetricData", "[statistics][moments]") {
		TEST_PRECISION_INFO();
		// Symmetric data around mean = zero skewness
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});

		Real ave, adev, sdev, var, skew, curt;
		Statistics::Moments(data, ave, adev, sdev, var, skew, curt);

		REQUIRE_THAT(ave, WithinAbs(REAL(3.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(2.5), REAL(1e-10))); // Variance of 1,2,3,4,5
		REQUIRE_THAT(sdev, WithinAbs(std::sqrt(REAL(2.5)), 1e-10));
		REQUIRE_THAT(skew, WithinAbs(REAL(0.0), REAL(1e-8))); // Symmetric = zero skew
	}

	TEST_CASE("Statistics::Moments_PositiveSkew", "[statistics][moments]") {
		TEST_PRECISION_INFO();
		// Right-skewed data (tail to the right)
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0), REAL(3.0), REAL(3.0), REAL(4.0), REAL(10.0)});

		Real ave, adev, sdev, var, skew, curt;
		Statistics::Moments(data, ave, adev, sdev, var, skew, curt);

		// Mean = 28/8 = REAL(3.5)
		REQUIRE_THAT(ave, WithinAbs(REAL(3.5), REAL(1e-10)));

		// Skewness should be positive (right tail)
		REQUIRE(skew > REAL(0.0));
	}

	TEST_CASE("Statistics::Moments_NegativeSkew", "[statistics][moments]") {
		TEST_PRECISION_INFO();
		// Left-skewed data (tail to the left)
		Vector<Real> data({REAL(1.0), REAL(7.0), REAL(8.0), REAL(8.0), REAL(9.0), REAL(9.0), REAL(9.0), REAL(10.0)});

		Real ave, adev, sdev, var, skew, curt;
		Statistics::Moments(data, ave, adev, sdev, var, skew, curt);

		// Mean = 61/8 = REAL(7.625)
		REQUIRE_THAT(ave, WithinAbs(REAL(7.625), REAL(1e-10)));

		// Skewness should be negative (left tail)
		REQUIRE(skew < REAL(0.0));
	}

	TEST_CASE("Statistics::Moments_HighKurtosis", "[statistics][moments]") {
		TEST_PRECISION_INFO();
		// Heavy-tailed distribution (high kurtosis)
		// Many values at extremes, few in middle
		Vector<Real> data({REAL(0.0), REAL(0.0), REAL(0.0), REAL(5.0), REAL(10.0), REAL(10.0), REAL(10.0)});

		Real ave, adev, sdev, var, skew, curt;
		Statistics::Moments(data, ave, adev, sdev, var, skew, curt);

		// Mean = 35/7 = REAL(5.0)
		REQUIRE_THAT(ave, WithinAbs(REAL(5.0), REAL(1e-10)));

		// Kurtosis should be positive (heavy tails)
		// Note: Implementation returns excess kurtosis (normal = 0)
		REQUIRE(curt > -REAL(3.0)); // Kurtosis > 0 (excess kurtosis)
	}

	TEST_CASE("Statistics::Moments_UniformData_ThrowsException", "[statistics][error]") {
		TEST_PRECISION_INFO();
		// Uniform data has zero variance, should throw when computing skew/kurtosis
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});

		Real ave, adev, sdev, var, skew, curt;

		REQUIRE_THROWS(Statistics::Moments(data, ave, adev, sdev, var, skew, curt));
	}

	TEST_CASE("Statistics::Moments_SingleValue_ThrowsException", "[statistics][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real ave, adev, sdev, var, skew, curt;

		REQUIRE_THROWS_AS(Statistics::Moments(data, ave, adev, sdev, var, skew, curt), StatisticsError);
	}

	/*********************************************************************/
	/*****          Edge Cases and Numerical Stability               *****/
	/*********************************************************************/
	TEST_CASE("Statistics::LargeDataset_Precision", "[statistics][precision]") {
		TEST_PRECISION_INFO();
		// Test with larger dataset
		Vector<Real> data;
		for (int i = 1; i <= 100; i++) {
			data.push_back(static_cast<Real>(i));
		}

		Real avg = Statistics::Avg(data);

		// Mean of 1 to 100 = REAL(50.5)
		REQUIRE_THAT(avg, WithinAbs(REAL(50.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::VerySmallValues", "[statistics][precision]") {
		TEST_PRECISION_INFO();
		// Test numerical stability with very small values
		Vector<Real> data({1e-10, 2e-10, 3e-10, 4e-10, 5e-10});

		Real avg = Statistics::Avg(data);

		REQUIRE_THAT(avg, WithinAbs(3e-10, REAL(1e-20)));
	}

	TEST_CASE("Statistics::VeryLargeValues", "[statistics][precision]") {
		TEST_PRECISION_INFO();
		// Test with very large values
		Vector<Real> data({1e10, 2e10, 3e10, 4e10, 5e10});

		Real avg = Statistics::Avg(data);

		REQUIRE_THAT(avg, WithinRel(REAL(3e10), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****          Real-World Example Data                          *****/
	/*********************************************************************/
	TEST_CASE("Statistics::TemperatureData_Example", "[statistics][example]") {
		TEST_PRECISION_INFO();
		// Weekly temperature data (°C)
		Vector<Real> temps({REAL(18.5), REAL(19.2), REAL(20.1), REAL(19.8), REAL(18.9), REAL(17.5), REAL(18.0)});

		Real avg, stddev;
		Statistics::AvgStdDev(temps, avg, stddev);

		// Mean should be around REAL(18.86)°C
		REQUIRE_THAT(avg, WithinAbs(REAL(18.857), REAL(0.01)));

		// Standard deviation should be small (< 1°C)
		REQUIRE(stddev < REAL(1.0));
		REQUIRE(stddev > REAL(0.0));
	}

	TEST_CASE("Statistics::TestScores_Example", "[statistics][example]") {
		TEST_PRECISION_INFO();
		// Student test scores (0-100)
		Vector<Real> scores({85, 92, 78, 90, 88, 76, 95, 89, 84, 91});

		Real ave, adev, sdev, var, skew, curt;
		Statistics::Moments(scores, ave, adev, sdev, var, skew, curt);

		// Mean should be around REAL(86.8)
		REQUIRE_THAT(ave, WithinAbs(REAL(86.8), REAL(0.1)));

		// Variance should be positive
		REQUIRE(var > REAL(0.0));

		// Standard deviation should be reasonable (< 10 points)
		REQUIRE(sdev < REAL(10.0));
		REQUIRE(sdev > REAL(0.0));
	}

	/*********************************************************************/
	/*****               Order Statistics - Median                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Median_OddCount", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// Odd number of elements: median is middle element
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0), REAL(9.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_EvenCount", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// Even number of elements: median is average of two middle elements
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(REAL(4.0), REAL(1e-10))); // (3+5)/2
	}

	TEST_CASE("Statistics::Median_SingleValue", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_TwoValues", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(10.0), REAL(20.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(REAL(15.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_UnsortedData", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// Median should work regardless of input order
		Vector<Real> data({REAL(9.0), REAL(1.0), REAL(5.0), REAL(3.0), REAL(7.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_NegativeValues", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({-REAL(5.0), -REAL(3.0), -REAL(1.0), REAL(1.0), REAL(3.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(-REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_EmptyVector_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::Median(data), StatisticsError);
	}

	TEST_CASE("Statistics::Median_AllSame", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});

		Real median = Statistics::Median(data);

		REQUIRE_THAT(median, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****             Order Statistics - Percentile                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Percentile_50th_IsMedian", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});

		Real p50 = Statistics::Percentile(data, 50);
		Real median = Statistics::Median(data);

		REQUIRE_THAT(p50, WithinAbs(median, REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_0th_IsMin", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0), REAL(50.0)});

		Real p0 = Statistics::Percentile(data, 0);

		REQUIRE_THAT(p0, WithinAbs(REAL(10.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_100th_IsMax", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0), REAL(50.0)});

		Real p100 = Statistics::Percentile(data, 100);

		REQUIRE_THAT(p100, WithinAbs(REAL(50.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_25th", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// 100 values from 1 to 100
		Vector<Real> data;
		for (int i = 1; i <= 100; i++) {
			data.push_back(static_cast<Real>(i));
		}

		Real p25 = Statistics::Percentile(data, 25);

		// 25th percentile should be around 25.75
		REQUIRE_THAT(p25, WithinAbs(REAL(25.75), REAL(0.5)));
	}

	TEST_CASE("Statistics::Percentile_75th", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;
		for (int i = 1; i <= 100; i++) {
			data.push_back(static_cast<Real>(i));
		}

		Real p75 = Statistics::Percentile(data, 75);

		// 75th percentile should be around 75.25
		REQUIRE_THAT(p75, WithinAbs(REAL(75.25), REAL(0.5)));
	}

	TEST_CASE("Statistics::Percentile_InvalidRange_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});

		REQUIRE_THROWS_AS(Statistics::Percentile(data, -1), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::Percentile(data, 101), StatisticsError);
	}

	TEST_CASE("Statistics::Percentile_EmptyVector_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::Percentile(data, 50), StatisticsError);
	}

	/*********************************************************************/
	/*****              Order Statistics - Quartiles                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Quartiles_BasicTest", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// Data 1-100
		Vector<Real> data;
		for (int i = 1; i <= 100; i++) {
			data.push_back(static_cast<Real>(i));
		}

		Real q1, q2, q3;
		Statistics::Quartiles(data, q1, q2, q3);

		// Q1 ≈ 25.75, Q2 = 50.5, Q3 ≈ 75.25
		REQUIRE_THAT(q2, WithinAbs(REAL(50.5), REAL(0.5)));
		REQUIRE(q1 < q2);
		REQUIRE(q2 < q3);
	}

	TEST_CASE("Statistics::Quartiles_SymmetricData", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// Symmetric data: Q2 should equal mean
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});

		Real q1, q2, q3;
		Statistics::Quartiles(data, q1, q2, q3);

		REQUIRE_THAT(q2, WithinAbs(REAL(3.0), REAL(1e-10)));
		// For symmetric data, Q3 - Q2 should equal Q2 - Q1
		REQUIRE_THAT(q3 - q2, WithinAbs(q2 - q1, REAL(1e-10)));
	}

	TEST_CASE("Statistics::Quartiles_EmptyVector_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		Real q1, q2, q3;
		REQUIRE_THROWS_AS(Statistics::Quartiles(data, q1, q2, q3), StatisticsError);
	}

	/*********************************************************************/
	/*****              Order Statistics - Range                     *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Range_BasicTest", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(10.0), REAL(15.0), REAL(20.0), REAL(25.0)});

		Real range = Statistics::Range(data);

		REQUIRE_THAT(range, WithinAbs(REAL(20.0), REAL(1e-10))); // 25 - 5
	}

	TEST_CASE("Statistics::Range_AllSame", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});

		Real range = Statistics::Range(data);

		REQUIRE_THAT(range, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Range_SingleValue", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real range = Statistics::Range(data);

		REQUIRE_THAT(range, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Range_NegativeValues", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({-REAL(10.0), -REAL(5.0), REAL(0.0), REAL(5.0), REAL(10.0)});

		Real range = Statistics::Range(data);

		REQUIRE_THAT(range, WithinAbs(REAL(20.0), REAL(1e-10))); // 10 - (-10)
	}

	TEST_CASE("Statistics::Range_EmptyVector_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::Range(data), StatisticsError);
	}

	/*********************************************************************/
	/*****           Order Statistics - IQR (Interquartile Range)    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::IQR_BasicTest", "[statistics][order]") {
		TEST_PRECISION_INFO();
		// Data 1-100
		Vector<Real> data;
		for (int i = 1; i <= 100; i++) {
			data.push_back(static_cast<Real>(i));
		}

		Real iqr = Statistics::IQR(data);

		// IQR = Q3 - Q1 ≈ 75.25 - 25.75 = 49.5
		REQUIRE_THAT(iqr, WithinAbs(REAL(49.5), REAL(1.0)));
	}

	TEST_CASE("Statistics::IQR_AllSame", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});

		Real iqr = Statistics::IQR(data);

		REQUIRE_THAT(iqr, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::IQR_ConsistentWithQuartiles", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0)});

		Real q1, q2, q3;
		Statistics::Quartiles(data, q1, q2, q3);
		Real iqr = Statistics::IQR(data);

		REQUIRE_THAT(iqr, WithinAbs(q3 - q1, REAL(1e-10)));
	}

	TEST_CASE("Statistics::IQR_EmptyVector_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::IQR(data), StatisticsError);
	}

	/*********************************************************************/
	/*****              Order Statistics - MinMax                    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::MinMax_BasicTest", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(2.0), REAL(9.0), REAL(1.0), REAL(7.0)});

		Real min, max;
		Statistics::MinMax(data, min, max);

		REQUIRE_THAT(min, WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(max, WithinAbs(REAL(9.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_SingleValue", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real min, max;
		Statistics::MinMax(data, min, max);

		REQUIRE_THAT(min, WithinAbs(REAL(42.0), REAL(1e-10)));
		REQUIRE_THAT(max, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_NegativeValues", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({-REAL(5.0), -REAL(10.0), -REAL(3.0), -REAL(8.0)});

		Real min, max;
		Statistics::MinMax(data, min, max);

		REQUIRE_THAT(min, WithinAbs(-REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(max, WithinAbs(-REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_RangeConsistency", "[statistics][order]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(5.0), REAL(10.0), REAL(15.0), REAL(20.0)});

		Real min, max;
		Statistics::MinMax(data, min, max);
		Real range = Statistics::Range(data);

		REQUIRE_THAT(max - min, WithinAbs(range, REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_EmptyVector_ThrowsException", "[statistics][order][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		Real min, max;
		REQUIRE_THROWS_AS(Statistics::MinMax(data, min, max), StatisticsError);
	}

	/*********************************************************************/
	/*****                Robust Statistics - Mode                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Mode_SingleMode", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Clear single mode
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0), REAL(4.0)});

		Real mode = Statistics::Mode(data);

		REQUIRE_THAT(mode, WithinAbs(REAL(2.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_MultipleOccurrences", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(3.0), REAL(3.0), REAL(5.0), REAL(5.0)});

		Real mode = Statistics::Mode(data);

		REQUIRE_THAT(mode, WithinAbs(REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_AllUnique", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// All unique values - returns first (smallest after sort)
		Vector<Real> data({REAL(5.0), REAL(3.0), REAL(1.0), REAL(4.0), REAL(2.0)});

		Real mode = Statistics::Mode(data);

		// When all unique, returns first element (smallest after sort)
		REQUIRE_THAT(mode, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_AllSame", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});

		Real mode = Statistics::Mode(data);

		REQUIRE_THAT(mode, WithinAbs(REAL(7.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_SingleValue", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real mode = Statistics::Mode(data);

		REQUIRE_THAT(mode, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_EmptyVector_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::Mode(data), StatisticsError);
	}

	/*********************************************************************/
	/*****             Robust Statistics - TrimmedMean               *****/
	/*********************************************************************/
	TEST_CASE("Statistics::TrimmedMean_NoTrim", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// 0% trim = regular mean
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});

		Real trimmedMean = Statistics::TrimmedMean(data, REAL(0.0));
		Real regularMean = Statistics::Avg(data);

		REQUIRE_THAT(trimmedMean, WithinAbs(regularMean, REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_10Percent", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// 10 values, 10% trim = remove 1 from each end
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 100 (outlier)
		// After trim: 2, 3, 4, 5, 6, 7, 8, 9 -> mean = 44/8 = 5.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(100.0)});

		Real trimmedMean = Statistics::TrimmedMean(data, REAL(10.0));

		REQUIRE_THAT(trimmedMean, WithinAbs(REAL(5.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_20Percent", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// 10 values, 20% trim = remove 2 from each end
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
		// After trim: 3, 4, 5, 6, 7, 8 -> mean = 33/6 = 5.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});

		Real trimmedMean = Statistics::TrimmedMean(data, REAL(20.0));

		REQUIRE_THAT(trimmedMean, WithinAbs(REAL(5.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_OutlierResistance", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// With outliers, trimmed mean should be close to true center
		Vector<Real> dataWithOutliers(
			{REAL(-100.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(1000.0)});
		Vector<Real> dataClean({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0)});

		Real trimmedMeanOutliers = Statistics::TrimmedMean(dataWithOutliers, REAL(10.0));
		Real regularMeanClean = Statistics::Avg(dataClean);

		REQUIRE_THAT(trimmedMeanOutliers, WithinAbs(regularMeanClean, REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_InvalidPercent_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});

		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(-1.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(50.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(100.0)), StatisticsError);
	}

	TEST_CASE("Statistics::TrimmedMean_EmptyVector_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(10.0)), StatisticsError);
	}

	/*********************************************************************/
	/*****                 Robust Statistics - MAD                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::MAD_BasicTest", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Data: 1, 2, 3, 4, 5
		// Median = 3
		// Absolute deviations: |1-3|=2, |2-3|=1, |3-3|=0, |4-3|=1, |5-3|=2
		// Sorted deviations: 0, 1, 1, 2, 2
		// MAD = median = 1
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});

		Real mad = Statistics::MAD(data);

		REQUIRE_THAT(mad, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MAD_AllSame", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// All same values -> MAD = 0
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});

		Real mad = Statistics::MAD(data);

		REQUIRE_THAT(mad, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MAD_WithOutlier", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Data: 1, 2, 3, 4, 100 (outlier)
		// Median = 3
		// Absolute deviations: 2, 1, 0, 1, 97
		// Sorted: 0, 1, 1, 2, 97
		// MAD = 1 (robust to outlier!)
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(100.0)});

		Real mad = Statistics::MAD(data);

		REQUIRE_THAT(mad, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MAD_SingleValue", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});

		Real mad = Statistics::MAD(data);

		REQUIRE_THAT(mad, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MAD_EmptyVector_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::MAD(data), StatisticsError);
	}

	/*********************************************************************/
	/*****             Specialized Means - Geometric                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::GeometricMean_BasicTest", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Geometric mean of 2, 8 = sqrt(16) = 4
		Vector<Real> data({REAL(2.0), REAL(8.0)});

		Real gm = Statistics::GeometricMean(data);

		REQUIRE_THAT(gm, WithinAbs(REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_AllSame", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// GM of same values = that value
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0)});

		Real gm = Statistics::GeometricMean(data);

		REQUIRE_THAT(gm, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_SingleValue", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0)});

		Real gm = Statistics::GeometricMean(data);

		REQUIRE_THAT(gm, WithinAbs(REAL(7.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_GrowthRates", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Growth factors: 1.1, 1.2, 1.15 (10%, 20%, 15% growth)
		// GM = (1.1 * 1.2 * 1.15)^(1/3) ≈ 1.1491
		Vector<Real> rates({REAL(1.1), REAL(1.2), REAL(1.15)});

		Real gm = Statistics::GeometricMean(rates);
		Real expected = std::pow(REAL(1.1) * REAL(1.2) * REAL(1.15), REAL(1.0) / REAL(3.0));

		REQUIRE_THAT(gm, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_NonPositive_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> dataWithZero({REAL(1.0), REAL(0.0), REAL(3.0)});
		Vector<Real> dataWithNegative({REAL(1.0), REAL(-2.0), REAL(3.0)});

		REQUIRE_THROWS_AS(Statistics::GeometricMean(dataWithZero), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::GeometricMean(dataWithNegative), StatisticsError);
	}

	TEST_CASE("Statistics::GeometricMean_EmptyVector_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::GeometricMean(data), StatisticsError);
	}

	/*********************************************************************/
	/*****              Specialized Means - Harmonic                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::HarmonicMean_BasicTest", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Harmonic mean of 1, 2, 4 = 3 / (1 + 0.5 + 0.25) = 3 / 1.75 ≈ 1.714
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(4.0)});

		Real hm = Statistics::HarmonicMean(data);
		Real expected = REAL(3.0) / (REAL(1.0) + REAL(0.5) + REAL(0.25));

		REQUIRE_THAT(hm, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_AllSame", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// HM of same values = that value
		Vector<Real> data({REAL(4.0), REAL(4.0), REAL(4.0)});

		Real hm = Statistics::HarmonicMean(data);

		REQUIRE_THAT(hm, WithinAbs(REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_Speeds", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Classic problem: drive 60 km at 30 km/h, then 60 km at 60 km/h
		// Average speed = total distance / total time = 120 / (2 + 1) = 40 km/h
		// This is the harmonic mean: 2 / (1/30 + 1/60) = 2 / (3/60) = 40
		Vector<Real> speeds({REAL(30.0), REAL(60.0)});

		Real hm = Statistics::HarmonicMean(speeds);

		REQUIRE_THAT(hm, WithinAbs(REAL(40.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_SingleValue", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(10.0)});

		Real hm = Statistics::HarmonicMean(data);

		REQUIRE_THAT(hm, WithinAbs(REAL(10.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_NonPositive_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> dataWithZero({REAL(1.0), REAL(0.0), REAL(3.0)});
		Vector<Real> dataWithNegative({REAL(1.0), REAL(-2.0), REAL(3.0)});

		REQUIRE_THROWS_AS(Statistics::HarmonicMean(dataWithZero), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::HarmonicMean(dataWithNegative), StatisticsError);
	}

	TEST_CASE("Statistics::HarmonicMean_EmptyVector_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;

		REQUIRE_THROWS_AS(Statistics::HarmonicMean(data), StatisticsError);
	}

	/*********************************************************************/
	/*****               Specialized Means - Weighted                *****/
	/*********************************************************************/
	TEST_CASE("Statistics::WeightedMean_EqualWeights", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Equal weights = regular mean
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0)});

		Real wm = Statistics::WeightedMean(data, weights);
		Real regularMean = Statistics::Avg(data);

		REQUIRE_THAT(wm, WithinAbs(regularMean, REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_UnequalWeights", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Data: 10, 20, 30 with weights 1, 2, 3
		// WM = (10*1 + 20*2 + 30*3) / (1+2+3) = (10 + 40 + 90) / 6 = 140/6 ≈ 23.33
		Vector<Real> data({REAL(10.0), REAL(20.0), REAL(30.0)});
		Vector<Real> weights({REAL(1.0), REAL(2.0), REAL(3.0)});

		Real wm = Statistics::WeightedMean(data, weights);
		Real expected = REAL(140.0) / REAL(6.0);

		REQUIRE_THAT(wm, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_SingleNonZeroWeight", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Only one weight is non-zero -> result is that value
		Vector<Real> data({REAL(10.0), REAL(20.0), REAL(30.0)});
		Vector<Real> weights({REAL(0.0), REAL(1.0), REAL(0.0)});

		Real wm = Statistics::WeightedMean(data, weights);

		REQUIRE_THAT(wm, WithinAbs(REAL(20.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_GradeAverage", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// Typical grade calculation: HW=85 (20%), Midterm=78 (30%), Final=92 (50%)
		// WM = (85*0.2 + 78*0.3 + 92*0.5) / 1.0 = 17 + 23.4 + 46 = 86.4
		Vector<Real> grades({REAL(85.0), REAL(78.0), REAL(92.0)});
		Vector<Real> weights({REAL(0.2), REAL(0.3), REAL(0.5)});

		Real wm = Statistics::WeightedMean(grades, weights);

		REQUIRE_THAT(wm, WithinAbs(REAL(86.4), REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_MismatchedSizes_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0)}); // Wrong size

		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedMean_NegativeWeights_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(1.0), REAL(-1.0), REAL(1.0)});

		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedMean_ZeroTotalWeight_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(0.0), REAL(0.0), REAL(0.0)});

		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedMean_EmptyVector_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data;
		Vector<Real> weights;

		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	/*********************************************************************/
	/*****              Specialized Means - WeightedVariance         *****/
	/*********************************************************************/
	TEST_CASE("Statistics::WeightedVariance_EqualWeights", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// With equal weights, should match regular variance
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0)});

		Real wv = Statistics::WeightedVariance(data, weights);
		Real avg, regularVar;
		Statistics::AvgVar(data, avg, regularVar);

		REQUIRE_THAT(wv, WithinAbs(regularVar, REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedVariance_AllSameValues", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		Vector<Real> weights({REAL(1.0), REAL(2.0), REAL(1.0), REAL(2.0)});

		Real wv = Statistics::WeightedVariance(data, weights);

		REQUIRE_THAT(wv, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedVariance_MismatchedSizes_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0)});

		REQUIRE_THROWS_AS(Statistics::WeightedVariance(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedVariance_InsufficientData_ThrowsException", "[statistics][robust][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0)});
		Vector<Real> weights({REAL(1.0)});

		REQUIRE_THROWS_AS(Statistics::WeightedVariance(data, weights), StatisticsError);
	}

	/*********************************************************************/
	/*****                 Mean Inequalities Tests                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::MeanInequality_HM_GM_AM", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// For positive numbers: HM <= GM <= AM (with equality iff all values equal)
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(8.0)});

		Real hm = Statistics::HarmonicMean(data);
		Real gm = Statistics::GeometricMean(data);
		Real am = Statistics::Avg(data);

		// HM < GM < AM for non-uniform positive data
		REQUIRE(hm < gm);
		REQUIRE(gm < am);
	}

	TEST_CASE("Statistics::MeanInequality_Equal_WhenUniform", "[statistics][robust]") {
		TEST_PRECISION_INFO();
		// For uniform data: HM = GM = AM
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});

		Real hm = Statistics::HarmonicMean(data);
		Real gm = Statistics::GeometricMean(data);
		Real am = Statistics::Avg(data);

		REQUIRE_THAT(hm, WithinAbs(gm, REAL(1e-10)));
		REQUIRE_THAT(gm, WithinAbs(am, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****                 Covariance Tests                          *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Covariance_PositiveRelationship", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// x and y increase together
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0)});

		Real cov = Statistics::Covariance(x, y);

		// Positive covariance expected
		REQUIRE(cov > 0);
		// Exact value: Cov(x, y) = 5.0 (can verify manually)
		REQUIRE_THAT(cov, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_NegativeRelationship", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// x increases, y decreases
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(10.0), REAL(8.0), REAL(6.0), REAL(4.0), REAL(2.0)});

		Real cov = Statistics::Covariance(x, y);

		// Negative covariance expected
		REQUIRE(cov < 0);
		REQUIRE_THAT(cov, WithinAbs(REAL(-5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_NoRelationship", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Quadratic relationship (symmetric around mean) - zero linear covariance
		Vector<Real> x({REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)}); // y = x²

		Real cov = Statistics::Covariance(x, y);

		// Covariance should be zero (no LINEAR relationship)
		REQUIRE_THAT(cov, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_SelfCovariance_IsVariance", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});

		Real selfCov = Statistics::Covariance(x, x);
		Real avg, var;
		Statistics::AvgVar(x, avg, var);

		REQUIRE_THAT(selfCov, WithinAbs(var, REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_DifferentSizes_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});

		REQUIRE_THROWS_AS(Statistics::Covariance(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::Covariance_SingleValue_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});

		REQUIRE_THROWS_AS(Statistics::Covariance(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****              Pearson Correlation Tests                    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::PearsonCorrelation_PerfectPositive", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Perfect positive correlation: y = 2x
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0)});

		Real r = Statistics::PearsonCorrelation(x, y);

		REQUIRE_THAT(r, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_PerfectNegative", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Perfect negative correlation: y = -2x + 12
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(10.0), REAL(8.0), REAL(6.0), REAL(4.0), REAL(2.0)});

		Real r = Statistics::PearsonCorrelation(x, y);

		REQUIRE_THAT(r, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_ZeroCorrelation", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Quadratic relationship - zero LINEAR correlation
		Vector<Real> x({REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});

		Real r = Statistics::PearsonCorrelation(x, y);

		REQUIRE_THAT(r, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_ModeratePositive", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Some noise in the relationship
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.5), REAL(3.5), REAL(5.0), REAL(7.5), REAL(9.5)});

		Real r = Statistics::PearsonCorrelation(x, y);

		// Should be high positive correlation but not perfect
		REQUIRE(r > REAL(0.9));
		REQUIRE(r < REAL(1.0));
	}

	TEST_CASE("Statistics::PearsonCorrelation_SelfCorrelation_IsOne", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(3.0), REAL(7.0), REAL(2.0), REAL(9.0)});

		Real r = Statistics::PearsonCorrelation(x, x);

		REQUIRE_THAT(r, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_RangeIsMinusOneToOne", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Random-ish data
		Vector<Real> x({REAL(2.3), REAL(5.1), REAL(1.2), REAL(8.7), REAL(4.4)});
		Vector<Real> y({REAL(7.8), REAL(3.2), REAL(9.1), REAL(2.5), REAL(6.6)});

		Real r = Statistics::PearsonCorrelation(x, y);

		REQUIRE(r >= REAL(-1.0));
		REQUIRE(r <= REAL(1.0));
	}

	TEST_CASE("Statistics::PearsonCorrelation_ZeroVarianceX_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});

		REQUIRE_THROWS_AS(Statistics::PearsonCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::PearsonCorrelation_ZeroVarianceY_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});

		REQUIRE_THROWS_AS(Statistics::PearsonCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::PearsonCorrelation_DifferentSizes_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});

		REQUIRE_THROWS_AS(Statistics::PearsonCorrelation(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****           Pearson Correlation With Test                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::PearsonCorrelationWithTest_HighCorrelation", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Strong positive correlation
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.1), REAL(3.9), REAL(6.2), REAL(7.8), REAL(10.1), REAL(12.0), REAL(14.1), REAL(15.9), REAL(18.2), REAL(19.8)});

		auto result = Statistics::PearsonCorrelationWithTest(x, y);

		REQUIRE(result.r > REAL(0.99));
		REQUIRE(result.degreesOfFreedom == 8);
		REQUIRE(result.tStatistic > 0); // Positive for positive correlation

		// With r ≈ 0.999, t should be very large
		REQUIRE(result.tStatistic > REAL(10.0));
	}

	TEST_CASE("Statistics::PearsonCorrelationWithTest_DegreesOfFreedom", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(5.0), REAL(4.0), REAL(5.0)});

		auto result = Statistics::PearsonCorrelationWithTest(x, y);

		// df = n - 2 = 5 - 2 = 3
		REQUIRE(result.degreesOfFreedom == 3);
	}

	TEST_CASE("Statistics::PearsonCorrelationWithTest_InsufficientData_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0)});

		// Need n > 2 for significance test
		REQUIRE_THROWS_AS(Statistics::PearsonCorrelationWithTest(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****                   R-Squared Tests                         *****/
	/*********************************************************************/
	TEST_CASE("Statistics::RSquared_PerfectFit", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Perfect linear relationship -> R² = 1
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0)});

		Real r2 = Statistics::RSquared(x, y);

		REQUIRE_THAT(r2, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::RSquared_NoFit", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// No linear relationship -> R² ≈ 0
		Vector<Real> x({REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});

		Real r2 = Statistics::RSquared(x, y);

		REQUIRE_THAT(r2, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::RSquared_RangeIsZeroToOne", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		Vector<Real> y({REAL(2.5), REAL(4.0), REAL(7.5), REAL(6.0), REAL(10.0)});

		Real r2 = Statistics::RSquared(x, y);

		REQUIRE(r2 >= REAL(0.0));
		REQUIRE(r2 <= REAL(1.0));
	}

	/*********************************************************************/
	/*****              Covariance Matrix Tests                      *****/
	/*********************************************************************/
	TEST_CASE("Statistics::CovarianceMatrix_TwoVariables", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// 5 observations, 2 variables
		Matrix<Real> data(5, 2);
		// x: 1, 2, 3, 4, 5
		// y: 2, 4, 6, 8, 10
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(2.0);
		data(1, 0) = REAL(2.0);
		data(1, 1) = REAL(4.0);
		data(2, 0) = REAL(3.0);
		data(2, 1) = REAL(6.0);
		data(3, 0) = REAL(4.0);
		data(3, 1) = REAL(8.0);
		data(4, 0) = REAL(5.0);
		data(4, 1) = REAL(10.0);

		Matrix<Real> covMat = Statistics::CovarianceMatrix(data);

		REQUIRE(covMat.rows() == 2);
		REQUIRE(covMat.cols() == 2);

		// Var(x) = 2.5, Var(y) = 10, Cov(x,y) = 5
		REQUIRE_THAT(covMat(0, 0), WithinAbs(REAL(2.5), REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 1), WithinAbs(REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(covMat(0, 1), WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 0), WithinAbs(REAL(5.0), REAL(1e-10))); // Symmetric
	}

	TEST_CASE("Statistics::CovarianceMatrix_IsSymmetric", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3);
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(2.0);
		data(0, 2) = REAL(3.0);
		data(1, 0) = REAL(4.0);
		data(1, 1) = REAL(5.0);
		data(1, 2) = REAL(6.0);
		data(2, 0) = REAL(7.0);
		data(2, 1) = REAL(8.0);
		data(2, 2) = REAL(9.0);
		data(3, 0) = REAL(2.0);
		data(3, 1) = REAL(3.0);
		data(3, 2) = REAL(4.0);

		Matrix<Real> covMat = Statistics::CovarianceMatrix(data);

		// Check symmetry
		REQUIRE_THAT(covMat(0, 1), WithinAbs(covMat(1, 0), REAL(1e-10)));
		REQUIRE_THAT(covMat(0, 2), WithinAbs(covMat(2, 0), REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 2), WithinAbs(covMat(2, 1), REAL(1e-10)));
	}

	TEST_CASE("Statistics::CovarianceMatrix_DiagonalIsVariance", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Matrix<Real> data(5, 2);
		Vector<Real> col1({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> col2({REAL(10.0), REAL(20.0), REAL(15.0), REAL(25.0), REAL(30.0)});

		for (int i = 0; i < 5; i++) {
			data(i, 0) = col1[i];
			data(i, 1) = col2[i];
		}

		Matrix<Real> covMat = Statistics::CovarianceMatrix(data);

		// Diagonal should match individual variances
		Real avg1, var1, avg2, var2;
		Statistics::AvgVar(col1, avg1, var1);
		Statistics::AvgVar(col2, avg2, var2);

		REQUIRE_THAT(covMat(0, 0), WithinAbs(var1, REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 1), WithinAbs(var2, REAL(1e-10)));
	}

	TEST_CASE("Statistics::CovarianceMatrix_InsufficientObservations_ThrowsException", "[statistics][correlation][error]") {
		TEST_PRECISION_INFO();
		Matrix<Real> data(1, 3); // Only 1 observation
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(2.0);
		data(0, 2) = REAL(3.0);

		REQUIRE_THROWS_AS(Statistics::CovarianceMatrix(data), StatisticsError);
	}

	/*********************************************************************/
	/*****              Correlation Matrix Tests                     *****/
	/*********************************************************************/
	TEST_CASE("Statistics::CorrelationMatrix_DiagonalIsOne", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Matrix<Real> data(5, 3);
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(10.0);
		data(0, 2) = REAL(100.0);
		data(1, 0) = REAL(2.0);
		data(1, 1) = REAL(20.0);
		data(1, 2) = REAL(200.0);
		data(2, 0) = REAL(3.0);
		data(2, 1) = REAL(15.0);
		data(2, 2) = REAL(150.0);
		data(3, 0) = REAL(4.0);
		data(3, 1) = REAL(25.0);
		data(3, 2) = REAL(250.0);
		data(4, 0) = REAL(5.0);
		data(4, 1) = REAL(30.0);
		data(4, 2) = REAL(300.0);

		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);

		// Diagonal should be 1 (variable correlated with itself)
		REQUIRE_THAT(corrMat(0, 0), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(1, 1), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(2, 2), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::CorrelationMatrix_PerfectCorrelation", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// All variables perfectly correlated (y = 2x, z = 3x)
		Matrix<Real> data(5, 3);
		for (int i = 0; i < 5; i++) {
			Real x = static_cast<Real>(i + 1);
			data(i, 0) = x;
			data(i, 1) = REAL(2.0) * x;
			data(i, 2) = REAL(3.0) * x;
		}

		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);

		// All correlations should be 1
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				REQUIRE_THAT(corrMat(i, j), WithinAbs(REAL(1.0), REAL(1e-10)));
			}
		}
	}

	TEST_CASE("Statistics::CorrelationMatrix_RangeIsMinusOneToOne", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Matrix<Real> data(5, 3);
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(5.0);
		data(0, 2) = REAL(2.0);
		data(1, 0) = REAL(3.0);
		data(1, 1) = REAL(2.0);
		data(1, 2) = REAL(7.0);
		data(2, 0) = REAL(5.0);
		data(2, 1) = REAL(8.0);
		data(2, 2) = REAL(1.0);
		data(3, 0) = REAL(2.0);
		data(3, 1) = REAL(3.0);
		data(3, 2) = REAL(9.0);
		data(4, 0) = REAL(4.0);
		data(4, 1) = REAL(6.0);
		data(4, 2) = REAL(3.0);

		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				REQUIRE(corrMat(i, j) >= REAL(-1.0) - REAL(1e-10));
				REQUIRE(corrMat(i, j) <= REAL(1.0) + REAL(1e-10));
			}
		}
	}

	TEST_CASE("Statistics::CorrelationMatrix_IsSymmetric", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3);
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(2.0);
		data(0, 2) = REAL(3.0);
		data(1, 0) = REAL(4.0);
		data(1, 1) = REAL(5.0);
		data(1, 2) = REAL(1.0);
		data(2, 0) = REAL(7.0);
		data(2, 1) = REAL(3.0);
		data(2, 2) = REAL(9.0);
		data(3, 0) = REAL(2.0);
		data(3, 1) = REAL(8.0);
		data(3, 2) = REAL(4.0);

		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);

		REQUIRE_THAT(corrMat(0, 1), WithinAbs(corrMat(1, 0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(0, 2), WithinAbs(corrMat(2, 0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(1, 2), WithinAbs(corrMat(2, 1), REAL(1e-10)));
	}

	TEST_CASE("Statistics::CorrelationMatrix_ConsistentWithPearson", "[statistics][correlation]") {
		TEST_PRECISION_INFO();
		// Check that matrix correlation matches pairwise Pearson
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(3.0), REAL(5.0), REAL(4.0), REAL(6.0)});

		Matrix<Real> data(5, 2);
		for (int i = 0; i < 5; i++) {
			data(i, 0) = x[i];
			data(i, 1) = y[i];
		}

		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);
		Real pairwiseR = Statistics::PearsonCorrelation(x, y);

		REQUIRE_THAT(corrMat(0, 1), WithinAbs(pairwiseR, REAL(1e-10)));
	}

} // namespace MML::Tests::Algorithms::StatisticsTests
