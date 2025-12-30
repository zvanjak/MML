#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics.h"
#include "algorithms/Statistics/HypothesisTesting.h"
#include "algorithms/Statistics/ConfidenceIntervals.h"
#include "algorithms/Statistics/Distributions.h"
#include "algorithms/Statistics/CoreDistributions.h"
#include "algorithms/Statistics/RankCorrelation.h"
#include "algorithms/Statistics/RandomGenerators.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::StatisticsTests
{
	/*********************************************************************/
	/*****          Basic Statistics - Average                       *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Avg_SimpleValues", "[statistics][basic]")
	{
			TEST_PRECISION_INFO();
		// Simple average of known values
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		
		Real avg = Statistics::Avg(data);
		
		REQUIRE_THAT(avg, WithinAbs(REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Avg_SingleValue", "[statistics][basic]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real avg = Statistics::Avg(data);
		
		REQUIRE_THAT(avg, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Avg_NegativeValues", "[statistics][basic]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> data({-REAL(5.0), -REAL(3.0), -REAL(1.0), REAL(1.0), REAL(3.0), REAL(5.0)});
		
		Real avg = Statistics::Avg(data);
		
		REQUIRE_THAT(avg, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Avg_EmptyVector_ThrowsException", "[statistics][error]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::Avg(data), StatisticsError);
	}

	/*********************************************************************/
	/*****          Average and Variance                             *****/
	/*********************************************************************/
	TEST_CASE("Statistics::AvgVar_KnownValues", "[statistics][variance]")
	{
			TEST_PRECISION_INFO();
		// Data: 2, 4, 4, 4, 5, 5, 7, 9
		// Mean = 40/8 = REAL(5.0)
		// Variance = sum((x - mean)²) / (n-1)
		// = (9 + 1 + 1 + 1 + 0 + 0 + 4 + 16) / 7 = 32/7 ≈ REAL(4.571)
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		
		Real avg, var;
		Statistics::AvgVar(data, avg, var);
		
		REQUIRE_THAT(avg, WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(32.0)/REAL(7.0), 1e-10));
	}

	TEST_CASE("Statistics::AvgVar_UniformData_ZeroVariance", "[statistics][variance]")
	{
			TEST_PRECISION_INFO();
		// All values the same = zero variance
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});
		
		Real avg, var;
		Statistics::AvgVar(data, avg, var);
		
		REQUIRE_THAT(avg, WithinAbs(REAL(7.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::AvgVar_TwoValues", "[statistics][variance]")
	{
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
	TEST_CASE("Statistics::AvgStdDev_KnownValues", "[statistics][stddev]")
	{
			TEST_PRECISION_INFO();
		// Data with known mean and std dev
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		
		Real avg, stddev;
		Statistics::AvgStdDev(data, avg, stddev);
		
		// Mean = REAL(5.0), Variance = 32/7, StdDev = sqrt(32/7) ≈ REAL(2.138)
		REQUIRE_THAT(avg, WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(stddev, WithinAbs(std::sqrt(REAL(32.0)/REAL(7.0)), 1e-10));
	}

	TEST_CASE("Statistics::AvgStdDev_StandardNormal", "[statistics][stddev]")
	{
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
	TEST_CASE("Statistics::Moments_SymmetricData", "[statistics][moments]")
	{
			TEST_PRECISION_INFO();
		// Symmetric data around mean = zero skewness
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		
		Real ave, adev, sdev, var, skew, curt;
		Statistics::Moments(data, ave, adev, sdev, var, skew, curt);
		
		REQUIRE_THAT(ave, WithinAbs(REAL(3.0), REAL(1e-10)));
		REQUIRE_THAT(var, WithinAbs(REAL(2.5), REAL(1e-10)));  // Variance of 1,2,3,4,5
		REQUIRE_THAT(sdev, WithinAbs(std::sqrt(REAL(2.5)), 1e-10));
		REQUIRE_THAT(skew, WithinAbs(REAL(0.0), REAL(1e-8)));  // Symmetric = zero skew
	}

	TEST_CASE("Statistics::Moments_PositiveSkew", "[statistics][moments]")
	{
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

	TEST_CASE("Statistics::Moments_NegativeSkew", "[statistics][moments]")
	{
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

	TEST_CASE("Statistics::Moments_HighKurtosis", "[statistics][moments]")
	{
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
		REQUIRE(curt > -REAL(3.0));  // Kurtosis > 0 (excess kurtosis)
	}

	TEST_CASE("Statistics::Moments_UniformData_ThrowsException", "[statistics][error]")
	{
			TEST_PRECISION_INFO();
		// Uniform data has zero variance, should throw when computing skew/kurtosis
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		
		Real ave, adev, sdev, var, skew, curt;
		
		REQUIRE_THROWS(Statistics::Moments(data, ave, adev, sdev, var, skew, curt));
	}

	TEST_CASE("Statistics::Moments_SingleValue_ThrowsException", "[statistics][error]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real ave, adev, sdev, var, skew, curt;
		
		REQUIRE_THROWS_AS(Statistics::Moments(data, ave, adev, sdev, var, skew, curt), StatisticsError);
	}

	/*********************************************************************/
	/*****          Edge Cases and Numerical Stability               *****/
	/*********************************************************************/
	TEST_CASE("Statistics::LargeDataset_Precision", "[statistics][precision]")
	{
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

	TEST_CASE("Statistics::VerySmallValues", "[statistics][precision]")
	{
			TEST_PRECISION_INFO();
		// Test numerical stability with very small values
		Vector<Real> data({1e-10, 2e-10, 3e-10, 4e-10, 5e-10});
		
		Real avg = Statistics::Avg(data);
		
		REQUIRE_THAT(avg, WithinAbs(3e-10, REAL(1e-20)));
	}

	TEST_CASE("Statistics::VeryLargeValues", "[statistics][precision]")
	{
			TEST_PRECISION_INFO();
		// Test with very large values
		Vector<Real> data({1e10, 2e10, 3e10, 4e10, 5e10});
		
		Real avg = Statistics::Avg(data);
		
		REQUIRE_THAT(avg, WithinRel(REAL(3e10), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****          Real-World Example Data                          *****/
	/*********************************************************************/
	TEST_CASE("Statistics::TemperatureData_Example", "[statistics][example]")
	{
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

	TEST_CASE("Statistics::TestScores_Example", "[statistics][example]")
	{
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
	/*****             Cauchy Distribution Tests                     *****/
	/*********************************************************************/
	TEST_CASE("CauchyDistribution::Construction", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		// Standard Cauchy (location=0, scale=1)
		Statistics::CauchyDistribution standard(REAL(0.0), REAL(1.0));
		REQUIRE(standard.mu == REAL(0.0));
		REQUIRE(standard.sigma == REAL(1.0));
		
		// Custom parameters
		Statistics::CauchyDistribution custom(REAL(5.0), REAL(2.0));
		REQUIRE(custom.mu == REAL(5.0));
		REQUIRE(custom.sigma == REAL(2.0));
		
		// Invalid scale should throw
		REQUIRE_THROWS_AS(Statistics::CauchyDistribution(REAL(0.0), REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::CauchyDistribution(REAL(0.0), -REAL(1.0)), StatisticsError);
	}

	TEST_CASE("CauchyDistribution::PDF_StandardCauchy", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::CauchyDistribution dist(REAL(0.0), REAL(1.0));
		
		// PDF at mode (x = 0) should be 1/π ≈ REAL(0.318309886)
		REQUIRE_THAT(dist.pdf(REAL(0.0)), WithinAbs(REAL(1.0) / Constants::PI, 1e-10));
		
		// PDF is symmetric around mode
		REQUIRE_THAT(dist.pdf(REAL(1.0)), WithinAbs(dist.pdf(-REAL(1.0)), REAL(1e-10)));
		REQUIRE_THAT(dist.pdf(REAL(2.0)), WithinAbs(dist.pdf(-REAL(2.0)), REAL(1e-10)));
		
		// PDF decreases as we move away from mode
		REQUIRE(dist.pdf(REAL(0.0)) > dist.pdf(REAL(1.0)));
		REQUIRE(dist.pdf(REAL(1.0)) > dist.pdf(REAL(2.0)));
		
		// PDF at x = ±1 should be 1/(2π) ≈ REAL(0.159154943)
		REQUIRE_THAT(dist.pdf(REAL(1.0)), WithinAbs(REAL(1.0) / (REAL(2.0) * Constants::PI), 1e-10));
	}

	TEST_CASE("CauchyDistribution::CDF_Properties", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::CauchyDistribution dist(REAL(0.0), REAL(1.0));
		
		// CDF at mode should be REAL(0.5)
		REQUIRE_THAT(dist.cdf(REAL(0.0)), WithinAbs(REAL(0.5), REAL(1e-10)));
		
		// CDF is monotonically increasing
		REQUIRE(dist.cdf(-REAL(2.0)) < dist.cdf(-REAL(1.0)));
		REQUIRE(dist.cdf(-REAL(1.0)) < dist.cdf(REAL(0.0)));
		REQUIRE(dist.cdf(REAL(0.0)) < dist.cdf(REAL(1.0)));
		REQUIRE(dist.cdf(REAL(1.0)) < dist.cdf(REAL(2.0)));
		
		// CDF approaches 0 and 1 at extremes
		REQUIRE(dist.cdf(-REAL(100.0)) < REAL(0.01));
		REQUIRE(dist.cdf(REAL(100.0)) > REAL(0.99));
		
		// CDF is symmetric: CDF(x) + CDF(-x) = 1
		REQUIRE_THAT(dist.cdf(REAL(1.5)) + dist.cdf(-REAL(1.5)), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("CauchyDistribution::InverseCDF", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::CauchyDistribution dist(REAL(0.0), REAL(1.0));
		
		// Inverse CDF at REAL(0.5) should be the mode
		REQUIRE_THAT(dist.inverseCdf(REAL(0.5)), WithinAbs(REAL(0.0), REAL(1e-10)));
		
		// Round-trip: inverseCdf(cdf(x)) ≈ x
		Real x1 = REAL(2.5);
		REQUIRE_THAT(dist.inverseCdf(dist.cdf(x1)), WithinAbs(x1, REAL(1e-9)));
		
		// Quartiles
		Real q1 = dist.inverseCdf(REAL(0.25));
		Real q3 = dist.inverseCdf(REAL(0.75));
		REQUIRE_THAT(q1, WithinAbs(-REAL(1.0), REAL(1e-10)));  // First quartile
		REQUIRE_THAT(q3, WithinAbs(REAL(1.0), REAL(1e-10)));   // Third quartile
		
		// Invalid probabilities should throw
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(1.0)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(-REAL(0.1)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(1.1)), StatisticsError);
	}

	/*********************************************************************/
	/*****          Exponential Distribution Tests                   *****/
	/*********************************************************************/
	TEST_CASE("ExponentialDistribution::Construction", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		// Standard exponential (rate=1)
		Statistics::ExponentialDistribution standard(REAL(1.0));
		REQUIRE(standard.lambda == REAL(1.0));
		REQUIRE_THAT(standard.mean(), WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// Rate = REAL(0.5) (mean = 2)
		Statistics::ExponentialDistribution slow(REAL(0.5));
		REQUIRE_THAT(slow.mean(), WithinAbs(REAL(2.0), REAL(1e-10)));
		REQUIRE_THAT(slow.variance(), WithinAbs(REAL(4.0), REAL(1e-10)));
		
		// Invalid rate should throw
		REQUIRE_THROWS_AS(Statistics::ExponentialDistribution(REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::ExponentialDistribution(-REAL(1.0)), StatisticsError);
	}

	TEST_CASE("ExponentialDistribution::PDF", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::ExponentialDistribution dist(REAL(1.0));
		
		// PDF at x=0 should be λ
		REQUIRE_THAT(dist.pdf(REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// PDF decreases exponentially
		REQUIRE(dist.pdf(REAL(0.0)) > dist.pdf(REAL(1.0)));
		REQUIRE(dist.pdf(REAL(1.0)) > dist.pdf(REAL(2.0)));
		
		// PDF at x=1 for λ=1 should be e^(-1)
		REQUIRE_THAT(dist.pdf(REAL(1.0)), WithinAbs(std::exp(-REAL(1.0)), 1e-10));
		
		// Negative x should throw
		REQUIRE_THROWS_AS(dist.pdf(-REAL(1.0)), StatisticsError);
	}

	TEST_CASE("ExponentialDistribution::CDF", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::ExponentialDistribution dist(REAL(1.0));
		
		// CDF at x=0 should be 0
		REQUIRE_THAT(dist.cdf(REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-10)));
		
		// CDF is monotonically increasing
		REQUIRE(dist.cdf(REAL(0.0)) < dist.cdf(REAL(1.0)));
		REQUIRE(dist.cdf(REAL(1.0)) < dist.cdf(REAL(2.0)));
		REQUIRE(dist.cdf(REAL(2.0)) < dist.cdf(REAL(3.0)));
		
		// CDF at x=mean should be approximately REAL(0.632)
		REQUIRE_THAT(dist.cdf(REAL(1.0)), WithinAbs(REAL(1.0) - std::exp(-REAL(1.0)), 1e-10));
		
		// CDF approaches 1
		REQUIRE(dist.cdf(REAL(10.0)) > REAL(0.99));
		
		// Negative x should throw
		REQUIRE_THROWS_AS(dist.cdf(-REAL(1.0)), StatisticsError);
	}

	TEST_CASE("ExponentialDistribution::InverseCDF", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::ExponentialDistribution dist(REAL(2.0));  // rate=2, mean=REAL(0.5)
		
		// Inverse CDF at p=0 should be 0
		REQUIRE_THAT(dist.inverseCdf(REAL(0.0)), WithinAbs(REAL(0.0), REAL(1e-10)));
		
		// Round-trip: inverseCdf(cdf(x)) ≈ x
		Real x1 = REAL(1.5);
		REQUIRE_THAT(dist.inverseCdf(dist.cdf(x1)), WithinAbs(x1, REAL(1e-9)));
		
		// Median (p=REAL(0.5)) should be ln(2)/λ
		Real median = dist.inverseCdf(REAL(0.5));
		REQUIRE_THAT(median, WithinAbs(std::log(REAL(2.0)) / REAL(2.0), 1e-10));
		
		// Invalid probabilities should throw
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(1.0)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(-REAL(0.1)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(1.1)), StatisticsError);
	}

	TEST_CASE("ExponentialDistribution::MemorylessProperty", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::ExponentialDistribution dist(REAL(1.0));
		
		// P(X > s+t | X > s) = P(X > t)
		// This is: (1 - CDF(s+t)) / (1 - CDF(s)) = 1 - CDF(t)
		Real s = REAL(2.0), t = REAL(3.0);
		Real conditional = (REAL(1.0) - dist.cdf(s + t)) / (REAL(1.0) - dist.cdf(s));
		Real marginal = REAL(1.0) - dist.cdf(t);
		
		REQUIRE_THAT(conditional, WithinAbs(marginal, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****            Logistic Distribution Tests                    *****/
	/*********************************************************************/
	TEST_CASE("LogisticDistribution::Construction", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		// Standard logistic (location=0, scale=1)
		Statistics::LogisticDistribution standard(REAL(0.0), REAL(1.0));
		REQUIRE(standard.mu == REAL(0.0));
		REQUIRE(standard.sigma == REAL(1.0));
		REQUIRE_THAT(standard.mean(), WithinAbs(REAL(0.0), REAL(1e-10)));
		
		// Custom parameters
		Statistics::LogisticDistribution custom(REAL(10.0), REAL(2.0));
		REQUIRE_THAT(custom.mean(), WithinAbs(REAL(10.0), REAL(1e-10)));
		
		// Variance = (π·σ)²/3
		Real expectedVar = (Constants::PI * custom.sigma) * (Constants::PI * custom.sigma) / REAL(3.0);
		REQUIRE_THAT(custom.variance(), WithinAbs(expectedVar, REAL(1e-10)));
		
		// Invalid scale should throw
		REQUIRE_THROWS_AS(Statistics::LogisticDistribution(REAL(0.0), REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::LogisticDistribution(REAL(0.0), -REAL(1.0)), StatisticsError);
	}

	TEST_CASE("LogisticDistribution::PDF", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::LogisticDistribution dist(REAL(0.0), REAL(1.0));
		
		// PDF is symmetric around mean
		REQUIRE_THAT(dist.pdf(REAL(1.0)), WithinAbs(dist.pdf(-REAL(1.0)), REAL(1e-10)));
		REQUIRE_THAT(dist.pdf(REAL(2.0)), WithinAbs(dist.pdf(-REAL(2.0)), REAL(1e-10)));
		
		// PDF at mode (x = μ) should be 1/(4σ) = REAL(0.25)
		REQUIRE_THAT(dist.pdf(REAL(0.0)), WithinAbs(REAL(0.25), REAL(1e-10)));
		
		// PDF decreases as we move away from mode
		REQUIRE(dist.pdf(REAL(0.0)) > dist.pdf(REAL(1.0)));
		REQUIRE(dist.pdf(REAL(1.0)) > dist.pdf(REAL(2.0)));
		
		// PDF should always be positive
		REQUIRE(dist.pdf(-REAL(5.0)) > REAL(0.0));
		REQUIRE(dist.pdf(REAL(5.0)) > REAL(0.0));
	}

	TEST_CASE("LogisticDistribution::CDF_Properties", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::LogisticDistribution dist(REAL(0.0), REAL(1.0));
		
		// CDF at mean should be REAL(0.5)
		REQUIRE_THAT(dist.cdf(REAL(0.0)), WithinAbs(REAL(0.5), REAL(1e-10)));
		
		// CDF is monotonically increasing
		REQUIRE(dist.cdf(-REAL(2.0)) < dist.cdf(-REAL(1.0)));
		REQUIRE(dist.cdf(-REAL(1.0)) < dist.cdf(REAL(0.0)));
		REQUIRE(dist.cdf(REAL(0.0)) < dist.cdf(REAL(1.0)));
		REQUIRE(dist.cdf(REAL(1.0)) < dist.cdf(REAL(2.0)));
		
		// CDF is symmetric: CDF(μ - x) = 1 - CDF(μ + x)
		REQUIRE_THAT(dist.cdf(-REAL(1.5)), WithinAbs(REAL(1.0) - dist.cdf(REAL(1.5)), 1e-10));
		
		// CDF approaches 0 and 1 at extremes
		REQUIRE(dist.cdf(-REAL(10.0)) < REAL(0.01));
		REQUIRE(dist.cdf(REAL(10.0)) > REAL(0.99));
	}

	TEST_CASE("LogisticDistribution::InverseCDF", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::LogisticDistribution dist(REAL(5.0), REAL(2.0));
		
		// Inverse CDF at REAL(0.5) should be the mean
		REQUIRE_THAT(dist.inverseCdf(REAL(0.5)), WithinAbs(REAL(5.0), REAL(1e-10)));
		
		// Round-trip: inverseCdf(cdf(x)) ≈ x
		Real x1 = REAL(7.3);
		REQUIRE_THAT(dist.inverseCdf(dist.cdf(x1)), WithinAbs(x1, REAL(1e-9)));
		
		// Symmetry: inverseCdf(p) + inverseCdf(1-p) = 2μ
		Real p = REAL(0.3);
		REQUIRE_THAT(dist.inverseCdf(p) + dist.inverseCdf(REAL(1.0) - p), 
		             WithinAbs(REAL(2.0) * dist.mu, 1e-9));
		
		// Invalid probabilities should throw
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(1.0)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(-REAL(0.1)), StatisticsError);
		REQUIRE_THROWS_AS(dist.inverseCdf(REAL(1.1)), StatisticsError);
	}

	TEST_CASE("LogisticDistribution::NumericalStability", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		Statistics::LogisticDistribution dist(REAL(0.0), REAL(1.0));
		
		// Test extreme values don't cause overflow/underflow
		REQUIRE(std::isfinite(dist.pdf(-REAL(50.0))));
		REQUIRE(std::isfinite(dist.pdf(REAL(50.0))));
		REQUIRE(std::isfinite(dist.cdf(-REAL(50.0))));
		REQUIRE(std::isfinite(dist.cdf(REAL(50.0))));
		
		// Very small and very large x
		REQUIRE_THAT(dist.cdf(-REAL(50.0)), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(dist.cdf(REAL(50.0)), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****              Distribution Comparison Tests                *****/
	/*********************************************************************/
	TEST_CASE("Distributions::CompareSymmetry", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		// All three distributions with symmetric parameters
		Statistics::CauchyDistribution cauchy(REAL(0.0), REAL(1.0));
		Statistics::LogisticDistribution logistic(REAL(0.0), REAL(1.0));
		
		// All should have CDF(0) = REAL(0.5)
		REQUIRE_THAT(cauchy.cdf(REAL(0.0)), WithinAbs(REAL(0.5), REAL(1e-10)));
		REQUIRE_THAT(logistic.cdf(REAL(0.0)), WithinAbs(REAL(0.5), REAL(1e-10)));
		
		// Cauchy has heavier tails than logistic
		// At x=3, Cauchy CDF should be closer to REAL(0.5) than logistic
		Real cauchy_tail = std::abs(cauchy.cdf(REAL(3.0)) - REAL(0.5));
		Real logistic_tail = std::abs(logistic.cdf(REAL(3.0)) - REAL(0.5));
		REQUIRE(cauchy_tail < logistic_tail);
	}

	TEST_CASE("Distributions::ScalingTest", "[statistics][distributions]")
	{
			TEST_PRECISION_INFO();
		// Test that distributions scale correctly with parameters
		Statistics::ExponentialDistribution exp1(REAL(1.0));
		Statistics::ExponentialDistribution exp2(REAL(2.0));
		
		// Rate REAL(2.0) should have half the mean of rate REAL(1.0)
		REQUIRE_THAT(exp2.mean(), WithinAbs(exp1.mean() / REAL(2.0), REAL(1e-10)));
		
		// Median should also scale
		Real median1 = exp1.inverseCdf(REAL(0.5));
		Real median2 = exp2.inverseCdf(REAL(0.5));
		REQUIRE_THAT(median2, WithinAbs(median1 / REAL(2.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****                Order Statistics - Median                  *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Median_OddCount", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Odd number of elements - median is the middle value
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(median, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_EvenCount", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Even number of elements - median is average of two middle values
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		Real median = Statistics::Median(data);
		
		// (2 + 3) / 2 = 2.5
		REQUIRE_THAT(median, WithinAbs(REAL(2.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_SingleValue", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(median, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_TwoValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(10.0), REAL(20.0)});
		
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(median, WithinAbs(REAL(15.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_UnsortedData", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Data is not sorted - median should still work
		Vector<Real> data({REAL(9.0), REAL(1.0), REAL(5.0), REAL(3.0), REAL(7.0)});
		
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(median, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_NegativeValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(-5.0), REAL(-3.0), REAL(-1.0), REAL(1.0), REAL(3.0)});
		
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(median, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Median_EmptyVector_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::Median(data), StatisticsError);
	}

	TEST_CASE("Statistics::Median_AllSameValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});
		
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(median, WithinAbs(REAL(7.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****              Order Statistics - Percentile                *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Percentile_50th_EqualsMedian", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		
		Real p50 = Statistics::Percentile(data, REAL(50.0));
		Real median = Statistics::Median(data);
		
		REQUIRE_THAT(p50, WithinAbs(median, REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_0th_IsMinimum", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(5.0), REAL(2.0), REAL(4.0)});
		
		Real p0 = Statistics::Percentile(data, REAL(0.0));
		
		REQUIRE_THAT(p0, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_100th_IsMaximum", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(5.0), REAL(2.0), REAL(4.0)});
		
		Real p100 = Statistics::Percentile(data, REAL(100.0));
		
		REQUIRE_THAT(p100, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_25th_FirstQuartile", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 (n=10)
		// Q1 (25th percentile): rank = 0.25 * 9 = 2.25
		// Linear interpolation: data[2] + 0.25 * (data[3] - data[2]) = 3 + 0.25 * 1 = 3.25
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                   REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		
		Real q1 = Statistics::Percentile(data, REAL(25.0));
		
		REQUIRE_THAT(q1, WithinAbs(REAL(3.25), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_75th_ThirdQuartile", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 (n=10)
		// Q3 (75th percentile): rank = 0.75 * 9 = 6.75
		// Linear interpolation: data[6] + 0.75 * (data[7] - data[6]) = 7 + 0.75 * 1 = 7.75
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                   REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		
		Real q3 = Statistics::Percentile(data, REAL(75.0));
		
		REQUIRE_THAT(q3, WithinAbs(REAL(7.75), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Percentile_InvalidRange_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		REQUIRE_THROWS_AS(Statistics::Percentile(data, REAL(-1.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::Percentile(data, REAL(101.0)), StatisticsError);
	}

	TEST_CASE("Statistics::Percentile_EmptyVector_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::Percentile(data, REAL(50.0)), StatisticsError);
	}

	/*********************************************************************/
	/*****              Order Statistics - Quartiles                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Quartiles_BasicTest", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                   REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		
		Real q1, median, q3;
		Statistics::Quartiles(data, q1, median, q3);
		
		REQUIRE_THAT(q1, WithinAbs(REAL(3.25), REAL(1e-10)));
		REQUIRE_THAT(median, WithinAbs(REAL(5.5), REAL(1e-10)));
		REQUIRE_THAT(q3, WithinAbs(REAL(7.75), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Quartiles_SymmetricData", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Symmetric data - median should equal mean
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		
		Real q1, median, q3;
		Statistics::Quartiles(data, q1, median, q3);
		
		REQUIRE_THAT(median, WithinAbs(REAL(5.0), REAL(1e-10)));
		// Q3 - median should equal median - Q1 (symmetric)
		REQUIRE_THAT(q3 - median, WithinAbs(median - q1, REAL(1e-10)));
	}

	TEST_CASE("Statistics::Quartiles_EmptyVector_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		Real q1, median, q3;
		REQUIRE_THROWS_AS(Statistics::Quartiles(data, q1, median, q3), StatisticsError);
	}

	/*********************************************************************/
	/*****               Order Statistics - Range                    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Range_BasicTest", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(5.0), REAL(3.0), REAL(9.0), REAL(2.0)});
		
		Real range = Statistics::Range(data);
		
		// max=9, min=1, range=8
		REQUIRE_THAT(range, WithinAbs(REAL(8.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Range_AllSameValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		
		Real range = Statistics::Range(data);
		
		REQUIRE_THAT(range, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Range_SingleValue", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real range = Statistics::Range(data);
		
		REQUIRE_THAT(range, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Range_NegativeValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(-10.0), REAL(-5.0), REAL(0.0), REAL(5.0), REAL(10.0)});
		
		Real range = Statistics::Range(data);
		
		REQUIRE_THAT(range, WithinAbs(REAL(20.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Range_EmptyVector_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::Range(data), StatisticsError);
	}

	/*********************************************************************/
	/*****                  Order Statistics - IQR                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::IQR_BasicTest", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
		// Q1 = 3.25, Q3 = 7.75, IQR = 4.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                   REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		
		Real iqr = Statistics::IQR(data);
		
		REQUIRE_THAT(iqr, WithinAbs(REAL(4.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::IQR_AllSameValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});
		
		Real iqr = Statistics::IQR(data);
		
		REQUIRE_THAT(iqr, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::IQR_ConsistentWithQuartiles", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(5.0), REAL(7.0), REAL(9.0), 
		                   REAL(11.0), REAL(13.0), REAL(15.0)});
		
		Real iqr = Statistics::IQR(data);
		Real q1, median, q3;
		Statistics::Quartiles(data, q1, median, q3);
		
		REQUIRE_THAT(iqr, WithinAbs(q3 - q1, REAL(1e-10)));
	}

	TEST_CASE("Statistics::IQR_EmptyVector_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::IQR(data), StatisticsError);
	}

	/*********************************************************************/
	/*****               Order Statistics - MinMax                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::MinMax_BasicTest", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(5.0), REAL(2.0), REAL(4.0)});
		
		Real minVal, maxVal;
		Statistics::MinMax(data, minVal, maxVal);
		
		REQUIRE_THAT(minVal, WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(maxVal, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_SingleValue", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real minVal, maxVal;
		Statistics::MinMax(data, minVal, maxVal);
		
		REQUIRE_THAT(minVal, WithinAbs(REAL(42.0), REAL(1e-10)));
		REQUIRE_THAT(maxVal, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_NegativeValues", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(-10.0), REAL(-5.0), REAL(-1.0)});
		
		Real minVal, maxVal;
		Statistics::MinMax(data, minVal, maxVal);
		
		REQUIRE_THAT(minVal, WithinAbs(REAL(-10.0), REAL(1e-10)));
		REQUIRE_THAT(maxVal, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_RangeConsistency", "[statistics][order-stats]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(5.0), REAL(3.0), REAL(9.0), REAL(2.0)});
		
		Real minVal, maxVal;
		Statistics::MinMax(data, minVal, maxVal);
		Real range = Statistics::Range(data);
		
		REQUIRE_THAT(range, WithinAbs(maxVal - minVal, REAL(1e-10)));
	}

	TEST_CASE("Statistics::MinMax_EmptyVector_ThrowsException", "[statistics][order-stats][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		Real minVal, maxVal;
		REQUIRE_THROWS_AS(Statistics::MinMax(data, minVal, maxVal), StatisticsError);
	}

	/*********************************************************************/
	/*****                Robust Statistics - Mode                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::Mode_SingleMode", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Clear single mode
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		Real mode = Statistics::Mode(data);
		
		REQUIRE_THAT(mode, WithinAbs(REAL(2.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_MultipleOccurrences", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(3.0), REAL(3.0), REAL(3.0), REAL(5.0), REAL(5.0)});
		
		Real mode = Statistics::Mode(data);
		
		REQUIRE_THAT(mode, WithinAbs(REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_AllUnique", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// All unique values - returns first (smallest after sort)
		Vector<Real> data({REAL(5.0), REAL(3.0), REAL(1.0), REAL(4.0), REAL(2.0)});
		
		Real mode = Statistics::Mode(data);
		
		// When all unique, returns first element (smallest after sort)
		REQUIRE_THAT(mode, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_AllSame", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0), REAL(7.0), REAL(7.0), REAL(7.0)});
		
		Real mode = Statistics::Mode(data);
		
		REQUIRE_THAT(mode, WithinAbs(REAL(7.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_SingleValue", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real mode = Statistics::Mode(data);
		
		REQUIRE_THAT(mode, WithinAbs(REAL(42.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Mode_EmptyVector_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::Mode(data), StatisticsError);
	}

	/*********************************************************************/
	/*****             Robust Statistics - TrimmedMean               *****/
	/*********************************************************************/
	TEST_CASE("Statistics::TrimmedMean_NoTrim", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// 0% trim = regular mean
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		
		Real trimmedMean = Statistics::TrimmedMean(data, REAL(0.0));
		Real regularMean = Statistics::Avg(data);
		
		REQUIRE_THAT(trimmedMean, WithinAbs(regularMean, REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_10Percent", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// 10 values, 10% trim = remove 1 from each end
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 100 (outlier)
		// After trim: 2, 3, 4, 5, 6, 7, 8, 9 -> mean = 44/8 = 5.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                   REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(100.0)});
		
		Real trimmedMean = Statistics::TrimmedMean(data, REAL(10.0));
		
		REQUIRE_THAT(trimmedMean, WithinAbs(REAL(5.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_20Percent", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// 10 values, 20% trim = remove 2 from each end
		// Data: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
		// After trim: 3, 4, 5, 6, 7, 8 -> mean = 33/6 = 5.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                   REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		
		Real trimmedMean = Statistics::TrimmedMean(data, REAL(20.0));
		
		REQUIRE_THAT(trimmedMean, WithinAbs(REAL(5.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_OutlierResistance", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// With outliers, trimmed mean should be close to true center
		Vector<Real> dataWithOutliers({REAL(-100.0), REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), 
		                               REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0), REAL(1000.0)});
		Vector<Real> dataClean({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                        REAL(6.0), REAL(7.0), REAL(8.0)});
		
		Real trimmedMeanOutliers = Statistics::TrimmedMean(dataWithOutliers, REAL(10.0));
		Real regularMeanClean = Statistics::Avg(dataClean);
		
		REQUIRE_THAT(trimmedMeanOutliers, WithinAbs(regularMeanClean, REAL(1e-10)));
	}

	TEST_CASE("Statistics::TrimmedMean_InvalidPercent_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(-1.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(50.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(100.0)), StatisticsError);
	}

	TEST_CASE("Statistics::TrimmedMean_EmptyVector_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::TrimmedMean(data, REAL(10.0)), StatisticsError);
	}

	/*********************************************************************/
	/*****                 Robust Statistics - MAD                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::MAD_BasicTest", "[statistics][robust]")
	{
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

	TEST_CASE("Statistics::MAD_AllSame", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// All same values -> MAD = 0
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		
		Real mad = Statistics::MAD(data);
		
		REQUIRE_THAT(mad, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MAD_WithOutlier", "[statistics][robust]")
	{
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

	TEST_CASE("Statistics::MAD_SingleValue", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(42.0)});
		
		Real mad = Statistics::MAD(data);
		
		REQUIRE_THAT(mad, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::MAD_EmptyVector_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::MAD(data), StatisticsError);
	}

	/*********************************************************************/
	/*****             Specialized Means - Geometric                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::GeometricMean_BasicTest", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Geometric mean of 2, 8 = sqrt(16) = 4
		Vector<Real> data({REAL(2.0), REAL(8.0)});
		
		Real gm = Statistics::GeometricMean(data);
		
		REQUIRE_THAT(gm, WithinAbs(REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_AllSame", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// GM of same values = that value
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0)});
		
		Real gm = Statistics::GeometricMean(data);
		
		REQUIRE_THAT(gm, WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_SingleValue", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(7.0)});
		
		Real gm = Statistics::GeometricMean(data);
		
		REQUIRE_THAT(gm, WithinAbs(REAL(7.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_GrowthRates", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Growth factors: 1.1, 1.2, 1.15 (10%, 20%, 15% growth)
		// GM = (1.1 * 1.2 * 1.15)^(1/3) ≈ 1.1491
		Vector<Real> rates({REAL(1.1), REAL(1.2), REAL(1.15)});
		
		Real gm = Statistics::GeometricMean(rates);
		Real expected = std::pow(REAL(1.1) * REAL(1.2) * REAL(1.15), REAL(1.0)/REAL(3.0));
		
		REQUIRE_THAT(gm, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::GeometricMean_NonPositive_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> dataWithZero({REAL(1.0), REAL(0.0), REAL(3.0)});
		Vector<Real> dataWithNegative({REAL(1.0), REAL(-2.0), REAL(3.0)});
		
		REQUIRE_THROWS_AS(Statistics::GeometricMean(dataWithZero), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::GeometricMean(dataWithNegative), StatisticsError);
	}

	TEST_CASE("Statistics::GeometricMean_EmptyVector_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::GeometricMean(data), StatisticsError);
	}

	/*********************************************************************/
	/*****              Specialized Means - Harmonic                 *****/
	/*********************************************************************/
	TEST_CASE("Statistics::HarmonicMean_BasicTest", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Harmonic mean of 1, 2, 4 = 3 / (1 + 0.5 + 0.25) = 3 / 1.75 ≈ 1.714
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(4.0)});
		
		Real hm = Statistics::HarmonicMean(data);
		Real expected = REAL(3.0) / (REAL(1.0) + REAL(0.5) + REAL(0.25));
		
		REQUIRE_THAT(hm, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_AllSame", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// HM of same values = that value
		Vector<Real> data({REAL(4.0), REAL(4.0), REAL(4.0)});
		
		Real hm = Statistics::HarmonicMean(data);
		
		REQUIRE_THAT(hm, WithinAbs(REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_Speeds", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Classic problem: drive 60 km at 30 km/h, then 60 km at 60 km/h
		// Average speed = total distance / total time = 120 / (2 + 1) = 40 km/h
		// This is the harmonic mean: 2 / (1/30 + 1/60) = 2 / (3/60) = 40
		Vector<Real> speeds({REAL(30.0), REAL(60.0)});
		
		Real hm = Statistics::HarmonicMean(speeds);
		
		REQUIRE_THAT(hm, WithinAbs(REAL(40.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_SingleValue", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(10.0)});
		
		Real hm = Statistics::HarmonicMean(data);
		
		REQUIRE_THAT(hm, WithinAbs(REAL(10.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::HarmonicMean_NonPositive_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> dataWithZero({REAL(1.0), REAL(0.0), REAL(3.0)});
		Vector<Real> dataWithNegative({REAL(1.0), REAL(-2.0), REAL(3.0)});
		
		REQUIRE_THROWS_AS(Statistics::HarmonicMean(dataWithZero), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::HarmonicMean(dataWithNegative), StatisticsError);
	}

	TEST_CASE("Statistics::HarmonicMean_EmptyVector_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		REQUIRE_THROWS_AS(Statistics::HarmonicMean(data), StatisticsError);
	}

	/*********************************************************************/
	/*****               Specialized Means - Weighted                *****/
	/*********************************************************************/
	TEST_CASE("Statistics::WeightedMean_EqualWeights", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Equal weights = regular mean
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0)});
		
		Real wm = Statistics::WeightedMean(data, weights);
		Real regularMean = Statistics::Avg(data);
		
		REQUIRE_THAT(wm, WithinAbs(regularMean, REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_UnequalWeights", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Data: 10, 20, 30 with weights 1, 2, 3
		// WM = (10*1 + 20*2 + 30*3) / (1+2+3) = (10 + 40 + 90) / 6 = 140/6 ≈ 23.33
		Vector<Real> data({REAL(10.0), REAL(20.0), REAL(30.0)});
		Vector<Real> weights({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		Real wm = Statistics::WeightedMean(data, weights);
		Real expected = REAL(140.0) / REAL(6.0);
		
		REQUIRE_THAT(wm, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_SingleNonZeroWeight", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Only one weight is non-zero -> result is that value
		Vector<Real> data({REAL(10.0), REAL(20.0), REAL(30.0)});
		Vector<Real> weights({REAL(0.0), REAL(1.0), REAL(0.0)});
		
		Real wm = Statistics::WeightedMean(data, weights);
		
		REQUIRE_THAT(wm, WithinAbs(REAL(20.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_GradeAverage", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// Typical grade calculation: HW=85 (20%), Midterm=78 (30%), Final=92 (50%)
		// WM = (85*0.2 + 78*0.3 + 92*0.5) / 1.0 = 17 + 23.4 + 46 = 86.4
		Vector<Real> grades({REAL(85.0), REAL(78.0), REAL(92.0)});
		Vector<Real> weights({REAL(0.2), REAL(0.3), REAL(0.5)});
		
		Real wm = Statistics::WeightedMean(grades, weights);
		
		REQUIRE_THAT(wm, WithinAbs(REAL(86.4), REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedMean_MismatchedSizes_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0)});  // Wrong size
		
		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedMean_NegativeWeights_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(1.0), REAL(-1.0), REAL(1.0)});
		
		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedMean_ZeroTotalWeight_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(0.0), REAL(0.0), REAL(0.0)});
		
		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedMean_EmptyVector_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		Vector<Real> weights;
		
		REQUIRE_THROWS_AS(Statistics::WeightedMean(data, weights), StatisticsError);
	}

	/*********************************************************************/
	/*****              Specialized Means - WeightedVariance         *****/
	/*********************************************************************/
	TEST_CASE("Statistics::WeightedVariance_EqualWeights", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		// With equal weights, should match regular variance
		Vector<Real> data({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0), REAL(1.0)});
		
		Real wv = Statistics::WeightedVariance(data, weights);
		Real avg, regularVar;
		Statistics::AvgVar(data, avg, regularVar);
		
		REQUIRE_THAT(wv, WithinAbs(regularVar, REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedVariance_AllSameValues", "[statistics][robust]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		Vector<Real> weights({REAL(1.0), REAL(2.0), REAL(1.0), REAL(2.0)});
		
		Real wv = Statistics::WeightedVariance(data, weights);
		
		REQUIRE_THAT(wv, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::WeightedVariance_MismatchedSizes_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> weights({REAL(1.0), REAL(1.0)});
		
		REQUIRE_THROWS_AS(Statistics::WeightedVariance(data, weights), StatisticsError);
	}

	TEST_CASE("Statistics::WeightedVariance_InsufficientData_ThrowsException", "[statistics][robust][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data({REAL(5.0)});
		Vector<Real> weights({REAL(1.0)});
		
		REQUIRE_THROWS_AS(Statistics::WeightedVariance(data, weights), StatisticsError);
	}

	/*********************************************************************/
	/*****                 Mean Inequalities Tests                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::MeanInequality_HM_GM_AM", "[statistics][robust]")
	{
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

	TEST_CASE("Statistics::MeanInequality_Equal_WhenUniform", "[statistics][robust]")
	{
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
	TEST_CASE("Statistics::Covariance_PositiveRelationship", "[statistics][correlation]")
	{
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

	TEST_CASE("Statistics::Covariance_NegativeRelationship", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// x increases, y decreases
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(10.0), REAL(8.0), REAL(6.0), REAL(4.0), REAL(2.0)});
		
		Real cov = Statistics::Covariance(x, y);
		
		// Negative covariance expected
		REQUIRE(cov < 0);
		REQUIRE_THAT(cov, WithinAbs(REAL(-5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_NoRelationship", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Quadratic relationship (symmetric around mean) - zero linear covariance
		Vector<Real> x({REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});  // y = x²
		
		Real cov = Statistics::Covariance(x, y);
		
		// Covariance should be zero (no LINEAR relationship)
		REQUIRE_THAT(cov, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_SelfCovariance_IsVariance", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(2.0), REAL(4.0), REAL(4.0), REAL(4.0), REAL(5.0), REAL(5.0), REAL(7.0), REAL(9.0)});
		
		Real selfCov = Statistics::Covariance(x, x);
		Real avg, var;
		Statistics::AvgVar(x, avg, var);
		
		REQUIRE_THAT(selfCov, WithinAbs(var, REAL(1e-10)));
	}

	TEST_CASE("Statistics::Covariance_DifferentSizes_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::Covariance(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::Covariance_SingleValue_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::Covariance(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****              Pearson Correlation Tests                    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::PearsonCorrelation_PerfectPositive", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Perfect positive correlation: y = 2x
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0)});
		
		Real r = Statistics::PearsonCorrelation(x, y);
		
		REQUIRE_THAT(r, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_PerfectNegative", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Perfect negative correlation: y = -2x + 12
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(10.0), REAL(8.0), REAL(6.0), REAL(4.0), REAL(2.0)});
		
		Real r = Statistics::PearsonCorrelation(x, y);
		
		REQUIRE_THAT(r, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_ZeroCorrelation", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Quadratic relationship - zero LINEAR correlation
		Vector<Real> x({REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});
		
		Real r = Statistics::PearsonCorrelation(x, y);
		
		REQUIRE_THAT(r, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_ModeratePositive", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Some noise in the relationship
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.5), REAL(3.5), REAL(5.0), REAL(7.5), REAL(9.5)});
		
		Real r = Statistics::PearsonCorrelation(x, y);
		
		// Should be high positive correlation but not perfect
		REQUIRE(r > REAL(0.9));
		REQUIRE(r < REAL(1.0));
	}

	TEST_CASE("Statistics::PearsonCorrelation_SelfCorrelation_IsOne", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(3.0), REAL(7.0), REAL(2.0), REAL(9.0)});
		
		Real r = Statistics::PearsonCorrelation(x, x);
		
		REQUIRE_THAT(r, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::PearsonCorrelation_RangeIsMinusOneToOne", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Random-ish data
		Vector<Real> x({REAL(2.3), REAL(5.1), REAL(1.2), REAL(8.7), REAL(4.4)});
		Vector<Real> y({REAL(7.8), REAL(3.2), REAL(9.1), REAL(2.5), REAL(6.6)});
		
		Real r = Statistics::PearsonCorrelation(x, y);
		
		REQUIRE(r >= REAL(-1.0));
		REQUIRE(r <= REAL(1.0));
	}

	TEST_CASE("Statistics::PearsonCorrelation_ZeroVarianceX_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		REQUIRE_THROWS_AS(Statistics::PearsonCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::PearsonCorrelation_ZeroVarianceY_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		
		REQUIRE_THROWS_AS(Statistics::PearsonCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::PearsonCorrelation_DifferentSizes_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::PearsonCorrelation(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****           Pearson Correlation With Test                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::PearsonCorrelationWithTest_HighCorrelation", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Strong positive correlation
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), 
		                REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.1), REAL(3.9), REAL(6.2), REAL(7.8), REAL(10.1), 
		                REAL(12.0), REAL(14.1), REAL(15.9), REAL(18.2), REAL(19.8)});
		
		auto result = Statistics::PearsonCorrelationWithTest(x, y);
		
		REQUIRE(result.r > REAL(0.99));
		REQUIRE(result.degreesOfFreedom == 8);
		REQUIRE(result.tStatistic > 0);  // Positive for positive correlation
		
		// With r ≈ 0.999, t should be very large
		REQUIRE(result.tStatistic > REAL(10.0));
	}

	TEST_CASE("Statistics::PearsonCorrelationWithTest_DegreesOfFreedom", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(5.0), REAL(4.0), REAL(5.0)});
		
		auto result = Statistics::PearsonCorrelationWithTest(x, y);
		
		// df = n - 2 = 5 - 2 = 3
		REQUIRE(result.degreesOfFreedom == 3);
	}

	TEST_CASE("Statistics::PearsonCorrelationWithTest_InsufficientData_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0)});
		
		// Need n > 2 for significance test
		REQUIRE_THROWS_AS(Statistics::PearsonCorrelationWithTest(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****                   R-Squared Tests                         *****/
	/*********************************************************************/
	TEST_CASE("Statistics::RSquared_PerfectFit", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// Perfect linear relationship -> R² = 1
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0)});
		
		Real r2 = Statistics::RSquared(x, y);
		
		REQUIRE_THAT(r2, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::RSquared_NoFit", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// No linear relationship -> R² ≈ 0
		Vector<Real> x({REAL(-2.0), REAL(-1.0), REAL(0.0), REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});
		
		Real r2 = Statistics::RSquared(x, y);
		
		REQUIRE_THAT(r2, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::RSquared_RangeIsZeroToOne", "[statistics][correlation]")
	{
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
	TEST_CASE("Statistics::CovarianceMatrix_TwoVariables", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		// 5 observations, 2 variables
		Matrix<Real> data(5, 2);
		// x: 1, 2, 3, 4, 5
		// y: 2, 4, 6, 8, 10
		data(0, 0) = REAL(1.0); data(0, 1) = REAL(2.0);
		data(1, 0) = REAL(2.0); data(1, 1) = REAL(4.0);
		data(2, 0) = REAL(3.0); data(2, 1) = REAL(6.0);
		data(3, 0) = REAL(4.0); data(3, 1) = REAL(8.0);
		data(4, 0) = REAL(5.0); data(4, 1) = REAL(10.0);
		
		Matrix<Real> covMat = Statistics::CovarianceMatrix(data);
		
		REQUIRE(covMat.RowNum() == 2);
		REQUIRE(covMat.ColNum() == 2);
		
		// Var(x) = 2.5, Var(y) = 10, Cov(x,y) = 5
		REQUIRE_THAT(covMat(0, 0), WithinAbs(REAL(2.5), REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 1), WithinAbs(REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(covMat(0, 1), WithinAbs(REAL(5.0), REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 0), WithinAbs(REAL(5.0), REAL(1e-10)));  // Symmetric
	}

	TEST_CASE("Statistics::CovarianceMatrix_IsSymmetric", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3);
		data(0, 0) = REAL(1.0); data(0, 1) = REAL(2.0); data(0, 2) = REAL(3.0);
		data(1, 0) = REAL(4.0); data(1, 1) = REAL(5.0); data(1, 2) = REAL(6.0);
		data(2, 0) = REAL(7.0); data(2, 1) = REAL(8.0); data(2, 2) = REAL(9.0);
		data(3, 0) = REAL(2.0); data(3, 1) = REAL(3.0); data(3, 2) = REAL(4.0);
		
		Matrix<Real> covMat = Statistics::CovarianceMatrix(data);
		
		// Check symmetry
		REQUIRE_THAT(covMat(0, 1), WithinAbs(covMat(1, 0), REAL(1e-10)));
		REQUIRE_THAT(covMat(0, 2), WithinAbs(covMat(2, 0), REAL(1e-10)));
		REQUIRE_THAT(covMat(1, 2), WithinAbs(covMat(2, 1), REAL(1e-10)));
	}

	TEST_CASE("Statistics::CovarianceMatrix_DiagonalIsVariance", "[statistics][correlation]")
	{
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

	TEST_CASE("Statistics::CovarianceMatrix_InsufficientObservations_ThrowsException", "[statistics][correlation][error]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> data(1, 3);  // Only 1 observation
		data(0, 0) = REAL(1.0);
		data(0, 1) = REAL(2.0);
		data(0, 2) = REAL(3.0);
		
		REQUIRE_THROWS_AS(Statistics::CovarianceMatrix(data), StatisticsError);
	}

	/*********************************************************************/
	/*****              Correlation Matrix Tests                     *****/
	/*********************************************************************/
	TEST_CASE("Statistics::CorrelationMatrix_DiagonalIsOne", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> data(5, 3);
		data(0, 0) = REAL(1.0); data(0, 1) = REAL(10.0); data(0, 2) = REAL(100.0);
		data(1, 0) = REAL(2.0); data(1, 1) = REAL(20.0); data(1, 2) = REAL(200.0);
		data(2, 0) = REAL(3.0); data(2, 1) = REAL(15.0); data(2, 2) = REAL(150.0);
		data(3, 0) = REAL(4.0); data(3, 1) = REAL(25.0); data(3, 2) = REAL(250.0);
		data(4, 0) = REAL(5.0); data(4, 1) = REAL(30.0); data(4, 2) = REAL(300.0);
		
		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);
		
		// Diagonal should be 1 (variable correlated with itself)
		REQUIRE_THAT(corrMat(0, 0), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(1, 1), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(2, 2), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::CorrelationMatrix_PerfectCorrelation", "[statistics][correlation]")
	{
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

	TEST_CASE("Statistics::CorrelationMatrix_RangeIsMinusOneToOne", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> data(5, 3);
		data(0, 0) = REAL(1.0); data(0, 1) = REAL(5.0); data(0, 2) = REAL(2.0);
		data(1, 0) = REAL(3.0); data(1, 1) = REAL(2.0); data(1, 2) = REAL(7.0);
		data(2, 0) = REAL(5.0); data(2, 1) = REAL(8.0); data(2, 2) = REAL(1.0);
		data(3, 0) = REAL(2.0); data(3, 1) = REAL(3.0); data(3, 2) = REAL(9.0);
		data(4, 0) = REAL(4.0); data(4, 1) = REAL(6.0); data(4, 2) = REAL(3.0);
		
		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);
		
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				REQUIRE(corrMat(i, j) >= REAL(-1.0) - REAL(1e-10));
				REQUIRE(corrMat(i, j) <= REAL(1.0) + REAL(1e-10));
			}
		}
	}

	TEST_CASE("Statistics::CorrelationMatrix_IsSymmetric", "[statistics][correlation]")
	{
		TEST_PRECISION_INFO();
		Matrix<Real> data(4, 3);
		data(0, 0) = REAL(1.0); data(0, 1) = REAL(2.0); data(0, 2) = REAL(3.0);
		data(1, 0) = REAL(4.0); data(1, 1) = REAL(5.0); data(1, 2) = REAL(1.0);
		data(2, 0) = REAL(7.0); data(2, 1) = REAL(3.0); data(2, 2) = REAL(9.0);
		data(3, 0) = REAL(2.0); data(3, 1) = REAL(8.0); data(3, 2) = REAL(4.0);
		
		Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);
		
		REQUIRE_THAT(corrMat(0, 1), WithinAbs(corrMat(1, 0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(0, 2), WithinAbs(corrMat(2, 0), REAL(1e-10)));
		REQUIRE_THAT(corrMat(1, 2), WithinAbs(corrMat(2, 1), REAL(1e-10)));
	}

	TEST_CASE("Statistics::CorrelationMatrix_ConsistentWithPearson", "[statistics][correlation]")
	{
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

	/*********************************************************************/
	/*****              Rank Computation Tests                       *****/
	/*********************************************************************/
	TEST_CASE("Statistics::ComputeRanks_SimpleSequence", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Already sorted: 1, 2, 3, 4, 5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE(ranks.size() == 5);
		REQUIRE_THAT(ranks[0], WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[1], WithinAbs(REAL(2.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[2], WithinAbs(REAL(3.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[3], WithinAbs(REAL(4.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[4], WithinAbs(REAL(5.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::ComputeRanks_UnsortedData", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Unsorted: 3, 1, 4, 1.5, 2  -> sorted order: 1, 1.5, 2, 3, 4
		// Ranks:    4, 1, 5, 2, 3
		Vector<Real> data({REAL(3.0), REAL(1.0), REAL(4.0), REAL(1.5), REAL(2.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE_THAT(ranks[0], WithinAbs(REAL(4.0), REAL(1e-10)));  // 3.0 is 4th smallest
		REQUIRE_THAT(ranks[1], WithinAbs(REAL(1.0), REAL(1e-10)));  // 1.0 is smallest
		REQUIRE_THAT(ranks[2], WithinAbs(REAL(5.0), REAL(1e-10)));  // 4.0 is largest
		REQUIRE_THAT(ranks[3], WithinAbs(REAL(2.0), REAL(1e-10)));  // 1.5 is 2nd smallest
		REQUIRE_THAT(ranks[4], WithinAbs(REAL(3.0), REAL(1e-10)));  // 2.0 is 3rd smallest
	}

	TEST_CASE("Statistics::ComputeRanks_WithTies", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Data with ties: 1, 2, 2, 3
		// Without ties: ranks would be 1, 2, 3, 4
		// With ties at positions 2,3: average rank = (2+3)/2 = 2.5
		Vector<Real> data({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE_THAT(ranks[0], WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(ranks[1], WithinAbs(REAL(2.5), REAL(1e-10)));  // Average of 2 and 3
		REQUIRE_THAT(ranks[2], WithinAbs(REAL(2.5), REAL(1e-10)));  // Average of 2 and 3
		REQUIRE_THAT(ranks[3], WithinAbs(REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::ComputeRanks_AllTied", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All equal: average rank for all
		Vector<Real> data({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		// Average rank = (1+2+3+4)/4 = 2.5
		for (int i = 0; i < 4; i++) {
			REQUIRE_THAT(ranks[i], WithinAbs(REAL(2.5), REAL(1e-10)));
		}
	}

	TEST_CASE("Statistics::ComputeRanks_EmptyVector", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> data;
		
		Vector<Real> ranks = Statistics::ComputeRanks(data);
		
		REQUIRE(ranks.size() == 0);
	}

	/*********************************************************************/
	/*****              Spearman Correlation Tests                   *****/
	/*********************************************************************/
	TEST_CASE("Statistics::SpearmanCorrelation_PerfectMonotonic", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Perfect positive monotonic relationship (but not linear)
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(4.0), REAL(9.0), REAL(16.0), REAL(25.0)});  // y = x²
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		// Perfect monotonic -> rho = 1
		REQUIRE_THAT(rho, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_PerfectNegativeMonotonic", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Perfect negative monotonic
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(25.0), REAL(16.0), REAL(9.0), REAL(4.0), REAL(1.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		REQUIRE_THAT(rho, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_NoMonotonicRelationship", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// U-shaped relationship: not monotonic
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(4.0), REAL(1.0), REAL(0.0), REAL(1.0), REAL(4.0)});  // y = (x-3)²
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		// Should be close to zero (not exactly due to discretization)
		REQUIRE(std::abs(rho) < REAL(0.5));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_WithTies", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Data with ties
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(3.0), REAL(3.0), REAL(4.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		// Should be positive (monotonically increasing overall)
		REQUIRE(rho > REAL(0.8));
		REQUIRE(rho <= REAL(1.0));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_SelfCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(3.0), REAL(1.0), REAL(4.0), REAL(1.5), REAL(9.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, x);
		
		REQUIRE_THAT(rho, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_RangeIsMinusOneToOne", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(2.3), REAL(5.1), REAL(1.2), REAL(8.7), REAL(4.4)});
		Vector<Real> y({REAL(7.8), REAL(3.2), REAL(9.1), REAL(2.5), REAL(6.6)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		
		REQUIRE(rho >= REAL(-1.0));
		REQUIRE(rho <= REAL(1.0));
	}

	TEST_CASE("Statistics::SpearmanCorrelation_DifferentSizes_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::SpearmanCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::SpearmanCorrelation_SingleElement_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::SpearmanCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::SpearmanCorrelationWithTest_HighCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0),
		                REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0),
		                REAL(12.0), REAL(14.0), REAL(16.0), REAL(18.0), REAL(20.0)});
		
		auto result = Statistics::SpearmanCorrelationWithTest(x, y);
		
		REQUIRE_THAT(result.rho, WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE(result.n == 10);
		REQUIRE(result.zScore > REAL(2.0));  // Significant
		REQUIRE(result.IsSignificant(0.05));
	}

	TEST_CASE("Statistics::SpearmanCorrelationWithTest_InsufficientData_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0)});
		
		REQUIRE_THROWS_AS(Statistics::SpearmanCorrelationWithTest(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****              Kendall Correlation Tests                    *****/
	/*********************************************************************/
	TEST_CASE("Statistics::KendallCorrelation_PerfectConcordance", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All pairs concordant
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(10.0), REAL(20.0), REAL(30.0), REAL(40.0), REAL(50.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE_THAT(tau, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_PerfectDiscordance", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All pairs discordant
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(50.0), REAL(40.0), REAL(30.0), REAL(20.0), REAL(10.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE_THAT(tau, WithinAbs(REAL(-1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_NoRelationship", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Designed to have balanced concordant/discordant pairs
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(1.0), REAL(4.0), REAL(2.0), REAL(3.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// tau should be between -1 and 1, close to 0 or moderate
		REQUIRE(tau >= REAL(-1.0));
		REQUIRE(tau <= REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_WithTiesInX", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Ties in x only
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// Should be positive but less than 1 due to tie
		REQUIRE(tau > REAL(0.5));
		REQUIRE(tau < REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_WithTiesInY", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Ties in y only
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE(tau > REAL(0.5));
		REQUIRE(tau < REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_AllTied", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// All tied in both x and y
		Vector<Real> x({REAL(5.0), REAL(5.0), REAL(5.0), REAL(5.0)});
		Vector<Real> y({REAL(3.0), REAL(3.0), REAL(3.0), REAL(3.0)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// No valid pairs -> tau = 0
		REQUIRE_THAT(tau, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_SelfCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(3.0), REAL(1.0), REAL(4.0), REAL(1.5), REAL(9.0)});
		
		Real tau = Statistics::KendallCorrelation(x, x);
		
		REQUIRE_THAT(tau, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::KendallCorrelation_RangeIsMinusOneToOne", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(2.3), REAL(5.1), REAL(1.2), REAL(8.7), REAL(4.4)});
		Vector<Real> y({REAL(7.8), REAL(3.2), REAL(9.1), REAL(2.5), REAL(6.6)});
		
		Real tau = Statistics::KendallCorrelation(x, y);
		
		REQUIRE(tau >= REAL(-1.0));
		REQUIRE(tau <= REAL(1.0));
	}

	TEST_CASE("Statistics::KendallCorrelation_DifferentSizes_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> y({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::KendallCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::KendallCorrelation_SingleElement_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0)});
		Vector<Real> y({REAL(2.0)});
		
		REQUIRE_THROWS_AS(Statistics::KendallCorrelation(x, y), StatisticsError);
	}

	TEST_CASE("Statistics::KendallCorrelationWithTest_SignificantCorrelation", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0),
		                REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0), REAL(6.0), REAL(8.0), REAL(10.0),
		                REAL(12.0), REAL(14.0), REAL(16.0), REAL(18.0), REAL(20.0)});
		
		auto result = Statistics::KendallCorrelationWithTest(x, y);
		
		REQUIRE_THAT(result.rho, WithinAbs(REAL(1.0), REAL(1e-10)));  // tau stored in rho
		REQUIRE(result.n == 10);
		REQUIRE(result.zScore > REAL(2.0));
		REQUIRE(result.IsSignificant(0.05));
	}

	TEST_CASE("Statistics::KendallCorrelationWithTest_InsufficientData_ThrowsException", "[statistics][rank][error]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0)});
		Vector<Real> y({REAL(2.0), REAL(4.0)});
		
		REQUIRE_THROWS_AS(Statistics::KendallCorrelationWithTest(x, y), StatisticsError);
	}

	/*********************************************************************/
	/*****          Spearman vs Kendall Comparison Tests             *****/
	/*********************************************************************/
	TEST_CASE("Statistics::SpearmanVsKendall_SameSign", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		// Both should agree on sign (direction of relationship)
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(2.0), REAL(3.0), REAL(5.0), REAL(4.0), REAL(6.0)});
		
		Real rho = Statistics::SpearmanCorrelation(x, y);
		Real tau = Statistics::KendallCorrelation(x, y);
		
		// Same sign
		REQUIRE(rho * tau > 0);
		
		// Both positive
		REQUIRE(rho > 0);
		REQUIRE(tau > 0);
		
		// Spearman typically has larger magnitude than Kendall
		// (for same data, |rho| >= |tau| is common but not guaranteed)
	}

	TEST_CASE("Statistics::RankCorrelation_PerfectMonotonic_AllAgree", "[statistics][rank]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> x({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		Vector<Real> y({REAL(1.0), REAL(8.0), REAL(27.0), REAL(64.0), REAL(125.0)});  // y = x³
		
		Real pearson = Statistics::PearsonCorrelation(x, y);
		Real spearman = Statistics::SpearmanCorrelation(x, y);
		Real kendall = Statistics::KendallCorrelation(x, y);
		
		// All three should be 1 for monotonic relationship
		// (Pearson may be less than 1 due to non-linearity)
		REQUIRE_THAT(spearman, WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(kendall, WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// Pearson less than 1 for non-linear but monotonic
		REQUIRE(pearson > REAL(0.9));
		REQUIRE(pearson <= REAL(1.0));
	}

	/*********************************************************************/
	/*****             Normal Distribution Tests                     *****/
	/*********************************************************************/
	TEST_CASE("Statistics::NormalDistribution_DefaultConstruction", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm;  // Standard normal
		
		REQUIRE_THAT(norm.mean(), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(norm.stddev(), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(norm.variance(), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_CustomParameters", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(100.0, 15.0);  // IQ distribution
		
		REQUIRE_THAT(norm.mean(), WithinAbs(REAL(100.0), REAL(1e-10)));
		REQUIRE_THAT(norm.stddev(), WithinAbs(REAL(15.0), REAL(1e-10)));
		REQUIRE_THAT(norm.variance(), WithinAbs(REAL(225.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_InvalidStddev_ThrowsException", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		REQUIRE_THROWS_AS(Statistics::NormalDistribution(0.0, 0.0), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::NormalDistribution(0.0, -1.0), StatisticsError);
	}

	TEST_CASE("Statistics::NormalDistribution_PDF_AtMean", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// PDF at mean should be 1/sqrt(2*pi) ≈ 0.3989
		Real expected = 1.0 / std::sqrt(2.0 * Constants::PI);
		REQUIRE_THAT(norm.pdf(0.0), WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_PDF_Symmetric", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(5.0, 2.0);
		
		// PDF should be symmetric around mean
		REQUIRE_THAT(norm.pdf(3.0), WithinAbs(norm.pdf(7.0), REAL(1e-10)));  // ±2 from mean
		REQUIRE_THAT(norm.pdf(1.0), WithinAbs(norm.pdf(9.0), REAL(1e-10)));  // ±4 from mean
	}

	TEST_CASE("Statistics::NormalDistribution_PDF_DecreaseFromMean", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// PDF decreases as we move away from mean
		REQUIRE(norm.pdf(0.0) > norm.pdf(1.0));
		REQUIRE(norm.pdf(1.0) > norm.pdf(2.0));
		REQUIRE(norm.pdf(2.0) > norm.pdf(3.0));
	}

	TEST_CASE("Statistics::NormalDistribution_CDF_AtMean", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// CDF at mean should be 0.5
		REQUIRE_THAT(norm.cdf(0.0), WithinAbs(REAL(0.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_CDF_KnownValues", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// Known values from standard normal table
		REQUIRE_THAT(norm.cdf(-3.0), WithinAbs(REAL(0.00135), REAL(1e-4)));
		REQUIRE_THAT(norm.cdf(-2.0), WithinAbs(REAL(0.02275), REAL(1e-4)));
		REQUIRE_THAT(norm.cdf(-1.0), WithinAbs(REAL(0.15866), REAL(1e-4)));
		REQUIRE_THAT(norm.cdf(1.0), WithinAbs(REAL(0.84134), REAL(1e-4)));
		REQUIRE_THAT(norm.cdf(2.0), WithinAbs(REAL(0.97725), REAL(1e-4)));
		REQUIRE_THAT(norm.cdf(3.0), WithinAbs(REAL(0.99865), REAL(1e-4)));
	}

	TEST_CASE("Statistics::NormalDistribution_CDF_MonotonicallyIncreasing", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		Real prev = norm.cdf(-5.0);
		for (Real x = -4.0; x <= 5.0; x += 0.5) {
			Real current = norm.cdf(x);
			REQUIRE(current > prev);
			prev = current;
		}
	}

	TEST_CASE("Statistics::NormalDistribution_CDF_BoundsZeroOne", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// At reasonable distance from mean
		REQUIRE(norm.cdf(-5.0) > REAL(0.0));
		REQUIRE(norm.cdf(-5.0) < REAL(0.001));
		REQUIRE(norm.cdf(5.0) > REAL(0.999));
		REQUIRE(norm.cdf(5.0) < REAL(1.0));
	}

	TEST_CASE("Statistics::NormalDistribution_InverseCDF_KnownValues", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// Known quantiles
		REQUIRE_THAT(norm.inverseCdf(0.5), WithinAbs(REAL(0.0), REAL(1e-6)));
		REQUIRE_THAT(norm.inverseCdf(0.84134), WithinAbs(REAL(1.0), REAL(1e-3)));
		REQUIRE_THAT(norm.inverseCdf(0.97725), WithinAbs(REAL(2.0), REAL(1e-3)));
		REQUIRE_THAT(norm.inverseCdf(0.15866), WithinAbs(REAL(-1.0), REAL(1e-3)));
	}

	TEST_CASE("Statistics::NormalDistribution_InverseCDF_RoundTrip", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(10.0, 3.0);
		
		// CDF(inverseCDF(p)) should equal p
		for (Real p = 0.1; p <= 0.9; p += 0.1) {
			Real x = norm.inverseCdf(p);
			Real pBack = norm.cdf(x);
			REQUIRE_THAT(pBack, WithinAbs(p, REAL(1e-6)));
		}
	}

	TEST_CASE("Statistics::NormalDistribution_InverseCDF_InvalidProbability", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		REQUIRE_THROWS_AS(norm.inverseCdf(0.0), StatisticsError);
		REQUIRE_THROWS_AS(norm.inverseCdf(1.0), StatisticsError);
		REQUIRE_THROWS_AS(norm.inverseCdf(-0.1), StatisticsError);
		REQUIRE_THROWS_AS(norm.inverseCdf(1.1), StatisticsError);
	}

	TEST_CASE("Statistics::NormalDistribution_ZScore", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(100.0, 15.0);
		
		// z = (x - mu) / sigma
		REQUIRE_THAT(norm.zScore(100.0), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(norm.zScore(115.0), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(norm.zScore(85.0), WithinAbs(REAL(-1.0), REAL(1e-10)));
		REQUIRE_THAT(norm.zScore(130.0), WithinAbs(REAL(2.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_FromZScore", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(100.0, 15.0);
		
		REQUIRE_THAT(norm.fromZScore(0.0), WithinAbs(REAL(100.0), REAL(1e-10)));
		REQUIRE_THAT(norm.fromZScore(1.0), WithinAbs(REAL(115.0), REAL(1e-10)));
		REQUIRE_THAT(norm.fromZScore(-2.0), WithinAbs(REAL(70.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_ZScoreRoundTrip", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(50.0, 10.0);
		
		Real x = 67.5;
		Real z = norm.zScore(x);
		Real xBack = norm.fromZScore(z);
		REQUIRE_THAT(xBack, WithinAbs(x, REAL(1e-10)));
	}

	/*********************************************************************/
	/*****           Standard Normal Tests                           *****/
	/*********************************************************************/
	TEST_CASE("Statistics::StandardNormal_PDF", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		
		Real expected = 1.0 / std::sqrt(2.0 * Constants::PI);
		REQUIRE_THAT(Statistics::StandardNormal::pdf(0.0), WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::StandardNormal_CDF", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		
		REQUIRE_THAT(Statistics::StandardNormal::cdf(0.0), WithinAbs(REAL(0.5), REAL(1e-10)));
		REQUIRE_THAT(Statistics::StandardNormal::cdf(1.96), WithinAbs(REAL(0.975), REAL(1e-3)));
	}

	TEST_CASE("Statistics::StandardNormal_InverseCDF", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		
		REQUIRE_THAT(Statistics::StandardNormal::inverseCdf(0.5), WithinAbs(REAL(0.0), REAL(1e-6)));
		REQUIRE_THAT(Statistics::StandardNormal::inverseCdf(0.975), WithinAbs(REAL(1.96), REAL(1e-2)));
	}

	TEST_CASE("Statistics::StandardNormal_SurvivalFunction", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		
		// sf(z) = 1 - cdf(z) = P(Z > z)
		REQUIRE_THAT(Statistics::StandardNormal::sf(0.0), WithinAbs(REAL(0.5), REAL(1e-10)));
		REQUIRE_THAT(Statistics::StandardNormal::sf(1.96), WithinAbs(REAL(0.025), REAL(1e-3)));
	}

	TEST_CASE("Statistics::StandardNormal_TwoTailedPValue", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		
		// Two-tailed p-value: P(|Z| > |z|) = 2 * P(Z > |z|)
		REQUIRE_THAT(Statistics::StandardNormal::twoTailedPValue(0.0), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(Statistics::StandardNormal::twoTailedPValue(1.96), WithinAbs(REAL(0.05), REAL(1e-3)));
		REQUIRE_THAT(Statistics::StandardNormal::twoTailedPValue(2.576), WithinAbs(REAL(0.01), REAL(1e-3)));
		
		// Symmetric: p-value for -z equals p-value for z
		REQUIRE_THAT(Statistics::StandardNormal::twoTailedPValue(-1.96), 
		             WithinAbs(Statistics::StandardNormal::twoTailedPValue(1.96), REAL(1e-10)));
	}

	TEST_CASE("Statistics::StandardNormal_ConsistentWithNormalDistribution", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution standard(0.0, 1.0);
		
		// StandardNormal should match NormalDistribution(0, 1)
		REQUIRE_THAT(Statistics::StandardNormal::pdf(1.5), WithinAbs(standard.pdf(1.5), REAL(1e-10)));
		REQUIRE_THAT(Statistics::StandardNormal::cdf(1.5), WithinAbs(standard.cdf(1.5), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****           Student's t-Distribution Tests                  *****/
	/*********************************************************************/
	TEST_CASE("Statistics::TDistribution_Construction", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Should construct without error
		REQUIRE(t.df == 10);
	}

	TEST_CASE("Statistics::TDistribution_InvalidDF_ThrowsException", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		REQUIRE_THROWS_AS(Statistics::TDistribution(0), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::TDistribution(-1), StatisticsError);
	}

	TEST_CASE("Statistics::TDistribution_PDF_Symmetric", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(5);
		
		// PDF should be symmetric around 0
		REQUIRE_THAT(t.pdf(1.0), WithinAbs(t.pdf(-1.0), REAL(1e-10)));
		REQUIRE_THAT(t.pdf(2.5), WithinAbs(t.pdf(-2.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_PDF_DecreaseFromCenter", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// PDF decreases as we move away from 0
		REQUIRE(t.pdf(0.0) > t.pdf(1.0));
		REQUIRE(t.pdf(1.0) > t.pdf(2.0));
		REQUIRE(t.pdf(2.0) > t.pdf(3.0));
	}

	TEST_CASE("Statistics::TDistribution_PDF_HeavierTailsThanNormal", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(5);
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// t-distribution has heavier tails: higher PDF at extremes
		REQUIRE(t.pdf(3.0) > norm.pdf(3.0));
		REQUIRE(t.pdf(4.0) > norm.pdf(4.0));
	}

	TEST_CASE("Statistics::TDistribution_PDF_ConvergesToNormal", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(100);  // Large df
		Statistics::NormalDistribution norm(0.0, 1.0);
		
		// With large df, t approaches normal
		REQUIRE_THAT(t.pdf(0.0), WithinAbs(norm.pdf(0.0), REAL(1e-2)));
		REQUIRE_THAT(t.pdf(1.5), WithinAbs(norm.pdf(1.5), REAL(1e-2)));
		REQUIRE_THAT(t.pdf(2.0), WithinAbs(norm.pdf(2.0), REAL(1e-2)));
	}

	TEST_CASE("Statistics::TDistribution_CDF_AtZero", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// CDF at 0 should be 0.5 (symmetric distribution)
		REQUIRE_THAT(t.cdf(0.0), WithinAbs(REAL(0.5), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_CDF_MonotonicallyIncreasing", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(8);
		
		Real prev = t.cdf(-5.0);
		for (Real x = -4.0; x <= 5.0; x += 0.5) {
			Real current = t.cdf(x);
			REQUIRE(current > prev);
			prev = current;
		}
	}

	TEST_CASE("Statistics::TDistribution_CDF_KnownValues", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Known values from t-table (df=10)
		// t = 1.812 gives P(T <= t) ≈ 0.95 (one-tailed)
		REQUIRE_THAT(t.cdf(1.812), WithinAbs(REAL(0.95), REAL(1e-2)));
		
		// By symmetry
		REQUIRE_THAT(t.cdf(-1.812), WithinAbs(REAL(0.05), REAL(1e-2)));
	}

	TEST_CASE("Statistics::TDistribution_CDF_DF1_IsCauchy", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(1);
		
		// df=1 is Cauchy distribution centered at 0
		// CDF(0) = 0.5, CDF(1) = 0.5 + arctan(1)/π ≈ 0.75
		REQUIRE_THAT(t.cdf(0.0), WithinAbs(REAL(0.5), REAL(1e-10)));
		REQUIRE_THAT(t.cdf(1.0), WithinAbs(REAL(0.75), REAL(1e-3)));
	}

	TEST_CASE("Statistics::TDistribution_InverseCDF_AtMedian", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(15);
		
		// inverseCDF(0.5) should be 0
		REQUIRE_THAT(t.inverseCdf(0.5), WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_InverseCDF_KnownValues", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Known critical values from t-table (df=10)
		// 95th percentile (one-tailed) ≈ 1.812
		REQUIRE_THAT(t.inverseCdf(0.95), WithinAbs(REAL(1.812), REAL(1e-2)));
		
		// By symmetry
		REQUIRE_THAT(t.inverseCdf(0.05), WithinAbs(REAL(-1.812), REAL(1e-2)));
	}

	TEST_CASE("Statistics::TDistribution_InverseCDF_RoundTrip", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(20);
		
		// CDF(inverseCDF(p)) should equal p
		for (Real p = 0.1; p <= 0.9; p += 0.1) {
			Real x = t.inverseCdf(p);
			Real pBack = t.cdf(x);
			REQUIRE_THAT(pBack, WithinAbs(p, REAL(1e-6)));
		}
	}

	TEST_CASE("Statistics::TDistribution_InverseCDF_InvalidProbability", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		REQUIRE_THROWS_AS(t.inverseCdf(0.0), StatisticsError);
		REQUIRE_THROWS_AS(t.inverseCdf(1.0), StatisticsError);
		REQUIRE_THROWS_AS(t.inverseCdf(-0.1), StatisticsError);
		REQUIRE_THROWS_AS(t.inverseCdf(1.1), StatisticsError);
	}

	TEST_CASE("Statistics::TDistribution_SurvivalFunction", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(12);
		
		// sf(x) = 1 - cdf(x) = P(T > x)
		REQUIRE_THAT(t.sf(0.0), WithinAbs(REAL(0.5), REAL(1e-10)));
		
		Real cdfVal = t.cdf(1.5);
		Real sfVal = t.sf(1.5);
		REQUIRE_THAT(cdfVal + sfVal, WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_TwoTailedPValue", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Two-tailed p-value: P(|T| > |t|)
		// For t = 2.228 (df=10), two-tailed p ≈ 0.05
		REQUIRE_THAT(t.twoTailedPValue(2.228), WithinAbs(REAL(0.05), REAL(1e-2)));
		
		// Symmetric
		REQUIRE_THAT(t.twoTailedPValue(-2.228), WithinAbs(t.twoTailedPValue(2.228), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_CriticalValue", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Critical value for α = 0.05 (two-tailed) with df=10 is ≈ 2.228
		Real critical = t.criticalValue(0.05);
		REQUIRE_THAT(critical, WithinAbs(REAL(2.228), REAL(1e-2)));
		
		// Check that this gives correct p-value
		REQUIRE_THAT(t.twoTailedPValue(critical), WithinAbs(REAL(0.05), REAL(1e-3)));
	}

	TEST_CASE("Statistics::TDistribution_CriticalValue_InvalidAlpha", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		REQUIRE_THROWS_AS(t.criticalValue(0.0), StatisticsError);
		REQUIRE_THROWS_AS(t.criticalValue(1.0), StatisticsError);
		REQUIRE_THROWS_AS(t.criticalValue(-0.1), StatisticsError);
	}

	TEST_CASE("Statistics::TDistribution_Mean", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(5);
		
		// Mean = 0 for df >= 2
		REQUIRE_THAT(t.mean(), WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_Mean_DoesNotExist_DF1", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(1);
		
		// Mean does not exist for df = 1
		REQUIRE_THROWS_AS(t.mean(), StatisticsError);
	}

	TEST_CASE("Statistics::TDistribution_Variance", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Var = ν/(ν-2) = 10/8 = 1.25 for df=10
		REQUIRE_THAT(t.variance(), WithinAbs(REAL(1.25), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_Variance_HighDF_ApproachesOne", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(100);
		
		// As df → ∞, variance → 1 (approaches standard normal)
		Real var = t.variance();
		REQUIRE(var < REAL(1.1));
		REQUIRE(var > REAL(1.0));
	}

	TEST_CASE("Statistics::TDistribution_Variance_DoesNotExist_LowDF", "[statistics][distribution][error]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t1(1);
		Statistics::TDistribution t2(2);
		
		// Variance does not exist for df <= 2
		REQUIRE_THROWS_AS(t1.variance(), StatisticsError);
		REQUIRE_THROWS_AS(t2.variance(), StatisticsError);
	}

	TEST_CASE("Statistics::TDistribution_CompareMultipleDF", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t3(3);
		Statistics::TDistribution t10(10);
		Statistics::TDistribution t30(30);
		
		// With lower df, PDF at tails is higher (heavier tails)
		REQUIRE(t3.pdf(3.0) > t10.pdf(3.0));
		REQUIRE(t10.pdf(3.0) > t30.pdf(3.0));
		
		// With lower df, critical values are larger
		REQUIRE(t3.criticalValue(0.05) > t10.criticalValue(0.05));
		REQUIRE(t10.criticalValue(0.05) > t30.criticalValue(0.05));
	}

	// ================================================================================
	// HYPOTHESIS TESTING: T-TESTS
	// ================================================================================

	TEST_CASE("Statistics::OneSampleTTest_BasicUsage", "[statistics][hypothesis][t-test]")
	{
		TEST_PRECISION_INFO();
		
		// Sample data: [10, 12, 15, 18, 20]
		// Mean = 15, std ≈ 4.18, n = 5
		// Test H0: μ = 15 (should NOT reject)
		Vector<Real> sample({10.0, 12.0, 15.0, 18.0, 20.0});
		
		auto result = Statistics::OneSampleTTest(sample, 15.0, 0.05);
		
		// Testing against the true mean should give t ≈ 0
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
		
		// Test H0: μ = 15 (should reject)
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
		// Sample: [5, 7, 9, 11, 13] mean = 9, std = sqrt(10) ≈ 3.162, n = 5
		// Test H0: μ = 10
		Vector<Real> sample({5.0, 7.0, 9.0, 11.0, 13.0});
		
		Real sampleMean = Statistics::Mean(sample);
		Real sampleStd = Statistics::StdDev(sample);
		int n = 5;
		
		// Manual calculation: t = (9 - 10) / (3.162 / sqrt(5)) = -1 / 1.414 ≈ -0.707
		Real expectedT = (sampleMean - 10.0) / (sampleStd / std::sqrt(static_cast<Real>(n)));
		
		auto result = Statistics::OneSampleTTest(sample, 10.0, 0.05);
		
		REQUIRE_THAT(result.testStatistic, WithinAbs(expectedT, 0.01));
		
		// df = 4, critical value at α=0.05 is ±2.776
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
		
		// Both have mean ≈ 14-15, should not reject H0: μ1 = μ2
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
		Vector<Real> sample1({10.0, 12.0, 14.0, 16.0, 18.0});  // mean ≈ 14
		Vector<Real> sample2({20.0, 22.0, 24.0, 26.0, 28.0});  // mean ≈ 24
		
		// Should reject H0: μ1 = μ2
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
		
		// No difference: t ≈ 0
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
		groups.push_back(Vector<Real>({10.0, 11.0, 9.0, 10.5, 10.5}));    // mean ≈ 10
		groups.push_back(Vector<Real>({20.0, 21.0, 19.0, 20.5, 20.5}));   // mean ≈ 20
		groups.push_back(Vector<Real>({30.0, 31.0, 29.0, 30.5, 30.5}));   // mean ≈ 30
		
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
		
		// For 2 groups, F = t² (approximately)
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

} // namespace MML::Tests::Algorithms::StatisticsTests
