#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Statistics.h"
#include "algorithms/Statistics/Distributions.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::DistributionsTests
{
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
		
		// PDF at mode (x = 0) should be 1/Ï€ â‰ˆ REAL(0.318309886)
		REQUIRE_THAT(dist.pdf(REAL(0.0)), WithinAbs(REAL(1.0) / Constants::PI, 1e-10));
		
		// PDF is symmetric around mode
		REQUIRE_THAT(dist.pdf(REAL(1.0)), WithinAbs(dist.pdf(-REAL(1.0)), REAL(1e-10)));
		REQUIRE_THAT(dist.pdf(REAL(2.0)), WithinAbs(dist.pdf(-REAL(2.0)), REAL(1e-10)));
		
		// PDF decreases as we move away from mode
		REQUIRE(dist.pdf(REAL(0.0)) > dist.pdf(REAL(1.0)));
		REQUIRE(dist.pdf(REAL(1.0)) > dist.pdf(REAL(2.0)));
		
		// PDF at x = Â±1 should be 1/(2Ï€) â‰ˆ REAL(0.159154943)
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
		
		// Round-trip: inverseCdf(cdf(x)) â‰ˆ x
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
		
		// PDF at x=0 should be Î»
		REQUIRE_THAT(dist.pdf(REAL(0.0)), WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// PDF decreases exponentially
		REQUIRE(dist.pdf(REAL(0.0)) > dist.pdf(REAL(1.0)));
		REQUIRE(dist.pdf(REAL(1.0)) > dist.pdf(REAL(2.0)));
		
		// PDF at x=1 for Î»=1 should be e^(-1)
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
		
		// Round-trip: inverseCdf(cdf(x)) â‰ˆ x
		Real x1 = REAL(1.5);
		REQUIRE_THAT(dist.inverseCdf(dist.cdf(x1)), WithinAbs(x1, REAL(1e-9)));
		
		// Median (p=REAL(0.5)) should be ln(2)/Î»
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
		
		// Variance = (Ï€Â·Ïƒ)Â²/3
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
		
		// PDF at mode (x = Î¼) should be 1/(4Ïƒ) = REAL(0.25)
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
		
		// CDF is symmetric: CDF(Î¼ - x) = 1 - CDF(Î¼ + x)
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
		
		// Round-trip: inverseCdf(cdf(x)) â‰ˆ x
		Real x1 = REAL(7.3);
		REQUIRE_THAT(dist.inverseCdf(dist.cdf(x1)), WithinAbs(x1, REAL(1e-9)));
		
		// Symmetry: inverseCdf(p) + inverseCdf(1-p) = 2Î¼
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
		
		// PDF at mean should be 1/sqrt(2*pi) â‰ˆ 0.3989
		Real expected = 1.0 / std::sqrt(2.0 * Constants::PI);
		REQUIRE_THAT(norm.pdf(0.0), WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("Statistics::NormalDistribution_PDF_Symmetric", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::NormalDistribution norm(5.0, 2.0);
		
		// PDF should be symmetric around mean
		REQUIRE_THAT(norm.pdf(3.0), WithinAbs(norm.pdf(7.0), REAL(1e-10)));  // Â±2 from mean
		REQUIRE_THAT(norm.pdf(1.0), WithinAbs(norm.pdf(9.0), REAL(1e-10)));  // Â±4 from mean
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
		// t = 1.812 gives P(T <= t) â‰ˆ 0.95 (one-tailed)
		REQUIRE_THAT(t.cdf(1.812), WithinAbs(REAL(0.95), REAL(1e-2)));
		
		// By symmetry
		REQUIRE_THAT(t.cdf(-1.812), WithinAbs(REAL(0.05), REAL(1e-2)));
	}

	TEST_CASE("Statistics::TDistribution_CDF_DF1_IsCauchy", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(1);
		
		// df=1 is Cauchy distribution centered at 0
		// CDF(0) = 0.5, CDF(1) = 0.5 + arctan(1)/Ï€ â‰ˆ 0.75
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
		// 95th percentile (one-tailed) â‰ˆ 1.812
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
		// For t = 2.228 (df=10), two-tailed p â‰ˆ 0.05
		REQUIRE_THAT(t.twoTailedPValue(2.228), WithinAbs(REAL(0.05), REAL(1e-2)));
		
		// Symmetric
		REQUIRE_THAT(t.twoTailedPValue(-2.228), WithinAbs(t.twoTailedPValue(2.228), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_CriticalValue", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(10);
		
		// Critical value for Î± = 0.05 (two-tailed) with df=10 is â‰ˆ 2.228
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
		
		// Var = Î½/(Î½-2) = 10/8 = 1.25 for df=10
		REQUIRE_THAT(t.variance(), WithinAbs(REAL(1.25), REAL(1e-10)));
	}

	TEST_CASE("Statistics::TDistribution_Variance_HighDF_ApproachesOne", "[statistics][distribution]")
	{
		TEST_PRECISION_INFO();
		Statistics::TDistribution t(100);
		
		// As df â†’ âˆž, variance â†’ 1 (approaches standard normal)
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

	/*********************************************************************/
	/*****     NEW DISTRIBUTIONS (Gamma, Beta, Weibull, etc.)        *****/
	/*********************************************************************/

	TEST_CASE("GammaDistribution_BasicProperties", "[statistics][distributions]")
	{
		TEST_PRECISION_INFO();
		
		// Gamma(2, 1) - shape=2, scale=1
		Statistics::GammaDistribution gamma(2.0, 1.0);
		
		// Mean = shape * scale = 2
		REQUIRE_THAT(gamma.mean(), WithinAbs(REAL(2.0), REAL(1e-10)));
		// Variance = shape * scale^2 = 2
		REQUIRE_THAT(gamma.variance(), WithinAbs(REAL(2.0), REAL(1e-10)));
		
		// PDF at x=1: f(1) = 1 * e^(-1) = 0.3679
		REQUIRE_THAT(gamma.pdf(1.0), WithinAbs(REAL(0.3678794), REAL(1e-5)));
		
		// CDF properties
		REQUIRE_THAT(gamma.cdf(0.0), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE(gamma.cdf(10.0) > 0.99);
	}

	TEST_CASE("BetaDistribution_BasicProperties", "[statistics][distributions]")
	{
		TEST_PRECISION_INFO();
		
		// Beta(2, 5)
		Statistics::BetaDistribution beta(2.0, 5.0);
		
		// Mean = α/(α+β) = 2/7
		REQUIRE_THAT(beta.mean(), WithinAbs(REAL(2.0/7.0), REAL(1e-10)));
		
		// CDF boundaries
		REQUIRE_THAT(beta.cdf(0.0), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(beta.cdf(1.0), WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// PDF integrates to 1 (spot check)
		REQUIRE(beta.pdf(0.3) > 0.0);
	}

	TEST_CASE("WeibullDistribution_BasicProperties", "[statistics][distributions]")
	{
		TEST_PRECISION_INFO();
		
		// Weibull(1, 1) = Exponential(1)
		Statistics::WeibullDistribution weibull1(1.0, 1.0);
		REQUIRE_THAT(weibull1.mean(), WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// Weibull(2, 1) - Rayleigh distribution
		Statistics::WeibullDistribution weibull2(2.0, 1.0);
		
		// CDF has closed form: 1 - exp(-(x/λ)^k)
		REQUIRE_THAT(weibull2.cdf(1.0), WithinAbs(REAL(1.0 - std::exp(-1.0)), REAL(1e-10)));
		
		// Inverse CDF round-trip
		Real q = weibull2.inverseCdf(0.5);
		REQUIRE_THAT(weibull2.cdf(q), WithinAbs(REAL(0.5), REAL(1e-8)));
	}

	TEST_CASE("ParetoDistribution_BasicProperties", "[statistics][distributions]")
	{
		TEST_PRECISION_INFO();
		
		// Pareto(3, 1) - 80/20 rule territory
		Statistics::ParetoDistribution pareto(3.0, 1.0);
		
		// Mean = α*xm/(α-1) = 3*1/2 = 1.5
		REQUIRE_THAT(pareto.mean(), WithinAbs(REAL(1.5), REAL(1e-10)));
		
		// CDF at x=1 is 0
		REQUIRE_THAT(pareto.cdf(1.0), WithinAbs(REAL(0.0), REAL(1e-10)));
		
		// CDF at x=2: 1 - (1/2)^3 = 0.875
		REQUIRE_THAT(pareto.cdf(2.0), WithinAbs(REAL(0.875), REAL(1e-10)));
	}

	TEST_CASE("LogNormalDistribution_BasicProperties", "[statistics][distributions]")
	{
		TEST_PRECISION_INFO();
		
		// LogNormal(0, 1) - standard log-normal
		Statistics::LogNormalDistribution lognorm(0.0, 1.0);
		
		// Median = exp(μ) = 1
		REQUIRE_THAT(lognorm.median(), WithinAbs(REAL(1.0), REAL(1e-10)));
		
		// Mean = exp(μ + σ²/2) = exp(0.5)
		REQUIRE_THAT(lognorm.mean(), WithinAbs(std::exp(0.5), REAL(1e-10)));
		
		// CDF at median is 0.5
		REQUIRE_THAT(lognorm.cdf(1.0), WithinAbs(REAL(0.5), REAL(1e-10)));
	}

	TEST_CASE("UniformDistribution_BasicProperties", "[statistics][distributions]")
	{
		TEST_PRECISION_INFO();
		
		// Uniform(0, 10)
		Statistics::UniformDistribution uniform(0.0, 10.0);
		
		// Mean = (a+b)/2 = 5
		REQUIRE_THAT(uniform.mean(), WithinAbs(REAL(5.0), REAL(1e-10)));
		
		// Variance = (b-a)^2/12 = 100/12
		REQUIRE_THAT(uniform.variance(), WithinAbs(REAL(100.0/12.0), REAL(1e-10)));
		
		// PDF = 1/(b-a) = 0.1
		REQUIRE_THAT(uniform.pdf(5.0), WithinAbs(REAL(0.1), REAL(1e-10)));
		
		// CDF at midpoint is 0.5
		REQUIRE_THAT(uniform.cdf(5.0), WithinAbs(REAL(0.5), REAL(1e-10)));
	}

}