#include <catch2/catch_all.hpp>
#include "TestPrecision.h"
#include "TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/Statistics/DiscreteDistributions.h>
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::DiscreteDistributionsTests
{
	TEST_CASE("BernoulliDistribution::PMF_CDF_Moments", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BernoulliDistribution bernoulli(static_cast<MML::Real>(0.3));

		REQUIRE_THAT(bernoulli.pmf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(bernoulli.pmf(0), WithinAbs(static_cast<MML::Real>(0.7), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(bernoulli.pmf(1), WithinAbs(static_cast<MML::Real>(0.3), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(bernoulli.pmf(2), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));

		REQUIRE_THAT(bernoulli.cdf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(bernoulli.cdf(0), WithinAbs(static_cast<MML::Real>(0.7), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(bernoulli.cdf(1), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));

		REQUIRE_THAT(bernoulli.mean(), WithinAbs(static_cast<MML::Real>(0.3), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(bernoulli.variance(), WithinAbs(static_cast<MML::Real>(0.21), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("BernoulliDistribution::InverseCDF", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BernoulliDistribution bernoulli(static_cast<MML::Real>(0.3));

		REQUIRE(bernoulli.inverseCdf(static_cast<MML::Real>(0.0)) == 0);
		REQUIRE(bernoulli.inverseCdf(static_cast<MML::Real>(0.69)) == 0);
		REQUIRE(bernoulli.inverseCdf(static_cast<MML::Real>(0.7)) == 0);
		REQUIRE(bernoulli.inverseCdf(static_cast<MML::Real>(0.71)) == 1);
		REQUIRE(bernoulli.inverseCdf(static_cast<MML::Real>(1.0)) == 1);
	}

	TEST_CASE("BernoulliDistribution::EdgeCases", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BernoulliDistribution alwaysFailure(static_cast<MML::Real>(0.0));
		REQUIRE_THAT(alwaysFailure.pmf(0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(alwaysFailure.pmf(1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE(alwaysFailure.inverseCdf(static_cast<MML::Real>(1.0)) == 0);

		Statistics::BernoulliDistribution alwaysSuccess(static_cast<MML::Real>(1.0));
		REQUIRE_THAT(alwaysSuccess.pmf(0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(alwaysSuccess.pmf(1), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE(alwaysSuccess.inverseCdf(static_cast<MML::Real>(0.0)) == 0);
		REQUIRE(alwaysSuccess.inverseCdf(static_cast<MML::Real>(0.01)) == 1);

		REQUIRE_THROWS_AS(Statistics::BernoulliDistribution(-static_cast<MML::Real>(0.1)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::BernoulliDistribution(static_cast<MML::Real>(1.1)), StatisticsError);
	}

	TEST_CASE("BinomialDistribution::PMF_CDF_Moments", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BinomialDistribution binomial(10, static_cast<MML::Real>(0.5));

		REQUIRE_THAT(binomial.mean(), WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(binomial.variance(), WithinAbs(static_cast<MML::Real>(2.5), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(binomial.stddev(), WithinAbs(std::sqrt(static_cast<MML::Real>(2.5)), TOL(1e-12, 1e-6)));

		REQUIRE(binomial.pmf(5) > binomial.pmf(3));
		REQUIRE(binomial.pmf(5) > binomial.pmf(7));

		Real pmfSum = 0.0;
		for (int k = 0; k <= 10; ++k)
			pmfSum += binomial.pmf(k);
		REQUIRE_THAT(pmfSum, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));

		REQUIRE_THAT(binomial.cdf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(binomial.cdf(10), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("BinomialDistribution::KnownValues", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BinomialDistribution binomial(10, static_cast<MML::Real>(0.3));

		REQUIRE_THAT(binomial.pmf(3), WithinAbs(static_cast<MML::Real>(0.266827932), TOL(1e-9, 1e-5)));
		REQUIRE_THAT(binomial.pmf(0), WithinAbs(static_cast<MML::Real>(0.0282475249), TOL(1e-9, 1e-5)));
		REQUIRE_THAT(binomial.cdf(2), WithinAbs(static_cast<MML::Real>(0.3827827864), TOL(1e-9, 1e-5)));
	}

	TEST_CASE("BinomialDistribution::InverseCDF", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BinomialDistribution binomial(20, static_cast<MML::Real>(0.4));

		for (Real probability : { static_cast<MML::Real>(0.1), static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.75), static_cast<MML::Real>(0.9) }) {
			int k = binomial.inverseCdf(probability);
			REQUIRE(binomial.cdf(k) >= probability - TOL(1e-12, 1e-6));
			if (k > 0)
				REQUIRE(binomial.cdf(k - 1) < probability + TOL(1e-12, 1e-6));
		}
	}

	TEST_CASE("BinomialDistribution::EdgeCases", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::BinomialDistribution noSuccesses(5, static_cast<MML::Real>(0.0));
		REQUIRE_THAT(noSuccesses.pmf(0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(noSuccesses.pmf(1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE(noSuccesses.inverseCdf(static_cast<MML::Real>(1.0)) == 5);

		Statistics::BinomialDistribution allSuccesses(5, static_cast<MML::Real>(1.0));
		REQUIRE_THAT(allSuccesses.pmf(4), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(allSuccesses.pmf(5), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE(allSuccesses.inverseCdf(static_cast<MML::Real>(0.5)) == 5);

		REQUIRE_THROWS_AS(Statistics::BinomialDistribution(-1, static_cast<MML::Real>(0.5)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::BinomialDistribution(5, -static_cast<MML::Real>(0.1)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::BinomialDistribution(5, static_cast<MML::Real>(1.1)), StatisticsError);
	}

	TEST_CASE("PoissonDistribution::PMF_CDF_Moments", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::PoissonDistribution poisson(static_cast<MML::Real>(5.0));

		REQUIRE_THAT(poisson.mean(), WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(poisson.variance(), WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(poisson.stddev(), WithinAbs(std::sqrt(static_cast<MML::Real>(5.0)), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(poisson.pmf(5), WithinAbs(static_cast<MML::Real>(0.1754673698), TOL(1e-9, 1e-5)));
		REQUIRE_THAT(poisson.pmf(0), WithinAbs(static_cast<MML::Real>(0.0067379470), TOL(1e-9, 1e-5)));
		REQUIRE_THAT(poisson.pmf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("PoissonDistribution::CDF_Monotonicity", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::PoissonDistribution poisson(static_cast<MML::Real>(3.0));

		Real pmfSum = 0.0;
		for (int k = 0; k <= 50; ++k)
			pmfSum += poisson.pmf(k);
		REQUIRE_THAT(pmfSum, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));

		Real previousCdf = poisson.cdf(-1);
		for (int k = 0; k <= 20; ++k) {
			Real currentCdf = poisson.cdf(k);
			REQUIRE(currentCdf >= previousCdf);
			previousCdf = currentCdf;
		}

		REQUIRE_THAT(poisson.cdf(2), WithinAbs(static_cast<MML::Real>(0.4231900811), TOL(1e-9, 1e-5)));
		REQUIRE(poisson.cdf(20) > static_cast<MML::Real>(0.999999));
	}

	TEST_CASE("PoissonDistribution::InverseCDF", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::PoissonDistribution poisson(static_cast<MML::Real>(7.5));

		for (Real probability : { static_cast<MML::Real>(0.1), static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.75), static_cast<MML::Real>(0.9) }) {
			int k = poisson.inverseCdf(probability);
			REQUIRE(poisson.cdf(k) >= probability - TOL(1e-12, 1e-6));
			if (k > 0)
				REQUIRE(poisson.cdf(k - 1) < probability + TOL(1e-12, 1e-6));
		}
	}

	TEST_CASE("PoissonDistribution::EdgeCases", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::PoissonDistribution poisson(static_cast<MML::Real>(0.5));
		REQUIRE(poisson.inverseCdf(static_cast<MML::Real>(0.0)) == 0);
		REQUIRE(poisson.inverseCdf(static_cast<MML::Real>(1.0)) == std::numeric_limits<int>::max());

		REQUIRE_THROWS_AS(Statistics::PoissonDistribution(static_cast<MML::Real>(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::PoissonDistribution(-static_cast<MML::Real>(1.0)), StatisticsError);
	}

	TEST_CASE("GeometricDistribution::PMF_CDF_Moments", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::GeometricDistribution geometric(static_cast<MML::Real>(0.2));

		REQUIRE_THAT(geometric.mean(), WithinAbs(static_cast<MML::Real>(4.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(geometric.variance(), WithinAbs(static_cast<MML::Real>(20.0), TOL(1e-12, 3e-6)));
		REQUIRE_THAT(geometric.pmf(0), WithinAbs(static_cast<MML::Real>(0.2), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(geometric.pmf(1), WithinAbs(static_cast<MML::Real>(0.16), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(geometric.pmf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(geometric.cdf(0), WithinAbs(static_cast<MML::Real>(0.2), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(geometric.cdf(1), WithinAbs(static_cast<MML::Real>(0.36), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("GeometricDistribution::InverseCDF_And_EdgeCases", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::GeometricDistribution geometric(static_cast<MML::Real>(0.2));
		for (Real probability : { static_cast<MML::Real>(0.1), static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.75), static_cast<MML::Real>(0.9) }) {
			int k = geometric.inverseCdf(probability);
			REQUIRE(geometric.cdf(k) >= probability - TOL(1e-12, 1e-6));
			if (k > 0)
				REQUIRE(geometric.cdf(k - 1) < probability + TOL(1e-12, 1e-6));
		}

		Statistics::GeometricDistribution certainSuccess(static_cast<MML::Real>(1.0));
		REQUIRE_THAT(certainSuccess.pmf(0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(certainSuccess.pmf(1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE(certainSuccess.inverseCdf(static_cast<MML::Real>(0.5)) == 0);

		REQUIRE_THROWS_AS(Statistics::GeometricDistribution(static_cast<MML::Real>(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::GeometricDistribution(-static_cast<MML::Real>(0.1)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::GeometricDistribution(static_cast<MML::Real>(1.1)), StatisticsError);
	}

	TEST_CASE("NegativeBinomialDistribution::PMF_CDF_Moments", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::NegativeBinomialDistribution negativeBinomial(3, static_cast<MML::Real>(0.5));

		REQUIRE_THAT(negativeBinomial.mean(), WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(negativeBinomial.variance(), WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(negativeBinomial.pmf(0), WithinAbs(static_cast<MML::Real>(0.125), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(negativeBinomial.pmf(1), WithinAbs(static_cast<MML::Real>(0.1875), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(negativeBinomial.cdf(1), WithinAbs(static_cast<MML::Real>(0.3125), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("NegativeBinomialDistribution::ReducesToGeometric", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Real probability = static_cast<MML::Real>(0.3);
		Statistics::NegativeBinomialDistribution negativeBinomial(1, probability);
		Statistics::GeometricDistribution geometric(probability);

		for (int k = 0; k <= 10; ++k)
			REQUIRE_THAT(negativeBinomial.pmf(k), WithinAbs(geometric.pmf(k), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("NegativeBinomialDistribution::InverseCDF_And_EdgeCases", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::NegativeBinomialDistribution negativeBinomial(4, static_cast<MML::Real>(0.4));
		for (Real probability : { static_cast<MML::Real>(0.1), static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.75), static_cast<MML::Real>(0.9) }) {
			int k = negativeBinomial.inverseCdf(probability);
			REQUIRE(negativeBinomial.cdf(k) >= probability - TOL(1e-12, 1e-6));
			if (k > 0)
				REQUIRE(negativeBinomial.cdf(k - 1) < probability + TOL(1e-12, 1e-6));
		}

		Statistics::NegativeBinomialDistribution certainSuccess(2, static_cast<MML::Real>(1.0));
		REQUIRE_THAT(certainSuccess.pmf(0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(certainSuccess.pmf(1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE(certainSuccess.inverseCdf(static_cast<MML::Real>(0.5)) == 0);

		REQUIRE_THROWS_AS(Statistics::NegativeBinomialDistribution(0, static_cast<MML::Real>(0.5)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::NegativeBinomialDistribution(2, static_cast<MML::Real>(0.0)), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::NegativeBinomialDistribution(2, static_cast<MML::Real>(1.1)), StatisticsError);
	}

	TEST_CASE("HypergeometricDistribution::PMF_CDF_Moments", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::HypergeometricDistribution hypergeometric(50, 10, 5);

		REQUIRE_THAT(hypergeometric.mean(), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(hypergeometric.variance(), WithinAbs(static_cast<MML::Real>(0.7346938776), TOL(1e-9, 1e-5)));

		Real pmfSum = 0.0;
		for (int k = 0; k <= 5; ++k)
			pmfSum += hypergeometric.pmf(k);
		REQUIRE_THAT(pmfSum, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));

		REQUIRE_THAT(hypergeometric.cdf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(hypergeometric.cdf(5), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("HypergeometricDistribution::KnownValues", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::HypergeometricDistribution hypergeometric(52, 13, 5);

		REQUIRE_THAT(hypergeometric.pmf(2), WithinAbs(static_cast<MML::Real>(0.2742797119), TOL(1e-9, 1e-5)));
		REQUIRE_THAT(hypergeometric.pmf(-1), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(hypergeometric.pmf(6), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
	}

	TEST_CASE("HypergeometricDistribution::SupportAndInverseCDF", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::HypergeometricDistribution constrained(10, 8, 5);
		REQUIRE_THAT(constrained.pmf(2), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE(constrained.inverseCdf(static_cast<MML::Real>(0.0)) == 3);
		REQUIRE(constrained.inverseCdf(static_cast<MML::Real>(1.0)) == 5);

		for (Real probability : { static_cast<MML::Real>(0.1), static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.75), static_cast<MML::Real>(0.9) }) {
			int k = constrained.inverseCdf(probability);
			REQUIRE(constrained.cdf(k) >= probability - TOL(1e-12, 1e-6));
			if (k > 3)
				REQUIRE(constrained.cdf(k - 1) < probability + TOL(1e-12, 1e-6));
		}
	}

	TEST_CASE("HypergeometricDistribution::EdgeCases", "[statistics][discrete]")
	{
		TEST_PRECISION_INFO();

		Statistics::HypergeometricDistribution emptyPopulation(0, 0, 0);
		REQUIRE_THAT(emptyPopulation.pmf(0), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(emptyPopulation.mean(), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(emptyPopulation.variance(), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-12, 1e-6)));

		REQUIRE_THROWS_AS(Statistics::HypergeometricDistribution(-1, 0, 0), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::HypergeometricDistribution(10, -1, 5), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::HypergeometricDistribution(10, 11, 5), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::HypergeometricDistribution(10, 5, -1), StatisticsError);
		REQUIRE_THROWS_AS(Statistics::HypergeometricDistribution(10, 5, 11), StatisticsError);
	}
} // namespace MML::Tests::Algorithms::DiscreteDistributionsTests
