#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Statistics/DiscreteDistributions.h>
#endif
#include <catch2/catch_all.hpp>

using namespace MML;
using namespace MML::Statistics;
using Catch::Matchers::WithinRel;

namespace {
	template<class Distribution>
	void CheckQuantile(const Distribution& distribution, Real probability) {
		const int k = distribution.inverseCdf(probability);
		CAPTURE(probability, k);
		REQUIRE(k >= 0);
		CHECK(distribution.cdf(k) >= probability);
		if (k > 0)
			CHECK(distribution.cdf(k - 1) < probability);
	}
}

TEST_CASE("Discrete quantiles repair the recorded support and degenerate endpoints", "[statistics][discrete-boundaries]") {
	CHECK(GeometricDistribution(Real(0.5)).inverseCdf(Real(1e-20)) == 0);
	CHECK(GeometricDistribution(1).inverseCdf(1) == 0);
	CHECK(NegativeBinomialDistribution(3, 1).inverseCdf(1) == 0);
	CHECK(BinomialDistribution(5, 0).inverseCdf(1) == 0);
	CHECK(BernoulliDistribution(Real(1e-20)).inverseCdf(1) == 1);
}

TEST_CASE("Discrete quantile NaN cannot silently choose an ordinary endpoint", "[statistics][discrete-boundaries]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS(BernoulliDistribution(Real(0.3)).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(BinomialDistribution(10, Real(0.3)).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(PoissonDistribution(3).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(GeometricDistribution(Real(0.5)).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(NegativeBinomialDistribution(3, Real(0.5)).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(HypergeometricDistribution(50, 10, 5).inverseCdf(nan), StatisticsError);
}

TEST_CASE("Geometric probabilities retain tiny success probabilities", "[statistics][discrete-boundaries]") {
	const Real p = Real(1e-20);
	CHECK_THAT(GeometricDistribution(p).cdf(0), WithinRel(p, Precision::EqualityTolerance));
	CHECK_THAT(GeometricDistribution(p).cdf(9), WithinRel(10 * p, Precision::EqualityTolerance));
	CHECK_THAT(NegativeBinomialDistribution(1, p).cdf(9), WithinRel(10 * p, Precision::EqualityTolerance));
	CHECK(GeometricDistribution(Real(0.5)).inverseCdf(std::numeric_limits<Real>::denorm_min()) == 0);
}

TEST_CASE("Discrete inverse probability clipping and unbounded endpoint sentinels remain intentional", "[statistics][discrete-boundaries]") {
	for (Real probability : {Real(-1), Real(0), -std::numeric_limits<Real>::infinity()}) {
		CHECK(BernoulliDistribution(1).inverseCdf(probability) == 0);
		CHECK(BinomialDistribution(5, 1).inverseCdf(probability) == 0);
		CHECK(PoissonDistribution(3).inverseCdf(probability) == 0);
		CHECK(GeometricDistribution(Real(0.5)).inverseCdf(probability) == 0);
		CHECK(NegativeBinomialDistribution(3, Real(0.5)).inverseCdf(probability) == 0);
		CHECK(HypergeometricDistribution(10, 8, 5).inverseCdf(probability) == 3);
	}
	for (Real probability : {Real(1), Real(2), std::numeric_limits<Real>::infinity()}) {
		CHECK(BernoulliDistribution(0).inverseCdf(probability) == 0);
		CHECK(BernoulliDistribution(1).inverseCdf(probability) == 1);
		CHECK(BinomialDistribution(5, 0).inverseCdf(probability) == 0);
		CHECK(BinomialDistribution(5, 1).inverseCdf(probability) == 5);
		CHECK(PoissonDistribution(3).inverseCdf(probability) == std::numeric_limits<int>::max());
		CHECK(GeometricDistribution(Real(0.5)).inverseCdf(probability) == std::numeric_limits<int>::max());
		CHECK(NegativeBinomialDistribution(3, Real(0.5)).inverseCdf(probability) == std::numeric_limits<int>::max());
		CHECK(GeometricDistribution(1).inverseCdf(probability) == 0);
		CHECK(NegativeBinomialDistribution(3, 1).inverseCdf(probability) == 0);
		CHECK(HypergeometricDistribution(10, 8, 5).inverseCdf(probability) == 5);
	}
}

TEST_CASE("Finite discrete CDF accumulation remains a bounded probability", "[statistics][discrete-boundaries]") {
	for (int n : {100, 500, 1000}) {
		Real binomial = BinomialDistribution(n, Real(0.01)).cdf(n - 1);
		Real hypergeometric = HypergeometricDistribution(10 * n, n, n).cdf(n - 1);
		CAPTURE(n, binomial, hypergeometric);
		CHECK(binomial >= 0);
		CHECK(binomial <= 1);
		CHECK(hypergeometric >= 0);
		CHECK(hypergeometric <= 1);
	}
}

TEST_CASE("Discrete quantiles certify the first integer across ordinary distributions", "[statistics][discrete-boundaries]") {
	for (Real probability : {Real(1e-20), Real(0.01), Real(0.1), Real(0.5), Real(0.9), Real(0.99)}) {
		CheckQuantile(BernoulliDistribution(Real(0.3)), probability);
		CheckQuantile(BinomialDistribution(20, Real(0.4)), probability);
		CheckQuantile(PoissonDistribution(Real(7.5)), probability);
		CheckQuantile(GeometricDistribution(Real(0.2)), probability);
		CheckQuantile(NegativeBinomialDistribution(4, Real(0.4)), probability);
		CheckQuantile(HypergeometricDistribution(10, 8, 5), probability);
	}
}

TEST_CASE("Discrete quantile jump neighbors retain their first-integer contract", "[statistics][discrete-boundaries]") {
	for (int k = 0; k < 12; ++k) {
		auto checkJump = [k](const auto& distribution) {
			const Real probability = distribution.cdf(k);
			if (probability <= 0 || probability >= 1) return;
			CAPTURE(k, probability);
			CheckQuantile(distribution, probability);
			CheckQuantile(distribution, std::nextafter(probability, Real(0)));
			CheckQuantile(distribution, std::nextafter(probability, Real(1)));
		};
		checkJump(BinomialDistribution(20, Real(0.4)));
		checkJump(PoissonDistribution(Real(7.5)));
		checkJump(GeometricDistribution(Real(0.5)));
		checkJump(NegativeBinomialDistribution(3, Real(0.5)));
	}
}

TEST_CASE("Interior integer-range failures never become unbounded endpoint sentinels", "[statistics][discrete-boundaries]") {
	CHECK_THROWS_AS(PoissonDistribution(2 * Real(std::numeric_limits<int>::max())).inverseCdf(Real(0.5)), ConvergenceError);
	CHECK_THROWS_AS(PoissonDistribution(std::numeric_limits<Real>::max() / 2).inverseCdf(Real(0.5)), ConvergenceError);
	CHECK_THROWS_AS(GeometricDistribution(Real(1e-20)).inverseCdf(Real(0.5)), ConvergenceError);
	CHECK_THROWS_AS(NegativeBinomialDistribution(3, Real(1e-20)).inverseCdf(Real(0.5)), ConvergenceError);
	const PoissonDistribution large(Real(1e6));
	const int median = large.inverseCdf(Real(0.5));
	CHECK(median == 1000000);
	CHECK(large.cdf(median - 1) < Real(0.5));
	CHECK(large.cdf(median) >= Real(0.5));
}

TEST_CASE("Discrete PMFs and CDFs form widened arguments before integer-boundary additions", "[statistics][discrete-boundaries]") {
	const int maximum = std::numeric_limits<int>::max();
	const Real p = 1 / (Real(maximum) + 1);
	const Real zeroMass = std::exp(Real(maximum) * std::log1p(-p));
	CHECK_THAT(BinomialDistribution(maximum, p).pmf(0), WithinRel(zeroMass, Precision::EqualityTolerance));
	CHECK_THAT(BinomialDistribution(maximum, p).pmf(1),
		WithinRel(Real(maximum) * p * std::exp(Real(maximum - 1) * std::log1p(-p)), Precision::EqualityTolerance));
	CHECK_THAT(BinomialDistribution(maximum, p).cdf(0), WithinRel(zeroMass, Precision::EqualityTolerance));
	CHECK(BinomialDistribution(maximum, 0).inverseCdf(Real(0.5)) == 0);
	CHECK(BinomialDistribution(maximum, 1).inverseCdf(Real(0.5)) == maximum);
	CHECK(PoissonDistribution(1).pmf(maximum) == 0);
	CHECK(PoissonDistribution(1).cdf(maximum) == 1);
	Real negativeMass = NegativeBinomialDistribution(maximum, Real(0.5)).pmf(maximum);
	CHECK(std::isfinite(negativeMass));
	CHECK(negativeMass > 0);
	CHECK(negativeMass <= 1);
	CHECK_THAT(GeometricDistribution(p).cdf(maximum),
		WithinRel(-std::expm1((Real(maximum) + 1) * std::log1p(-p)), Precision::EqualityTolerance));
	CHECK_THAT(NegativeBinomialDistribution(1, p).cdf(maximum),
		WithinRel(GeometricDistribution(p).cdf(maximum), Precision::EqualityTolerance));
	CheckQuantile(GeometricDistribution(p), Real(0.5));
	CheckQuantile(NegativeBinomialDistribution(1, p), Real(0.5));
	CHECK_THAT(HypergeometricDistribution(maximum, 1, 1).pmf(1), WithinRel(1 / Real(maximum), Precision::EqualityTolerance));
	CHECK(HypergeometricDistribution(maximum, maximum, maximum).pmf(maximum) == 1);
	CHECK(HypergeometricDistribution(maximum, maximum, maximum).cdf(maximum - 1) == 0);
	CHECK(HypergeometricDistribution(maximum, maximum, maximum).inverseCdf(Real(0.5)) == maximum);
}

TEST_CASE("Negative binomial masses resolve large integer modes without gamma cancellation",
	"[statistics][discrete-boundaries][float-contract]") {
	for (int count : {1000, 1000000, std::numeric_limits<int>::max()}) {
		CAPTURE(count);
		const Real n = Real(count);
		const Real inverse = 1 / n;
		// Central-binomial Stirling expansion, independent of the Poisson evaluation.
		const Real correction = 1 - inverse / 8 + inverse * inverse / 128 +
			Real(5) * inverse * inverse * inverse / 1024 -
			Real(21) * inverse * inverse * inverse * inverse / 32768;
		const Real expected = correction / (2 * std::sqrt(Constants::PI * n));
		NegativeBinomialDistribution distribution(count, Real(0.5));
		const Real mass = distribution.pmf(count);
		REQUIRE(std::isfinite(mass));
		REQUIRE(mass > 0);
		CHECK(IsWithinRel(mass, expected, Precision::EqualityTolerance));
		const Real previous = distribution.pmf(count - 1);
		CHECK(IsWithinRel(mass / previous, (n + Real(count - 1)) / n / 2, Precision::EqualityTolerance));
	}
	for (Real probability : {Real(0.125), Real(0.5), Real(0.875)}) {
		const int successes = 1000000;
		const int failures = int(Real(successes) * (1 - probability) / probability);
		NegativeBinomialDistribution distribution(successes, probability);
		const Real mass = distribution.pmf(failures), next = distribution.pmf(failures + 1);
		REQUIRE(std::isfinite(mass));
		REQUIRE(mass > 0);
		CHECK(IsWithinRel(next / mass, (Real(failures) + successes) / Real(failures + 1) * (1 - probability),
			Precision::EqualityTolerance));
	}
	CHECK(IsWithinRel(NegativeBinomialDistribution(2, Real(0.25)).pmf(3), Real(0.10546875), Precision::EqualityTolerance));
	CHECK(NegativeBinomialDistribution(2, Real(0.25)).pmf(-1) == 0);
	CHECK(NegativeBinomialDistribution(2, 1).pmf(0) == 1);
	CHECK(NegativeBinomialDistribution(2, 1).pmf(1) == 0);
}

TEST_CASE("Discrete search and accumulation budgets report exact failures", "[statistics][discrete-boundaries]") {
	try {
		StatisticsDetail::DiscreteQuantile(Real(0.5), [](int k) { return Real(k) / 100; }, 0, 100, 0);
		FAIL("Expected exhausted integer-search budget");
	} catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 0);
		CHECK(error.residual() == 100);
	}
	try {
		StatisticsDetail::SumDiscreteProbabilities(0, Precision::DiscreteCdfMaxTerms, [](int) { return Real(0); });
		FAIL("Expected accumulation budget failure");
	} catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 0);
		CHECK(error.residual() == Real(Precision::DiscreteCdfMaxTerms) + 1);
	}
	const int maximum = std::numeric_limits<int>::max();
	CHECK_THAT(StatisticsDetail::SumDiscreteProbabilities(maximum - 2, maximum, [](int) { return Real(0.1); }),
		WithinRel(Real(0.3), Precision::EqualityTolerance));
	if constexpr (std::numeric_limits<Real>::digits >= std::numeric_limits<int>::digits)
		CHECK(StatisticsDetail::DiscreteQuantile(1, [](int k) { return Real(k) / Real(maximum); }, 0, maximum) == maximum);
	CHECK_THROWS_AS(StatisticsDetail::DiscreteQuantile(Real(0.5),
		[](int) { return std::numeric_limits<Real>::quiet_NaN(); }, 0, 1), NumericalMethodError);
	CHECK_THROWS_AS(HypergeometricDistribution(2000000, 1000000, 1000000).cdf(500000), ConvergenceError);
}

TEST_CASE("Discrete inverse entry points revalidate public parameter state", "[statistics][discrete-boundaries]") {
	for (Real invalid : {Real(-1), Real(2), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		BernoulliDistribution bernoulli(Real(0.5));
		BinomialDistribution binomial(10, Real(0.5));
		GeometricDistribution geometric(Real(0.5));
		NegativeBinomialDistribution negative(3, Real(0.5));
		bernoulli.p = binomial.p = geometric.p = negative.p = invalid;
		CHECK_THROWS_AS(bernoulli.inverseCdf(Real(0.5)), StatisticsError);
		CHECK_THROWS_AS(binomial.inverseCdf(Real(0.5)), StatisticsError);
		CHECK_THROWS_AS(geometric.inverseCdf(Real(0.5)), StatisticsError);
		CHECK_THROWS_AS(negative.inverseCdf(Real(0.5)), StatisticsError);
	}
	for (Real invalid : {Real(-1), Real(0), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		PoissonDistribution poisson(3);
		poisson.lambda = invalid;
		CHECK_THROWS_AS(poisson.inverseCdf(Real(0.5)), StatisticsError);
	}
	BinomialDistribution binomial(10, Real(0.5));
	NegativeBinomialDistribution negative(3, Real(0.5));
	HypergeometricDistribution hypergeometric(50, 10, 5);
	binomial.n = -1;
	negative.r = 0;
	hypergeometric.K = 51;
	CHECK_THROWS_AS(binomial.inverseCdf(Real(0.5)), StatisticsError);
	CHECK_THROWS_AS(negative.inverseCdf(Real(0.5)), StatisticsError);
	CHECK_THROWS_AS(hypergeometric.inverseCdf(Real(0.5)), StatisticsError);
}

TEST_CASE("Discrete integer quantiles match independent high-precision recurrence references", "[statistics][discrete-boundaries]") {
	// Local mpmath PMF recurrences at 100 decimal digits, independent of the beta/gamma CDF identities.
	CHECK(PoissonDistribution(40).inverseCdf(Real(1e-15)) == 2);
	CHECK(PoissonDistribution(40).inverseCdf(Real(0.99)) == 55);
	if (Real(1e-100) > 0) {
		CHECK(PoissonDistribution(1000).inverseCdf(Real(1e-100)) == 408);
		CheckQuantile(PoissonDistribution(1000), Real(1e-100));
	}
	CHECK(PoissonDistribution(1000).inverseCdf(Real(0.9)) == 1041);
	CHECK(BinomialDistribution(1000, Real(0.001)).inverseCdf(Real(0.99)) == 4);
	CHECK(NegativeBinomialDistribution(5, Real(0.2)).inverseCdf(Real(0.99)) == 50);
	CHECK(HypergeometricDistribution(1000, 100, 100).inverseCdf(Real(0.99)) == 17);
	CheckQuantile(PoissonDistribution(40), Real(1e-15));
	CheckQuantile(PoissonDistribution(40), Real(0.99));
	CheckQuantile(PoissonDistribution(1000), Real(0.9));
	CheckQuantile(BinomialDistribution(1000, Real(0.001)), Real(0.99));
	CheckQuantile(NegativeBinomialDistribution(5, Real(0.2)), Real(0.99));
	CheckQuantile(HypergeometricDistribution(1000, 100, 100), Real(0.99));
}
