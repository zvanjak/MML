#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Statistics/Distributions.h>
#include <mml/algorithms/Statistics/DiscreteDistributions.h>
#endif
#include <catch2/catch_all.hpp>
#include <limits>
#include <random>
#include <stdexcept>

using namespace MML;
using namespace MML::Statistics;

namespace {
	struct RejectingUrbg {
		using result_type = unsigned int;
		static constexpr result_type min() { return 0; }
		static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }
		int calls = 0;
		result_type operator()() {
			++calls;
			throw std::runtime_error("Unexpected RNG access in invalid-parameter test");
		}
	};

	template<class Distribution, class Value>
	void CheckInvalidState(Distribution distribution, Value Distribution::* field, Value invalid) {
		distribution.*field = invalid;
		RejectingUrbg gen;
		REQUIRE_THROWS_AS(distribution.sample(gen), StatisticsError);
		CHECK(gen.calls == 0);
		const auto before = Random::Engine();
		CHECK_THROWS_AS(distribution.sample(), StatisticsError);
		CHECK(Random::Engine() == before);
		if constexpr (requires { distribution.pmf(0); }) {
			CHECK_THROWS_AS(distribution.pmf(-1), StatisticsError);
			CHECK_THROWS_AS(distribution.pmf(0), StatisticsError);
			CHECK_THROWS_AS(distribution.cdf(-1), StatisticsError);
			CHECK_THROWS_AS(distribution.cdf(std::numeric_limits<int>::max()), StatisticsError);
		} else {
			CHECK_THROWS_AS(distribution.pdf(-1), StatisticsError);
			CHECK_THROWS_AS(distribution.pdf(1), StatisticsError);
			CHECK_THROWS_AS(distribution.cdf(-1), StatisticsError);
			CHECK_THROWS_AS(distribution.cdf(1), StatisticsError);
		}
		for (Real probability : {-std::numeric_limits<Real>::infinity(), Real(0), Real(0.5), Real(1),
			std::numeric_limits<Real>::infinity()})
			CHECK_THROWS_AS(distribution.inverseCdf(probability), StatisticsError);
		if constexpr (requires { distribution.mean(); })
			CHECK_THROWS_AS(distribution.mean(), StatisticsError);
		if constexpr (requires { distribution.variance(); })
			CHECK_THROWS_AS(distribution.variance(), StatisticsError);
		if constexpr (requires { distribution.stddev(); })
			CHECK_THROWS_AS(distribution.stddev(), StatisticsError);
		if constexpr (requires { distribution.median(); })
			CHECK_THROWS_AS(distribution.median(), StatisticsError);
		if constexpr (requires { distribution.zScore(1); }) {
			CHECK_THROWS_AS(distribution.zScore(1), StatisticsError);
			CHECK_THROWS_AS(distribution.fromZScore(1), StatisticsError);
		}
		if constexpr (requires { distribution.sf(1); }) {
			CHECK_THROWS_AS(distribution.sf(1), StatisticsError);
			CHECK_THROWS_AS(distribution.twoTailedPValue(1), StatisticsError);
			CHECK_THROWS_AS(distribution.criticalValue(Real(0.05)), StatisticsError);
		}
	}

	template<class Distribution>
	void CheckSamplingEngineAgreement(const Distribution& distribution) {
		const auto saved = Random::Engine();
		std::mt19937 explicitEngine(17019);
		Random::SetSeed(17019);
		CHECK(distribution.sample(explicitEngine) == distribution.sample());
		CHECK(explicitEngine == Random::Engine());
		Random::Engine() = saved;
	}
}

TEST_CASE("Location and scale distribution constructors reject nonfinite parameters", "[statistics][distribution-parameters]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(NormalDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(NormalDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(CauchyDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(CauchyDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(LogisticDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(LogisticDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(LogNormalDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(LogNormalDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(UniformDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(UniformDistribution(0, invalid), StatisticsError);
	}
}

TEST_CASE("Shape scale and rate constructors reject nonfinite parameters", "[statistics][distribution-parameters]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(GammaDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(GammaDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(BetaDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(BetaDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(WeibullDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(WeibullDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(ParetoDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(ParetoDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(ExponentialDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(PoissonDistribution(invalid), StatisticsError);
	}
}

TEST_CASE("Discrete probability constructors reject nonfinite success probabilities", "[statistics][distribution-parameters]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(BernoulliDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(BinomialDistribution(5, invalid), StatisticsError);
		CHECK_THROWS_AS(GeometricDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(NegativeBinomialDistribution(3, invalid), StatisticsError);
	}
}

TEST_CASE("Direct gamma sampling rejects invalid parameters before touching the URBG", "[statistics][distribution-parameters]") {
	for (Real invalid : {Real(-1), Real(0), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		RejectingUrbg gen;
		CHECK_THROWS_AS(SampleGamma(gen, invalid, 1), StatisticsError);
		CHECK(gen.calls == 0);
		CHECK_THROWS_AS(SampleGamma(gen, 1, invalid), StatisticsError);
		CHECK(gen.calls == 0);
	}
}

TEST_CASE("Closed-form continuous inverse probability guards reject NaN explicitly", "[statistics][distribution-parameters]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS(NormalDistribution().inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(CauchyDistribution().inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(ExponentialDistribution(1).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(LogisticDistribution().inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(UniformDistribution().inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(WeibullDistribution(1).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(ParetoDistribution(1).inverseCdf(nan), StatisticsError);
	CHECK_THROWS_AS(LogNormalDistribution().inverseCdf(nan), StatisticsError);
}

TEST_CASE("Mutable location and interval parameters are checked by every operation", "[statistics][distribution-parameters]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		CheckInvalidState(NormalDistribution(), &NormalDistribution::mu, invalid);
		CheckInvalidState(CauchyDistribution(), &CauchyDistribution::mu, invalid);
		CheckInvalidState(LogisticDistribution(), &LogisticDistribution::mu, invalid);
		CheckInvalidState(LogNormalDistribution(), &LogNormalDistribution::mu, invalid);
		CheckInvalidState(UniformDistribution(), &UniformDistribution::a, invalid);
		CheckInvalidState(UniformDistribution(), &UniformDistribution::b, invalid);
	}
	CheckInvalidState(UniformDistribution(), &UniformDistribution::a, Real(1));
	CheckInvalidState(UniformDistribution(), &UniformDistribution::a, Real(2));
	CheckInvalidState(UniformDistribution(), &UniformDistribution::b, Real(0));
	CheckInvalidState(UniformDistribution(), &UniformDistribution::b, Real(-1));
}

TEST_CASE("Mutable positive parameters are checked by every operation", "[statistics][distribution-parameters]") {
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CheckInvalidState(NormalDistribution(), &NormalDistribution::sigma, invalid);
		CheckInvalidState(CauchyDistribution(), &CauchyDistribution::sigma, invalid);
		CheckInvalidState(LogisticDistribution(), &LogisticDistribution::sigma, invalid);
		CheckInvalidState(LogNormalDistribution(), &LogNormalDistribution::sigma, invalid);
		CheckInvalidState(ExponentialDistribution(1), &ExponentialDistribution::lambda, invalid);
		CheckInvalidState(PoissonDistribution(1), &PoissonDistribution::lambda, invalid);
		CheckInvalidState(GammaDistribution(2, 3), &GammaDistribution::shape, invalid);
		CheckInvalidState(GammaDistribution(2, 3), &GammaDistribution::scale, invalid);
		CheckInvalidState(BetaDistribution(2, 3), &BetaDistribution::alpha, invalid);
		CheckInvalidState(BetaDistribution(2, 3), &BetaDistribution::beta, invalid);
		CheckInvalidState(WeibullDistribution(2, 3), &WeibullDistribution::shape, invalid);
		CheckInvalidState(WeibullDistribution(2, 3), &WeibullDistribution::scale, invalid);
		CheckInvalidState(ParetoDistribution(3, 2), &ParetoDistribution::alpha, invalid);
		CheckInvalidState(ParetoDistribution(3, 2), &ParetoDistribution::xm, invalid);
		CHECK_THROWS_AS(NormalDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(CauchyDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(LogisticDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(LogNormalDistribution(0, invalid), StatisticsError);
		CHECK_THROWS_AS(ExponentialDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(PoissonDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(GammaDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(GammaDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(BetaDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(BetaDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(WeibullDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(WeibullDistribution(1, invalid), StatisticsError);
		CHECK_THROWS_AS(ParetoDistribution(invalid, 1), StatisticsError);
		CHECK_THROWS_AS(ParetoDistribution(1, invalid), StatisticsError);
	}
}

TEST_CASE("Mutable success probabilities are checked by every operation", "[statistics][distribution-parameters]") {
	for (Real invalid : {Real(-1), Real(2), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CheckInvalidState(BernoulliDistribution(Real(0.5)), &BernoulliDistribution::p, invalid);
		CheckInvalidState(BinomialDistribution(5, Real(0.5)), &BinomialDistribution::p, invalid);
		CheckInvalidState(GeometricDistribution(Real(0.5)), &GeometricDistribution::p, invalid);
		CheckInvalidState(NegativeBinomialDistribution(3, Real(0.5)), &NegativeBinomialDistribution::p, invalid);
		CHECK_THROWS_AS(BernoulliDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(BinomialDistribution(5, invalid), StatisticsError);
		CHECK_THROWS_AS(GeometricDistribution(invalid), StatisticsError);
		CHECK_THROWS_AS(NegativeBinomialDistribution(3, invalid), StatisticsError);
	}
	CheckInvalidState(GeometricDistribution(Real(0.5)), &GeometricDistribution::p, Real(0));
	CheckInvalidState(NegativeBinomialDistribution(3, Real(0.5)), &NegativeBinomialDistribution::p, Real(0));
	CHECK_THROWS_AS(GeometricDistribution(0), StatisticsError);
	CHECK_THROWS_AS(NegativeBinomialDistribution(3, 0), StatisticsError);
}

TEST_CASE("Mutable integer distribution invariants are checked by every operation", "[statistics][distribution-parameters]") {
	CheckInvalidState(TDistribution(5), &TDistribution::df, 0);
	CheckInvalidState(TDistribution(5), &TDistribution::df, -1);
	CheckInvalidState(BinomialDistribution(5, Real(0.5)), &BinomialDistribution::n, -1);
	CheckInvalidState(NegativeBinomialDistribution(3, Real(0.5)), &NegativeBinomialDistribution::r, 0);
	CheckInvalidState(NegativeBinomialDistribution(3, Real(0.5)), &NegativeBinomialDistribution::r, -1);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::N, -1);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::N, 9);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::N, 4);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::K, -1);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::K, 51);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::n, -1);
	CheckInvalidState(HypergeometricDistribution(50, 10, 5), &HypergeometricDistribution::n, 51);
	CHECK_THROWS_AS(TDistribution(0), StatisticsError);
	CHECK_THROWS_AS(ChiSquareDistribution(0), StatisticsError);
	CHECK_THROWS_AS(FDistribution(0, 1), StatisticsError);
	CHECK_THROWS_AS(FDistribution(1, 0), StatisticsError);
	CHECK_THROWS_AS(BinomialDistribution(-1, Real(0.5)), StatisticsError);
	CHECK_THROWS_AS(NegativeBinomialDistribution(0, Real(0.5)), StatisticsError);
	CHECK_THROWS_AS(HypergeometricDistribution(-1, 0, 0), StatisticsError);
	CHECK_THROWS_AS(HypergeometricDistribution(50, -1, 5), StatisticsError);
	CHECK_THROWS_AS(HypergeometricDistribution(50, 51, 5), StatisticsError);
	CHECK_THROWS_AS(HypergeometricDistribution(50, 10, -1), StatisticsError);
	CHECK_THROWS_AS(HypergeometricDistribution(50, 10, 51), StatisticsError);
	CHECK_THROWS_AS(UniformDistribution(1, 1), StatisticsError);
	CHECK_THROWS_AS(UniformDistribution(2, 1), StatisticsError);
}

TEST_CASE("Finite valid parameters preserve endpoints clipping and undefined moments", "[statistics][distribution-parameters]") {
	const Real inf = std::numeric_limits<Real>::infinity();
	const Real maximum = std::numeric_limits<Real>::max();
	CHECK(GammaDistribution(Real(0.5)).pdf(0) == inf);
	CHECK(BetaDistribution(Real(0.5), Real(0.5)).pdf(0) == inf);
	CHECK(BetaDistribution(Real(0.5), Real(0.5)).pdf(1) == inf);
	CHECK(WeibullDistribution(Real(0.5)).pdf(0) == inf);
	CHECK_THROWS_AS(TDistribution(1).mean(), StatisticsError);
	CHECK_THROWS_AS(TDistribution(2).variance(), StatisticsError);
	CHECK_THROWS_AS(ParetoDistribution(1).mean(), StatisticsError);
	CHECK_THROWS_AS(ParetoDistribution(2).variance(), StatisticsError);
	CHECK_THROWS_AS(FDistribution(1, 2).mean(), StatisticsError);
	CHECK_THROWS_AS(FDistribution(1, 4).variance(), StatisticsError);
	CHECK(NormalDistribution(maximum, maximum).mean() == maximum);
	CHECK(NormalDistribution(0, maximum).variance() == inf);
	CHECK(LogNormalDistribution(maximum, 1).mean() == inf);
	CHECK(UniformDistribution(-maximum, maximum).variance() == inf);
	CHECK_NOTHROW(GammaDistribution(std::numeric_limits<Real>::min(), 1));
	CHECK(NormalDistribution().cdf(inf) == 1);
	CHECK(NormalDistribution().cdf(-inf) == 0);
	CHECK(std::isnan(NormalDistribution().cdf(std::numeric_limits<Real>::quiet_NaN())));
	CHECK_THROWS_AS(ExponentialDistribution(1).cdf(-1), StatisticsError);
	CHECK(ExponentialDistribution(1).inverseCdf(0) == 0);
	CHECK_THROWS_AS(ExponentialDistribution(1).inverseCdf(1), StatisticsError);
	CHECK(UniformDistribution(2, 3).inverseCdf(0) == 2);
	CHECK(UniformDistribution(2, 3).inverseCdf(1) == 3);
	CHECK(GammaDistribution(1).inverseCdf(-inf) == 0);
	CHECK(GammaDistribution(1).inverseCdf(inf) == inf);
	CHECK(BetaDistribution(1, 1).inverseCdf(-inf) == 0);
	CHECK(BetaDistribution(1, 1).inverseCdf(inf) == 1);
	CHECK(WeibullDistribution(1).inverseCdf(-inf) == 0);
	CHECK(WeibullDistribution(1).inverseCdf(inf) == inf);
	CHECK(ParetoDistribution(1, 2).inverseCdf(-inf) == 2);
	CHECK(ParetoDistribution(1, 2).inverseCdf(inf) == inf);
	CHECK(LogNormalDistribution().inverseCdf(-inf) == 0);
	CHECK(LogNormalDistribution().inverseCdf(inf) == inf);
	CHECK(BernoulliDistribution(0).pmf(0) == 1);
	CHECK(BernoulliDistribution(1).pmf(1) == 1);
	CHECK(BinomialDistribution(0, 0).pmf(0) == 1);
	CHECK(BinomialDistribution(5, 1).pmf(5) == 1);
	CHECK(GeometricDistribution(1).pmf(0) == 1);
	CHECK(NegativeBinomialDistribution(3, 1).pmf(0) == 1);
	CHECK(HypergeometricDistribution(0, 0, 0).pmf(0) == 1);
}

TEST_CASE("Valid distribution sampling retains explicit and default engine sequences", "[statistics][distribution-parameters]") {
	CheckSamplingEngineAgreement(NormalDistribution(2, 3));
	CheckSamplingEngineAgreement(TDistribution(5));
	CheckSamplingEngineAgreement(ChiSquareDistribution(5));
	CheckSamplingEngineAgreement(FDistribution(5, 9));
	CheckSamplingEngineAgreement(CauchyDistribution(2, 3));
	CheckSamplingEngineAgreement(ExponentialDistribution(2));
	CheckSamplingEngineAgreement(LogisticDistribution(2, 3));
	CheckSamplingEngineAgreement(UniformDistribution(2, 3));
	CheckSamplingEngineAgreement(GammaDistribution(Real(0.5), 3));
	CheckSamplingEngineAgreement(GammaDistribution(2, 3));
	CheckSamplingEngineAgreement(BetaDistribution(2, 3));
	CheckSamplingEngineAgreement(WeibullDistribution(2, 3));
	CheckSamplingEngineAgreement(ParetoDistribution(3, 2));
	CheckSamplingEngineAgreement(LogNormalDistribution(2, 3));
	CheckSamplingEngineAgreement(BernoulliDistribution(Real(0.5)));
	CheckSamplingEngineAgreement(BinomialDistribution(5, Real(0.5)));
	CheckSamplingEngineAgreement(PoissonDistribution(2));
	CheckSamplingEngineAgreement(GeometricDistribution(Real(0.5)));
	CheckSamplingEngineAgreement(NegativeBinomialDistribution(3, Real(0.5)));
	CheckSamplingEngineAgreement(HypergeometricDistribution(50, 10, 5));
}
