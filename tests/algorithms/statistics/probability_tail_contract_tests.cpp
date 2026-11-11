#include "../../TestPrecision.h"
#include "../../ProbabilityTailTestUtils.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Statistics/Distributions.h>
#include <mml/algorithms/Statistics/DiscreteDistributions.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <limits>
#include <type_traits>

using namespace MML;
using namespace MML::Statistics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;
using MML::Testing::GammaTailTolerance;

namespace {
	void CheckProbability(Real actual, Real expected, Real tolerance = Precision::EqualityTolerance) {
		REQUIRE(std::isfinite(expected));
		REQUIRE(expected > 0);
		INFO("actual/expected=" << actual / expected << ", relative error=" << std::abs(actual / expected - 1)
			<< ", log(expected)=" << std::log(expected));
		CHECK(std::isfinite(actual));
		CHECK(actual > 0);
		CHECK_THAT(actual, WithinRel(expected, tolerance));
	}

	[[maybe_unused]] bool FloatQuantileBracketsProbability(Real z, Real p) {
		Real lower = z, upper = z;
		// One ULP for the inverse-erfc iterate plus rounding in sqrt(2) and its final product.
		for (int i = 0; i < 2; ++i) {
			lower = std::nextafter(lower, -std::numeric_limits<Real>::infinity());
			upper = std::nextafter(upper, std::numeric_limits<Real>::infinity());
		}
		return Real(0) < p && p < Real(1) &&
			0.5 * std::erfc(-static_cast<double>(lower) / std::sqrt(2.0)) <= p &&
			p <= 0.5 * std::erfc(-static_cast<double>(upper) / std::sqrt(2.0));
	}

	void CheckNormalRoundtrip(Real z, Real p) {
		if constexpr (std::is_same_v<Real, float>) {
			// Quantile spacing and z/sqrt(2) rounding amplify in the tail; check each stage independently.
			CHECK(FloatQuantileBracketsProbability(z, p));
			const Real tail = std::min(p, 1 - p);
			const Real inverseErfc = SpecialFunctions::ErfcInv(2 * tail);
			const double lower = std::nextafter(inverseErfc, -std::numeric_limits<Real>::infinity());
			const double upper = std::nextafter(inverseErfc, std::numeric_limits<Real>::infinity());
			CHECK(0.5 * std::erfc(upper) <= tail);
			CHECK(tail <= 0.5 * std::erfc(lower));
			const Real argument = -z / std::sqrt(Real(2));
			const Real expected = static_cast<Real>(0.5 * std::erfc(static_cast<double>(argument)));
			CheckProbability(StandardNormal::cdf(z), expected);
		} else CheckProbability(StandardNormal::cdf(z), p);
	}

	long double PoissonFiniteSum(Real rate, int k) {
		long double term = std::exp(-static_cast<long double>(rate)), sum = term;
		for (int i = 1; i <= k; ++i) {
			term *= static_cast<long double>(rate) / i;
			sum += term;
		}
		return sum;
	}
}

TEST_CASE("Standard normal CDF survival and two-sided values preserve small tails", "[statistics][probability-tails]") {
	for (Real z : {Real(0), Real(2), Real(8), Real(10), Real(20), Real(30), Real(37)}) {
		Real expected = Real(0.5) * std::erfc(z / std::sqrt(Real(2)));
		if (expected < std::numeric_limits<Real>::min())
			continue;
		CAPTURE(z);
		CheckProbability(StandardNormal::cdf(-z), expected);
		CheckProbability(StandardNormal::sf(z), expected);
		CheckProbability(StandardNormal::twoTailedPValue(z), 2 * expected);
		CheckProbability(StandardNormal::twoTailedPValue(-z), 2 * expected);
		CheckProbability(NormalDistribution().cdf(-z), expected);
	}
	CHECK(StandardNormal::cdf(-std::numeric_limits<Real>::infinity()) == 0);
	CHECK(StandardNormal::cdf(std::numeric_limits<Real>::infinity()) == 1);
	CHECK(StandardNormal::sf(std::numeric_limits<Real>::infinity()) == 0);
	CHECK(StandardNormal::sf(-std::numeric_limits<Real>::infinity()) == 1);
}

TEST_CASE("Normal CDF scaling preserves valid large-scale positive controls", "[statistics][probability-tails]") {
	const Real large = std::pow(Real(10), MML::Real(std::numeric_limits<Real>::max_exponent10));
	CheckProbability(NormalDistribution(0, large).cdf(large), Real(0.5) * std::erfc(-Real(1) / std::sqrt(Real(2))));
	const Real maximum = std::numeric_limits<Real>::max();
	CheckProbability(NormalDistribution(0, maximum).cdf(maximum), StandardNormal::cdf(1));
	CheckProbability(NormalDistribution(-maximum, maximum).cdf(maximum), StandardNormal::cdf(2));
	CheckProbability(NormalDistribution(maximum, maximum).cdf(-maximum), StandardNormal::cdf(-2));
	for (Real unit : {Real(0.125), Real(1), Real(1e3)})
		CheckProbability(NormalDistribution(3 * unit, unit).cdf(-7 * unit), StandardNormal::cdf(-10));
}

TEST_CASE("Inverse standard normal retains extreme probabilities and shared wrappers", "[statistics][probability-tails]") {
	// Independently generated with mpmath at 400 decimal digits.
	const Real probabilities[] = {Real(1e-20L), Real(1e-100L), Real(1e-300L)};
	const Real expected[] = {Real(-9.2623400897984075737L), Real(-21.273453560965324295L), Real(-37.047096299361199237L)};
	for (int i = 0; i < 3; ++i) {
		Real p = probabilities[i];
		if (p <= 0)
			continue;
		CAPTURE(p);
		const Real z = InverseStandardNormalCdf(p);
		CHECK(std::isfinite(z));
		CHECK_THAT(z, WithinRel(expected[i], Precision::EqualityTolerance));
		CHECK_THAT(StandardNormal::inverseCdf(p), WithinRel(expected[i], Precision::EqualityTolerance));
		CHECK_THAT(NormalDistribution(3, 2).inverseCdf(p), WithinRel(3 + 2 * expected[i], Precision::EqualityTolerance));
		CheckNormalRoundtrip(z, p);
		CheckProbability(LogNormalDistribution().inverseCdf(p), std::exp(expected[i]));
	}
}

TEST_CASE("Inverse normal preserves central upper and smallest-double probabilities", "[statistics][probability-tails]") {
	CHECK(InverseStandardNormalCdf(Real(0.5)) == 0);
	for (Real p : {Real(0.001), Real(0.1), Real(0.25), Real(0.75), Real(0.9), Real(0.999)})
		CheckProbability(StandardNormal::cdf(InverseStandardNormalCdf(p)), p);
	if constexpr (std::numeric_limits<Real>::digits == 53 && std::numeric_limits<Real>::max_exponent == 1024) {
		CHECK_THAT(InverseStandardNormalCdf(std::numeric_limits<Real>::denorm_min()),
			WithinRel(Real(-38.467405617144346251), Precision::EqualityTolerance));
		CHECK_THAT(InverseStandardNormalCdf(std::nextafter(Real(1), Real(0))),
			WithinRel(Real(8.2095361516013868556), Precision::EqualityTolerance));
	}
}

TEST_CASE("Inverse normal rejects probabilities outside its existing open interval", "[statistics][probability-tails]") {
	for (Real p : {Real(-1), Real(0), Real(1), Real(2), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(InverseStandardNormalCdf(p), StatisticsError);
		CHECK_THROWS_AS(StandardNormal::inverseCdf(p), StatisticsError);
		CHECK_THROWS_AS(NormalDistribution().inverseCdf(p), StatisticsError);
	}
}

TEST_CASE("Inverse normal stays finite and monotone across the representable binary range", "[statistics][probability-tails]") {
	Real previous = std::numeric_limits<Real>::infinity();
	for (int exponent = -1; exponent >= std::numeric_limits<Real>::min_exponent - std::numeric_limits<Real>::digits; exponent -= 17) {
		Real p = std::scalbn(Real(1), exponent);
		Real z = InverseStandardNormalCdf(p);
		CAPTURE(exponent, p);
		CHECK(std::isfinite(z));
		CHECK(z < previous);
		if (p >= std::numeric_limits<Real>::min())
			CheckNormalRoundtrip(z, p);
		previous = z;
	}
}

TEST_CASE("Lognormal lower CDF uses the same representable normal tail", "[statistics][probability-tails]") {
	for (Real z : {Real(-5), Real(-8), Real(-10)})
		CheckProbability(LogNormalDistribution().cdf(std::exp(z)), Real(0.5) * std::erfc(-z / std::sqrt(Real(2))));
	const Real maximum = std::numeric_limits<Real>::max();
	CheckProbability(LogNormalDistribution(0, maximum).cdf(std::exp(Real(1))), StandardNormal::cdf(Real(1) / maximum));
	CheckProbability(LogNormalDistribution(-maximum, maximum).cdf(1), StandardNormal::cdf(1));
}

TEST_CASE("Poisson lower CDF retains independently summed finite integer tails", "[statistics][probability-tails]") {
	for (Real rate : {Real(1), Real(5), Real(40), Real(100), Real(500), Real(700)}) {
		Real term = std::exp(-rate);
		if (term < std::numeric_limits<Real>::min())
			continue;
		Real sum = term;
		PoissonDistribution poisson(rate);
		for (int k = 0; k <= 5; ++k) {
			if (k > 0) {
				term *= rate / k;
				sum += term;
			}
			CAPTURE(rate, k);
			const Real expected = std::is_same_v<Real, float> ? static_cast<Real>(PoissonFiniteSum(rate, k)) : sum;
			CheckProbability(poisson.cdf(k), expected, GammaTailTolerance(Real(k) + 1, rate));
		}
		CHECK(poisson.cdf(-1) == 0);
		CHECK(poisson.cdf(std::numeric_limits<int>::max()) == 1);
	}
}

TEST_CASE("Float inverse normal quantiles bracket independent probabilities through the normal range",
	"[statistics][probability-tails]") {
	if constexpr (std::is_same_v<Real, float>) {
		for (int exponent = -1; exponent >= std::numeric_limits<Real>::min_exponent - 1; --exponent) {
			const Real p = std::scalbn(Real(1), exponent);
			const Real z = InverseStandardNormalCdf(p);
			CAPTURE(exponent, p, z);
			REQUIRE(std::isfinite(z));
			CheckNormalRoundtrip(z, p);
			if (exponent < -1) {
				Real displaced = z;
				for (int i = 0; i < 4; ++i)
					displaced = std::nextafter(displaced, -std::numeric_limits<Real>::infinity());
				CHECK_FALSE(FloatQuantileBracketsProbability(displaced, p));
				displaced = z;
				for (int i = 0; i < 4; ++i)
					displaced = std::nextafter(displaced, std::numeric_limits<Real>::infinity());
				CHECK_FALSE(FloatQuantileBracketsProbability(displaced, p));
			}
		}
	}
}

TEST_CASE("Float inverse normal endpoint controls remain finite without demanding subnormal relative accuracy",
	"[statistics][probability-tails]") {
	if constexpr (std::is_same_v<Real, float>) {
		for (Real p : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(),
			std::nextafter(Real(1), Real(0))}) {
			const Real z = InverseStandardNormalCdf(p);
			CAPTURE(p, z);
			CHECK(std::isfinite(z));
			CHECK(FloatQuantileBracketsProbability(z, p));
			CHECK(StandardNormal::inverseCdf(p) == z);
			CHECK(NormalDistribution().inverseCdf(p) == z);
			const Real actual = StandardNormal::cdf(z);
			CHECK(std::isfinite(actual));
			CHECK(actual > 0);
			if (p == std::numeric_limits<Real>::denorm_min()) CHECK(actual == p);
			else CheckNormalRoundtrip(z, p);
		}
	}
}

TEST_CASE("Float Poisson finite tails survive an underflowed initial exponential",
	"[statistics][probability-tails]") {
	if constexpr (std::is_same_v<Real, float>) {
		const Real rate = 110;
		REQUIRE(std::exp(-rate) == 0);
		for (int k : {5, 10, 20}) {
			const Real expected = static_cast<Real>(PoissonFiniteSum(rate, k));
			CAPTURE(k, expected);
			REQUIRE(expected > 0);
			const Real tolerance = GammaTailTolerance(Real(k) + 1, rate);
			const Real actual = PoissonDistribution(rate).cdf(k);
			if (std::isnormal(expected)) CheckProbability(actual, expected, tolerance);
			else {
				CHECK(std::isfinite(actual));
				CHECK(actual > 0);
				CHECK_THAT(actual, WithinAbs(expected, std::numeric_limits<Real>::denorm_min() + tolerance * expected));
			}
			CHECK_FALSE(WithinRel(expected, tolerance).match(Real(0)));
			CHECK_FALSE(WithinRel(expected, tolerance).match(expected * (1 + 4 * tolerance)));
		}
	}
}

TEST_CASE("Chi-square right tails consume the shared upper gamma directly", "[statistics][probability-tails]") {
	for (Real x : {Real(1), Real(10), Real(80), Real(200), Real(1000)}) {
		Real expected = std::exp(-x / 2);
		if (expected >= std::numeric_limits<Real>::min())
			CheckProbability(ChiSquareDistribution(2).rightTailPValue(x), expected);
	}
	CHECK(ChiSquareDistribution(2).rightTailPValue(0) == 1);
	CHECK(ChiSquareDistribution(2).rightTailPValue(-1) == 1);
	CHECK(ChiSquareDistribution(2).rightTailPValue(std::numeric_limits<Real>::infinity()) == 0);
}
