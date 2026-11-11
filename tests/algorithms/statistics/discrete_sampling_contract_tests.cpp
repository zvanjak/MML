#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Statistics/DiscreteDistributions.h>
#endif
#include <catch2/catch_all.hpp>
#include <cstdint>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace MML;
using namespace MML::Statistics;

namespace {
	struct FixedUrbg {
		using result_type = std::uint64_t;
		static constexpr result_type min() { return 0; }
		static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }
		result_type value = result_type(1) << 63;
		std::size_t calls = 0;
		std::size_t limit = 2048;
		std::vector<result_type> sequence;
		result_type operator()() {
			const auto index = calls++;
			if (index >= limit || (!sequence.empty() && index >= sequence.size()))
				throw std::runtime_error("Controlled URBG draw limit exceeded");
			return sequence.empty() ? value : sequence[index];
		}
	};

	std::uint64_t RawFraction(std::uint64_t numerator, unsigned denominatorBits) {
		return numerator << (64 - denominatorBits);
	}

	std::uint64_t RawBelowOne() {
		// Use a fraction exactly representable in float, not a raw draw that rounds to 1 in libc++.
		return std::is_same_v<Real, float> ? RawFraction(16777215, 24) : FixedUrbg::max() - 2047;
	}

	Real LargestSupportedPoissonRate() {
		Real rate = Real(std::numeric_limits<int>::max());
		if (StatisticsDetail::DiscreteSamplingReal(rate) > std::numeric_limits<int>::max())
			rate = std::nextafter(rate, Real(0));
		return rate;
	}
}

TEST_CASE("Large Poisson samples use bounded central transformed proposals", "[statistics][discrete-sampling-contract]") {
	for (Real rate : {Real(1000), Real(1000000)}) {
		FixedUrbg gen;
		CHECK(PoissonDistribution(rate).sample(gen) == static_cast<int>(rate));
		CHECK(gen.calls == 2);
	}
}

TEST_CASE("Unrepresentable Poisson sampling rates fail before RNG use", "[statistics][discrete-sampling-contract]") {
	for (Real rate : {Real(std::numeric_limits<int>::max()) + 1, std::numeric_limits<Real>::max()}) {
		FixedUrbg gen;
		gen.limit = 0;
		CHECK_THROWS_AS(PoissonDistribution(rate).sample(gen), NumericalMethodError);
		CHECK(gen.calls == 0);
		const auto before = Random::Engine();
		CHECK_THROWS_AS(PoissonDistribution(rate).sample(), NumericalMethodError);
		CHECK(Random::Engine() == before);
		CHECK(PoissonDistribution(rate).mean() == rate);
		CHECK(PoissonDistribution(rate).variance() == rate);
		CHECK(PoissonDistribution(rate).cdf(-1) == 0);
	}
}

TEST_CASE("Discrete sample conversion checks exact integer limits before casting", "[statistics][discrete-sampling-contract]") {
	using S = StatisticsDetail::DiscreteSamplingReal;
	const S maximum = std::numeric_limits<int>::max();
	CHECK(StatisticsDetail::CheckedDiscreteSample(0) == 0);
	CHECK(StatisticsDetail::CheckedDiscreteSample(maximum) == std::numeric_limits<int>::max());
	CHECK(StatisticsDetail::CheckedDiscreteSample(maximum + S(0.5)) == std::numeric_limits<int>::max());
	CHECK_THROWS_AS(StatisticsDetail::CheckedDiscreteSample(maximum + 1), NumericalMethodError);
	CHECK_THROWS_AS(StatisticsDetail::CheckedDiscreteSample(-S(0.5)), NumericalMethodError);
	CHECK_THROWS_AS(StatisticsDetail::CheckedDiscreteSample(std::numeric_limits<S>::infinity()), NumericalMethodError);
	CHECK_THROWS_AS(StatisticsDetail::CheckedDiscreteSample(std::numeric_limits<S>::quiet_NaN()), NumericalMethodError);
}

TEST_CASE("Geometric sampling retains tiny probabilities without unchecked casts", "[statistics][discrete-sampling-contract]") {
	FixedUrbg gen;
	gen.value = RawBelowOne();
	const Real p = TOL(1e-20, std::ldexp(Real(1), -37));
	REQUIRE(p > Real(0));
	REQUIRE(Real(1) - p == Real(1));
	CHECK(GeometricDistribution(p).sample(gen) == (std::is_same_v<Real, float> ? 8192 : 11102));
	CHECK(gen.calls == 1);
	if constexpr (std::is_same_v<Real, float>) {
		FixedUrbg quantized;
		quantized.value = gen.value;
		CHECK_THROWS_AS(GeometricDistribution(Real(1e-20)).sample(quantized), NumericalMethodError);
		CHECK(quantized.calls == 1);
	}
	FixedUrbg overflow;
	CHECK_THROWS_AS(GeometricDistribution(Real(1e-20)).sample(overflow), NumericalMethodError);
	CHECK(overflow.calls == 1);
	FixedUrbg extreme;
	CHECK_THROWS_AS(GeometricDistribution(std::numeric_limits<Real>::min()).sample(extreme), NumericalMethodError);
	FixedUrbg zero;
	zero.value = 0;
	CHECK(GeometricDistribution(Real(0.5)).sample(zero) == std::numeric_limits<Real>::min_exponent * -1 + 1);
	FixedUrbg certain;
	certain.limit = 0;
	CHECK(GeometricDistribution(1).sample(certain) == 0);
	CHECK(NegativeBinomialDistribution(std::numeric_limits<int>::max(), 1).sample(certain) == 0);
	CHECK(certain.calls == 0);
}

TEST_CASE("Geometric and negative-binomial samples check conversion and accumulation separately", "[statistics][discrete-sampling-contract]") {
	using S = StatisticsDetail::DiscreteSamplingReal;
	const int maximum = std::numeric_limits<int>::max();
	const Real p = Real(-std::expm1(std::log(S(0.5)) / (S(maximum) + S(0.25))));
	const int expectedLimit = std::is_same_v<Real, float> ? 2147483641 : maximum;
	FixedUrbg atLimit;
	CHECK(GeometricDistribution(p).sample(atLimit) == expectedLimit);
	CHECK(atLimit.calls == 1);
	FixedUrbg single;
	CHECK(NegativeBinomialDistribution(1, p).sample(single) == expectedLimit);
	CHECK(single.calls == 1);
	if constexpr (std::is_same_v<Real, float>) {
		FixedUrbg lowerProbability;
		CHECK_THROWS_AS(GeometricDistribution(std::nextafter(p, Real(0))).sample(lowerProbability), NumericalMethodError);
		CHECK(lowerProbability.calls == 1);
		FixedUrbg higherProbability;
		CHECK(GeometricDistribution(std::nextafter(p, Real(1))).sample(higherProbability) == 2147483457);
		CHECK(higherProbability.calls == 1);
	}
	FixedUrbg checkedSum;
	checkedSum.sequence = {std::is_same_v<Real, float> ? RawFraction(8388609, 24) : RawFraction(1, 1),
		RawBelowOne()};
	CHECK(NegativeBinomialDistribution(2, p).sample(checkedSum)
		== (std::is_same_v<Real, float> ? 2147483456 : maximum));
	CHECK(checkedSum.calls == 2);
	FixedUrbg sumOverflow;
	CHECK_THROWS_AS(NegativeBinomialDistribution(2, p).sample(sumOverflow), NumericalMethodError);
	CHECK(sumOverflow.calls == 2);
	FixedUrbg individualOverflow;
	CHECK_THROWS_AS(NegativeBinomialDistribution(2, Real(1e-20)).sample(individualOverflow), NumericalMethodError);
	CHECK(individualOverflow.calls == 1);
	FixedUrbg fits;
	CHECK(NegativeBinomialDistribution(3, Real(1e-9)).sample(fits)
		== (std::is_same_v<Real, float> ? 2079441597 : 2079441540));
	CHECK(fits.calls == 3);
	FixedUrbg tooMany;
	CHECK_THROWS_AS(NegativeBinomialDistribution(4, Real(1e-9)).sample(tooMany), NumericalMethodError);
	CHECK(tooMany.calls == 4);
}

TEST_CASE("Poisson rejection handles boundaries and checks range only after acceptance", "[statistics][discrete-sampling-contract]") {
	for (auto first : {RawFraction(0, 1), RawFraction(1, 10), RawFraction(1023, 10)}) {
		FixedUrbg gen;
		gen.sequence = {first, RawFraction(1, 1), RawFraction(1, 1), RawFraction(1, 1)};
		CHECK(PoissonDistribution(1000).sample(gen) == 1000);
		CHECK(gen.calls == 4);
	}
	FixedUrbg fullAcceptance;
	fullAcceptance.sequence = {RawFraction(1, 1), RawFraction(15, 4)};
	CHECK(PoissonDistribution(1000).sample(fullAcceptance) == 1000);
	CHECK(fullAcceptance.calls == 2);
	FixedUrbg logRejected;
	logRejected.sequence = {RawFraction(127, 7), RawFraction(1, 8), RawFraction(1, 1), RawFraction(1, 1)};
	CHECK(PoissonDistribution(1000).sample(logRejected) == 1000);
	CHECK(logRejected.calls == 4);
	FixedUrbg rejectedOutOfRange;
	rejectedOutOfRange.sequence = logRejected.sequence;
	const Real rate = LargestSupportedPoissonRate();
	REQUIRE(StatisticsDetail::DiscreteSamplingReal(rate) <= std::numeric_limits<int>::max());
	CHECK(PoissonDistribution(rate).sample(rejectedOutOfRange) == static_cast<int>(rate));
	CHECK(rejectedOutOfRange.calls == 4);
	FixedUrbg acceptedOutOfRange;
	acceptedOutOfRange.sequence = {RawFraction(3, 2), RawFraction(1, 1)};
	CHECK_THROWS_AS(PoissonDistribution(rate).sample(acceptedOutOfRange), NumericalMethodError);
	CHECK(acceptedOutOfRange.calls == 2);
}

TEST_CASE("Poisson rate support ends between adjacent representable Real values", "[statistics][discrete-sampling-contract]") {
	using S = StatisticsDetail::DiscreteSamplingReal;
	const Real supported = LargestSupportedPoissonRate();
	const Real unsupported = std::nextafter(supported, std::numeric_limits<Real>::infinity());
	REQUIRE(S(supported) <= std::numeric_limits<int>::max());
	REQUIRE(S(unsupported) > std::numeric_limits<int>::max());
	FixedUrbg center;
	CHECK(PoissonDistribution(supported).sample(center) == static_cast<int>(supported));
	CHECK(center.calls == 2);
	FixedUrbg invalid;
	invalid.limit = 0;
	CHECK_THROWS_AS(PoissonDistribution(unsupported).sample(invalid), NumericalMethodError);
	CHECK(invalid.calls == 0);
	const auto saved = Random::Engine();
	CHECK_THROWS_AS(PoissonDistribution(unsupported).sample(), NumericalMethodError);
	CHECK(Random::Engine() == saved);
}

TEST_CASE("Both Poisson algorithms enforce the approved exact work budget", "[statistics][discrete-sampling-contract]") {
	CHECK(Precision::PoissonSamplingMaxIterations == 100000);
	FixedUrbg small;
	small.value = FixedUrbg::max() - 2047;
	small.limit = Precision::PoissonSamplingMaxIterations + 1;
	CHECK_THROWS_AS(PoissonDistribution(1).sample(small), ConvergenceError);
	CHECK(small.calls == Precision::PoissonSamplingMaxIterations);
	FixedUrbg large;
	large.value = 0;
	large.limit = 2 * Precision::PoissonSamplingMaxIterations + 1;
	CHECK_THROWS_AS(PoissonDistribution(1000).sample(large), ConvergenceError);
	CHECK(large.calls == 2 * Precision::PoissonSamplingMaxIterations);
}

TEST_CASE("Poisson log acceptance remains accurate near large integer modes", "[statistics][discrete-sampling-contract]") {
	using S = StatisticsDetail::DiscreteSamplingReal;
	struct Reference { S count, rate, logMass; };
	for (const auto& ref : {
		Reference{0, 30, -30},
		Reference{31, 30, S(-2.655104721788493993606470618906736567L)},
		Reference{32, 30, S(-2.719643242926065165279394534590729495L)},
		Reference{1000, 1000, S(-4.372899506026296824156522672401259668L)},
		Reference{1000000, 1000000, S(-7.826693895520143127164859656014266612L)},
		Reference{2147483647, 2147483647, S(-11.66271983168979950133306285582726763L)},
		Reference{2147483000, 2147483647, S(-11.66281714606109317389718138222258841L)},
		Reference{1286, 1000, S(-41.97475008465682264147797950556017698L)},
		Reference{S(2147908000LL), 2147483647, S(-53.58714451510012342388963977863944549L)}
	}) {
		CHECK(std::abs(StatisticsDetail::PoissonSamplingLogPmf(ref.count, ref.rate) - ref.logMass)
			<= Precision::EqualityTolerance);
	}
}

TEST_CASE("Small-rate Poisson and ordinary geometric sampling preserve legacy engine sequences", "[statistics][discrete-sampling-contract]") {
	for (Real rate : {Real(0.125), Real(3), std::nextafter(Precision::PoissonSamplingSwitchRate, Real(0))}) {
		std::mt19937 actual(34019), expected(34019);
		std::uniform_real_distribution<Real> uniform(0, 1);
		for (int i = 0; i < 100; ++i) {
			Real product = 1;
			int k = 0;
			do { ++k; product *= uniform(expected); } while (product > std::exp(-rate));
			CHECK(PoissonDistribution(rate).sample(actual) == k - 1);
			CHECK(actual == expected);
		}
	}
	std::mt19937 actual(34020), expected(34020);
	std::uniform_real_distribution<Real> uniform(0, 1);
	for (int i = 0; i < 100; ++i) {
		const int geometric = static_cast<int>(std::floor(std::log(std::max(uniform(expected), std::numeric_limits<Real>::min()))
			/ std::log(Real(0.5))));
		CHECK(GeometricDistribution(Real(0.5)).sample(actual) == geometric);
		int negative = 0;
		for (int j = 0; j < 3; ++j)
			negative += static_cast<int>(std::floor(std::log(std::max(uniform(expected), std::numeric_limits<Real>::min()))
				/ std::log(Real(0.5))));
		CHECK(NegativeBinomialDistribution(3, Real(0.5)).sample(actual) == negative);
		CHECK(actual == expected);
	}
	const auto saved = Random::Engine();
	Random::SetSeed(34021);
	std::mt19937 explicitEngine(34021);
	CHECK(PoissonDistribution(1000000).sample() == PoissonDistribution(1000000).sample(explicitEngine));
	CHECK(Random::Engine() == explicitEngine);
	Random::Engine() = saved;
}

TEST_CASE("Poisson samples pass independent moment and CDF sanity checks across regimes", "[statistics][discrete-sampling-contract]") {
	using S = StatisticsDetail::DiscreteSamplingReal;
	constexpr int samples = 20000;
	for (Real rate : {Real(0.5), std::nextafter(Precision::PoissonSamplingSwitchRate, Real(0)),
		Precision::PoissonSamplingSwitchRate, Real(1000), Real(1000000), Real(100000000)}) {
		CAPTURE(rate);
		std::mt19937 gen(34022);
		S mean = 0, m2 = 0;
		int lower = 0;
		for (int i = 0; i < samples; ++i) {
			const int value = PoissonDistribution(rate).sample(gen);
			CHECK(value >= 0);
			if (value <= static_cast<int>(rate)) ++lower;
			const S delta = S(value) - mean;
			mean += delta / (i + 1);
			m2 += delta * (S(value) - mean);
		}
		// Eight standard errors make these seeded sanity checks conservative, not accuracy benchmarks.
		const S confidence = 8;
		CHECK(std::abs(mean - rate) <= confidence * std::sqrt(S(rate) / samples));
		CHECK(std::abs(m2 / (samples - 1) - rate)
			<= confidence * std::sqrt((S(rate) + 2 * S(rate) * rate) / (samples - 1)));
		if (rate == 30 || rate == 1000) {
			const S expected = rate == 30 ? S(0.548351512577911426150573781465540606L)
				: S(0.508409367168505991214259092872459958L);
			CHECK(std::abs(S(lower) / samples - expected)
				<= confidence * std::sqrt(expected * (1 - expected) / samples));
		}
	}
}
