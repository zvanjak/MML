#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Fourier/FourierWindowing.h>
#include <mml/algorithms/Fourier/FourierSpectrum.h>
#endif
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>

using namespace MML;
using namespace MML::Fourier;

namespace {
	long double BesselI0IntegralReference(long double x) {
		// DLMF 10.32.1: I0 is the mean of exp(x*cos(theta)) over a full period.
		// A 64-point periodic trapezoidal rule resolves this fixture's |x| <= 8.6.
		constexpr int samples = 64;
		const long double pi = std::acos(-1.0L);
		long double sum = 0;
		for (int i = 0; i < samples; ++i)
			sum += std::exp(x * std::cos(2 * pi * i / samples));
		return sum / samples;
	}

	Real PhaseRoundoffBound(Real originalPhase, int operations = 8) {
		// Reducing rounded many-turn inputs cannot recover phase bits already
		// lost when the fixture was constructed; include that absolute error.
		return std::max(Precision::EqualityTolerance,
			Real(operations) * std::numeric_limits<Real>::epsilon() * std::max(Real(1), std::abs(originalPhase)));
	}
}

TEST_CASE("Independent Kaiser oracle matches known modified Bessel values", "[Fourier][FourierHelperContract]")
{
	CHECK(BesselI0IntegralReference(0) == 1);
	CHECK(std::abs(BesselI0IntegralReference(1) - 1.266065877752008335598244625214717538L) <=
		Precision::EqualityTolerance);
	CHECK(std::abs(BesselI0IntegralReference(8.6L) / 750.461159563166341717596245680364480L - 1) <=
		Precision::EqualityTolerance);
}

TEST_CASE("Parameterized singleton windows have unit weight after parameter validation", "[Fourier][FourierHelperContract]")
{
	for (const Real beta : {Real{0}, Real{1}, Real{8.6}, std::numeric_limits<Real>::max()}) {
		const auto window = Windows::Kaiser(1, beta);
		REQUIRE(window.size() == 1);
		CHECK(window[0] == 1);
	}
	for (const Real sigma : {std::numeric_limits<Real>::denorm_min(), Real{0.4}, Real{1}, std::numeric_limits<Real>::max()}) {
		const auto window = Windows::Gaussian(1, sigma);
		REQUIRE(window.size() == 1);
		CHECK(window[0] == 1);
	}
	for (const Real invalid : {Real{-1}, std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(Windows::Kaiser(1, invalid), FourierError);
		CHECK_THROWS_AS(Windows::Gaussian(1, invalid), FourierError);
	}
	CHECK_THROWS_AS(Windows::Gaussian(1, 0), FourierError);
	CHECK_THROWS_AS(Windows::Kaiser(0, 1), FourierError);
	CHECK_THROWS_AS(Windows::Gaussian(0, 1), FourierError);
}

TEST_CASE("Ordinary parameterized windows retain independent formula controls", "[Fourier][FourierHelperContract]")
{
	const Real beta = Real{8.6}, sigma = Real{0.4};
	for (int n : {2, 3, 8, 17}) {
		const auto kaiser = Windows::Kaiser(n, beta);
		const auto gaussian = Windows::Gaussian(n, sigma);
		const auto rectangular = Windows::Kaiser(n, 0);
		for (int i = 0; i < n; ++i) {
			const Real position = Real{2} * i / Real(n - 1) - Real{1};
			const Real expectedKaiser = static_cast<Real>(
				BesselI0IntegralReference(beta * std::sqrt(Real{1} - position * position)) /
				BesselI0IntegralReference(beta));
			const Real expectedGaussian = std::exp(-Real{0.5} * (position / sigma) * (position / sigma));
			CHECK(IsWithinAbs(kaiser[i], expectedKaiser, Precision::EqualityTolerance));
			CHECK(IsWithinAbs(gaussian[i], expectedGaussian, Precision::EqualityTolerance));
			CHECK(rectangular[i] == 1);
		}
	}
}

TEST_CASE("DCT verification rejects invalid input instead of claiming success", "[Fourier][FourierHelperContract]")
{
	for (int n : {1, 2, 3, 7, 32, 33}) {
		for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
			for (int slot : {0, n - 1}) {
				Vector<Real> data(n, Real{1});
				data[slot] = invalid;
				CHECK_THROWS_AS(DCT::VerifyRoundTrip(data), NumericInputError);
				CHECK_THROWS_AS(DCT::VerifyFastMatchesReference(data), NumericInputError);
			}
		}
	}
	const Vector<Real> data{1, 2, 3};
	for (const Real invalid : {Real{-1}, std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(DCT::VerifyRoundTrip(data, invalid), FourierError);
		CHECK_THROWS_AS(DCT::VerifyFastMatchesReference(data, invalid), FourierError);
	}
	CHECK_THROWS_AS(DCT::VerifyRoundTrip(Vector<Real>{}), FourierError);
	CHECK_THROWS_AS(DCT::VerifyFastMatchesReference(Vector<Real>{}), FourierError);
}

TEST_CASE("DCT verification actually compares all supported fast size families", "[Fourier][FourierHelperContract]")
{
	for (int n = 1; n <= 97; ++n) {
		Vector<Real> data(n);
		for (int i = 0; i < n; ++i)
			data[i] = std::sin(Real(i + 1)) + Real(i % 5) / Real{8};
		const auto fast = DCT::ForwardII_Fast(data), reference = DCT::ForwardII_Reference(data);
		const auto recovered = DCT::InverseII(DCT::ForwardII(data));
		Real inputL1 = 0;
		for (const Real value : data) inputL1 += std::abs(value);
		// Direct cosine angles grow with n; the normalized transform divides
		// the O(n*eps*L1) phase/sum bound by n, while inverse errors grow with n.
		const Real toleranceBound = std::max(Precision::EqualityTolerance,
			Real(8) * n * std::numeric_limits<Real>::epsilon() * inputL1);
		for (const Real tolerance : {Real{0}, Precision::EqualityTolerance, toleranceBound}) {
			bool matches = true, roundTrip = true;
			for (int i = 0; i < n; ++i) {
				matches = matches && std::isfinite(fast[i]) && std::isfinite(reference[i]) &&
					std::abs(fast[i] - reference[i]) <= tolerance;
				roundTrip = roundTrip && std::isfinite(recovered[i]) && std::abs(recovered[i] - data[i]) <= tolerance;
			}
			CHECK(DCT::VerifyFastMatchesReference(data, tolerance) == matches);
			CHECK(DCT::VerifyRoundTrip(data, tolerance) == roundTrip);
		}
		CHECK(DCT::VerifyFastMatchesReference(data, toleranceBound));
		CHECK(DCT::VerifyRoundTrip(data, toleranceBound));
	}
	const Vector<Real> large{std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()};
	CHECK_FALSE(DCT::VerifyRoundTrip(large, std::numeric_limits<Real>::max()));
	CHECK_FALSE(DCT::VerifyFastMatchesReference(large, std::numeric_limits<Real>::max()));
}

TEST_CASE("Phase unwrapping removes multiple turns and preserves signed half-turn ties", "[Fourier][Spectrum][FourierHelperContract]")
{
	for (const Real sign : {Real{-1}, Real{1}})
		for (int turns : {0, 1, 2, 17, 1000}) {
			const Real jump = sign * (Real{2} * turns + Real{0.25}) * Constants::PI;
			const auto unwrapped = UnwrapPhase(Vector<Real>{0, jump});
			CHECK(IsWithinAbs(unwrapped[1], sign * Constants::PI / Real{4}, PhaseRoundoffBound(jump)));
		}
	for (const Real sign : {Real{-1}, Real{1}}) {
		// 5*pi is not generally representable: rounding can put it on either
		// side of the branch cut. Probe both sides; the exact tie is tested below.
		const Real belowTie = std::nextafter(sign * Real{5} * Constants::PI, Real{0});
		const Real aboveTie = std::nextafter(sign * Real{5} * Constants::PI, sign * Constants::PosInf);
		const auto tie = UnwrapPhase(Vector<Real>{0, belowTie});
		CHECK(IsWithinAbs(tie[1], sign * Constants::PI, PhaseRoundoffBound(sign * Real{5} * Constants::PI)));
		CHECK(IsWithinAbs(UnwrapPhase(Vector<Real>{0, aboveTie})[1], -sign * Constants::PI,
			PhaseRoundoffBound(aboveTie)));
		const auto singleTie = UnwrapPhase(Vector<Real>{0, sign * Constants::PI});
		CHECK(singleTie[1] == sign * Constants::PI);
	}
	const Vector<Real> expected{Real{0.2}, Real{0.4}, Real{0.8}, Real{1.2}, Real{1.7}};
	Vector<Real> shifted = expected;
	for (int i = 1; i < shifted.size(); ++i)
		shifted[i] += Real((i % 2 == 0 ? 1 : -1) * 20 * i) * Constants::PI;
	const auto unwrapped = UnwrapPhase(shifted);
	for (int i = 0; i < expected.size(); ++i)
		CHECK(IsWithinAbs(unwrapped[i], expected[i], PhaseRoundoffBound(shifted[i], 8 * (i + 1))));
	CHECK(shifted[1] != unwrapped[1]);
}

TEST_CASE("Phase unwrapping retains custom threshold behavior and validates inputs", "[Fourier][Spectrum][FourierHelperContract]")
{
	const Real pi = Constants::PI;
	const auto narrow = UnwrapPhase(Vector<Real>{0, pi / Real{2}}, pi / Real{4});
	CHECK(IsWithinAbs(narrow[1], -Real{1.5} * pi, Precision::EqualityTolerance));
	const auto wide = UnwrapPhase(Vector<Real>{0, Real{5} * pi}, Real{4} * pi);
	CHECK(IsWithinAbs(wide[1], Real{3} * pi, Precision::EqualityTolerance));
	const auto unchanged = UnwrapPhase(Vector<Real>{0, Real{3} * pi}, Real{4} * pi);
	CHECK(unchanged[1] == Real{3} * pi);
	CHECK(UnwrapPhase(Vector<Real>{2})[0] == 2);
	for (const Real invalid : {Real{0}, Real{-1}, std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
		CHECK_THROWS_AS(UnwrapPhase(Vector<Real>{0}, invalid), FourierError);
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(UnwrapPhase(Vector<Real>{invalid}), NumericInputError);
		CHECK_THROWS_AS(UnwrapPhase(Vector<Real>{0, invalid}), NumericInputError);
	}
	CHECK_THROWS_AS(UnwrapPhase(Vector<Real>{}), FourierError);
	const Real huge = std::numeric_limits<Real>::max();
	const Vector<Real> extremeInput{0, huge, -huge, 0};
	const auto extreme = UnwrapPhase(extremeInput);
	for (int i = 0; i < extreme.size(); ++i)
		CHECK(std::isfinite(extreme[i]));
	for (int i = 1; i < extreme.size(); ++i)
		CHECK(std::abs(extreme[i] - extreme[i - 1]) <= pi + Precision::EqualityTolerance);
	for (int i = 0; i < extreme.size(); ++i)
		CHECK(std::abs(std::remainder(extreme[i] - std::fmod(extremeInput[i], Real{2} * pi), Real{2} * pi)) <=
			Precision::EqualityTolerance);
	CHECK_THROWS_AS(UnwrapPhase(Vector<Real>{0, huge, -huge / Real{2}, 0, huge, -huge / Real{2}, 0, huge},
		huge / Real{2}), NumericalMethodError);
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
		CHECK_THROWS_AS(FrequencyAxis(2, invalid), FourierError);
}

TEST_CASE("General phase unwrapping agrees with repeated whole-turn corrections", "[Fourier][Spectrum][FourierHelperContract]")
{
	const Real period = Real{2} * Constants::PI;
	for (const Real threshold : {Constants::PI / Real{4}, Constants::PI, Real{4} * Constants::PI})
		for (int turns = -30; turns <= 30; ++turns) {
			const Real jump = (Real{2} * turns + Real{0.3}) * Constants::PI;
			Real expected = jump;
			if (jump > threshold)
				while (expected > threshold)
					expected -= period;
			else if (jump < -threshold)
				while (expected < -threshold)
					expected += period;
			CHECK(IsWithinAbs(UnwrapPhase(Vector<Real>{0, jump}, threshold)[1], expected,
				PhaseRoundoffBound(jump, 8 + 2 * std::abs(turns))));
		}
	for (const Real direction : {Real{-1}, Real{1}}) {
		Vector<Complex> spectrum(64);
		for (int i = 0; i < spectrum.size(); ++i) {
			const Real phase = Real{0.2} + direction * Real(i) * Constants::PI / Real{3};
			spectrum[i] = Complex(std::cos(phase), std::sin(phase));
		}
		const auto unwrapped = UnwrapPhase(Phase(spectrum));
		for (int i = 0; i < spectrum.size(); ++i)
			CHECK(IsWithinAbs(unwrapped[i], Real{0.2} + direction * Real(i) * Constants::PI / Real{3},
				Precision::EqualityTolerance));
	}
}
