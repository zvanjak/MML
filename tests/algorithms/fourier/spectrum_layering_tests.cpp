#include <catch2/catch_all.hpp>
#include <mml/algorithms/Fourier/FourierSpectrum.h>

using namespace MML;
using namespace MML::Fourier;

TEST_CASE("Spectrum frequency bins preserve even and odd layouts", "[Fourier][Spectrum][layering]")
{
	for (int n : {1, 5, 8}) {
		const auto one = FrequencyAxis(n, Real(40), FrequencyAxisLayout::OneSided);
		const auto two = FrequencyAxis(n, Real(40));
		const auto shifted = FrequencyAxis(n, Real(40), FrequencyAxisLayout::TwoSided, true);
		REQUIRE(one.size() == n / 2 + 1);
		REQUIRE(two.size() == n);
		REQUIRE(shifted.size() == n);
		const auto shiftedBins = FFTShift(two);
		for (int i = 0; i < one.size(); ++i)
			REQUIRE(IsNearlyEqual(one[i], Real(i) * Real(40) / n));
		for (int i = 0; i < n; ++i)
			REQUIRE(IsNearlyEqual(shifted[i], shiftedBins[i]));
		const auto restored = IFFTShift(shiftedBins);
		for (int i = 0; i < n; ++i)
			REQUIRE(restored[i] == two[i]);
	}
}

TEST_CASE("Squared coefficient magnitude has no implicit signal normalization", "[Fourier][Spectrum][layering]")
{
	const auto squared = MagnitudeSquared(Vector<Complex>{Complex(3, 4), Complex(-2, 0), Complex(0, 0)});
	REQUIRE(squared.size() == 3);
	REQUIRE(squared[0] == 25);
	REQUIRE(squared[1] == 4);
	REQUIRE(squared[2] == 0);
	REQUIRE_THROWS_AS(MagnitudeSquared(Vector<Complex>{}), FourierError);
}
