#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Fourier/Fourier.h>
#endif
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <type_traits>

using namespace MML;
using namespace MML::Fourier;

TEST_CASE("DCT reference resolves high modes and its inverse across input sizes", "[Fourier][DCT][float-contract]") {
	using Wide = std::common_type_t<Real, double>;
	const Wide pi = std::acos(Wide(-1));
	for (int n : {17, 64, 257, 1024})
		for (int mode : {0, 1, n - 1}) {
			CAPTURE(n, mode);
			Vector<Real> signal(n);
			for (int i = 0; i < n; ++i)
				signal[i] = Real(std::cos(pi * mode * (Wide(i) + Wide(0.5)) / n));
			const auto coefficients = DCT::ForwardII_Reference(signal);
			const auto restored = DCT::InverseII_Reference(coefficients);
			const Real budget = std::is_same_v<Real, float> ? Real(16) * Precision::Eps :
				Real(16 * n) * Precision::Eps;
			REQUIRE(coefficients.size() == n);
			REQUIRE(restored.size() == n);
			for (int i = 0; i < n; ++i) {
				const Real expected = i == mode ? (mode == 0 ? Real(2) : Real(1)) : Real(0);
				CHECK(IsWithinAbs(coefficients[i], expected, budget));
				CHECK(IsWithinAbs(restored[i], signal[i], budget));
			}
		}
}
