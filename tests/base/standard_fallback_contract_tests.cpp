#include <cmath>

namespace {
#if defined(MML_HAS_STD_SPECIAL_FUNCTIONS)
	constexpr int expectedSpecialFunctionSelection = MML_HAS_STD_SPECIAL_FUNCTIONS;
#elif defined(__cpp_lib_math_special_functions) && __cpp_lib_math_special_functions >= 201603L
	constexpr int expectedSpecialFunctionSelection = 1;
#else
	constexpr int expectedSpecialFunctionSelection = 0;
#endif
}

#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/StandardFunctions.h>
#endif
#include <catch2/catch_all.hpp>

using namespace MML;

namespace {
	void CheckValue(Real actual, Real expected) {
		CAPTURE(actual, expected, actual - expected, Precision::RoundoffZeroThreshold);
		CHECK(std::isfinite(actual));
		CHECK(IsNearlyEqual(actual, expected, Precision::RoundoffZeroThreshold));
	}
}

TEST_CASE("Standard function selection follows library capability or an explicit override", "[fallback-parity]") {
	STATIC_REQUIRE(MML_HAS_STD_SPECIAL_FUNCTIONS == 0 || MML_HAS_STD_SPECIAL_FUNCTIONS == 1);
	STATIC_REQUIRE(MML_HAS_STD_SPECIAL_FUNCTIONS == expectedSpecialFunctionSelection);
#ifdef MML_EXPECT_STANDARD_FALLBACK
	STATIC_REQUIRE(MML_HAS_STD_SPECIAL_FUNCTIONS == 0);
#endif
}

TEST_CASE("Associated Legendre omits the Condon Shortley phase", "[fallback-parity]") {
	for (Real x : {Real(-1), Real(-0.75), Real(-0.25), Real(0), Real(0.25), Real(0.75), Real(1)}) {
		CAPTURE(x);
		const Real root = std::sqrt((1 - x) * (1 + x));
		CheckValue(Functions::AssocLegendre(1, 1, x), root);
		CheckValue(Functions::AssocLegendre(2, 1, x), 3 * x * root);
		CheckValue(Functions::AssocLegendre(2, 2, x), 3 * (1 - x) * (1 + x));
		CheckValue(Functions::AssocLegendre(3, 1, x), Real(1.5) * (5 * x * x - 1) * root);
	}
}

TEST_CASE("Spherical Legendre retains its Condon Shortley phase and normalization", "[fallback-parity]") {
	for (Real theta : {Real(0), Real(0.25), Real(0.75), Real(1.25), Constants::PI / 2, Real(2.25), Real(2.75), Constants::PI}) {
		CAPTURE(theta);
		const Real s = std::sin(theta), c = std::cos(theta);
		CheckValue(Functions::SphLegendre(1, 0, theta), std::sqrt(3 / (4 * Constants::PI)) * c);
		CheckValue(Functions::SphLegendre(1, 1, theta), -std::sqrt(3 / (8 * Constants::PI)) * s);
		CheckValue(Functions::SphLegendre(2, 1, theta), -std::sqrt(15 / (8 * Constants::PI)) * s * c);
		CheckValue(Functions::SphLegendre(2, 2, theta), std::sqrt(15 / (32 * Constants::PI)) * s * s);
	}
}

TEST_CASE("Legendre wrappers match standard library conventions across degrees and orders", "[fallback-parity]") {
#if defined(__cpp_lib_math_special_functions) && __cpp_lib_math_special_functions >= 201603L
	for (unsigned int l = 0; l <= 12; ++l) {
		for (unsigned int m = 0; m <= l; ++m) {
			CAPTURE(l, m);
			for (Real x : {Real(-1), Real(-0.75), Real(-0.25), Real(0), Real(0.25), Real(0.75), Real(1)}) {
				CAPTURE(x);
				CheckValue(Functions::AssocLegendre(l, m, x), std::assoc_legendre(l, m, x));
			}
			for (Real theta : {Real(0), Real(0.25), Real(0.75), Real(1.25), Constants::PI / 2, Real(2.25), Real(2.75), Constants::PI}) {
				CAPTURE(theta);
				CheckValue(Functions::SphLegendre(l, m, theta), std::sph_legendre(l, m, theta));
			}
		}
	}
#else
	SUCCEED("Independent analytic and reference-value contracts cover platforms without standard special functions");
#endif
}

TEST_CASE("Complete elliptic E matches independent values and exact endpoints", "[fallback-parity]") {
	CheckValue(Functions::Comp_ellint_2(0), Constants::PI / 2);
	CheckValue(Functions::Comp_ellint_2(Real(0.5)), Real(1.467462209339427155459795266990916L));
	CheckValue(Functions::Comp_ellint_2(Real(0.75)), Real(1.318472107994620997371842794497931L));
	CheckValue(Functions::Comp_ellint_2(Real(0.9L)), Real(1.171697052781614141185913957957410L));
	CHECK(Functions::Comp_ellint_2(Real(1)) == 1);
	CHECK(Functions::Comp_ellint_2(Real(-1)) == 1);
	for (Real k : {Real(0.25), Real(0.5), Real(0.75), Real(0.9), std::nextafter(Real(1), Real(0))}) {
		CAPTURE(k);
		CHECK(Functions::Comp_ellint_2(-k) == Functions::Comp_ellint_2(k));
	}
}

TEST_CASE("Complete elliptic E matches independent near-endpoint references in either selection", "[fallback-parity]") {
	// libstdc++ loses accuracy near |k|=1; use independent 70-digit mpmath E(k^2) references.
	const Real references[][2] = {
		{0, Constants::PI / 2},
		{Real(0.25L), Real(1.545957256105465034950412439920611L)},
		{Real(0.5L), Real(1.467462209339427155459795266990916L)},
		{Real(0.75L), Real(1.318472107994620997371842794497931L)},
		{Real(0.875L), Real(1.2011106307369146977945187107974920523039457729777985L)},
		{Real(0.9L), Real(1.171697052781614141185913957957410L)},
		{Real(0.9375L), Real(1.1215931702323738421144018082563277876389259783330498L)},
		{Real(0.96875L), Real(1.0714005291435245968440387023217754976918342419680415L)},
		{Real(0.99L), Real(1.028475809028804000983887138518022L)},
		{Real(0.9999L), Real(1.000514500083781179564348852768442L)},
		{1 - std::ldexp(Real(1), -20), Real(1.0000071250860012986600886383993150546059388215752588L)},
		{1, 1}
	};
	for (const auto& reference : references) {
		CAPTURE(reference[0]);
		CheckValue(Functions::Comp_ellint_2(reference[0]), reference[1]);
		CheckValue(Functions::Comp_ellint_2(-reference[0]), reference[1]);
	}
	const Real neighbor = std::nextafter(Real(1), Real(0));
	const Real complementarySquared = (1 - neighbor) * (1 + neighbor);
	// DLMF 19.12.2: omitted O(k'^4 log(1/k')) terms are below roundoff at these neighbors.
	const Real expected = 1 + complementarySquared / 2 *
		(std::log(4 / std::sqrt(complementarySquared)) - Real(0.5));
	CheckValue(Functions::Comp_ellint_2(neighbor), expected);
	CheckValue(Functions::Comp_ellint_2(-neighbor), expected);
}

TEST_CASE("Modified special function domains propagate NaN and reject invalid real arguments", "[fallback-parity]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real inf = std::numeric_limits<Real>::infinity();
	for (Real x : {std::nextafter(Real(1), inf), std::nextafter(Real(-1), -inf), inf, -inf}) {
		CAPTURE(x);
		CHECK_THROWS_AS(Functions::AssocLegendre(2, 1, x), DomainError);
		CHECK_THROWS_AS(Functions::Comp_ellint_2(x), DomainError);
	}
	CHECK(std::isnan(Functions::AssocLegendre(0, 0, nan)));
	CHECK(std::isnan(Functions::AssocLegendre(2, 1, nan)));
	CHECK(std::isnan(Functions::AssocLegendre(1, 2, nan)));
	CHECK(std::isnan(Functions::Comp_ellint_2(nan)));
	CHECK(std::isnan(Functions::SphLegendre(2, 1, nan)));
	CHECK(Functions::AssocLegendre(2, 3, Real(0.5)) == 0);
	CHECK(Functions::SphLegendre(2, 3, Real(0.5)) == 0);
}
