// TOLERANCE_POLICY regression tests (task MinimalMathLibrary-wewi.4.1): the root vocabulary in
// MMLPrecision.h / MMLBase.h behaves as the policy states, tiers drive the comparison families,
// and the deprecated lowerCamel helpers forward (including the (0,0) fix for the relative test).
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#endif

#include <cmath>
#include <limits>

using namespace MML;

namespace MML::Tests::Base::TolerancePolicyTests
{
	TEST_CASE("TolerancePolicy_tiers_are_ordered_and_drive_the_families", "[tolerance][policy]")
	{
		static_assert(Precision::StrictTolerance < Precision::ZeroTolerance);
		static_assert(Precision::ZeroTolerance <= Precision::EqualityTolerance);
		static_assert(Precision::EqualityTolerance <= Precision::ModerateTolerance);
		static_assert(Precision::ModerateTolerance <= Precision::LooseTolerance);
		static_assert(Precision::Eps == std::numeric_limits<Real>::epsilon());

		// comparison families are tiers, not independent literals
		static_assert(Precision::GeometryEpsilon == Precision::EqualityTolerance);
		static_assert(Precision::VectorIsEqualTolerance == Precision::EqualityTolerance);
		static_assert(Precision::MatrixIsEqualTolerance == Precision::EqualityTolerance);
		static_assert(Precision::IsMatrixSymmetricTolerance == Precision::EqualityTolerance);
		static_assert(Precision::NumericalZeroThreshold == Precision::ZeroTolerance);
		static_assert(Precision::DeterminantZeroThreshold == Precision::StrictTolerance);
		static_assert(Precision::DefaultTolerance == Precision::LooseTolerance);
		static_assert(Precision::DefaultToleranceStrict == Precision::ModerateTolerance);

		// finite-difference step scales are the stated roots of eps (to 3 significant digits)
		REQUIRE(Precision::SqrtEps == Catch::Approx(std::sqrt(Precision::Eps)).epsilon(1e-3));
		REQUIRE(Precision::CbrtEps == Catch::Approx(std::cbrt(Precision::Eps)).epsilon(1e-3));
		REQUIRE(Precision::FourthRootEps == Catch::Approx(std::pow(Precision::Eps, 0.25)).epsilon(1e-3));
		REQUIRE(Precision::FifthRootEps == Catch::Approx(std::pow(Precision::Eps, 0.2)).epsilon(1e-3));

		// the runtime-adjustable integration defaults come from Precision, not from literals
		REQUIRE(Real(Defaults::TrapezoidIntegrationEPS) == Precision::TrapezoidIntegrationTolerance);
		REQUIRE(Real(Defaults::RombergIntegrationEPS) == Precision::RombergIntegrationTolerance);
	}

	TEST_CASE("TolerancePolicy_comparison_vocabulary", "[tolerance][policy]")
	{
		REQUIRE(IsNearlyEqual(1.0, 1.0 + 1e-11));
		REQUIRE_FALSE(IsNearlyEqual(1.0, 1.0 + 1e-9));
		REQUIRE(IsNearlyEqual(1e12, 1e12 + 1.0, 0.0, 1e-10));      // relative part carries large values
		REQUIRE_FALSE(IsNearlyEqual(1e12, 1e12 + 1.0, 1e-10, 0.0)); // absolute part alone does not
		REQUIRE(IsNearlyEqual(0.0, 1e-12, 1e-10, 0.0));            // absolute part carries near-zero values

		REQUIRE(IsWithinAbs(0.0, 1e-13, 1e-12));
		REQUIRE_FALSE(IsWithinAbs(0.0, 1e-11, 1e-12));

		REQUIRE(IsWithinRel(0.0, 0.0, 1e-10));                      // the historical (0,0) == false defect
		REQUIRE(IsWithinRel(5.0, 5.0, 1e-10));
		REQUIRE(IsWithinRel(1e-300, 1e-300, 1e-10));
		REQUIRE(IsWithinRel(100.0, 100.0 + 1e-9, 1e-10));
		REQUIRE_FALSE(IsWithinRel(0.0, 1e-300, 1e-10));            // one side zero: nothing is relatively close to 0

		REQUIRE(IsNearlyZero(1e-13));
		REQUIRE_FALSE(IsNearlyZero(1e-11));
		REQUIRE(IsNearlyZero(1e-7, Precision::LooseTolerance));

		// scale-aware negligibility: a pivot against the matrix norm
		REQUIRE(IsNegligibleRelativeTo(1e-13, 1e3 * 10));           // 1e-13 <= eps * 1e4 ≈ 2.2e-12
		REQUIRE_FALSE(IsNegligibleRelativeTo(1e-13, 1.0));
		REQUIRE(IsNegligibleRelativeTo(0.0, 0.0));                  // 0 <= 0: exact zero is negligible at any scale

		REQUIRE(IsUnsafeDenominator(0.0));
		REQUIRE(IsUnsafeDenominator(1e-31));
		REQUIRE_FALSE(IsUnsafeDenominator(1e-20));

		REQUIRE(AnglesAreEqual(Constants::PI, -Constants::PI, 1e-12));
	}

	TEST_CASE("TolerancePolicy_geometry_epsilon_unified", "[tolerance][policy]")
	{
		static_assert(Precision::GeometryEpsilon == PrecisionValues<Real>::GeometryEpsilon);
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
		// deprecated lowerCamel helpers forward
		REQUIRE(isNearlyEqual(1.0, 1.0 + 1e-11, 1e-10) == IsNearlyEqual(1.0, 1.0 + 1e-11, 1e-10));
		REQUIRE(isNearlyEqual(1e12, 1e12 + 1.0, 0.0, 1e-10) == IsNearlyEqual(1e12, 1e12 + 1.0, 0.0, 1e-10));
		REQUIRE(isWithinAbsPrec(0.0, 1e-13, 1e-12) == IsWithinAbs(0.0, 1e-13, 1e-12));
		REQUIRE(isWithinRelPrec(0.0, 0.0, 1e-10));                  // fixed through the forwarder too
		REQUIRE(isNearlyZero(1e-13) == IsNearlyZero(1e-13));
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
	}
}
