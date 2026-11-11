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
		static_assert(Precision::StrictTolerance <= Precision::ZeroTolerance);
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
		const Real near = Precision::EqualityTolerance / Real(4);
		const Real far = Real(4) * Precision::EqualityTolerance;
		REQUIRE(1 + near != 1);
		REQUIRE(IsNearlyEqual(Real(1), 1 + near));
		REQUIRE_FALSE(IsNearlyEqual(Real(1), 1 + far));
		const Real large = TOL(1e12, 1e6), increment = TOL(1, 0.25);
		REQUIRE(large + increment != large);
		REQUIRE(IsNearlyEqual(large, large + increment, 0, TOL(1e-10, Precision::EqualityTolerance)));
		REQUIRE_FALSE(IsNearlyEqual(large, large + increment, Real(1e-10), 0));
		REQUIRE(IsNearlyEqual(MML::Real(0.0), MML::Real(1e-12), MML::Real(1e-10), MML::Real(0.0)));            // absolute part carries near-zero values

		REQUIRE(IsWithinAbs(MML::Real(0.0), MML::Real(1e-13), MML::Real(1e-12)));
		REQUIRE_FALSE(IsWithinAbs(MML::Real(0.0), MML::Real(1e-11), MML::Real(1e-12)));

		REQUIRE(IsWithinRel(MML::Real(0.0), MML::Real(0.0), MML::Real(1e-10)));                      // the historical (0,0) == false defect
		REQUIRE(IsWithinRel(MML::Real(5.0), MML::Real(5.0), MML::Real(1e-10)));
		const Real tiny = TOL(1e-300, std::numeric_limits<Real>::min());
		REQUIRE(tiny > 0);
		REQUIRE(IsWithinRel(tiny, tiny, Real(1e-10)));
		REQUIRE(IsWithinRel(Real(100), Real(100) + TOL(1e-9, Real(400) * Precision::Eps),
			TOL(1e-10, Precision::EqualityTolerance)));
		REQUIRE_FALSE(IsWithinRel(Real(0), tiny, Real(1e-10)));

		REQUIRE(IsNearlyZero(Precision::ZeroTolerance / Real(4)));
		REQUIRE_FALSE(IsNearlyZero(Real(4) * Precision::ZeroTolerance));
		REQUIRE(IsNearlyZero(MML::Real(1e-7), Precision::LooseTolerance));

		// scale-aware negligibility: a pivot against the matrix norm
		const Real pivot = Precision::RoundoffZeroThreshold;
		REQUIRE(IsNegligibleRelativeTo(pivot, Real(1e3) * 10));
		REQUIRE_FALSE(IsNegligibleRelativeTo(pivot, Real(1)));
		REQUIRE(IsNegligibleRelativeTo(0.0, 0.0));                  // 0 <= 0: exact zero is negligible at any scale

		REQUIRE(IsUnsafeDenominator(0.0));
		REQUIRE(IsUnsafeDenominator(Precision::DivisionSafetyThreshold / 4));
		REQUIRE_FALSE(IsUnsafeDenominator(4 * Precision::DivisionSafetyThreshold));

		REQUIRE(AnglesAreEqual(Constants::PI, -Constants::PI, MML::Real(1e-12)));
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
		REQUIRE(isNearlyEqual(MML::Real(1.0), MML::Real(1.0) + MML::Real(1e-11), MML::Real(1e-10)) == IsNearlyEqual(MML::Real(1.0), MML::Real(1.0) + MML::Real(1e-11), MML::Real(1e-10)));
		REQUIRE(isNearlyEqual(MML::Real(1e12), MML::Real(1e12) + MML::Real(1.0), MML::Real(0.0), MML::Real(1e-10)) == IsNearlyEqual(MML::Real(1e12), MML::Real(1e12) + MML::Real(1.0), MML::Real(0.0), MML::Real(1e-10)));
		REQUIRE(isWithinAbsPrec(MML::Real(0.0), MML::Real(1e-13), MML::Real(1e-12)) == IsWithinAbs(MML::Real(0.0), MML::Real(1e-13), MML::Real(1e-12)));
		REQUIRE(isWithinRelPrec(MML::Real(0.0), MML::Real(0.0), MML::Real(1e-10)));                  // fixed through the forwarder too
		REQUIRE(isNearlyZero(MML::Real(1e-13)) == IsNearlyZero(MML::Real(1e-13)));
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
	}
}
