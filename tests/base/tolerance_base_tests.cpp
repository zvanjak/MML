// TOLERANCE_POLICY in mml/base (task MinimalMathLibrary-wewi.4.2): the literals that used to sit in
// comparisons are now the named Precision:: constants, and behaviour at the thresholds is pinned.
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Intervals.h>
#include <mml/base/Polynom.h>
#include <mml/base/Matrix/MatrixBandDiag.h>
#include <mml/base/SparseMatrix/SparseMatrixCOO.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#include <mml/base/Quaternions.h>
#include <mml/base/Combinatorics.h>
#include <mml/base/AlgorithmTypes.h>
#endif

using namespace MML;

namespace MML::Tests::Base::ToleranceBaseTests
{
	TEST_CASE("ToleranceBase_normalization_threshold_is_RoundoffZeroThreshold", "[tolerance][base]")
	{
		const Real thr = Precision::RoundoffZeroThreshold;
		REQUIRE(thr == 100 * std::numeric_limits<Real>::epsilon());

		Vector3Cartesian tiny(thr * 0.5, 0.0, 0.0);
		Vector3Cartesian small(thr * 2.0, 0.0, 0.0);
		REQUIRE_THROWS_AS(tiny.Normalized(), VectorDimensionError);
		REQUIRE(small.Normalized()[0] == Catch::Approx(1.0));
		REQUIRE(tiny.NormalizedOrZero().IsZero());

		VectorN<Real, 2> v{ thr * 0.5, 0.0 };
		REQUIRE_THROWS_AS(v.Normalized(), VectorDimensionError);
	}

	TEST_CASE("ToleranceBase_named_defaults_in_signatures", "[tolerance][base]")
	{
		// MatrixBandDiag::IsSymmetric and SparseMatrixCSR::IsSymmetric default to the same family constant as Matrix
		Matrix<Real> data(2, 3, { 0.0, 1.0, 2.0, 2.0 + 0.5 * Precision::IsMatrixSymmetricTolerance, 3.0, 0.0 });
		BandDiagonalMatrix band(2, 1, 1, data);
		REQUIRE(band.IsSymmetric());
		REQUIRE_FALSE(band.IsSymmetric(Precision::IsMatrixSymmetricTolerance * 0.1));

		SparseMatrix::SparseMatrixCOO<Real> coo(2, 2);
		coo.addEntry(0, 1, 2.0);
		coo.addEntry(1, 0, 2.0 + 0.5 * Precision::IsMatrixSymmetricTolerance);
		SparseMatrix::SparseMatrixCSR<Real> csr(coo);
		REQUIRE(csr.IsSymmetric());

		ConfigBase cfg;
		REQUIRE(cfg.tolerance == Precision::EqualityTolerance);
	}

	TEST_CASE("ToleranceBase_interval_integer_detection_and_polynomial_reduce", "[tolerance][base]")
	{
		// Polynom::Reduce() drops leading coefficients that are zero up to roundoff (100 ulp), not more
		Polynom<Real> p({ 1.0, 2.0, 0.5 * Precision::RoundoffZeroThreshold });
		p.Reduce();
		REQUIRE(p.degree() == 1);
		Polynom<Real> q({ 1.0, 2.0, 2.0 * Precision::RoundoffZeroThreshold });
		q.Reduce();
		REQUIRE(q.degree() == 2);
	}

	TEST_CASE("ToleranceBase_gimbal_lock_and_slerp_thresholds_are_named", "[tolerance][base]")
	{
		// Euler extraction near the gimbal-lock pole must not produce NaN
		Quaternion q = Quaternion::FromEulerZYX(0.3, Constants::PI / 2, 0.0);
		Vec3Cart e = q.ToEulerZYX();
		REQUIRE(std::isfinite(e[0]));
		REQUIRE(std::isfinite(e[1]));
		REQUIRE(std::isfinite(e[2]));
		REQUIRE(std::abs(e[1]) == Catch::Approx(Constants::PI / 2).margin(1e-9));

		// slerp of two nearly identical rotations falls back to lerp and stays unit
		Quaternion a = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 0.5);
		Quaternion b = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 0.5 + 1e-9);
		REQUIRE(Quaternion::Slerp(a, b, 0.5).IsUnit());
	}

	TEST_CASE("ToleranceBase_bernoulli_zeta_tail_converges_to_roundoff", "[tolerance][base]")
	{
		// The zeta-series Bernoulli path used an unreachable 1e-18 relative cutoff; with Eps it still agrees
		REQUIRE(Combinatorics::BernoulliReal(40) == Catch::Approx(-261082718496449122051.0 / 13530.0).epsilon(1e-6));
		REQUIRE(Combinatorics::BernoulliReal(60) == Catch::Approx(-2.1399949257225335e+34).epsilon(1e-9));   // zeta path (B_60)
		REQUIRE(std::isfinite(Combinatorics::BernoulliReal(100)));
	}
}
