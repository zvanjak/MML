///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        bugfix_epic_i5u3_tests.cpp                                          ///
///  Description: Regression tests for P0 bugs fixed in epic i5u3 (AI improvements)   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include <atomic>
#include <cmath>
#include <sstream>
#include <thread>

#include <mml/MMLBase.h>

#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/algorithms/MatrixAlg.h>
#include <mml/base/Tensor/Tensor2.h>
#include <mml/base/Tensor/Tensor3.h>
#include <mml/base/Function.h>
#include <mml/base/Intervals.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#include <mml/base/ODESystem.h>
#include <mml/base/DAESystem.h>
#include <mml/base/StandardFunctions.h>
#include <mml/base/Algebra/ModInt.h>
#include <mml/base/ChebyshevApproximation.h>
#include <mml/core/Derivation.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/Curves.h>
#include <mml/core/OrthogonalBasis/ChebyshevBasis.h>
#include <mml/systems/LinearSystem.h>
#include <mml/interfaces/IODESystemStepCalculator.h>
#include <mml/interfaces/IODESystemStepper.h>
#include <mml/algorithms/EigenSystemSolvers.h>
#include <mml/algorithms/ODESolvers/ODESteppers.h>
#include <mml/algorithms/Fourier/Fourier.h>
#include <mml/algorithms/RootFinding/RootFindingBracketing.h>
#include <mml/algorithms/Optimization/LP/LPSolvers.h>
#include <mml/algorithms/DAESolvers.h>
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#include <mml/systems/DynamicalSystem.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/core/Integration/Integration2DAdaptive.h>
#include <mml/core/Derivation/ADJacobians.h>
#include <mml/tools/Timer.h>
#include <mml/tools/ThreadPool.h>
#include <mml/tools/CsvUtils.h>
#include <mml/tools/ConsolePrinter.h>
#include <mml/tools/serializer/SerializerBase.h>

using namespace MML;

namespace
{
	constexpr Real LinearSolveTolerance = MML::Testing::Tol(REAL(1e-10), REAL(1e-6));
	constexpr Real OneSidedDerivativeTolerance = MML::Testing::Tol(REAL(1e-3), REAL(3e-2));
	constexpr Real EllipticIntegralTolerance = MML::Testing::Tol(REAL(1e-9), REAL(1e-6));
}

namespace MML::Tests::BugfixI5u3 {

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.1 — Vector3Cylindrical::operator-() must negate Z                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.1 - Cylindrical negation matches Cartesian negation", "[bugfix][i5u3][vector]")
{
	Vector3Cylindrical v(2.0, Constants::PI / 4, 3.0);
	Vector3Cylindrical neg = -v;

	// Convert both to Cartesian; -v in cylindrical must equal the negated Cartesian vector
	Real x, y, z, nx, ny, nz;
	Utils::CylindricalToCartesian(v.R(), v.Phi(), v.Z(), x, y, z);
	Utils::CylindricalToCartesian(neg.R(), neg.Phi(), neg.Z(), nx, ny, nz);

	REQUIRE(nx == Catch::Approx(-x).margin(1e-12));
	REQUIRE(ny == Catch::Approx(-y).margin(1e-12));
	REQUIRE(nz == Catch::Approx(-z).margin(1e-12));
}

TEST_CASE("BugFix_i5u3.1.1 - Cylindrical double negation is identity", "[bugfix][i5u3][vector]")
{
	Vector3Cylindrical v(1.5, 0.7, -2.5);
	Vector3Cylindrical vv = -(-v);

	Real x, y, z, xx, yy, zz;
	Utils::CylindricalToCartesian(v.R(), v.Phi(), v.Z(), x, y, z);
	Utils::CylindricalToCartesian(vv.R(), vv.Phi(), vv.Z(), xx, yy, zz);

	REQUIRE(xx == Catch::Approx(x).margin(1e-12));
	REQUIRE(yy == Catch::Approx(y).margin(1e-12));
	REQUIRE(zz == Catch::Approx(z).margin(1e-12));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.2 — MatrixNM Identity/diagonal ctor must not write past min(N,M)      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.2 - Non-square MatrixNM Identity sets min(N,M) diagonal", "[bugfix][i5u3][matrix]")
{
	// N > M: previously wrote _vals[2][2] out of bounds
	auto id32 = MatrixNM<Real, 3, 2>::Identity();
	REQUIRE(id32(0, 0) == 1.0);
	REQUIRE(id32(1, 1) == 1.0);
	REQUIRE(id32(0, 1) == 0.0);
	REQUIRE(id32(2, 0) == 0.0);
	REQUIRE(id32(2, 1) == 0.0);

	auto id23 = MatrixNM<Real, 2, 3>::Identity();
	REQUIRE(id23(0, 0) == 1.0);
	REQUIRE(id23(1, 1) == 1.0);
	REQUIRE(id23(0, 2) == 0.0);
	REQUIRE(id23(1, 2) == 0.0);
}

TEST_CASE("BugFix_i5u3.1.2 - Non-square MatrixNM scalar diagonal ctor is safe", "[bugfix][i5u3][matrix]")
{
	MatrixNM<Real, 3, 2> d(5.0);
	REQUIRE(d(0, 0) == 5.0);
	REQUIRE(d(1, 1) == 5.0);
	REQUIRE(d(2, 0) == 0.0);
	REQUIRE(d(2, 1) == 0.0);

	// Square path unchanged
	MatrixNM<Real, 3, 3> sq(2.0);
	REQUIRE(sq(0, 0) == 2.0);
	REQUIRE(sq(1, 1) == 2.0);
	REQUIRE(sq(2, 2) == 2.0);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.3 — Defaults::* must resolve to the calling thread's context          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.3 - Defaults read/write round trip on main thread", "[bugfix][i5u3][defaults]")
{
	int original = Defaults::VectorPrintWidth;
	Defaults::VectorPrintWidth = original + 7;
	REQUIRE(static_cast<int>(Defaults::VectorPrintWidth) == original + 7);
	Defaults::VectorPrintWidth = original;
	REQUIRE(static_cast<int>(Defaults::VectorPrintWidth) == original);
}

TEST_CASE("BugFix_i5u3.1.3 - Defaults writes in worker thread do not leak to main thread", "[bugfix][i5u3][defaults]")
{
	int mainBefore = Defaults::VectorPrintWidth;
	int workerSeen = -1;

	std::thread worker([&]() {
		Defaults::VectorPrintWidth = 99;
		workerSeen = Defaults::VectorPrintWidth;  // must see its own thread's value
	});
	worker.join();

	REQUIRE(workerSeen == 99);
	// Main thread's per-thread context must be unaffected
	REQUIRE(static_cast<int>(Defaults::VectorPrintWidth) == mainBefore);
}

TEST_CASE("BugFix_i5u3.1.3 - Algorithm parameter Defaults are per-thread too", "[bugfix][i5u3][defaults]")
{
	Real mainBefore = Defaults::TrapezoidIntegrationEPS;
	Real workerSeen = 0.0;

	std::thread worker([&]() {
		Defaults::TrapezoidIntegrationEPS = REAL(1.0e-9);
		workerSeen = Defaults::TrapezoidIntegrationEPS;
	});
	worker.join();

	REQUIRE(workerSeen == Catch::Approx(1.0e-9));
	REQUIRE(static_cast<Real>(Defaults::TrapezoidIntegrationEPS) == Catch::Approx(mainBefore));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.4 — transfTensor2 must preserve tensor variance (was inverted)        ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.4 - transfTensor2 preserves (0,2) covariant variance", "[bugfix][i5u3][tensor]")
{
	CoordTransfSphericalToCartesian transf;
	Vector3Spherical pos(2.0, Constants::PI / 3, Constants::PI / 5);

	Tensor2<3> covariant(2, 0);  // (nCovar=2, nContravar=0) - fully covariant
	covariant(0, 0) = 1.0; covariant(1, 1) = 2.0; covariant(2, 2) = 3.0;

	Tensor2<3> result = transf.transfTensor2(covariant, pos);

	REQUIRE(result.NumCovar() == 2);
	REQUIRE(result.NumContravar() == 0);
	REQUIRE_FALSE(result.IsContravar(0));
	REQUIRE_FALSE(result.IsContravar(1));
}

TEST_CASE("BugFix_i5u3.1.4 - transfTensor2 preserves (2,0) contravariant variance", "[bugfix][i5u3][tensor]")
{
	CoordTransfSphericalToCartesian transf;
	Vector3Spherical pos(2.0, Constants::PI / 3, Constants::PI / 5);

	Tensor2<3> contravariant(0, 2);  // (nCovar=0, nContravar=2) - fully contravariant
	contravariant(0, 0) = 1.0; contravariant(1, 1) = 1.0; contravariant(2, 2) = 1.0;

	Tensor2<3> result = transf.transfTensor2(contravariant, pos);

	REQUIRE(result.NumCovar() == 0);
	REQUIRE(result.NumContravar() == 2);
	REQUIRE(result.IsContravar(0));
	REQUIRE(result.IsContravar(1));
}

TEST_CASE("BugFix_i5u3.1.4 - transfTensor2 preserves mixed per-slot variance order", "[bugfix][i5u3][tensor]")
{
	CoordTransfSphericalToCartesian transf;
	Vector3Spherical pos(1.5, Constants::PI / 4, Constants::PI / 6);

	Tensor2<3> mixed(1, 1);        // ctor default: slot 0 covariant, slot 1 contravariant
	mixed(0, 0) = 1.0; mixed(1, 2) = 0.5;
	REQUIRE_FALSE(mixed.IsContravar(0));
	REQUIRE(mixed.IsContravar(1));

	Tensor2<3> result = transf.transfTensor2(mixed, pos);

	REQUIRE(result.NumCovar() == 1);
	REQUIRE(result.NumContravar() == 1);
	REQUIRE_FALSE(result.IsContravar(0));  // per-slot order preserved
	REQUIRE(result.IsContravar(1));
}

TEST_CASE("BugFix_i5u3.1.4 - transfTensor3 preserves (3,0) covariant variance", "[bugfix][i5u3][tensor]")
{
	CoordTransfSphericalToCartesian transf;
	Vector3Spherical pos(2.0, Constants::PI / 3, Constants::PI / 5);

	Tensor3<3> covariant(3, 0);
	covariant(0, 0, 0) = 1.0;

	Tensor3<3> result = transf.transfTensor3(covariant, pos);

	REQUIRE(result.NumCovar() == 3);
	REQUIRE(result.NumContravar() == 0);
	REQUIRE_FALSE(result.IsContravar(0));
	REQUIRE_FALSE(result.IsContravar(1));
	REQUIRE_FALSE(result.IsContravar(2));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.5 — LinearSystem::GetLU must return a real decomposition              ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.5 - GetLU returns filled L, U and P*A = L*U holds", "[bugfix][i5u3][linearsystem]")
{
	Matrix<Real> A(3, 3, {2.0, 1.0, 1.0,
	                      4.0, -6.0, 0.0,
	                      -2.0, 7.0, 2.0});
	Systems::LinearSystem<Real> sys(A);

	const auto& lu = sys.LUDecompose();

	int n = 3;
	// L must be unit lower triangular, U upper triangular - and NOT all zeros (the old stub)
	Real lSum = 0.0, uSum = 0.0;
	for (int i = 0; i < n; ++i) {
		REQUIRE(lu.L(i, i) == Catch::Approx(1.0));
		for (int j = i + 1; j < n; ++j)
			REQUIRE(lu.L(i, j) == Catch::Approx(0.0));
		for (int j = 0; j < i; ++j)
			REQUIRE(lu.U(i, j) == Catch::Approx(0.0));
		for (int j = 0; j < n; ++j) { lSum += std::abs(lu.L(i, j)); uSum += std::abs(lu.U(i, j)); }
	}
	REQUIRE(lSum > 0.0);
	REQUIRE(uSum > 0.0);

	// Permutation must be a valid permutation of 0..n-1
	REQUIRE(static_cast<int>(lu.permutation.size()) == n);
	std::vector<bool> seen(n, false);
	for (int i = 0; i < n; ++i) {
		REQUIRE(lu.permutation[i] >= 0);
		REQUIRE(lu.permutation[i] < n);
		seen[lu.permutation[i]] = true;
	}
	for (int i = 0; i < n; ++i)
		REQUIRE(seen[i]);

	// Reconstruction: row k of P*A is row permutation[k] of A; P*A = L*U
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			Real pa = A(lu.permutation[i], j);
			Real luVal = 0.0;
			for (int k = 0; k < n; ++k)
				luVal += lu.L(i, k) * lu.U(k, j);
			REQUIRE(luVal == Catch::Approx(pa).margin(1e-10));
		}

	// Determinant must match direct computation: 2*(-12-0) - 1*(8-0) + 1*(28-12) = -24 - 8 + 16 = -16
	REQUIRE(lu.determinant == Catch::Approx(-16.0).margin(1e-10));
}

TEST_CASE("BugFix_i5u3.1.5 - SolveMultiple solves all columns correctly", "[bugfix][i5u3][linearsystem]")
{
	Matrix<Real> A(3, 3, {4.0, 1.0, 0.0,
	                      1.0, 3.0, 1.0,
	                      0.0, 1.0, 2.0});
	// Two right-hand sides as columns
	Matrix<Real> B(3, 2, {1.0, 2.0,
	                      0.0, 1.0,
	                      1.0, 0.0});
	Systems::LinearSystem<Real> sys(A, B);

	Matrix<Real> X = sys.SolveMultiple();
	REQUIRE(X.rows() == 3);
	REQUIRE(X.cols() == 2);

	// Verify A*X = B column by column
	for (int j = 0; j < 2; ++j)
		for (int i = 0; i < 3; ++i) {
			Real ax = 0.0;
			for (int k = 0; k < 3; ++k)
				ax += A(i, k) * X(k, j);
			REQUIRE(ax == Catch::Approx(B(i, j)).margin(LinearSolveTolerance));
		}
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.6 — NDer2Left/Right overloads must use the same offset (3h)           ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.6 - NDer2Left auto-step and explicit-step overloads agree", "[bugfix][i5u3][derivation]")
{
	RealFunctionFromStdFunc f([](Real x) { return std::sin(x); });
	Real x = 1.0;
	Real h = Derivation::ScaleStep(Derivation::NDer2_h, x);

	REQUIRE(Derivation::NDer2Left(f, x) == Catch::Approx(Derivation::NDer2Left(f, x, h)));
	REQUIRE(Derivation::NDer2Right(f, x) == Catch::Approx(Derivation::NDer2Right(f, x, h)));
}

TEST_CASE("BugFix_i5u3.1.6 - NDer2Left never evaluates f at or beyond x", "[bugfix][i5u3][derivation]")
{
	// One-sided derivative for a boundary: the entire stencil (incl. error estimate)
	// must stay strictly left of x, matching the NDer4/6/8 family pattern.
	Real x = 2.0;
	Real maxArgSeen = -1e30;
	RealFunctionFromStdFunc f([&maxArgSeen](Real arg) {
		if (arg > maxArgSeen) maxArgSeen = arg;
		return arg * arg;
	});

	Real error = 0.0;
	Real deriv = Derivation::NDer2Left(f, x, &error);

	REQUIRE(maxArgSeen < x);                          // strictly one-sided
	REQUIRE(deriv == Catch::Approx(2.0 * x).epsilon(OneSidedDerivativeTolerance));  // approximates f'(x)=2x from the left
}

TEST_CASE("BugFix_i5u3.1.6 - NDer2Right never evaluates f at or below x", "[bugfix][i5u3][derivation]")
{
	Real x = 2.0;
	Real minArgSeen = 1e30;
	RealFunctionFromStdFunc f([&minArgSeen](Real arg) {
		if (arg < minArgSeen) minArgSeen = arg;
		return arg * arg;
	});

	Real error = 0.0;
	Real deriv = Derivation::NDer2Right(f, x, &error);

	REQUIRE(minArgSeen > x);
	REQUIRE(deriv == Catch::Approx(2.0 * x).epsilon(OneSidedDerivativeTolerance));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.7 — Interval infinite bounds must use Real limits (float-build UB)    ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.7 - Infinite-interval bounds are finite Real values", "[bugfix][i5u3][intervals]")
{
	// Previously stored double::max() into Real members: out-of-range conversion UB when Real=float.
	CompleteRInterval all;
	REQUIRE(std::isfinite(all.getLowerBound()));
	REQUIRE(std::isfinite(all.getUpperBound()));
	REQUIRE(all.getLowerBound() == -std::numeric_limits<Real>::max());
	REQUIRE(all.getUpperBound() == std::numeric_limits<Real>::max());

	NegInfToOpenInterval leftInf(5.0);
	REQUIRE(std::isfinite(leftInf.getLowerBound()));
	REQUIRE(leftInf.getLowerBound() == -std::numeric_limits<Real>::max());

	ClosedToInfInterval rightInf(-3.0);
	REQUIRE(std::isfinite(rightInf.getUpperBound()));
	REQUIRE(rightInf.getUpperBound() == std::numeric_limits<Real>::max());
}

TEST_CASE("BugFix_i5u3.1.7 - Infinite-interval containment still correct", "[bugfix][i5u3][intervals]")
{
	CompleteRInterval all;
	REQUIRE(all.contains(0.0));
	REQUIRE(all.contains(-1e30));
	REQUIRE(all.contains(1e30));

	NegInfToClosedInterval upTo(2.0);
	REQUIRE(upTo.contains(2.0));
	REQUIRE(upTo.contains(-1e30));
	REQUIRE_FALSE(upTo.contains(2.5));

	OpenToInfInterval from(1.0);
	REQUIRE(from.contains(1.5));
	REQUIRE_FALSE(from.contains(1.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.8 — SparseMatrixCSR norms must be complex-correct and real-valued      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.8 - Complex CSR Frobenius norm uses |v|^2", "[bugfix][i5u3][sparse]")
{
	// Entry 3+4i has |v| = 5; previously sum += v*v produced the complex (3+4i)^2 = -7+24i
	SparseMatrix::SparseMatrixCOO<Complex> coo(2, 2);
	coo.addEntry(0, 0, Complex(3.0, 4.0));
	coo.addEntry(1, 1, Complex(0.0, 2.0));   // |v| = 2
	SparseMatrix::SparseMatrixCSR<Complex> csr(coo);

	Real frob = csr.normFrobenius();   // sqrt(25 + 4) = sqrt(29)
	REQUIRE(frob == Catch::Approx(std::sqrt(29.0)));

	Real inf = csr.normInf();          // max(5, 2) = 5
	REQUIRE(inf == Catch::Approx(5.0));

	Real one = csr.norm1();            // max col sums: max(5, 2) = 5
	REQUIRE(one == Catch::Approx(5.0));
}

TEST_CASE("BugFix_i5u3.1.8 - Real CSR norms unchanged", "[bugfix][i5u3][sparse]")
{
	SparseMatrix::SparseMatrixCOO<Real> coo(2, 2);
	coo.addEntry(0, 0, 3.0);
	coo.addEntry(0, 1, -4.0);
	coo.addEntry(1, 1, 12.0);
	SparseMatrix::SparseMatrixCSR<Real> csr(coo);

	REQUIRE(csr.normFrobenius() == Catch::Approx(13.0));  // sqrt(9+16+144)
	REQUIRE(csr.normInf() == Catch::Approx(12.0));         // max(3+4, 12)
	REQUIRE(csr.norm1() == Catch::Approx(16.0));           // max(3, 4+12)
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.9 - diagonal dominance must reject non-square matrices                ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.9 - isDiagonallyDominant returns false for non-square", "[bugfix][i5u3][matrix]")
{
	// rows > cols: (2,2) read was out of bounds
	Matrix<Real> tall(3, 2, {5.0, 1.0,
	                         1.0, 5.0,
	                         1.0, 1.0});
	REQUIRE_FALSE(MatrixAlg::IsDiagonallyDominant(tall));

	Matrix<Real> wide(2, 3, {5.0, 1.0, 1.0,
	                         1.0, 5.0, 1.0});
	REQUIRE_FALSE(MatrixAlg::IsDiagonallyDominant(wide));
}

TEST_CASE("BugFix_i5u3.1.9 - isDiagonallyDominant square semantics unchanged", "[bugfix][i5u3][matrix]")
{
	Matrix<Real> dominant(3, 3, {4.0, 1.0, 1.0,
	                             1.0, 5.0, 2.0,
	                             0.0, 1.0, 3.0});
	REQUIRE(MatrixAlg::IsDiagonallyDominant(dominant));

	Matrix<Real> notDominant(2, 2, {1.0, 5.0,
	                                5.0, 1.0});
	REQUIRE_FALSE(MatrixAlg::IsDiagonallyDominant(notDominant));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.10 — Line3D IsPerpendicular abs() and parallel-line distance          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.10 - Anti-parallel lines are NOT perpendicular", "[bugfix][i5u3][geometry]")
{
	// dot = -1: the old `dot < eps` check classified these as perpendicular
	Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
	Line3D xAxisReversed(Pnt3Cart(0, 5, 0), Vec3Cart(-1, 0, 0));
	REQUIRE_FALSE(xAxis.IsPerpendicular(xAxisReversed));

	// 135 degrees: dot = -sqrt(2)/2 < 0, also wrongly classified before
	Line3D diag(Pnt3Cart(0, 0, 0), Vec3Cart(-1, 1, 0));
	REQUIRE_FALSE(xAxis.IsPerpendicular(diag));

	// Genuinely perpendicular still detected
	Line3D yAxis(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0));
	REQUIRE(xAxis.IsPerpendicular(yAxis));
}

TEST_CASE("BugFix_i5u3.1.10 - Parallel line distance is perpendicular separation", "[bugfix][i5u3][geometry]")
{
	// Two parallel lines along x, separated by 5 in y.
	// Old formula projected p1->p2 onto the line direction: gave 0 here instead of 5.
	Line3D l1(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
	Line3D l2(Pnt3Cart(3, 5, 0), Vec3Cart(1, 0, 0));

	REQUIRE(l1.Dist(l2) == Catch::Approx(5.0));

	// Offset with both along-line and perpendicular components: distance is still 5
	Line3D l3(Pnt3Cart(100.0, 3.0, 4.0), Vec3Cart(1, 0, 0));  // perp offset (0,3,4), |.| = 5
	REQUIRE(l1.Dist(l3) == Catch::Approx(5.0));

	// Anti-parallel direction: same geometric lines family, same distance
	Line3D l4(Pnt3Cart(-7.0, 0.0, 5.0), Vec3Cart(-1, 0, 0));
	REQUIRE(l1.Dist(l4) == Catch::Approx(5.0));

	// Coincident lines: distance 0
	Line3D l5(Pnt3Cart(42.0, 0.0, 0.0), Vec3Cart(1, 0, 0));
	REQUIRE(l1.Dist(l5) == Catch::Approx(0.0).margin(1e-12));

	// The nearest-points overload's parallel branch must agree
	Real dist = 0.0;
	Pnt3Cart pa, pb;
	bool ok = l1.Dist(l3, dist, pa, pb);
	REQUIRE_FALSE(ok);  // parallel: nearest points undefined
	REQUIRE(dist == Catch::Approx(5.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.11 — Timer init-before-Start and ThreadPool enqueue-after-stop        ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.11 - Timer MarkTime before Start measures from construction", "[bugfix][i5u3][tools]")
{
	// Previously _startTime was default-constructed (epoch): intervals were garbage (~years)
	Timer t;
	t.MarkTime("no explicit Start");
	double interval = t.GetIntervalTime(0);
	REQUIRE(interval >= 0.0);
	REQUIRE(interval < 60.0);  // sane: seconds since construction, not since clock epoch
}

TEST_CASE("BugFix_i5u3.1.11 - Timer normal Start/MarkTime flow unchanged", "[bugfix][i5u3][tools]")
{
	Timer t;
	t.Start();
	t.MarkTime("phase 1");
	REQUIRE(t.GetMarkCount() == 1);
	REQUIRE(t.GetIntervalTime(0) >= 0.0);
}

TEST_CASE("BugFix_i5u3.1.11 - ThreadPool executes tasks normally", "[bugfix][i5u3][tools]")
{
	std::atomic<int> counter{0};
	{
		ThreadPool pool(2);
		for (int i = 0; i < 10; ++i)
			pool.enqueue([&counter] { counter.fetch_add(1); });
		pool.wait_for_tasks();
		REQUIRE(counter.load() == 10);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.12 — LinearSystem property caches ignored the tolerance argument      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.12 - isSymmetric cache is keyed on tolerance", "[bugfix][i5u3][systems]")
{
	// Nearly-symmetric matrix: off-diagonal mismatch of 1e-6
	Matrix<Real> A(2, 2, { 2.0, 1.0 + 1e-6,
	                       1.0, 2.0 });
	Systems::LinearSystem<Real> sys(A, Vector<Real>{1.0, 1.0});

	// Strict tolerance first: not symmetric
	REQUIRE(sys.IsSymmetric({REAL(1e-9), REAL(0.0)}) == false);
	// Loose tolerance afterwards must recompute, not return the cached strict result
	REQUIRE(sys.IsSymmetric({REAL(1e-3), REAL(0.0)}) == true);
	// And back again
	REQUIRE(sys.IsSymmetric({REAL(1e-9), REAL(0.0)}) == false);
}

TEST_CASE("BugFix_i5u3.1.12 - IsPositiveDefinite respects the requested tolerance", "[bugfix][i5u3][systems]")
{
	// SPD matrix with a tiny symmetry perturbation
	Matrix<Real> A(2, 2, { 4.0, 1.0 + 1e-6,
	                       1.0, 3.0 });
	Systems::LinearSystem<Real> sys(A, Vector<Real>{1.0, 1.0});

	REQUIRE_THROWS_AS(sys.IsPositiveDefinite(REAL(1e-9)), MatrixDimensionError);
	REQUIRE(sys.IsPositiveDefinite(REAL(1e-3)) == true);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.14 — Stepper interfaces lacked virtual destructors                    ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.14 - stepper interfaces have virtual destructors", "[bugfix][i5u3][interfaces]")
{
	STATIC_REQUIRE(std::has_virtual_destructor<IODESystemStepCalculator>::value);
	STATIC_REQUIRE(std::has_virtual_destructor<IODESystemStepper>::value);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.15 — Low-severity base/core cleanup batch                             ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.15 - Tensor2 setContravar keeps counts in sync", "[bugfix][i5u3][tensor]")
{
	Tensor2<3> t(2, 0);   // both indices covariant
	REQUIRE(t.NumCovar() == 2);
	REQUIRE(t.NumContravar() == 0);

	t.setContravar(0, true);
	REQUIRE(t.IsContravar(0) == true);
	REQUIRE(t.NumContravar() == 1);   // previously stayed 0 - desynced
	REQUIRE(t.NumCovar() == 1);

	t.setContravar(0, false);
	REQUIRE(t.NumContravar() == 0);
	REQUIRE(t.NumCovar() == 2);
}

TEST_CASE("BugFix_i5u3.1.15 - Tensor3 setContravar keeps counts in sync", "[bugfix][i5u3][tensor]")
{
	Tensor3<3> t(3, 0);
	t.setContravar(1, true);
	REQUIRE(t.NumContravar() == 1);
	REQUIRE(t.NumCovar() == 2);
}

TEST_CASE("BugFix_i5u3.1.15 - MatrixNM unary minus works on const object", "[bugfix][i5u3][matrix]")
{
	const MatrixNM<Real, 2, 2> m({ 1.0, -2.0, 3.0, -4.0 });
	MatrixNM<Real, 2, 2> neg = -m;   // previously did not compile: operator-() was non-const
	REQUIRE(neg(0, 0) == REAL(-1.0));
	REQUIRE(neg(0, 1) == REAL(2.0));
	REQUIRE(neg(1, 0) == REAL(-3.0));
	REQUIRE(neg(1, 1) == REAL(4.0));
}

TEST_CASE("BugFix_i5u3.1.15 - Vector::erase validates position and range", "[bugfix][i5u3][vector]")
{
	Vector<Real> v({ 1.0, 2.0, 3.0 });

	REQUIRE_THROWS_AS(v.erase(-1), VectorAccessBoundsError);
	REQUIRE_THROWS_AS(v.erase(3), VectorAccessBoundsError);
	REQUIRE_THROWS_AS(v.erase(0, 5), VectorAccessBoundsError);
	REQUIRE_THROWS_AS(v.erase(2, 1), VectorAccessBoundsError);

	v.erase(1);   // valid erase still works
	REQUIRE(v.size() == 2);
	REQUIRE(v[0] == REAL(1.0));
	REQUIRE(v[1] == REAL(3.0));

	v.erase(0, 2);
	REQUIRE(v.size() == 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug i5u3.1.16 — Low-severity algorithms/tools cleanup batch                      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_i5u3.1.16 - DAESolution fillValues grows safely beyond preallocation", "[bugfix][i5u3][algorithms]")
{
	// num_steps preallocation can underestimate the accepted-step count (bounded by
	// max_steps) - fillValues past the preallocated size must extend, not corrupt.
	DAESolution sol(0.0, 1.0, 2, 1, 4);   // preallocates 4+1 slots
	Vector<Real> x({ 1.0, 2.0 });
	Vector<Real> y({ 3.0 });

	for (int i = 0; i < 20; i++)
		REQUIRE_NOTHROW(sol.fillValues(i, 0.05 * i, x, y));

	REQUIRE(sol.capacity() >= 20);
	Vector<Real> t = sol.getTValues();
	REQUIRE(t[19] == Catch::Approx(0.95));
}

TEST_CASE("BugFix_i5u3.1.16 - Jacobi eigensolver validates symmetry", "[bugfix][i5u3][algorithms]")
{
	// Non-symmetric matrix: previously silently averaged with its transpose
	Matrix<Real> nonSym(2, 2, { 1.0, 5.0,
	                            0.0, 2.0 });
	REQUIRE_THROWS_AS(SymmMatEigenSolverJacobi::Solve(nonSym), MatrixDimensionError);

	// Explicit symmetric-part solve is available and works
	auto result = SymmMatEigenSolverJacobi::SolveSymmetricPart(nonSym);
	REQUIRE(result.converged);

	// Symmetric matrices still solve directly
	Matrix<Real> sym(2, 2, { 2.0, 1.0,
	                         1.0, 2.0 });
	auto res2 = SymmMatEigenSolverJacobi::Solve(sym);
	REQUIRE(res2.converged);
	REQUIRE(res2.eigenvalues[0] == Catch::Approx(1.0));
	REQUIRE(res2.eigenvalues[1] == Catch::Approx(3.0));
}

TEST_CASE("BugFix_i5u3.1.16 - stepper interpolate clamps outside last step", "[bugfix][i5u3][algorithms]")
{
	// dx/dt = x, x(0) = 1
	ODESystem sys(1, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) { dxdt[0] = x[0]; });

	DormandPrince5_Stepper stepper(sys);
	Vector<Real> x({ 1.0 });
	Vector<Real> dxdt({ 1.0 });
	auto res = stepper.doStep(0.0, x, dxdt, 0.1, 1e-8);
	REQUIRE(res.accepted);

	// Interpolating far outside the accepted step must clamp to the endpoints,
	// not silently Hermite-extrapolate
	Vector<Real> atStart = stepper.interpolate(-100.0);
	Vector<Real> atEnd = stepper.interpolate(100.0);
	Vector<Real> t0 = stepper.interpolate(0.0);
	Vector<Real> t1 = stepper.interpolate(res.hDone);

	REQUIRE(atStart[0] == Catch::Approx(t0[0]));
	REQUIRE(atEnd[0] == Catch::Approx(t1[0]));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.3.2 — Shared ODE step controller and Hermite interpolation machinery  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Task_i5u3.3.2 - PI step controller is configurable and resettable", "[ODE][i5u3][step-controller]")
{
	StepSizeControllerConfig config;
	config.safety = REAL(0.8);
	config.alpha = REAL(0.25);
	config.beta = REAL(0.1);
	config.minFactor = REAL(0.1);
	config.maxFactor = REAL(4.0);
	config.errorFloor = REAL(1e-6);

	StepSizeController controller(config);
	Real first = controller.acceptedStep(REAL(2.0), REAL(0.25));
	REQUIRE(first == Catch::Approx(REAL(2.0) * REAL(0.8) * std::pow(REAL(0.25), REAL(-0.25))));
	REQUIRE(controller.previousError() == Catch::Approx(REAL(0.25)));

	Real second = controller.acceptedStep(REAL(2.0), REAL(0.25));
	REQUIRE(second == Catch::Approx(first * std::pow(REAL(0.25), REAL(0.1))));

	controller.reset();
	REQUIRE(controller.acceptedStep(REAL(2.0), REAL(0.25)) == Catch::Approx(first));
	REQUIRE(controller.rejectedStep(REAL(-2.0), REAL(16.0)) < 0);

	controller.reset();
	Real orderSpecific = controller.acceptedStep(REAL(2.0), REAL(0.25), REAL(0.5));
	REQUIRE(orderSpecific == Catch::Approx(REAL(2.0) * REAL(0.8) * std::pow(REAL(0.25), REAL(-0.5))));
	REQUIRE_THROWS_AS(controller.rejectedStep(REAL(2.0), REAL(2.0), REAL(0.0)), ODESolverError);
}

TEST_CASE("Task_i5u3.3.2 - Hermite interpolator validates state and clamps both step directions", "[ODE][i5u3][interpolation]")
{
	HermiteInterpolator interpolator;
	REQUIRE_THROWS_AS(interpolator.interpolate(0), ODESolverError);

	Vector<Real> x0({ REAL(0.0) });
	Vector<Real> x1({ REAL(2.0) });
	Vector<Real> dx0({ REAL(1.0) });
	Vector<Real> dx1({ REAL(1.0) });
	interpolator.setStep(REAL(0.0), REAL(2.0), x0, x1, dx0, dx1);

	REQUIRE(interpolator.interpolate(REAL(1.0))[0] == Catch::Approx(REAL(1.0)));
	REQUIRE(interpolator.interpolate(REAL(-10.0))[0] == Catch::Approx(REAL(0.0)));
	REQUIRE(interpolator.interpolate(REAL(10.0))[0] == Catch::Approx(REAL(2.0)));

	interpolator.setStep(REAL(2.0), REAL(-2.0), x1, x0, dx1, dx0);
	REQUIRE(interpolator.interpolate(REAL(1.0))[0] == Catch::Approx(REAL(1.0)));
	REQUIRE(interpolator.interpolate(REAL(10.0))[0] == Catch::Approx(REAL(2.0)));
	REQUIRE(interpolator.interpolate(REAL(-10.0))[0] == Catch::Approx(REAL(0.0)));

	interpolator.reset();
	REQUIRE_THROWS_AS(interpolator.interpolate(REAL(1.0)), ODESolverError);
}

TEST_CASE("BugFix_i5u3.1.16 - serializer sanitizes newlines in titles", "[bugfix][i5u3][tools]")
{
	std::ostringstream out;
	auto res = Serializer::WriteRealFuncHeader(out, "REAL_FUNCTION", "evil\ntitle\rhere", 0.0, 1.0, 10);
	REQUIRE(res.success);

	// Header must still be exactly 6 lines - newlines in the title became spaces
	std::istringstream in(out.str());
	std::string line;
	int lines = 0;
	bool foundSanitized = false;
	while (std::getline(in, line)) {
		lines++;
		if (line == "evil title here")
			foundSanitized = true;
	}
	REQUIRE(lines == 6);
	REQUIRE(foundSanitized);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.2 — Contract policy conformance: base layer                          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Contract_i5u3.2.2 - Vector operator== returns false on size mismatch", "[bugfix][i5u3][contract]")
{
	Vector<Real> a({ 1.0, 2.0 });
	Vector<Real> b({ 1.0, 2.0, 3.0 });
	REQUIRE_NOTHROW(a == b);   // previously threw VectorDimensionError
	REQUIRE_FALSE(a == b);
	REQUIRE(a != b);
	REQUIRE(a == Vector<Real>({ 1.0, 2.0 }));
}

TEST_CASE("Contract_i5u3.2.2 - VectorN isZero is exact, isNearZero has tolerance", "[bugfix][i5u3][contract]")
{
	VectorN<Real, 3> tiny({ 1e-14, 0.0, 0.0 });
	REQUIRE_FALSE(tiny.isZero());       // previously true (hidden epsilon*100 tolerance)
	REQUIRE(tiny.isNearZero(REAL(1e-13)));

	VectorN<Real, 3> zero({ 0.0, 0.0, 0.0 });
	REQUIRE(zero.isZero());
	REQUIRE(zero.isNearZero());
}

TEST_CASE("Contract_i5u3.2.2 - Vec2/Vec3 normalization throws on near-zero", "[bugfix][i5u3][contract]")
{
	Vector2Cartesian z2(0.0, 0.0);
	Vector3Cartesian z3(0.0, 0.0, 0.0);

	REQUIRE_THROWS_AS(z2.Normalized(), VectorDimensionError);       // previously returned zero vector
	REQUIRE_THROWS_AS(z2.GetAsUnitVector(), VectorDimensionError);
	REQUIRE_THROWS_AS(z3.Normalized(), VectorDimensionError);
	REQUIRE_THROWS_AS(z3.GetAsUnitVector(), VectorDimensionError);

	// Explicit fallback variant preserves the old behavior under an honest name
	REQUIRE(z3.GetAsUnitVectorOrZero().isZero());

	Vector3Cartesian v(3.0, 0.0, 4.0);
	REQUIRE(v.Normalized().NormL2() == Catch::Approx(1.0));
	REQUIRE(v.GetAsUnitVectorOrZero().NormL2() == Catch::Approx(1.0));
}

TEST_CASE("Contract_i5u3.2.2 - DAESystem null callbacks throw NotImplementedError", "[bugfix][i5u3][contract]")
{
	DAESystem sys;   // default ctor leaves both callbacks null
	Vector<Real> x(2), y(1), out(2);
	REQUIRE_THROWS_AS(sys.diffEqs(0.0, x, y, out), NotImplementedError);      // previously silent no-op
	REQUIRE_THROWS_AS(sys.algConstraints(0.0, x, y, out), NotImplementedError);
}

TEST_CASE("Contract_i5u3.2.2 - SparseMatrix throws MML exceptions", "[bugfix][i5u3][contract]")
{
	REQUIRE_THROWS_AS(SparseMatrix::SparseMatrixCSR<Real>(-1, 3), MatrixDimensionError);

	SparseMatrix::SparseMatrixCSR<Real> m(2, 2);
	REQUIRE_THROWS_AS(m(5, 0), MatrixAccessBoundsError);

	std::vector<Real> wrongSize(3), out;
	REQUIRE_THROWS_AS(m.multiply(wrongSize, out), VectorDimensionError);
}

TEST_CASE("Contract_i5u3.2.2 - StandardFunctions throw MML DomainError", "[bugfix][i5u3][contract]")
{
	REQUIRE_THROWS_AS(Functions::Factorial(-1), DomainError);
	REQUIRE_THROWS_AS(Functions::Sec(Constants::PI / 2), DomainError);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.3 — Contract policy conformance: core and algorithms layers          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Contract_i5u3.2.3 - Fourier validation and arbitrary fast DCT", "[bugfix][i5u3][contract]")
{
	Vector<Real> empty(0);
	REQUIRE_THROWS_AS(Fourier::DCT::ForwardII(empty), FourierError);

	Vector<Real> notPow2({ 1.0, 2.0, 3.0 });
	auto fast = Fourier::DCT::ForwardII_Fast(notPow2);
	auto reference = Fourier::DCT::ForwardII_Reference(notPow2);
	for (int i = 0; i < notPow2.size(); i++)
		REQUIRE(std::abs(fast[i] - reference[i]) < Real(1e-5));
}

TEST_CASE("Contract_i5u3.2.3 - GetValues throws MML ArgumentError", "[bugfix][i5u3][contract]")
{
	RealFunctionFromStdFunc f([](Real x) { return x; });
	Vector<Real> outX, outY;
	REQUIRE_THROWS_AS(f.GetValues(0.0, 1.0, 1, outX, outY), ArgumentError);
}

TEST_CASE("Contract_i5u3.2.3 - LUSolver::Solve throws on dimension mismatch", "[bugfix][i5u3][contract]")
{
	Matrix<Real> A(2, 2, { 2.0, 1.0,
	                       1.0, 3.0 });
	LUSolver<Real> solver(A);

	Vector<Real> bWrong(3), xWrong(3);
	REQUIRE_THROWS_AS(solver.Solve(bWrong, xWrong), VectorDimensionError);   // previously returned false

	Vector<Real> b({ 3.0, 4.0 }), x(2);
	REQUIRE_NOTHROW(solver.Solve(b, x));
	REQUIRE(x[0] == Catch::Approx(1.0));
	REQUIRE(x[1] == Catch::Approx(1.0));
}

TEST_CASE("Contract_i5u3.2.3 - LP Steepest pivot rule throws NotImplementedError", "[bugfix][i5u3][contract]")
{
	Optimization::LinearProgram lp;
	lp.SetObjective(Vector<Real>({ 1.0, 1.0 }), Optimization::LPObjective::Maximize);
	lp.AddConstraint(Vector<Real>({ 1.0, 1.0 }), Optimization::LPConstraintType::LessEqual, 4.0);

	Optimization::LPConfig config;
	config.pivotRule = Optimization::LPPivotRule::Steepest;

	Optimization::SimplexSolver solver(config);
	REQUIRE_THROWS_AS(solver.Solve(lp), NotImplementedError);   // previously silently fell back to Dantzig
}

TEST_CASE("Contract_i5u3.2.3 - root finding verbose output goes to config stream", "[bugfix][i5u3][contract]")
{
	RealFunctionFromStdFunc f([](Real x) -> Real { return x * x - REAL(2.0); });

	std::ostringstream log;
	RootFinding::RootFindingConfig config;
	config.verbose = true;
	config.verboseStream = &log;

	auto result = RootFinding::FindRootBisection(f, 0.0, 2.0, config);
	REQUIRE(result.converged);
	REQUIRE(log.str().find("Bisection iter") != std::string::npos);
}

TEST_CASE("Contract_i5u3.2.3 - 2D adaptive integration rejects unimplemented GK rules", "[bugfix][i5u3][contract]")
{
	auto f = [](Real x, Real y) { return x + y; };
	Integration::AdaptiveConfig2D config;
	config.rule = Integration::GKRule::GK21;   // only tensor-GK15 is implemented
	REQUIRE_THROWS_AS(Integration::IntegrateAdaptive2D(f, 0.0, 1.0, 0.0, 1.0, config), NotImplementedError);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.4 — Reverse-time integration policy: validate and reject loudly      ///
///////////////////////////////////////////////////////////////////////////////////////////

namespace ReverseTimeDetail {
	// Index-1 linear DAE: x' = -x + y, constraint x + y = 1
	class SimpleDAE : public IODESystemDAEWithJacobian
	{
	public:
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& dxdt) const override { dxdt[0] = -x[0] + y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = x[0] + y[0] - 1.0; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1.0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1.0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1.0; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1.0; }
	};
}

TEST_CASE("Contract_i5u3.2.4 - DAE solvers reject reverse time and bad step size", "[bugfix][i5u3][contract]")
{
	ReverseTimeDetail::SimpleDAE sys;
	Vector<Real> x0({ 1.0 }), y0({ 0.0 });

	// Reverse time: previously silently returned an empty/undefined result
	REQUIRE_THROWS_AS(SolveDAEBDF2(sys, 1.0, x0, y0, 0.0), ArgumentError);
	REQUIRE_THROWS_AS(SolveDAEBackwardEuler(sys, 1.0, x0, y0, 0.0), ArgumentError);
	REQUIRE_THROWS_AS(SolveDAERODAS(sys, 1.0, x0, y0, 0.0), ArgumentError);

	// Non-positive step size
	DAESolverConfig badStep;
	badStep.step_size = 0.0;
	REQUIRE_THROWS_AS(SolveDAEBDF2(sys, 0.0, x0, y0, 1.0, badStep), ArgumentError);

	// Forward time still works
	auto result = SolveDAEBDF2(sys, 0.0, x0, y0, 0.5);
	REQUIRE(result.solution.getTotalSavedSteps() > 0);
}

TEST_CASE("Contract_i5u3.2.4 - stiff ODE solvers reject reverse time", "[bugfix][i5u3][contract]")
{
	// y' = -y with analytic Jacobian
	ODESystemWithJacobianFromStdFunc sys(1,
		[](Real, const Vector<Real>& y, Vector<Real>& dydt) { dydt[0] = -y[0]; },
		[](Real, const Vector<Real>& y, Vector<Real>& dydt, Matrix<Real>& J) { dydt[0] = -y[0]; J(0, 0) = -1.0; });
	Vector<Real> y0({ 1.0 });

	REQUIRE_THROWS_AS(SolveBackwardEuler(sys, 1.0, y0, 0.0, 0.1), ODESolverError);
	REQUIRE_THROWS_AS(SolveBDF2(sys, 1.0, y0, 0.0, 0.1), ODESolverError);
	REQUIRE_THROWS_AS(SolveBackwardEuler(sys, 0.0, y0, 1.0, -0.1), ODESolverError);
}

TEST_CASE("Contract_i5u3.2.4 - systems analyzers reject reverse time and bad steps", "[bugfix][i5u3][contract]")
{
	Systems::VanDerPolSystem vdp(1.0);
	Vector<Real> x({ 0.5, 0.5 });

	REQUIRE_THROWS_AS(Systems::PhaseSpaceAnalyzer::IntegrateTrajectory(vdp, x, -1.0, 0.1), ArgumentError);
	REQUIRE_THROWS_AS(Systems::PhaseSpaceAnalyzer::IntegrateTrajectory(vdp, x, 1.0, 0.1, -0.01), ArgumentError);
	REQUIRE_THROWS_AS(Systems::PhaseSpaceAnalyzer::ComputePoincareSection(vdp, x, Systems::PoincareSection<Real>(0, 0.0, +1), 10, -0.01), ArgumentError);

	// Forward time still works
	auto traj = Systems::PhaseSpaceAnalyzer::IntegrateTrajectory(vdp, x, 1.0, 0.1);
	REQUIRE(traj.size() > 1);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.5 — Precision-literal sweep                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Contract_i5u3.2.5 - named precision constants exist and are sane", "[bugfix][i5u3][contract]")
{
	STATIC_REQUIRE(Precision::PolynomialRootTolerance > 0);
	STATIC_REQUIRE(Precision::AGMConvergenceTolerance > 0);
	STATIC_REQUIRE(Precision::LyapunovChaosThreshold > 0);
	STATIC_REQUIRE(Precision::FixedPointTolerance > 0);
	STATIC_REQUIRE(Precision::FixedPointUniquenessTolerance > Precision::FixedPointTolerance);
	STATIC_REQUIRE(Systems::DefaultTransientTime > Systems::DefaultRecordTime);
}

TEST_CASE("Contract_i5u3.2.5 - singular condition numbers report infinity", "[bugfix][i5u3][contract]")
{
	Matrix<Real> singular(2, 2, { 1.0, 2.0,
	                              2.0, 4.0 });
	Systems::LinearSystem<Real> sys(singular);

	// previously returned the magic sentinel 1e30
	REQUIRE(std::isinf(sys.ConditionNumber1()));
	REQUIRE(std::isinf(sys.ConditionNumberInfinity()));
}

TEST_CASE("Contract_i5u3.2.5 - elliptic integrals still converge with named tolerance", "[bugfix][i5u3][contract]")
{
	// K(0) = pi/2, E(0) = pi/2
	REQUIRE(Functions::Comp_ellint_1(0.0) == Catch::Approx(Constants::PI / 2));
	REQUIRE(Functions::Comp_ellint_2(0.0) == Catch::Approx(Constants::PI / 2));
	// K(0.5) reference value
	REQUIRE(Functions::Comp_ellint_1(0.5) == Catch::Approx(1.6857503548).epsilon(EllipticIntegralTolerance));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.6 — Real-policy sweep: AD on Real, MML::AD namespace                 ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Contract_i5u3.2.6 - AD lives in MML::AD on Real", "[bugfix][i5u3][contract]")
{
	// Dual defaults to Real; ADVar/Tape store Real
	STATIC_REQUIRE(std::is_same_v<decltype(AD::Dual<>{}.value), Real>);
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<AD::ADVar>().value()), Real>);
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<AD::TapeEntry>().value), Real>);

	// Reverse-mode gradient of f(x,y) = x*y + sin(x) at (2, 3)
	auto grad = AD::gradientReverse(
		[](const std::vector<AD::ADVar>& v) { return v[0] * v[1] + AD::sin(v[0]); },
		std::vector<Real>{ 2.0, 3.0 });
	REQUIRE(grad[0] == Catch::Approx(3.0 + std::cos(2.0)));
	REQUIRE(grad[1] == Catch::Approx(2.0));

	// Backward-compat alias still works
	MML::AD::Dual<> d(1.0, 1.0);
	REQUIRE(d.value == REAL(1.0));
}

TEST_CASE("Contract_i5u3.2.6 - exception payloads carry Real", "[bugfix][i5u3][contract]")
{
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<SingularMatrixError>().determinant()), Real>);
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<ConvergenceError>().residual()), Real>);
	STATIC_REQUIRE(std::is_same_v<decltype(std::declval<IntegrationTooManySteps>().achieved_precision()), Real>);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.7 — Concepts on core containers and real-only solvers                ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Contract_i5u3.2.7 - containers constrained with Field, solvers with MMLReal", "[bugfix][i5u3][contract]")
{
	// Field admits the types the library actually uses...
	STATIC_REQUIRE(Field<Real>);
	STATIC_REQUIRE(Field<Complex>);
	STATIC_REQUIRE(Field<int>);
	STATIC_REQUIRE(Field<Algebra::PrimeFieldElement<5>>);
	// ...and rejects non-arithmetic types (Matrix/MatrixNM/MatrixSym now require Field,
	// so Matrix<std::string> fails with a one-line concept error instead of a template wall)
	STATIC_REQUIRE_FALSE(Field<std::string>);

	// QR/SVD are real-only per their pivoting/ordering semantics (MMLReal);
	// complex QR/SVD would silently miscompute - now rejected at instantiation
	STATIC_REQUIRE(MMLReal<Real>);
	STATIC_REQUIRE_FALSE(MMLReal<Complex>);

	// Positive instantiations still work
	STATIC_REQUIRE(requires { typename Matrix<Real>; });
	STATIC_REQUIRE(requires { typename Matrix<Complex>; });
	STATIC_REQUIRE(requires { typename MatrixNM<Algebra::PrimeFieldElement<5>, 2, 2>; });
	STATIC_REQUIRE(requires { typename QRSolver<Real>; });
	STATIC_REQUIRE(requires { typename SVDecompositionSolver<Real>; });
	// LU stays generic - complex LU is legitimate and used
	STATIC_REQUIRE(requires { typename LUSolver<Complex>; });
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.2.8 — All raw std:: exceptions converted to MML taxonomy               ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Contract_i5u3.2.8 - converted throws are MML exceptions with std bases intact", "[bugfix][i5u3][contract]")
{
	// Spot-checks across the converted areas; catch as MMLException proves taxonomy membership,
	// catch as the std base proves existing std-typed handlers keep working
	REQUIRE_THROWS_AS(Curves::LogSpiralCurve(1.0), MMLException);
	REQUIRE_THROWS_AS(Curves::LogSpiralCurve(1.0), std::invalid_argument);

	ChebyshevBasis basis;
	REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.5), MMLException);

	ChebyshevApproximation cheb([](Real x) { return x * x; }, 0.0, 1.0, 8);
	REQUIRE_THROWS_AS(cheb(2.0), DomainError);

	// InvalidStateError: new taxonomy member for state misuse
	Timer t;
	REQUIRE_THROWS_AS(t.GetIntervalTime(0), InvalidStateError);
	REQUIRE_THROWS_AS(t.GetIntervalTime(0), std::runtime_error);   // std base preserved

	// ODESystem null callback now NotImplementedError per policy section 5
	ODESystem emptySys;
	Vector<Real> x(1), dxdt(1);
	REQUIRE_THROWS_AS(emptySys.derivs(0.0, x, dxdt), NotImplementedError);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Task i5u3.3.7 — One shared RFC 4180 CSV module (CsvUtils)                        ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Dedup_i5u3.3.7 - CsvUtils escape/unescape/split roundtrip", "[bugfix][i5u3][dedup]")
{
	using namespace CsvUtils;

	REQUIRE(Escape("plain") == "plain");
	REQUIRE(Escape("a,b") == "\"a,b\"");
	REQUIRE(Escape("He said \"Hi\"") == "\"He said \"\"Hi\"\"\"");
	REQUIRE(Unescape(Escape("quote\"comma,nl\n")) == "quote\"comma,nl\n");

	auto fields = SplitLine("a,\"b,c\",\"d\"\"e\"", ',');
	REQUIRE(fields.size() == 3);
	REQUIRE(fields[0] == "a");
	REQUIRE(fields[1] == "b,c");
	REQUIRE(fields[2] == "d\"e");

	// The three former call sites all route here now
	REQUIRE(Csv::escape("a,b") == Escape("a,b"));
}

} // namespace MML::Tests::BugfixI5u3
