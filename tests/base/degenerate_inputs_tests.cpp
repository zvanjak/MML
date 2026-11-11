#ifndef MML_USE_SINGLE_HEADER
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#endif
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        degenerate_inputs_tests.cpp                                         ///
///  Description: Degenerate-input sweep across base: empty / zero-size containers,   ///
///               N = 0/1/2 point counts, zero vectors in normalisation, OOB indices  ///
///               in every matrix variant, null-polynomial corner cases.              ///
///                                                                                   ///
///  Analysis v4, MinimalMathLibrary-wewi.2.2. Each case either pins already-correct  ///
///  behaviour or guards a defect found by the sweep (noted inline).                  ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/BaseUtils.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/base/Matrix/Matrix3D.h>
#include <mml/base/Polynom.h>
#include <mml/base/InterpolatedFunction.h>
#include <mml/base/ChebyshevApproximation.h>
#include <mml/base/Intervals.h>

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <cmath>
#include <functional>

using namespace MML;

namespace MML::Tests::Base::DegenerateInputsTests {

///////////////////////////////////////////////////////////////////////////////////////////
//                               Empty / zero-size containers
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Empty Vector: norms, predicates and arithmetic are total", "[degenerate][vector][wewi.2.2]")
{
    Vector<Real> e;
    REQUIRE(e.size() == 0);
    REQUIRE(e.NormL2() == static_cast<MML::Real>(0.0));
    REQUIRE(e.NormL1() == static_cast<MML::Real>(0.0));
    REQUIRE(e.NormLInf() == static_cast<MML::Real>(0.0));
    REQUIRE(e.IsZero());
    REQUIRE(e.IsEqualTo(e));
    REQUIRE((e + e).size() == 0);
    REQUIRE((e * static_cast<MML::Real>(2.0)).size() == 0);
    REQUIRE(Utils::ScalarProduct(e, e) == static_cast<MML::Real>(0.0));
    REQUIRE(e.to_string(8, 4) == "[]");
    REQUIRE(Vector<Real>(0).size() == 0);
    REQUIRE_THROWS_AS(Vector<Real>(-1), VectorInitializationError);
}

TEST_CASE("Empty / zero-size matrices behave consistently", "[degenerate][matrix][wewi.2.2]")
{
    Matrix<Real> m0;
    REQUIRE(m0.rows() == 0);
    REQUIRE(m0.IsSquare());
    REQUIRE((m0 * m0).rows() == 0);
    REQUIRE(m0.Transposed().rows() == 0);
    REQUIRE(m0.Trace() == static_cast<MML::Real>(0.0));
    REQUIRE(m0.IsEqualTo(m0));
    REQUIRE_FALSE(m0.to_string(8, 4).empty());

    REQUIRE(Matrix<Real>(0, 5).rows() == 0);
    REQUIRE((Matrix<Real>(0, 5) * Vector<Real>(5)).size() == 0);
    Matrix<Real> p = Matrix<Real>(2, 0) * Matrix<Real>(0, 3);   // 2x0 * 0x3 = 2x3 zero matrix
    REQUIRE(p.rows() == 2);
    REQUIRE(p.cols() == 3);
    REQUIRE(p.IsZero());

    REQUIRE_THROWS_AS(Matrix<Real>(-1, 2), MatrixDimensionError);
    REQUIRE_THROWS_AS(([] { Matrix<Real> t; t.resize(-1, -1); }()), MatrixDimensionError);
    REQUIRE_THROWS_AS(Matrix<Real>::Identity(0), MatrixDimensionError);
    REQUIRE_THROWS_AS(Matrix<Real>::Identity(-1), MatrixDimensionError);

    REQUIRE(MatrixSym<Real>(0).rows() == 0);
    REQUIRE(MatrixSym<Real>(0).Trace() == static_cast<MML::Real>(0.0));
    REQUIRE(MatrixSym<Real>(0).NormFrobenius() == static_cast<MML::Real>(0.0));
    REQUIRE_THROWS_AS(MatrixSym<Real>(-1), MatrixDimensionError);

    REQUIRE_THROWS_AS(Matrix3D<Real>(0, 0, 0), ArgumentError);
    REQUIRE_THROWS_AS(Matrix3D<Real>(-1, 1, 1), ArgumentError);
}

TEST_CASE("Singular matrices are rejected by inverse() including the all-zero case (regression)", "[degenerate][matrix][wewi.2.2]")
{
    // The pivot test was |pivot| < eps*norm*n; for an all-zero matrix norm = 0, so 0 < 0 was
    // false and the inverse silently became Inf/NaN.
    REQUIRE_THROWS_AS((MatrixAlg::Inverse(Matrix<Real>(1, 1, {static_cast<MML::Real>(0.0)}))), SingularMatrixError);
    REQUIRE_THROWS_AS(MatrixAlg::Inverse(Matrix<Real>(2, 2)), SingularMatrixError);
    REQUIRE_THROWS_AS(([] { MatrixNM<Real, 2, 2> z; return MatrixAlg::Inverse(z); }()), SingularMatrixError);
    REQUIRE(MatrixAlg::Inverse(Matrix<Real>(1, 1, {static_cast<MML::Real>(4.0)}))(0, 0) == Catch::Approx(0.25));
    REQUIRE_THROWS_AS(MatrixAlg::Inverse(Matrix<Real>(2, 3)), MatrixDimensionError);
    REQUIRE_THROWS_AS(Matrix<Real>(2, 3).Trace(), MatrixDimensionError);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Out-of-bounds indices
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Checked element access throws for OOB indices in every container variant", "[degenerate][oob][wewi.2.2]")
{
    Vector<Real> v{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};
		REQUIRE_THROWS_AS(v.at(-1), VectorAccessBoundsError);
		REQUIRE_THROWS_AS(v.at(2), VectorAccessBoundsError);
		REQUIRE_THROWS_AS(Vector<Real>::UnitVector(0, 0), VectorDimensionError);
    REQUIRE_THROWS_AS(Vector<Real>::UnitVector(3, -1), VectorDimensionError);

    VectorN<Real, 3> vn;
		REQUIRE_THROWS_AS(vn.at(3), VectorAccessBoundsError);
		REQUIRE_THROWS_AS(vn.at(-1), VectorAccessBoundsError);

		Matrix<Real> m(2, 2, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)});
    REQUIRE_THROWS_AS(m.at(2, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(m.at(0, -1), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(m.at(-1, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(m.GetSubmatrix(5, 5, 1, 1), MatrixDimensionError);
    REQUIRE_THROWS_AS(m * Vector<Real>(3), MatrixDimensionError);
    REQUIRE_THROWS_AS(m + Matrix<Real>(3, 3), MatrixDimensionError);
    REQUIRE_THROWS_AS(m * Matrix<Real>(3, 3), MatrixDimensionError);

    MatrixNM<Real, 2, 2> nm;
    REQUIRE_THROWS_AS(nm.at(2, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(nm.at(0, -1), MatrixAccessBoundsError);

    MatrixSym<Real> sy(2);
    REQUIRE_THROWS_AS(sy.at(2, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(sy.at(0, -1), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(sy.at(-1, 1), MatrixAccessBoundsError);

    Matrix3D<Real> m3(2, 2, 2);
    REQUIRE_THROWS_AS(m3.at(2, 0, 0), IndexError);   // Matrix3D spells it At() (naming drift: wewi.5)
    REQUIRE_THROWS_AS(m3.at(0, -1, 0), IndexError);

    // TridiagonalMatrix::operator() is the checked accessor; negative / too-large indices used
    // to read memory before the band vectors (const overload) instead of throwing.
    TridiagonalMatrix<Real> t3(3, {static_cast<MML::Real>(2), static_cast<MML::Real>(1), static_cast<MML::Real>(1), static_cast<MML::Real>(2), static_cast<MML::Real>(1), static_cast<MML::Real>(1), static_cast<MML::Real>(2)});
    const TridiagonalMatrix<Real>& ct3 = t3;
    REQUIRE_THROWS_AS(t3(3, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(t3(-1, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(t3(0, 2), MatrixAccessBoundsError);       // off-band not writable
    REQUIRE_THROWS_AS(ct3(-1, 0), MatrixAccessBoundsError);
    REQUIRE_THROWS_AS(ct3(0, 3), MatrixAccessBoundsError);
    REQUIRE(ct3(0, 2) == static_cast<MML::Real>(0.0));                             // off-band readable as 0
    REQUIRE(ct3(2, 1) == static_cast<MML::Real>(1.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               N = 1 / 2 sized objects
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("TridiagonalMatrix supports dimensions 1 and 2 (regression)", "[degenerate][tridiag][wewi.2.2]")
{
    // The initializer-list constructor read past the list for dim == 1 and then overwrote the
    // diagonal, so TridiagonalMatrix(1, {5}) had (0,0) == 0 and Solve reported a zero pivot.
    TridiagonalMatrix<Real> t1(1, {static_cast<MML::Real>(5.0)});
    REQUIRE(t1(0, 0) == static_cast<MML::Real>(5.0));
    REQUIRE(t1.Solve(Vector<Real>{static_cast<MML::Real>(10.0)})[0] == Catch::Approx(2.0));

    TridiagonalMatrix<Real> t2(2, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
    Vector<Real> s = t2.Solve(Vector<Real>{static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0)});
    REQUIRE(s[0] == Catch::Approx(1.0));
    REQUIRE(s[1] == Catch::Approx(1.0));

    Vector<Real> a{static_cast<MML::Real>(0.0)}, d{static_cast<MML::Real>(3.0)}, c{static_cast<MML::Real>(0.0)};
    TridiagonalMatrix<Real> tv(1, a, d, c);
    REQUIRE(tv.Solve(Vector<Real>{static_cast<MML::Real>(6.0)})[0] == Catch::Approx(2.0));

    REQUIRE_THROWS_AS(TridiagonalMatrix<Real>(0, {}), MatrixDimensionError);
    REQUIRE_THROWS_AS(TridiagonalMatrix<Real>(-1, {}), MatrixDimensionError);
    REQUIRE_THROWS_AS((TridiagonalMatrix<Real>(2, {static_cast<MML::Real>(1.0)})), MatrixDimensionError);
    REQUIRE_THROWS_AS((TridiagonalMatrix<Real>(2, {static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}).Solve(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)})),
                      SingularMatrixError);
}

TEST_CASE("Interpolation rejects 0/1 points, size mismatch and non-monotonic abscissas", "[degenerate][interpolation][wewi.2.2]")
{
    Vector<Real> x0, y0, x1{static_cast<MML::Real>(1.0)}, y1{static_cast<MML::Real>(2.0)}, x2{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)}, y2{static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0)};

    REQUIRE_THROWS_AS(LinearInterpRealFunc(x0, y0), RealFuncInterpInitError);
    REQUIRE_THROWS_AS(LinearInterpRealFunc(x1, y1), RealFuncInterpInitError);
    REQUIRE(LinearInterpRealFunc(x2, y2)(static_cast<MML::Real>(1.5)) == Catch::Approx(3.0));

    // size mismatch used to be accepted and read y out of bounds
    REQUIRE_THROWS_AS(LinearInterpRealFunc(x2, y1), RealFuncInterpInitError);
    REQUIRE_THROWS_AS(SplineInterpRealFunc(x2, y1), RealFuncInterpInitError);

    // non-monotonic x used to be accepted and produce silent garbage (Linear returned 0)
    Vector<Real> xu{static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0)}, yu{static_cast<MML::Real>(1.0), static_cast<MML::Real>(9.0), static_cast<MML::Real>(4.0)};
    REQUIRE_THROWS_AS(LinearInterpRealFunc(xu, yu), RealFuncInterpInitError);
    REQUIRE_THROWS_AS(SplineInterpRealFunc(xu, yu), RealFuncInterpInitError);
    REQUIRE_THROWS_AS(PolynomInterpRealFunc(xu, yu, 3), RealFuncInterpInitError);
    // duplicates still reported
    REQUIRE_THROWS_AS((LinearInterpRealFunc(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)}, y2)), RealFuncInterpInitError);
    // strictly descending is allowed (NR locate/hunt handle both directions)
    LinearInterpRealFunc desc(Vector<Real>{static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0)}, Vector<Real>{static_cast<MML::Real>(9.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(1.0)});
    REQUIRE(desc(static_cast<MML::Real>(2.5)) == Catch::Approx(6.5));

    REQUIRE_THROWS_AS(PolynomInterpRealFunc(x2, y2, 3), RealFuncInterpInitError);   // m > n
    REQUIRE_THROWS_AS(PolynomInterpRealFunc(x2, y2, 0), RealFuncInterpInitError);
    REQUIRE(SplineInterpRealFunc(x2, y2)(static_cast<MML::Real>(1.5)) == Catch::Approx(3.0));         // 2 points = linear
    REQUIRE_THROWS_AS(BarycentricPolynomialInterp(x1, y1), RealFuncInterpInitError);
    REQUIRE(BarycentricPolynomialInterp(x2, y2)(static_cast<MML::Real>(1.5)) == Catch::Approx(3.0));
    REQUIRE_THROWS_AS(BarycentricRationalInterp(x2, y2, 3), RealFuncInterpInitError);
    REQUIRE_THROWS_AS(MonotoneCubicInterpRealFunc(x1, y1), RealFuncInterpInitError);
    REQUIRE(MonotoneCubicInterpRealFunc(x2, y2)(static_cast<MML::Real>(1.5)) == Catch::Approx(3.0));
}

TEST_CASE("ChebyshevApproximation validates n and the interval before allocating", "[degenerate][chebyshev][wewi.2.2]")
{
    std::function<Real(Real)> f = [](Real v) { return std::sin(v); };
    REQUIRE_THROWS_AS(ChebyshevApproximation(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 0), ArgumentError);
    REQUIRE_THROWS_AS(ChebyshevApproximation(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), -1), ArgumentError);   // used to leak VectorInitializationError
    REQUIRE_THROWS_AS(ChebyshevApproximation(f, static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), 10), ArgumentError);
    REQUIRE_THROWS_AS(ChebyshevApproximation(f, static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), 10), ArgumentError);
    // n = 1 is the constant approximation f(midpoint)
    REQUIRE(ChebyshevApproximation(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 1)(static_cast<MML::Real>(0.5)) == Catch::Approx(std::sin(0.5)));
}

TEST_CASE("Intervals: point coverings with 0/1/2 points and reversed bounds", "[degenerate][intervals][wewi.2.2]")
{
    ClosedInterval ci(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
    std::vector<Real> pts;
    ci.getEquidistantCovering(0, pts);
    REQUIRE(pts.empty());
    ci.getEquidistantCovering(1, pts);
    REQUIRE(pts.size() == 1);
    REQUIRE(pts[0] == Catch::Approx(0.5));
    ci.getEquidistantCovering(2, pts);
    REQUIRE(pts.size() == 2);
    REQUIRE(pts.front() == Catch::Approx(0.0));
    REQUIRE(pts.back() == Catch::Approx(1.0));

    ClosedInterval degenerate(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0));
    degenerate.getEquidistantCovering(3, pts);
    REQUIRE(pts.size() == 3);
    for (Real p : pts) REQUIRE(p == static_cast<MML::Real>(1.0));

    REQUIRE_THROWS_AS(ClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)), ArgumentError);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Zero vectors in normalisation / angles
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Zero vectors: normalisation, angle and projection throw instead of producing NaN", "[degenerate][vector][wewi.2.2]")
{
    Vector<Real> z{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, x{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};
    REQUIRE_THROWS_AS(Utils::VectorsAngle(z, x), ArgumentError);          // used to return NaN
    REQUIRE_THROWS_AS(Utils::VectorsAngle(x, z), ArgumentError);
    REQUIRE_THROWS_AS(Utils::VectorProjectionParallelTo(x, z), ArgumentError);
    REQUIRE_THROWS_AS(Utils::VectorProjectionPerpendicularTo(x, z), ArgumentError);
    REQUIRE_THROWS_AS((Utils::VectorProjectionParallelTo(x, Vector<Real>{static_cast<MML::Real>(1.0)})), VectorDimensionError);

    VectorN<Real, 3> zn, en{static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)};
    REQUIRE_THROWS_AS(zn.Normalized(), VectorDimensionError);
    REQUIRE_THROWS_AS(Utils::VectorsAngle(zn, en), ArgumentError);
    REQUIRE_THROWS_AS(Vector3Cartesian(0, 0, 0).Normalized(), VectorDimensionError);
    REQUIRE_THROWS_AS(Vector3Cartesian(0, 0, 0).Normalized(), VectorDimensionError);
}

TEST_CASE("VectorsAngle is exact for parallel, antiparallel and orthogonal vectors", "[degenerate][vector][wewi.2.2]")
{
    // acos(dot/(|a||b|)) returned ~2e-8 for identical vectors; Kahan's atan2 formula is exact
    Vector<Real> a{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};
    REQUIRE(Utils::VectorsAngle(a, a) == static_cast<MML::Real>(0.0));
    REQUIRE(Utils::VectorsAngle(a, a * static_cast<MML::Real>(3.0)) == static_cast<MML::Real>(0.0));
    REQUIRE(Utils::VectorsAngle(a, a * static_cast<MML::Real>(-1.0)) == Catch::Approx(Constants::PI));
    REQUIRE(Utils::VectorsAngle(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}, Vector<Real>{static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)}) == Catch::Approx(Constants::PI / 2));
    REQUIRE(Utils::VectorsAngle(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}, Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)}) == Catch::Approx(Constants::PI / 4).epsilon(1e-15));

    REQUIRE(Utils::VectorsAngle(Vector3Cartesian(1, 2, 3), Vector3Cartesian(2, 4, 6)) == static_cast<MML::Real>(0.0));
    REQUIRE(Utils::VectorsAngle(Vector3Cartesian(1, 2, 3), Vector3Cartesian(-2, -4, -6)) == Catch::Approx(Constants::PI));
    REQUIRE(Utils::VectorsAngle(Vector3Cartesian(1, 0, 0), Vector3Cartesian(0, 1, 0)) == Catch::Approx(Constants::PI / 2));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Null / constant polynomials
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Null and constant polynomial corner cases are total (regression)", "[degenerate][polynom][wewi.2.2]")
{
    Polynom<Real> null, c{static_cast<MML::Real>(3.0)}, lin{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};
    REQUIRE(null(static_cast<MML::Real>(2.0)) == static_cast<MML::Real>(0.0));
    REQUIRE(null.derivative().isNull());
    REQUIRE(c.derivative()(static_cast<MML::Real>(1.0)) == static_cast<MML::Real>(0.0));

    // poldiv(null, lin) wrote r[0] of an empty remainder -> access violation
    Polynom<Real> q, r;
    REQUIRE_NOTHROW(Polynom<Real>::poldiv(null, lin, q, r));
    REQUIRE(q.isNull());
    REQUIRE(r.isNull());
    REQUIRE_NOTHROW(Polynom<Real>::poldiv(c, lin, q, r));       // deg(u) < deg(v)
    REQUIRE(q.isNull());
    REQUIRE(r.degree() == 0);
    REQUIRE(r[0] == static_cast<MML::Real>(3.0));
    REQUIRE_THROWS_AS(Polynom<Real>::poldiv(lin, null, q, r), DivisionByZeroError);

    // Derive with an empty output vector wrote pd[0] out of bounds
    Vector<Real> empty;
    REQUIRE_THROWS_AS(lin.Derive(static_cast<MML::Real>(1.0), empty), VectorDimensionError);
    Vector<Real> jet(3);
    REQUIRE_NOTHROW(null.Derive(static_cast<MML::Real>(1.0), jet));
    REQUIRE(jet[0] == static_cast<MML::Real>(0.0));
    Vector<Real> wide(5);                                        // more derivatives than the degree
    lin.Derive(static_cast<MML::Real>(1.0), wide);
    REQUIRE(wide[0] == Catch::Approx(3.0));
    REQUIRE(wide[1] == Catch::Approx(2.0));
    REQUIRE(wide[2] == static_cast<MML::Real>(0.0));
    REQUIRE(wide[4] == static_cast<MML::Real>(0.0));

    // SetDegree(-5) leaked std::length_error
    Polynom<Real> p{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};
    REQUIRE_NOTHROW(p.SetDegree(-1));
    REQUIRE(p.isNull());
    REQUIRE_THROWS_AS(p.SetDegree(-5), ArgumentError);

    REQUIRE(Polynom<Real>::FromValues({static_cast<MML::Real>(1.0)}, {static_cast<MML::Real>(2.0)}).degree() == 0);
    REQUIRE(Polynom<Real>::FromValues({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)}, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0)})(static_cast<MML::Real>(3.0)) == Catch::Approx(6.0));
}

} // namespace MML::Tests::Base::DegenerateInputsTests
