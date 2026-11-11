#ifndef MML_USE_SINGLE_HEADER
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#endif
///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        complex_instantiation_tests.cpp                                     ///
///  Description: Instantiate public templates with Complex and exercise the members  ///
///               that previously only compiled for Real (analysis v4, BASE D5-D7).  ///
///                                                                                   ///
///  Started with MinimalMathLibrary-wewi.1.6; the full sweep is wewi.2.1.           ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Polynom.h>

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include <sstream>

using namespace MML;

namespace MML::Tests::Base::ComplexInstantiationTests {

TEST_CASE("VectorN<Complex>::IsEqualTo compiles and compares magnitudes", "[complex-instantiation][vectorn][wewi.1.6]")
{
    // eps used to be `Type eps`, so Abs(difference) > eps was Real > Complex and did not compile
    VectorN<Complex, 2> a{Complex(1.0, 2.0), Complex(-3.0, 0.5)};
    VectorN<Complex, 2> b{Complex(1.0, 2.0 + 1e-12), Complex(-3.0, 0.5)};
    VectorN<Complex, 2> c{Complex(1.0, 2.5), Complex(-3.0, 0.5)};

    REQUIRE(a.IsEqualTo(b));
    REQUIRE(a.IsEqualTo(b, static_cast<MML::Real>(1e-9)));
    REQUIRE_FALSE(a.IsEqualTo(c));
    REQUIRE(a.IsEqualTo(c, static_cast<MML::Real>(1.0)));
}

TEST_CASE("Matrix<Complex>::Print with threshold compiles and zeroes small entries by magnitude", "[complex-instantiation][matrix][wewi.1.6]")
{
    // static_cast<Real>(complex) was instantiated because is_MML_simple_numeric<Complex> is true
    Matrix<Complex> m(2, 2);
    m(0, 0) = Complex(1.0, 0.0);
    m(0, 1) = Complex(1e-14, -1e-14);
    m(1, 0) = Complex(0.0, 2.0);
    m(1, 1) = Complex(-1e-13, 0.0);

    std::ostringstream out;
    m.Print(out, 8, 3, static_cast<MML::Real>(1e-10));
    const std::string text = out.str();
    INFO(text);
    REQUIRE_FALSE(text.empty());
    REQUIRE(text.find("1e-14") == std::string::npos);
    REQUIRE(text.find("1e-13") == std::string::npos);
}

TEST_CASE("PolynomComplex::Reduce compiles and strips trailing near-zero coefficients", "[complex-instantiation][polynom][wewi.1.6]")
{
    // std::numeric_limits<Complex>::epsilon() is 0 and std::abs(c) < Complex did not compile
    PolynomComplex p{Complex(1.0, 1.0), Complex(2.0, 0.0), Complex(1e-15, -1e-15), Complex(0.0, 1e-16)};
    REQUIRE(p.degree() == 3);
    p.Reduce();
    REQUIRE(p.degree() == 1);
    REQUIRE(p[1] == Complex(2.0, 0.0));

    PolynomComplex q{Complex(0.5, 0.5), Complex(0.0, 0.01)};
    q.Reduce(static_cast<MML::Real>(0.1));
    REQUIRE(q.degree() == 0);
}

TEST_CASE("Polynom: arithmetic with the null polynomial is total", "[polynom][null][wewi.1.6]")
{
    // size() + size() - 1 underflowed for two empty coefficient vectors -> std::length_error
    Polynom<Real> null1, null2;
    Polynom<Real> p{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};

    REQUIRE(null1.isNull());
    REQUIRE_NOTHROW(null1 * null2);
    REQUIRE((null1 * null2).isNull());
    REQUIRE((null1 * p).isNull());
    REQUIRE((p * null1).isNull());

    Polynom<Real> sum = null1 + p;
    REQUIRE(sum.degree() == 1);
    REQUIRE(sum[1] == static_cast<MML::Real>(2.0));
    REQUIRE((null1 + null2).isNull());

    PolynomComplex cnull;
    REQUIRE_NOTHROW(cnull * cnull);
    REQUIRE((cnull * cnull).isNull());
}



///////////////////////////////////////////////////////////////////////////////////////////
///  wewi.2.1 - full Complex instantiation sweep
///
///  Explicit instantiation forces every non-template member of each public class template
///  to compile with Complex; the TEST_CASEs below exercise the members that were fixed by
///  this sweep plus a representative set of free operators. Serializers and VectorTypes2D/3D
///  are Real-only by design and are intentionally not covered here.
///////////////////////////////////////////////////////////////////////////////////////////

} // namespace MML::Tests::Base::ComplexInstantiationTests

#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/base/Matrix/Matrix3D.h>
#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/base/Tensor/MatrixNDim.h>

template class MML::Vector<Complex>;
template class MML::VectorN<Complex, 3>;
template class MML::Matrix<Complex>;
template class MML::MatrixNM<Complex, 3, 3>;
template class MML::MatrixSym<Complex>;
template class MML::TridiagonalMatrix<Complex>;
template class MML::Matrix3D<Complex>;
template class MML::Polynom<Complex>;
template class MML::Polynom<Real, Complex>;
template class MML::SparseMatrix::SparseMatrixCOO<Complex>;
template class MML::SparseMatrix::SparseMatrixCSR<Complex>;
template class MML::SparseMatrix::SparseMatrixCSC<Complex>;
template class MML::MatrixNDim<Complex, 3, 2>;

namespace MML::Tests::Base::ComplexInstantiationTests {

namespace {
    const Complex I(0.0, 1.0);
}

TEST_CASE("VectorN<Complex>::AreEqual compiles with a real tolerance", "[complex-instantiation][vectorn][wewi.2.1]")
{
    VectorN<Complex, 2> a{Complex(1.0, 2.0), Complex(-3.0, 0.5)};
    VectorN<Complex, 2> b{Complex(1.0, 2.0 + 1e-12), Complex(-3.0, 0.5)};
    REQUIRE(VectorN<Complex, 2>::AreEqual(a, b));
    REQUIRE_FALSE(VectorN<Complex, 2>::AreEqual(a, b, static_cast<MML::Real>(1e-15)));
}

TEST_CASE("MatrixNM<Complex>::Invert / SolveGaussJordan pivot on magnitude", "[complex-instantiation][matrixnm][wewi.2.1]")
{
    // `Type big` made std::abs(..) >= big a Complex comparison; now Real
    MatrixNM<Complex, 2, 2> m;
    m(0, 0) = Complex(0.0, 0.0); m(0, 1) = Complex(0.0, 1.0);
    m(1, 0) = Complex(2.0, 0.0); m(1, 1) = Complex(1.0, 1.0);
    // determinant = 0*(1+i) - i*2 = -2i  -> invertible, but needs a row swap (zero pivot)
    auto inv = MatrixAlg::Inverse(m);
    auto prod = m * inv;
    REQUIRE(std::abs(prod(0, 0) - Complex(1.0)) < 1e-14);
    REQUIRE(std::abs(prod(0, 1)) < 1e-14);
    REQUIRE(std::abs(prod(1, 0)) < 1e-14);
    REQUIRE(std::abs(prod(1, 1) - Complex(1.0)) < 1e-14);
}

TEST_CASE("TridiagonalMatrix<Complex>: IsEqualTo and Solve", "[complex-instantiation][tridiag][wewi.2.1]")
{
    TridiagonalMatrix<Complex> t(3, {Complex(2.0, 1.0), Complex(-1.0, 0.0),
                                     Complex(-1.0, 0.0), Complex(2.0, -1.0), Complex(-1.0, 0.0),
                                     Complex(-1.0, 0.0), Complex(2.0, 0.5)});
    TridiagonalMatrix<Complex> u = t;
    REQUIRE(t.IsEqualTo(u));
    REQUIRE(t.IsEqualTo(u, static_cast<MML::Real>(1e-15)));

    Vector<Complex> xExact{Complex(1.0, 1.0), Complex(-2.0, 0.5), Complex(0.0, 3.0)};
    Vector<Complex> rhs(3);
    for (int i = 0; i < 3; ++i) {
        rhs[i] = Complex(0.0);
        for (int j = 0; j < 3; ++j)
            if (std::abs(i - j) <= 1) rhs[i] += t(i, j) * xExact[j];
    }
    Vector<Complex> x = t.Solve(rhs);
    REQUIRE(x.IsEqualTo(xExact, static_cast<MML::Real>(1e-12)));
}

TEST_CASE("PolynomComplex: IsEqualTo exists and Derive produces complex jets", "[complex-instantiation][polynom][wewi.2.1]")
{
    // IsEqualTo was SFINAE-restricted to floating-point CoefT -> no such member on PolynomComplex
    PolynomComplex p{Complex(1.0, 1.0), Complex(0.0, 2.0), Complex(3.0, 0.0)};   // (1+i) + 2i x + 3 x^2
    PolynomComplex q{Complex(1.0, 1.0 + 1e-13), Complex(0.0, 2.0), Complex(3.0, 0.0)};
    REQUIRE(p.IsEqualTo(q));
    REQUIRE_FALSE(p.IsEqualTo(q, static_cast<MML::Real>(1e-15)));

    // Derive took Vector<Real> and could not store complex coefficients
    Vector<Complex> jet(3);
    const Complex x(0.5, -1.0);
    p.Derive(x, jet);
    REQUIRE(std::abs(jet[0] - (Complex(1.0, 1.0) + Complex(0.0, 2.0) * x + Complex(3.0) * x * x)) < 1e-14);
    REQUIRE(std::abs(jet[1] - (Complex(0.0, 2.0) + Complex(6.0) * x)) < 1e-14);
    REQUIRE(std::abs(jet[2] - Complex(6.0)) < 1e-14);

    // Real-coefficient path unchanged
    Polynom<Real> r{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)};
    Vector<Real> rjet(3);
    r.Derive(static_cast<MML::Real>(2.0), rjet);
    REQUIRE(rjet[0] == Catch::Approx(17.0));
    REQUIRE(rjet[1] == Catch::Approx(14.0));
    REQUIRE(rjet[2] == Catch::Approx(6.0));
}

TEST_CASE("SparseMatrixCSR<Complex>: isSymmetric, norms, norm2 of complex vectors, triangular solves",
          "[complex-instantiation][sparse][wewi.2.1]")
{
    using SparseMatrix::SparseMatrixCOO;
    using SparseMatrix::SparseMatrixCSR;

    SparseMatrixCOO<Complex> coo(3, 3);
    coo.addEntry(0, 0, Complex(2.0, 0.0));
    coo.addEntry(0, 1, Complex(1.0, 1.0));
    coo.addEntry(1, 0, Complex(1.0, 1.0));    // symmetric (not Hermitian) off-diagonal
    coo.addEntry(1, 1, Complex(0.0, 3.0));
    coo.addEntry(2, 2, Complex(-1.0, 0.0));
    SparseMatrixCSR<Complex> A(coo);

    SECTION("isSymmetric takes a real tolerance")
    {
        REQUIRE(A.IsSymmetric());
        REQUIRE(A.IsSymmetric(static_cast<MML::Real>(1e-14)));
        SparseMatrixCOO<Complex> coo2(2, 2);
        coo2.addEntry(0, 1, Complex(1.0, 1.0));
        coo2.addEntry(1, 0, Complex(1.0, -1.0));   // Hermitian but not symmetric
        REQUIRE_FALSE(SparseMatrixCSR<Complex>(coo2).IsSymmetric());
    }

    SECTION("norm2 is sqrt(sum |x|^2) and returns a Real")
    {
        std::vector<Complex> v{Complex(3.0, 4.0), Complex(0.0, 0.0)};
        static_assert(std::is_same_v<decltype(SparseMatrix::norm2(v)), Real>);
        REQUIRE(SparseMatrix::norm2(v) == Catch::Approx(5.0));
        std::vector<Complex> w{Complex(1.0, 1.0)};   // old formula: sqrt((1+i)^2) = sqrt(2i) != sqrt(2)
        REQUIRE(SparseMatrix::norm2(w) == Catch::Approx(std::sqrt(2.0)));
        REQUIRE(A.NormFrobenius() == Catch::Approx(std::sqrt(4.0 + 2.0 + 2.0 + 9.0 + 1.0)));
    }

    SECTION("mat-vec, diagonal and lower solve with complex entries")
    {
        std::vector<Complex> x{Complex(1.0, 0.0), I, Complex(2.0, -1.0)};
        auto y = A * x;
        REQUIRE(std::abs(y[0] - (Complex(2.0) + Complex(1.0, 1.0) * I)) < 1e-15);
        REQUIRE(std::abs(y[1] - (Complex(1.0, 1.0) + Complex(0.0, 3.0) * I)) < 1e-15);
        REQUIRE(std::abs(y[2] - Complex(-2.0, 1.0)) < 1e-15);

        auto d = A.getDiagonal();
        REQUIRE(d[1] == Complex(0.0, 3.0));

        SparseMatrixCOO<Complex> lcoo(3, 3);
        lcoo.addEntry(0, 0, Complex(2.0, 0.0));
        lcoo.addEntry(1, 0, Complex(1.0, 1.0));
        lcoo.addEntry(1, 1, Complex(0.0, 3.0));
        lcoo.addEntry(2, 1, Complex(0.0, -1.0));
        lcoo.addEntry(2, 2, Complex(-1.0, 0.0));
        SparseMatrixCSR<Complex> L(lcoo);
        auto b = L * x;
        std::vector<Complex> z;
        L.SolveLower(b, z);
        for (int i = 0; i < 3; ++i) REQUIRE(std::abs(z[i] - x[i]) < 1e-14);
        std::vector<Complex> bu = A * x;   // A's strict lower part is ignored by solveUpper
        std::vector<Complex> zu;
        A.SolveUpper(bu, zu);
        REQUIRE(std::abs(zu[2] - x[2]) < 1e-14);
    }

    SECTION("COO Laplacian stencils accept complex scale factors")
    {
        SparseMatrixCOO<Complex> lap(3, 3);
        for (int i = 0; i < 3; ++i) lap.addLaplacian1D(i, i, Complex(1.0, 0.0));
        SparseMatrixCSR<Complex> L(lap);
        REQUIRE(L(1, 1) == Complex(-2.0, 0.0));
        SparseMatrixCOO<Complex> lap2(4, 4);
        lap2.addLaplacian2D_5pt(0, 0, 0, 2, 2, I);
        SparseMatrixCSR<Complex> L2(lap2);
        REQUIRE(L2(0, 0) == Complex(-4.0) * I);
    }
}

TEST_CASE("MatrixNDim<Complex>::IsEqualTo and arithmetic", "[complex-instantiation][tensor][wewi.2.1]")
{
    MatrixNDim<Complex, 2, 2> a, b;
    a(0, 0) = Complex(1.0, 1.0); a(1, 1) = Complex(0.0, -2.0);
    b = a;
    b(1, 1) += Complex(1e-13, 0.0);
    REQUIRE(a.IsEqualTo(b, static_cast<MML::Real>(1e-12)));
    REQUIRE_FALSE(a.IsEqualTo(b, static_cast<MML::Real>(1e-15)));
    auto s = a + b;
    REQUIRE(std::abs(s(0, 0) - Complex(2.0, 2.0)) < 1e-15);
}

TEST_CASE("Representative Complex free operators and norms across Vector/Matrix/MatrixSym", "[complex-instantiation][wewi.2.1]")
{
    Vector<Complex> v{Complex(3.0, 4.0), Complex(0.0, 0.0)};
    REQUIRE(v.NormL2() == Catch::Approx(5.0));
    REQUIRE(v.NormL1() == Catch::Approx(5.0));
    REQUIRE(v.NormLInf() == Catch::Approx(5.0));
    Vector<Complex> w = I * v;
    REQUIRE(std::abs(w[0] - Complex(-4.0, 3.0)) < 1e-15);
    Vector<Complex> u = v * I;
    REQUIRE(u.IsEqualTo(w));
    REQUIRE((v - v).IsZero());

    Matrix<Complex> m(2, 2, {Complex(1.0, 1.0), Complex(0.0, 1.0), Complex(2.0, 0.0), Complex(1.0, -1.0)});
    REQUIRE(m.Trace() == Complex(2.0, 0.0));
    Matrix<Complex> mi = MatrixAlg::Inverse(m);
    Matrix<Complex> id = m * mi;
    REQUIRE(id.IsEqualTo(Matrix<Complex>::Identity(2), static_cast<MML::Real>(1e-13)));
    Vector<Complex> mv = m * v;
    REQUIRE(std::abs(mv[0] - Complex(1.0, 1.0) * Complex(3.0, 4.0)) < 1e-14);
    std::ostringstream os;
    os << m;
    REQUIRE_FALSE(os.str().empty());

    MatrixSym<Complex> s(2);
    s(0, 0) = Complex(1.0, 0.0); s(0, 1) = Complex(0.0, 1.0); s(1, 1) = Complex(2.0, 0.0);
    REQUIRE(s(1, 0) == s(0, 1));
    REQUIRE(s.Trace() == Complex(3.0, 0.0));
    REQUIRE(s.NormFrobenius() == Catch::Approx(std::sqrt(1.0 + 1.0 + 1.0 + 4.0)));
    MatrixSym<Complex> s2 = s + s;
    REQUIRE(s2.IsEqualTo(s * Complex(2.0)));
}

} // namespace MML::Tests::Base::ComplexInstantiationTests
