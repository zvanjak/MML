#ifndef MML_USE_SINGLE_HEADER
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#endif
/**
 * @file eigensolver_adversarial_tests.cpp
 * @brief Adversarial tests for the general (nonsymmetric) real EigenSolver.
 *
 * Regression tests for MinimalMathLibrary-5372 (analysis v4, ALGORITHMS A1): the QR driver
 * abandoned the lower block on an interior deflation, took shifts from already-converged blocks,
 * and the Schur back-substitution ignored 2x2 blocks - producing wrong eigenvalues with
 * converged == true and eigenvector residuals of O(1). The pre-existing tests checked only the
 * trace of a 4x4 and eigenvectors of a symmetric 3x3, so none of this was red.
 *
 * Every test here checks the thing that actually matters: ||A v - lambda v|| / ||v|| for every
 * returned eigenpair (complex pairs included), and the full spectrum against known values.
 */

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/algorithms/Eigen/detail/RealSchurAnalysis.h>
#include <mml/algorithms/Eigen/detail/SchurEigenvectors.h>
#include <mml/algorithms/Eigen/EigenSolver.h>
#endif

#include <catch2/catch_all.hpp>
#include <algorithm>
#include <complex>
#include <random>
#include <vector>
#include "../TestPrecision.h"

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Algorithms::EigensolverAdversarialTests
{

namespace
{
    using CReal = std::complex<Real>;

    // Largest ||A v - lambda v|| / ||v|| over all eigenpairs (complex pairs use both columns).
    Real MaxEigenpairResidual(const Matrix<Real>& A, const EigenSolver::Result& r)
    {
        const int n = A.rows();
        Real worst = static_cast<MML::Real>(0.0);
        Vector<Real> vr(n), vi(n);
        std::size_t i = 0;
        while (i < r.eigenvalues.size())
        {
            const bool pair = r.isComplexPair[i] && i + 1 < r.eigenvalues.size();
            for (int j = 0; j < n; j++)
            {
                vr[j] = r.eigenvectors(j, static_cast<int>(i));
                vi[j] = pair ? r.eigenvectors(j, static_cast<int>(i + 1)) : static_cast<MML::Real>(0.0);
            }
            const Real res = MML::detail::EigenpairResidual(A, vr, vi, r.eigenvalues[i].real,
                                                            pair ? r.eigenvalues[i].imag : static_cast<MML::Real>(0.0));
            worst = std::max(worst, res);
            i += pair ? 2 : 1;
        }
        return worst;
    }

    // Every expected eigenvalue must be matched by a computed one within tol.
    void RequireSpectrum(const EigenSolver::Result& r, std::vector<CReal> expected, Real tol)
    {
        REQUIRE(r.eigenvalues.size() == expected.size());
        std::vector<CReal> got;
        for (const auto& e : r.eigenvalues)
            got.emplace_back(e.real, e.imag);
        for (const auto& ev : expected)
        {
            Real best = std::numeric_limits<Real>::max();
            for (const auto& g : got)
                best = std::min(best, std::abs(g - ev));
            INFO("expected eigenvalue " << ev.real() << " + " << ev.imag() << "i, closest distance " << best);
            REQUIRE(best < tol);
        }
    }

    Real Trace(const Matrix<Real>& A)
    {
        Real t = static_cast<MML::Real>(0.0);
        for (int i = 0; i < A.rows(); i++)
            t += A(i, i);
        return t;
    }

    Real SumOfEigenvalues(const EigenSolver::Result& r)
    {
        Real s = static_cast<MML::Real>(0.0);
        for (const auto& e : r.eigenvalues)
            s += e.real;
        return s;
    }

    // S diag(1,2,3) S^-1 with S = [[1,1,0],[0,1,1],[1,0,1]]: dense 3x3 with spectrum {1,2,3}
    Matrix<Real> DenseBlock123()
    {
        Matrix<Real> S{3, 3, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0),
                              static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0),
                              static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)}};
        Matrix<Real> Sinv{3, 3, {static_cast<MML::Real>(0.5), -static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.5),
                                 static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.5), -static_cast<MML::Real>(0.5),
                                 -static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.5)}};
        Matrix<Real> D{3, 3, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0),
                              static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0),
                              static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0)}};
        return S * D * Sinv;
    }
} // namespace

TEST_CASE("EigenSolver - block upper triangular, dense block at the bottom (interior deflation)",
          "[eigensolver][adversarial][deflation]")
{
    // Top-left 3x3 is upper triangular with 10, 20, 30; bottom-right 3x3 is dense with 1, 2, 3.
    // The Hessenberg form splits at an INTERIOR subdiagonal; the old driver then iterated only
    // the upper block and reported 2.881, 1.319, 1.800 for the lower one with converged == true.
    Matrix<Real> A(6, 6);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            A(i, j) = static_cast<MML::Real>(0.0);
    A(0, 0) = static_cast<MML::Real>(10.0); A(1, 1) = static_cast<MML::Real>(20.0); A(2, 2) = static_cast<MML::Real>(30.0);
    A(0, 1) = static_cast<MML::Real>(1.0);  A(0, 2) = static_cast<MML::Real>(2.0);  A(1, 2) = static_cast<MML::Real>(3.0);
    for (int i = 0; i < 3; i++)
        for (int j = 3; j < 6; j++)
            A(i, j) = static_cast<MML::Real>(1.0);
    const Matrix<Real> B = DenseBlock123();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A(3 + i, 3 + j) = B(i, j);

    auto r = EigenSolver::Solve(A);

    INFO("iterations " << r.iterations << " maxResidual " << r.maxResidual);
    REQUIRE(r.converged);
    REQUIRE(r.status == AlgorithmStatus::Success);
    RequireSpectrum(r, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(20.0), static_cast<MML::Real>(30.0)}, TOL(1e-9, 1e-3));
    REQUIRE(MaxEigenpairResidual(A, r) < TOL(1e-9, 1e-3));
}

TEST_CASE("EigenSolver - block upper triangular, dense block at the top (shift from converged block)",
          "[eigensolver][adversarial][deflation]")
{
    // Mirror image: the old driver deflated 30, 20, 10 from the bottom, then kept taking Wilkinson
    // shifts from the bottom 2x2 of the FULL matrix (already converged) and failed after 1000 iterations.
    Matrix<Real> A(6, 6);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            A(i, j) = static_cast<MML::Real>(0.0);
    const Matrix<Real> B = DenseBlock123();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A(i, j) = B(i, j);
    A(3, 3) = static_cast<MML::Real>(10.0); A(4, 4) = static_cast<MML::Real>(20.0); A(5, 5) = static_cast<MML::Real>(30.0);
    A(3, 4) = static_cast<MML::Real>(1.0);  A(3, 5) = static_cast<MML::Real>(2.0);  A(4, 5) = static_cast<MML::Real>(3.0);
    for (int i = 0; i < 3; i++)
        for (int j = 3; j < 6; j++)
            A(i, j) = static_cast<MML::Real>(1.0);

    auto r = EigenSolver::Solve(A);

    INFO("iterations " << r.iterations << " maxResidual " << r.maxResidual);
    REQUIRE(r.converged);
    REQUIRE(r.iterations < 100);
    RequireSpectrum(r, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(10.0), static_cast<MML::Real>(20.0), static_cast<MML::Real>(30.0)}, TOL(1e-9, 1e-3));
    REQUIRE(MaxEigenpairResidual(A, r) < TOL(1e-9, 1e-3));
}

TEST_CASE("EigenSolver - dense nonsymmetric 6x6: eigenvector residuals, not just the trace",
          "[eigensolver][adversarial][eigenvector]")
{
    // The old back-substitution returned eigenvectors with ||Av - lv||/||v|| = 0.51 here.
    Matrix<Real> A{6, 6, {static_cast<MML::Real>(4.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0),
                          static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0),
                          static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(6.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(7.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0),
                          static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(8.0), static_cast<MML::Real>(2.0),
                          static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(9.0)}};

    auto r = EigenSolver::Solve(A);

    REQUIRE(r.converged);
    REQUIRE(r.eigenvalues.size() == 6);
    REQUIRE(std::abs(SumOfEigenvalues(r) - Trace(A)) < TOL(1e-9, 1e-3));
    const Real residual = MaxEigenpairResidual(A, r);
    INFO("max eigenpair residual " << residual);
    REQUIRE(residual < TOL(1e-9, 1e-3));
    REQUIRE(r.maxResidual < TOL(1e-9, 1e-3)); // the solver's own report must agree
}

TEST_CASE("EigenSolver - complex conjugate pairs: eigenvalues and complex eigenvectors",
          "[eigensolver][adversarial][complex]")
{
    // Block-diagonal with rotations (+-i, 2+-3i) and 7, then a dense similarity transform of it.
    Matrix<Real> A{5, 5, {static_cast<MML::Real>(0.0), -static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0),
                          static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), -static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(7.0)}};
    const std::vector<CReal> expected{{static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)}, {static_cast<MML::Real>(0.0), -static_cast<MML::Real>(1.0)},
                                      {static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)}, {static_cast<MML::Real>(2.0), -static_cast<MML::Real>(3.0)},
                                      {static_cast<MML::Real>(7.0), static_cast<MML::Real>(0.0)}};

    SECTION("block diagonal")
    {
        auto r = EigenSolver::Solve(A);
        REQUIRE(r.converged);
        RequireSpectrum(r, expected, TOL(1e-9, 1e-3));
        REQUIRE(MaxEigenpairResidual(A, r) < TOL(1e-9, 1e-3));
    }

    SECTION("dense similarity transform")
    {
        std::mt19937 gen(7);
        std::uniform_real_distribution<Real> u(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0));
        Matrix<Real> S(5, 5);
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                S(i, j) = u(gen) + (i == j ? static_cast<MML::Real>(3.0) : static_cast<MML::Real>(0.0));
        const Matrix<Real> M = S * A * MatrixAlg::Inverse(S);

        auto r = EigenSolver::Solve(M);
        REQUIRE(r.converged);
        RequireSpectrum(r, expected, TOL(1e-8, 1e-2));
        REQUIRE(MaxEigenpairResidual(M, r) < TOL(1e-8, 1e-2));

        // Pair convention relied upon by MatrixAlg/Eigensystems.h: +beta first, columns k / k+1 = Re / Im
        for (std::size_t i = 0; i + 1 < r.eigenvalues.size(); i++)
            if (r.isComplexPair[i] && (i == 0 || !r.isComplexPair[i - 1]))
            {
                REQUIRE(r.eigenvalues[i].imag > static_cast<MML::Real>(0.0));
                REQUIRE(r.eigenvalues[i + 1].imag == Catch::Approx(-r.eigenvalues[i].imag));
                REQUIRE(r.eigenvalues[i + 1].real == Catch::Approx(r.eigenvalues[i].real));
            }
    }
}

TEST_CASE("EigenSolver - companion matrix (polynomial roots path)", "[eigensolver][adversarial][companion]")
{
    // (x-1)(x-2)(x+3)(x^2+1) = x^5 + 0x^4 - 6x^3 + 0x^2 - 7x + 6  ->  check by expansion below
    std::vector<Real> p{static_cast<MML::Real>(1.0)};
    auto mul = [&](std::vector<Real> q) {
        std::vector<Real> out(p.size() + q.size() - 1, static_cast<MML::Real>(0.0));
        for (std::size_t i = 0; i < p.size(); i++)
            for (std::size_t j = 0; j < q.size(); j++)
                out[i + j] += p[i] * q[j];
        p = out;
    };
    mul({static_cast<MML::Real>(1.0), -static_cast<MML::Real>(1.0)});
    mul({static_cast<MML::Real>(1.0), -static_cast<MML::Real>(2.0)});
    mul({static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0)});
    mul({static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});

    const int n = static_cast<int>(p.size()) - 1;
    Matrix<Real> C(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C(i, j) = static_cast<MML::Real>(0.0);
    for (int j = 0; j < n; j++)
        C(0, j) = -p[j + 1] / p[0];
    for (int i = 1; i < n; i++)
        C(i, i - 1) = static_cast<MML::Real>(1.0);

    auto r = EigenSolver::Solve(C);
    REQUIRE(r.converged);
    RequireSpectrum(r, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), -static_cast<MML::Real>(3.0), {static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)}, {static_cast<MML::Real>(0.0), -static_cast<MML::Real>(1.0)}},
                    TOL(1e-8, 1e-2));
    REQUIRE(MaxEigenpairResidual(C, r) < TOL(1e-8, 1e-2));
}

TEST_CASE("EigenSolver - random dense matrices: residual, trace, honest convergence report",
          "[eigensolver][adversarial][random]")
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<Real> u(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0));

    for (int n : {8, 12, 20})
    {
        Matrix<Real> A(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A(i, j) = u(gen);

        auto r = EigenSolver::Solve(A);
        INFO("n = " << n << " iterations " << r.iterations << " maxResidual " << r.maxResidual);
        REQUIRE(r.converged);
        REQUIRE(r.status == AlgorithmStatus::Success);
        REQUIRE(r.eigenvalues.size() == static_cast<std::size_t>(n));
        REQUIRE(std::abs(SumOfEigenvalues(r) - Trace(A)) < TOL(1e-9, 1e-3));
        REQUIRE(MaxEigenpairResidual(A, r) < TOL(1e-8, 1e-2));
    }

    SECTION("maxIter exhausted is reported, never converged == true with iterations == 0 on an unreduced matrix")
    {
        Matrix<Real> A(8, 8);
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                A(i, j) = u(gen);
        auto r = EigenSolver::Solve(A, PrecisionValues<Real>::EigenSolverConvergenceTolerance, 1);
        REQUIRE_FALSE(r.converged);
        REQUIRE(r.status == AlgorithmStatus::MaxIterationsExceeded);
        REQUIRE_FALSE(r.error_message.empty());
    }
}

TEST_CASE("EigenSolver - repeated and defective eigenvalues do not break the solver",
          "[eigensolver][adversarial][repeated]")
{
    SECTION("symmetric with repeated eigenvalue 3, 3, 6")
    {
        Matrix<Real> A{3, 3, {static_cast<MML::Real>(4.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0),
                              static_cast<MML::Real>(1.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(1.0),
                              static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(4.0)}};
        auto r = EigenSolver::Solve(A);
        REQUIRE(r.converged);
        RequireSpectrum(r, {static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(6.0)}, TOL(1e-9, 1e-3));
        REQUIRE(MaxEigenpairResidual(A, r) < TOL(1e-9, 1e-3));
    }

    SECTION("Jordan block (defective): eigenvalues exact, finite eigenvectors, no division by zero")
    {
        Matrix<Real> J{4, 4, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0),
                              static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0),
                              static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0),
                              static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0)}};
        auto r = EigenSolver::Solve(J);
        REQUIRE(r.converged);
        RequireSpectrum(r, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0)}, TOL(1e-9, 1e-3));
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                REQUIRE(std::isfinite(r.eigenvectors(i, j)));
    }
}

TEST_CASE("EigenSolver - input validation", "[eigensolver][adversarial][validation]")
{
    Matrix<Real> rect(2, 3);
    REQUIRE_THROWS_AS(EigenSolver::Solve(rect), MatrixDimensionError);

    Matrix<Real> A{2, 2, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)}};
    REQUIRE_THROWS_AS(EigenSolver::Solve(A, -static_cast<MML::Real>(1.0)), DomainError);
}

// ---------------------------------------------------------------------------------------------
// Building blocks added with the fix
// ---------------------------------------------------------------------------------------------

TEST_CASE("Schur building blocks - FindUnreducedBlockStart returns the LOWEST split, not the first from the bottom",
          "[eigensolver][building-block][deflation]")
{
    // Splits at (2,1) and (4,3): the active block ending at iHi = 5 starts at 4.
    Matrix<Real> H(6, 6);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            H(i, j) = (j >= i - 1) ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(0.0);
    H(2, 1) = static_cast<MML::Real>(0.0);
    H(4, 3) = static_cast<MML::Real>(1e-14);

    const Real anorm = MML::detail::HessenbergAbsNorm(H);
    REQUIRE(MML::detail::FindUnreducedBlockStart(H, 0, 5, static_cast<MML::Real>(1e-10), anorm) == 4);
    REQUIRE(MML::detail::FindUnreducedBlockStart(H, 0, 3, static_cast<MML::Real>(1e-10), anorm) == 2);
    REQUIRE(MML::detail::FindUnreducedBlockStart(H, 0, 1, static_cast<MML::Real>(1e-10), anorm) == 0);
}

TEST_CASE("Schur building blocks - StandardizeSchurBlock triangularises a real 2x2 block and keeps the similarity",
          "[eigensolver][building-block][schur]")
{
    // Embed [[4,1],[2,3]] (eigenvalues 5, 2) as the trailing block of a 4x4 quasi-triangular matrix.
    Matrix<Real> H{4, 4, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(6.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(1.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)}};
    const Matrix<Real> H0 = H;
    Matrix<Real> Q = Matrix<Real>::Identity(4);

    REQUIRE(MML::detail::StandardizeSchurBlock(H, Q, 2));
    REQUIRE(H(3, 2) == static_cast<MML::Real>(0.0));
    std::vector<Real> diag{H(2, 2), H(3, 3)};
    std::sort(diag.begin(), diag.end());
    REQUIRE(diag[0] == Catch::Approx(static_cast<MML::Real>(2.0)).margin(TOL(1e-12, 1e-5)));
    REQUIRE(diag[1] == Catch::Approx(static_cast<MML::Real>(5.0)).margin(TOL(1e-12, 1e-5)));

    // H == Q^T H0 Q  and Q orthogonal
    const Matrix<Real> QtH0Q = Q.Transposed() * H0 * Q;
    const Matrix<Real> QtQ = Q.Transposed() * Q;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            REQUIRE(QtH0Q(i, j) == Catch::Approx(H(i, j)).margin(TOL(1e-12, 1e-5)));
            REQUIRE(QtQ(i, j) == Catch::Approx(i == j ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(0.0)).margin(TOL(1e-12, 1e-5)));
        }

    // A genuine complex block is left alone
    Matrix<Real> R{2, 2, {static_cast<MML::Real>(0.0), -static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}};
    Matrix<Real> Q2 = Matrix<Real>::Identity(2);
    REQUIRE_FALSE(MML::detail::StandardizeSchurBlock(R, Q2, 0));
    REQUIRE(R(1, 0) == static_cast<MML::Real>(1.0));
}

TEST_CASE("Schur building blocks - back-substitution honours a complex 2x2 block ABOVE the target eigenvalue",
          "[eigensolver][building-block][eigenvector]")
{
    // T = [[0,-1, 1],[1, 0, 2],[0, 0, 3]]: complex block (+-i) on top, real eigenvalue 3 below.
    Matrix<Real> T{3, 3, {static_cast<MML::Real>(0.0), -static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0),
                          static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0),
                          static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0)}};

    Vector<Real> v = MML::detail::ComputeRealEigenvector(T, 2);
    Vector<Real> zero(3, static_cast<MML::Real>(0.0));
    REQUIRE(MML::detail::EigenpairResidual(T, v, zero, static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0)) < TOL(1e-12, 1e-5));

    auto [vr, vi] = MML::detail::ComputeComplexEigenvectors(T, 0);
    REQUIRE(MML::detail::EigenpairResidual(T, vr, vi, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)) < TOL(1e-12, 1e-5));

    // Full pipeline on this T with Q = I
    auto vecs = MML::detail::ComputeEigenvectorsFromSchur(T, Matrix<Real>::Identity(3), static_cast<MML::Real>(0.0));
    REQUIRE(vecs.isComplexPair[0]);
    REQUIRE(vecs.isComplexPair[1]);
    REQUIRE_FALSE(vecs.isComplexPair[2]);
}

} // namespace MML::Tests::Algorithms::EigensolverAdversarialTests
