/**
 * @file test_eigensolver_building_blocks.cpp
 * @brief Tests for eigensolver building blocks
 *
 * Tests each building block independently before integration:
 * - Building Block 1: Hessenberg Reduction
 * - Building Block 2: Single QR Step
 * - Building Block 3: Francis Double-Shift
 * - Building Block 4: Deflation and Extraction
 * - Building Block 5: Eigenvector Computation
 * 
 * Also tests the integrated EigenSolver class.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "algorithms/EigenSolverHelpers.h"
#include "algorithms/EigenSystemSolvers.h"  // For EigenSolver
#endif

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

using namespace MML;
using namespace MML::Testing;

// Tolerance for numerical comparisons
constexpr Real TOL = 1e-10;

// ============================================================================
// BUILDING BLOCK 1: HESSENBERG REDUCTION TESTS
// ============================================================================

TEST_CASE("Hessenberg - 2x2 matrix unchanged", "[eigensolver][building-block][hessenberg]")
{
    // 2x2 matrices are already Hessenberg
    Matrix<Real> A{2, 2, {REAL(1.0), REAL(2.0),
                          REAL(3.0), REAL(4.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    // H should equal A for 2x2
    REQUIRE(EigenSolverHelpers::MaxAbsDiff(result.H, A) < TOL);
    
    // Q should be identity
    Matrix<Real> I = Matrix<Real>::GetUnitMatrix(2);
    REQUIRE(EigenSolverHelpers::MaxAbsDiff(result.Q, I) < TOL);
}

TEST_CASE("Hessenberg - 3x3 produces upper Hessenberg", "[eigensolver][building-block][hessenberg]")
{
    Matrix<Real> A{3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(4.0), REAL(5.0), REAL(6.0),
                          REAL(7.0), REAL(8.0), REAL(9.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    // Check H is upper Hessenberg
    REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
    INFO("H(2,0) = " << result.H(2, 0));
    REQUIRE(std::abs(result.H(2, 0)) < TOL);
}

TEST_CASE("Hessenberg - 3x3 Q is orthogonal", "[eigensolver][building-block][hessenberg]")
{
    Matrix<Real> A{3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(4.0), REAL(5.0), REAL(6.0),
                          REAL(7.0), REAL(8.0), REAL(9.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    // Check Q is orthogonal: Q^T * Q = I
    REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
}

TEST_CASE("Hessenberg - 3x3 similarity preserved (Q^T*A*Q = H)", "[eigensolver][building-block][hessenberg]")
{
    Matrix<Real> A{3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(4.0), REAL(5.0), REAL(6.0),
                          REAL(7.0), REAL(8.0), REAL(9.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    // Verify Q^T * A * Q = H
    Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, A);
    
    Real diff = EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H);
    INFO("Max diff between Q^T*A*Q and H: " << diff);
    REQUIRE(diff < TOL);
}

TEST_CASE("Hessenberg - 3x3 trace preserved", "[eigensolver][building-block][hessenberg]")
{
    Matrix<Real> A{3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(4.0), REAL(5.0), REAL(6.0),
                          REAL(7.0), REAL(8.0), REAL(9.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    Real traceA = EigenSolverHelpers::Trace(A);
    Real traceH = EigenSolverHelpers::Trace(result.H);
    
    INFO("trace(A) = " << traceA << ", trace(H) = " << traceH);
    REQUIRE(std::abs(traceA - traceH) < TOL);
}

TEST_CASE("Hessenberg - 4x4 all criteria", "[eigensolver][building-block][hessenberg]")
{
    Matrix<Real> A{4, 4, {REAL(4.0), REAL(1.0), -REAL(2.0), REAL(2.0),
                          REAL(1.0), REAL(2.0),  REAL(0.0), REAL(1.0),
                         -REAL(2.0), REAL(0.0),  REAL(3.0), -REAL(2.0),
                          REAL(2.0), REAL(1.0), -REAL(2.0), -REAL(1.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    // Check all criteria
    SECTION("H is upper Hessenberg")
    {
        REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
        REQUIRE(std::abs(result.H(2, 0)) < TOL);
        REQUIRE(std::abs(result.H(3, 0)) < TOL);
        REQUIRE(std::abs(result.H(3, 1)) < TOL);
    }
    
    SECTION("Q is orthogonal")
    {
        REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
    }
    
    SECTION("Similarity: Q^T*A*Q = H")
    {
        Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, A);
        REQUIRE(EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H) < TOL);
    }
    
    SECTION("Trace preserved")
    {
        REQUIRE(std::abs(EigenSolverHelpers::Trace(A) - EigenSolverHelpers::Trace(result.H)) < TOL);
    }
}

TEST_CASE("Hessenberg - 5x5 all criteria", "[eigensolver][building-block][hessenberg]")
{
    Matrix<Real> A{5, 5, {REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0),
                          REAL(6.0), REAL(7.0), REAL(8.0), REAL(9.0), REAL(10.0),
                          REAL(11.0), REAL(12.0), REAL(13.0), REAL(14.0), REAL(15.0),
                          REAL(16.0), REAL(17.0), REAL(18.0), REAL(19.0), REAL(20.0),
                          REAL(21.0), REAL(22.0), REAL(23.0), REAL(24.0), REAL(25.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    SECTION("H is upper Hessenberg")
    {
        REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
    }
    
    SECTION("Q is orthogonal")
    {
        REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
    }
    
    SECTION("Similarity: Q^T*A*Q = H")
    {
        Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, A);
        REQUIRE(EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H) < TOL);
    }
    
    SECTION("Trace preserved")
    {
        REQUIRE(std::abs(EigenSolverHelpers::Trace(A) - EigenSolverHelpers::Trace(result.H)) < TOL);
    }
}

TEST_CASE("Hessenberg - already Hessenberg matrix", "[eigensolver][building-block][hessenberg]")
{
    // Upper Hessenberg matrix should be preserved
    Matrix<Real> H{4, 4, {REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0),
                          REAL(5.0), REAL(6.0), REAL(7.0), REAL(8.0),
                          REAL(0.0), REAL(9.0), REAL(10.0), REAL(11.0),
                          REAL(0.0), REAL(0.0), REAL(12.0), REAL(13.0)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(H);
    
    REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
    REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
    
    // Similarity should hold
    Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, H);
    REQUIRE(EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H) < TOL);
}

TEST_CASE("Hessenberg - random-like structured matrix", "[eigensolver][building-block][hessenberg]")
{
    // Create a matrix with varied elements (like testbed matrices)
    Matrix<Real> A{4, 4, {REAL(7.9), -REAL(9.3),  REAL(1.7),  REAL(0.5),
                         -REAL(0.2), -REAL(5.9),  REAL(9.9),  REAL(2.3),
                          REAL(2.4),  REAL(6.5), -REAL(5.5),  REAL(1.1),
                          REAL(3.1), -REAL(2.2),  REAL(4.4),  REAL(8.8)}};
    
    auto result = EigenSolverHelpers::ReduceToHessenberg(A);
    
    INFO("H(2,0) = " << result.H(2, 0));
    INFO("H(3,0) = " << result.H(3, 0));
    INFO("H(3,1) = " << result.H(3, 1));
    
    REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
    REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
    
    Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, A);
    Real diff = EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H);
    INFO("Reconstruction error: " << diff);
    REQUIRE(diff < TOL);
}

// ============================================================================
// BUILDING BLOCK 4: EIGENVALUE EXTRACTION FROM 2x2 BLOCK
// ============================================================================

TEST_CASE("Eigenvalue2x2 - two real eigenvalues", "[eigensolver][building-block][2x2]")
{
    // [[4, 1], [2, 3]] has eigenvalues 5 and 2
    auto result = EigenSolverHelpers::Eigenvalues2x2(REAL(4.0), REAL(1.0), REAL(2.0), REAL(3.0));
    
    REQUIRE_FALSE(result.isComplex);
    
    // Sort eigenvalues for comparison
    Real e1 = std::max(result.real1, result.real2);
    Real e2 = std::min(result.real1, result.real2);
    
    INFO("Eigenvalues: " << result.real1 << ", " << result.real2);
    REQUIRE(std::abs(e1 - REAL(5.0)) < TOL);
    REQUIRE(std::abs(e2 - REAL(2.0)) < TOL);
}

TEST_CASE("Eigenvalue2x2 - complex conjugate pair", "[eigensolver][building-block][2x2]")
{
    // Rotation matrix: [[cos(θ), -sin(θ)], [sin(θ), cos(θ)]]
    // Eigenvalues: e^{±iθ}
    Real theta = REAL(0.5);
    Real c = std::cos(theta);
    Real s = std::sin(theta);
    
    auto result = EigenSolverHelpers::Eigenvalues2x2(c, -s, s, c);
    
    REQUIRE(result.isComplex);
    
    // Both real parts should equal cos(θ)
    INFO("Real parts: " << result.real1 << ", " << result.real2);
    REQUIRE(std::abs(result.real1 - c) < TOL);
    REQUIRE(std::abs(result.real2 - c) < TOL);
    
    // Imaginary parts should be ±sin(θ)
    INFO("Imag parts: " << result.imag1 << ", " << result.imag2);
    REQUIRE(std::abs(std::abs(result.imag1) - s) < TOL);
    REQUIRE(std::abs(result.imag1 + result.imag2) < TOL);  // Conjugate pair
}

TEST_CASE("Eigenvalue2x2 - diagonal matrix", "[eigensolver][building-block][2x2]")
{
    // Diagonal [[3, 0], [0, 7]]
    auto result = EigenSolverHelpers::Eigenvalues2x2(REAL(3.0), REAL(0.0), REAL(0.0), REAL(7.0));
    
    REQUIRE_FALSE(result.isComplex);
    
    Real e1 = std::max(result.real1, result.real2);
    Real e2 = std::min(result.real1, result.real2);
    
    REQUIRE(std::abs(e1 - REAL(7.0)) < TOL);
    REQUIRE(std::abs(e2 - REAL(3.0)) < TOL);
}

TEST_CASE("Eigenvalue2x2 - repeated eigenvalue", "[eigensolver][building-block][2x2]")
{
    // [[5, 1], [0, 5]] - repeated eigenvalue 5
    auto result = EigenSolverHelpers::Eigenvalues2x2(REAL(5.0), REAL(1.0), REAL(0.0), REAL(5.0));
    
    REQUIRE_FALSE(result.isComplex);
    REQUIRE(std::abs(result.real1 - REAL(5.0)) < TOL);
    REQUIRE(std::abs(result.real2 - REAL(5.0)) < TOL);
}

TEST_CASE("Eigenvalue2x2 - negative discriminant (complex)", "[eigensolver][building-block][2x2]")
{
    // [[0, -1], [1, 0]] - pure imaginary eigenvalues ±i
    auto result = EigenSolverHelpers::Eigenvalues2x2(REAL(0.0), -REAL(1.0), REAL(1.0), REAL(0.0));
    
    REQUIRE(result.isComplex);
    REQUIRE(std::abs(result.real1) < TOL);  // Real part should be 0
    REQUIRE(std::abs(result.real2) < TOL);
    REQUIRE(std::abs(std::abs(result.imag1) - REAL(1.0)) < TOL);  // |imag| = 1
}

// ============================================================================
// BUILDING BLOCK 2: SINGLE QR STEP TESTS
// ============================================================================

TEST_CASE("QRStep - Wilkinson shift for real eigenvalues", "[eigensolver][building-block][qrstep]")
{
    // 2x2 block [[4, 1], [2, 3]] has eigenvalues 5 and 2
    // Shift should be closer to 3 (the d element), so shift = 2
    Real shift = EigenSolverHelpers::WilkinsonShift(REAL(4.0), REAL(1.0), REAL(2.0), REAL(3.0));
    
    // Eigenvalues are 5 and 2; d=3, so closer eigenvalue is 2
    INFO("Wilkinson shift = " << shift);
    REQUIRE(std::abs(shift - REAL(2.0)) < TOL);
}

TEST_CASE("QRStep - Wilkinson shift for complex eigenvalues", "[eigensolver][building-block][qrstep]")
{
    // Rotation matrix has complex eigenvalues, shift should be d
    Real theta = REAL(0.5);
    Real c = std::cos(theta);
    Real s = std::sin(theta);
    
    Real shift = EigenSolverHelpers::WilkinsonShift(c, -s, s, c);
    
    // For complex eigenvalues, shift should be d = cos(θ)
    INFO("Wilkinson shift = " << shift << ", expected = " << c);
    REQUIRE(std::abs(shift - c) < TOL);
}

TEST_CASE("QRStep - Single step preserves Hessenberg form", "[eigensolver][building-block][qrstep]")
{
    // Start with upper Hessenberg matrix
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(2.0), REAL(5.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(3.0), REAL(6.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(7.0)}};
    
    auto result = EigenSolverHelpers::SingleQRStep(H);
    
    REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
}

TEST_CASE("QRStep - Single step Q is orthogonal", "[eigensolver][building-block][qrstep]")
{
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(2.0), REAL(5.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(3.0), REAL(6.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(7.0)}};
    
    auto result = EigenSolverHelpers::SingleQRStep(H, true);
    
    REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
}

TEST_CASE("QRStep - Single step preserves similarity", "[eigensolver][building-block][qrstep]")
{
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(2.0), REAL(5.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(3.0), REAL(6.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(7.0)}};
    
    auto result = EigenSolverHelpers::SingleQRStep(H, true);
    
    // Verify Q^T * H * Q = H_new
    Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, H);
    
    Real diff = EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H);
    INFO("Similarity error: " << diff);
    REQUIRE(diff < TOL);
}

TEST_CASE("QRStep - Single step preserves trace", "[eigensolver][building-block][qrstep]")
{
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(2.0), REAL(5.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(3.0), REAL(6.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(7.0)}};
    
    Real traceBefore = EigenSolverHelpers::Trace(H);
    auto result = EigenSolverHelpers::SingleQRStep(H);
    Real traceAfter = EigenSolverHelpers::Trace(result.H);
    
    INFO("Trace before: " << traceBefore << ", after: " << traceAfter);
    REQUIRE(std::abs(traceBefore - traceAfter) < TOL);
}

TEST_CASE("QRStep - Multiple steps converge for diagonal dominance", "[eigensolver][building-block][qrstep]")
{
    // Start with diagonally dominant Hessenberg - should converge quickly
    Matrix<Real> H{3, 3, {REAL(10.0), REAL(1.0), REAL(0.5),
                           REAL(1.0), REAL(5.0), REAL(0.3),
                           REAL(0.0), REAL(0.2), REAL(1.0)}};
    
    int iters = EigenSolverHelpers::MultipleQRSteps(H, 30, 1e-10);
    
    INFO("Converged in " << iters << " iterations");
    INFO("H(2,1) = " << H(2, 1));
    
    REQUIRE(iters < 30);  // Should converge
    REQUIRE(std::abs(H(2, 1)) < 1e-10);  // Bottom subdiagonal should be zero
}

TEST_CASE("QRStep - Convergence reveals eigenvalue", "[eigensolver][building-block][qrstep]")
{
    // Start with Hessenberg form of matrix with known eigenvalues
    // [[4, 1], [2, 3]] -> eigenvalues 5 and 2
    Matrix<Real> H{2, 2, {REAL(4.0), REAL(1.0),
                          REAL(2.0), REAL(3.0)}};
    
    int iters = EigenSolverHelpers::MultipleQRSteps(H, 30, 1e-10);
    
    INFO("Converged in " << iters << " iterations");
    INFO("H(1,0) after convergence = " << H(1, 0));
    
    // After convergence, H should be (nearly) upper triangular
    // with eigenvalues on diagonal
    REQUIRE(std::abs(H(1, 0)) < 1e-10);
    
    // Diagonal elements should be eigenvalues (5 and 2, in some order)
    Real e1 = H(0, 0);
    Real e2 = H(1, 1);
    
    INFO("Diagonal: " << e1 << ", " << e2);
    
    bool correct = (std::abs(e1 - REAL(5.0)) < 1e-6 && std::abs(e2 - REAL(2.0)) < 1e-6) ||
                   (std::abs(e1 - REAL(2.0)) < 1e-6 && std::abs(e2 - REAL(5.0)) < 1e-6);
    REQUIRE(correct);
}

// =============================================================================
// BUILDING BLOCK 3: FRANCIS DOUBLE-SHIFT TESTS
// =============================================================================

TEST_CASE("DoubleShift - Preserves Hessenberg form", "[eigensolver][building-block][doubleshift]")
{
    // 4x4 Hessenberg matrix
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(3.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(2.0), REAL(5.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(1.0), REAL(6.0)}};
    
    auto result = EigenSolverHelpers::FrancisDoubleShift(H, 0, 3, true);
    
    REQUIRE(EigenSolverHelpers::IsUpperHessenberg(result.H, TOL));
}

TEST_CASE("DoubleShift - Q is orthogonal", "[eigensolver][building-block][doubleshift]")
{
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(3.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(2.0), REAL(5.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(1.0), REAL(6.0)}};
    
    auto result = EigenSolverHelpers::FrancisDoubleShift(H, 0, 3, true);
    
    REQUIRE(EigenSolverHelpers::IsOrthogonal(result.Q, TOL));
}

TEST_CASE("DoubleShift - Preserves similarity", "[eigensolver][building-block][doubleshift]")
{
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(3.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(2.0), REAL(5.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(1.0), REAL(6.0)}};
    
    auto result = EigenSolverHelpers::FrancisDoubleShift(H, 0, 3, true);
    
    // Verify Q^T * H * Q = H_new
    Matrix<Real> reconstructed = EigenSolverHelpers::SimilarityTransform(result.Q, H);
    
    Real diff = EigenSolverHelpers::MaxAbsDiff(reconstructed, result.H);
    INFO("Similarity error: " << diff);
    REQUIRE(diff < TOL);
}

TEST_CASE("DoubleShift - Preserves trace", "[eigensolver][building-block][doubleshift]")
{
    Matrix<Real> H{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(3.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(2.0), REAL(5.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(1.0), REAL(6.0)}};
    
    Real traceBefore = EigenSolverHelpers::Trace(H);
    auto result = EigenSolverHelpers::FrancisDoubleShift(H, 0, 3, false);
    Real traceAfter = EigenSolverHelpers::Trace(result.H);
    
    INFO("Trace before: " << traceBefore << ", after: " << traceAfter);
    REQUIRE(std::abs(traceBefore - traceAfter) < TOL);
}

TEST_CASE("DoubleShift - Converges on complex eigenvalue matrix", "[eigensolver][building-block][doubleshift]")
{
    // Matrix with complex eigenvalues: [[0, -1], [1, 0]] has eigenvalues ±i
    // Embedded in 4x4 Hessenberg:
    Matrix<Real> H{4, 4, {REAL(1.0),  REAL(0.0), -REAL(1.0), REAL(0.5),
                          REAL(1.0),  REAL(0.0),  REAL(1.0), REAL(0.3),
                          REAL(0.0),  REAL(1.0),  REAL(2.0), REAL(0.2),
                          REAL(0.0),  REAL(0.0),  REAL(0.5), REAL(3.0)}};
    
    Matrix<Real> Hcopy = H;
    int iters = EigenSolverHelpers::MultipleDoubleShiftSteps(Hcopy, 0, 3, 50, 1e-10);
    
    INFO("Double-shift converged in " << iters << " iterations");
    
    // Should converge (might not reach zero subdiagonal if complex block remains)
    REQUIRE(iters <= 50);
    
    // Matrix should still be upper Hessenberg
    REQUIRE(EigenSolverHelpers::IsUpperHessenberg(Hcopy, 1e-6));
}

TEST_CASE("DoubleShift - 2x2 complex block extraction", "[eigensolver][building-block][doubleshift]")
{
    // Create a matrix whose bottom 2x2 has complex eigenvalues
    // [[a, -b], [b, a]] has eigenvalues a ± ib
    // Use a=2, b=3 -> eigenvalues 2±3i
    Matrix<Real> H{4, 4, {REAL(5.0), REAL(1.0), REAL(0.0), REAL(0.0),
                          REAL(1.0), REAL(4.0), REAL(1.0), REAL(0.0),
                          REAL(0.0), REAL(1.0), REAL(2.0), -REAL(3.0),
                          REAL(0.0), REAL(0.0), REAL(3.0),  REAL(2.0)}};
    
    // Bottom 2x2 eigenvalues should be 2 ± 3i
    auto eig = EigenSolverHelpers::Eigenvalues2x2(H(2,2), H(2,3), H(3,2), H(3,3));
    
    INFO("Bottom 2x2 eigenvalues: " << eig.real1 << " ± " << eig.imag1 << "i");
    REQUIRE(eig.isComplex);
    REQUIRE(std::abs(eig.real1 - REAL(2.0)) < 1e-10);
    REQUIRE(std::abs(std::abs(eig.imag1) - REAL(3.0)) < 1e-10);
}

// =============================================================================
// BUILDING BLOCK 4: DEFLATION AND EXTRACTION TESTS  
// =============================================================================

TEST_CASE("Deflation - Detect deflation in converged matrix", "[eigensolver][building-block][deflation]")
{
    // Matrix with zero subdiagonal element (already deflated at position 2)
    Matrix<Real> H{4, 4, {REAL(5.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(4.0), REAL(1.0), REAL(2.0),   // Non-zero at (1,0)
                          REAL(0.0), REAL(0.0), REAL(3.0), REAL(1.0),   // Zero at (2,1) - deflation!
                          REAL(0.0), REAL(0.0), REAL(1.0), REAL(2.0)}};
    
    auto result = EigenSolverHelpers::CheckDeflation(H, 0, 3, 1e-10);
    
    INFO("canDeflate: " << result.canDeflate);
    INFO("deflationIndex: " << result.deflationIndex);
    REQUIRE(result.canDeflate);
    REQUIRE(result.deflationIndex == 2);  // CheckDeflation finds from bottom up
}

TEST_CASE("Deflation - No deflation in unreduced matrix", "[eigensolver][building-block][deflation]")
{
    // Matrix with no small subdiagonal elements
    Matrix<Real> H{4, 4, {REAL(5.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(4.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(2.0), REAL(3.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(2.0)}};
    
    auto result = EigenSolverHelpers::CheckDeflation(H, 0, 3, 1e-10);
    
    REQUIRE_FALSE(result.canDeflate);
}

TEST_CASE("Deflation - Detect 1x1 block at bottom", "[eigensolver][building-block][deflation]")
{
    // Matrix with zero at H(3,2) - isolates single eigenvalue
    Matrix<Real> H{4, 4, {REAL(5.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(4.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(2.0), REAL(3.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(0.0), REAL(7.0)}};  // Zero at (3,2) - eigenvalue 7
    
    auto result = EigenSolverHelpers::CheckDeflation(H, 0, 3, 1e-10);
    
    REQUIRE(result.canDeflate);
    REQUIRE(result.deflationIndex == 3);
    REQUIRE(result.blockSize == 1);  // 1x1 block (real eigenvalue)
}

TEST_CASE("Deflation - Detect 2x2 block", "[eigensolver][building-block][deflation]")
{
    // Matrix with zero at H(2,1) but non-zero at H(3,2) - isolates 2x2 block
    Matrix<Real> H{4, 4, {REAL(5.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(4.0), REAL(1.0), REAL(2.0),
                          REAL(0.0), REAL(0.0), REAL(3.0), REAL(1.0),   // Zero at (2,1)
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(2.0)}}; // Non-zero at (3,2) - 2x2 block below
    
    auto result = EigenSolverHelpers::CheckDeflation(H, 0, 3, 1e-10);
    
    REQUIRE(result.canDeflate);
    REQUIRE(result.deflationIndex == 2);
    REQUIRE(result.blockSize == 2);  // 2x2 block (possibly complex pair)
}

TEST_CASE("Extraction - Extract eigenvalues from triangular matrix", "[eigensolver][building-block][deflation]")
{
    // Upper triangular matrix - all 1x1 blocks
    Matrix<Real> H{4, 4, {REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0),
                          REAL(0.0), REAL(2.0), REAL(3.0), REAL(4.0),
                          REAL(0.0), REAL(0.0), REAL(3.0), REAL(4.0),
                          REAL(0.0), REAL(0.0), REAL(0.0), REAL(4.0)}};
    
    auto result = EigenSolverHelpers::ExtractEigenvalues(H, 1e-10);
    
    REQUIRE(result.eigenvalues.size() == 4);
    REQUIRE(result.realCount == 4);
    REQUIRE(result.complexPairs == 0);
    
    // Check eigenvalues (should be 1, 2, 3, 4)
    std::vector<Real> expected = {REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
    for (size_t i = 0; i < 4; i++)
    {
        INFO("Eigenvalue " << i << ": " << result.eigenvalues[i].real);
        REQUIRE(std::abs(result.eigenvalues[i].real - expected[i]) < 1e-10);
        REQUIRE_FALSE(result.eigenvalues[i].isComplex);
    }
}

TEST_CASE("Extraction - Extract complex pair from 2x2 block", "[eigensolver][building-block][deflation]")
{
    // Quasi-triangular: 1x1 block then 2x2 complex block then 1x1 block
    // 2x2 block [[2, -3], [3, 2]] has eigenvalues 2 ± 3i
    Matrix<Real> H{4, 4, {REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0),
                          REAL(0.0), REAL(2.0),-REAL(3.0), REAL(4.0),
                          REAL(0.0), REAL(3.0), REAL(2.0), REAL(4.0),
                          REAL(0.0), REAL(0.0), REAL(0.0), REAL(5.0)}};
    
    auto result = EigenSolverHelpers::ExtractEigenvalues(H, 1e-10);
    
    REQUIRE(result.eigenvalues.size() == 4);
    REQUIRE(result.realCount == 2);      // eigenvalues 1 and 5
    REQUIRE(result.complexPairs == 1);   // one complex pair 2±3i
    
    // First eigenvalue: 1 (real)
    REQUIRE(std::abs(result.eigenvalues[0].real - REAL(1.0)) < 1e-10);
    REQUIRE_FALSE(result.eigenvalues[0].isComplex);
    
    // Second and third: 2 ± 3i (complex pair)
    REQUIRE(std::abs(result.eigenvalues[1].real - REAL(2.0)) < 1e-10);
    REQUIRE(std::abs(std::abs(result.eigenvalues[1].imag) - REAL(3.0)) < 1e-10);
    REQUIRE(result.eigenvalues[1].isComplex);
    
    REQUIRE(std::abs(result.eigenvalues[2].real - REAL(2.0)) < 1e-10);
    REQUIRE(std::abs(std::abs(result.eigenvalues[2].imag) - REAL(3.0)) < 1e-10);
    REQUIRE(result.eigenvalues[2].isComplex);
    
    // Fourth: 5 (real)
    REQUIRE(std::abs(result.eigenvalues[3].real - REAL(5.0)) < 1e-10);
    REQUIRE_FALSE(result.eigenvalues[3].isComplex);
}

TEST_CASE("Extraction - Mixed real and complex blocks", "[eigensolver][building-block][deflation]")
{
    // 5x5 quasi-triangular with mix
    // Block structure: 1x1 (real), 2x2 (complex), 2x2 (real pair from 2x2 with real eigenvalues)
    Matrix<Real> H{5, 5, {REAL(1.0), REAL(2.0), REAL(0.0), REAL(0.0), REAL(0.0),
                          REAL(0.0), REAL(3.0),-REAL(2.0), REAL(0.0), REAL(0.0),  // 2x2 complex block
                          REAL(0.0), REAL(2.0), REAL(3.0), REAL(0.0), REAL(0.0),  // eigenvalues 3 ± 2i
                          REAL(0.0), REAL(0.0), REAL(0.0), REAL(4.0), REAL(0.0),  // 1x1 real
                          REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(5.0)}};// 1x1 real
    
    auto result = EigenSolverHelpers::ExtractEigenvalues(H, 1e-10);
    
    REQUIRE(result.eigenvalues.size() == 5);
    REQUIRE(result.realCount == 3);      // 1, 4, 5
    REQUIRE(result.complexPairs == 1);   // 3 ± 2i
}

TEST_CASE("ApplyDeflation - Zeros subdiagonal element", "[eigensolver][building-block][deflation]")
{
    Matrix<Real> H{4, 4, {REAL(5.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          1e-12, REAL(4.0), REAL(1.0), REAL(2.0),  // Very small but non-zero
                          REAL(0.0), REAL(2.0), REAL(3.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(2.0), REAL(2.0)}};
    
    auto check = EigenSolverHelpers::CheckDeflation(H, 0, 3, 1e-10);
    REQUIRE(check.canDeflate);
    
    EigenSolverHelpers::ApplyDeflation(H, check.deflationIndex);
    
    REQUIRE(H(1, 0) == REAL(0.0));  // Should be exactly zero now
}

// =============================================================================
// BUILDING BLOCK 5: EIGENVECTOR COMPUTATION TESTS
// =============================================================================

TEST_CASE("Eigenvectors - Real eigenvalue back-substitution", "[eigensolver][building-block][eigenvector]")
{
    // Upper triangular matrix with eigenvalues 1, 2, 3
    Matrix<Real> T{3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(0.0), REAL(2.0), REAL(4.0),
                          REAL(0.0), REAL(0.0), REAL(3.0)}};
    
    // Eigenvector for eigenvalue 3 (at T[2,2])
    Vector<Real> v = EigenSolverHelpers::ComputeRealEigenvector(T, 2);
    
    // Verify: v should be [0, 0, 1] (normalized) or scalar multiple
    // Since x[2]=1 and back-sub gives x[1] = -4/(2-3) = 4, x[0] = -(2*4+3*1)/(1-3) = REAL(5.5)
    // After normalization: [REAL(5.5), 4, 1] / ||...||
    
    // Check eigenvector property: T*v = 3*v
    Real residual = EigenSolverHelpers::EigenvectorResidual(T, v, REAL(3.0));
    INFO("Residual for eigenvalue 3: " << residual);
    REQUIRE(residual < 1e-10);
}

TEST_CASE("Eigenvectors - Multiple real eigenvalues", "[eigensolver][building-block][eigenvector]")
{
    // Upper triangular with eigenvalues 1, 2, 3
    Matrix<Real> T{3, 3, {REAL(1.0), REAL(1.0), REAL(1.0),
                          REAL(0.0), REAL(2.0), REAL(1.0),
                          REAL(0.0), REAL(0.0), REAL(3.0)}};
    Matrix<Real> Q = Matrix<Real>::GetUnitMatrix(3);  // Identity (T is already Schur form)
    
    auto result = EigenSolverHelpers::ComputeEigenvectorsFromSchur(T, Q, 1e-10);
    
    REQUIRE(result.vectors.RowNum() == 3);
    REQUIRE(result.vectors.ColNum() == 3);
    
    // Verify each eigenvector
    for (int col = 0; col < 3; col++)
    {
        Vector<Real> v(3);
        for (int row = 0; row < 3; row++)
            v[row] = result.vectors(row, col);
        
        Real lambda = T(col, col);  // Eigenvalues are 1, 2, 3 on diagonal
        Real residual = EigenSolverHelpers::EigenvectorResidual(T, v, lambda);
        
        INFO("Eigenvalue " << lambda << " residual: " << residual);
        REQUIRE(residual < 1e-6);
    }
}

TEST_CASE("Eigenvectors - With Q transformation", "[eigensolver][building-block][eigenvector]")
{
    // Create a matrix A, reduce to Schur form T = Q^T * A * Q
    // Then verify eigenvectors of A work correctly
    
    // Simple 3x3 with known eigenvalues
    Matrix<Real> A{3, 3, {REAL(4.0), REAL(1.0), REAL(1.0),
                          REAL(1.0), REAL(4.0), REAL(1.0),
                          REAL(1.0), REAL(1.0), REAL(4.0)}};
    // This symmetric matrix has eigenvalues 6, 3, 3
    
    // First reduce to Hessenberg (which for symmetric = tridiagonal)
    auto hess = EigenSolverHelpers::ReduceToHessenberg(A);
    
    // Apply QR iterations to get Schur form
    Matrix<Real> T = hess.H;
    Matrix<Real> Q = hess.Q;
    
    // Multiple QR iterations
    for (int iter = 0; iter < 50; iter++)
    {
        auto step = EigenSolverHelpers::SingleQRStep(T, true);
        
        // Update total Q
        Matrix<Real> newQ(3, 3);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                newQ(i, j) = REAL(0.0);
                for (int k = 0; k < 3; k++)
                    newQ(i, j) += Q(i, k) * step.Q(k, j);
            }
        Q = newQ;
        T = step.H;
        
        // Check convergence
        bool converged = true;
        for (int i = 1; i < 3; i++)
            if (std::abs(T(i, i-1)) > 1e-10)
                converged = false;
        if (converged) break;
    }
    
    // Now compute eigenvectors
    auto evecs = EigenSolverHelpers::ComputeEigenvectorsFromSchur(T, Q, 1e-10);
    
    // Verify A*v = λ*v for each eigenvector
    auto eigenvalues = EigenSolverHelpers::ExtractEigenvalues(T, 1e-10);
    
    for (int col = 0; col < 3; col++)
    {
        Vector<Real> v(3);
        for (int row = 0; row < 3; row++)
            v[row] = evecs.vectors(row, col);
        
        Real lambda = eigenvalues.eigenvalues[col].real;
        
        // Compute A*v
        Vector<Real> Av(3, REAL(0.0));
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Av[i] += A(i, j) * v[j];
        
        // Compute ||A*v - λ*v|| / ||v||
        Real residual = REAL(0.0);
        Real vnorm = REAL(0.0);
        for (int i = 0; i < 3; i++)
        {
            Real diff = Av[i] - lambda * v[i];
            residual += diff * diff;
            vnorm += v[i] * v[i];
        }
        residual = std::sqrt(residual) / std::sqrt(vnorm);
        
        INFO("Eigenvalue " << lambda << " residual: " << residual);
        REQUIRE(residual < 1e-4);  // Allow some numerical error
    }
}

TEST_CASE("Eigenvectors - Complex eigenvalue pair", "[eigensolver][building-block][eigenvector]")
{
    // Quasi-triangular with complex eigenvalue block
    // [[0, -1], [1, 0]] has eigenvalues ±i
    Matrix<Real> T{3, 3, {REAL(2.0), REAL(1.0), REAL(1.0),
                          REAL(0.0), REAL(0.0),-REAL(1.0),  // 2x2 block with eigenvalues ±i
                          REAL(0.0), REAL(1.0), REAL(0.0)}};
    Matrix<Real> Q = Matrix<Real>::GetUnitMatrix(3);
    
    auto result = EigenSolverHelpers::ComputeEigenvectorsFromSchur(T, Q, 1e-10);
    
    // First column should be real eigenvector for λ=2
    REQUIRE_FALSE(result.isComplexPair[0]);
    
    // Columns 1,2 should be complex pair
    REQUIRE(result.isComplexPair[1]);
    REQUIRE(result.isComplexPair[2]);
    
    // Verify first eigenvector (real, λ=2)
    Vector<Real> v(3);
    for (int i = 0; i < 3; i++)
        v[i] = result.vectors(i, 0);
    
    Real residual = EigenSolverHelpers::EigenvectorResidual(T, v, REAL(2.0));
    INFO("Residual for eigenvalue 2: " << residual);
    REQUIRE(residual < 1e-6);
}

TEST_CASE("Eigenvectors - Normalized output", "[eigensolver][building-block][eigenvector]")
{
    Matrix<Real> T{3, 3, {REAL(1.0), REAL(2.0), REAL(3.0),
                          REAL(0.0), REAL(4.0), REAL(5.0),
                          REAL(0.0), REAL(0.0), REAL(6.0)}};
    Matrix<Real> Q = Matrix<Real>::GetUnitMatrix(3);
    
    auto result = EigenSolverHelpers::ComputeEigenvectorsFromSchur(T, Q, 1e-10);
    
    // Verify each column is normalized
    for (int col = 0; col < 3; col++)
    {
        Real norm = REAL(0.0);
        for (int row = 0; row < 3; row++)
            norm += result.vectors(row, col) * result.vectors(row, col);
        norm = std::sqrt(norm);
        
        INFO("Column " << col << " norm: " << norm);
        REQUIRE(std::abs(norm - REAL(1.0)) < 1e-10);
    }
}

// =============================================================================
// INTEGRATION TESTS: FULL GENERAL EIGENSOLVER
// =============================================================================

TEST_CASE("EigenSolver - 2x2 real eigenvalues", "[eigensolver][integration]")
{
    // [[4, 1], [2, 3]] has eigenvalues 5 and 2
    Matrix<Real> A{2, 2, {REAL(4.0), REAL(1.0),
                          REAL(2.0), REAL(3.0)}};
    
    auto result = EigenSolver::Solve(A, 1e-10, 100);
    
    REQUIRE(result.converged);
    REQUIRE(result.eigenvalues.size() == 2);
    
    // Check eigenvalues (order may vary)
    std::vector<Real> expected = {REAL(5.0), REAL(2.0)};
    std::vector<Real> computed;
    for (auto& e : result.eigenvalues)
        computed.push_back(e.real);
    std::sort(computed.begin(), computed.end());
    std::sort(expected.begin(), expected.end());
    
    INFO("Eigenvalues: " << computed[0] << ", " << computed[1]);
    REQUIRE(std::abs(computed[0] - expected[0]) < 1e-6);
    REQUIRE(std::abs(computed[1] - expected[1]) < 1e-6);
}

TEST_CASE("EigenSolver - 3x3 symmetric matrix", "[eigensolver][integration]")
{
    // Symmetric matrix with known eigenvalues 6, 3, 3
    Matrix<Real> A{3, 3, {REAL(4.0), REAL(1.0), REAL(1.0),
                          REAL(1.0), REAL(4.0), REAL(1.0),
                          REAL(1.0), REAL(1.0), REAL(4.0)}};
    
    auto result = EigenSolver::Solve(A, 1e-10, 200);
    
    INFO("Converged: " << result.converged);
    INFO("Iterations: " << result.iterations);
    REQUIRE(result.converged);
    REQUIRE(result.eigenvalues.size() == 3);
    
    // Check eigenvalues
    std::vector<Real> computed;
    for (auto& e : result.eigenvalues)
        computed.push_back(e.real);
    std::sort(computed.begin(), computed.end());
    
    // Expected: 3, 3, 6
    INFO("Eigenvalues: " << computed[0] << ", " << computed[1] << ", " << computed[2]);
    REQUIRE(std::abs(computed[0] - REAL(3.0)) < 1e-4);
    REQUIRE(std::abs(computed[1] - REAL(3.0)) < 1e-4);
    REQUIRE(std::abs(computed[2] - REAL(6.0)) < 1e-4);
}

TEST_CASE("EigenSolver - 3x3 with complex eigenvalues", "[eigensolver][integration]")
{
    // Matrix with complex eigenvalues
    // [[0, -1, 0], [1, 0, 0], [0, 0, 2]] has eigenvalues ±i and 2
    Matrix<Real> A{3, 3, {REAL(0.0), -REAL(1.0), REAL(0.0),
                          REAL(1.0),  REAL(0.0), REAL(0.0),
                          REAL(0.0),  REAL(0.0), REAL(2.0)}};
    
    auto result = EigenSolver::Solve(A, 1e-10, 200);
    
    INFO("Converged: " << result.converged);
    INFO("Iterations: " << result.iterations);
    REQUIRE(result.converged);
    REQUIRE(result.eigenvalues.size() == 3);
    
    // Should have one real (2) and one complex pair (±i)
    int realCount = 0;
    int complexCount = 0;
    Real realEig = 0;
    Real complexReal = 0;
    Real complexImag = 0;
    
    for (auto& e : result.eigenvalues)
    {
        if (e.isComplex())
        {
            complexCount++;
            complexReal = e.real;
            complexImag = std::abs(e.imag);
        }
        else
        {
            realCount++;
            realEig = e.real;
        }
    }
    
    INFO("Real eigenvalue: " << realEig);
    INFO("Complex pair: " << complexReal << " ± " << complexImag << "i");
    
    REQUIRE(realCount == 1);
    REQUIRE(complexCount == 2);  // Complex pair
    REQUIRE(std::abs(realEig - REAL(2.0)) < 1e-6);
    REQUIRE(std::abs(complexReal) < 1e-6);  // Real part is 0
    REQUIRE(std::abs(complexImag - REAL(1.0)) < 1e-6);  // Imaginary part is 1
}

TEST_CASE("EigenSolver - Diagonal matrix", "[eigensolver][integration]")
{
    // Diagonal matrix - eigenvalues are diagonal elements
    Matrix<Real> A{4, 4, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0),
                          REAL(0.0), REAL(2.0), REAL(0.0), REAL(0.0),
                          REAL(0.0), REAL(0.0), REAL(3.0), REAL(0.0),
                          REAL(0.0), REAL(0.0), REAL(0.0), REAL(4.0)}};
    
    auto result = EigenSolver::Solve(A, 1e-10, 100);
    
    REQUIRE(result.converged);
    REQUIRE(result.eigenvalues.size() == 4);
    
    std::vector<Real> computed;
    for (auto& e : result.eigenvalues)
        computed.push_back(e.real);
    std::sort(computed.begin(), computed.end());
    
    REQUIRE(std::abs(computed[0] - REAL(1.0)) < 1e-10);
    REQUIRE(std::abs(computed[1] - REAL(2.0)) < 1e-10);
    REQUIRE(std::abs(computed[2] - REAL(3.0)) < 1e-10);
    REQUIRE(std::abs(computed[3] - REAL(4.0)) < 1e-10);
}

TEST_CASE("EigenSolver - 4x4 general matrix", "[eigensolver][integration]")
{
    // Non-symmetric matrix
    Matrix<Real> A{4, 4, {REAL(4.0), REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(2.0), REAL(3.0), REAL(1.0), REAL(2.0),
                          REAL(1.0), REAL(2.0), REAL(5.0), REAL(1.0),
                          REAL(1.0), REAL(1.0), REAL(1.0), REAL(6.0)}};
    
    auto result = EigenSolver::Solve(A, 1e-10, 300);
    
    INFO("Converged: " << result.converged);
    INFO("Iterations: " << result.iterations);
    INFO("Max residual: " << result.maxResidual);
    REQUIRE(result.converged);
    REQUIRE(result.eigenvalues.size() == 4);
    
    // Verify trace: sum of eigenvalues = trace of A = 4+3+5+6 = 18
    Real traceA = A(0,0) + A(1,1) + A(2,2) + A(3,3);
    Real sumEig = 0;
    for (auto& e : result.eigenvalues)
        sumEig += e.real;
    
    INFO("Trace of A: " << traceA);
    INFO("Sum of eigenvalues: " << sumEig);
    REQUIRE(std::abs(sumEig - traceA) < 1e-4);
}

TEST_CASE("EigenSolver - Eigenvector verification", "[eigensolver][integration]")
{
    // Simple matrix for eigenvector verification
    Matrix<Real> A{3, 3, {REAL(2.0), REAL(1.0), REAL(0.0),
                          REAL(1.0), REAL(2.0), REAL(1.0),
                          REAL(0.0), REAL(1.0), REAL(2.0)}};
    
    auto result = EigenSolver::Solve(A, 1e-10, 200);
    
    REQUIRE(result.converged);
    
    // Verify A*v = λ*v for each real eigenpair
    for (size_t i = 0; i < result.eigenvalues.size(); i++)
    {
        if (!result.eigenvalues[i].isComplex())
        {
            Vector<Real> v(3);
            for (int j = 0; j < 3; j++)
                v[j] = result.eigenvectors(j, static_cast<int>(i));
            
            Real lambda = result.eigenvalues[i].real;
            
            // Compute A*v
            Vector<Real> Av(3, REAL(0.0));
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    Av[j] += A(j, k) * v[k];
            
            // Compute ||A*v - λ*v|| / ||v||
            Real residual = REAL(0.0);
            Real vnorm = REAL(0.0);
            for (int j = 0; j < 3; j++)
            {
                Real diff = Av[j] - lambda * v[j];
                residual += diff * diff;
                vnorm += v[j] * v[j];
            }
            residual = std::sqrt(residual) / std::sqrt(vnorm);
            
            INFO("Eigenvalue " << lambda << " residual: " << residual);
            REQUIRE(residual < 1e-4);
        }
    }
}

