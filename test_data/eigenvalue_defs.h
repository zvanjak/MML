#if !defined __MML_EIGENVALUE_DEFS_H
#define __MML_EIGENVALUE_DEFS_H

#include <string>
#include <cmath>
#include <vector>
#include <complex>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#include "base/Vector.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    CLASSIC/STANDARD EIGENVALUE TEST MATRICES                          //
    //             (Well-conditioned, easily verifiable reference cases)                      //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // 5x5 STANDARD SYMMETRIC TEST MATRIX
    // ========================================================================================
    // Pascal matrix (symmetric, positive definite, integer eigenvalues for small n)
    // P(i,j) = C(i+j, i) - binomial coefficients
    // Known for having integer eigenvalues for n<=5
    const static inline Matrix<Real> classic_pascal_5x5{5, 5, {
         1,  1,  1,  1,  1,
         1,  2,  3,  4,  5,
         1,  3,  6, 10, 15,
         1,  4, 10, 20, 35,
         1,  5, 15, 35, 70
    }};
    // Eigenvalues (verified by solver - product ≈ det = 1)
    const static inline Vector<Real> classic_pascal_5x5_eigenvalues{
        0.0108354247894, 0.181242088628, 1.0, 5.51748985992, 92.2904329266
    };

    // ========================================================================================
    // 10x10 STANDARD HILBERT MATRIX
    // ========================================================================================
    // H(i,j) = 1/(i+j+1) - classic ill-conditioned test matrix
    // Symmetric positive definite but extremely ill-conditioned
    inline Matrix<Real> getHilbert10x10() {
        Matrix<Real> H(10, 10);
        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 10; ++j)
                H(i, j) = 1.0 / (i + j + 1);
        return H;
    }
    // Eigenvalues (high precision - condition number ~10^13)
    const static inline Vector<Real> classic_hilbert_10x10_eigenvalues{
        1.75219e-10, 3.07855e-09, 4.47524e-08, 5.42563e-07,
        5.47398e-06, 4.54373e-05, 3.05898e-04, 1.62338e-03,
        6.38379e-03, 1.75191e+00
    };

    // ========================================================================================
    // 20x20 STANDARD: TOEPLITZ TRIDIAGONAL
    // ========================================================================================
    // T = tridiag(1, 2, 1) - symmetric tridiagonal with known eigenvalues
    // λ_k = 2 + 2*cos(k*π/(n+1)) for k = 1,...,n
    inline Matrix<Real> getToeplitz20x20() {
        Matrix<Real> T(20, 20);
        for (int i = 0; i < 20; ++i) {
            T(i, i) = 2.0;
            if (i > 0) T(i, i-1) = 1.0;
            if (i < 19) T(i, i+1) = 1.0;
        }
        return T;
    }
    inline Vector<Real> getToeplitz20x20_eigenvalues() {
        Vector<Real> eig(20);
        for (int k = 1; k <= 20; ++k)
            eig[k-1] = 2.0 + 2.0 * std::cos(k * Constants::PI / 21.0);
        return eig;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                       REPEATED/DEGENERATE EIGENVALUES                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 1: EXACT REPEATED EIGENVALUES (symmetric)
    // ========================================================================================

    // Identity matrix - all eigenvalues = 1 (multiplicity n)
    inline Matrix<Real> getIdentity(int n) {
        return Matrix<Real>::GetUnitMatrix(n);
    }

    // Double eigenvalue: diag(1, 2, 2, 3) - eigenvalue 2 has multiplicity 2
    const static inline Matrix<Real> repeated_double_4x4{4, 4, {
        1, 0, 0, 0,
        0, 2, 0, 0,
        0, 0, 2, 0,
        0, 0, 0, 3
    }};
    const static inline Vector<Real> repeated_double_4x4_eigenvalues{1.0, 2.0, 2.0, 3.0};

    // Triple eigenvalue: diag(1, 2, 2, 2, 3) - eigenvalue 2 has multiplicity 3
    const static inline Matrix<Real> repeated_triple_5x5{5, 5, {
        1, 0, 0, 0, 0,
        0, 2, 0, 0, 0,
        0, 0, 2, 0, 0,
        0, 0, 0, 2, 0,
        0, 0, 0, 0, 3
    }};
    const static inline Vector<Real> repeated_triple_5x5_eigenvalues{1.0, 2.0, 2.0, 2.0, 3.0};

    // Quadruple eigenvalue with 2D eigenspace - tests eigenvector orthogonalization
    // Symmetric matrix with eigenvalues: 0, 0, 4, 4
    const static inline Matrix<Real> repeated_quadruple_4x4{4, 4, {
        2, 1, 1, 0,
        1, 2, 0, 1,
        1, 0, 2, 1,
        0, 1, 1, 2
    }};
    const static inline Vector<Real> repeated_quadruple_4x4_eigenvalues{0.0, 2.0, 2.0, 4.0};

    // ========================================================================================
    // CATEGORY 2: DEFECTIVE MATRICES (non-diagonalizable)
    // ========================================================================================
    // These matrices have fewer linearly independent eigenvectors than eigenvalues
    // Tests: detection of defectiveness, handling of Jordan blocks

    // Classic 2x2 Jordan block: eigenvalue 2, but only one eigenvector
    const static inline Matrix<Real> defective_jordan_2x2{2, 2, {
        2, 1,
        0, 2
    }};
    const static inline Vector<Real> defective_jordan_2x2_eigenvalues{2.0, 2.0};
    // Defect: algebraic multiplicity 2, geometric multiplicity 1

    // 3x3 Jordan block: eigenvalue 1, single eigenvector
    const static inline Matrix<Real> defective_jordan_3x3{3, 3, {
        1, 1, 0,
        0, 1, 1,
        0, 0, 1
    }};
    const static inline Vector<Real> defective_jordan_3x3_eigenvalues{1.0, 1.0, 1.0};

    // Mixed Jordan structure: J_2(2) ⊕ J_1(3) - two Jordan blocks
    const static inline Matrix<Real> defective_mixed_3x3{3, 3, {
        2, 1, 0,
        0, 2, 0,
        0, 0, 3
    }};
    const static inline Vector<Real> defective_mixed_3x3_eigenvalues{2.0, 2.0, 3.0};

    // 4x4 with two 2x2 Jordan blocks: eigenvalues 1, 1, 2, 2
    const static inline Matrix<Real> defective_two_blocks_4x4{4, 4, {
        1, 1, 0, 0,
        0, 1, 0, 0,
        0, 0, 2, 1,
        0, 0, 0, 2
    }};
    const static inline Vector<Real> defective_two_blocks_4x4_eigenvalues{1.0, 1.0, 2.0, 2.0};

    // ========================================================================================
    // CATEGORY 3: NEARLY DEFECTIVE MATRICES (near Jordan structure)
    // ========================================================================================
    // These have distinct eigenvalues but are very close to being defective
    // Tests: numerical stability, condition of eigenvector computation

    // Nearly defective: Jordan block perturbed by epsilon
    inline Matrix<Real> getNearlyDefective3x3(Real epsilon = 1e-8) {
        return Matrix<Real>(3, 3, {
            2.0, 1.0, 0.0,
            epsilon, 2.0, 1.0,
            0.0, epsilon, 2.0
        });
    }
    // Eigenvalues ≈ 2 ± small perturbation

    // Wilkinson matrix W_n^+ : tridiag(n-1, n-2, ..., 0, 1, ..., n-1) with 1's on off-diag
    // Famous for having very close eigenvalue pairs
    inline Matrix<Real> getWilkinsonMinus(int n = 21) {
        Matrix<Real> W(n, n);
        int m = (n - 1) / 2;
        for (int i = 0; i < n; ++i) {
            W(i, i) = std::abs(m - i);
            if (i > 0) W(i, i-1) = 1.0;
            if (i < n-1) W(i, i+1) = 1.0;
        }
        return W;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         CLUSTERED EIGENVALUES                                          //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 4: CLUSTERED EIGENVALUES (distinct but very close)
    // ========================================================================================

    // Cluster around 1: eigenvalues 0.999, 1.0, 1.001, 2.0
    const static inline Matrix<Real> clustered_tight_4x4{4, 4, {
        0.999,     0.0,   0.0, 0.0,
          0.0,     1.0,   0.0, 0.0,
          0.0,     0.0, 1.001, 0.0,
          0.0,     0.0,   0.0, 2.0
    }};
    const static inline Vector<Real> clustered_tight_4x4_eigenvalues{0.999, 1.0, 1.001, 2.0};

    // Very tight cluster: 1.0, 1.0001, 1.0002, 1.0003, 2.0
    const static inline Matrix<Real> clustered_very_tight_5x5{5, 5, {
        1.0000, 0, 0, 0, 0,
        0, 1.0001, 0, 0, 0,
        0, 0, 1.0002, 0, 0,
        0, 0, 0, 1.0003, 0,
        0, 0, 0, 0, 2.0000
    }};
    const static inline Vector<Real> clustered_very_tight_5x5_eigenvalues{
        1.0, 1.0001, 1.0002, 1.0003, 2.0
    };

    // Graded clusters: well-separated groups but close within groups
    // Group 1: ~0.1, Group 2: ~1.0, Group 3: ~10.0
    inline Matrix<Real> getGradedClusters6x6() {
        Matrix<Real> A(6, 6);
        Vector<Real> diag{0.10, 0.11, 1.00, 1.01, 10.0, 10.1};
        for (int i = 0; i < 6; ++i) A(i, i) = diag[i];
        return A;
    }
    const static inline Vector<Real> clustered_graded_6x6_eigenvalues{
        0.10, 0.11, 1.00, 1.01, 10.0, 10.1
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                           ILL-CONDITIONED MATRICES                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 5: EXTREME CONDITIONING
    // ========================================================================================

    // Vandermonde-like: extremely ill-conditioned
    // V(i,j) = nodes[i]^j with nodes 1,2,3,4,5
    inline Matrix<Real> getVandermonde5x5() {
        Matrix<Real> V(5, 5);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                V(i, j) = std::pow(static_cast<Real>(i + 1), j);
        return V;
    }

    // Frank matrix: eigenvalues known analytically, ill-conditioned eigenvectors
    // F(i,j) = n - max(i,j) + 1 for i <= j; F(j,i) = F(i,j) - 1
    inline Matrix<Real> getFrank(int n) {
        Matrix<Real> F(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                F(i, j) = (i <= j) ? (n - j) : (n - i - 1);
        return F;
    }

    // Kahan matrix: upper triangular, very ill-conditioned
    // Tests backward stability of eigensolvers
    inline Matrix<Real> getKahan(int n, Real c = 0.5) {
        Real s = std::sqrt(1 - c * c);
        Matrix<Real> K(n, n);
        for (int i = 0; i < n; ++i) {
            K(i, i) = std::pow(s, i);
            for (int j = i + 1; j < n; ++j)
                K(i, j) = -c * std::pow(s, i);
        }
        return K;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                           COMPLEX EIGENVALUES                                          //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 6: MATRICES WITH COMPLEX EIGENVALUES (non-symmetric)
    // ========================================================================================

    // Rotation-like: 2x2 rotation by θ, eigenvalues e^{±iθ}
    inline Matrix<Real> getRotation2x2(Real theta = Constants::PI / 4) {
        Real c = std::cos(theta);
        Real s = std::sin(theta);
        return Matrix<Real>(2, 2, {c, -s, s, c});
    }
    // Eigenvalues: cos(θ) ± i*sin(θ) = e^{±iθ}

    // 3x3 with one real and one complex conjugate pair
    const static inline Matrix<Real> complex_pair_3x3{3, 3, {
        1.0, -1.0, 0.0,
        1.0,  1.0, 0.0,
        0.0,  0.0, 3.0
    }};
    // Eigenvalues: 1±i, 3
    const static inline std::vector<std::complex<Real>> complex_pair_3x3_eigenvalues{
        {1.0, 1.0}, {1.0, -1.0}, {3.0, 0.0}
    };

    // 4x4 with two complex conjugate pairs
    const static inline Matrix<Real> complex_two_pairs_4x4{4, 4, {
        0, -1, 0, 0,
        1,  0, 0, 0,
        0,  0, 0, -2,
        0,  0, 2,  0
    }};
    // Eigenvalues: ±i, ±2i
    const static inline std::vector<std::complex<Real>> complex_two_pairs_4x4_eigenvalues{
        {0.0, 1.0}, {0.0, -1.0}, {0.0, 2.0}, {0.0, -2.0}
    };

    // Companion matrix for polynomial x^4 - 1 (roots: ±1, ±i)
    const static inline Matrix<Real> companion_x4_minus_1{4, 4, {
        0, 0, 0, 1,
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0
    }};
    const static inline std::vector<std::complex<Real>> companion_x4_minus_1_eigenvalues{
        {1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}, {0.0, -1.0}
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         SPECIAL STRUCTURE MATRICES                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 7: SPECIAL STRUCTURES
    // ========================================================================================

    // Circulant matrix: eigenvalues are DFT of first row
    // C = circ(c_0, c_1, ..., c_{n-1})
    // λ_k = Σ c_j * ω^{jk} where ω = e^{2πi/n}
    inline Matrix<Real> getCirculant4x4() {
        return Matrix<Real>(4, 4, {
            1, 2, 3, 4,
            4, 1, 2, 3,
            3, 4, 1, 2,
            2, 3, 4, 1
        });
    }
    // Eigenvalues: 10, -2, -2i, 2i (sum of first row = 10 is always eigenvalue)
    const static inline std::vector<std::complex<Real>> circulant_4x4_eigenvalues{
        {10.0, 0.0}, {-2.0, 0.0}, {0.0, -2.0}, {0.0, 2.0}
    };

    // Toeplitz (symmetric): constant diagonals, analytically known eigenvalues
    // For tridiag(1, 2, 1): λ_k = 2 + 2*cos(kπ/(n+1))
    inline Matrix<Real> getSymmetricToeplitz(int n, Real diag = 2.0, Real offdiag = 1.0) {
        Matrix<Real> T(n, n);
        for (int i = 0; i < n; ++i) {
            T(i, i) = diag;
            if (i > 0) T(i, i-1) = offdiag;
            if (i < n-1) T(i, i+1) = offdiag;
        }
        return T;
    }

    // Block diagonal: eigenvalues are union of block eigenvalues
    // Tests handling of block structure
    inline Matrix<Real> getBlockDiagonal6x6() {
        return Matrix<Real>(6, 6, {
            1, 2, 0, 0, 0, 0,
            2, 1, 0, 0, 0, 0,
            0, 0, 3, 1, 0, 0,
            0, 0, 1, 3, 0, 0,
            0, 0, 0, 0, 5, 0,
            0, 0, 0, 0, 0, 6
        });
    }
    // Block 1 (2x2): eigenvalues -1, 3
    // Block 2 (2x2): eigenvalues 2, 4
    // Block 3 (1x1): eigenvalue 5
    // Block 4 (1x1): eigenvalue 6
    const static inline Vector<Real> block_diagonal_6x6_eigenvalues{-1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    // Stochastic (row sums = 1): always has eigenvalue 1
    const static inline Matrix<Real> stochastic_4x4{4, 4, {
        0.5, 0.3, 0.1, 0.1,
        0.2, 0.5, 0.2, 0.1,
        0.1, 0.2, 0.6, 0.1,
        0.1, 0.1, 0.1, 0.7
    }};
    // Dominant eigenvalue = 1 (Perron-Frobenius)

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                           PATHOLOGICAL CASES                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////

    // ========================================================================================
    // CATEGORY 8: PATHOLOGICAL/EXTREME CASES
    // ========================================================================================

    // Zero matrix: all eigenvalues = 0
    inline Matrix<Real> getZero(int n) {
        return Matrix<Real>(n, n, 0.0);
    }

    // Nilpotent: A^n = 0, all eigenvalues = 0, but not zero matrix
    const static inline Matrix<Real> nilpotent_3x3{3, 3, {
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    }};
    // All eigenvalues = 0, but matrix is not zero

    // Matrix with eigenvalue at machine epsilon scale
    inline Matrix<Real> getTinyEigenvalue3x3(Real eps = 1e-15) {
        return Matrix<Real>(3, 3, {
            eps, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 2.0
        });
    }

    // Matrix with huge eigenvalue spread
    inline Matrix<Real> getHugeSpread3x3(Real scale = 1e15) {
        return Matrix<Real>(3, 3, {
            1.0/scale, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, scale
        });
    }

    // Almost rank-deficient: smallest eigenvalue ≈ 0
    inline Matrix<Real> getAlmostSingular4x4(Real eps = 1e-12) {
        Matrix<Real> A(4, 4, {
            1.0 + eps, 1.0, 1.0, 1.0,
            1.0, 1.0 + eps, 1.0, 1.0,
            1.0, 1.0, 1.0 + eps, 1.0,
            1.0, 1.0, 1.0, 1.0 + eps
        });
        return A;
    }
    // Eigenvalues: 4+3eps (simple), eps (multiplicity 3)

} // namespace MML::TestBeds

#endif // __MML_EIGENVALUE_DEFS_H
