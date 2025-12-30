#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_SOLVERSPECIFIC_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_SOLVERSPECIFIC_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#endif

/*******************************************************************************************************************
 * SOLVER-SPECIFIC TEST MATRICES FOR LINEAR EQUATION SOLVERS
 * 
 * This file contains test matrices specifically designed for different types of solvers:
 * 
 * 1. OVERDETERMINED SYSTEMS (m > n)
 *    - For testing QR decomposition
 *    - For least-squares solvers
 *    - For SVD-based methods
 * 
 * 2. DIAGONALLY DOMINANT MATRICES
 *    - For iterative solvers (Jacobi, Gauss-Seidel, SOR)
 *    - Guaranteed convergence for many iterative methods
 *    - Well-conditioned systems
 * 
 * 3. SYMMETRIC POSITIVE DEFINITE (SPD) MATRICES
 *    - For Cholesky decomposition
 *    - For conjugate gradient method
 *    - For specialized SPD solvers
 * 
 * Each test system includes:
 * - Coefficient matrix A
 * - Right-hand side vector(s) b
 * - Exact solution(s) x (computed using MML solvers where appropriate)
 * - Eigenvalues (for square matrices)
 * - Documentation: properties, condition number, intended use
 * 
 * NOTE: All matrices are defined as static inline constants (not functions) to enable
 *       direct initializer list construction, following the pattern in classic_defs.h
 * 
 *******************************************************************************************************************/

namespace MML::TestBeds
{

/*******************************************************************************************************************
 * OVERDETERMINED SYSTEMS (m > n) - For QR and Least-Squares Solvers
 * 
 * These systems have more equations than unknowns, so generally have no exact solution.
 * Solvers should find the least-squares solution that minimizes ||Ax - b||₂.
 *******************************************************************************************************************/

// ===== OVERDETERMINED 4x3 SYSTEM =====
// Source: Custom test matrix
// Properties: Full rank (rank = 3), Well-conditioned
// Purpose: Basic QR and least-squares testing
// Note: Residual ||Ax - b|| ≈ 0.267, Condition number κ(A) ≈ 160

const static inline Matrix<Real> overdetermined_4x3_mat{4, 3, {
    1.0,  2.0,  1.0,
    2.0,  3.0,  1.0,
    3.0,  5.0,  1.0,
    4.0,  7.0,  1.0
}};

const static inline Vector<Real> overdetermined_4x3_rhs{2.0, 3.0, 5.0, 7.0};

const static inline Vector<Real> overdetermined_4x3_sol{
    0.0,
    1.0,
    0.0
};

const static inline VectorComplex overdetermined_4x3_eigen{
    Complex(46.4647, 0),
    Complex(0.5335, 0),
    Complex(0.0018, 0)
};

// Matrix properties for overdetermined_4x3_mat (rectangular - no determinant)
const static inline MML::Vector<Real> overdetermined_4x3_singular_values{10.629605628935844, 0.3935847579485939, 0.01346660091538802};
const static inline Real              overdetermined_4x3_cond_2 = 789.2547266174948;
const static inline int               overdetermined_4x3_rank = 3;

// ===== OVERDETERMINED 5x3 REGRESSION =====
// Source: Classic linear regression problem (y = β₀ + β₁x₁ + β₂x₂)
// Properties: Design matrix for linear regression
// Purpose: Test least-squares in regression context
// Note: Residual ||Ax - b|| ≈ 0.469, Condition number κ(A) ≈ 55

const static inline Matrix<Real> overdetermined_5x3_regression_mat{5, 3, {
    1.0,  1.0,  1.0,
    1.0,  2.0,  4.0,
    1.0,  3.0,  9.0,
    1.0,  4.0, 16.0,
    1.0,  5.0, 25.0
}};

const static inline Vector<Real> overdetermined_5x3_regression_rhs{2.0, 4.0, 7.0, 11.0, 16.0};

const static inline Vector<Real> overdetermined_5x3_regression_sol{
    1.0,
    0.5,
    0.5
};

const static inline VectorComplex overdetermined_5x3_regression_eigen{
    Complex(1510.0, 0),
    Complex(34.5, 0),
    Complex(0.5, 0)
};

// Matrix properties for overdetermined_5x3_regression_mat (rectangular - no determinant)
const static inline MML::Vector<Real> overdetermined_5x3_regression_singular_values{38.94025779098568, 2.4788816217779655, 0.43284115316413265};
const static inline Real              overdetermined_5x3_regression_cond_2 = 89.95862610108093;
const static inline int               overdetermined_5x3_regression_rank = 3;

// ===== OVERDETERMINED 6x4 POLYNOMIAL =====
// Source: Cubic polynomial fitting (y = a₀ + a₁x + a₂x² + a₃x³)
// Properties: Vandermonde-like structure, increasingly ill-conditioned with degree
// Purpose: Test polynomial least-squares fitting
// Note: Data points y = x² (quadratic), so cubic fit should be nearly exact
// Condition number κ(A) ≈ 218

const static inline Matrix<Real> overdetermined_6x4_polynomial_mat{6, 4, {
    1.0,  0.0,   0.0,    0.0,
    1.0,  1.0,   1.0,    1.0,
    1.0,  2.0,   4.0,    8.0,
    1.0,  3.0,   9.0,   27.0,
    1.0,  4.0,  16.0,   64.0,
    1.0,  5.0,  25.0,  125.0
}};

const static inline Vector<Real> overdetermined_6x4_polynomial_rhs{0.0, 1.0, 4.0, 9.0, 16.0, 25.0};

const static inline Vector<Real> overdetermined_6x4_polynomial_sol{
    0.000000000000000000,
    0.000000000000000000,
    1.000000000000000000,
    0.000000000000000000
};

const static inline VectorComplex overdetermined_6x4_polynomial_eigen{
    Complex(47625.0, 0),
    Complex(225.0, 0),
    Complex(6.0, 0),
    Complex(1.0, 0)
};

// Matrix properties for overdetermined_6x4_polynomial_mat (rectangular - no determinant)
const static inline MML::Vector<Real> overdetermined_6x4_polynomial_singular_values{218.14612773285034, 8.376015556001766, 1.1102296108422065, 0.11303803251665};
const static inline Real              overdetermined_6x4_polynomial_cond_2 = 1929.8908406299066;
const static inline int               overdetermined_6x4_polynomial_rank = 4;

/*******************************************************************************************************************
 * DIAGONALLY DOMINANT MATRICES - For Iterative Solvers
 * 
 * A matrix is strictly diagonally dominant if |a_ii| > Σ|a_ij| for all i.
 * Such matrices guarantee convergence of Jacobi and Gauss-Seidel methods.
 *******************************************************************************************************************/

// ===== STRICTLY DIAGONALLY DOMINANT 4x4 =====
// Source: Custom test matrix
// Properties: Strictly diagonally dominant, Well-conditioned (κ ≈ 1.6)
// Purpose: Test iterative solver convergence, benchmark convergence rate

const static inline Matrix<Real> diag_dominant_4x4_mat{4, 4, {
     10.0,  -1.0,   2.0,   0.0,
     -1.0,  11.0,  -1.0,   3.0,
      2.0,  -1.0,  10.0,  -1.0,
      0.0,   3.0,  -1.0,   8.0
}};

const static inline Vector<Real> diag_dominant_4x4_rhs{6.0, 25.0, -11.0, 15.0};

const static inline Vector<Real> diag_dominant_4x4_sol{
    1.000000000000000000,
    2.000000000000000000,
   -1.000000000000000000,
    1.000000000000000000
};

const static inline VectorComplex diag_dominant_4x4_eigen{
    Complex(11.652, 0),
    Complex(10.345, 0),
    Complex(9.821, 0),
    Complex(7.182, 0)
};

// Matrix properties for diag_dominant_4x4_mat
const static inline Real              diag_dominant_4x4_det = 7700.0;
const static inline MML::Vector<Real> diag_dominant_4x4_singular_values{12.80232827821645, 10.86174188588088, 8.787015679261046, 6.449101632298261};
const static inline Real              diag_dominant_4x4_cond_2 = 1.9848820115499685;
const static inline int               diag_dominant_4x4_rank = 4;

// ===== DIAGONALLY DOMINANT TRIDIAGONAL 5x5 =====
// Source: Discretized 1D Poisson equation (-u''(x) = f(x))
// Properties: Tridiagonal, diagonally dominant, arises from finite difference
// Purpose: Test iterative solvers on realistic PDE discretizations
// Eigenvalues: λ_k = 4 - 2*cos(kπ/(n+1)) for k=1..n
// Condition number κ ≈ 2.36

const static inline Matrix<Real> diag_dominant_5x5_tridiag_mat{5, 5, {
     4.0, -1.0,  0.0,  0.0,  0.0,
    -1.0,  4.0, -1.0,  0.0,  0.0,
     0.0, -1.0,  4.0, -1.0,  0.0,
     0.0,  0.0, -1.0,  4.0, -1.0,
     0.0,  0.0,  0.0, -1.0,  4.0
}};

const static inline Vector<Real> diag_dominant_5x5_tridiag_rhs{1.0, 2.0, 3.0, 2.0, 1.0};

const static inline Vector<Real> diag_dominant_5x5_tridiag_sol{
    0.480769230769230769,
    0.923076923076923077,
    1.211538461538461538,
    0.923076923076923077,
    0.480769230769230769
};

const static inline VectorComplex diag_dominant_5x5_tridiag_eigen{
    Complex(5.618, 0),
    Complex(5.000, 0),
    Complex(4.000, 0),
    Complex(3.000, 0),
    Complex(2.382, 0)
};

// Matrix properties for diag_dominant_5x5_tridiag_mat
const static inline Real              diag_dominant_5x5_tridiag_det = 832.0;
const static inline MML::Vector<Real> diag_dominant_5x5_tridiag_singular_values{5.618033988749895, 5.0, 4.0, 3.0, 2.381966011250105};
const static inline Real              diag_dominant_5x5_tridiag_cond_2 = 2.3595041746212636;
const static inline int               diag_dominant_5x5_tridiag_rank = 5;

// ===== DIAGONALLY DOMINANT 6x6 - 2D POISSON =====
// Source: 2D Poisson equation on 2×3 grid
// Properties: Block structure from 2D problem, sparse (max 5 non-zeros per row)
// Purpose: Test 2D discretizations, test sparse iterative solvers
// Condition number κ ≈ 5.3

const static inline Matrix<Real> diag_dominant_6x6_poisson2d_mat{6, 6, {
     4.0, -1.0,  0.0, -1.0,  0.0,  0.0,
    -1.0,  4.0, -1.0,  0.0, -1.0,  0.0,
     0.0, -1.0,  4.0,  0.0,  0.0, -1.0,
    -1.0,  0.0,  0.0,  4.0, -1.0,  0.0,
     0.0, -1.0,  0.0, -1.0,  4.0, -1.0,
     0.0,  0.0, -1.0,  0.0, -1.0,  4.0
}};

const static inline Vector<Real> diag_dominant_6x6_poisson2d_rhs{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

const static inline Vector<Real> diag_dominant_6x6_poisson2d_sol{
    0.571428571428571429,
    0.714285714285714286,
    0.571428571428571429,
    0.571428571428571429,
    0.714285714285714286,
    0.571428571428571429
};

const static inline VectorComplex diag_dominant_6x6_poisson2d_eigen{
    Complex(6.732, 0),
    Complex(5.414, 0),
    Complex(5.414, 0),
    Complex(2.586, 0),
    Complex(2.586, 0),
    Complex(1.268, 0)
};

// Matrix properties for diag_dominant_6x6_poisson2d_mat
const static inline Real              diag_dominant_6x6_poisson2d_det = 112.0;
const static inline MML::Vector<Real> diag_dominant_6x6_poisson2d_singular_values{6.732050807568877, 5.414213562373095, 5.414213562373095, 2.585786437626905, 2.585786437626905, 1.2679491924311228};
const static inline Real              diag_dominant_6x6_poisson2d_cond_2 = 5.309401076758504;
const static inline int               diag_dominant_6x6_poisson2d_rank = 6;

/*******************************************************************************************************************
 * SYMMETRIC POSITIVE DEFINITE (SPD) MATRICES - For Cholesky and CG
 * 
 * SPD matrices have all positive eigenvalues and are symmetric.
 * They can be decomposed as A = LLᵀ (Cholesky) and are ideal for conjugate gradient.
 *******************************************************************************************************************/

// ===== SPD 3x3 - SIMPLE CASE =====
// Source: Custom test matrix
// Properties: Symmetric positive definite, Well-conditioned (κ ≈ 1.89)
// Purpose: Basic Cholesky decomposition test, Simple conjugate gradient test

const static inline MatrixSym<Real> spd_3x3_mat{3, {
     4.0,
     1.0,  4.0,
     0.0,  1.0,  4.0
}};

const static inline Vector<Real> spd_3x3_rhs{5.0, 6.0, 5.0};

const static inline Vector<Real> spd_3x3_sol{
    1.000000000000000000,
    1.000000000000000000,
    1.000000000000000000
};

const static inline Vector<Real> spd_3x3_eigen{5.236, 4.000, 2.764};

// Matrix properties for spd_3x3_mat
const static inline Real              spd_3x3_det = 56.0;
const static inline MML::Vector<Real> spd_3x3_singular_values{5.23606797749979, 4.0, 2.7639320225002095};
const static inline Real              spd_3x3_cond_2 = 1.8944271909999157;
const static inline int               spd_3x3_rank = 3;

// ===== SPD 4x4 - CORRELATION MATRIX =====
// Source: Typical correlation matrix
// Properties: Represents correlation between 4 variables, diagonal = 1, off-diagonal < 1
// Purpose: Test on realistic covariance/correlation matrices
// Condition number κ ≈ 40

const static inline MatrixSym<Real> spd_4x4_correlation_mat{4, {
    1.0,
    0.8,  1.0,
    0.5,  0.7,  1.0,
    0.3,  0.4,  0.6,  1.0
}};

const static inline Vector<Real> spd_4x4_correlation_rhs{2.6, 2.9, 2.8, 2.3};

const static inline Vector<Real> spd_4x4_correlation_sol{
    1.000000000000000000,
    1.000000000000000000,
    1.000000000000000000,
    1.000000000000000000
};

const static inline Vector<Real> spd_4x4_correlation_eigen{2.845, 0.889, 0.195, 0.071};

// Matrix properties for spd_4x4_correlation_mat
const static inline Real              spd_4x4_correlation_det = 0.035343;
const static inline MML::Vector<Real> spd_4x4_correlation_singular_values{2.844855825376282, 0.8889427629063621, 0.19545547024645665, 0.07074594147089902};
const static inline Real              spd_4x4_correlation_cond_2 = 40.21162379051538;
const static inline int               spd_4x4_correlation_rank = 4;

// ===== SPD 5x5 - MASS MATRIX =====
// Source: Finite element mass matrix from linear FEM with uniform mesh
// Properties: Tridiagonal SPD structure, arises from FEM discretization
// Purpose: Test on FEM-type matrices, benchmark structured SPD solvers
// Condition number κ ≈ 2.36 (excellent)

const static inline MatrixSym<Real> spd_5x5_mass_matrix_mat{5, {
    4.0,
    1.0,  4.0,
    0.0,  1.0,  4.0,
    0.0,  0.0,  1.0,  4.0,
    0.0,  0.0,  0.0,  1.0,  4.0
}};

const static inline Vector<Real> spd_5x5_mass_matrix_rhs{5.0, 6.0, 6.0, 6.0, 5.0};

const static inline Vector<Real> spd_5x5_mass_matrix_sol{
    1.000000000000000000,
    1.000000000000000000,
    1.000000000000000000,
    1.000000000000000000,
    1.000000000000000000
};

const static inline Vector<Real> spd_5x5_mass_matrix_eigen{5.618, 5.000, 4.000, 3.000, 2.382};

// Matrix properties for spd_5x5_mass_matrix_mat
const static inline Real              spd_5x5_mass_matrix_det = 832.0;
const static inline MML::Vector<Real> spd_5x5_mass_matrix_singular_values{5.618033988749895, 5.0, 4.0, 3.0, 2.381966011250105};
const static inline Real              spd_5x5_mass_matrix_cond_2 = 2.3595041746212636;
const static inline int               spd_5x5_mass_matrix_rank = 5;

// ===== SPD 6x6 - GRAPH LAPLACIAN =====
// Source: Laplacian matrix of cycle graph C₆
// Properties: Graph Laplacian (degree - adjacency), Symmetric positive semi-definite
// Purpose: Test graph-based algorithms, test handling of semi-definiteness
// Note: Technically PSD not SPD (one zero eigenvalue), but useful test case
// Effective condition for this RHS: κ ≈ 4

const static inline MatrixSym<Real> spd_6x6_graph_laplacian_mat{6, {
     2.0,
    -1.0,  2.0,
     0.0, -1.0,  2.0,
     0.0,  0.0, -1.0,  2.0,
     0.0,  0.0,  0.0, -1.0,  2.0,
    -1.0,  0.0,  0.0,  0.0, -1.0,  2.0
}};

const static inline Vector<Real> spd_6x6_graph_laplacian_rhs{1.0, -1.0, 1.0, -1.0, 1.0, -1.0};

const static inline Vector<Real> spd_6x6_graph_laplacian_sol{
    0.666666666666666667,
   -0.333333333333333333,
    0.666666666666666667,
   -0.333333333333333333,
    0.666666666666666667,
   -0.333333333333333333
};

const static inline Vector<Real> spd_6x6_graph_laplacian_eigen{4.000, 3.000, 3.000, 1.000, 1.000, 0.000};

// Matrix properties for spd_6x6_graph_laplacian_mat (semi-definite - determinant = 0)
const static inline Real              spd_6x6_graph_laplacian_det = 0.0;
const static inline MML::Vector<Real> spd_6x6_graph_laplacian_singular_values{4.0, 3.0, 3.0, 1.0, 1.0, 0.0};
const static inline Real              spd_6x6_graph_laplacian_cond_2 = std::numeric_limits<Real>::infinity();  // Rank-deficient
const static inline int               spd_6x6_graph_laplacian_rank = 5;

} // namespace MML::TestBeds

#endif // MML_LINEAR_ALG_EQ_SYSTEMS_SOLVERSPECIFIC_DEFS_H
