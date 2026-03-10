#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_TRIDIAG_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_TRIDIAG_DEFS_H

/*******************************************************************************************************************
 * TRIDIAGONAL TEST MATRICES FOR LINEAR EQUATION SOLVERS
 * 
 * This file contains tridiagonal matrices specifically designed for:
 * 
 * 1. CUBIC SPLINE INTERPOLATION
 *    - Natural cubic splines produce tridiagonal systems
 *    - Clamped splines also produce tridiagonal systems
 * 
 * 2. FINITE DIFFERENCE METHODS
 *    - 1D Poisson equation: -u''(x) = f(x) discretization
 *    - Heat equation implicit schemes
 *    - Wave equation implicit schemes
 * 
 * 3. THOMAS ALGORITHM TESTING
 *    - Direct O(n) solver for tridiagonal systems
 *    - Numerical stability verification
 * 
 * 4. ITERATIVE SOLVER BENCHMARKS
 *    - Tridiagonal matrices are well-suited for SOR, Jacobi
 *    - Known spectral properties for convergence analysis
 * 
 * Matrix Types:
 * - General tridiagonal (non-symmetric)
 * - Symmetric positive definite (SPD) tridiagonal
 * - Toeplitz tridiagonal (constant diagonals)
 * 
 * Each test system includes:
 * - Coefficient matrix A (stored as full matrix for compatibility)
 * - Right-hand side vector b
 * - Exact solution x
 * - Eigenvalues (real for SPD, complex for general)
 * - Documentation: source, properties, condition number, intended use
 * 
 * SIOF-SAFE PATTERN: All test data uses inline functions with static local variables
 * to avoid Static Initialization Order Fiasco issues across translation units.
 * 
 *******************************************************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#endif

namespace MML::TestBeds
{

/*******************************************************************************************************************
 * SYMMETRIC POSITIVE DEFINITE TRIDIAGONAL MATRICES
 * 
 * These arise from finite difference discretization of -u''(x) = f(x)
 * Eigenvalues: λ_k = 2 - 2*cos(k*π/(n+1)) = 4*sin²(k*π/(2(n+1))) for k=1..n
 * Condition number: κ ≈ (n+1)²/π² for large n
 *******************************************************************************************************************/

// ===== SPD TRIDIAGONAL 3x3 =====
// Source: 1D Poisson equation with Dirichlet BCs, h=0.25
// Properties: Symmetric positive definite, well-conditioned
// Eigenvalues: λ_k = 2 - 2*cos(k*π/4) for k=1,2,3
// Condition number κ ≈ 5.83

inline const MML::Matrix<Real>& tridiag_spd_3x3() {
    static const MML::Matrix<Real> m{3, 3, {
         2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_spd_3x3_rhs() {
    static const MML::Vector<Real> v{1.0, 2.0, 1.0};
    return v;
}

// Solution computed analytically: x = A^(-1) * b
inline const MML::Vector<Real>& tridiag_spd_3x3_sol() {
    static const MML::Vector<Real> v{
        2.0,
        3.0,
        2.0
    };
    return v;
}

// Eigenvalues: 2 - 2*cos(kπ/4) for k=1,2,3
inline const VectorComplex& tridiag_spd_3x3_eigen() {
    static const VectorComplex v{
        Complex(2.0 - 2.0*0.707106781186548, 0),   // 2 - √2 ≈ 0.5858
        Complex(2.0, 0),                            // 2
        Complex(2.0 + 2.0*0.707106781186548, 0)    // 2 + √2 ≈ 3.4142
    };
    return v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ===== SPD TRIDIAGONAL 5x5 =====
// Source: 1D Poisson equation with Dirichlet BCs, n=5 interior points
// Properties: Symmetric positive definite, Toeplitz structure
// Condition number κ ≈ 14.93

inline const MML::Matrix<Real>& tridiag_spd_5x5() {
    static const MML::Matrix<Real> m{5, 5, {
         2.0, -1.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0, -1.0,  2.0
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_spd_5x5_rhs() {
    static const MML::Vector<Real> v{1.0, 0.0, 0.0, 0.0, 1.0};
    return v;
}

inline const MML::Vector<Real>& tridiag_spd_5x5_sol() {
    static const MML::Vector<Real> v{
        0.833333333333333333,
        0.666666666666666667,
        0.5,
        0.666666666666666667,
        0.833333333333333333
    };
    return v;
}

// Eigenvalues: 2 - 2*cos(kπ/6) for k=1..5
inline const VectorComplex& tridiag_spd_5x5_eigen() {
    static const VectorComplex v{
        Complex(0.267949192431123, 0),   // 2 - 2*cos(π/6) = 2 - √3
        Complex(1.0, 0),                  // 2 - 2*cos(2π/6) = 2 - 1 = 1
        Complex(2.0, 0),                  // 2 - 2*cos(3π/6) = 2 - 0 = 2
        Complex(3.0, 0),                  // 2 - 2*cos(4π/6) = 2 + 1 = 3
        Complex(3.732050807568877, 0)    // 2 - 2*cos(5π/6) = 2 + √3
    };
    return v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ===== SPD TRIDIAGONAL 10x10 =====
// Source: 1D Poisson discretization
// Properties: SPD, Toeplitz, well-conditioned for moderate size
// Condition number κ ≈ 55.9

inline const MML::Matrix<Real>& tridiag_spd_10x10() {
    static const MML::Matrix<Real> m{10, 10, {
         2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0
    }};
    return m;
}

// RHS corresponding to f(x) = sin(πx) on [0,1] with h=1/11
inline const MML::Vector<Real>& tridiag_spd_10x10_rhs() {
    static const MML::Vector<Real> v{
        0.281732556841429696,   // sin(π/11) * h²
        0.540640817455597582,   // sin(2π/11) * h²
        0.755749574354258273,   // sin(3π/11) * h²
        0.909631995354518436,   // sin(4π/11) * h²
        0.989821441880932732,   // sin(5π/11) * h²
        0.989821441880932732,   // sin(6π/11) * h²
        0.909631995354518436,   // sin(7π/11) * h²
        0.755749574354258273,   // sin(8π/11) * h²
        0.540640817455597582,   // sin(9π/11) * h²
        0.281732556841429696    // sin(10π/11) * h²
    };
    return v;
}

// Solution: u(x) = sin(πx)/π² (scaled by h²)
inline const MML::Vector<Real>& tridiag_spd_10x10_sol() {
    static const MML::Vector<Real> v{
        0.028559933214452666,
        0.054823100567116405,
        0.076616599662254800,
        0.092219329846383850,
        0.100352567143534780,
        0.100352567143534780,
        0.092219329846383850,
        0.076616599662254800,
        0.054823100567116405,
        0.028559933214452666
    };
    return v;
}

// Eigenvalues: 2 - 2*cos(kπ/11) for k=1..10
inline const VectorComplex& tridiag_spd_10x10_eigen() {
    static const VectorComplex v{
        Complex(0.081014258856924, 0),
        Complex(0.317493087989917, 0),
        Complex(0.690279942783750, 0),
        Complex(1.169769746792498, 0),
        Complex(1.715616038353807, 0),
        Complex(2.284383961646193, 0),
        Complex(2.830230253207502, 0),
        Complex(3.309720057216250, 0),
        Complex(3.682506912010083, 0),
        Complex(3.918985741143076, 0)
    };
    return v;
}

/*******************************************************************************************************************
 * GENERAL (NON-SYMMETRIC) TRIDIAGONAL MATRICES
 * 
 * These test non-symmetric tridiagonal solvers and stability
 *******************************************************************************************************************/

// ===== GENERAL TRIDIAGONAL 4x4 =====
// Source: Custom test matrix with varying off-diagonals
// Properties: Non-symmetric, diagonally dominant
// Purpose: Test Thomas algorithm on non-symmetric case

inline const MML::Matrix<Real>& tridiag_general_4x4() {
    static const MML::Matrix<Real> m{4, 4, {
         4.0, -1.0,  0.0,  0.0,
        -2.0,  5.0, -1.0,  0.0,
         0.0, -2.0,  6.0, -1.0,
         0.0,  0.0, -2.0,  7.0
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_general_4x4_rhs() {
    static const MML::Vector<Real> v{3.0, 2.0, 3.0, 5.0};
    return v;
}

inline const MML::Vector<Real>& tridiag_general_4x4_sol() {
    static const MML::Vector<Real> v{
        0.895833333333333333,
        0.583333333333333333,
        0.666666666666666667,
        0.904761904761904762
    };
    return v;
}

inline const VectorComplex& tridiag_general_4x4_eigen() {
    static const VectorComplex v{
        Complex(3.267949192431123, 0),
        Complex(4.585786437626905, 0),
        Complex(6.414213562373095, 0),
        Complex(7.732050807568877, 0)
    };
    return v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ===== GENERAL TRIDIAGONAL 6x6 =====
// Source: Convection-diffusion discretization (upwind scheme)
// Properties: Non-symmetric, diagonally dominant, stable for Thomas algorithm
// Purpose: Test solvers on physically-motivated non-symmetric system

inline const MML::Matrix<Real>& tridiag_general_6x6() {
    static const MML::Matrix<Real> m{6, 6, {
         3.0, -1.0,  0.0,  0.0,  0.0,  0.0,
        -1.5,  3.0, -1.0,  0.0,  0.0,  0.0,
         0.0, -1.5,  3.0, -1.0,  0.0,  0.0,
         0.0,  0.0, -1.5,  3.0, -1.0,  0.0,
         0.0,  0.0,  0.0, -1.5,  3.0, -1.0,
         0.0,  0.0,  0.0,  0.0, -1.5,  3.0
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_general_6x6_rhs() {
    static const MML::Vector<Real> v{2.0, 1.5, 1.5, 1.5, 1.5, 1.5};
    return v;
}

inline const MML::Vector<Real>& tridiag_general_6x6_sol() {
    static const MML::Vector<Real> v{
        0.932515337423313,
        0.797546012269939,
        0.768404907975460,
        0.777914110429448,
        0.805521472392638,
        0.902760736196319
    };
    return v;
}

inline const VectorComplex& tridiag_general_6x6_eigen() {
    static const VectorComplex v{
        Complex(1.340142287498996, 0),
        Complex(1.933050714633240, 0),
        Complex(2.709430584957905, 0),
        Complex(3.290569415042095, 0),
        Complex(4.066949285366760, 0),
        Complex(4.659857712501004, 0)
    };
    return v;
}

/*******************************************************************************************************************
 * CUBIC SPLINE TRIDIAGONAL MATRICES
 * 
 * These are specifically designed for cubic spline interpolation testing
 * Natural spline: second derivative = 0 at endpoints
 * The system has the form: [1, 4, 1] on tridiagonal
 *******************************************************************************************************************/

// ===== NATURAL CUBIC SPLINE 5x5 =====
// Source: Natural cubic spline second derivative system
// Properties: Symmetric positive definite, [1, 4, 1] structure
// Purpose: Test spline interpolation solver
// Note: Boundary rows modified for natural spline conditions

inline const MML::Matrix<Real>& tridiag_spline_5x5() {
    static const MML::Matrix<Real> m{5, 5, {
         1.0,  0.0,  0.0,  0.0,  0.0,
         1.0,  4.0,  1.0,  0.0,  0.0,
         0.0,  1.0,  4.0,  1.0,  0.0,
         0.0,  0.0,  1.0,  4.0,  1.0,
         0.0,  0.0,  0.0,  0.0,  1.0
    }};
    return m;
}

// RHS for spline fitting y = x² at x = {0, 1, 2, 3, 4}
// Natural spline: M_0 = M_4 = 0, interior from 6(y_{i+1} - 2y_i + y_{i-1})/h²
inline const MML::Vector<Real>& tridiag_spline_5x5_rhs() {
    static const MML::Vector<Real> v{0.0, 12.0, 12.0, 12.0, 0.0};
    return v;
}

// Second derivatives (M values) for y = x²
inline const MML::Vector<Real>& tridiag_spline_5x5_sol() {
    static const MML::Vector<Real> v{
        0.0,
        2.571428571428571,
        2.0,
        2.571428571428571,
        0.0
    };
    return v;
}

inline const VectorComplex& tridiag_spline_5x5_eigen() {
    static const VectorComplex v{
        Complex(1.0, 0),
        Complex(1.0, 0),
        Complex(2.585786437626905, 0),
        Complex(4.0, 0),
        Complex(5.414213562373095, 0)
    };
    return v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ===== INTERIOR SPLINE SYSTEM 4x4 =====
// Source: Pure interior spline system (4x4 for 6 data points)
// Properties: SPD, strictly diagonally dominant
// Purpose: Test standard [1,4,1] spline solver

inline const MML::Matrix<Real>& tridiag_spline_interior_4x4() {
    static const MML::Matrix<Real> m{4, 4, {
         4.0,  1.0,  0.0,  0.0,
         1.0,  4.0,  1.0,  0.0,
         0.0,  1.0,  4.0,  1.0,
         0.0,  0.0,  1.0,  4.0
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_spline_interior_4x4_rhs() {
    static const MML::Vector<Real> v{6.0, 0.0, 0.0, 6.0};
    return v;
}

inline const MML::Vector<Real>& tridiag_spline_interior_4x4_sol() {
    static const MML::Vector<Real> v{
        1.565217391304348,
       -0.260869565217391,
       -0.260869565217391,
        1.565217391304348
    };
    return v;
}

inline const VectorComplex& tridiag_spline_interior_4x4_eigen() {
    static const VectorComplex v{
        Complex(2.381966011250105, 0),
        Complex(3.381966011250105, 0),
        Complex(4.618033988749895, 0),
        Complex(5.618033988749895, 0)
    };
    return v;
}

/*******************************************************************************************************************
 * VARIABLE COEFFICIENT TRIDIAGONAL MATRICES
 * 
 * These test solvers on non-constant coefficient problems
 *******************************************************************************************************************/

// ===== VARIABLE COEFFICIENT 5x5 =====
// Source: -d/dx(a(x)*du/dx) = f(x) with a(x) = 1+x
// Properties: SPD but non-Toeplitz
// Purpose: Test on realistic variable-coefficient BVP

inline const MML::Matrix<Real>& tridiag_variable_5x5() {
    static const MML::Matrix<Real> m{5, 5, {
         2.5, -1.5,  0.0,  0.0,  0.0,
        -1.5,  3.5, -2.0,  0.0,  0.0,
         0.0, -2.0,  4.5, -2.5,  0.0,
         0.0,  0.0, -2.5,  5.5, -3.0,
         0.0,  0.0,  0.0, -3.0,  6.5
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_variable_5x5_rhs() {
    static const MML::Vector<Real> v{1.0, 1.0, 1.0, 1.0, 1.0};
    return v;
}

inline const MML::Vector<Real>& tridiag_variable_5x5_sol() {
    static const MML::Vector<Real> v{
        0.759259259259259259,
        0.518518518518518519,
        0.388888888888888889,
        0.333333333333333333,
        0.307692307692307692
    };
    return v;
}

inline const VectorComplex& tridiag_variable_5x5_eigen() {
    static const VectorComplex v{
        Complex(0.744429004028256, 0),
        Complex(2.069903290878389, 0),
        Complex(3.884762024618692, 0),
        Complex(6.096135665096665, 0),
        Complex(9.704770015377998, 0)
    };
    return v;
}

/*******************************************************************************************************************
 * NEARLY SINGULAR TRIDIAGONAL
 * 
 * For testing solver robustness near singularity
 *******************************************************************************************************************/

// ===== NEARLY SINGULAR TRIDIAGONAL 4x4 =====
// Source: Modified Poisson with small diagonal perturbation
// Properties: SPD but poorly conditioned
// Condition number κ ≈ 1000

inline const MML::Matrix<Real>& tridiag_nearly_singular_4x4() {
    static const MML::Matrix<Real> m{4, 4, {
         0.002, -0.001,  0.0,    0.0,
        -0.001,  0.002, -0.001,  0.0,
         0.0,   -0.001,  0.002, -0.001,
         0.0,    0.0,   -0.001,  0.002
    }};
    return m;
}

inline const MML::Vector<Real>& tridiag_nearly_singular_4x4_rhs() {
    static const MML::Vector<Real> v{0.001, 0.0, 0.0, 0.001};
    return v;
}

inline const MML::Vector<Real>& tridiag_nearly_singular_4x4_sol() {
    static const MML::Vector<Real> v{
        0.8,
        0.6,
        0.6,
        0.8
    };
    return v;
}

inline const VectorComplex& tridiag_nearly_singular_4x4_eigen() {
    static const VectorComplex v{
        Complex(0.000381966011250105, 0),
        Complex(0.001, 0),
        Complex(0.002, 0),
        Complex(0.003618033988749895, 0)
    };
    return v;
}

/*******************************************************************************************************************
 * LARGER TRIDIAGONAL FOR PERFORMANCE TESTING
 *******************************************************************************************************************/

// ===== SPD TRIDIAGONAL 20x20 =====
// Source: Standard 1D Poisson discretization
// Properties: SPD, Toeplitz, moderate condition number
// Condition number κ ≈ 211.2

inline const MML::Matrix<Real>& tridiag_spd_20x20() {
    static const MML::Matrix<Real> m{20, 20, {
         2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  2.0
    }};
    return m;
}

// Unit vector RHS for simple testing
inline const MML::Vector<Real>& tridiag_spd_20x20_rhs() {
    static const MML::Vector<Real> v{
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
    };
    return v;
}

inline const MML::Vector<Real>& tridiag_spd_20x20_sol() {
    static const MML::Vector<Real> v{
        0.952380952380952381,
        0.904761904761904762,
        0.857142857142857143,
        0.809523809523809524,
        0.761904761904761905,
        0.714285714285714286,
        0.666666666666666667,
        0.619047619047619048,
        0.571428571428571429,
        0.523809523809523810,
        0.523809523809523810,
        0.571428571428571429,
        0.619047619047619048,
        0.666666666666666667,
        0.714285714285714286,
        0.761904761904761905,
        0.809523809523809524,
        0.857142857142857143,
        0.904761904761904762,
        0.952380952380952381
    };
    return v;
}

// Eigenvalues: 2 - 2*cos(kπ/21) for k=1..20
inline const VectorComplex& tridiag_spd_20x20_eigen() {
    static const VectorComplex v{
        Complex(0.022193688576881, 0),
        Complex(0.088390450088978, 0),
        Complex(0.197477208423344, 0),
        Complex(0.347296355333861, 0),
        Complex(0.534565752092256, 0),
        Complex(0.754915028125263, 0),
        Complex(1.002923048454133, 0),
        Complex(1.273096188178695, 0),
        Complex(1.559816222183262, 0),
        Complex(1.857372825564638, 0),
        Complex(2.159975698706451, 0),
        Complex(2.461780721227131, 0),
        Complex(2.756918851094295, 0),
        Complex(3.039543269700609, 0),
        Complex(3.303863589714455, 0),
        Complex(3.544187046831166, 0),
        Complex(3.754949771723356, 0),
        Complex(3.930757232315044, 0),
        Complex(4.066419952737040, 0),
        Complex(4.156986339905760, 0)
    };
    return v;
}

}  // namespace MML::TestBeds

#endif  // MML_LINEAR_ALG_EQ_SYSTEMS_TRIDIAG_DEFS_H