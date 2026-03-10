///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgEqSolvers.h                                                   ///
///  Description: Umbrella header for all linear algebra equation solvers             ///
///               Includes Direct (LU, Cholesky, Gauss), QR, SVD, and iterative       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
///
/// This is an umbrella header that includes all linear algebra equation solvers.
/// For finer-grained includes, use the individual headers:
///
///   - LinAlgDirect.h           : GaussJordanSolver, LUSolver, LUSolverInPlace,
///                                BandDiagonalSolver, CholeskySolver
///   - LinAlgQR.h               : QRSolver (Householder reflections)
///   - LinAlgSVD.h              : SVDecompositionSolver (pseudoinverse, least squares)
///   - LinAlgEqSolvers_iterative.h : Jacobi, Gauss-Seidel, SOR, Conjugate Gradient
///
/// Thread Safety:
///   All linear solvers are REENTRANT - safe to call from multiple threads with
///   different matrices. Solvers are stateless pure functions.
///   See docs/THREADING.md for details.
///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_LINEAR_ALG_EQ_SOLVERS_H
#define MML_LINEAR_ALG_EQ_SOLVERS_H

// Direct solvers: Gauss-Jordan, LU decomposition, Band diagonal, Cholesky
#include "core/LinAlgEqSolvers/LinAlgDirect.h"

// QR decomposition solver using Householder reflections
#include "core/LinAlgEqSolvers/LinAlgQR.h"

// Singular Value Decomposition solver
#include "core/LinAlgEqSolvers/LinAlgSVD.h"

// Iterative solvers: Jacobi, Gauss-Seidel, SOR, Conjugate Gradient
#include "core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h"

#endif // MML_LINEAR_ALG_EQ_SOLVERS_H
