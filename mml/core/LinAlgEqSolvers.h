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

/// @file LinAlgEqSolvers.h
/// @brief Aggregate header for linear algebra equation solvers.
/// This header includes direct solvers, QR solvers, real and complex SVD solvers,
/// and iterative solvers such as Jacobi, Gauss-Seidel, SOR, and conjugate gradient.
///
/// Main types included here:
/// - GaussJordanSolver - direct dense solve and matrix inverse helper
/// - LUSolver and LUSolverInPlace - LU factorization based direct solvers
/// - BandDiagonalSolver - direct solver for banded systems
/// - CholeskySolver - symmetric positive-definite direct solver
/// - QRSolver - Householder QR solver
/// - SVDecompositionSolver and ComplexSVDecompositionSolver - SVD-based least-squares and pseudoinverse solvers
/// - JacobiSolver, GaussSeidelSolver, and SORSolver - stationary iterative solvers
///
/// Linear solvers are reentrant when separate matrices or solver instances are used.
/// For more focused includes, use the individual headers under LinAlgEqSolvers/.

#if !defined  MML_LINEAR_ALG_EQ_SOLVERS_H
#define MML_LINEAR_ALG_EQ_SOLVERS_H

#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/core/LinAlgEqSolvers/LinAlgQR.h>
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>
#include <mml/core/LinAlgEqSolvers/LinAlgComplexSVD.h>
#include <mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h>

#endif // MML_LINEAR_ALG_EQ_SOLVERS_H
