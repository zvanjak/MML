///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        EigenSystemSolvers.h                                                ///
///  Description: Aggregate header for eigensystem solver implementations             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file EigenSystemSolvers.h
/// @brief Aggregate header for eigenvalue and eigenvector solvers.
/// This header includes symmetric, Hermitian, QR-based, complex, and general
/// eigensystem solver implementations.
///
/// Main solver types included here:
/// - SymmMatEigenSolverJacobi - Jacobi solver for real symmetric matrices
/// - HermitianMatEigenSolverJacobi - Jacobi solver for complex Hermitian matrices
/// - SymmMatEigenSolverQR - QR-based solver for real symmetric matrices
/// - ComplexEigenSolver - nonsymmetric complex eigensystem solver
/// - EigenSolver - general real eigensystem front end
///
/// For more focused includes, use the individual headers under Eigen/.

#ifndef MML_EIGENSYSTEM_SOLVERS_H
#define MML_EIGENSYSTEM_SOLVERS_H

#include <mml/algorithms/Eigen/SymmMatEigenSolverJacobi.h>
#include <mml/algorithms/Eigen/HermitianMatEigenSolverJacobi.h>
#include <mml/algorithms/Eigen/SymmMatEigenSolverQR.h>
#include <mml/algorithms/Eigen/ComplexEigenSolver.h>
#include <mml/algorithms/Eigen/EigenSolver.h>

#endif // MML_EIGENSYSTEM_SOLVERS_H
