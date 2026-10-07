///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFinding.h                                                       ///
///  Description: Aggregate header for root-finding algorithms                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file RootFinding.h
/// @brief Aggregate header for scalar and system root-finding algorithms.
/// This header includes all root-finding functionality:
/// - Root solver configuration and result types
/// - Root isolation and bracketing utilities
/// - Bisection, false-position, Newton, secant, Ridders, and Brent methods
/// - All-root, polynomial, complex, and nonlinear-system solvers
///
/// Main types and functions included here:
/// - RootFindingConfig, RootFindingResult, and RootFindingResultFinalizer
/// - RootIsolationConfig, RootInterval, RootCandidate, and RootCandidateType
/// - BracketRoot and FindRootBrackets - bracketing utilities for real functions
/// - FindRootBisection, FindRootFalsePosition, FindRootNewton, FindRootSecant, FindRootRidders, and FindRootBrent
/// - FindAllRealRootsConfig, RealRootResult, and FindAllRealRootsResult
/// - NonlinearSystemConfig, NonlinearSystemIteration, and BasicNonlinearSystemResult
/// - ComplexRootFindingResult, FindRootNewtonComplex, FindRootMuller, and polynomial root helpers
///
/// For more focused includes, use the individual headers under RootFinding/.

#if !defined MML_ROOTFINDING_H
#define MML_ROOTFINDING_H

#include <mml/algorithms/RootFinding/RootFindingBase.h>
#include <mml/algorithms/RootFinding/RootIsolation.h>
#include <mml/algorithms/RootFinding/RootFindingBracketing.h>
#include <mml/algorithms/RootFinding/RootFindingMethods.h>
#include <mml/algorithms/RootFinding/RootFindingAllRoots.h>
#include <mml/algorithms/RootFinding/NonlinearSystemSolvers.h>
#include <mml/algorithms/RootFinding/RootFindingPolynoms.h>
#include <mml/algorithms/RootFinding/RootFindingComplex.h>

#endif // MML_ROOTFINDING_H
