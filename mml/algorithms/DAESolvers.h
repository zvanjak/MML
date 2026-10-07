///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAESolvers.h                                                        ///
///  Description: Aggregate header for DAE solvers and support types                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file DAESolvers.h
/// @brief Aggregate header for DAE solver functionality.
/// This header includes all differential-algebraic equation solver components:
/// - DAE interfaces and system data structures
/// - Solver configuration, result, and initial-condition helpers
/// - Numerical Jacobian adapters
/// - Backward Euler, BDF2, BDF4, RODAS, Radau IIA, adaptive, and event-aware solvers
///
/// Main types and functions included here:
/// - IODESystemDAE and IODESystemDAEWithEvents - DAE system interfaces
/// - DAESystem - reusable DAE system container
/// - DAESolverConfig, DAESolverResult, DAEFailureReason, and DAENewtonResult
/// - ComputeConsistentIC and VerifyConsistentIC - consistent initial-condition helpers
/// - DAESystemNumericalJacobian - Jacobian adapter for systems without analytical Jacobians
/// - SolveDAEBackwardEuler, SolveDAEBDF2, SolveDAEBDF4, SolveDAERODAS, and SolveDAERadauIIA
/// - DAEEventConfig, DAEEventInfo, and DAEEventResult - event detection support
///
/// For more focused includes, use the individual headers under DAESolvers/.

#if !defined MML_DAE_SOLVERS_H
#define MML_DAE_SOLVERS_H

#include <mml/interfaces/IODESystemDAE.h>
#include <mml/interfaces/IODESystemDAEWithEvents.h>
#include <mml/base/DAESystem.h>
#include <mml/algorithms/DAESolvers/DAESolverBase.h>
#include <mml/algorithms/DAESolvers/DAENumericalJacobian.h>
#include <mml/algorithms/DAESolvers/DAEBackwardEuler.h>
#include <mml/algorithms/DAESolvers/DAEBDF2.h>
#include <mml/algorithms/DAESolvers/DAEBDF4.h>
#include <mml/algorithms/DAESolvers/DAERODAS.h>
#include <mml/algorithms/DAESolvers/DAERadauIIA.h>
#include <mml/algorithms/DAESolvers/DAEAdaptive.h>
#include <mml/algorithms/DAESolvers/DAEEventDetection.h>

#endif // MML_DAE_SOLVERS_H
