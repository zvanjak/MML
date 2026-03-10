///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAESolvers.h                                                        ///
///  Description: Umbrella header for DAE (Differential-Algebraic Equation) solvers  ///
///                                                                                   ///
///  Includes:                                                                        ///
///  - DAE interfaces (IODESystemDAE, IODESystemDAEWithJacobian)                     ///
///  - DAE data structures (DAESolution, DAESystem)                                   ///
///  - DAE Index-1 solvers:                                                           ///
///      - SolveDAEBackwardEuler()  - 1st order, A-stable                             ///
///      - SolveDAEBDF2()           - 2nd order, A-stable                             ///
///      - SolveDAEBDF4()           - 4th order, A(α)-stable                          ///
///      - SolveDAERODAS()          - Rosenbrock, L-stable, no Newton                 ///
///      - SolveDAERadauIIA()       - Order 5, L-stable, gold standard IRK            ///
///                                                                                   ///
///  Common types:                                                                    ///
///  - DAESolverConfig   - Configuration for all solvers                              ///
///  - DAESolverResult   - Result structure with solution and diagnostics             ///
///  - ComputeConsistentIC()  - Compute consistent initial conditions                 ///
///  - VerifyConsistentIC()   - Verify initial conditions satisfy constraints         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_SOLVERS_H
#define MML_DAE_SOLVERS_H

// DAE interfaces
#include "mml/interfaces/IODESystemDAE.h"

// DAE data structures
#include "mml/base/DAESystem.h"

// DAE solver base types (DAESolverConfig, DAESolverResult, IC functions)
#include "mml/algorithms/DAESolvers/DAESolverBase.h"

// Individual DAE solver implementations
#include "mml/algorithms/DAESolvers/DAEBackwardEuler.h"
#include "mml/algorithms/DAESolvers/DAEBDF2.h"
#include "mml/algorithms/DAESolvers/DAEBDF4.h"
#include "mml/algorithms/DAESolvers/DAERODAS.h"
#include "mml/algorithms/DAESolvers/DAERadauIIA.h"

#endif // MML_DAE_SOLVERS_H
