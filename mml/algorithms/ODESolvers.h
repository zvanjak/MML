///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESolvers.h                                                        ///
///  Description: Aggregate header for ODE solver functionality                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file ODESolvers.h
/// @brief Aggregate header for ordinary differential equation solvers.
/// This header includes all ODE solver functionality:
/// - Base stepper interfaces and Runge-Kutta step calculators
/// - Fixed-step, adaptive-step, and stiff ODE integrators
/// - Boundary-value shooting helpers
///
/// Main types included here:
/// - ODESystemFixedStepSolver - fixed-step integration using explicit step calculators
/// - EulerStep_Calculator, Midpoint_StepCalculator, and RungeKutta4_StepCalculator
/// - DormandPrince5_Stepper, CashKarp_Stepper, DormandPrince8_Stepper, and BulirschStoer_Stepper
/// - ODEAdaptiveIntegrator - adaptive integration with step-size control and solution statistics
/// - Rosenbrock23Solver - stiff ODE integration
/// - BVPShootingSolver and BVPShootingSolverND - shooting-method solvers for boundary-value problems
///
/// All ODE solvers are reentrant when separate solver instances are used per thread.
/// For more focused includes, use the individual headers under ODESolvers/.

#ifndef MML_ODE_SOLVERS_H
#define MML_ODE_SOLVERS_H

#include <mml/algorithms/ODESolvers/ODESteppers.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>
#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#include <mml/algorithms/ODESolvers/BVPShootingMethod.h>

#endif // MML_ODE_SOLVERS_H
