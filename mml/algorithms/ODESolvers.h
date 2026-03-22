///////////////////////////////////////////////////////////////////////////////////////////
// ODESolvers.h - Umbrella header for all ODE solver functionality
///////////////////////////////////////////////////////////////////////////////////////////
//
// This header provides a convenient way to include all ODE (Ordinary Differential Equation)
// solver functionality in the MML library.
//
// Components:
// - ODESteppers.h      - Base stepper classes and interfaces
// - ODEStepCalculators.h - Runge-Kutta and other step calculators
// - ODESolverFixedStep.h  - Fixed step size integrators
// - ODESolverAdaptive.h    - Adaptive step size integrators
// - ODESolverStiff.h          - Solvers for stiff ODEs (implicit methods)
// - BVPShootingMethod.h        - Boundary value problem shooting method
//
// Thread Safety:
//   All ODE solvers are REENTRANT - safe to call from multiple threads with different
//   solver instances. Create one solver per thread for parallel ODE solving.
//   See docs/THREADING.md for parallel usage patterns.
//
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_ODE_SOLVERS_H
#define MML_ODE_SOLVERS_H

#include "mml/algorithms/ODESolvers/ODESteppers.h"
#include "mml/algorithms/ODESolvers/ODEStepCalculators.h"
#include "mml/algorithms/ODESolvers/ODESolverFixedStep.h"
#include "mml/algorithms/ODESolvers/ODESolverAdaptive.h"
#include "mml/algorithms/ODESolvers/ODESolverStiff.h"
#include "mml/algorithms/ODESolvers/BVPShootingMethod.h"

#endif // MML_ODE_SOLVERS_H
