///////////////////////////////////////////////////////////////////////////////////////////
// ODESolvers.h - Umbrella header for all ODE solver functionality
///////////////////////////////////////////////////////////////////////////////////////////
//
// This header provides a convenient way to include all ODE (Ordinary Differential Equation)
// solver functionality in the MML library.
//
// Components:
// - ODESystemSteppers.h      - Base stepper classes and interfaces
// - ODESystemStepCalculators.h - Runge-Kutta and other step calculators
// - ODEFixedStepIntegrators.h  - Fixed step size integrators
// - ODEAdaptiveIntegrator.h    - Adaptive step size integrators
// - ODEStiffSolvers.h          - Solvers for stiff ODEs (implicit methods)
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

#include "mml/algorithms/ODESolvers/ODESystemSteppers.h"
#include "mml/algorithms/ODESolvers/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "mml/algorithms/ODESolvers/ODEStiffSolvers.h"
#include "mml/algorithms/ODESolvers/BVPShootingMethod.h"

#endif // MML_ODE_SOLVERS_H
