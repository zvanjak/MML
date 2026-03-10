///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme07_ode_solvers.cpp                                            ///
///  Description: README Differential Equations section demo                          ///
///               Demonstrates ODE system solver with RK4                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"
#include "algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "algorithms/ODESolvers/ODESystemStepCalculators.h"
#endif

using namespace MML;

// Define ODE system: Simple harmonic oscillator
// dx/dt = v, dv/dt = -ω²x, with ω = 2.0
static void harmonicOscillator(Real t, const Vector<Real>& y, Vector<Real>& dydt)
{
    const Real omega = 2.0;
    dydt[0] = y[1];                    // dx/dt = v
    dydt[1] = -omega*omega * y[0];     // dv/dt = -ω²x
}

void Readme_ODESolvers()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            README - Differential Equations                    ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    const Real omega = 2.0;
    ODESystem system(2, harmonicOscillator);

    Vector<Real> initial_cond{1.0, 0.0};  // x(0) = 1, v(0) = 0

    // Solve with RK4 stepper
    RungeKutta4_StepCalculator stepper;
    ODESystemFixedStepSolver solver(system, stepper);
    ODESystemSolution solution = solver.integrate(initial_cond, 0.0, 10.0, 1000);

    // Access final values (index 999 for 1000 points)
    int last = solution.size() - 1;
    std::cout << "Final position: " << solution.getXValue(last, 0) << std::endl;
    std::cout << "Final velocity: " << solution.getXValue(last, 1) << std::endl;

    // Analytical solution at t=10: x(t) = cos(ωt), v(t) = -ω sin(ωt)
    Real t_final = 10.0;
    std::cout << "Analytical x:   " << std::cos(omega * t_final) << std::endl;
    std::cout << "Analytical v:   " << -omega * std::sin(omega * t_final) << std::endl;

/* Expected OUTPUT:
    Final position: 0.4080820860  (RK4 numerical)
    Final velocity: -1.8258904789 (RK4 numerical)
    Analytical x:   0.4080820618  (cos(20))
    Analytical v:   -1.8258905015 (-2*sin(20))
*/
}
