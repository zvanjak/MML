///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme15_lorenz_attractor.cpp                                       ///
///  Description: README example - Lorenz Strange Attractor simulation                ///
///               Demonstrates ODE system solving with fixed-step RK4                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/ODESystem.h"
#include "algorithms/ODESolvers/ODESystemStepCalculators.h"
#include "algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace MML;

// Lorenz parameters (global for simplicity in function pointer interface)
static const Real SIGMA = 10.0;
static const Real RHO = 28.0;
static const Real BETA = 8.0 / 3.0;

// Lorenz system derivatives
void lorenz_derivs(Real t, const Vector<Real>& state, Vector<Real>& deriv)
{
    deriv[0] = SIGMA * (state[1] - state[0]);           // dx/dt = σ(y - x)
    deriv[1] = state[0] * (RHO - state[2]) - state[1];  // dy/dt = x(ρ - z) - y
    deriv[2] = state[0] * state[1] - BETA * state[2];   // dz/dt = xy - βz
}

void Readme_LorenzAttractor()
{
    std::cout << std::endl;
    std::cout << "=== Lorenz Strange Attractor ===" << std::endl;
    std::cout << "Simulating the Lorenz system with fixed-step RK4" << std::endl;

    std::cout << "Parameters: σ = " << SIGMA << ", ρ = " << RHO << ", β = " << BETA << std::endl;

    // Create ODE system and step calculator
    ODESystem lorenz(3, lorenz_derivs);
    RK4_StepCalculator rk4;
    ODESystemFixedStepSolver solver(lorenz, rk4);

    // Initial conditions and solve
    Vector<Real> initCond{1.0, 1.0, 1.0};
    Real t_end = 50.0;
    int numSteps = 5000;

    std::cout << "Initial state: (" << initCond[0] << ", " << initCond[1] << ", " << initCond[2] << ")" << std::endl;
    std::cout << "Integrating from t=0 to t=" << t_end << " with " << numSteps << " steps" << std::endl;
    std::cout << std::endl;

    // Integrate
    ODESystemSolution sol = solver.integrate(initCond, 0.0, t_end, numSteps);

    // Display sample points
    std::cout << "Sample points along trajectory:" << std::endl;
    std::cout << std::setw(10) << "t" << std::setw(14) << "x" 
              << std::setw(14) << "y" << std::setw(14) << "z" << std::endl;
    std::cout << std::string(52, '-') << std::endl;

    int sample_interval = numSteps / 10;  // 10 samples
    for (int i = 0; i <= numSteps; i += sample_interval)
    {
        std::cout << std::fixed << std::setprecision(4)
                  << std::setw(10) << sol.t[i]
                  << std::setw(14) << sol.y[0][i]
                  << std::setw(14) << sol.y[1][i]
                  << std::setw(14) << sol.y[2][i] << std::endl;
    }

    // Final state
    std::cout << std::string(52, '-') << std::endl;
    std::cout << "Final state at t=" << std::fixed << std::setprecision(4) << sol.t[numSteps] 
              << ": (" << sol.y[0][numSteps] << ", " << sol.y[1][numSteps] << ", " << sol.y[2][numSteps] << ")" << std::endl;

    // Demonstrate sensitivity to initial conditions (butterfly effect)
    std::cout << std::endl << "--- Demonstrating Chaos (Butterfly Effect) ---" << std::endl;
    
    Vector<Real> initCond1{1.0, 1.0, 1.0};
    Vector<Real> initCond2{1.0 + 1e-10, 1.0, 1.0};  // Tiny perturbation

    std::cout << "Initial difference: " << std::scientific << std::setprecision(2) 
              << 1e-10 << std::endl;

    ODESystemSolution sol1 = solver.integrate(initCond1, 0.0, 30.0, 3000);
    ODESystemSolution sol2 = solver.integrate(initCond2, 0.0, 30.0, 3000);

    Real diff = std::sqrt(std::pow(sol1.y[0][3000]-sol2.y[0][3000], 2) +
                          std::pow(sol1.y[1][3000]-sol2.y[1][3000], 2) +
                          std::pow(sol1.y[2][3000]-sol2.y[2][3000], 2));
    std::cout << "Difference at t=30: " << std::fixed << std::setprecision(4) << diff << std::endl;
    std::cout << "Amplification factor: ~" << std::scientific << std::setprecision(1) 
              << diff / 1e-10 << std::endl;

    std::cout << std::endl;
}
