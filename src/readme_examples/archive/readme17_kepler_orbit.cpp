///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme17_kepler_orbit.cpp                                           ///
///  Description: README example - Kepler Orbit Simulation                            ///
///               Demonstrates orbital mechanics with energy conservation check       ///
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
#include <cmath>

using namespace MML;

// Gravitational parameter (GM = 1 for normalized units)
static const Real GM = 1.0;

// Kepler two-body equations of motion
// State vector: [x, y, vx, vy]
// Equations: dx/dt = vx, dy/dt = vy, dvx/dt = -GM*x/r³, dvy/dt = -GM*y/r³
static void kepler_derivs(Real t, const Vector<Real>& state, Vector<Real>& deriv)
{
    (void)t;  // unused
    Real x = state[0];
    Real y = state[1];
    Real vx = state[2];
    Real vy = state[3];

    Real r = std::sqrt(x*x + y*y);
    Real r3 = r * r * r;

    deriv[0] = vx;                    // dx/dt = vx
    deriv[1] = vy;                    // dy/dt = vy
    deriv[2] = -GM * x / r3;          // dvx/dt = -GM*x/r³
    deriv[3] = -GM * y / r3;          // dvy/dt = -GM*y/r³
}

void Readme_KeplerOrbit()
{
    std::cout << std::endl;
    std::cout << "=== Kepler Orbit Simulation ===" << std::endl;
    std::cout << "Simulating planetary motion with energy conservation verification" << std::endl;

    // Function to compute total energy (should be conserved)
    auto computeEnergy = [](Real x, Real y, Real vx, Real vy) -> Real {
        Real r = std::sqrt(x*x + y*y);
        Real KE = 0.5 * (vx*vx + vy*vy);  // Kinetic energy
        Real PE = -GM / r;                 // Potential energy
        return KE + PE;                    // Total energy
    };

    // Function to compute angular momentum (should also be conserved)
    auto computeAngularMomentum = [](Real x, Real y, Real vx, Real vy) -> Real {
        return x * vy - y * vx;  // L = x*vy - y*vx
    };

    // Initial conditions for elliptical orbit (e ≈ 0.5)
    // Start at perihelion with velocity giving eccentricity 0.5
    Real r_peri = 1.0;                           // Perihelion distance
    Real e = 0.5;                                 // Eccentricity
    Real v_peri = std::sqrt(GM * (1 + e) / r_peri);  // Velocity at perihelion

    Vector<Real> initCond({r_peri, 0.0, 0.0, v_peri});

    // Orbital period from Kepler's third law
    Real a = r_peri / (1 - e);                   // Semi-major axis
    Real period = 2 * Constants::PI * std::sqrt(a*a*a / GM);

    std::cout << std::endl << "Orbital Parameters:" << std::endl;
    std::cout << "  Eccentricity e = " << e << std::endl;
    std::cout << "  Semi-major axis a = " << std::setprecision(4) << a << std::endl;
    std::cout << "  Perihelion distance = " << r_peri << std::endl;
    std::cout << "  Aphelion distance = " << a * (1 + e) << std::endl;
    std::cout << "  Orbital period T = " << std::setprecision(6) << period << std::endl;

    // Initial conserved quantities
    Real E0 = computeEnergy(initCond[0], initCond[1], initCond[2], initCond[3]);
    Real L0 = computeAngularMomentum(initCond[0], initCond[1], initCond[2], initCond[3]);

    std::cout << std::endl << "Initial conserved quantities:" << std::endl;
    std::cout << "  Total energy E = " << std::setprecision(10) << E0 << std::endl;
    std::cout << "  Angular momentum L = " << L0 << std::endl;

    // Create ODE system and solver
    ODESystem kepler(4, kepler_derivs);
    RK4_StepCalculator rk4;
    ODESystemFixedStepSolver solver(kepler, rk4);

    // Integrate for one complete orbit
    std::cout << std::endl << "Integrating one complete orbit..." << std::endl;
    std::cout << std::setw(10) << "t/T" << std::setw(12) << "x" 
              << std::setw(12) << "y" << std::setw(12) << "r"
              << std::setw(14) << "ΔE/E0" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    int numSteps = 1000;  // 1000 steps per orbit
    ODESystemSolution sol = solver.integrate(initCond, 0.0, period, numSteps);

    // Sample output and track conservation
    int output_interval = 100;  // Output every 100 steps (10 points per orbit)
    
    for (int i = 0; i <= numSteps; i += output_interval)
    {
        Real t = sol.t[i];
        Real x = sol.y[0][i];
        Real y = sol.y[1][i];
        Real vx = sol.y[2][i];
        Real vy = sol.y[3][i];
        
        Real r = std::sqrt(x*x + y*y);
        Real E = computeEnergy(x, y, vx, vy);
        Real dE_rel = (E - E0) / std::abs(E0);

        std::cout << std::fixed
                  << std::setw(10) << std::setprecision(2) << t/period
                  << std::setw(12) << std::setprecision(4) << x
                  << std::setw(12) << y
                  << std::setw(12) << r
                  << std::scientific << std::setw(14) << std::setprecision(2) << dE_rel
                  << std::endl;
    }

    // Final conservation check
    Real E_final = computeEnergy(sol.y[0][numSteps], sol.y[1][numSteps], 
                                  sol.y[2][numSteps], sol.y[3][numSteps]);
    Real L_final = computeAngularMomentum(sol.y[0][numSteps], sol.y[1][numSteps], 
                                           sol.y[2][numSteps], sol.y[3][numSteps]);

    std::cout << std::string(60, '-') << std::endl;
    std::cout << std::endl << "Conservation check after one orbit:" << std::endl;
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  Energy error ΔE/E0 = " << (E_final - E0) / std::abs(E0) << std::endl;
    std::cout << "  Angular momentum error ΔL/L0 = " << (L_final - L0) / std::abs(L0) << std::endl;

    // Check return to initial position
    Real dx = sol.y[0][numSteps] - r_peri;
    Real dy = sol.y[1][numSteps] - 0.0;
    Real position_error = std::sqrt(dx*dx + dy*dy);
    std::cout << "  Position error after one orbit = " << position_error << std::endl;

    // Demonstrate multiple orbits
    std::cout << std::endl << "--- Long-term stability (10 orbits) ---" << std::endl;

    int steps_10_orbits = 5000;  // Coarser for speed
    ODESystemSolution sol10 = solver.integrate(initCond, 0.0, 10.0 * period, steps_10_orbits);

    int steps_per_orbit = steps_10_orbits / 10;
    for (int orbit = 1; orbit <= 10; orbit++)
    {
        int idx = orbit * steps_per_orbit;
        Real E = computeEnergy(sol10.y[0][idx], sol10.y[1][idx], 
                               sol10.y[2][idx], sol10.y[3][idx]);
        std::cout << "  Orbit " << std::setw(2) << orbit 
                  << ": ΔE/E0 = " << std::scientific << std::setprecision(4) 
                  << (E - E0) / std::abs(E0) << std::endl;
    }

    std::cout << std::endl;
}
