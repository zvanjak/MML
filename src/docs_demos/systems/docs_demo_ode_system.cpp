///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_ode_system.cpp                                            ///
///  Description: Documentation examples for ODESystem.h and ODESystemSolution.h     ///
///               ODE system representations and solution containers                  ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                              ODE SYSTEM EXAMPLES                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

// Example 1: Simple harmonic oscillator
// dx/dt = v
// dv/dt = -omega^2 * x
void harmonicOscillator(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    Real omega = 2.0;  // angular frequency
    dxdt[0] = x[1];              // dx/dt = v
    dxdt[1] = -omega*omega * x[0]; // dv/dt = -omega^2 * x
}

// Example 2: Lorenz system (chaotic attractor)
// dx/dt = sigma * (y - x)
// dy/dt = x * (rho - z) - y
// dz/dt = x * y - beta * z
void lorenzSystem(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    Real sigma = 10.0;
    Real rho = 28.0;
    Real beta = 8.0 / 3.0;
    
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
}

// Example 3: Van der Pol oscillator (nonlinear)
// dx/dt = v
// dv/dt = mu * (1 - x^2) * v - x
void vanDerPolOscillator(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    Real mu = 1.0;  // nonlinearity parameter
    dxdt[0] = x[1];
    dxdt[1] = mu * (1.0 - x[0]*x[0]) * x[1] - x[0];
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              ODE SYSTEM DEFINITION                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_ODESystem()
{
    std::cout << "\n=== ODE System Definition ===\n\n";

    std::cout << "ODESystem wraps function pointers for ODE solvers.\n";
    std::cout << "It defines dx/dt = f(t, x) where x is a state vector.\n\n";

    // Create ODE system for harmonic oscillator
    ODESystem harmonic(2, harmonicOscillator);
    
    std::cout << "Created harmonic oscillator ODE system:\n";
    std::cout << "  Dimension: " << harmonic.getDim() << "\n";
    std::cout << "  Equations: dx/dt = v, dv/dt = -omega^2 * x\n\n";

    // Test computing derivatives
    Real t = 0.0;
    Vector<Real> x({1.0, 0.0});  // Initial: x=1, v=0
    Vector<Real> dxdt(2);
    
    harmonic.derivs(t, x, dxdt);
    
    std::cout << "At t=0, x=[1, 0]:\n";
    std::cout << "  dx/dt = " << dxdt << "\n";
    std::cout << "  (v=0, acceleration=-4 since omega=2)\n\n";

    // Alternative: use operator()
    Vector<Real> x2({0.5, 1.0});
    harmonic(t, x2, dxdt);
    std::cout << "At t=0, x=[0.5, 1.0]:\n";
    std::cout << "  dx/dt = " << dxdt << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              LORENZ SYSTEM                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_LorenzSystem()
{
    std::cout << "\n=== Lorenz System (Chaotic) ===\n\n";

    std::cout << "The Lorenz system is a classic example of chaos:\n";
    std::cout << "  dx/dt = sigma * (y - x)\n";
    std::cout << "  dy/dt = x * (rho - z) - y\n";
    std::cout << "  dz/dt = x * y - beta * z\n\n";

    ODESystem lorenz(3, lorenzSystem);
    
    std::cout << "Standard parameters: sigma=10, rho=28, beta=8/3\n";
    std::cout << "System dimension: " << lorenz.getDim() << "\n\n";

    // Compute derivatives at a point
    Vector<Real> x({1.0, 1.0, 1.0});
    Vector<Real> dxdt(3);
    lorenz.derivs(0.0, x, dxdt);
    
    std::cout << "At x=[1, 1, 1]:\n";
    std::cout << "  dx/dt = " << dxdt << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              ODE SOLUTION CONTAINER                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_ODESystemSolution()
{
    std::cout << "\n=== ODE Solution Container ===\n\n";

    std::cout << "ODESystemSolution stores integration results:\n";
    std::cout << "  - Time points t\n";
    std::cout << "  - State vectors x(t)\n";
    std::cout << "  - Integration statistics\n\n";

    // Create solution container
    Real t1 = 0.0, t2 = 10.0;
    int dim = 2;
    
    ODESystemSolution solution(t1, t2, dim);
    
    std::cout << "Created ODESystemSolution:\n";
    std::cout << "  Time interval: [" << t1 << ", " << t2 << "]\n";
    std::cout << "  Dimension: " << solution.getSysDim() << "\n";
    std::cout << "  Initial capacity: " << solution.capacity() << " points\n\n";

    // Simulate filling with harmonic oscillator solution
    int numPoints = 11;
    Real omega = 2.0;
    for (int i = 0; i < numPoints; ++i) {
        Real t = t1 + i * (t2 - t1) / (numPoints - 1);
        Vector<Real> x(dim);
        x[0] = std::cos(omega * t);   // x(t) = cos(omega*t)
        x[1] = -omega * std::sin(omega * t);  // v(t) = -omega*sin(omega*t)
        
        solution.setTVal(i, t);
        solution.setXVal(i, 0, x[0]);
        solution.setXVal(i, 1, x[1]);
    }
    solution.setFinalSize(numPoints);
    
    std::cout << "Filled with " << numPoints << " harmonic oscillator points.\n\n";

    // Access solution data
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Time points: ";
    for (int i = 0; i < solution.size(); ++i) {
        std::cout << solution.getTValue(i) << " ";
    }
    std::cout << "\n\n";

    std::cout << "Position x(t): ";
    for (int i = 0; i < solution.size(); ++i) {
        std::cout << solution.getXValue(i, 0) << " ";
    }
    std::cout << "\n\n";

    // Get final values
    Vector<Real> finalState = solution.getXValuesAtEnd();
    std::cout << "Final state at t=" << t2 << ": " << finalState << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              ODE WITH JACOBIAN                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

// Jacobian for harmonic oscillator
void harmonicJacobian(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& J)
{
    Real omega = 2.0;
    
    // Derivatives
    dxdt[0] = x[1];
    dxdt[1] = -omega*omega * x[0];
    
    // Jacobian: J[i][j] = df_i/dx_j
    // f0 = v      -> df0/dx = 0, df0/dv = 1
    // f1 = -w^2*x -> df1/dx = -w^2, df1/dv = 0
    J(0, 0) = 0.0;           J(0, 1) = 1.0;
    J(1, 0) = -omega*omega;  J(1, 1) = 0.0;
}

void Demo_ODEWithJacobian()
{
    std::cout << "\n=== ODE System with Jacobian ===\n\n";

    std::cout << "ODESystemWithJacobian provides analytical Jacobian:\n";
    std::cout << "  - Required for implicit/stiff solvers (BDF, Rosenbrock)\n";
    std::cout << "  - J[i][j] = df_i/dx_j\n\n";

    ODESystemWithJacobian system(2, harmonicOscillator, harmonicJacobian);
    
    std::cout << "Harmonic oscillator Jacobian:\n";
    std::cout << "  f0 = v          -> df0/dx = 0,     df0/dv = 1\n";
    std::cout << "  f1 = -omega^2*x -> df1/dx = -w^2,  df1/dv = 0\n\n";

    // Compute Jacobian at a point
    Vector<Real> x({1.0, 0.5});
    Vector<Real> dxdt(2);
    Matrix<Real> J(2, 2);
    
    system.jacobian(0.0, x, dxdt, J);
    
    std::cout << "At x=[1, 0.5], Jacobian matrix:\n";
    J.Print(std::cout, 8, 3);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              INTERPOLATION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_SolutionInterpolation()
{
    std::cout << "\n=== Solution Interpolation ===\n\n";

    std::cout << "ODESystemSolution provides interpolation between saved points:\n";
    std::cout << "  - Linear interpolation (fast)\n";
    std::cout << "  - Spline interpolation (smooth)\n\n";

    // Create and fill a solution
    Real t1 = 0.0, t2 = Constants::PI;
    int dim = 1;
    int numPoints = 5;
    
    ODESystemSolution solution(t1, t2, dim, numPoints);
    
    for (int i = 0; i < numPoints; ++i) {
        Real t = t1 + i * (t2 - t1) / (numPoints - 1);
        solution.setTVal(i, t);
        solution.setXVal(i, 0, std::sin(t));
    }
    solution.setFinalSize(numPoints);
    
    // Get linear interpolation
    auto linInterp = solution.getSolAsLinInterp(0);
    
    std::cout << "sin(t) interpolation on [0, pi] with " << numPoints << " points:\n\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "t\t\tExact\t\tLinear Interp\n";
    std::cout << std::string(45, '-') << "\n";
    
    for (Real t = 0.0; t <= Constants::PI + 0.01; t += Constants::PI / 8) {
        Real exact = std::sin(t);
        Real interp = linInterp(t);
        std::cout << t << "\t\t" << exact << "\t\t" << interp << "\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              MAIN DEMO FUNCTION                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_ODESystem()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               ODE SYSTEM DEMO\n";
    std::cout << "               ODESystem.h & ODESystemSolution.h\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "This demo covers:\n";
    std::cout << "  - ODESystem class for defining dx/dt = f(t,x)\n";
    std::cout << "  - Classic examples: harmonic oscillator, Lorenz, Van der Pol\n";
    std::cout << "  - ODESystemSolution for storing integration results\n";
    std::cout << "  - ODESystemWithJacobian for stiff solvers\n";
    std::cout << "  - Solution interpolation\n";
    std::cout << std::string(70, '=') << "\n";

    Demo_ODESystem();
    Demo_LorenzSystem();
    Demo_ODESystemSolution();
    Demo_ODEWithJacobian();
    Demo_SolutionInterpolation();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               END OF ODE SYSTEM DEMO\n";
    std::cout << std::string(70, '=') << "\n";
}
