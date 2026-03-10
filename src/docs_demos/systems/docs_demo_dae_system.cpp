///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_dae_system.cpp                                            ///
///  Description: Documentation examples for DAESystem.h                              ///
///               Differential-Algebraic Equation system representations              ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/DAESystem.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                              DAE SOLUTION CONTAINER                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_DAESolution()
{
    std::cout << "\n=== DAE Solution Container ===\n\n";

    std::cout << "DAESolution stores results from DAE integration:\n";
    std::cout << "  - Time points t\n";
    std::cout << "  - Differential state vectors x(t)\n";
    std::cout << "  - Algebraic state vectors y(t)\n\n";

    // Create solution container: t in [0, 10], 2 differential vars, 1 algebraic var
    Real t1 = 0.0, t2 = 10.0;
    int diffDim = 2;  // e.g., position and velocity
    int algDim = 1;   // e.g., constraint force
    
    DAESolution solution(t1, t2, diffDim, algDim);
    
    std::cout << "Created DAESolution:\n";
    std::cout << "  Time interval: [" << t1 << ", " << t2 << "]\n";
    std::cout << "  Differential dimensions: " << diffDim << "\n";
    std::cout << "  Algebraic dimensions: " << algDim << "\n\n";

    // Simulate filling with some data points
    int numPoints = 5;
    for (int i = 0; i < numPoints; ++i) {
        Real t = t1 + i * (t2 - t1) / (numPoints - 1);
        Vector<Real> x(diffDim);
        x[0] = std::cos(t);      // position
        x[1] = -std::sin(t);     // velocity
        Vector<Real> y(algDim);
        y[0] = std::cos(t);      // constraint force (example)
        
        solution.fillValues(i, t, x, y);
    }
    solution.setFinalSize(numPoints - 1);
    
    std::cout << "Filled with " << numPoints << " data points.\n";
    std::cout << "  getTotalSavedSteps(): " << solution.getTotalSavedSteps() << "\n\n";

    // Access data
    std::cout << "Time points: ";
    for (int i = 0; i < solution.getTotalSavedSteps(); ++i) {
        std::cout << std::fixed << std::setprecision(2) << solution.getTValue(i) << " ";
    }
    std::cout << "\n";

    std::cout << "x[0] values (position): ";
    for (int i = 0; i < solution.getTotalSavedSteps(); ++i) {
        std::cout << std::fixed << std::setprecision(4) << solution.getXValue(i, 0) << " ";
    }
    std::cout << "\n";

    std::cout << "y[0] values (constraint): ";
    for (int i = 0; i < solution.getTotalSavedSteps(); ++i) {
        std::cout << std::fixed << std::setprecision(4) << solution.getYValue(i, 0) << " ";
    }
    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              DAE SYSTEM DEFINITION                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

// Example: Simple pendulum as DAE
// Differential equations: dx/dt = v, dv/dt = -g*sin(theta) + lambda*x/L
// Algebraic constraint: x^2 + y^2 = L^2

void pendulumDiffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& dxdt)
{
    // x[0] = position x, x[1] = position y, x[2] = velocity vx, x[3] = velocity vy
    // y[0] = constraint force lambda
    Real L = 1.0;  // pendulum length
    Real g = 9.81; // gravity
    
    dxdt[0] = x[2];  // dx/dt = vx
    dxdt[1] = x[3];  // dy/dt = vy
    dxdt[2] = y[0] * x[0] / L;           // dvx/dt = lambda * x / L
    dxdt[3] = -g + y[0] * x[1] / L;      // dvy/dt = -g + lambda * y / L
}

void pendulumAlgConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g)
{
    // Constraint: x^2 + y^2 - L^2 = 0
    Real L = 1.0;
    g[0] = x[0]*x[0] + x[1]*x[1] - L*L;
}

void Demo_DAESystem()
{
    std::cout << "\n=== DAE System Definition ===\n\n";

    std::cout << "DAESystem wraps function pointers for DAE solvers.\n";
    std::cout << "It defines:\n";
    std::cout << "  - Differential equations: dx/dt = f(t, x, y)\n";
    std::cout << "  - Algebraic constraints: 0 = g(t, x, y)\n\n";

    // Create DAE system for pendulum
    int diffDim = 4;  // x, y, vx, vy
    int algDim = 1;   // lambda (constraint force)
    
    DAESystem pendulum(diffDim, algDim, pendulumDiffEqs, pendulumAlgConstraints);
    
    std::cout << "Created pendulum DAE system:\n";
    std::cout << "  Differential dimension: " << pendulum.getDiffDim() << "\n";
    std::cout << "  Algebraic dimension: " << pendulum.getAlgDim() << "\n\n";

    // Test the system
    Real t = 0.0;
    Vector<Real> x({0.0, -1.0, 0.1, 0.0});  // Initial state: hanging down, small velocity
    Vector<Real> y({9.81});                  // Initial constraint force (approximate)
    
    Vector<Real> dxdt(diffDim);
    pendulum.diffEqs(t, x, y, dxdt);
    
    std::cout << "At t = 0, x = " << x << ", y = " << y << ":\n";
    std::cout << "  dx/dt = " << dxdt << "\n";

    Vector<Real> constraint(algDim);
    pendulum.algConstraints(t, x, y, constraint);
    std::cout << "  Constraint residual: " << constraint << "\n";
    std::cout << "  (Should be ~0 if constraint is satisfied)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DAE SYSTEM WITH JACOBIAN                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_DAESystemWithJacobian()
{
    std::cout << "\n=== DAE System with Jacobian ===\n\n";

    std::cout << "DAESystemWithJacobian extends DAESystem with Jacobian matrices:\n";
    std::cout << "  - df/dx: Jacobian of differential eqs w.r.t. diff variables\n";
    std::cout << "  - df/dy: Jacobian of differential eqs w.r.t. alg variables\n";
    std::cout << "  - dg/dx: Jacobian of constraints w.r.t. diff variables\n";
    std::cout << "  - dg/dy: Jacobian of constraints w.r.t. alg variables\n\n";

    std::cout << "Jacobians are required for implicit solvers (BDF, Radau, etc.).\n";
    std::cout << "They enable faster convergence of Newton iterations.\n\n";

    std::cout << "Example structure for pendulum:\n";
    std::cout << "  df/dx is a 4x4 matrix\n";
    std::cout << "  df/dy is a 4x1 matrix\n";
    std::cout << "  dg/dx is a 1x4 matrix\n";
    std::cout << "  dg/dy is a 1x1 matrix (typically zero for index-1 DAEs)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              DAE CONCEPTS                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_DAEConcepts()
{
    std::cout << "\n=== DAE Concepts ===\n\n";

    std::cout << "DIFFERENTIAL-ALGEBRAIC EQUATIONS (DAEs):\n\n";
    
    std::cout << "General form:\n";
    std::cout << "  dx/dt = f(t, x, y)    (differential equations)\n";
    std::cout << "  0 = g(t, x, y)        (algebraic constraints)\n\n";

    std::cout << "Where:\n";
    std::cout << "  x(t) = differential variables (have derivatives)\n";
    std::cout << "  y(t) = algebraic variables (determined by constraints)\n\n";

    std::cout << "DAE INDEX:\n";
    std::cout << "  Index 0: Pure ODE (no algebraic constraints)\n";
    std::cout << "  Index 1: dg/dy is invertible (most common)\n";
    std::cout << "  Index 2: Need to differentiate g once\n";
    std::cout << "  Index 3: Need to differentiate g twice (e.g., position constraints)\n\n";

    std::cout << "COMMON EXAMPLES:\n";
    std::cout << "  - Constrained mechanical systems (pendulum, robots)\n";
    std::cout << "  - Electrical circuits with ideal components\n";
    std::cout << "  - Chemical reaction networks at equilibrium\n";
    std::cout << "  - Power system dynamics\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              MAIN DEMO FUNCTION                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_DAESystem()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               DAE SYSTEM DEMO\n";
    std::cout << "               DAESystem.h\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "This demo covers:\n";
    std::cout << "  - DAESolution container for storing results\n";
    std::cout << "  - DAESystem class for defining DAE problems\n";
    std::cout << "  - DAESystemWithJacobian for implicit solvers\n";
    std::cout << "  - DAE concepts and terminology\n";
    std::cout << std::string(70, '=') << "\n";

    Demo_DAESolution();
    Demo_DAESystem();
    Demo_DAESystemWithJacobian();
    Demo_DAEConcepts();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "               END OF DAE SYSTEM DEMO\n";
    std::cout << std::string(70, '=') << "\n";
}
