///////////////////////////////////////////////////////////////////////////////
// Test suite for ODESystem.h and ODESystemSolution.h
///////////////////////////////////////////////////////////////////////////////

#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include <string>

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

///////////////////////////////////////////////////////////////////////////////
// Test functions for ODE systems

// Simple 1D exponential growth: dy/dt = y
void expGrowth(Real t, const MML::Vector<Real>& y, MML::Vector<Real>& dydt)
{
    dydt[0] = y[0];
}

// 2D harmonic oscillator: dx/dt = v, dv/dt = -x
void harmonicOscillator(Real t, const MML::Vector<Real>& y, MML::Vector<Real>& dydt)
{
    dydt[0] = y[1];      // dx/dt = v
    dydt[1] = -y[0];     // dv/dt = -x
}

// 3D Lorenz system (simplified parameters)
void lorenzSystem(Real t, const MML::Vector<Real>& y, MML::Vector<Real>& dydt)
{
    Real sigma = REAL(10.0), rho = REAL(28.0), beta = REAL(8.0)/REAL(3.0);
    dydt[0] = sigma * (y[1] - y[0]);
    dydt[1] = y[0] * (rho - y[2]) - y[1];
    dydt[2] = y[0] * y[1] - beta * y[2];
}

// Jacobian for harmonic oscillator
void harmonicJacobian(Real t, const MML::Vector<Real>& y, MML::Vector<Real>& dydt, MML::Matrix<Real>& J)
{
    harmonicOscillator(t, y, dydt);
    J.Resize(2, 2);
    J[0][0] = REAL(0.0);  J[0][1] = REAL(1.0);
    J[1][0] = -REAL(1.0); J[1][1] = REAL(0.0);
}

///////////////////////////////////////////////////////////////////////////////
// ODESystem Tests

TEST_CASE("ODESystem - Default Constructor", "[ode_system]")
{
    ODESystem sys;
    
    REQUIRE(sys.getDim() == 0);
}

TEST_CASE("ODESystem - Constructor with Function", "[ode_system]")
{
    ODESystem sys(1, expGrowth);
    
    REQUIRE(sys.getDim() == 1);
}

TEST_CASE("ODESystem - Exponential Growth Derivatives", "[ode_system]")
{
    ODESystem sys(1, expGrowth);
    
    MML::Vector<Real> y(1);
    y[0] = REAL(2.0);
    
    MML::Vector<Real> dydt(1);
    sys.derivs(REAL(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(REAL(2.0), REAL(1e-10)));
}

TEST_CASE("ODESystem - Harmonic Oscillator Derivatives", "[ode_system]")
{
    ODESystem sys(2, harmonicOscillator);
    
    REQUIRE(sys.getDim() == 2);
    
    MML::Vector<Real> y(2);
    y[0] = REAL(1.0);  // position
    y[1] = REAL(0.0);  // velocity
    
    MML::Vector<Real> dydt(2);
    sys.derivs(REAL(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(REAL(0.0), REAL(1e-10)));  // dx/dt = v = 0
    REQUIRE_THAT(dydt[1], WithinAbs(-REAL(1.0), REAL(1e-10))); // dv/dt = -x = -1
}

TEST_CASE("ODESystem - Lorenz System Derivatives", "[ode_system]")
{
    ODESystem sys(3, lorenzSystem);
    
    REQUIRE(sys.getDim() == 3);
    
    MML::Vector<Real> y(3);
    y[0] = REAL(1.0);
    y[1] = REAL(1.0);
    y[2] = REAL(1.0);
    
    MML::Vector<Real> dydt(3);
    sys.derivs(REAL(0.0), y, dydt);
    
    // sigma * (y - x) = 10 * (1 - 1) = 0
    REQUIRE_THAT(dydt[0], WithinAbs(REAL(0.0), REAL(1e-10)));
    // x * (rho - z) - y = 1 * (28 - 1) - 1 = 26
    REQUIRE_THAT(dydt[1], WithinAbs(REAL(26.0), REAL(1e-10)));
    // x * y - beta * z = 1 * 1 - (8/3) * 1 ≈ -REAL(1.6667)
    REQUIRE_THAT(dydt[2], WithinAbs(-REAL(5.0)/REAL(3.0), REAL(1e-10)));
}

TEST_CASE("ODESystem - Operator() Call", "[ode_system]")
{
    ODESystem sys(1, expGrowth);
    
    MML::Vector<Real> y(1);
    y[0] = REAL(3.0);
    
    MML::Vector<Real> dydt(1);
    sys(REAL(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(REAL(3.0), REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// ODESystemWithJacobian Tests

TEST_CASE("ODESystemWithJacobian - Default Constructor", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys;
    
    REQUIRE(sys.getDim() == 0);
}

TEST_CASE("ODESystemWithJacobian - Constructor with Functions", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    
    REQUIRE(sys.getDim() == 2);
}

TEST_CASE("ODESystemWithJacobian - Harmonic Oscillator Derivatives", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    
    MML::Vector<Real> y(2);
    y[0] = REAL(1.0);
    y[1] = REAL(0.5);
    
    MML::Vector<Real> dydt(2);
    sys.derivs(REAL(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(REAL(0.5), REAL(1e-10)));
    REQUIRE_THAT(dydt[1], WithinAbs(-REAL(1.0), REAL(1e-10)));
}

TEST_CASE("ODESystemWithJacobian - Jacobian Computation", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    
    MML::Vector<Real> y(2);
    y[0] = REAL(1.0);
    y[1] = REAL(0.5);
    
    MML::Vector<Real> dydt(2);
    MML::Matrix<Real> J(2, 2);
    
    sys.jacobian(REAL(0.0), y, dydt, J);
    
    // Check Jacobian elements
    REQUIRE_THAT(J[0][0], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(J[0][1], WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(J[1][0], WithinAbs(-REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(J[1][1], WithinAbs(REAL(0.0), REAL(1e-10)));
    
    // Check that derivatives were also computed
    REQUIRE_THAT(dydt[0], WithinAbs(REAL(0.5), REAL(1e-10)));
    REQUIRE_THAT(dydt[1], WithinAbs(-REAL(1.0), REAL(1e-10)));
}

///////////////////////////////////////////////////////////////////////////////
// ODESystemSolution Tests

TEST_CASE("ODESystemSolution - Constructor with Max Steps", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(10.0), 2, 500);
    
    REQUIRE(sol.getSysDim() == 2);
    REQUIRE_THAT(sol.getT1(), WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(sol.getT2(), WithinAbs(REAL(10.0), REAL(1e-10)));
    REQUIRE(sol.getNumStepsOK() == 0);
    REQUIRE(sol.getNumStepsBad() == 0);
    REQUIRE(sol.getTotalSavedSteps() == 501);  // maxSteps + 1
}

TEST_CASE("ODESystemSolution - Constructor with Default Max Steps", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 3);
    
    REQUIRE(sol.getSysDim() == 3);
    REQUIRE(sol.getTotalSavedSteps() == 1001);  // default 1000 + 1
}

TEST_CASE("ODESystemSolution - Increment Step Counters", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 1);
    
    sol.incrementSuccessfulSteps();
    sol.incrementSuccessfulSteps();
    sol.incrementRejectedSteps();
    
    REQUIRE(sol.getNumStepsOK() == 2);
    REQUIRE(sol.getNumStepsBad() == 1);
    REQUIRE(sol.getTotalNumSteps() == 3);
}

TEST_CASE("ODESystemSolution - Set and Get T Values", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(10.0), 2);
    
    sol.setTVal(0, REAL(0.0));
    sol.setTVal(1, REAL(0.5));
    sol.setTVal(2, REAL(1.0));
    
    MML::Vector<Real> tvals = sol.getTValues();
    
    REQUIRE_THAT(tvals[0], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(tvals[1], WithinAbs(REAL(0.5), REAL(1e-10)));
    REQUIRE_THAT(tvals[2], WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("ODESystemSolution - Set and Get X Values", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 2);
    
    sol.setXVal(0, 0, REAL(1.0));
    sol.setXVal(0, 1, REAL(2.0));
    sol.setXVal(1, 0, REAL(1.5));
    sol.setXVal(1, 1, REAL(2.5));
    
    MML::Matrix<Real> xvals = sol.getXValues();
    
    REQUIRE_THAT(xvals[0][0], WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(xvals[1][0], WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(xvals[0][1], WithinAbs(REAL(1.5), REAL(1e-10)));
    REQUIRE_THAT(xvals[1][1], WithinAbs(REAL(2.5), REAL(1e-10)));
}

TEST_CASE("ODESystemSolution - Get X Values by Component", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 3);
    
    sol.setXVal(0, 0, REAL(1.0));
    sol.setXVal(1, 0, REAL(2.0));
    sol.setXVal(2, 0, REAL(3.0));
    
    sol.setXVal(0, 1, REAL(4.0));
    sol.setXVal(1, 1, REAL(5.0));
    sol.setXVal(2, 1, REAL(6.0));
    
    MML::Vector<Real> comp0 = sol.getXValues(0);
    MML::Vector<Real> comp1 = sol.getXValues(1);
    
    REQUIRE_THAT(comp0[0], WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(comp0[1], WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(comp0[2], WithinAbs(REAL(3.0), REAL(1e-10)));
    
    REQUIRE_THAT(comp1[0], WithinAbs(REAL(4.0), REAL(1e-10)));
    REQUIRE_THAT(comp1[1], WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(comp1[2], WithinAbs(REAL(6.0), REAL(1e-10)));
}

TEST_CASE("ODESystemSolution - Get X Values by Component - Out of Range", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 2);
    
    REQUIRE_THROWS(sol.getXValues(-1));
    REQUIRE_THROWS(sol.getXValues(2));
    REQUIRE_THROWS(sol.getXValues(3));
}

TEST_CASE("ODESystemSolution - Fill Values", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 3);
    
    MML::Vector<Real> y(3);
    y[0] = REAL(1.0);
    y[1] = REAL(2.0);
    y[2] = REAL(3.0);
    
    sol.fillValues(0, REAL(0.0), y);
    
    y[0] = REAL(4.0);
    y[1] = REAL(5.0);
    y[2] = REAL(6.0);
    
    sol.fillValues(1, REAL(0.5), y);
    
    MML::Vector<Real> tvals = sol.getTValues();
    MML::Matrix<Real> xvals = sol.getXValues();
    
    REQUIRE_THAT(tvals[0], WithinAbs(REAL(0.0), REAL(1e-10)));
    REQUIRE_THAT(tvals[1], WithinAbs(REAL(0.5), REAL(1e-10)));
    
    REQUIRE_THAT(xvals[0][0], WithinAbs(REAL(1.0), REAL(1e-10)));
    REQUIRE_THAT(xvals[1][0], WithinAbs(REAL(2.0), REAL(1e-10)));
    REQUIRE_THAT(xvals[2][0], WithinAbs(REAL(3.0), REAL(1e-10)));
    
    REQUIRE_THAT(xvals[0][1], WithinAbs(REAL(4.0), REAL(1e-10)));
    REQUIRE_THAT(xvals[1][1], WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(xvals[2][1], WithinAbs(REAL(6.0), REAL(1e-10)));
}

TEST_CASE("ODESystemSolution - Fill Values - Size Mismatch", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 3);
    
    MML::Vector<Real> y(2);  // Wrong size
    y[0] = REAL(1.0);
    y[1] = REAL(2.0);
    
    REQUIRE_THROWS(sol.fillValues(0, REAL(0.0), y));
}

TEST_CASE("ODESystemSolution - Get X Values at End", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 2, 10);
    
    // Fill the last step (index 10 for maxSteps=10)
    sol.setXVal(10, 0, REAL(5.0));
    sol.setXVal(10, 1, REAL(6.0));
    
    MML::Vector<Real> endVals = sol.getXValuesAtEnd();
    
    REQUIRE_THAT(endVals[0], WithinAbs(REAL(5.0), REAL(1e-10)));
    REQUIRE_THAT(endVals[1], WithinAbs(REAL(6.0), REAL(1e-10)));
}

TEST_CASE("ODESystemSolution - Set Final Size", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 2, 1000);
    
    REQUIRE(sol.getTotalSavedSteps() == 1001);
    
    sol.setFinalSize(50);
    
    REQUIRE(sol.getTotalSavedSteps() == 51);
}

TEST_CASE("ODESystemSolution - Auto-Extend Storage", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 2, 10);
    
    REQUIRE(sol.getTotalSavedSteps() == 11);
    
    // Try to set value beyond initial capacity
    sol.setTVal(15, REAL(1.5));
    
    // Should have auto-extended by 100 steps
    REQUIRE(sol.getTotalSavedSteps() >= 15);
}

TEST_CASE("ODESystemSolution - Linear Interpolation", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(2.0), 1);
    
    // Set up simple linear data: y = 2*t
    sol.setTVal(0, REAL(0.0));
    sol.setXVal(0, 0, REAL(0.0));
    
    sol.setTVal(1, REAL(1.0));
    sol.setXVal(1, 0, REAL(2.0));
    
    sol.setTVal(2, REAL(2.0));
    sol.setXVal(2, 0, REAL(4.0));
    
    sol.setFinalSize(2);
    
    LinearInterpRealFunc interp = sol.getSolAsLinInterp(0);
    
    REQUIRE_THAT(interp(REAL(0.5)), WithinAbs(REAL(1.0), REAL(1e-9)));
    REQUIRE_THAT(interp(REAL(1.5)), WithinAbs(REAL(3.0), REAL(1e-9)));
}

TEST_CASE("ODESystemSolution - 2D Parametric Curve", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 3, 10);
    
    // Set up circle: x = cos(t), y = sin(t)
    for (int i = 0; i <= 10; i++)
    {
        Real t = i * REAL(0.1);
        sol.setTVal(i, t);
        sol.setXVal(i, 0, std::cos(t));
        sol.setXVal(i, 1, std::sin(t));
        sol.setXVal(i, 2, REAL(0.0));
    }
    
    sol.setFinalSize(10);
    
    SplineInterpParametricCurve<2> curve = sol.getSolAsParamCurve2D(0, 1);
    
    // Verify curve returns points near the circle
    VectorN<Real, 2> point = curve(REAL(0.0));
    REQUIRE_THAT(point[0], WithinAbs(REAL(1.0), REAL(1e-1)));  // cos(0) = 1
    REQUIRE_THAT(point[1], WithinAbs(REAL(0.0), REAL(1e-1)));  // sin(0) = 0
}

TEST_CASE("ODESystemSolution - 3D Parametric Curve", "[ode_system_solution]")
{
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 3, 10);
    
    // Set up helix: x = cos(t), y = sin(t), z = t
    for (int i = 0; i <= 10; i++)
    {
        Real t = i * REAL(0.1);
        sol.setTVal(i, t);
        sol.setXVal(i, 0, std::cos(t));
        sol.setXVal(i, 1, std::sin(t));
        sol.setXVal(i, 2, t);
    }
    
    sol.setFinalSize(10);
    
    SplineInterpParametricCurve<3> curve = sol.getSolAsParamCurve3D(0, 1, 2);
    
    // Verify curve returns points near the helix
    VectorN<Real, 3> point = curve(REAL(0.0));
    REQUIRE_THAT(point[0], WithinAbs(REAL(1.0), REAL(1e-1)));  // cos(0) = 1
    REQUIRE_THAT(point[1], WithinAbs(REAL(0.0), REAL(1e-1)));  // sin(0) = 0
    REQUIRE_THAT(point[2], WithinAbs(REAL(0.0), REAL(1e-1)));  // z = 0
}

///////////////////////////////////////////////////////////////////////////////
// Integration Tests

TEST_CASE("ODESystem Integration - Exponential Growth Manual Step", "[ode_system_integration]")
{
    ODESystem sys(1, expGrowth);
    ODESystemSolution sol(REAL(0.0), REAL(1.0), 1);
    
    MML::Vector<Real> y(1);
    y[0] = REAL(1.0);
    
    Real t = REAL(0.0);
    Real dt = REAL(0.01);
    
    sol.fillValues(0, t, y);
    
    // Simple Euler step: y(t+dt) = y(t) + dt * f(t, y)
    for (int i = 1; i <= 10; i++)
    {
        MML::Vector<Real> dydt(1);
        sys.derivs(t, y, dydt);
        
        y[0] += dt * dydt[0];
        t += dt;
        
        sol.fillValues(i, t, y);
    }
    
    sol.setFinalSize(10);
    
    // After small time step, should approximate exp(REAL(0.1)) ≈ REAL(1.1052)
    REQUIRE_THAT(y[0], WithinAbs(REAL(1.1052), REAL(0.01)));
}

TEST_CASE("ODESystem Integration - Harmonic Oscillator Conservation", "[ode_system_integration]")
{
    ODESystem sys(2, harmonicOscillator);
    
    MML::Vector<Real> y(2);
    y[0] = REAL(1.0);  // Initial position
    y[1] = REAL(0.0);  // Initial velocity
    
    // Energy should be conserved: E = REAL(0.5) * (x^2 + v^2)
    Real initialEnergy = REAL(0.5) * (y[0]*y[0] + y[1]*y[1]);
    
    // Take a few Euler steps
    Real t = REAL(0.0);
    Real dt = REAL(0.01);
    
    for (int i = 0; i < 10; i++)
    {
        MML::Vector<Real> dydt(2);
        sys.derivs(t, y, dydt);
        
        y[0] += dt * dydt[0];
        y[1] += dt * dydt[1];
        t += dt;
    }
    
    Real finalEnergy = REAL(0.5) * (y[0]*y[0] + y[1]*y[1]);
    
    // Energy should be approximately conserved (within numerical error)
    REQUIRE_THAT(finalEnergy, WithinAbs(initialEnergy, REAL(0.05)));
}
