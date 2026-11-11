///////////////////////////////////////////////////////////////////////////////
// Test suite for ODESystem.h and ODESystemSolution.h
///////////////////////////////////////////////////////////////////////////////

#include <mml/base/ODESystem.h>
#include <mml/base/ODESystemSolution.h>
#include <mml/base/DAESystem.h>
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include <cmath>
#include <string>
#include <type_traits>

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

namespace MML::Tests::Base::ODESystemTests {

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
    Real sigma = static_cast<MML::Real>(10.0), rho = static_cast<MML::Real>(28.0), beta = static_cast<MML::Real>(8.0)/static_cast<MML::Real>(3.0);
    dydt[0] = sigma * (y[1] - y[0]);
    dydt[1] = y[0] * (rho - y[2]) - y[1];
    dydt[2] = y[0] * y[1] - beta * y[2];
}

// Jacobian for harmonic oscillator
void harmonicJacobian(Real t, const MML::Vector<Real>& y, MML::Vector<Real>& dydt, MML::Matrix<Real>& J)
{
    harmonicOscillator(t, y, dydt);
    J.resize(2, 2);
    J[0][0] = static_cast<MML::Real>(0.0);  J[0][1] = static_cast<MML::Real>(1.0);
    J[1][0] = -static_cast<MML::Real>(1.0); J[1][1] = static_cast<MML::Real>(0.0);
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
    y[0] = static_cast<MML::Real>(2.0);
    
    MML::Vector<Real> dydt(1);
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystem - Harmonic Oscillator Derivatives", "[ode_system]")
{
    ODESystem sys(2, harmonicOscillator);
    
    REQUIRE(sys.getDim() == 2);
    
    MML::Vector<Real> y(2);
    y[0] = static_cast<MML::Real>(1.0);  // position
    y[1] = static_cast<MML::Real>(0.0);  // velocity
    
    MML::Vector<Real> dydt(2);
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));  // dx/dt = v = 0
    REQUIRE_THAT(dydt[1], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5))); // dv/dt = -x = -1
}

TEST_CASE("ODESystem - Lorenz System Derivatives", "[ode_system]")
{
    ODESystem sys(3, lorenzSystem);
    
    REQUIRE(sys.getDim() == 3);
    
    MML::Vector<Real> y(3);
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(1.0);
    y[2] = static_cast<MML::Real>(1.0);
    
    MML::Vector<Real> dydt(3);
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    // sigma * (y - x) = 10 * (1 - 1) = 0
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    // x * (rho - z) - y = 1 * (28 - 1) - 1 = 26
    REQUIRE_THAT(dydt[1], WithinAbs(static_cast<MML::Real>(26.0), TOL(1e-10, 1e-5)));
    // x * y - beta * z = 1 * 1 - (8/3) * 1 ≈ -static_cast<MML::Real>(1.6667)
    REQUIRE_THAT(dydt[2], WithinAbs(-static_cast<MML::Real>(5.0)/static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystem - Operator() Call", "[ode_system]")
{
    ODESystem sys(1, expGrowth);
    
    MML::Vector<Real> y(1);
    y[0] = static_cast<MML::Real>(3.0);
    
    MML::Vector<Real> dydt(1);
    sys(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
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
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(0.5);
    
    MML::Vector<Real> dydt(2);
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(dydt[1], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemWithJacobian - Jacobian Computation", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    
    MML::Vector<Real> y(2);
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(0.5);
    
    MML::Vector<Real> dydt(2);
    MML::Matrix<Real> J(2, 2);
    
    sys.derivsAndJacobian(static_cast<MML::Real>(0.0), y, dydt, J);
    
    // Check Jacobian elements
    REQUIRE_THAT(J[0][0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(J[0][1], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(J[1][0], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(J[1][1], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    
    // Check that derivatives were also computed
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(dydt[1], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
}

///////////////////////////////////////////////////////////////////////////////
// ODESystemSolution Tests

TEST_CASE("ODESystemSolution - Constructor with Max Steps", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 2, 500);
    
    REQUIRE(sol.getSysDim() == 2);
    REQUIRE_THAT(sol.getT1(), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getT2(), WithinAbs(static_cast<MML::Real>(10.0), TOL(1e-10, 1e-5)));
    REQUIRE(sol.getNumStepsOK() == 0);
    REQUIRE(sol.getNumStepsBad() == 0);
    REQUIRE(sol.getTotalSavedSteps() == 501);  // maxSteps + 1
}

TEST_CASE("ODESystemSolution - Constructor with Default Max Steps", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
    
    REQUIRE(sol.getSysDim() == 3);
    REQUIRE(sol.getTotalSavedSteps() == 1001);  // default 1000 + 1
}

TEST_CASE("ODESystemSolution - Increment Step Counters", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 1);
    
    sol.incrementSuccessfulSteps();
    sol.incrementSuccessfulSteps();
    sol.incrementRejectedSteps();
    
    REQUIRE(sol.getNumStepsOK() == 2);
    REQUIRE(sol.getNumStepsBad() == 1);
    REQUIRE(sol.getTotalNumSteps() == 3);
}

TEST_CASE("ODESystemSolution - Set and Get T Values", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 2);
    
    sol.setTVal(0, static_cast<MML::Real>(0.0));
    sol.setTVal(1, static_cast<MML::Real>(0.5));
    sol.setTVal(2, static_cast<MML::Real>(1.0));
    
    MML::Vector<Real> tvals = sol.getTValues();
    
    REQUIRE_THAT(tvals[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(tvals[1], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(tvals[2], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - Set and Get X Values", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2);
    
    sol.setXVal(0, 0, static_cast<MML::Real>(1.0));
    sol.setXVal(0, 1, static_cast<MML::Real>(2.0));
    sol.setXVal(1, 0, static_cast<MML::Real>(1.5));
    sol.setXVal(1, 1, static_cast<MML::Real>(2.5));
    
    MML::Matrix<Real> xvals = sol.getXValues();
    
    REQUIRE_THAT(xvals[0][0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[1][0], WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[0][1], WithinAbs(static_cast<MML::Real>(1.5), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[1][1], WithinAbs(static_cast<MML::Real>(2.5), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - Get X Values by Component", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
    
    sol.setXVal(0, 0, static_cast<MML::Real>(1.0));
    sol.setXVal(1, 0, static_cast<MML::Real>(2.0));
    sol.setXVal(2, 0, static_cast<MML::Real>(3.0));
    
    sol.setXVal(0, 1, static_cast<MML::Real>(4.0));
    sol.setXVal(1, 1, static_cast<MML::Real>(5.0));
    sol.setXVal(2, 1, static_cast<MML::Real>(6.0));
    
    MML::Vector<Real> comp0 = sol.getXValues(0);
    MML::Vector<Real> comp1 = sol.getXValues(1);
    
    REQUIRE_THAT(comp0[0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(comp0[1], WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(comp0[2], WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
    
    REQUIRE_THAT(comp1[0], WithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(comp1[1], WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(comp1[2], WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - Get X Values by Component - Out of Range", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2);
    
    REQUIRE_THROWS(sol.getXValues(-1));
    REQUIRE_THROWS(sol.getXValues(2));
    REQUIRE_THROWS(sol.getXValues(3));
}

TEST_CASE("ODESystemSolution - Fill Values", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
    
    MML::Vector<Real> y(3);
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(2.0);
    y[2] = static_cast<MML::Real>(3.0);
    
    sol.fillValues(0, static_cast<MML::Real>(0.0), y);
    
    y[0] = static_cast<MML::Real>(4.0);
    y[1] = static_cast<MML::Real>(5.0);
    y[2] = static_cast<MML::Real>(6.0);
    
    sol.fillValues(1, static_cast<MML::Real>(0.5), y);
    
    MML::Vector<Real> tvals = sol.getTValues();
    MML::Matrix<Real> xvals = sol.getXValues();
    
    REQUIRE_THAT(tvals[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(tvals[1], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-10, 1e-5)));
    
    REQUIRE_THAT(xvals[0][0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[1][0], WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[2][0], WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
    
    REQUIRE_THAT(xvals[0][1], WithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[1][1], WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(xvals[2][1], WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - Fill Values - Size Mismatch", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
    
    MML::Vector<Real> y(2);  // Wrong size
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(2.0);
    
    REQUIRE_THROWS(sol.fillValues(0, static_cast<MML::Real>(0.0), y));
}

TEST_CASE("ODESystemSolution - Get X Values at End", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    // Fill the last step (index 10 for maxSteps=10)
    sol.setXVal(10, 0, static_cast<MML::Real>(5.0));
    sol.setXVal(10, 1, static_cast<MML::Real>(6.0));
    
    MML::Vector<Real> endVals = sol.getXValuesAtEnd();
    
    REQUIRE_THAT(endVals[0], WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(endVals[1], WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - Set Final Size", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 1000);
    
    REQUIRE(sol.getTotalSavedSteps() == 1001);
    
    sol.setFinalSize(50);
    
    REQUIRE(sol.getTotalSavedSteps() == 51);
}

TEST_CASE("ODESystemSolution - Auto-Extend Storage", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    REQUIRE(sol.getTotalSavedSteps() == 11);
    
    // Try to set value beyond initial capacity
    sol.setTVal(15, static_cast<MML::Real>(1.5));
    
    // Should have auto-extended by 100 steps
    REQUIRE(sol.getTotalSavedSteps() >= 15);
}

TEST_CASE("ODESystemSolution - Linear Interpolation", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), 1);
    
    // Set up simple linear data: y = 2*t
    sol.setTVal(0, static_cast<MML::Real>(0.0));
    sol.setXVal(0, 0, static_cast<MML::Real>(0.0));
    
    sol.setTVal(1, static_cast<MML::Real>(1.0));
    sol.setXVal(1, 0, static_cast<MML::Real>(2.0));
    
    sol.setTVal(2, static_cast<MML::Real>(2.0));
    sol.setXVal(2, 0, static_cast<MML::Real>(4.0));
    
    sol.setFinalSize(2);
    
    LinearInterpRealFunc interp = sol.getSolAsLinInterp(0);
    
    REQUIRE_THAT(interp(static_cast<MML::Real>(0.5)), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-9, 1e-4)));
    REQUIRE_THAT(interp(static_cast<MML::Real>(1.5)), WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-9, 1e-4)));
}

TEST_CASE("ODESystemSolution - 2D Parametric Curve", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3, 10);
    
    // Set up circle: x = cos(t), y = sin(t)
    for (int i = 0; i <= 10; i++)
    {
        Real t = i * static_cast<MML::Real>(0.1);
        sol.setTVal(i, t);
        sol.setXVal(i, 0, std::cos(t));
        sol.setXVal(i, 1, std::sin(t));
        sol.setXVal(i, 2, static_cast<MML::Real>(0.0));
    }
    
    sol.setFinalSize(10);
    
    SplineInterpParametricCurve<2> curve = sol.getSolAsParamCurve2D(0, 1);
    
    // Verify curve returns points near the circle
    VectorN<Real, 2> point = curve(static_cast<MML::Real>(0.0));
    REQUIRE_THAT(point[0], WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-1)));  // cos(0) = 1
    REQUIRE_THAT(point[1], WithinAbs(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-1)));  // sin(0) = 0
}

TEST_CASE("ODESystemSolution - 3D Parametric Curve", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3, 10);
    
    // Set up helix: x = cos(t), y = sin(t), z = t
    for (int i = 0; i <= 10; i++)
    {
        Real t = i * static_cast<MML::Real>(0.1);
        sol.setTVal(i, t);
        sol.setXVal(i, 0, std::cos(t));
        sol.setXVal(i, 1, std::sin(t));
        sol.setXVal(i, 2, t);
    }
    
    sol.setFinalSize(10);
    
    SplineInterpParametricCurve<3> curve = sol.getSolAsParamCurve3D(0, 1, 2);
    
    // Verify curve returns points near the helix
    VectorN<Real, 3> point = curve(static_cast<MML::Real>(0.0));
    REQUIRE_THAT(point[0], WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-1)));  // cos(0) = 1
    REQUIRE_THAT(point[1], WithinAbs(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-1)));  // sin(0) = 0
    REQUIRE_THAT(point[2], WithinAbs(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-1)));  // z = 0
}

///////////////////////////////////////////////////////////////////////////////
// Integration Tests

TEST_CASE("ODESystem Integration - Exponential Growth Manual Step", "[ode_system_integration]")
{
    ODESystem sys(1, expGrowth);
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 1);
    
    MML::Vector<Real> y(1);
    y[0] = static_cast<MML::Real>(1.0);
    
    Real t = static_cast<MML::Real>(0.0);
    Real dt = static_cast<MML::Real>(0.01);
    
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
    
    // After small time step, should approximate exp(static_cast<MML::Real>(0.1)) ≈ static_cast<MML::Real>(1.1052)
    REQUIRE_THAT(y[0], WithinAbs(static_cast<MML::Real>(1.1052), static_cast<MML::Real>(0.01)));
}

TEST_CASE("ODESystem Integration - Harmonic Oscillator Conservation", "[ode_system_integration]")
{
    ODESystem sys(2, harmonicOscillator);
    
    MML::Vector<Real> y(2);
    y[0] = static_cast<MML::Real>(1.0);  // Initial position
    y[1] = static_cast<MML::Real>(0.0);  // Initial velocity
    
    // Energy should be conserved: E = static_cast<MML::Real>(0.5) * (x^2 + v^2)
    Real initialEnergy = static_cast<MML::Real>(0.5) * (y[0]*y[0] + y[1]*y[1]);
    
    // Take a few Euler steps
    Real t = static_cast<MML::Real>(0.0);
    Real dt = static_cast<MML::Real>(0.01);
    
    for (int i = 0; i < 10; i++)
    {
        MML::Vector<Real> dydt(2);
        sys.derivs(t, y, dydt);
        
        y[0] += dt * dydt[0];
        y[1] += dt * dydt[1];
        t += dt;
    }
    
    Real finalEnergy = static_cast<MML::Real>(0.5) * (y[0]*y[0] + y[1]*y[1]);
    
    // Energy should be approximately conserved (within numerical error)
    REQUIRE_THAT(finalEnergy, WithinAbs(initialEnergy, static_cast<MML::Real>(0.05)));
}

///////////////////////////////////////////////////////////////////////////////
// Additional ODESystem Tests

TEST_CASE("ODESystem - Copy and Move Semantics", "[ode_system]")
{
    ODESystem sys1(1, expGrowth);
    
    // Copy constructor
    ODESystem sys2(sys1);
    REQUIRE(sys2.getDim() == 1);
    
    MML::Vector<Real> y(1), dydt(1);
    y[0] = static_cast<MML::Real>(5.0);
    sys2.derivs(static_cast<MML::Real>(0.0), y, dydt);
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
    
    // Copy assignment
    ODESystem sys3;
    sys3 = sys1;
    REQUIRE(sys3.getDim() == 1);
    sys3.derivs(static_cast<MML::Real>(0.0), y, dydt);
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystem - Different Time Values", "[ode_system]")
{
    // Time-dependent ODE: dy/dt = t (linear growth)
    auto linearGrowth = [](Real t, const MML::Vector<Real>& y, MML::Vector<Real>& dydt) {
        dydt[0] = t;
    };
    
    ODESystem sys(1, +linearGrowth);  // Unary + to convert lambda to function pointer
    
    MML::Vector<Real> y(1), dydt(1);
    y[0] = static_cast<MML::Real>(0.0);
    
    sys.derivs(static_cast<MML::Real>(1.0), y, dydt);
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    
    sys.derivs(static_cast<MML::Real>(5.0), y, dydt);
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemWithJacobian - Copy Semantics", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys1(2, harmonicOscillator, harmonicJacobian);
    
    // Copy constructor
    ODESystemWithJacobian sys2(sys1);
    REQUIRE(sys2.getDim() == 2);
    
    MML::Vector<Real> y(2), dydt(2);
    MML::Matrix<Real> J(2, 2);
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(0.0);
    
    sys2.jacobian(static_cast<MML::Real>(0.0), y, J);
    REQUIRE_THAT(J[0][1], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(J[1][0], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemWithJacobian - Operator() Call", "[ode_system_jacobian]")
{
    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    
    MML::Vector<Real> y(2), dydt(2);
    y[0] = static_cast<MML::Real>(2.0);
    y[1] = static_cast<MML::Real>(3.0);
    
    sys(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));  // dx/dt = v = 3
    REQUIRE_THAT(dydt[1], WithinAbs(-static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5))); // dv/dt = -x = -2
}

///////////////////////////////////////////////////////////////////////////////
// Additional ODESystemSolution Tests

TEST_CASE("ODESystemSolution - Size and Capacity Methods", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 50);
    
    REQUIRE(sol.size() == 0);       // no data stored yet
    REQUIRE(sol.capacity() == 51);  // capacity is maxSteps + 1
    REQUIRE(sol.isEmpty());
    REQUIRE(sol.getNumSteps() == 0);
}

TEST_CASE("ODESystemSolution - Reserve Method", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    REQUIRE(sol.capacity() == 11);
    
    sol.reserve(100);
    
    REQUIRE(sol.capacity() == 100);
    
    // Reserve smaller shouldn't shrink
    sol.reserve(50);
    REQUIRE(sol.capacity() == 100);
}

TEST_CASE("ODESystemSolution - GetTValue and GetXValue", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    sol.setTVal(5, static_cast<MML::Real>(0.5));
    sol.setXVal(5, 0, static_cast<MML::Real>(1.5));
    sol.setXVal(5, 1, static_cast<MML::Real>(2.5));
    
    REQUIRE_THAT(sol.getTValue(5), WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getXValue(5, 0), WithinAbs(static_cast<MML::Real>(1.5), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getXValue(5, 1), WithinAbs(static_cast<MML::Real>(2.5), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - GetTValue Out of Range", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    REQUIRE_THROWS_AS(sol.getTValue(-1), IndexError);
    REQUIRE_THROWS_AS(sol.getTValue(100), IndexError);
}

TEST_CASE("ODESystemSolution - GetXValue Out of Range", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    // Index out of range
    REQUIRE_THROWS_AS(sol.getXValue(-1, 0), IndexError);
    REQUIRE_THROWS_AS(sol.getXValue(100, 0), IndexError);
    
    // Component out of range
    REQUIRE_THROWS_AS(sol.getXValue(0, -1), IndexError);
    REQUIRE_THROWS_AS(sol.getXValue(0, 2), IndexError);
}

TEST_CASE("ODESystemSolution - SetTVal Negative Index", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    REQUIRE_THROWS_AS(sol.setTVal(-1, static_cast<MML::Real>(0.5)), IndexError);
}

TEST_CASE("ODESystemSolution - SetXVal Invalid Indices", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    // Negative index
    REQUIRE_THROWS_AS(sol.setXVal(-1, 0, static_cast<MML::Real>(1.0)), IndexError);
    
    // Component out of range
    REQUIRE_THROWS_AS(sol.setXVal(0, -1, static_cast<MML::Real>(1.0)), IndexError);
    REQUIRE_THROWS_AS(sol.setXVal(0, 2, static_cast<MML::Real>(1.0)), IndexError);
}

TEST_CASE("ODESystemSolution - FillValues Negative Index", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 10);
    
    MML::Vector<Real> y(2);
    y[0] = static_cast<MML::Real>(1.0);
    y[1] = static_cast<MML::Real>(2.0);
    
    REQUIRE_THROWS_AS(sol.fillValues(-1, static_cast<MML::Real>(0.0), y), IndexError);
}

TEST_CASE("ODESystemSolution - Polynomial Interpolation", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), 1);
    
    // Set up quadratic data: y = t^2
    sol.setTVal(0, static_cast<MML::Real>(0.0));
    sol.setXVal(0, 0, static_cast<MML::Real>(0.0));
    
    sol.setTVal(1, static_cast<MML::Real>(1.0));
    sol.setXVal(1, 0, static_cast<MML::Real>(1.0));
    
    sol.setTVal(2, static_cast<MML::Real>(2.0));
    sol.setXVal(2, 0, static_cast<MML::Real>(4.0));
    
    sol.setFinalSize(2);
    
    PolynomInterpRealFunc interp = sol.getSolAsPolyInterp(0, 2);  // Quadratic
    
    REQUIRE_THAT(interp(static_cast<MML::Real>(0.5)), WithinAbs(static_cast<MML::Real>(0.25), TOL(1e-8, 1e-4)));
    REQUIRE_THAT(interp(static_cast<MML::Real>(1.5)), WithinAbs(static_cast<MML::Real>(2.25), TOL(1e-8, 1e-4)));
}

TEST_CASE("ODESystemSolution - Spline Interpolation", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), 1);
    
    // Set up sinusoidal data
    for (int i = 0; i <= 10; i++)
    {
        Real t = i * static_cast<MML::Real>(0.2);
        sol.setTVal(i, t);
        sol.setXVal(i, 0, std::sin(t));
    }
    
    sol.setFinalSize(10);
    
    SplineInterpRealFunc interp = sol.getSolAsSplineInterp(0);
    
    // Spline should approximate sin well
    REQUIRE_THAT(interp(static_cast<MML::Real>(0.5)), WithinAbs(std::sin(static_cast<MML::Real>(0.5)), static_cast<MML::Real>(0.01)));
    REQUIRE_THAT(interp(static_cast<MML::Real>(1.5)), WithinAbs(std::sin(static_cast<MML::Real>(1.5)), static_cast<MML::Real>(0.01)));
}

TEST_CASE("ODESystemSolution - Multiple Extensions", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(100.0), 1, 5);
    
    REQUIRE(sol.getTotalSavedSteps() == 6);
    
    // Force multiple extensions
    for (int i = 0; i < 50; i++)
    {
        sol.setTVal(i, static_cast<Real>(i));
        sol.setXVal(i, 0, static_cast<Real>(i * 2));
    }
    
    // Verify data integrity after extensions
    REQUIRE_THAT(sol.getTValue(0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getTValue(25), WithinAbs(static_cast<MML::Real>(25.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getTValue(49), WithinAbs(static_cast<MML::Real>(49.0), TOL(1e-10, 1e-5)));
    
    REQUIRE_THAT(sol.getXValue(0, 0), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getXValue(25, 0), WithinAbs(static_cast<MML::Real>(50.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(sol.getXValue(49, 0), WithinAbs(static_cast<MML::Real>(98.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemSolution - High Dimension System", "[ode_system_solution]")
{
    int dim = 10;
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), dim, 20);
    
    REQUIRE(sol.getSysDim() == dim);
    
    // Fill all components
    MML::Vector<Real> y(dim);
    for (int i = 0; i < dim; i++)
        y[i] = static_cast<Real>(i + 1);
    
    sol.fillValues(0, static_cast<MML::Real>(0.0), y);
    
    // Verify all components
    for (int i = 0; i < dim; i++)
    {
        REQUIRE_THAT(sol.getXValue(0, i), WithinAbs(static_cast<Real>(i + 1), TOL(1e-10, 1e-5)));
    }
    
    // Get end values
    MML::Vector<Real> endVals = sol.getXValuesAtEnd();
    REQUIRE(endVals.size() == dim);
}

TEST_CASE("ODESystemSolution - Step Counter Reset", "[ode_system_solution]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 1);
    
    // Initial state
    REQUIRE(sol.getNumStepsOK() == 0);
    REQUIRE(sol.getNumStepsBad() == 0);
    REQUIRE(sol.getTotalNumSteps() == 0);
    
    // Increment counters
    for (int i = 0; i < 100; i++)
        sol.incrementSuccessfulSteps();
    for (int i = 0; i < 5; i++)
        sol.incrementRejectedSteps();
    
    REQUIRE(sol.getNumStepsOK() == 100);
    REQUIRE(sol.getNumStepsBad() == 5);
    REQUIRE(sol.getTotalNumSteps() == 105);
}

///////////////////////////////////////////////////////////////////////////////
// Edge Case Tests

TEST_CASE("ODESystem - Zero State Vector", "[ode_system]")
{
    ODESystem sys(2, harmonicOscillator);
    
    MML::Vector<Real> y(2), dydt(2);
    y[0] = static_cast<MML::Real>(0.0);
    y[1] = static_cast<MML::Real>(0.0);
    
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(dydt[1], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystem - Large State Values", "[ode_system]")
{
    ODESystem sys(1, expGrowth);
    
    MML::Vector<Real> y(1), dydt(1);
    y[0] = static_cast<MML::Real>(1e10);
    
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(1e10), static_cast<MML::Real>(1.0)));
}

TEST_CASE("ODESystem - Negative State Values", "[ode_system]")
{
    ODESystem sys(2, harmonicOscillator);
    
    MML::Vector<Real> y(2), dydt(2);
    y[0] = -static_cast<MML::Real>(3.0);
    y[1] = -static_cast<MML::Real>(2.0);
    
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    
    REQUIRE_THAT(dydt[0], WithinAbs(-static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));  // dx/dt = v = -2
    REQUIRE_THAT(dydt[1], WithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));   // dv/dt = -x = 3
}

///////////////////////////////////////////////////////////////////////////////
// ODESystemFromStdFunc Tests

TEST_CASE("ODESystemFromStdFunc - Lambda with capture", "[ode_system_stdfunc]")
{
    Real k = static_cast<MML::Real>(2.0);  // captured parameter
    ODESystemFromStdFunc sys(1, [k](Real, const Vector<Real>& y, Vector<Real>& dydt) {
        dydt[0] = k * y[0];
    });

    REQUIRE(sys.getDim() == 1);

    Vector<Real> y(1), dydt(1);
    y[0] = static_cast<MML::Real>(3.0);
    sys.derivs(static_cast<MML::Real>(0.0), y, dydt);
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemFromStdFunc - Operator() call", "[ode_system_stdfunc]")
{
    ODESystemFromStdFunc sys(2, [](Real, const Vector<Real>& y, Vector<Real>& dydt) {
        dydt[0] = y[1];
        dydt[1] = -y[0];
    });

    Vector<Real> y(2), dydt(2);
    y[0] = static_cast<MML::Real>(1.0); y[1] = static_cast<MML::Real>(0.0);
    sys(static_cast<MML::Real>(0.0), y, dydt);
    REQUIRE_THAT(dydt[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(dydt[1], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemWithJacobianFromStdFunc - Jacobian", "[ode_system_stdfunc]")
{
    auto derivFunc = [](Real, const Vector<Real>& y, Vector<Real>& dydt) {
        dydt[0] = y[1];
        dydt[1] = -y[0];
    };
    auto jacFunc = [](const Real, const Vector<Real>&, Vector<Real>& dydt, Matrix<Real>& J) {
        J.resize(2, 2);
        J[0][0] = static_cast<MML::Real>(0.0);  J[0][1] = static_cast<MML::Real>(1.0);
        J[1][0] = -static_cast<MML::Real>(1.0); J[1][1] = static_cast<MML::Real>(0.0);
    };

    ODESystemWithJacobianFromStdFunc sys(2, derivFunc, jacFunc);
    REQUIRE(sys.getDim() == 2);

    Vector<Real> y(2), dydt(2);
    Matrix<Real> J(2, 2);
    y[0] = static_cast<MML::Real>(1.0); y[1] = static_cast<MML::Real>(0.0);
    sys.jacobian(static_cast<MML::Real>(0.0), y, J);
    REQUIRE_THAT(J[0][1], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(J[1][0], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ODESystemWithJacobian - Null jacobian throws NotImplementedError", "[ode_system_jacobian]")
{
    // Default constructor leaves _funcJac as nullptr
    ODESystemWithJacobian sys(2, harmonicOscillator, static_cast<ODESystemWithJacobian::DerivsJacobianFunc>(nullptr));

    Vector<Real> y(2), dydt(2);
    Matrix<Real> J(2, 2);
    y[0] = static_cast<MML::Real>(1.0); y[1] = static_cast<MML::Real>(0.0);

    REQUIRE_THROWS_AS(sys.jacobian(static_cast<MML::Real>(0.0), y, J), NotImplementedError);
    REQUIRE_THROWS_AS(sys.derivsAndJacobian(static_cast<MML::Real>(0.0), y, dydt, J), NotImplementedError);
}

TEST_CASE("ODESystemSolution - size/isEmpty reflect actual data, not capacity", "[ode_system_solution]")
{
    // Create with capacity 500 (internally 501), but no data stored yet
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 2, 500);

    REQUIRE(sol.isEmpty() == true);
    REQUIRE(sol.size() == 0);
    REQUIRE(sol.getNumSteps() == 0);
    REQUIRE(sol.capacity() == 501);  // capacity is maxSteps + 1

    // Store 3 data points using fillValues
    Vector<Real> y(2);
    y[0] = static_cast<MML::Real>(1.0); y[1] = static_cast<MML::Real>(2.0);
    sol.fillValues(0, static_cast<MML::Real>(0.0), y);
    REQUIRE(sol.isEmpty() == false);
    REQUIRE(sol.size() == 1);

    y[0] = static_cast<MML::Real>(3.0); y[1] = static_cast<MML::Real>(4.0);
    sol.fillValues(1, static_cast<MML::Real>(5.0), y);
    REQUIRE(sol.size() == 2);
    REQUIRE(sol.getNumSteps() == 1);  // 1 interval between 2 points

    y[0] = static_cast<MML::Real>(5.0); y[1] = static_cast<MML::Real>(6.0);
    sol.fillValues(2, static_cast<MML::Real>(10.0), y);
    REQUIRE(sol.size() == 3);
    REQUIRE(sol.getNumSteps() == 2);  // 2 intervals between 3 points

    // After setFinalSize, size and capacity should match
    sol.setFinalSize(2);
    REQUIRE(sol.size() == 3);
    REQUIRE(sol.capacity() == 3);
}

TEST_CASE("ODESystemSolution - writes far beyond capacity grow storage enough (no grow-once OOB)", "[ode_system_solution][wewi.1.2]")
{
    // Capacity 11 (maxSteps 10). The old ExtendSavedSteps grew exactly once (to 17), so a write at
    // index 500 landed 480 slots past the end of the storage (analysis v4, BASE D2 [EXEC]).
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0), 2, 10);
    REQUIRE(sol.capacity() == 11);

    sol.setTVal(500, static_cast<MML::Real>(7.5));
    REQUIRE(sol.capacity() >= 501);
    REQUIRE(sol.size() == 501);
    REQUIRE(sol.getTValue(500) == static_cast<MML::Real>(7.5));

    sol.setXVal(5000, 1, static_cast<MML::Real>(-2.0));
    REQUIRE(sol.capacity() >= 5001);
    REQUIRE(sol.size() == 5001);
    REQUIRE(sol.getXValue(5000, 1) == static_cast<MML::Real>(-2.0));
    REQUIRE(sol.getTValue(500) == static_cast<MML::Real>(7.5)); // earlier data preserved across the regrowth

    Vector<Real> y(2);
    y[0] = static_cast<MML::Real>(1.0); y[1] = static_cast<MML::Real>(2.0);
    sol.fillValues(20000, static_cast<MML::Real>(3.0), y);
    REQUIRE(sol.size() == 20001);
    REQUIRE(sol.getTValue(20000) == static_cast<MML::Real>(3.0));
    REQUIRE(sol.getXValue(20000, 0) == static_cast<MML::Real>(1.0));

    // Bulk accessors report the stored points, not the (larger) capacity
    REQUIRE(sol.getTValues().size() == sol.size());
    REQUIRE(sol.getXValues().cols() == sol.size());
    REQUIRE(sol.getXValues().rows() == 2);
}

TEST_CASE("ODESystemSolution - bulk accessors are trimmed to size()", "[ode_system_solution][wewi.1.2]")
{
    ODESystemSolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);   // default capacity 1001
    REQUIRE(sol.getTValues().size() == 0);
    REQUIRE(sol.getXValues().cols() == 0);

    Vector<Real> y({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
    sol.fillValues(0, static_cast<MML::Real>(0.0), y);
    sol.fillValues(1, static_cast<MML::Real>(0.5), y);
    REQUIRE(sol.getTValues().size() == 2);
    REQUIRE(sol.getXValues().cols() == 2);
    REQUIRE(sol.getXValues()[2][1] == static_cast<MML::Real>(3.0));
}

TEST_CASE("ODESystemWithJacobian - is an IODESystemWithJacobian and an IODESystem (single base chain)", "[ode_system_jacobian][wewi.1.1]")
{
    // Before the fix the class derived from ODESystem only, so it could not be passed to any stiff
    // solver (they take IODESystemWithJacobian&) despite carrying a Jacobian (analysis v4, BASE D1).
    static_assert(std::is_base_of_v<IODESystemWithJacobian, ODESystemWithJacobian>);
    static_assert(std::is_convertible_v<ODESystemWithJacobian*, IODESystem*>);               // unambiguous
    static_assert(std::is_convertible_v<ODESystemWithJacobian*, IODESystemWithJacobian*>);
    static_assert(std::is_convertible_v<ODESystemWithJacobianFromStdFunc*, IODESystem*>);    // was ambiguous (diamond)
    static_assert(std::is_convertible_v<ODESystemWithJacobianFromStdFunc*, IODESystemWithJacobian*>);

    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    IODESystemWithJacobian& stiffView = sys;
    IODESystem& plainView = sys;

    Vector<Real> y({static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}), dydt(2);
    Matrix<Real> J(2, 2);
    plainView.derivs(static_cast<MML::Real>(0.0), y, dydt);
    REQUIRE(dydt[0] == static_cast<MML::Real>(0.0));
    REQUIRE(dydt[1] == -static_cast<MML::Real>(1.0));

    stiffView.jacobian(static_cast<MML::Real>(0.0), y, J);
    REQUIRE(J[0][1] == static_cast<MML::Real>(1.0));
    REQUIRE(J[1][0] == -static_cast<MML::Real>(1.0));
    REQUIRE(stiffView.getDim() == 2);
    REQUIRE(plainView.getDim() == 2);

    ODESystemWithJacobian empty;
    REQUIRE_THROWS_AS(empty.derivs(static_cast<MML::Real>(0.0), y, dydt), NotImplementedError);
}

TEST_CASE("ODESystemWithJacobian - accepted by a stiff solver end to end", "[ode_system_jacobian][wewi.1.1]")
{
    ODESystemWithJacobian sys(2, harmonicOscillator, harmonicJacobian);
    Vector<Real> y0({static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)});

    // Backward Euler is only first order and dissipative; a loose check is all that is needed here -
    // the point is that the call compiles and runs with the function-pointer system.
    ODESystemSolution sol = SolveBackwardEuler(sys, static_cast<MML::Real>(0.0), y0, static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.001));
    REQUIRE(sol.size() > 2);
    REQUIRE_THAT(sol.getXValuesAtEnd()[0], WithinAbs(std::cos(static_cast<MML::Real>(1.0)), static_cast<MML::Real>(0.01)));
    REQUIRE_THAT(sol.getXValuesAtEnd()[1], WithinAbs(-std::sin(static_cast<MML::Real>(1.0)), static_cast<MML::Real>(0.01)));
}

TEST_CASE("DAESolution - writes far beyond capacity grow storage enough; setXVal/setYVal update size()", "[dae_solution][wewi.1.2]")
{
    DAESolution sol(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 2, 1, 10);
    REQUIRE(sol.capacity() == 11);

    sol.setTVal(500, static_cast<MML::Real>(1.5));
    REQUIRE(sol.capacity() >= 501);
    REQUIRE(sol.size() == 501);

    sol.setXVal(4000, 1, static_cast<MML::Real>(2.5));
    REQUIRE(sol.size() == 4001);
    REQUIRE(sol.getXValue(4000, 1) == static_cast<MML::Real>(2.5));

    sol.setYVal(6000, 0, static_cast<MML::Real>(3.5));
    REQUIRE(sol.size() == 6001);
    REQUIRE(sol.getYValue(6000, 0) == static_cast<MML::Real>(3.5));
    REQUIRE(sol.getTValue(500) == static_cast<MML::Real>(1.5));

    REQUIRE(sol.getTValues().size() == sol.size());
    REQUIRE(sol.getXValues().cols() == sol.size());
    REQUIRE(sol.getYValues().cols() == sol.size());
    REQUIRE(sol.getYValues().rows() == 1);
}

} // namespace MML::Tests::Base::ODESystemTests
