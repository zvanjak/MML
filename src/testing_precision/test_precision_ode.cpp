///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        test_precision_ode.cpp                                              ///
///  Description: Comprehensive ODE solver precision tests                            ///
///               Compares fixed-step and adaptive methods                            ///
///                                                                                   ///
///  Tests:       - Euler, Midpoint, RK4 (fixed-step)                                 ///
///               - Cash-Karp, Dormand-Prince 5/8 (adaptive)                          ///
///               - Leapfrog/Verlet for Hamiltonian systems                           ///
///               - Energy conservation tests                                         ///
///               - Long-time integration stability                                   ///
///               - Step size vs accuracy tradeoffs                                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "PrecisionTestFramework.h"

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODEAdaptiveIntegrator.h"

using namespace MML;
using namespace MML::PrecisionTesting;

///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST ODE SYSTEMS (local definitions)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Exponential decay: y' = -ky, exact: y(t) = y0*exp(-kt)
 */
class TestExpDecay : public IODESystem
{
    Real _k;
public:
    TestExpDecay(Real k = 1.0) : _k(k) {}
    int getDim() const override { return 1; }
    void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = -_k * y[0];
    }
    Real exact(Real y0, Real t) const { return y0 * std::exp(-_k * t); }
};

/**
 * @brief Simple harmonic oscillator: y'' = -w^2*y
 * State: [x, v], exact: x(t) = A*cos(wt + phi)
 */
class TestSHO : public IODESystem
{
    Real _omega;
public:
    TestSHO(Real omega = 1.0) : _omega(omega) {}
    int getDim() const override { return 2; }
    void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = y[1];           // dx/dt = v
        dydt[1] = -_omega * _omega * y[0];  // dv/dt = -w^2*x
    }
    Real exactX(Real x0, Real v0, Real t) const {
        return x0 * std::cos(_omega * t) + (v0 / _omega) * std::sin(_omega * t);
    }
    Real exactV(Real x0, Real v0, Real t) const {
        return -_omega * x0 * std::sin(_omega * t) + v0 * std::cos(_omega * t);
    }
    Real energy(Real x, Real v) const {
        return 0.5 * v * v + 0.5 * _omega * _omega * x * x;  // KE + PE
    }
    Real omega() const { return _omega; }
};

/**
 * @brief Damped harmonic oscillator: y'' + 2*zeta*w*y' + w^2*y = 0
 */
class TestDampedSHO : public IODESystem
{
    Real _omega, _zeta;
public:
    TestDampedSHO(Real omega = 1.0, Real zeta = 0.1) : _omega(omega), _zeta(zeta) {}
    int getDim() const override { return 2; }
    void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = y[1];
        dydt[1] = -2.0 * _zeta * _omega * y[1] - _omega * _omega * y[0];
    }
    // Underdamped solution (zeta < 1)
    Real exactX(Real x0, Real v0, Real t) const {
        Real omega_d = _omega * std::sqrt(1.0 - _zeta * _zeta);
        Real exp_term = std::exp(-_zeta * _omega * t);
        Real A = x0;
        Real B = (v0 + _zeta * _omega * x0) / omega_d;
        return exp_term * (A * std::cos(omega_d * t) + B * std::sin(omega_d * t));
    }
};

/**
 * @brief Logistic growth: y' = r*y*(1 - y/K)
 * Exact: y(t) = K / (1 + (K/y0 - 1)*exp(-r*t))
 */
class TestLogistic : public IODESystem
{
    Real _r, _K;
public:
    TestLogistic(Real r = 1.0, Real K = 10.0) : _r(r), _K(K) {}
    int getDim() const override { return 1; }
    void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = _r * y[0] * (1.0 - y[0] / _K);
    }
    Real exact(Real y0, Real t) const {
        return _K / (1.0 + (_K / y0 - 1.0) * std::exp(-_r * t));
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
//                     FIXED-STEP SOLVER COMPARISON
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compare Euler, Midpoint, RK4 on exponential decay
 */
void Test_ODE_FixedStep_ExpDecay()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Fixed-Step ODE Solvers\n";
    std::cout << "  Exponential decay y' = -y with varying step counts\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Fixed-Step ODE Solvers", "Exponential decay comparison");
    
    TestExpDecay system(1.0);  // k = 1
    Real y0 = 1.0;
    Real t_end = 5.0;
    Real exact = system.exact(y0, t_end);
    
    EulerStep_Calculator euler;
    Midpoint_StepCalculator midpoint;
    RungeKutta4_StepCalculator rk4;
    
    std::vector<int> stepCounts = {10, 20, 50, 100, 200, 500, 1000, 2000};
    
    for (int n : stepCounts)
    {
        std::string nStr = "n=" + std::to_string(n);
        
        // Euler
        {
            ODESystemFixedStepSolver solver(system, euler);
            auto sol = solver.integrate(Vector<Real>{y0}, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("Euler", nStr, exact, computed);
        }
        
        // Midpoint
        {
            ODESystemFixedStepSolver solver(system, midpoint);
            auto sol = solver.integrate(Vector<Real>{y0}, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("Midpoint", nStr, exact, computed);
        }
        
        // RK4
        {
            ODESystemFixedStepSolver solver(system, rk4);
            auto sol = solver.integrate(Vector<Real>{y0}, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("RK4", nStr, exact, computed);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Compare step calculators on harmonic oscillator (1 period)
 */
void Test_ODE_FixedStep_SHO()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Fixed-Step ODE - Harmonic Oscillator\n";
    std::cout << "  One period integration, comparing position accuracy\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Fixed-Step SHO", "One period integration");
    
    TestSHO system(1.0);  // omega = 1
    Real x0 = 1.0, v0 = 0.0;
    Real t_end = 2.0 * Constants::PI;  // One period
    Real exact_x = system.exactX(x0, v0, t_end);
    
    EulerStep_Calculator euler;
    Midpoint_StepCalculator midpoint;
    RungeKutta4_StepCalculator rk4;
    Leapfrog_StepCalculator leapfrog;
    VelocityVerlet_StepCalculator verlet;
    
    std::vector<int> stepCounts = {50, 100, 200, 500, 1000, 2000};
    
    for (int n : stepCounts)
    {
        std::string nStr = "n=" + std::to_string(n);
        Vector<Real> ic{x0, v0};
        
        // Euler
        {
            ODESystemFixedStepSolver solver(system, euler);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("Euler", nStr, exact_x, computed);
        }
        
        // Midpoint
        {
            ODESystemFixedStepSolver solver(system, midpoint);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("Midpoint", nStr, exact_x, computed);
        }
        
        // RK4
        {
            ODESystemFixedStepSolver solver(system, rk4);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("RK4", nStr, exact_x, computed);
        }
        
        // Leapfrog
        {
            ODESystemFixedStepSolver solver(system, leapfrog);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("Leapfrog", nStr, exact_x, computed);
        }
        
        // Velocity Verlet
        {
            ODESystemFixedStepSolver solver(system, verlet);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("Verlet", nStr, exact_x, computed);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     ADAPTIVE SOLVER TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compare adaptive solvers (Cash-Karp, DP5) on exponential decay
 */
void Test_ODE_Adaptive_ExpDecay()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Adaptive ODE Solvers\n";
    std::cout << "  Exponential decay with varying tolerances\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Adaptive ODE Solvers", "Tolerance vs accuracy");
    
    TestExpDecay system(1.0);
    Real y0 = 1.0;
    Real t_end = 5.0;
    Real exact = system.exact(y0, t_end);
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        // Cash-Karp
        {
            ODEAdaptiveIntegrator<CashKarp_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("CashKarp", tolStr.str(), exact, computed);
        }
        
        // Dormand-Prince 5
        {
            ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP5", tolStr.str(), exact, computed);
        }
        
        // Dormand-Prince 8
        {
            ODEAdaptiveIntegrator<DormandPrince8_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP8", tolStr.str(), exact, computed);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Adaptive solvers on harmonic oscillator
 */
void Test_ODE_Adaptive_SHO()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Adaptive ODE - Harmonic Oscillator\n";
    std::cout << "  10 periods with varying tolerances\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Adaptive SHO", "Long-time oscillator");
    
    TestSHO system(1.0);
    Real x0 = 1.0, v0 = 0.0;
    Real t_end = 20.0 * Constants::PI;  // 10 periods
    Real exact_x = system.exactX(x0, v0, t_end);
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        // Cash-Karp
        {
            ODEAdaptiveIntegrator<CashKarp_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("CashKarp", tolStr.str(), exact_x, computed);
        }
        
        // Dormand-Prince 5
        {
            ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP5", tolStr.str(), exact_x, computed);
        }
        
        // Dormand-Prince 8
        {
            ODEAdaptiveIntegrator<DormandPrince8_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP8", tolStr.str(), exact_x, computed);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     ENERGY CONSERVATION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test energy conservation for symplectic vs non-symplectic methods
 */
void Test_ODE_EnergyConservation()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Energy Conservation\n";
    std::cout << "  Harmonic oscillator - energy drift over 100 periods\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Energy Conservation", "100 period integration");
    
    TestSHO system(1.0);
    Real x0 = 1.0, v0 = 0.0;
    Real E0 = system.energy(x0, v0);  // Initial energy
    Real t_end = 200.0 * Constants::PI;  // 100 periods
    int n_steps = 10000;  // Fixed steps for fair comparison
    
    EulerStep_Calculator euler;
    RungeKutta4_StepCalculator rk4;
    Leapfrog_StepCalculator leapfrog;
    VelocityVerlet_StepCalculator verlet;
    
    Vector<Real> ic{x0, v0};
    
    // Euler (expected: energy grows exponentially)
    {
        ODESystemFixedStepSolver solver(system, euler);
        auto sol = solver.integrate(ic, 0.0, t_end, n_steps);
        auto final_state = sol.getXValuesAtEnd();
        Real E_final = system.energy(final_state[0], final_state[1]);
        suite.addResult("Euler", "E_100T", E0, E_final);
    }
    
    // RK4 (expected: small secular drift)
    {
        ODESystemFixedStepSolver solver(system, rk4);
        auto sol = solver.integrate(ic, 0.0, t_end, n_steps);
        auto final_state = sol.getXValuesAtEnd();
        Real E_final = system.energy(final_state[0], final_state[1]);
        suite.addResult("RK4", "E_100T", E0, E_final);
    }
    
    // Leapfrog (expected: bounded energy oscillations)
    {
        ODESystemFixedStepSolver solver(system, leapfrog);
        auto sol = solver.integrate(ic, 0.0, t_end, n_steps);
        auto final_state = sol.getXValuesAtEnd();
        Real E_final = system.energy(final_state[0], final_state[1]);
        suite.addResult("Leapfrog", "E_100T", E0, E_final);
    }
    
    // Velocity Verlet (expected: bounded energy, symplectic)
    {
        ODESystemFixedStepSolver solver(system, verlet);
        auto sol = solver.integrate(ic, 0.0, t_end, n_steps);
        auto final_state = sol.getXValuesAtEnd();
        Real E_final = system.energy(final_state[0], final_state[1]);
        suite.addResult("Verlet", "E_100T", E0, E_final);
    }
    
    suite.printDetailedTable();
    suite.printSummary();
    
    std::cout << "\n  NOTE: Symplectic methods (Leapfrog, Verlet) should preserve energy better\n";
    std::cout << "        than non-symplectic methods (Euler, RK4) over long integration.\n\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     CONVERGENCE ORDER VERIFICATION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Verify theoretical convergence orders of methods
 */
void Test_ODE_ConvergenceOrder()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Convergence Order Verification\n";
    std::cout << "  Expected: Euler=1, Midpoint=2, RK4=4\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Convergence Order", "Error reduction with step halving");
    
    TestExpDecay system(1.0);
    Real y0 = 1.0;
    Real t_end = 1.0;  // Short integration for clean order verification
    Real exact = system.exact(y0, t_end);
    
    EulerStep_Calculator euler;
    Midpoint_StepCalculator midpoint;
    RungeKutta4_StepCalculator rk4;
    
    // Powers of 2 for clean halving
    std::vector<int> stepCounts = {16, 32, 64, 128, 256, 512, 1024};
    
    std::cout << "Step count    Euler Error      Midpoint Error   RK4 Error\n";
    std::cout << "------------------------------------------------------------\n";
    
    for (int n : stepCounts)
    {
        std::string nStr = "n=" + std::to_string(n);
        Vector<Real> ic{y0};
        
        Real euler_err, midpoint_err, rk4_err;
        
        {
            ODESystemFixedStepSolver solver(system, euler);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            euler_err = std::abs(computed - exact);
            suite.addResult("Euler", nStr, exact, computed);
        }
        
        {
            ODESystemFixedStepSolver solver(system, midpoint);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            midpoint_err = std::abs(computed - exact);
            suite.addResult("Midpoint", nStr, exact, computed);
        }
        
        {
            ODESystemFixedStepSolver solver(system, rk4);
            auto sol = solver.integrate(ic, 0.0, t_end, n);
            Real computed = sol.getXValuesAtEnd()[0];
            rk4_err = std::abs(computed - exact);
            suite.addResult("RK4", nStr, exact, computed);
        }
        
        std::cout << std::setw(8) << n << "    "
                  << std::scientific << std::setprecision(3) << euler_err << "     "
                  << midpoint_err << "     "
                  << rk4_err << "\n";
    }
    
    std::cout << "\n";
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    std::cout << "\n  Order verification: Error should decrease by factor of:\n";
    std::cout << "    - Euler:    2^1 = 2   (halving step -> halving error)\n";
    std::cout << "    - Midpoint: 2^2 = 4   (halving step -> quarter error)\n";
    std::cout << "    - RK4:      2^4 = 16  (halving step -> 1/16 error)\n\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     DAMPED OSCILLATOR (DECAY + OSCILLATION)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test on damped harmonic oscillator
 */
void Test_ODE_DampedOscillator()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Damped Harmonic Oscillator\n";
    std::cout << "  Combined decay and oscillation behavior\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Damped Oscillator", "w=1, zeta=0.1");
    
    TestDampedSHO system(1.0, 0.1);  // omega=1, zeta=0.1 (underdamped)
    Real x0 = 1.0, v0 = 0.0;
    Real t_end = 20.0;  // Several damped oscillations
    Real exact = system.exactX(x0, v0, t_end);
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        // Cash-Karp
        {
            ODEAdaptiveIntegrator<CashKarp_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("CashKarp", tolStr.str(), exact, computed);
        }
        
        // Dormand-Prince 5
        {
            ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP5", tolStr.str(), exact, computed);
        }
        
        // Dormand-Prince 8
        {
            ODEAdaptiveIntegrator<DormandPrince8_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP8", tolStr.str(), exact, computed);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     LOGISTIC EQUATION (NONLINEAR)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test on logistic growth equation
 */
void Test_ODE_LogisticGrowth()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Logistic Growth Equation\n";
    std::cout << "  Nonlinear ODE with exact solution\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Logistic Growth", "r=1, K=10");
    
    TestLogistic system(1.0, 10.0);  // r=1, K=10
    Real y0 = 0.5;  // Start below carrying capacity
    Real t_end = 10.0;  // Long enough to approach K
    Real exact = system.exact(y0, t_end);
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        // Cash-Karp
        {
            ODEAdaptiveIntegrator<CashKarp_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("CashKarp", tolStr.str(), exact, computed);
        }
        
        // Dormand-Prince 5
        {
            ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP5", tolStr.str(), exact, computed);
        }
        
        // Dormand-Prince 8
        {
            ODEAdaptiveIntegrator<DormandPrince8_Stepper> integrator(system);
            auto sol = integrator.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
            Real computed = sol.getXValuesAtEnd()[0];
            suite.addResult("DP8", tolStr.str(), exact, computed);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     COMPREHENSIVE SUMMARY
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Summary comparison across all test problems
 */
void Test_ODE_ComprehensiveSummary()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "  COMPREHENSIVE ODE SOLVER ACCURACY SUMMARY\n";
    std::cout << "================================================================================\n\n";
    
    PrecisionTestSuite suite("All ODE Solvers", "Complete error comparison");
    
    // Test problems
    TestExpDecay expDecay(1.0);
    TestSHO sho(1.0);
    TestLogistic logistic(1.0, 10.0);
    
    Real tol = 1e-8;  // Standard tolerance for comparison
    
    // Exponential decay
    {
        Real y0 = 1.0, t_end = 5.0;
        Real exact = expDecay.exact(y0, t_end);
        
        ODEAdaptiveIntegrator<CashKarp_Stepper> ck(expDecay);
        ODEAdaptiveIntegrator<DormandPrince5_Stepper> dp5(expDecay);
        ODEAdaptiveIntegrator<DormandPrince8_Stepper> dp8(expDecay);
        
        auto sol_ck = ck.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
        auto sol_dp5 = dp5.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
        auto sol_dp8 = dp8.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
        
        suite.addResult("CashKarp", "ExpDecay", exact, sol_ck.getXValuesAtEnd()[0]);
        suite.addResult("DP5", "ExpDecay", exact, sol_dp5.getXValuesAtEnd()[0]);
        suite.addResult("DP8", "ExpDecay", exact, sol_dp8.getXValuesAtEnd()[0]);
    }
    
    // Harmonic oscillator (1 period)
    {
        Real x0 = 1.0, v0 = 0.0, t_end = 2.0 * Constants::PI;
        Real exact = sho.exactX(x0, v0, t_end);
        
        ODEAdaptiveIntegrator<CashKarp_Stepper> ck(sho);
        ODEAdaptiveIntegrator<DormandPrince5_Stepper> dp5(sho);
        ODEAdaptiveIntegrator<DormandPrince8_Stepper> dp8(sho);
        
        auto sol_ck = ck.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
        auto sol_dp5 = dp5.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
        auto sol_dp8 = dp8.integrate(Vector<Real>{x0, v0}, 0.0, t_end, 0.1, tol);
        
        suite.addResult("CashKarp", "SHO", exact, sol_ck.getXValuesAtEnd()[0]);
        suite.addResult("DP5", "SHO", exact, sol_dp5.getXValuesAtEnd()[0]);
        suite.addResult("DP8", "SHO", exact, sol_dp8.getXValuesAtEnd()[0]);
    }
    
    // Logistic growth
    {
        Real y0 = 0.5, t_end = 10.0;
        Real exact = logistic.exact(y0, t_end);
        
        ODEAdaptiveIntegrator<CashKarp_Stepper> ck(logistic);
        ODEAdaptiveIntegrator<DormandPrince5_Stepper> dp5(logistic);
        ODEAdaptiveIntegrator<DormandPrince8_Stepper> dp8(logistic);
        
        auto sol_ck = ck.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
        auto sol_dp5 = dp5.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
        auto sol_dp8 = dp8.integrate(Vector<Real>{y0}, 0.0, t_end, 0.1, tol);
        
        suite.addResult("CashKarp", "Logistic", exact, sol_ck.getXValuesAtEnd()[0]);
        suite.addResult("DP5", "Logistic", exact, sol_dp5.getXValuesAtEnd()[0]);
        suite.addResult("DP8", "Logistic", exact, sol_dp8.getXValuesAtEnd()[0]);
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Export results
    suite.exportCSV("results/ode_precision.csv");
    suite.exportMarkdown("results/ode_precision.md");
    
    std::cout << "\n  Results exported to:\n";
    std::cout << "    - results/ode_precision.csv\n";
    std::cout << "    - results/ode_precision.md\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     MASTER TEST FUNCTION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Run all ODE precision tests
 */
void Test_Precision_ODE()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "               ODE SOLVER PRECISION TEST SUITE\n";
    std::cout << "================================================================================\n";
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  1. FIXED-STEP METHOD COMPARISON\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ODE_FixedStep_ExpDecay();
    Test_ODE_FixedStep_SHO();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  2. ADAPTIVE SOLVER COMPARISON\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ODE_Adaptive_ExpDecay();
    Test_ODE_Adaptive_SHO();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  3. ENERGY CONSERVATION (Symplectic vs Non-Symplectic)\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ODE_EnergyConservation();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  4. CONVERGENCE ORDER VERIFICATION\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ODE_ConvergenceOrder();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  5. ADDITIONAL TEST PROBLEMS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ODE_DampedOscillator();
    Test_ODE_LogisticGrowth();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  6. COMPREHENSIVE SUMMARY\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ODE_ComprehensiveSummary();
    
    std::cout << "\n================================================================================\n";
    std::cout << "            ODE PRECISION TESTS COMPLETE\n";
    std::cout << "================================================================================\n";
}
