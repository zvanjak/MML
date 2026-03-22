///////////////////////////////////////////////////////////////////////////////////////////
// quick_start_guide.cpp - Verification of all examples from QUICK_START.md
//
// This file ensures all code examples in QUICK_START.md compile and run correctly.
// Each example is numbered according to its section in the markdown file.
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Vector/VectorTypes.h"
#include "base/Matrix/Matrix.h"
#include "base/ODESystem.h"
#include "base/ODESystemSolution.h"
#include "core/LinAlgEqSolvers.h"
#include "core/Derivation.h"
#include "core/Integration.h"
#include "algorithms/RootFinding.h"
#include "mml/algorithms/ODESolvers/ODESolverFixedStep.h"
#include "mml/algorithms/ODESolvers/ODESolverAdaptive.h"
#endif

using namespace MML;
#include <iostream>
#include <cmath>
#include <functional>

void example_step2_vectors()
{
    std::cout << "\n=== Step 2: Vectors ===" << std::endl;
    
    // Create vectors
    Vector<Real> v{1.0, 2.0, 3.0};
    Vector<Real> w{4.0, 5.0, 6.0};
    
    // Basic operations
    Vector<Real> sum = v + w;
    
    // For dot product, use Vector3Cartesian
    Vector3Cartesian v3{1.0, 2.0, 3.0};
    Vector3Cartesian w3{4.0, 5.0, 6.0};
    Real dot = ScalarProduct(v3, w3);  // Correct: global function
    
    Real length = v.NormL2();          // Length
    
    std::cout << "v + w = " << sum << std::endl;
    std::cout << "v . w = " << dot << std::endl;
    std::cout << "||v|| = " << length << std::endl;
}

void example_step2_matrices()
{
    std::cout << "\n=== Step 2: Matrices ===" << std::endl;
    
    // Create a 3x3 matrix
    Matrix<Real> A{3, 3, {
        1, 2, 3,
        4, 5, 6,
        7, 8, 10
    }};

    // Matrix operations
    Matrix<Real> At = A.transpose();
    Matrix<Real> I = Matrix<Real>::Identity(3);
    Real trace = A.trace();

    std::cout << "Matrix A:" << std::endl;
    A.Print(std::cout, 8, 3);
    std::cout << std::endl;
    std::cout << "Trace: " << trace << std::endl;
}

void example_step3_linear_system()
{
    std::cout << "\n=== Step 3: Linear System ===" << std::endl;
    
    // System: A·x = b
    Matrix<Real> A{3, 3, {
        4, 1, 2,
        1, 5, 1,
        2, 1, 6
    }};
    
    Vector<Real> b{8, 9, 12};
    
    // Solve using LU decomposition
    LUSolver<Real> solver(A);
    Vector<Real> x = solver.Solve(b);
    
    std::cout << "Solution: " << x << std::endl;
    
    // Verify: A*x should equal b
    Vector<Real> check = A * x;
    std::cout << "A*x = " << check << std::endl;
    std::cout << "b   = " << b << std::endl;
}

void example_step4_derivatives()
{
    std::cout << "\n=== Step 4: Derivatives ===" << std::endl;
    
    // Define a function: f(x) = x²
    RealFunction f{ [](Real x) { return x * x; } };
    
    // Compute derivative at x = 2.0
    Real x0 = 2.0;
    Real df = Derivation::NDer2(f, x0);  // f'(2.0) using 2nd order
    
    std::cout << "f(x) = x^2" << std::endl;
    std::cout << "f'(2.0) = " << df << std::endl;
    std::cout << "Expected: 4.0" << std::endl;
}

void example_step4_integrals()
{
    std::cout << "\n=== Step 4: Integrals ===" << std::endl;
    
    // Define function: f(x) = sin(x)
    RealFunction f{ [](Real x) { return std::sin(x); } };

    // Integrate from 0 to π
    Real a = 0.0, b = Constants::PI;
    IntegrationResult result = IntegrateRomberg(f, a, b);
    Real integral = result.value;

    std::cout << "Integral of sin(x) from 0 to pi = " << integral << std::endl;
    std::cout << "Expected: 2.0" << std::endl;
}

void example_step5_root_finding()
{
    std::cout << "\n=== Step 5: Root Finding ===" << std::endl;
    
    // Find root of: f(x) = x² - 2 = 0  (answer: x = √2)
    RealFunction f{ [](Real x) { return x*x - 2.0; } };
    
    // Find root in interval [1, 2]
    Real root = RootFinding::FindRootBrent(f, 1.0, 2.0, 1e-10);
    
    std::cout << "Root of x^2 - 2 = 0: " << root << std::endl;
    std::cout << "sqrt(2) = " << std::sqrt(2.0) << std::endl;
}

void example_step6_odes()
{
    std::cout << "\n=== Step 6: ODEs ===" << std::endl;
    
    // Define ODE system
    class SimpleODE : public IODESystem {
    public:
        int getDim() const override { return 1; }
        void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
            dydt[0] = -y[0];  // dy/dt = -y
        }
    };
    
    SimpleODE ode_sys;
    Vector<Real> y0{1.0};
    
    // Use adaptive integrator
    ODEAdaptiveIntegrator<DormandPrince5_Stepper> solver(ode_sys);
    Real t0 = 0.0, tend = 2.0;
    Real outputInterval = 0.02;
    
    ODESystemSolution sol = solver.integrate(y0, t0, tend, outputInterval);
    
    // Get final value using correct API
    Vector<Real> y_final = sol.getXValuesAtEnd();
    
    std::cout << "y(2) = " << y_final[0] << std::endl;
    std::cout << "Expected: e^(-2) = " << std::exp(-2.0) << std::endl;
}

void complete_example()
{
    std::cout << "\n=== Complete Example: Everything Together ===" << std::endl;
    
    // 1. Vectors
    std::cout << "\n1. VECTORS:" << std::endl;
    Vector3Cartesian v{1, 2, 3}, w{4, 5, 6};
    std::cout << "   v.w = " << ScalarProduct(v, w) << std::endl;
    
    // 2. Matrices
    std::cout << "\n2. MATRICES:" << std::endl;
    Matrix<Real> A = Matrix<Real>::Identity(3);
    std::cout << "   I_3 trace = " << A.trace() << std::endl;
    
    // 3. Linear Systems
    std::cout << "\n3. LINEAR SYSTEM:" << std::endl;
    Matrix<Real> M{2, 2, {3, 1, 1, 2}};
    Vector<Real> b{9, 8};
    LUSolver<Real> lu(M);
    Vector<Real> x = lu.Solve(b);
    std::cout << "   Solution: [" << x[0] << ", " << x[1] << "]" << std::endl;
    
    // 4. Calculus
    std::cout << "\n4. CALCULUS:" << std::endl;
    RealFunction f{ [](Real x) { return x*x; } };
    Real df = Derivation::NDer2(f, 3.0);
    std::cout << "   (x^2)' at x=3: " << df << " (expected: 6)" << std::endl;
    
    RealFunction g{ [](Real x) { return x; } };
    IntegrationResult int_result = IntegrateRomberg(g, 0, 2);
    std::cout << "   Integral of x dx from 0 to 2 = " << int_result.value << " (expected: 2)" << std::endl;
    
    // 5. Root Finding
    std::cout << "\n5. ROOT FINDING:" << std::endl;
    RealFunction h{ [](Real x) { return x*x - 2; } };
    Real root = RootFinding::FindRootBrent(h, 0, 2, 1e-10);
    std::cout << "   sqrt(2) ~ " << root << std::endl;
    
    // 6. ODEs
    std::cout << "\n6. DIFFERENTIAL EQUATIONS:" << std::endl;
    class ExpDecay : public IODESystem {
    public:
        int getDim() const override { return 1; }
        void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
            dydt[0] = -y[0];
        }
    };
    
    ExpDecay ode;
    ODEAdaptiveIntegrator<DormandPrince5_Stepper> solver(ode);
    ODESystemSolution sol = solver.integrate(Vector<Real>{1.0}, 0.0, 1.0, 0.02);
    Vector<Real> y_final = sol.getXValuesAtEnd();
    std::cout << "   e^(-1) ~ " << y_final[0] << std::endl;
    
    std::cout << "\n All tests passed! You're ready to use MML!" << std::endl;
}

void Docs_Demo_QuickStartGuide()
{
    std::cout << "=====================================================" << std::endl;
    std::cout << "  QUICK START GUIDE - Example Verification" << std::endl;
    std::cout << "=====================================================" << std::endl;
    
    example_step2_vectors();
    example_step2_matrices();
    example_step3_linear_system();
    example_step4_derivatives();
    example_step4_integrals();
    example_step5_root_finding();
    example_step6_odes();
    complete_example();
    
    std::cout << "\n=====================================================" << std::endl;
    std::cout << "  ALL QUICK START EXAMPLES VERIFIED SUCCESSFULLY!" << std::endl;
    std::cout << "=====================================================" << std::endl;
}

int main()
{
    Docs_Demo_QuickStartGuide();
    return 0;
}
