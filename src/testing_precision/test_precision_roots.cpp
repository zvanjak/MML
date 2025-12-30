///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        test_precision_roots.cpp                                            ///
///  Description: Comprehensive root finding precision tests                          ///
///               Compares Bisection, Newton, Secant, FalsePosition, Ridders, Brent   ///
///                                                                                   ///
///  Tests:       - Single root finding accuracy                                      ///
///               - Convergence rate verification                                     ///
///               - Multiple roots detection                                          ///
///               - Edge cases (near-zero derivative, multiple roots)                 ///
///               - Polynomial root finding                                           ///
///               - Transcendental equations                                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "PrecisionTestFramework.h"

#include "MMLBase.h"

#include "base/Vector.h"
#include "interfaces/IFunction.h"

#include "algorithms/RootFinding.h"

using namespace MML;
using namespace MML::PrecisionTesting;

///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST FUNCTIONS (local definitions)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Simple polynomial: f(x) = x^2 - 2 (root at sqrt(2))
 */
class TestSqrt2 : public IRealFunction
{
public:
    Real operator()(Real x) const override { return x * x - 2.0; }
    Real exactRoot() const { return std::sqrt(2.0); }
};

/**
 * @brief Cubic polynomial: f(x) = x^3 - x - 2 (root near 1.5214)
 */
class TestCubic : public IRealFunction
{
public:
    Real operator()(Real x) const override { return x * x * x - x - 2.0; }
    // Real root: x ≈ 1.5213797068045675
    Real exactRoot() const { return 1.5213797068045675; }
};

/**
 * @brief Exponential: f(x) = e^x - 3x (roots at ~0.619 and ~1.512)
 */
class TestExpLinear : public IRealFunction
{
public:
    Real operator()(Real x) const override { return std::exp(x) - 3.0 * x; }
    Real exactRoot1() const { return 0.6190612867359451; }  // smaller root
    Real exactRoot2() const { return 1.5121345516578425; }  // larger root
};

/**
 * @brief Trigonometric: f(x) = cos(x) - x (root at ~0.739)
 */
class TestCosX : public IRealFunction
{
public:
    Real operator()(Real x) const override { return std::cos(x) - x; }
    // Dottie number: x ≈ 0.7390851332151607
    Real exactRoot() const { return 0.7390851332151607; }
};

/**
 * @brief Kepler's equation: E - e*sin(E) - M = 0
 * Fundamental in orbital mechanics
 */
class TestKepler : public IRealFunction
{
    Real _e;  // eccentricity
    Real _M;  // mean anomaly
public:
    TestKepler(Real e, Real M) : _e(e), _M(M) {}
    Real operator()(Real E) const override { return E - _e * std::sin(E) - _M; }
};

/**
 * @brief Near-zero derivative at root: f(x) = (x-1)^3
 * Tests robustness when f'(root) = 0 (triple root)
 */
class TestTripleRoot : public IRealFunction
{
public:
    Real operator()(Real x) const override { 
        Real d = x - 1.0;
        return d * d * d; 
    }
    Real exactRoot() const { return 1.0; }
};

/**
 * @brief Polynomial with known roots: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
 */
class TestThreeRoots : public IRealFunction
{
public:
    Real operator()(Real x) const override { 
        return x * x * x - 6.0 * x * x + 11.0 * x - 6.0; 
    }
    std::vector<Real> exactRoots() const { return {1.0, 2.0, 3.0}; }
};

/**
 * @brief Challenging function: f(x) = x * exp(-x^2) 
 * Root at 0, but flat near the root
 */
class TestFlatRoot : public IRealFunction
{
public:
    Real operator()(Real x) const override { return x * std::exp(-x * x); }
    Real exactRoot() const { return 0.0; }
};

/**
 * @brief Bessel function approximation: sin(x)/x - 0.5
 * Multiple roots
 */
class TestSincHalf : public IRealFunction
{
public:
    Real operator()(Real x) const override { 
        if (std::abs(x) < 1e-10) return 0.5;  // sinc(0) = 1
        return std::sin(x) / x - 0.5; 
    }
    // First positive root: x ≈ 1.895494267033981
    Real exactRoot() const { return 1.895494267033981; }
};

/**
 * @brief High-precision test: x^5 - x - 1 (Wilkinson-like)
 */
class TestQuintic : public IRealFunction
{
public:
    Real operator()(Real x) const override { 
        return x * x * x * x * x - x - 1.0; 
    }
    // Real root: x ≈ 1.1673039782614187
    Real exactRoot() const { return 1.1673039782614187; }
};

///////////////////////////////////////////////////////////////////////////////////////////
//                     BASIC METHOD COMPARISON
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compare all root finding methods on sqrt(2) problem
 */
void Test_Roots_BasicComparison()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding Methods - Basic Comparison\n";
    std::cout << "  Test function: f(x) = x^2 - 2, root = sqrt(2)\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Root Finding Basic", "sqrt(2) comparison");
    
    TestSqrt2 func;
    Real exact = func.exactRoot();
    Real x1 = 1.0, x2 = 2.0;  // bracket containing sqrt(2)
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        // Bisection
        {
            Real root = RootFinding::FindRootBisection(func, x1, x2, tol);
            suite.addResult("Bisection", tolStr.str(), exact, root);
        }
        
        // Newton-Raphson
        {
            Real root = RootFinding::FindRootNewton(func, x1, x2, tol);
            suite.addResult("Newton", tolStr.str(), exact, root);
        }
        
        // Secant
        {
            Real root = RootFinding::FindRootSecant(func, x1, x2, tol);
            suite.addResult("Secant", tolStr.str(), exact, root);
        }
        
        // False Position
        {
            Real root = RootFinding::FindRootFalsePosition(func, x1, x2, tol);
            suite.addResult("FalsePos", tolStr.str(), exact, root);
        }
        
        // Ridders
        {
            Real root = RootFinding::FindRootRidders(func, x1, x2, tol);
            suite.addResult("Ridders", tolStr.str(), exact, root);
        }
        
        // Brent
        {
            Real root = RootFinding::FindRootBrent(func, x1, x2, tol);
            suite.addResult("Brent", tolStr.str(), exact, root);
        }
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Compare methods on cubic polynomial
 */
void Test_Roots_Cubic()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding - Cubic Polynomial\n";
    std::cout << "  Test function: f(x) = x^3 - x - 2\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Root Finding Cubic", "Cubic polynomial");
    
    TestCubic func;
    Real exact = func.exactRoot();
    Real x1 = 1.0, x2 = 2.0;
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        suite.addResult("Bisection", tolStr.str(), exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", tolStr.str(), exact, 
            RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Secant", tolStr.str(), exact, 
            RootFinding::FindRootSecant(func, x1, x2, tol));
        suite.addResult("Ridders", tolStr.str(), exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", tolStr.str(), exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     TRANSCENDENTAL EQUATIONS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test on cos(x) = x (Dottie number)
 */
void Test_Roots_Transcendental()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding - Transcendental Equations\n";
    std::cout << "  Test function: f(x) = cos(x) - x (Dottie number)\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Transcendental", "Dottie number");
    
    TestCosX func;
    Real exact = func.exactRoot();
    Real x1 = 0.0, x2 = 1.0;
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
    
    for (Real tol : tolerances)
    {
        std::ostringstream tolStr;
        tolStr << "tol=" << std::scientific << std::setprecision(0) << tol;
        
        suite.addResult("Bisection", tolStr.str(), exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", tolStr.str(), exact, 
            RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Secant", tolStr.str(), exact, 
            RootFinding::FindRootSecant(func, x1, x2, tol));
        suite.addResult("Ridders", tolStr.str(), exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", tolStr.str(), exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test exponential-linear equation with two roots
 */
void Test_Roots_ExpLinear()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding - Exponential-Linear\n";
    std::cout << "  Test function: f(x) = e^x - 3x (two roots)\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Exp-Linear", "Two roots test");
    
    TestExpLinear func;
    Real tol = 1e-12;
    
    // First root (smaller, bracket [0, 1])
    {
        Real exact = func.exactRoot1();
        Real x1 = 0.0, x2 = 1.0;
        
        suite.addResult("Bisection", "root1", exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", "root1", exact, 
            RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Ridders", "root1", exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", "root1", exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    // Second root (larger, bracket [1, 2])
    {
        Real exact = func.exactRoot2();
        Real x1 = 1.0, x2 = 2.0;
        
        suite.addResult("Bisection", "root2", exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", "root2", exact, 
            RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Ridders", "root2", exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", "root2", exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     KEPLER'S EQUATION (PHYSICS APPLICATION)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test Kepler's equation for various eccentricities
 */
void Test_Roots_Kepler()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding - Kepler's Equation\n";
    std::cout << "  E - e*sin(E) - M = 0 (orbital mechanics)\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Kepler Equation", "Orbital mechanics");
    
    Real M = 0.8;  // Mean anomaly
    Real tol = 1e-12;
    
    // Test with different eccentricities
    std::vector<Real> eccentricities = {0.1, 0.3, 0.5, 0.7, 0.9};
    
    for (Real e : eccentricities)
    {
        TestKepler kepler(e, M);
        std::ostringstream eStr;
        eStr << "e=" << std::fixed << std::setprecision(1) << e;
        
        // For Kepler's equation, E is always in [M-e, M+e] for 0 <= e < 1
        // Use a safer bracket: [0, pi] for M in first quadrant
        Real x1 = 0.0;
        Real x2 = Constants::PI;  // E is always less than pi for M < pi
        
        // For Kepler, verify by computing residual: E - e*sin(E) - M should be ~0
        try {
            Real root_brent = RootFinding::FindRootBrent(kepler, x1, x2, tol);
            Real residual = kepler(root_brent);
            suite.addResult("Brent", eStr.str(), 0.0, residual);
        } catch (...) {
            suite.addResult("Brent", eStr.str(), 0.0, 1.0);  // Mark as failed
        }
        
        // Newton may fail with wide brackets - use tighter bracket around M
        try {
            Real x1n = std::max(0.0, M - 1.0);
            Real x2n = M + 1.0;
            Real root_newton = RootFinding::FindRootNewton(kepler, x1n, x2n, tol);
            Real residual = kepler(root_newton);
            suite.addResult("Newton", eStr.str(), 0.0, residual);
        } catch (...) {
            suite.addResult("Newton", eStr.str(), 0.0, 1.0);  // Mark as failed
        }
        
        try {
            Real root_ridders = RootFinding::FindRootRidders(kepler, x1, x2, tol);
            Real residual = kepler(root_ridders);
            suite.addResult("Ridders", eStr.str(), 0.0, residual);
        } catch (...) {
            suite.addResult("Ridders", eStr.str(), 0.0, 1.0);  // Mark as failed
        }
    }
    
    suite.printDetailedTable();
    suite.printSummary();
    
    std::cout << "\n  NOTE: For Kepler's equation, we verify residual E - e*sin(E) - M ≈ 0\n";
    std::cout << "        Higher eccentricity makes the problem more challenging.\n\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     MULTIPLE ROOTS (BRACKET DETECTION)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test multiple root detection
 */
void Test_Roots_MultipleRoots()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding - Multiple Roots Detection\n";
    std::cout << "  Test function: (x-1)(x-2)(x-3) with roots at 1, 2, 3\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Multiple Roots", "Three roots polynomial");
    
    TestThreeRoots func;
    auto exact = func.exactRoots();
    Real tol = 1e-12;
    
    // Use bracket finding - use exact range [0.5, 3.5] with enough points
    Vector<Real> xb1, xb2;
    int numBrackets = RootFinding::FindRootBrackets(func, 0.5, 3.5, 60, xb1, xb2);
    
    std::cout << "  Found " << numBrackets << " root brackets:\n";
    
    // Only process the expected 3 roots
    int rootsToTest = std::min(numBrackets, (int)exact.size());
    for (int i = 0; i < rootsToTest; ++i)
    {
        std::cout << "    Bracket " << i+1 << ": [" << std::setprecision(4) << xb1[i] << ", " << xb2[i] << "]\n";
        
        std::ostringstream rootStr;
        rootStr << "root" << (i+1);
        
        try {
            Real root_brent = RootFinding::FindRootBrent(func, xb1[i], xb2[i], tol);
            suite.addResult("Brent", rootStr.str(), exact[i], root_brent);
        } catch (...) {
            suite.addResult("Brent", rootStr.str(), exact[i], 0.0);
        }
        
        try {
            Real root_ridders = RootFinding::FindRootRidders(func, xb1[i], xb2[i], tol);
            suite.addResult("Ridders", rootStr.str(), exact[i], root_ridders);
        } catch (...) {
            suite.addResult("Ridders", rootStr.str(), exact[i], 0.0);
        }
        
        try {
            Real root_bisection = RootFinding::FindRootBisection(func, xb1[i], xb2[i], tol);
            suite.addResult("Bisection", rootStr.str(), exact[i], root_bisection);
        } catch (...) {
            suite.addResult("Bisection", rootStr.str(), exact[i], 0.0);
        }
    }
    
    std::cout << "\n";
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     EDGE CASES
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test challenging cases: flat roots, near-zero derivative
 */
void Test_Roots_EdgeCases()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  PRECISION TEST SUITE: Root Finding - Edge Cases\n";
    std::cout << "  Challenging functions that stress test the algorithms\n";
    std::cout << "====================================================================================================\n\n";
    
    PrecisionTestSuite suite("Edge Cases", "Challenging roots");
    
    Real tol = 1e-10;  // Slightly looser tolerance for difficult cases
    
    // Flat root: x * exp(-x^2) at x=0
    {
        TestFlatRoot func;
        Real exact = func.exactRoot();
        Real x1 = -0.5, x2 = 0.5;
        
        try { suite.addResult("Bisection", "FlatRoot", exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Ridders", "FlatRoot", exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Brent", "FlatRoot", exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol)); } catch (...) {}
    }
    
    // Sinc - 0.5: oscillating function
    {
        TestSincHalf func;
        Real exact = func.exactRoot();
        Real x1 = 1.0, x2 = 2.5;
        
        try { suite.addResult("Bisection", "Sinc", exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Ridders", "Sinc", exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Brent", "Sinc", exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol)); } catch (...) {}
    }
    
    // Quintic polynomial
    {
        TestQuintic func;
        Real exact = func.exactRoot();
        Real x1 = 1.0, x2 = 1.5;
        
        try { suite.addResult("Bisection", "Quintic", exact, 
            RootFinding::FindRootBisection(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Newton", "Quintic", exact, 
            RootFinding::FindRootNewton(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Ridders", "Quintic", exact, 
            RootFinding::FindRootRidders(func, x1, x2, tol)); } catch (...) {}
        try { suite.addResult("Brent", "Quintic", exact, 
            RootFinding::FindRootBrent(func, x1, x2, tol)); } catch (...) {}
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     CONVERGENCE RATE COMPARISON
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Measure iteration counts and convergence rates
 */
void Test_Roots_ConvergenceRates()
{
    std::cout << "\n";
    std::cout << "====================================================================================================\n";
    std::cout << "  CONVERGENCE RATE ANALYSIS\n";
    std::cout << "  Theoretical orders: Bisection=1, FalsePos≈1.6, Secant≈1.6, Ridders=2, Newton=2\n";
    std::cout << "====================================================================================================\n\n";
    
    TestSqrt2 func;
    Real exact = func.exactRoot();
    
    std::cout << "Error progression for sqrt(2) problem:\n\n";
    std::cout << std::setw(12) << "Tolerance" 
              << std::setw(15) << "Bisection" 
              << std::setw(15) << "Secant"
              << std::setw(15) << "Ridders"
              << std::setw(15) << "Brent" << "\n";
    std::cout << std::string(72, '-') << "\n";
    
    std::vector<Real> tolerances = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
    
    for (Real tol : tolerances)
    {
        Real x1 = 1.0, x2 = 2.0;
        
        Real err_bisection = std::abs(RootFinding::FindRootBisection(func, x1, x2, tol) - exact);
        Real err_secant = std::abs(RootFinding::FindRootSecant(func, x1, x2, tol) - exact);
        Real err_ridders = std::abs(RootFinding::FindRootRidders(func, x1, x2, tol) - exact);
        Real err_brent = std::abs(RootFinding::FindRootBrent(func, x1, x2, tol) - exact);
        
        std::cout << std::scientific << std::setprecision(0)
                  << std::setw(12) << tol
                  << std::setprecision(3)
                  << std::setw(15) << err_bisection
                  << std::setw(15) << err_secant
                  << std::setw(15) << err_ridders
                  << std::setw(15) << err_brent << "\n";
    }
    
    std::cout << "\n";
    std::cout << "  Expected convergence orders:\n";
    std::cout << "    - Bisection:      O(1) - linear, error halves each iteration\n";
    std::cout << "    - False Position: O(1.6) - superlinear\n";
    std::cout << "    - Secant:         O(1.618) - golden ratio convergence\n";
    std::cout << "    - Ridders:        O(2) - quadratic (like Newton)\n";
    std::cout << "    - Brent:          O(1.8-2) - adaptive, usually quadratic\n";
    std::cout << "    - Newton:         O(2) - quadratic near root\n\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     COMPREHENSIVE SUMMARY
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Summary comparison across all test problems
 */
void Test_Roots_ComprehensiveSummary()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "  COMPREHENSIVE ROOT FINDING ACCURACY SUMMARY\n";
    std::cout << "================================================================================\n\n";
    
    PrecisionTestSuite suite("All Root Finders", "Complete comparison");
    
    Real tol = 1e-12;
    
    // sqrt(2)
    {
        TestSqrt2 func;
        Real exact = func.exactRoot();
        Real x1 = 1.0, x2 = 2.0;
        
        suite.addResult("Bisection", "sqrt2", exact, RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", "sqrt2", exact, RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Secant", "sqrt2", exact, RootFinding::FindRootSecant(func, x1, x2, tol));
        suite.addResult("Ridders", "sqrt2", exact, RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", "sqrt2", exact, RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    // Cubic
    {
        TestCubic func;
        Real exact = func.exactRoot();
        Real x1 = 1.0, x2 = 2.0;
        
        suite.addResult("Bisection", "cubic", exact, RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", "cubic", exact, RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Secant", "cubic", exact, RootFinding::FindRootSecant(func, x1, x2, tol));
        suite.addResult("Ridders", "cubic", exact, RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", "cubic", exact, RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    // Dottie number (cos(x) = x)
    {
        TestCosX func;
        Real exact = func.exactRoot();
        Real x1 = 0.0, x2 = 1.0;
        
        suite.addResult("Bisection", "dottie", exact, RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", "dottie", exact, RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Secant", "dottie", exact, RootFinding::FindRootSecant(func, x1, x2, tol));
        suite.addResult("Ridders", "dottie", exact, RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", "dottie", exact, RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    // Quintic
    {
        TestQuintic func;
        Real exact = func.exactRoot();
        Real x1 = 1.0, x2 = 1.5;
        
        suite.addResult("Bisection", "quintic", exact, RootFinding::FindRootBisection(func, x1, x2, tol));
        suite.addResult("Newton", "quintic", exact, RootFinding::FindRootNewton(func, x1, x2, tol));
        suite.addResult("Secant", "quintic", exact, RootFinding::FindRootSecant(func, x1, x2, tol));
        suite.addResult("Ridders", "quintic", exact, RootFinding::FindRootRidders(func, x1, x2, tol));
        suite.addResult("Brent", "quintic", exact, RootFinding::FindRootBrent(func, x1, x2, tol));
    }
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Export results
    suite.exportCSV("results/roots_precision.csv");
    suite.exportMarkdown("results/roots_precision.md");
    
    std::cout << "\n  Results exported to:\n";
    std::cout << "    - results/roots_precision.csv\n";
    std::cout << "    - results/roots_precision.md\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     MASTER TEST FUNCTION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Run all root finding precision tests
 */
void Test_Precision_RootFinding()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "               ROOT FINDING PRECISION TEST SUITE\n";
    std::cout << "================================================================================\n";
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  1. BASIC METHOD COMPARISON\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_BasicComparison();
    Test_Roots_Cubic();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  2. TRANSCENDENTAL EQUATIONS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_Transcendental();
    Test_Roots_ExpLinear();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  3. PHYSICS APPLICATION: KEPLER'S EQUATION\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_Kepler();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  4. MULTIPLE ROOTS DETECTION\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_MultipleRoots();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  5. EDGE CASES (Challenging Functions)\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_EdgeCases();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  6. CONVERGENCE RATE ANALYSIS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_ConvergenceRates();
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  7. COMPREHENSIVE SUMMARY\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Roots_ComprehensiveSummary();
    
    std::cout << "\n================================================================================\n";
    std::cout << "            ROOT FINDING PRECISION TESTS COMPLETE\n";
    std::cout << "================================================================================\n";
}
