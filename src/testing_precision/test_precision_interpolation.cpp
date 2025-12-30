///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        test_precision_interpolation.cpp                                    ///
///  Description: Comprehensive interpolation precision tests                          ///
///               Compares Linear, Polynomial, Spline, Rational, Barycentric          ///
///                                                                                   ///
///  Tests:       - Interpolation at known points                                     ///
///               - Interpolation between points                                      ///
///               - Extrapolation behavior                                            ///
///               - Derivative estimation (for spline)                                ///
///               - Convergence with increasing data points                           ///
///               - Runge phenomenon demonstration                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "PrecisionTestFramework.h"

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/InterpolatedFunction.h"

using namespace MML;
using namespace MML::PrecisionTesting;

///////////////////////////////////////////////////////////////////////////////////////////
//                     BASIC INTERPOLATION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test interpolation at known data points (should be exact or near-exact)
 */
void Test_Interpolation_AtKnownPoints()
{
    PrecisionTestSuite suite("Interpolation at Known Points",
                             "Testing interpolation accuracy at data nodes");
    
    // Create test data: sin(x) at 11 points from 0 to pi
    const int n = 11;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++) {
        x[i] = i * Constants::PI / (n - 1);
        y[i] = sin(x[i]);
    }
    
    // Create interpolators
    LinearInterpRealFunc linear(x, y);
    PolynomInterpRealFunc poly(x, y, 4);       // 4-point polynomial
    SplineInterpRealFunc spline(x, y);          // Cubic spline
    RationalInterpRealFunc rational(x, y, 4);   // 4-point rational
    
    // Test at each data point (should recover exact values)
    for (int i = 0; i < n; i++) {
        double exact = y[i];
        
        suite.addResult("Linear", "node " + std::to_string(i), exact, linear(x[i]));
        suite.addResult("Poly4", "node " + std::to_string(i), exact, poly(x[i]));
        suite.addResult("Spline", "node " + std::to_string(i), exact, spline(x[i]));
        suite.addResult("Rational", "node " + std::to_string(i), exact, rational(x[i]));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test interpolation between data points against known function
 */
void Test_Interpolation_BetweenPoints()
{
    PrecisionTestSuite suite("Interpolation Between Points",
                             "Testing interpolation accuracy at midpoints");
    
    // Create test data: sin(x) at 21 points from 0 to 2*pi
    const int n = 21;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++) {
        x[i] = i * 2.0 * Constants::PI / (n - 1);
        y[i] = sin(x[i]);
    }
    
    // Create interpolators
    LinearInterpRealFunc linear(x, y);
    PolynomInterpRealFunc poly4(x, y, 4);
    PolynomInterpRealFunc poly6(x, y, 6);
    SplineInterpRealFunc spline(x, y);
    RationalInterpRealFunc rational(x, y, 4);
    
    // Test at midpoints between data nodes
    for (int i = 0; i < n - 1; i++) {
        double test_x = (x[i] + x[i + 1]) / 2.0;
        double exact = sin(test_x);
        
        suite.addResult("Linear", "mid " + std::to_string(i), exact, linear(test_x));
        suite.addResult("Poly4", "mid " + std::to_string(i), exact, poly4(test_x));
        suite.addResult("Poly6", "mid " + std::to_string(i), exact, poly6(test_x));
        suite.addResult("Spline", "mid " + std::to_string(i), exact, spline(test_x));
        suite.addResult("Rational", "mid " + std::to_string(i), exact, rational(test_x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     CONVERGENCE TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test convergence with increasing number of data points
 */
void Test_Interpolation_Convergence()
{
    PrecisionTestSuite suite("Interpolation Convergence",
                             "Error reduction with increasing data points");
    
    // Test function: exp(x) on [0, 1]
    auto f = [](Real x) { return exp(x); };
    
    // Test point in middle of domain
    double test_x = 0.5;
    double exact = f(test_x);
    
    std::vector<int> point_counts = {5, 11, 21, 41, 81};
    
    for (int n : point_counts) {
        Vector<Real> x(n), y(n);
        for (int i = 0; i < n; i++) {
            x[i] = (double)i / (n - 1);
            y[i] = f(x[i]);
        }
        
        std::string nStr = "n=" + std::to_string(n);
        
        LinearInterpRealFunc linear(x, y);
        SplineInterpRealFunc spline(x, y);
        
        suite.addResult("Linear", nStr, exact, linear(test_x));
        suite.addResult("Spline", nStr, exact, spline(test_x));
        
        // Polynomial - use min of n and reasonable order
        int polyOrder = std::min(n, 8);
        PolynomInterpRealFunc poly(x, y, polyOrder);
        suite.addResult("Poly" + std::to_string(polyOrder), nStr, exact, poly(test_x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     RUNGE PHENOMENON
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Demonstrate Runge phenomenon with high-degree polynomial interpolation
 */
void Test_Interpolation_RungePhenom()
{
    PrecisionTestSuite suite("Runge Phenomenon",
                             "Polynomial vs Spline on uniform grid");
    
    // Runge function: 1/(1+25x^2) on [-1, 1]
    auto runge = [](Real x) { return 1.0 / (1.0 + 25.0 * x * x); };
    
    // Create data with uniform spacing
    const int n = 21;  // 21 points
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++) {
        x[i] = -1.0 + 2.0 * i / (n - 1);
        y[i] = runge(x[i]);
    }
    
    // Create interpolators
    SplineInterpRealFunc spline(x, y);
    PolynomInterpRealFunc poly(x, y, n);  // Full degree polynomial
    
    // Test at points near the boundary (where Runge phenomenon is worst)
    std::vector<double> test_points = {-0.95, -0.85, -0.75, 0.75, 0.85, 0.95};
    
    for (double tx : test_points) {
        double exact = runge(tx);
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << tx;
        
        suite.addResult("Spline", "x=" + ss.str(), exact, spline(tx));
        suite.addResult("PolyFull", "x=" + ss.str(), exact, poly(tx));
    }
    
    // Also test interior points
    std::vector<double> interior_points = {-0.3, -0.1, 0.0, 0.1, 0.3};
    for (double tx : interior_points) {
        double exact = runge(tx);
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << tx;
        
        suite.addResult("Spline", "x=" + ss.str(), exact, spline(tx));
        suite.addResult("PolyFull", "x=" + ss.str(), exact, poly(tx));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    std::cout << "\n  NOTE: High-order polynomial shows oscillations near boundaries (Runge phenomenon)\n";
    std::cout << "        Spline remains stable throughout the domain.\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     SPLINE DERIVATIVE TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test spline derivative accuracy
 */
void Test_Spline_Derivatives()
{
    PrecisionTestSuite suite("Spline Derivatives",
                             "Testing derivative and second derivative accuracy");
    
    // Test function: sin(x), f'(x) = cos(x), f''(x) = -sin(x)
    const int n = 21;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++) {
        x[i] = i * Constants::PI / (n - 1);  // [0, pi]
        y[i] = sin(x[i]);
    }
    
    SplineInterpRealFunc spline(x, y);
    
    // Test derivatives at interior points
    for (int i = 1; i < n - 1; i++) {
        double test_x = x[i];
        
        // First derivative
        double exact_deriv = cos(test_x);
        double computed_deriv = spline.Derivative(test_x);
        suite.addResult("Deriv1", "x=" + std::to_string(i), exact_deriv, computed_deriv);
        
        // Second derivative
        double exact_sec = -sin(test_x);
        double computed_sec = spline.SecondDerivative(test_x);
        suite.addResult("Deriv2", "x=" + std::to_string(i), exact_sec, computed_sec);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test spline integration accuracy
 */
void Test_Spline_Integration()
{
    PrecisionTestSuite suite("Spline Integration",
                             "Testing definite integral accuracy");
    
    // Test function: sin(x), integral from 0 to x is 1 - cos(x)
    const int n = 21;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++) {
        x[i] = i * Constants::PI / (n - 1);  // [0, pi]
        y[i] = sin(x[i]);
    }
    
    SplineInterpRealFunc spline(x, y);
    
    // Test integral from 0 to various endpoints
    std::vector<double> endpoints = {
        Constants::PI / 4,    // integral = 1 - cos(pi/4) ≈ 0.293
        Constants::PI / 2,    // integral = 1 - cos(pi/2) = 1
        3 * Constants::PI / 4, // integral = 1 - cos(3pi/4) ≈ 1.707
        Constants::PI         // integral = 1 - cos(pi) = 2
    };
    
    for (double b : endpoints) {
        double exact = 1.0 - cos(b);  // ∫sin(x)dx from 0 to b = -cos(b) + cos(0)
        double computed = spline.Integrate(0.0, b);
        
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(3) << b;
        suite.addResult("Integral", "[0," + ss.str() + "]", exact, computed);
    }
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     DIFFERENT FUNCTION TYPES
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test interpolation on various function types
 */
void Test_Interpolation_VariousFunctions()
{
    PrecisionTestSuite suite("Various Functions",
                             "Interpolation accuracy for different function types");
    
    // Test functions and their properties
    struct TestFunc {
        std::string name;
        std::function<double(double)> f;
        double a, b;  // domain
        double test_x;  // test point
    };
    
    std::vector<TestFunc> funcs = {
        {"sin", [](double x) { return sin(x); }, 0.0, Constants::PI, Constants::PI/3},
        {"exp", [](double x) { return exp(x); }, 0.0, 2.0, 1.0},
        {"poly", [](double x) { return x*x*x - 2*x*x + x; }, -1.0, 2.0, 0.5},
        {"sqrt", [](double x) { return sqrt(x); }, 0.1, 4.0, 1.5},
        {"log", [](double x) { return log(x); }, 0.5, 5.0, 2.0},
        {"atan", [](double x) { return atan(x); }, -5.0, 5.0, 1.0}
    };
    
    const int n = 21;  // number of data points
    
    for (const auto& tf : funcs) {
        Vector<Real> x(n), y(n);
        for (int i = 0; i < n; i++) {
            x[i] = tf.a + (tf.b - tf.a) * i / (n - 1);
            y[i] = tf.f(x[i]);
        }
        
        LinearInterpRealFunc linear(x, y);
        SplineInterpRealFunc spline(x, y);
        PolynomInterpRealFunc poly(x, y, 6);
        RationalInterpRealFunc rational(x, y, 4);
        
        double exact = tf.f(tf.test_x);
        
        suite.addResult("Linear", tf.name, exact, linear(tf.test_x));
        suite.addResult("Spline", tf.name, exact, spline(tf.test_x));
        suite.addResult("Poly6", tf.name, exact, poly(tf.test_x));
        suite.addResult("Rational", tf.name, exact, rational(tf.test_x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     BARYCENTRIC INTERPOLATION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test Barycentric Rational Interpolation
 */
void Test_BarycentricInterpolation()
{
    PrecisionTestSuite suite("Barycentric Rational Interpolation",
                             "Testing BarycentricRationalInterp class");
    
    // Test with sin(x) on [0, 2*pi]
    const int n = 21;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; i++) {
        x[i] = i * 2.0 * Constants::PI / (n - 1);
        y[i] = sin(x[i]);
    }
    
    // Create Barycentric interpolator with different orders
    BarycentricRationalInterp bary3(x, y, 3);
    BarycentricRationalInterp bary5(x, y, 5);
    SplineInterpRealFunc spline(x, y);  // For comparison
    
    // Test at random interior points
    std::vector<double> test_points = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    
    for (double tx : test_points) {
        double exact = sin(tx);
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << tx;
        
        suite.addResult("Bary3", "x=" + ss.str(), exact, bary3(tx));
        suite.addResult("Bary5", "x=" + ss.str(), exact, bary5(tx));
        suite.addResult("Spline", "x=" + ss.str(), exact, spline(tx));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     COMPREHENSIVE SUMMARY
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Generate comprehensive interpolation accuracy comparison
 */
void Test_Interpolation_ComprehensiveSummary()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "  COMPREHENSIVE INTERPOLATION ACCURACY SUMMARY\n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    
    PrecisionTestSuite suite("All Interpolators x Test Functions",
                             "Complete error comparison");
    
    // Standard test functions
    struct TestFunc {
        std::string name;
        std::function<double(double)> f;
        double a, b;
    };
    
    std::vector<TestFunc> funcs = {
        {"sin", [](double x) { return sin(x); }, 0.0, 2*Constants::PI},
        {"exp", [](double x) { return exp(x); }, 0.0, 3.0},
        {"sqrt", [](double x) { return sqrt(x); }, 0.1, 4.0},
        {"log", [](double x) { return log(x); }, 0.5, 5.0},
        {"1/(1+x^2)", [](double x) { return 1.0/(1+x*x); }, -3.0, 3.0}
    };
    
    const int n = 21;  // data points
    
    for (const auto& tf : funcs) {
        Vector<Real> x(n), y(n);
        for (int i = 0; i < n; i++) {
            x[i] = tf.a + (tf.b - tf.a) * i / (n - 1);
            y[i] = tf.f(x[i]);
        }
        
        // Test at several interior points and average error
        double err_linear = 0, err_spline = 0, err_poly = 0, err_rational = 0;
        int num_tests = 50;
        
        LinearInterpRealFunc linear(x, y);
        SplineInterpRealFunc spline(x, y);
        PolynomInterpRealFunc poly(x, y, 6);
        RationalInterpRealFunc rational(x, y, 4);
        
        for (int t = 0; t < num_tests; t++) {
            double test_x = tf.a + (tf.b - tf.a) * (t + 0.5) / num_tests;
            double exact = tf.f(test_x);
            
            err_linear += std::abs(linear(test_x) - exact);
            err_spline += std::abs(spline(test_x) - exact);
            err_poly += std::abs(poly(test_x) - exact);
            err_rational += std::abs(rational(test_x) - exact);
        }
        
        // Add average errors (use 0.0 as "exact" to show absolute error)
        suite.addResult("Linear", tf.name, 0.0, err_linear / num_tests);
        suite.addResult("Spline", tf.name, 0.0, err_spline / num_tests);
        suite.addResult("Poly6", tf.name, 0.0, err_poly / num_tests);
        suite.addResult("Rational", tf.name, 0.0, err_rational / num_tests);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Export results
    suite.exportCSV("results/interpolation_precision.csv");
    suite.exportMarkdown("results/interpolation_precision.md");
    
    std::cout << "\n  Results exported to:\n";
    std::cout << "    - results/interpolation_precision.csv\n";
    std::cout << "    - results/interpolation_precision.md\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     MAIN TEST ENTRY POINT
///////////////////////////////////////////////////////////////////////////////////////////

void Test_Precision_Interpolation()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "               INTERPOLATION PRECISION TEST SUITE                              \n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    
    //-----------------------------------------------------------------------------------
    // BASIC INTERPOLATION TESTS
    //-----------------------------------------------------------------------------------
    
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "  1. BASIC INTERPOLATION TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Interpolation_AtKnownPoints();
    Test_Interpolation_BetweenPoints();
    
    //-----------------------------------------------------------------------------------
    // CONVERGENCE AND RUNGE PHENOMENON
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  2. CONVERGENCE AND RUNGE PHENOMENON\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Interpolation_Convergence();
    Test_Interpolation_RungePhenom();
    
    //-----------------------------------------------------------------------------------
    // SPLINE SPECIFIC TESTS
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  3. SPLINE SPECIFIC TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Spline_Derivatives();
    Test_Spline_Integration();
    
    //-----------------------------------------------------------------------------------
    // VARIOUS FUNCTIONS AND BARYCENTRIC
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  4. VARIOUS FUNCTIONS AND BARYCENTRIC\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_Interpolation_VariousFunctions();
    Test_BarycentricInterpolation();
    
    //-----------------------------------------------------------------------------------
    // COMPREHENSIVE SUMMARY
    //-----------------------------------------------------------------------------------
    
    Test_Interpolation_ComprehensiveSummary();
    
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "            INTERPOLATION PRECISION TESTS COMPLETE                             \n";
    std::cout << "================================================================================\n";
}