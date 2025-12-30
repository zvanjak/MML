///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        test_precision_integration.cpp                                      ///
///  Description: Comprehensive precision testing for numerical integration           ///
///                                                                                   ///
///  Tests:       1D: Trapezoidal, Simpson, Romberg, Gauss-Legendre (10-100 pts)      ///
///               2D: Surface integration                                              ///
///               3D: Volume integration                                               ///
///               Improper integrals                                                   ///
///                                                                                   ///
///  Uses:        PrecisionTestFramework for unified output                           ///
///               RealFunctionsTestBed for test functions with known integrals        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Function.h"
#include "core/Integration.h"
#endif

#include "PrecisionTestFramework.h"
#include "../test_data/real_functions_test_bed.h"

using namespace MML;
using namespace MML::TestBeds;
using namespace MML::PrecisionTesting;

///////////////////////////////////////////////////////////////////////////////////////////
//                     1D INTEGRATION PRECISION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compare all 1D integrators on a single function
 */
void Test_1D_AllIntegrators_SingleFunc(const std::string& funcName,
                                        const TestFunctionRealWithIntegral& testFunc)
{
    PrecisionTestSuite suite("1D Integration: " + funcName,
                             "Compare Trap, Simpson, Romberg, Gauss on " + funcName);
    
    const IRealFunction& f = testFunc._func;
    const IRealFunction& f_int = testFunc._funcIntegrated;
    
    double a = testFunc._intervalTest->getLowerBound();
    double b = testFunc._intervalTest->getUpperBound();
    
    // Exact integral
    double exact = f_int(b) - f_int(a);
    
    // Test with default tolerance
    auto result_trap = IntegrateTrap(f, a, b);
    auto result_simp = IntegrateSimpson(f, a, b);
    auto result_romb = IntegrateRomberg(f, a, b);
    auto result_gauss = IntegrateGauss10(f, a, b);
    
    PrecisionTestResult r_trap("Trapezoidal", funcName, exact, result_trap.value);
    r_trap.iterations = result_trap.iterations;
    r_trap.converged = result_trap.converged;
    suite.addResult(r_trap);
    
    PrecisionTestResult r_simp("Simpson", funcName, exact, result_simp.value);
    r_simp.iterations = result_simp.iterations;
    r_simp.converged = result_simp.converged;
    suite.addResult(r_simp);
    
    PrecisionTestResult r_romb("Romberg", funcName, exact, result_romb.value);
    r_romb.iterations = result_romb.iterations;
    r_romb.converged = result_romb.converged;
    suite.addResult(r_romb);
    
    PrecisionTestResult r_gauss("Gauss10", funcName, exact, result_gauss.value);
    r_gauss.converged = true;
    suite.addResult(r_gauss);
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

/**
 * @brief Test all standard functions with all integrators
 */
void Test_1D_AllFunctions_AllIntegrators()
{
    PrecisionTestSuite suite("1D Integration: All Functions",
                             "Complete comparison across test bed");
    
    int numFuncs = RealFunctionsTestBed::getNumFuncWithIntegral();
    
    for (int i = 0; i < numFuncs; i++) {
        const auto& testFunc = RealFunctionsTestBed::getFuncWithIntegral(i);
        
        const IRealFunction& f = testFunc._func;
        const IRealFunction& f_int = testFunc._funcIntegrated;
        
        double a = testFunc._intervalTest->getLowerBound();
        double b = testFunc._intervalTest->getUpperBound();
        double exact = f_int(b) - f_int(a);
        
        // Skip if interval is problematic
        if (!std::isfinite(exact) || std::abs(b - a) > 1e6) continue;
        
        auto result_trap = IntegrateTrap(f, a, b);
        auto result_simp = IntegrateSimpson(f, a, b);
        auto result_romb = IntegrateRomberg(f, a, b);
        auto result_gauss = IntegrateGauss10(f, a, b);
        
        suite.addResult("Trap", testFunc._funcName, exact, result_trap.value);
        suite.addResult("Simpson", testFunc._funcName, exact, result_simp.value);
        suite.addResult("Romberg", testFunc._funcName, exact, result_romb.value);
        suite.addResult("Gauss10", testFunc._funcName, exact, result_gauss.value);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test accuracy vs tolerance setting
 */
void Test_1D_ToleranceSweep()
{
    PrecisionTestSuite suite("1D Integration: Tolerance Sweep",
                             "How accuracy varies with requested tolerance");
    
    // Use sin(x) from 0 to pi (integral = 2)
    RealFunction f = [](Real x) { return sin(x); };
    double a = 0.0;
    double b = Constants::PI;
    double exact = 2.0;
    
    std::vector<double> tolerances = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    
    for (double tol : tolerances) {
        std::ostringstream tolStr;
        tolStr << std::scientific << std::setprecision(0) << tol;
        
        auto result_simp = IntegrateSimpson(f, a, b, tol);
        auto result_romb = IntegrateRomberg(f, a, b, tol);
        
        PrecisionTestResult r_simp("Simpson", "tol=" + tolStr.str(), exact, result_simp.value);
        r_simp.iterations = result_simp.iterations;
        r_simp.converged = result_simp.converged;
        suite.addResult(r_simp);
        
        PrecisionTestResult r_romb("Romberg", "tol=" + tolStr.str(), exact, result_romb.value);
        r_romb.iterations = result_romb.iterations;
        r_romb.converged = result_romb.converged;
        suite.addResult(r_romb);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test Gauss-Legendre with varying number of points
 */
void Test_1D_GaussLegendre_PointSweep()
{
    PrecisionTestSuite suite("Gauss-Legendre: Point Count Sweep",
                             "Accuracy vs number of quadrature points");
    
    // Test with polynomial: x^5 on [0, 1] (integral = 1/6)
    RealFunction f_poly = [](Real x) { return x*x*x*x*x; };
    double exact_poly = 1.0 / 6.0;
    
    // Test with sin(x) on [0, pi] (integral = 2)
    RealFunction f_sin = [](Real x) { return sin(x); };
    double exact_sin = 2.0;
    
    // Test with exp(x) on [0, 1] (integral = e - 1)
    RealFunction f_exp = [](Real x) { return exp(x); };
    double exact_exp = exp(1.0) - 1.0;
    
    std::vector<int> n_points = {5, 10, 15, 20, 30, 50};
    
    for (int n : n_points) {
        std::string nStr = "n=" + std::to_string(n);
        
        // Compute using Gauss-Legendre with n points
        std::vector<Real> nodes(n), weights(n);
        
        // For polynomial: [0, 1]
        GaussLegendre(nodes, weights, 0.0, 1.0);
        double computed_poly = 0.0;
        for (int i = 0; i < n; i++)
            computed_poly += weights[i] * f_poly(nodes[i]);
        suite.addResult("GaussLeg", "x^5 " + nStr, exact_poly, computed_poly);
        
        // For sin: [0, pi]
        GaussLegendre(nodes, weights, 0.0, Constants::PI);
        double computed_sin = 0.0;
        for (int i = 0; i < n; i++)
            computed_sin += weights[i] * f_sin(nodes[i]);
        suite.addResult("GaussLeg", "sin " + nStr, exact_sin, computed_sin);
        
        // For exp: [0, 1]
        GaussLegendre(nodes, weights, 0.0, 1.0);
        double computed_exp = 0.0;
        for (int i = 0; i < n; i++)
            computed_exp += weights[i] * f_exp(nodes[i]);
        suite.addResult("GaussLeg", "exp " + nStr, exact_exp, computed_exp);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     CHALLENGING INTEGRATION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test with numerically challenging functions
 */
void Test_1D_ChallengingFunctions()
{
    PrecisionTestSuite suite("1D Integration: Challenging Functions",
                             "Functions that stress numerical methods");
    
    // 1. Highly oscillatory: sin(50x) on [0, 2pi] (integral = 0)
    RealFunction f_osc = [](Real x) { return sin(50*x); };
    double exact_osc = 0.0;
    
    auto r_osc_trap = IntegrateTrap(f_osc, 0.0, 2*Constants::PI);
    auto r_osc_simp = IntegrateSimpson(f_osc, 0.0, 2*Constants::PI);
    auto r_osc_romb = IntegrateRomberg(f_osc, 0.0, 2*Constants::PI);
    
    suite.addResult("Trap", "sin(50x)", exact_osc, r_osc_trap.value);
    suite.addResult("Simpson", "sin(50x)", exact_osc, r_osc_simp.value);
    suite.addResult("Romberg", "sin(50x)", exact_osc, r_osc_romb.value);
    
    // 2. Peaked function: 1/(1+100x^2) on [-1, 1] (integral = 2*atan(10)/10)
    RealFunction f_peak = [](Real x) { return 1.0 / (1.0 + 100*x*x); };
    double exact_peak = 2.0 * atan(10.0) / 10.0;
    
    auto r_peak_trap = IntegrateTrap(f_peak, -1.0, 1.0);
    auto r_peak_simp = IntegrateSimpson(f_peak, -1.0, 1.0);
    auto r_peak_romb = IntegrateRomberg(f_peak, -1.0, 1.0);
    
    suite.addResult("Trap", "peaked", exact_peak, r_peak_trap.value);
    suite.addResult("Simpson", "peaked", exact_peak, r_peak_simp.value);
    suite.addResult("Romberg", "peaked", exact_peak, r_peak_romb.value);
    
    // 3. Near-singular: sqrt(x) on [0, 1] (integral = 2/3)
    RealFunction f_sqrt = [](Real x) { return (x > 1e-15) ? sqrt(x) : 0.0; };
    double exact_sqrt = 2.0 / 3.0;
    
    auto r_sqrt_trap = IntegrateTrap(f_sqrt, 0.0, 1.0);
    auto r_sqrt_simp = IntegrateSimpson(f_sqrt, 0.0, 1.0);
    auto r_sqrt_romb = IntegrateRomberg(f_sqrt, 0.0, 1.0);
    
    suite.addResult("Trap", "sqrt(x)", exact_sqrt, r_sqrt_trap.value);
    suite.addResult("Simpson", "sqrt(x)", exact_sqrt, r_sqrt_simp.value);
    suite.addResult("Romberg", "sqrt(x)", exact_sqrt, r_sqrt_romb.value);
    
    // 4. Exponentially decaying: exp(-x^2) on [0, 5] (approx sqrt(pi)/2)
    RealFunction f_gauss = [](Real x) { return exp(-x*x); };
    double exact_gauss = 0.5 * sqrt(Constants::PI) * std::erf(5.0);
    
    auto r_gauss_trap = IntegrateTrap(f_gauss, 0.0, 5.0);
    auto r_gauss_simp = IntegrateSimpson(f_gauss, 0.0, 5.0);
    auto r_gauss_romb = IntegrateRomberg(f_gauss, 0.0, 5.0);
    
    suite.addResult("Trap", "exp(-x^2)", exact_gauss, r_gauss_trap.value);
    suite.addResult("Simpson", "exp(-x^2)", exact_gauss, r_gauss_simp.value);
    suite.addResult("Romberg", "exp(-x^2)", exact_gauss, r_gauss_romb.value);
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test with extended function set from test bed
 */
void Test_1D_ExtendedTestBed()
{
    PrecisionTestSuite suite("1D Integration: Extended Test Bed",
                             "Challenging functions from extended test bed");
    
    int numExtended = RealFunctionsTestBed::getNumFuncExtendedWithIntegral();
    
    for (int i = 0; i < numExtended; i++) {
        const auto& testFunc = RealFunctionsTestBed::getFuncExtendedWithIntegral(i);
        
        const IRealFunction& f = testFunc._func;
        const IRealFunction& f_int = testFunc._funcIntegrated;
        
        double a = testFunc._intervalTest->getLowerBound();
        double b = testFunc._intervalTest->getUpperBound();
        double exact = f_int(b) - f_int(a);
        
        if (!std::isfinite(exact)) continue;
        
        auto result_simp = IntegrateSimpson(f, a, b);
        auto result_romb = IntegrateRomberg(f, a, b);
        
        suite.addResult("Simpson", testFunc._funcName, exact, result_simp.value);
        suite.addResult("Romberg", testFunc._funcName, exact, result_romb.value);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     2D AND 3D INTEGRATION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

// Helper functions for constant bounds in 2D integration
Real ConstY0(Real) { return 0.0; }
Real ConstY1(Real) { return 1.0; }
Real ConstY2(Real) { return 2.0; }
Real ConstY3(Real) { return 3.0; }
Real ConstYPi(Real) { return Constants::PI; }
Real ConstYPiHalf(Real) { return Constants::PI / 2.0; }
Real ConstYNeg1(Real) { return -1.0; }

// Helper functions for constant bounds in 3D integration
Real ConstZ0(Real, Real) { return 0.0; }
Real ConstZ1(Real, Real) { return 1.0; }
Real ConstZ2(Real, Real) { return 2.0; }
Real ConstZ3(Real, Real) { return 3.0; }
Real ConstZ4(Real, Real) { return 4.0; }

/**
 * @brief Test 2D surface integration
 */
void Test_2D_SurfaceIntegration()
{
    PrecisionTestSuite suite("2D Surface Integration",
                             "Integrate over rectangular regions");
    
    // 1. f(x,y) = x*y over [0,1]x[0,1] (integral = 0.25)
    ScalarFunction<2> f_xy = [](const VectorN<Real, 2>& p) { return p[0] * p[1]; };
    double exact_xy = 0.25;
    
    auto result_xy = Integrate2D(f_xy, IntegrationMethod::ROMBERG, 0.0, 1.0, ConstY0, ConstY1);
    PrecisionTestResult r_xy("Surface", "x*y [0,1]^2", exact_xy, result_xy.value);
    r_xy.converged = result_xy.converged;
    suite.addResult(r_xy);
    
    // 2. f(x,y) = sin(x)*cos(y) over [0,pi]x[0,pi/2] (integral = 2)
    ScalarFunction<2> f_sincos = [](const VectorN<Real, 2>& p) { return sin(p[0]) * cos(p[1]); };
    double exact_sincos = 2.0;
    
    auto result_sincos = Integrate2D(f_sincos, IntegrationMethod::ROMBERG, 0.0, Constants::PI, ConstY0, ConstYPiHalf);
    PrecisionTestResult r_sincos("Surface", "sin*cos", exact_sincos, result_sincos.value);
    r_sincos.converged = result_sincos.converged;
    suite.addResult(r_sincos);
    
    // 3. f(x,y) = x^2 + y^2 over [0,1]x[0,1] (integral = 2/3)
    ScalarFunction<2> f_circ = [](const VectorN<Real, 2>& p) { return p[0]*p[0] + p[1]*p[1]; };
    double exact_circ = 2.0 / 3.0;
    
    auto result_circ = Integrate2D(f_circ, IntegrationMethod::ROMBERG, 0.0, 1.0, ConstY0, ConstY1);
    PrecisionTestResult r_circ("Surface", "x^2+y^2", exact_circ, result_circ.value);
    r_circ.converged = result_circ.converged;
    suite.addResult(r_circ);
    
    // 4. f(x,y) = exp(-(x^2+y^2)) over [-1,1]x[-1,1] (approx pi*(1-exp(-2)))
    ScalarFunction<2> f_gauss2d = [](const VectorN<Real, 2>& p) { return exp(-(p[0]*p[0] + p[1]*p[1])); };
    // Exact value for [-1,1]x[-1,1]: (sqrt(pi)*erf(1))^2
    double erf1 = std::erf(1.0);
    double exact_gauss2d = Constants::PI * erf1 * erf1;
    
    auto result_gauss2d = Integrate2D(f_gauss2d, IntegrationMethod::ROMBERG, -1.0, 1.0, ConstYNeg1, ConstY1);
    PrecisionTestResult r_gauss2d("Surface", "exp(-r^2)", exact_gauss2d, result_gauss2d.value);
    r_gauss2d.converged = result_gauss2d.converged;
    suite.addResult(r_gauss2d);
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

/**
 * @brief Test 3D volume integration
 */
void Test_3D_VolumeIntegration()
{
    PrecisionTestSuite suite("3D Volume Integration",
                             "Integrate over rectangular boxes");
    
    // Helper bound functions for 3D
    auto y0 = [](Real) { return 0.0; };
    auto y1 = [](Real) { return 1.0; };
    auto y2 = [](Real) { return 2.0; };
    auto y3 = [](Real) { return 3.0; };
    auto z0 = [](Real, Real) { return 0.0; };
    auto z1 = [](Real, Real) { return 1.0; };
    auto z4 = [](Real, Real) { return 4.0; };
    
    // 1. f(x,y,z) = x*y*z over [0,1]^3 (integral = 1/8)
    ScalarFunction<3> f_xyz = [](const VectorN<Real, 3>& p) { return p[0] * p[1] * p[2]; };
    double exact_xyz = 1.0 / 8.0;
    
    auto result_xyz = Integrate3D(f_xyz, IntegrationMethod::ROMBERG, 0.0, 1.0, ConstY0, ConstY1, ConstZ0, ConstZ1);
    PrecisionTestResult r_xyz("Volume", "xyz [0,1]^3", exact_xyz, result_xyz.value);
    r_xyz.converged = result_xyz.converged;
    suite.addResult(r_xyz);
    
    // 2. f(x,y,z) = x^2 + y^2 + z^2 over [0,1]^3 (integral = 1)
    ScalarFunction<3> f_r2 = [](const VectorN<Real, 3>& p) { return p[0]*p[0] + p[1]*p[1] + p[2]*p[2]; };
    double exact_r2 = 1.0;
    
    auto result_r2 = Integrate3D(f_r2, IntegrationMethod::ROMBERG, 0.0, 1.0, ConstY0, ConstY1, ConstZ0, ConstZ1);
    PrecisionTestResult r_r2("Volume", "r^2 [0,1]^3", exact_r2, result_r2.value);
    r_r2.converged = result_r2.converged;
    suite.addResult(r_r2);
    
    // 3. f(x,y,z) = 1 over [0,2]x[0,3]x[0,4] (integral = 24)
    ScalarFunction<3> f_const = [](const VectorN<Real, 3>& p) { return 1.0; };
    double exact_const = 24.0;  // 2 * 3 * 4
    
    auto result_const = Integrate3D(f_const, IntegrationMethod::ROMBERG, 0.0, 2.0, ConstY0, ConstY3, ConstZ0, ConstZ4);
    PrecisionTestResult r_const("Volume", "1 [0,2]x[0,3]x[0,4]", exact_const, result_const.value);
    r_const.converged = result_const.converged;
    suite.addResult(r_const);
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     COMPREHENSIVE SUMMARY
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Generate comprehensive integration accuracy tables
 */
void Test_Integration_ComprehensiveSummary()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "  COMPREHENSIVE INTEGRATION ACCURACY SUMMARY\n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    
    PrecisionTestSuite suite("All Integrators x All Functions",
                             "Complete error order matrix");
    
    int numFuncs = RealFunctionsTestBed::getNumFuncWithIntegral();
    
    for (int i = 0; i < std::min(numFuncs, 10); i++) {  // Limit to first 10 for readable output
        const auto& testFunc = RealFunctionsTestBed::getFuncWithIntegral(i);
        
        const IRealFunction& f = testFunc._func;
        const IRealFunction& f_int = testFunc._funcIntegrated;
        
        double a = testFunc._intervalTest->getLowerBound();
        double b = testFunc._intervalTest->getUpperBound();
        double exact = f_int(b) - f_int(a);
        
        if (!std::isfinite(exact) || std::abs(b - a) > 1e6) continue;
        
        auto result_trap = IntegrateTrap(f, a, b, 1e-10);
        auto result_simp = IntegrateSimpson(f, a, b, 1e-10);
        auto result_romb = IntegrateRomberg(f, a, b, 1e-10);
        auto result_gauss = IntegrateGauss10(f, a, b);
        
        suite.addResult("Trap", testFunc._funcName, exact, result_trap.value);
        suite.addResult("Simpson", testFunc._funcName, exact, result_simp.value);
        suite.addResult("Romberg", testFunc._funcName, exact, result_romb.value);
        suite.addResult("Gauss10", testFunc._funcName, exact, result_gauss.value);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Export to files
    suite.exportCSV("results/integration_precision.csv");
    suite.exportMarkdown("results/integration_precision.md");
    
    std::cout << "\n  Results exported to:\n";
    std::cout << "    - results/integration_precision.csv\n";
    std::cout << "    - results/integration_precision.md\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     MAIN TEST ENTRY POINT
///////////////////////////////////////////////////////////////////////////////////////////

void Test_Precision_Integration()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "                INTEGRATION PRECISION TEST SUITE                               \n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    
    //-----------------------------------------------------------------------------------
    // 1D INTEGRATION TESTS
    //-----------------------------------------------------------------------------------
    
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "  1. 1D INTEGRATION TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    // 1a. Single function comparison (use TestInt2 which is 1/(4+x^2))
    Test_1D_AllIntegrators_SingleFunc("TestInt2", RealFunctionsTestBed::getFuncWithIntegral("TestInt2"));
    
    // 1b. All functions comparison
    Test_1D_AllFunctions_AllIntegrators();
    
    // 1c. Tolerance sweep
    Test_1D_ToleranceSweep();
    
    // 1d. Gauss-Legendre point sweep
    Test_1D_GaussLegendre_PointSweep();
    
    //-----------------------------------------------------------------------------------
    // CHALLENGING FUNCTIONS
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  2. CHALLENGING INTEGRATION TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_1D_ChallengingFunctions();
    Test_1D_ExtendedTestBed();
    
    //-----------------------------------------------------------------------------------
    // 2D AND 3D INTEGRATION
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  3. MULTI-DIMENSIONAL INTEGRATION\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_2D_SurfaceIntegration();
    Test_3D_VolumeIntegration();
    
    //-----------------------------------------------------------------------------------
    // COMPREHENSIVE SUMMARY
    //-----------------------------------------------------------------------------------
    
    Test_Integration_ComprehensiveSummary();
    
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "             INTEGRATION PRECISION TESTS COMPLETE                              \n";
    std::cout << "================================================================================\n";
}
