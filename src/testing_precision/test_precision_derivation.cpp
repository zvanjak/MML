///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        test_precision_derivation.cpp                                       ///
///  Description: Comprehensive precision testing for numerical derivatives           ///
///                                                                                   ///
///  Tests:       NDer1, NDer2, NDer4, NDer6, NDer8 (first derivatives)               ///
///               NSecDer2, NSecDer4 (second derivatives)                             ///
///               NThirdDer2 (third derivatives)                                      ///
///               Gradient computation                                                ///
///               Partial derivatives                                                 ///
///                                                                                   ///
///  Uses:        PrecisionTestFramework for unified output                           ///
///               RealFunctionsTestBed for comprehensive test functions               ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Function.h"
#include "core/Derivation.h"
#endif

#include "PrecisionTestFramework.h"
#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"

using namespace MML;
using namespace MML::TestBeds;
using namespace MML::PrecisionTesting;

///////////////////////////////////////////////////////////////////////////////////////////
//                     FIRST DERIVATIVE PRECISION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test all first derivative orders (NDer1-8) against a single function
 * 
 * Produces an error order matrix showing precision of each derivative order
 * at different test points.
 */
void Test_FirstDerivative_AllOrders_SingleFunc(const std::string& funcName, 
                                                const TestFunctionReal& testFunc,
                                                const std::vector<Real>& testPoints)
{
    PrecisionTestSuite suite("First Derivative: " + funcName, 
                             "Compare NDer1, NDer2, NDer4, NDer6, NDer8 accuracy");
    
    const RealFunction& f = testFunc._func;
    const RealFunction& f_der = testFunc._funcDerived;
    
    for (double x : testPoints) {
        double exact = f_der(x);
        std::string xStr = "x=" + std::to_string(x).substr(0, 5);
        
        // Test all derivative orders
        suite.addResult("NDer1", xStr, exact, Derivation::NDer1(f, x));
        suite.addResult("NDer2", xStr, exact, Derivation::NDer2(f, x));
        suite.addResult("NDer4", xStr, exact, Derivation::NDer4(f, x));
        suite.addResult("NDer6", xStr, exact, Derivation::NDer6(f, x));
        suite.addResult("NDer8", xStr, exact, Derivation::NDer8(f, x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test first derivatives across ALL standard test functions
 * 
 * Shows which derivative order works best for different function types
 */
void Test_FirstDerivative_AllFunctions()
{
    PrecisionTestSuite suite("First Derivative: All Functions", 
                             "NDer8 accuracy across standard test bed");
    
    int numFuncs = RealFunctionsTestBed::getNumFunc();
    
    for (int i = 0; i < numFuncs; i++) {
        const TestFunctionReal& testFunc = RealFunctionsTestBed::getFunc(i);
        
        const RealFunction& f = testFunc._func;
        const RealFunction& f_der = testFunc._funcDerived;
        
        // Test at center of test interval
        double x1 = testFunc._intervalTest->getLowerBound();
        double x2 = testFunc._intervalTest->getUpperBound();
        double x_mid = (x1 + x2) / 2.0;
        
        // Avoid division by zero for some functions
        if (std::abs(x_mid) < 1e-10) x_mid = 1.0;
        
        double exact = f_der(x_mid);
        
        // Test with NDer8 (highest accuracy)
        double computed = Derivation::NDer8(f, x_mid);
        
        PrecisionTestResult result("NDer8", testFunc._funcName, exact, computed);
        result.parameters = "x=" + std::to_string(x_mid);
        suite.addResult(result);
    }
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

/**
 * @brief Test derivative accuracy vs step size h
 * 
 * Shows optimal h for each derivative order
 */
void Test_FirstDerivative_StepSizeAnalysis(int derivOrder)
{
    std::string orderStr = "NDer" + std::to_string(derivOrder);
    PrecisionTestSuite suite(orderStr + " Step Size Analysis", 
                             "Finding optimal step size for " + orderStr);
    
    const TestFunctionReal& testFunc = RealFunctionsTestBed::getFunc("Sin");
    const RealFunction& f = testFunc._func;
    const RealFunction& f_der = testFunc._funcDerived;
    
    double x = 1.0;  // Test point
    double exact = f_der(x);
    
    // Test step sizes from 1e-2 to 1e-12
    std::vector<double> h_values = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
    
    for (double h : h_values) {
        double computed = 0.0;
        
        switch (derivOrder) {
            case 1: computed = Derivation::NDer1(f, x, h); break;
            case 2: computed = Derivation::NDer2(f, x, h); break;
            case 4: computed = Derivation::NDer4(f, x, h); break;
            case 6: computed = Derivation::NDer6(f, x, h); break;
            case 8: computed = Derivation::NDer8(f, x, h); break;
        }
        
        // Format h nicely
        std::ostringstream hStr;
        hStr << std::scientific << std::setprecision(0) << h;
        
        PrecisionTestResult result(orderStr, "h=" + hStr.str(), exact, computed);
        suite.addResult(result);
    }
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

/**
 * @brief Multi-function comparison showing best derivative order per function type
 */
void Test_FirstDerivative_MultiFunction_Comparison()
{
    PrecisionTestSuite suite("First Derivative: Multi-Function Comparison", 
                             "Best derivative order for each function type");
    
    // Select representative functions from different categories
    std::vector<std::string> funcNames = {"Sin", "Cos", "Exp", "x^2", "x^3", "Sqrt", "Tanh"};
    
    for (const auto& funcName : funcNames) {
        try {
            const TestFunctionReal& testFunc = RealFunctionsTestBed::getFunc(funcName);
            const RealFunction& f = testFunc._func;
            const RealFunction& f_der = testFunc._funcDerived;
            
            // Test at a safe point within the interval
            double x1 = testFunc._intervalTest->getLowerBound();
            double x2 = testFunc._intervalTest->getUpperBound();
            double x = (x1 + x2) / 2.0;
            if (funcName == "Sqrt") x = 1.0;  // Avoid 0 for sqrt
            
            double exact = f_der(x);
            
            suite.addResult("NDer1", funcName, exact, Derivation::NDer1(f, x));
            suite.addResult("NDer2", funcName, exact, Derivation::NDer2(f, x));
            suite.addResult("NDer4", funcName, exact, Derivation::NDer4(f, x));
            suite.addResult("NDer6", funcName, exact, Derivation::NDer6(f, x));
            suite.addResult("NDer8", funcName, exact, Derivation::NDer8(f, x));
        }
        catch (...) {
            // Skip functions not in test bed
        }
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     SECOND DERIVATIVE PRECISION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test second derivative accuracy (NSecDer2, NSecDer4)
 */
void Test_SecondDerivative_AllFunctions()
{
    PrecisionTestSuite suite("Second Derivative: All Functions", 
                             "NSecDer2 and NSecDer4 accuracy comparison");
    
    int numFuncs = RealFunctionsTestBed::getNumFunc();
    
    for (int i = 0; i < numFuncs; i++) {
        const TestFunctionReal& testFunc = RealFunctionsTestBed::getFunc(i);
        
        const RealFunction& f = testFunc._func;
        const RealFunction& f_sec_der = testFunc._funcSecDer;
        
        double x1 = testFunc._intervalTest->getLowerBound();
        double x2 = testFunc._intervalTest->getUpperBound();
        double x = (x1 + x2) / 2.0;
        
        // Avoid problematic points
        if (std::abs(x) < 1e-10) x = 1.0;
        if (testFunc._funcName == "Sqrt") x = 4.0;
        
        double exact = f_sec_der(x);
        
        suite.addResult("NSecDer2", testFunc._funcName, exact, Derivation::NSecDer2(f, x));
        suite.addResult("NSecDer4", testFunc._funcName, exact, Derivation::NSecDer4(f, x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

/**
 * @brief Test second derivative with varying step sizes
 */
void Test_SecondDerivative_StepSize()
{
    PrecisionTestSuite suite("Second Derivative: Step Size Analysis", 
                             "NSecDer4 with sin(x) at x=1");
    
    RealFunctionFromStdFunc f([](Real x) { return sin(x); });
    RealFunctionFromStdFunc f_sec([](Real x) { return -sin(x); });
    
    double x = 1.0;
    double exact = f_sec(x);
    
    std::vector<double> h_values = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    
    for (double h : h_values) {
        std::ostringstream hStr;
        hStr << std::scientific << std::setprecision(0) << h;
        
        suite.addResult("NSecDer2", "h=" + hStr.str(), exact, Derivation::NSecDer2(f, x, h));
        suite.addResult("NSecDer4", "h=" + hStr.str(), exact, Derivation::NSecDer4(f, x, h));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     THIRD DERIVATIVE PRECISION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test third derivative accuracy (NThirdDer2)
 */
void Test_ThirdDerivative()
{
    PrecisionTestSuite suite("Third Derivative: NThirdDer2", 
                             "Third derivative accuracy for standard functions");
    
    // Test with sin (third derivative = -cos)
    RealFunctionFromStdFunc f_sin([](Real x) { return sin(x); });
    RealFunctionFromStdFunc f_sin_3rd([](Real x) { return -cos(x); });
    
    // Test with x^4 (third derivative = 24x)
    RealFunctionFromStdFunc f_x4([](Real x) { return x*x*x*x; });
    RealFunctionFromStdFunc f_x4_3rd([](Real x) { return 24*x; });
    
    // Test with exp (third derivative = exp)
    RealFunctionFromStdFunc f_exp([](Real x) { return exp(x); });
    
    std::vector<double> test_points = {0.5, 1.0, 1.5, 2.0};
    
    for (double x : test_points) {
        std::string ptStr = "x=" + std::to_string(x).substr(0, 4);
        
        suite.addResult("NThirdDer2", "sin@" + ptStr, f_sin_3rd(x), Derivation::NThirdDer2(f_sin, x));
        suite.addResult("NThirdDer2", "x^4@" + ptStr, f_x4_3rd(x), Derivation::NThirdDer2(f_x4, x));
        suite.addResult("NThirdDer2", "exp@" + ptStr, exp(x), Derivation::NThirdDer2(f_exp, x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     EDGE CASES AND CHALLENGING FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test derivatives near discontinuities and steep gradients
 */
void Test_EdgeCases()
{
    PrecisionTestSuite suite("Edge Cases", 
                             "Derivatives near challenging points");
    
    // 1. Near-zero: f(x) = x^3, test near x=0
    RealFunctionFromStdFunc f_cubic([](Real x) { return x*x*x; });
    RealFunctionFromStdFunc f_cubic_der([](Real x) { return 3*x*x; });
    
    std::vector<double> near_zero = {0.001, 0.01, 0.1};
    for (double x : near_zero) {
        std::ostringstream xStr;
        xStr << "x^3@" << x;
        suite.addResult("NDer8", xStr.str(), f_cubic_der(x), Derivation::NDer8(f_cubic, x));
    }
    
    // 2. Steep gradient: f(x) = exp(x), large x
    RealFunctionFromStdFunc f_exp([](Real x) { return exp(x); });
    std::vector<double> steep_points = {5.0, 10.0, 15.0};
    for (double x : steep_points) {
        std::ostringstream xStr;
        xStr << "exp@" << x;
        suite.addResult("NDer8", xStr.str(), exp(x), Derivation::NDer8(f_exp, x));
    }
    
    // 3. Oscillatory: f(x) = sin(10x), high frequency
    RealFunctionFromStdFunc f_osc([](Real x) { return sin(10*x); });
    RealFunctionFromStdFunc f_osc_der([](Real x) { return 10*cos(10*x); });
    std::vector<double> osc_points = {0.1, 0.5, 1.0};
    for (double x : osc_points) {
        std::ostringstream xStr;
        xStr << "sin10x@" << x;
        suite.addResult("NDer8", xStr.str(), f_osc_der(x), Derivation::NDer8(f_osc, x));
    }
    
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

/**
 * @brief Test with numerically challenging functions from extended test bed
 */
void Test_NumericallyChallengingFunctions()
{
    PrecisionTestSuite suite("Numerically Challenging Functions", 
                             "Testing extended function set with known difficulties");
    
    int numExtended = RealFunctionsTestBed::getNumFuncExtended();
    
    for (int i = 0; i < numExtended; i++) {
        const auto& testFunc = RealFunctionsTestBed::getFuncExtended(i);
        
        const auto& f = testFunc._func;
        const auto& f_der = testFunc._funcDerived;
        
        double x1 = testFunc._intervalTest->getLowerBound();
        double x2 = testFunc._intervalTest->getUpperBound();
        double x = (x1 + x2) / 2.0;
        
        // Avoid problematic points
        if (std::abs(x) < 0.1) x = 1.0;
        
        double exact = f_der(x);
        
        suite.addResult("NDer4", testFunc._funcName, exact, Derivation::NDer4(f, x));
        suite.addResult("NDer8", testFunc._funcName, exact, Derivation::NDer8(f, x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     COMPREHENSIVE SUMMARY TABLES
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Generate comprehensive derivative accuracy tables
 * 
 * Shows error orders for all derivative methods across all standard functions
 */
void Test_ComprehensiveSummary()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "  COMPREHENSIVE DERIVATIVE ACCURACY SUMMARY\n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    
    // Test all standard functions with all derivative orders at their midpoint
    PrecisionTestSuite suite("All Derivatives x All Functions", 
                             "Complete error order matrix");
    
    int numFuncs = RealFunctionsTestBed::getNumFunc();
    
    for (int i = 0; i < numFuncs; i++) {
        const TestFunctionReal& testFunc = RealFunctionsTestBed::getFunc(i);
        
        const RealFunction& f = testFunc._func;
        const RealFunction& f_der = testFunc._funcDerived;
        
        double x1 = testFunc._intervalTest->getLowerBound();
        double x2 = testFunc._intervalTest->getUpperBound();
        double x = (x1 + x2) / 2.0;
        
        // Adjust for problematic functions
        if (testFunc._funcName == "Sqrt") x = 4.0;
        if (std::abs(x) < 0.1) x = 1.0;
        
        double exact = f_der(x);
        
        suite.addResult("NDer1", testFunc._funcName, exact, Derivation::NDer1(f, x));
        suite.addResult("NDer2", testFunc._funcName, exact, Derivation::NDer2(f, x));
        suite.addResult("NDer4", testFunc._funcName, exact, Derivation::NDer4(f, x));
        suite.addResult("NDer6", testFunc._funcName, exact, Derivation::NDer6(f, x));
        suite.addResult("NDer8", testFunc._funcName, exact, Derivation::NDer8(f, x));
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Export to files
    suite.exportCSV("results/derivative_precision.csv");
    suite.exportMarkdown("results/derivative_precision.md");
    
    std::cout << "\n  Results exported to:\n";
    std::cout << "    - results/derivative_precision.csv\n";
    std::cout << "    - results/derivative_precision.md\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     LEGACY COMPATIBILITY FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////

// Keep the original function signatures for backward compatibility
void NDer_Error_Order_Diff_Der_Orders_Single_Func(std::string funcName, 
                                                   const TestFunctionReal& inFunc, 
                                                   std::vector<Real> intervalPoints)
{
    Test_FirstDerivative_AllOrders_SingleFunc(funcName, inFunc, intervalPoints);
}

void NDer_Average_Error_Diff_Orders_Single_Func(std::string funcName, const TestFunctionReal& inFunc)
{
    // Legacy output format
    const RealFunction& f = inFunc._func;
    const RealFunction& f_der = inFunc._funcDerived;
    double x1 = inFunc._intervalTest->getLowerBound();
    double x2 = inFunc._intervalTest->getUpperBound();

    const int numPntForEval = 20;

    double err_sum1 = 0.0, err_sum2 = 0.0, err_sum4 = 0.0, err_sum6 = 0.0, err_sum8 = 0.0;

    std::cout << "\nAVERAGE DERIVATION ERROR FOR DIFFERENT ORDERS" << std::endl;
    std::cout << "    X    Exact der.       Nder1        NDer1 err.          Nder2         NDer2 err.          Nder4         NDer4 err.          Nder6         NDer6 err.          Nder8         NDer8 err.         " << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

    for (int i = 0; i < numPntForEval; i++) {
        double x = x1 + (x2 - x1) * i / (numPntForEval - 1);
        double exact_der = f_der(x);

        double num_der1 = Derivation::NDer1(f, x);
        double num_der2 = Derivation::NDer2(f, x);
        double num_der4 = Derivation::NDer4(f, x);
        double num_der6 = Derivation::NDer6(f, x);
        double num_der8 = Derivation::NDer8(f, x);

        double err1 = num_der1 - exact_der;
        double err2 = num_der2 - exact_der;
        double err4 = num_der4 - exact_der;
        double err6 = num_der6 - exact_der;
        double err8 = num_der8 - exact_der;

        std::cout << std::setw(6) << std::setprecision(3) << x
            << std::setw(13) << std::setprecision(8) << exact_der << " "
            << std::setw(13) << std::setprecision(8) << num_der1 << "   "
            << std::scientific << std::setw(15) << err1 << "   " << std::fixed
            << std::setw(13) << std::setprecision(8) << num_der2 << "   "
            << std::scientific << std::setw(15) << err2 << "   " << std::fixed
            << std::setw(13) << std::setprecision(8) << num_der4 << "   "
            << std::scientific << std::setw(15) << err4 << "   " << std::fixed
            << std::setw(13) << std::setprecision(8) << num_der6 << "   "
            << std::scientific << std::setw(15) << err6 << "   " << std::fixed
            << std::setw(13) << std::setprecision(8) << num_der8 << "   "
            << std::scientific << std::setw(15) << err8 << "   " << std::fixed
            << std::endl;

        err_sum1 += std::abs(err1);
        err_sum2 += std::abs(err2);
        err_sum4 += std::abs(err4);
        err_sum6 += std::abs(err6);
        err_sum8 += std::abs(err8);
    }

    std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Total abs error =                    " << std::scientific << err_sum1 << "                    "
        << err_sum2 << "                    " << err_sum4 << "                    "
        << err_sum6 << "                    " << err_sum8 << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
//                     MAIN TEST ENTRY POINT
///////////////////////////////////////////////////////////////////////////////////////////

void Test_Precision_Derivation()
{
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "                 DERIVATION PRECISION TEST SUITE                               \n";
    std::cout << "================================================================================\n";
    std::cout << "\n";
    
    //-----------------------------------------------------------------------------------
    // FIRST DERIVATIVE TESTS
    //-----------------------------------------------------------------------------------
    
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "  1. FIRST DERIVATIVE TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    // 1a. All orders for sin(x)
    std::vector<Real> testPoints = {0.0, 0.5, 1.0, 2.0, 3.0, 5.0};
    Test_FirstDerivative_AllOrders_SingleFunc("Sin", RealFunctionsTestBed::getFunc("Sin"), testPoints);
    
    // 1b. Legacy average error display
    NDer_Average_Error_Diff_Orders_Single_Func("Sin", RealFunctionsTestBed::getFunc(0));
    
    // 1c. Multi-function comparison
    Test_FirstDerivative_MultiFunction_Comparison();
    
    // 1d. Step size analysis for NDer4
    Test_FirstDerivative_StepSizeAnalysis(4);
    
    //-----------------------------------------------------------------------------------
    // SECOND DERIVATIVE TESTS  
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  2. SECOND DERIVATIVE TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_SecondDerivative_AllFunctions();
    Test_SecondDerivative_StepSize();
    
    //-----------------------------------------------------------------------------------
    // THIRD DERIVATIVE TESTS
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  3. THIRD DERIVATIVE TESTS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_ThirdDerivative();
    
    //-----------------------------------------------------------------------------------
    // EDGE CASES AND CHALLENGING FUNCTIONS
    //-----------------------------------------------------------------------------------
    
    std::cout << "\n--------------------------------------------------------------------------------\n";
    std::cout << "  4. EDGE CASES AND CHALLENGING FUNCTIONS\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    Test_EdgeCases();
    Test_NumericallyChallengingFunctions();
    
    //-----------------------------------------------------------------------------------
    // COMPREHENSIVE SUMMARY
    //-----------------------------------------------------------------------------------
    
    Test_ComprehensiveSummary();
    
    std::cout << "\n";
    std::cout << "================================================================================\n";
    std::cout << "              DERIVATION PRECISION TESTS COMPLETE                              \n";
    std::cout << "================================================================================\n";
}
