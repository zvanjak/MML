///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme14_function_analysis.cpp                                      ///
///  Description: README example - Function Analysis                                  ///
///               Demonstrates RealFunctionAnalyzer for interval analysis             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "algorithms/FunctionsAnalyzer.h"
#include "algorithms/RootFinding.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_FunctionAnalysis()
{
    std::cout << std::endl;
    std::cout << "=== Function Analysis ===" << std::endl;

    // Analyze function behavior over an interval
    RealFunction f{[](Real x) { return x*x*x - 3*x + 1; }};

    std::cout << "Analyzing f(x) = x³ - 3x + 1" << std::endl;
    std::cout << std::endl;

    RealFunctionAnalyzer analyzer(f, "x³ - 3x + 1");
    analyzer.PrintIntervalAnalysis(-3.0, 3.0, 100, 1e-6);

    // Find roots using root bracket search + bisection
    Vector<Real> xb1, xb2;
    int numBrackets = RootFinding::FindRootBrackets(f, -3.0, 3.0, 100, xb1, xb2);
    std::cout << std::endl << "Found " << numBrackets << " root brackets:" << std::endl;
    std::cout << std::setprecision(10);
    for (int i = 0; i < numBrackets; i++) {
        Real root = RootFinding::FindRootBisection(f, xb1[i], xb2[i], 1e-10);
        std::cout << "  x = " << root << ", f(x) = " << f(root) << std::endl;
    }

    // Analyze a function with discontinuity
    std::cout << std::endl << "--- Analyzing step function ---" << std::endl;
    RealFunctionFromStdFunc step([](Real x) -> Real { 
        if (x < 0) return 0.0;
        else if (x > 0) return 1.0;
        else return 0.5;
    });
    RealFunctionAnalyzer step_analyzer(step, "step(x)");
    step_analyzer.PrintIntervalAnalysis(-2.0, 2.0, 100, 1e-6);

    // Analyze a function with singularity
    std::cout << std::endl << "--- Analyzing 1/(x-1) ---" << std::endl;
    RealFunction singular{[](Real x) { return 1.0 / (x - 1.0); }};
    RealFunctionAnalyzer sing_analyzer(singular, "1/(x-1)");
    sing_analyzer.PrintIntervalAnalysis(-2.0, 3.0, 100, 1e-6);

    // Analyze oscillating function
    std::cout << std::endl << "--- Analyzing sin(x) ---" << std::endl;
    RealFunction sinf{[](Real x) { return std::sin(x); }};
    RealFunctionAnalyzer sin_analyzer(sinf, "sin(x)");
    sin_analyzer.PrintIntervalAnalysis(0.0, 4*Constants::PI, 200, 1e-6);

    // Find roots of sin(x) using root bracketing
    Vector<Real> sin_xb1, sin_xb2;
    int sin_numBrackets = RootFinding::FindRootBrackets(sinf, 0.1, 4*Constants::PI - 0.1, 200, sin_xb1, sin_xb2);
    std::cout << std::endl << "Roots of sin(x) in [0, 4π]: " << sin_numBrackets << " roots" << std::endl;
    for (int i = 0; i < sin_numBrackets; i++) {
        Real root = RootFinding::FindRootBisection(sinf, sin_xb1[i], sin_xb2[i], 1e-10);
        std::cout << "  x = " << root << " ≈ " << root/Constants::PI << "π" << std::endl;
    }

    std::cout << std::endl;
}
