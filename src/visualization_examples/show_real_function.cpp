/**
 * @file show_real_function.cpp
 * @brief Demonstrates visualization of single Real functions
 * 
 * Uses the MML Visualizer infrastructure to display real-valued functions.
 * Cross-platform: works on Windows, Linux, and macOS.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Function.h"
#include "core/FunctionHelpers.h"
#include "tools/Visualizer.h"
#endif

using namespace MML;

void Show_Real_Function_Examples()
{
    std::cout << "\n=== Real Function Visualization Examples ===\n\n";

    // Example 1: Simple trigonometric function
    std::cout << "1. Visualizing f(x) = sin(x) * cos(2x)\n";
    RealFunction sinCosFunc{[](Real x) { return std::sin(x) * std::cos(2*x); }};
    Visualizer::VisualizeRealFunction(sinCosFunc, "sin(x)*cos(2x)", 
                                      -2*Constants::PI, 2*Constants::PI, 500, 
                                      "viz_real_func_sincos.txt");

    // Example 2: Polynomial with interesting features
    std::cout << "2. Visualizing f(x) = x^3 - 3x^2 - x + 3 (cubic with roots at -1, 1, 3)\n";
    RealFunction cubicFunc{[](Real x) { return x*x*x - 3*x*x - x + 3; }};
    Visualizer::VisualizeRealFunction(cubicFunc, "x^3 - 3x^2 - x + 3", 
                                      -2.0, 4.0, 300, 
                                      "viz_real_func_cubic.txt");

    // Example 3: Function with singularity (avoiding the singular point)
    std::cout << "3. Visualizing f(x) = sin(x)/(x-2) with singularity at x=2\n";
    RealFunction singularFunc{[](Real x) { 
        Real denom = x - 2.0;
        if (std::abs(denom) < 1e-10) return Real(0.0);
        return std::sin(x) / denom; 
    }};
    Visualizer::VisualizeRealFunction(singularFunc, "sin(x)/(x-2)", 
                                      -5.0, 7.0, 500, 
                                      "viz_real_func_singular.txt");

    // Example 4: Gaussian (bell curve)
    std::cout << "4. Visualizing Gaussian: f(x) = exp(-x^2/2)\n";
    RealFunction gaussianFunc{[](Real x) { return std::exp(-x*x/2.0); }};
    Visualizer::VisualizeRealFunction(gaussianFunc, "Gaussian exp(-x^2/2)", 
                                      -4.0, 4.0, 200, 
                                      "viz_real_func_gaussian.txt");

    // Example 5: Derivative visualization
    std::cout << "5. Visualizing derivative of sin(x)\n";
    RealFunction sinFunc{[](Real x) { return std::sin(x); }};
    RealFuncDerived4 sinDeriv(sinFunc);
    Visualizer::VisualizeRealFunction(sinDeriv, "d/dx sin(x) = cos(x)", 
                                      -2*Constants::PI, 2*Constants::PI, 300, 
                                      "viz_real_func_derivative.txt");

    std::cout << "\nReal function visualization examples complete!\n";
}
