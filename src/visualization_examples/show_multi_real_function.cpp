/**
 * @file show_multi_real_function.cpp
 * @brief Demonstrates visualization of multiple Real functions together
 * 
 * Uses the MML Visualizer infrastructure to display multiple real-valued 
 * functions on the same plot for comparison.
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

void Show_Multi_Real_Function_Examples()
{
    std::cout << "\n=== Multi Real Function Visualization Examples ===\n\n";

    // Example 1: Trigonometric functions comparison
    std::cout << "1. Comparing sin(x), cos(x), and tan(x)\n";
    RealFunction sinFunc{[](Real x) { return std::sin(x); }};
    RealFunction cosFunc{[](Real x) { return std::cos(x); }};
    RealFunction tanFunc{[](Real x) { 
        Real c = std::cos(x);
        if (std::abs(c) < 1e-10) return Real(0.0);
        return std::sin(x) / c; 
    }};
    
    Visualizer::VisualizeMultiRealFunction(
        {&sinFunc, &cosFunc, &tanFunc}, 
        "Trigonometric Functions",
        {"sin(x)", "cos(x)", "tan(x)"},
        -Constants::PI, Constants::PI, 500, 
        "viz_multi_func_trig.txt");

    // Example 2: Function and its derivatives
    std::cout << "2. Function and its 1st and 2nd derivatives: x^3 - 2x^2 + x\n";
    RealFunction polyFunc{[](Real x) { return x*x*x - 2*x*x + x; }};
    RealFunction polyDeriv1{[](Real x) { return 3*x*x - 4*x + 1; }};  // Analytical first derivative
    RealFunction polyDeriv2{[](Real x) { return 6*x - 4; }};          // Analytical second derivative
    
    Visualizer::VisualizeMultiRealFunction(
        {&polyFunc, &polyDeriv1, &polyDeriv2}, 
        "Polynomial and Derivatives",
        {"f(x) = x^3 - 2x^2 + x", "f'(x) = 3x^2 - 4x + 1", "f''(x) = 6x - 4"},
        -1.0, 2.0, 300, 
        "viz_multi_func_poly_derivs.txt");

    // Example 3: Damped oscillations with different damping
    std::cout << "3. Damped oscillations with different damping coefficients\n";
    RealFunction dampedLow{[](Real x) { return std::exp(-0.1*x) * std::cos(2*x); }};
    RealFunction dampedMed{[](Real x) { return std::exp(-0.3*x) * std::cos(2*x); }};
    RealFunction dampedHigh{[](Real x) { return std::exp(-0.5*x) * std::cos(2*x); }};
    
    Visualizer::VisualizeMultiRealFunction(
        {&dampedLow, &dampedMed, &dampedHigh}, 
        "Damped Oscillations",
        {"Low damping (0.1)", "Medium damping (0.3)", "High damping (0.5)"},
        0.0, 20.0, 500, 
        "viz_multi_func_damped.txt");

    // Example 4: Bessel-like functions (approximations)
    std::cout << "4. Comparing different wave functions\n";
    RealFunction wave1{[](Real x) { return std::sin(x) / (x + 0.001); }};  // sinc-like
    RealFunction wave2{[](Real x) { return std::cos(x) * std::exp(-x*x/50); }};  // Gaussian modulated
    RealFunction wave3{[](Real x) { return std::sin(2*x) * std::cos(x/2); }};  // Beat pattern
    
    Visualizer::VisualizeMultiRealFunction(
        {&wave1, &wave2, &wave3}, 
        "Wave Functions",
        {"sinc-like", "Gaussian modulated", "Beat pattern"},
        -10.0, 10.0, 500, 
        "viz_multi_func_waves.txt");

    // Example 5: Exponential growth comparison
    std::cout << "5. Exponential growth rates comparison\n";
    RealFunction exp1{[](Real x) { return std::exp(0.5*x); }};
    RealFunction exp2{[](Real x) { return std::exp(x); }};
    RealFunction exp3{[](Real x) { return std::exp(1.5*x); }};
    
    Visualizer::VisualizeMultiRealFunction(
        {&exp1, &exp2, &exp3}, 
        "Exponential Growth",
        {"exp(0.5x)", "exp(x)", "exp(1.5x)"},
        0.0, 3.0, 200, 
        "viz_multi_func_exp.txt");

    std::cout << "\nMulti real function visualization examples complete!\n";
}
