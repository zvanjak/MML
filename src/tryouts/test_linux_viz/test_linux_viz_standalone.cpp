#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Function.h"
#include "core/FunctionHelpers.h"
#include "tools/Visualizer.h"
#endif

#include <iostream>

using namespace MML;

int main()
{
    std::cout << "\n***********************************************************\n";
    std::cout << "**********  Testing Linux RealFunction Visualization  ********\n";
    std::cout << "***********************************************************\n\n";

    std::cout << "Current backend: " << (GetVisualizerBackend() == VisualizerBackend::FLTK ? "FLTK" : "Qt") << "\n\n";

    // Test 1: Simple function
    std::cout << "Test 1: Visualizing single function f(x) = sin(x)*(x-3)*(x+5)/sqrt(|2-x|)\n";
    RealFunction f1{[](Real x) { return sin(x) * (x-3)*(x+5) / sqrt(std::abs(2 - x)); } };
    
    try {
        Visualizer::VisualizeRealFunction(f1, "Linux Test - Simple Function", -10.0, 10.0, 500, "linux_test_simple.txt");
        std::cout << "✓ Test 1 completed successfully\n\n";
    } catch (const std::exception& e) {
        std::cout << "✗ Test 1 failed: " << e.what() << "\n\n";
        return 1;
    }

    // Test 2: Derivative visualization
    std::cout << "Test 2: Visualizing derivative of f(x)\n";
    RealFuncDerived4 f2(f1);
    
    try {
        Visualizer::VisualizeRealFunction(f2, "Linux Test - Derivative", -10.0, 10.0, 500, "linux_test_derivative.txt");
        std::cout << "✓ Test 2 completed successfully\n\n";
    } catch (const std::exception& e) {
        std::cout << "✗ Test 2 failed: " << e.what() << "\n\n";
        return 1;
    }

    // Test 3: Multiple functions
    std::cout << "Test 3: Visualizing both functions together\n";
    
    try {
        Visualizer::VisualizeMultiRealFunction({&f1, &f2}, "Linux Test - Multi Function", 
                                              {"Original", "Derivative"},
                                              -10.0, 10.0, 500, "linux_test_multi.txt");
        std::cout << "✓ Test 3 completed successfully\n\n";
    } catch (const std::exception& e) {
        std::cout << "✗ Test 3 failed: " << e.what() << "\n\n";
        return 1;
    }

    std::cout << "***********************************************************\n";
    std::cout << "All Linux RealFunction visualization tests passed! ✓\n";
    std::cout << "Backend used: " << (GetVisualizerBackend() == VisualizerBackend::FLTK ? "FLTK" : "Qt") << "\n";
    std::cout << "***********************************************************\n";

    return 0;
}
