/**
 * @file test_linux_realfunc_viz.cpp
 * @brief Test Linux RealFunction visualization with FLTK and Qt backends
 * 
 * This example demonstrates the new cross-platform visualization support
 * for RealFunction on Linux with both FLTK and Qt backends.
 * 
 * Usage:
 *   1. Default (FLTK):
 *      ./test_linux_realfunc_viz
 * 
 *   2. Using Qt backend (environment variable):
 *      export MML_LINUX_VISUALIZER=Qt
 *      ./test_linux_realfunc_viz
 * 
 *   3. Using Qt backend (programmatic):
 *      Uncomment the SetLinuxVisualizerBackend line below
 */

#include "MML.h"

using namespace MML;

int main()
{
    std::cout << "=== Linux RealFunction Visualization Test ===" << std::endl;
    
    // OPTION 1: Programmatically select Qt backend (uncomment to use)
    // SetLinuxVisualizerBackend(LinuxVisualizerBackend::Qt);
    
    // OPTION 2: Use environment variable (see usage above)
    // This takes precedence over programmatic setting
    
    // Show which backend is active
    auto backend = GetVisualizerBackend();
    std::cout << "Active visualizer backend: " 
              << (backend == LinuxVisualizerBackend::FLTK ? "FLTK" : "Qt") 
              << std::endl;
    
    // Create a simple sine function
    auto sineFunc = RealFunctionFromStdFunc([](Real x) { 
        return std::sin(x); 
    });
    
    std::cout << "\n1. Visualizing single function: sin(x)" << std::endl;
    Visualizer::VisualizeRealFunction(sineFunc, "Sine Function", 
                                       -Constants::PI, Constants::PI, 
                                       100, "sine_test.txt");
    
    // Create multiple functions for comparison
    auto cosFunc = RealFunctionFromStdFunc([](Real x) { 
        return std::cos(x); 
    });
    auto tanFunc = RealFunctionFromStdFunc([](Real x) { 
        return std::tan(x); 
    });
    
    std::cout << "2. Visualizing multiple functions: sin, cos, tan" << std::endl;
    std::vector<IRealFunction*> funcs = {&sineFunc, &cosFunc, &tanFunc};
    std::vector<std::string> legends = {"sin(x)", "cos(x)", "tan(x)"};
    
    Visualizer::VisualizeMultiRealFunction(funcs, "Trigonometric Functions",
                                            legends, -Constants::PI, Constants::PI,
                                            100, "trig_multi_test.txt");
    
    // Test with polynomial
    auto polyFunc = RealFunctionFromStdFunc([](Real x) { 
        return x*x*x - 3*x*x + 2*x + 1; 
    });
    
    std::cout << "3. Visualizing polynomial: x³ - 3x² + 2x + 1" << std::endl;
    Visualizer::VisualizeRealFunction(polyFunc, "Cubic Polynomial", 
                                       -2.0, 4.0, 150, "poly_test.txt");
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    std::cout << "Visualizer executables called:" << std::endl;
    std::cout << "  " << GetRealFuncVisualizerPath() << std::endl;
    
    return 0;
}
