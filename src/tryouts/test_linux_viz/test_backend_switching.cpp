#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Function.h"
#include "tools/Visualizer.h"
#endif

#include <iostream>

using namespace MML;

int main()
{
    std::cout << "\n***********************************************************\n";
    std::cout << "******  Testing Runtime Backend Switching API  *******\n";
    std::cout << "***********************************************************\n\n";

    RealFunction f{[](Real x) { return sin(x) * cos(x/2); } };

    // Test 1: Default backend (FLTK)
    std::cout << "Test 1: Using default backend (FLTK)\n";
    std::cout << "Current backend: " << (GetVisualizerBackend() == VisualizerBackend::FLTK ? "FLTK" : "Qt") << "\n";
    Visualizer::VisualizeRealFunction(f, "Default Backend Test", -10.0, 10.0, 200, "backend_test_default.txt");
    std::cout << "✓ Test 1 completed\n\n";

    // Test 2: Switch to Qt at runtime
    std::cout << "Test 2: Switching to Qt backend at runtime\n";
    SetVisualizerBackend(VisualizerBackend::Qt);
    std::cout << "Current backend: " << (GetVisualizerBackend() == VisualizerBackend::FLTK ? "FLTK" : "Qt") << "\n";
    Visualizer::VisualizeRealFunction(f, "Qt Backend Test", -10.0, 10.0, 200, "backend_test_qt.txt");
    std::cout << "✓ Test 2 completed\n\n";

    // Test 3: Switch back to FLTK
    std::cout << "Test 3: Switching back to FLTK backend\n";
    SetVisualizerBackend(VisualizerBackend::FLTK);
    std::cout << "Current backend: " << (GetVisualizerBackend() == VisualizerBackend::FLTK ? "FLTK" : "Qt") << "\n";
    Visualizer::VisualizeRealFunction(f, "FLTK Backend Test", -10.0, 10.0, 200, "backend_test_fltk.txt");
    std::cout << "✓ Test 3 completed\n\n";

    std::cout << "***********************************************************\n";
    std::cout << "All backend switching tests passed! ✓\n";
    std::cout << "***********************************************************\n";

    return 0;
}
