#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Function.h"
#include "base/VectorN.h"
#include "core/Fields.h"
#include "core/Curves.h"
#include "tools/Visualizer.h"
#include "tools/Serializer.h"
#endif

#include <iostream>
#include <vector>

using namespace MML;

void TestRealFunction()
{
    std::cout << "\n=== Testing RealFunction Visualization ===\n";
    RealFunction f{[](Real x) { return sin(x) * cos(x/2); }};
    
    // Test FLTK
    SetVisualizerBackend(VisualizerBackend::FLTK);
    std::cout << "FLTK: ";
    Visualizer::VisualizeRealFunction(f, "RealFunction - FLTK", -10.0, 10.0, 200, "test_realfunc_fltk.txt");
    std::cout << "✓\n";
    
    // Test Qt
    SetVisualizerBackend(VisualizerBackend::Qt);
    std::cout << "Qt:   ";
    Visualizer::VisualizeRealFunction(f, "RealFunction - Qt", -10.0, 10.0, 200, "test_realfunc_qt.txt");
    std::cout << "✓\n";
}

void TestScalarFunction2D()
{
    std::cout << "\n=== Testing ScalarFunction2D Visualization (Qt) ===\n";
    SetVisualizerBackend(VisualizerBackend::Qt);
    
    ScalarFunction<2> func{[](const VectorN<Real, 2> &x) -> Real { 
        return REAL(sin(x[0]) * cos(x[1]) * exp(-(x[0]*x[0] + x[1]*x[1])/REAL(10.0))); 
    }};
    
    std::cout << "Qt:   ";
    Visualizer::VisualizeScalarFunc2DCartesian(func, "ScalarFunction2D - Qt", 
                                                -5.0, 5.0, 30, -5.0, 5.0, 30, 
                                                "test_scalar2d_qt.txt");
    std::cout << "✓\n";
}

void TestVectorField2D()
{
    std::cout << "\n=== Testing VectorField2D Visualization (FLTK) ===\n";
    SetVisualizerBackend(VisualizerBackend::FLTK);
    
    // Simple rotation field
    VectorFunction<2> field{[](const VectorN<Real, 2> &pos) {
        VectorN<Real, 2> result;
        result[0] = -pos[1];  // x component
        result[1] =  pos[0];  // y component
        return result;
    }};
    
    std::cout << "FLTK: ";
    Visualizer::VisualizeVectorField2DCartesian(field, "VectorField2D - FLTK", 
                                                 -5.0, 5.0, 15, -5.0, 5.0, 15, 
                                                 "test_vectorfield2d_fltk.txt");
    std::cout << "✓\n";
}

void TestVectorField3D()
{
    std::cout << "\n=== Testing VectorField3D Visualization (Qt) ===\n";
    SetVisualizerBackend(VisualizerBackend::Qt);
    
    // Simple 3D field
    VectorFunction<3> field{[](const VectorN<Real, 3> &pos) {
        VectorN<Real, 3> result;
        result[0] = -pos[1];
        result[1] =  pos[0];
        result[2] = -pos[2] * 0.1;
        return result;
    }};
    
    std::cout << "Qt:   ";
    Visualizer::VisualizeVectorField3DCartesian(field, "VectorField3D - Qt", 
                                                 -3.0, 3.0, 8, -3.0, 3.0, 8, -3.0, 3.0, 8,
                                                 "test_vectorfield3d_qt.txt");
    std::cout << "✓\n";
}

void TestParametricCurve2D()
{
    std::cout << "\n=== Testing ParametricCurve2D Visualization (FLTK) ===\n";
    SetVisualizerBackend(VisualizerBackend::FLTK);
    
    // Use a simple circle curve
    Curves::Circle2DCurve curve(2.0);
    
    std::cout << "FLTK: ";
    Visualizer::VisualizeParamCurve2D(curve, "ParametricCurve2D - FLTK", 
                                      0.0, 2*Constants::PI, 200, 
                                      "test_paramcurve2d_fltk.txt");
    std::cout << "✓\n";
}

void TestParametricCurve3D()
{
    std::cout << "\n=== Testing ParametricCurve3D Visualization (Qt) ===\n";
    SetVisualizerBackend(VisualizerBackend::Qt);
    
    // 3D helix
    Curves::HelixCurve curve(2.0, 0.5);
    
    std::cout << "Qt:   ";
    Visualizer::VisualizeParamCurve3D(curve, "ParametricCurve3D - Qt", 
                                      0.0, 4*Constants::PI, 300, 
                                      "test_paramcurve3d_qt.txt");
    std::cout << "✓\n";
}

void TestParticleSimulation2D()
{
    std::cout << "\n=== Testing ParticleSimulation2D Visualization (FLTK) ===\n";
    SetVisualizerBackend(VisualizerBackend::FLTK);
    
    // Create a simple particle simulation file
    std::string filename = "results/test_particles2d.txt";
    {
        std::ofstream out(filename);
        out << "# 2D Particle Simulation Test Data\n";
        out << "# Format: time x1 y1 x2 y2 ...\n";
        
        // Simulate 3 particles moving in circles
        for (int frame = 0; frame < 50; frame++) {
            Real t = frame * 0.1;
            out << t;
            for (int p = 0; p < 3; p++) {
                Real angle = t + p * 2.0 * Constants::PI / 3.0;
                Real radius = 2.0 + p * 0.5;
                out << " " << radius * cos(angle) << " " << radius * sin(angle);
            }
            out << "\n";
        }
    }
    
    std::cout << "FLTK: ";
    Visualizer::VisualizeParticleSimulation2D(filename);
    std::cout << "✓\n";
}

void TestParticleSimulation3D()
{
    std::cout << "\n=== Testing ParticleSimulation3D Visualization (Qt) ===\n";
    SetVisualizerBackend(VisualizerBackend::Qt);
    
    // Create a simple 3D particle simulation file
    std::string filename = "results/test_particles3d.txt";
    {
        std::ofstream out(filename);
        out << "# 3D Particle Simulation Test Data\n";
        out << "# Format: time x1 y1 z1 x2 y2 z2 ...\n";
        
        // Simulate 2 particles in 3D helical paths
        for (int frame = 0; frame < 50; frame++) {
            Real t = frame * 0.1;
            out << t;
            for (int p = 0; p < 2; p++) {
                Real angle = t + p * Constants::PI;
                Real radius = 2.0 + p * 0.5;
                out << " " << radius * cos(angle) 
                    << " " << radius * sin(angle) 
                    << " " << t;
            }
            out << "\n";
        }
    }
    
    std::cout << "Qt:   ";
    Visualizer::VisualizeParticleSimulation3D(filename);
    std::cout << "✓\n";
}

int main()
{
    std::cout << "\n***********************************************************\n";
    std::cout << "******  Testing ALL Linux Visualizers  *******\n";
    std::cout << "***********************************************************\n";

    try {
        TestRealFunction();
        TestScalarFunction2D();
        TestVectorField2D();
        TestVectorField3D();
        TestParametricCurve2D();
        TestParametricCurve3D();
        TestParticleSimulation2D();
        TestParticleSimulation3D();
        
        std::cout << "\n***********************************************************\n";
        std::cout << "All visualizer tests completed successfully! ✓\n";
        std::cout << "***********************************************************\n";
        
        return 0;
    } 
    catch (const std::exception& e) {
        std::cout << "\n✗ Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}
