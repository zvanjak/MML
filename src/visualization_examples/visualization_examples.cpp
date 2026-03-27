/**
 * @file visualization_examples.cpp
 * @brief Main entry point for MML visualization examples
 * 
 * This file demonstrates all visualization capabilities of the MinimalMathLibrary.
 * Each visualization type has its own dedicated file:
 * 
 * - show_real_function.cpp       : Single real-valued functions y = f(x)
 * - show_multi_real_function.cpp : Multiple functions on same plot
 * - show_scalar_function.cpp     : 2D scalar fields as 3D surfaces
 * - show_scalar_function_3d.cpp  : 3D scalar fields (volumetric/isosurface)
 * - show_vector_field_2d.cpp     : 2D vector fields
 * - show_vector_field_3d.cpp     : 3D vector fields  
 * - show_parametric_curve_2d.cpp : Parametric curves in 2D
 * - show_parametric_curve_3d.cpp : Parametric curves in 3D
 * - show_parametric_surface.cpp  : Parametric surfaces in 3D
 * - show_particle_visualizer_2d.cpp : 2D particle simulations
 * - show_particle_visualizer_3d.cpp : 3D particle simulations
 * 
 * The visualization system is cross-platform:
 * - Windows: Uses Windows-native visualizers from MML_Visualizers project
 * - Linux/macOS: Uses FLTK or Qt-based visualizers (platform-dependent)
 * 
 * All paths are resolved using the new visualization infrastructure:
 * - GetResultFilesPath()        : Path for output data files
 * - GetRealFuncVisualizerPath() : Path to real function visualizer
 * - GetSurfaceVisualizerPath()  : Path to surface visualizer
 * - etc.
 * 
 * Command-line usage:
 *   MML_VisualizationApp [visualization_type]
 * 
 * Where visualization_type is one of:
 *   real_function   - Real-valued functions y = f(x)
 *   multi_function  - Multiple real functions on same plot
 *   scalar_function - 2D scalar fields as 3D surfaces
 *   curve_2d        - Parametric curves in 2D
 *   curve_3d        - Parametric curves in 3D
 *   field_2d        - 2D vector fields
 *   field_3d        - 3D vector fields
 *   particle_2d     - 2D particle simulations
 *   particle_3d     - 3D particle simulations
 * 
 * If no argument provided, runs demo with selected visualizations.
 * 
 * @author MinimalMathLibrary
 * @date December 2025
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "tools/Visualizer.h"
#endif

#include <cstring>

using namespace MML;

// Forward declarations for all visualization functions
void Show_Real_Function_Examples();
void Show_Multi_Real_Function_Examples();
void Show_Scalar_Function_Examples();
void Show_Scalar_Function_3D_Examples();
void Show_Vector_Field_2D_Examples();
void Show_Vector_Field_3D_Examples();
void Show_Parametric_Curve_2D_Examples();
void Show_Parametric_Curve_3D_Examples();
void Show_Parametric_Surface_Examples();
void Show_Particle_Visualizer_2D_Examples();
void Show_Particle_Visualizer_3D_Examples();
void Show_Readme_Files();
void Show_Backend_Switching_Demo();

/**
 * @brief Print usage information
 */
void Print_Usage(const char* program_name)
{
    std::cout << "\nUsage: " << program_name << " [visualization_type]\n\n";
    std::cout << "Available visualization types:\n";
    std::cout << "  real_function   - Real-valued functions y = f(x)\n";
    std::cout << "  multi_function  - Multiple real functions on same plot\n";
    std::cout << "  scalar_function - 2D scalar fields as 3D surfaces\n";
    std::cout << "  scalar_3d       - 3D scalar fields (volumetric/isosurface)\n";
    std::cout << "  curve_2d        - Parametric curves in 2D\n";
    std::cout << "  curve_3d        - Parametric curves in 3D\n";
    std::cout << "  surface         - Parametric surfaces in 3D\n";
    std::cout << "  field_2d        - 2D vector fields\n";
    std::cout << "  field_3d        - 3D vector fields\n";
    std::cout << "  particle_2d     - 2D particle simulations\n";
    std::cout << "  particle_3d     - 3D particle simulations\n";
    std::cout << "  readme          - One sample of each visualization type\n";
    std::cout << "  backend         - Runtime backend switching demo\n";
    std::cout << "  all             - Run all visualizations\n";
    std::cout << "\nIf no argument is provided, runs the default demo.\n\n";
}

/**
 * @brief Print information about the visualization environment
 */
void Print_Visualization_Info()
{
    std::cout << "\n";
    std::cout << "+================================================================+\n";
    std::cout << "|        MinimalMathLibrary - Visualization Examples            |\n";
    std::cout << "+================================================================+\n";
    std::cout << "| Cross-platform visualization infrastructure                   |\n";
    std::cout << "+================================================================+\n\n";

    std::cout << "Configuration:\n";
    std::cout << "  Results folder: " << GetResultFilesPath() << "\n";
    
    auto backend = GetVisualizerBackend();
    std::string backendName;
    switch (backend) {
        case VisualizerBackend::WPF:  backendName = "WPF"; break;
        case VisualizerBackend::Qt:   backendName = "Qt"; break;
        case VisualizerBackend::FLTK: backendName = "FLTK"; break;
        case VisualizerBackend::Auto: backendName = "Auto"; break;
    }
    std::cout << "  Backend: " << backendName << "\n";
    
    // Show the actual path to the visualizer
    std::cout << "  RealFunc Visualizer: " << GetRealFuncVisualizerPath() << "\n";
    
#ifdef _WIN32
    std::cout << "  Platform: Windows\n";
#elif defined(__APPLE__)
    std::cout << "  Platform: macOS\n";
#else
    std::cout << "  Platform: Linux\n";
#endif

    std::cout << "\n";
}

/**
 * @brief Run a specific visualization based on the type string
 * @param viz_type The visualization type to run
 * @return true if valid type, false otherwise
 */
bool Run_Visualization(const char* viz_type)
{
    if (strcmp(viz_type, "real_function") == 0) {
        Show_Real_Function_Examples();
        return true;
    }
    else if (strcmp(viz_type, "multi_function") == 0) {
        Show_Multi_Real_Function_Examples();
        return true;
    }
    else if (strcmp(viz_type, "scalar_function") == 0) {
        Show_Scalar_Function_Examples();
        return true;
    }
    else if (strcmp(viz_type, "scalar_3d") == 0) {
        Show_Scalar_Function_3D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "curve_2d") == 0) {
        Show_Parametric_Curve_2D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "curve_3d") == 0) {
        Show_Parametric_Curve_3D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "surface") == 0) {
        Show_Parametric_Surface_Examples();
        return true;
    }
    else if (strcmp(viz_type, "field_2d") == 0) {
        Show_Vector_Field_2D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "field_3d") == 0) {
        Show_Vector_Field_3D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "particle_2d") == 0) {
        Show_Particle_Visualizer_2D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "particle_3d") == 0) {
        Show_Particle_Visualizer_3D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "readme") == 0) {
        Show_Readme_Files();
        return true;
    }
    else if (strcmp(viz_type, "backend") == 0) {
        Show_Backend_Switching_Demo();
        return true;
    }
    else if (strcmp(viz_type, "all") == 0) {
        Show_Real_Function_Examples();
        Show_Multi_Real_Function_Examples();
        Show_Scalar_Function_Examples();
        Show_Scalar_Function_3D_Examples();
        Show_Parametric_Curve_2D_Examples();
        Show_Parametric_Curve_3D_Examples();
        Show_Parametric_Surface_Examples();
        Show_Vector_Field_2D_Examples();
        Show_Vector_Field_3D_Examples();
        Show_Particle_Visualizer_2D_Examples();
        Show_Particle_Visualizer_3D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "--help") == 0 || strcmp(viz_type, "-h") == 0) {
        return false;  // Will print usage
    }
    
    std::cout << "Unknown visualization type: " << viz_type << "\n";
    return false;
}

/**
 * @brief Default demo - runs selected visualizations
 */
void Demo_Visualization_Examples()
{
    // ═══════════════════════════════════════════════════════════════════
    // REAL FUNCTIONS (1D)
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Real_Function_Examples();           // Single functions
    // Show_Multi_Real_Function_Examples();     // Multiple functions together
    
    // ═══════════════════════════════════════════════════════════════════
    // SCALAR FUNCTIONS (2D → Surface)
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Scalar_Function_Examples();         // 2D scalar fields as surfaces
    
    // ═══════════════════════════════════════════════════════════════════
    // SCALAR FUNCTIONS (3D → Isosurfaces/Volumetric)
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Scalar_Function_3D_Examples();      // 3D scalar fields (Gaussian, Gyroid, SDFs)
    
    // ═══════════════════════════════════════════════════════════════════
    // PARAMETRIC SURFACES
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Parametric_Surface_Examples();      // 3D parametric surfaces (Sphere, Torus, Möbius)
    
    // ═══════════════════════════════════════════════════════════════════
    // VECTOR FIELDS
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Vector_Field_2D_Examples();         // 2D vector fields
    // Show_Vector_Field_3D_Examples();         // 3D vector fields
    
    // ═══════════════════════════════════════════════════════════════════
    // PARAMETRIC CURVES
    // ═══════════════════════════════════════════════════════════════════
    
    Show_Parametric_Curve_2D_Examples();     // 2D parametric curves
    // Show_Parametric_Curve_3D_Examples();     // 3D parametric curves
    
    // ═══════════════════════════════════════════════════════════════════

    // PARTICLE SIMULATIONS
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Particle_Visualizer_2D_Examples();  // 2D particle simulations
    // Show_Particle_Visualizer_3D_Examples();  // 3D particle simulations
    
    // ═══════════════════════════════════════════════════════════════════
    // QUICK TEST - Uncomment ONE of these to test specific visualizer
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Real_Function_Examples();  // Default: start with simplest visualizer
}

/**
 * @brief Demonstrate runtime backend switching
 * 
 * Shows how to switch between visualization backends (FLTK, Qt, WPF)
 * at runtime using SetVisualizerBackend() / GetVisualizerBackend().
 */
void Show_Backend_Switching_Demo()
{
    std::cout << "\n***********************************************************\n";
    std::cout << "******  Runtime Backend Switching Demo  *******\n";
    std::cout << "***********************************************************\n\n";

    RealFunction f{[](Real x) { return sin(x) * cos(x / 2); }};

    auto backendName = [](VisualizerBackend b) -> const char* {
        switch (b) {
            case VisualizerBackend::WPF:  return "WPF";
            case VisualizerBackend::Qt:   return "Qt";
            case VisualizerBackend::FLTK: return "FLTK";
            case VisualizerBackend::Auto: return "Auto";
        }
        return "Unknown";
    };

    // Save original backend to restore later
    auto original = GetVisualizerBackend();

    // Step 1: Visualize with current backend
    std::cout << "Step 1: Using current backend (" << backendName(original) << ")\n";
    Visualizer::VisualizeRealFunction(f, "Current Backend", -10.0, 10.0, 200, "backend_demo_current.mml");

    // Step 2: Switch to Qt
    std::cout << "\nStep 2: Switching to Qt backend\n";
    SetVisualizerBackend(VisualizerBackend::Qt);
    std::cout << "  Backend is now: " << backendName(GetVisualizerBackend()) << "\n";
    Visualizer::VisualizeRealFunction(f, "Qt Backend", -10.0, 10.0, 200, "backend_demo_qt.mml");

    // Step 3: Switch to FLTK
    std::cout << "\nStep 3: Switching to FLTK backend\n";
    SetVisualizerBackend(VisualizerBackend::FLTK);
    std::cout << "  Backend is now: " << backendName(GetVisualizerBackend()) << "\n";
    Visualizer::VisualizeRealFunction(f, "FLTK Backend", -10.0, 10.0, 200, "backend_demo_fltk.mml");

    // Restore original backend
    SetVisualizerBackend(original);
    std::cout << "\nRestored original backend: " << backendName(original) << "\n";

    std::cout << "\n***********************************************************\n";
    std::cout << "Backend switching demo complete!\n";
    std::cout << "***********************************************************\n";
}

/**
 * @brief Main entry point for the visualization examples application
 */
int main(int argc, char* argv[])
{
    // Uncomment ONE of these to override the default WPF backend on Windows:
    // SetVisualizerBackend(VisualizerBackend::Qt);
    // SetVisualizerBackend(VisualizerBackend::FLTK);
    SetVisualizerBackend(VisualizerBackend::WPF);
    
    Print_Visualization_Info();
    
    if (argc > 1) {
        // Run specific visualization from command line
        if (!Run_Visualization(argv[1])) {
            Print_Usage(argv[0]);
            return 1;
        }
    }
    else {
        // Run default demo
        std::cout << "Running default demo. Use --help for options.\n\n";
        Demo_Visualization_Examples();
    }
    
    std::cout << "\n================================================================\n";
    std::cout << "Visualization examples complete!\n";
    std::cout << "================================================================\n\n";
    
    return 0;
}
