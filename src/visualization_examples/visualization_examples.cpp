/**
 * @file visualization_examples.cpp
 * @brief Main entry point for MML visualization examples
 * 
 * This file demonstrates all visualization capabilities of the MinimalMathLibrary.
 * Each visualization type has its own dedicated file:
 * 
 * - show_real_function.cpp      : Single real-valued functions y = f(x)
 * - show_multi_real_function.cpp : Multiple functions on same plot
 * - show_scalar_function.cpp     : 2D scalar fields as 3D surfaces
 * - show_vector_field_2d.cpp     : 2D vector fields
 * - show_vector_field_3d.cpp     : 3D vector fields  
 * - show_parametric_curve_2d.cpp : Parametric curves in 2D
 * - show_parametric_curve_3d.cpp : Parametric curves in 3D
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
void Show_Vector_Field_2D_Examples();
void Show_Vector_Field_3D_Examples();
void Show_Parametric_Curve_2D_Examples();
void Show_Parametric_Curve_3D_Examples();
void Show_Particle_Visualizer_2D_Examples();
void Show_Particle_Visualizer_3D_Examples();

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
    std::cout << "  curve_2d        - Parametric curves in 2D\n";
    std::cout << "  curve_3d        - Parametric curves in 3D\n";
    std::cout << "  field_2d        - 2D vector fields\n";
    std::cout << "  field_3d        - 3D vector fields\n";
    std::cout << "  particle_2d     - 2D particle simulations\n";
    std::cout << "  particle_3d     - 3D particle simulations\n";
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
    std::cout << "  Backend: " << (GetVisualizerBackend() == VisualizerBackend::FLTK ? "FLTK" : "Qt") << "\n";
    
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
    else if (strcmp(viz_type, "curve_2d") == 0) {
        Show_Parametric_Curve_2D_Examples();
        return true;
    }
    else if (strcmp(viz_type, "curve_3d") == 0) {
        Show_Parametric_Curve_3D_Examples();
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
    else if (strcmp(viz_type, "all") == 0) {
        Show_Real_Function_Examples();
        Show_Multi_Real_Function_Examples();
        Show_Scalar_Function_Examples();
        Show_Parametric_Curve_2D_Examples();
        Show_Parametric_Curve_3D_Examples();
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
    // VECTOR FIELDS
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Vector_Field_2D_Examples();         // 2D vector fields
    Show_Vector_Field_3D_Examples();         // 3D vector fields
    
    // ═══════════════════════════════════════════════════════════════════
    // PARAMETRIC CURVES
    // ═══════════════════════════════════════════════════════════════════
    
    // Show_Parametric_Curve_2D_Examples();     // 2D parametric curves
    Show_Parametric_Curve_3D_Examples();     // 3D parametric curves
    
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
 * @brief Main entry point for the visualization examples application
 */
int main(int argc, char* argv[])
{
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