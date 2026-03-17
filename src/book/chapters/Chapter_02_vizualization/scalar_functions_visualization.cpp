///////////////////////////////////////////////////////////////////////////////////////////
/// @file scalar_functions_visualization.cpp
/// @brief Chapter 02: Scalar Function Visualization Demos
/// @details Demonstrates visualization of scalar functions:
///          - 2D scalar functions f(x,y) → R displayed as surfaces
///          - 3D scalar functions f(x,y,z) → R displayed as isosurfaces/volumes
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/VectorN.h"
#include "mml/core/Fields.h"
#include "mml/tools/Visualizer.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 2D Scalar Functions (Surfaces)
/// @details Scalar functions f: R² → R visualized as 3D surfaces z = f(x,y)
///          Shows mathematical surfaces with interesting geometric properties.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_ScalarFunction2D()
{
    std::cout << "\n=== Demo: 2D Scalar Functions (Surfaces) ===\n";
    
    //-------------------------------------------------------------------------
    // Monkey Saddle - a surface with three "legs" meeting at the origin
    // Unlike a regular saddle (2 up, 2 down), monkey saddle has 3 alternating
    // regions - perfect for a monkey with two legs and a tail!
    //-------------------------------------------------------------------------
    ScalarFunction<2> monkey_saddle{ [](const VectorN<Real, 2>& v) {
        Real x = v[0], y = v[1];
        return x*x*x - 3.0 * x * y*y;  // z = x³ - 3xy²
    }};
    
    Visualizer::VisualizeScalarFunc2DCartesian(
        monkey_saddle, "Monkey Saddle: z = x³ - 3xy²",
        -3.0, 3.0, 50, -3.0, 3.0, 50,
        "ch02_monkey_saddle.mml");
    
    std::cout << "   Visualized: Monkey Saddle - 3 alternating up/down regions\n";
    
    //-------------------------------------------------------------------------
    // 2D Sinc Function (Mexican Hat / Sombrero)
    // The radial sinc function creates beautiful concentric ripples
    // Central peak with decaying oscillations - appears in optics, signal processing
    //-------------------------------------------------------------------------
    ScalarFunction<2> sinc_2d{ [](const VectorN<Real, 2>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        if (r < 1e-10) return Real(50.0);  // sinc(0) = 1, scaled by 50
        return 50.0 * std::sin(r) / r;     // Scale by 50 for better visualization
    }};
    
    Visualizer::VisualizeScalarFunc2DCartesian(
        sinc_2d, "2D Sinc (Sombrero): z = 50·sin(r)/r",
        -15.0, 15.0, 80, -15.0, 15.0, 80,
        "ch02_sinc_2d.mml");
    
    std::cout << "   Visualized: 2D Sinc / Mexican Hat - concentric ripples\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 3D Scalar Functions (Volumetric/Isosurfaces)
/// @details Scalar functions f: R³ → R visualized as isosurfaces or volumes.
///          Shows famous surfaces from differential geometry and materials science.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_ScalarFunction3D()
{
    std::cout << "\n=== Demo: 3D Scalar Functions (Isosurfaces) ===\n";
    
    //-------------------------------------------------------------------------
    // Gyroid - a triply periodic minimal surface
    // Famous in materials science, 3D printing, and architecture
    // Has zero mean curvature everywhere - soap film-like geometry
    // Discovered by Alan Schoen in 1970
    //-------------------------------------------------------------------------
    ScalarFunction<3> gyroid{ [](const VectorN<Real, 3>& v) {
        Real scale = Constants::PI / 2.0;  // One full period
        Real x = v[0] * scale;
        Real y = v[1] * scale;
        Real z = v[2] * scale;
        return std::sin(x) * std::cos(y) + 
               std::sin(y) * std::cos(z) + 
               std::sin(z) * std::cos(x);
    }};
    
    Visualizer::VisualizeScalarFunc3DCartesian(
        gyroid, "Gyroid (Triply Periodic Minimal Surface)",
        -2.0, 2.0, 40,
        -2.0, 2.0, 40,
        -2.0, 2.0, 40,
        "ch02_gyroid.mml");
    
    std::cout << "   Visualized: Gyroid - minimal surface used in 3D printing\n";
    
    //-------------------------------------------------------------------------
    // Torus Signed Distance Function (SDF)
    // The zero-level set defines a perfect torus surface
    // Negative inside, positive outside - fundamental in computer graphics
    //-------------------------------------------------------------------------
    ScalarFunction<3> torus_sdf{ [](const VectorN<Real, 3>& v) {
        const Real R = 0.7;  // Major radius (center of tube to center of torus)
        const Real r = 0.25; // Minor radius (tube radius)
        
        // Distance to torus: d = ||(||(x,y)|| - R, z)|| - r
        Real xy_dist = std::sqrt(v[0]*v[0] + v[1]*v[1]) - R;
        return std::sqrt(xy_dist*xy_dist + v[2]*v[2]) - r;
    }};
    
    Visualizer::VisualizeScalarFunc3DCartesian(
        torus_sdf, "Torus (Signed Distance Function)",
        -1.2, 1.2, 40,
        -1.2, 1.2, 40,
        -0.6, 0.6, 25,
        "ch02_torus_sdf.mml");
    
    std::cout << "   Visualized: Torus SDF - zero-level set is the torus surface\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Entry point for Scalar Function visualization demos
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_ScalarFunctionVisualization()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         CHAPTER 02: Scalar Function Visualization                 ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    Chapter02_Demo_ScalarFunction2D();
    Chapter02_Demo_ScalarFunction3D();
}
