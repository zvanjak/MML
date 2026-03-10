/**
 * @file show_parametric_surface.cpp
 * @brief Demonstrates visualization of 3D parametric surfaces
 * 
 * Uses the MML Visualizer infrastructure to display parametric surfaces in 3D.
 * Parametric surfaces are defined as mappings (u, w) → (x, y, z).
 * Uses predefined surfaces from MML::Surfaces namespace.
 * Cross-platform: works on Windows, Linux, and macOS.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorN.h"
#include "core/Surfaces.h"
#include "tools/Visualizer.h"
#endif

using namespace MML;
using namespace MML::Surfaces;

///////////////////////////////////////////////////////////////////////////////
// Main Example Function
///////////////////////////////////////////////////////////////////////////////

void Show_Parametric_Surface_Examples()
{
    std::cout << "\n=== 3D Parametric Surface Visualization Examples ===\n\n";

    // Example 1: Sphere
    std::cout << "1. Sphere (radius = 50)\n";
    std::cout << "   Classic parametric: x = r*sin(u)*cos(w), y = r*sin(u)*sin(w), z = r*cos(u)\n";
    Sphere sphere(50.0);
    Visualizer::VisualizeParametricSurface(sphere, "Sphere",
                                           30, 50,
                                           "viz_surface_sphere.mml");

    // Example 2: Torus
    std::cout << "2. Torus (R=40, r=15)\n";
    std::cout << "   Donut shape: (R + r*cos(w))*cos(u), (R + r*cos(w))*sin(u), r*sin(w)\n";
    Torus torus(40.0, 15.0);
    Visualizer::VisualizeParametricSurface(torus, "Torus (Donut)",
                                           50, 30,
                                           "viz_surface_torus.mml");

    // Example 3: Möbius Strip
    std::cout << "3. Mobius Strip\n";
    std::cout << "   Non-orientable surface with only one side!\n";
    MobiusStrip mobius(30.0);  // Scale = 30 for visibility
    Visualizer::VisualizeParametricSurface(mobius, "Mobius Strip",
                                           60, 15,
                                           "viz_surface_mobius.mml");

    // Example 4: Helicoid
    std::cout << "4. Helicoid (Spiral Ramp)\n";
    std::cout << "   Minimal surface - like a spiral staircase\n";
    Helicoid helicoid(5.0);  // pitch = 5
    Visualizer::VisualizeParametricSurface(helicoid, "Helicoid",
                                           25, 80,
                                           "viz_surface_helicoid.mml");

    // Example 5: Enneper Surface
    std::cout << "5. Enneper Surface\n";
    std::cout << "   Self-intersecting minimal surface with elegant shape\n";
    EnneperSurface enneper(15.0);  // scale = 15 for visibility
    Visualizer::VisualizeParametricSurface(enneper, "Enneper Surface",
                                           50, 50,
                                           "viz_surface_enneper.mml");

    std::cout << "\nParametric surface examples complete!\n";
}
