/**
 * @file show_scalar_function_3d.cpp
 * @brief Demonstrates visualization of 3D scalar functions (volumetric data)
 * 
 * Uses the MML Visualizer infrastructure to display 3D scalar functions
 * as isosurfaces or volume renderings. These are functions f: R³ → R.
 * Cross-platform: works on Windows, Linux, and macOS.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorN.h"
#include "core/Fields.h"
#include "tools/Visualizer.h"
#endif

using namespace MML;

void Show_Scalar_Function_3D_Examples()
{
    std::cout << "\n=== 3D Scalar Function (Volumetric) Visualization Examples ===\n\n";

    // Example 1: Gaussian Blob
    std::cout << "1. Gaussian Blob: f(x,y,z) = exp(-(x^2 + y^2 + z^2) / sigma^2)\n";
    std::cout << "   A smooth peak centered at the origin\n";
    ScalarFunction<3> gaussianBlob{[](const VectorN<Real, 3>& v) {
        Real sigma = 0.5;
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        return std::exp(-r2 / (sigma * sigma));
    }};
    Visualizer::VisualizeScalarFunc3DCartesian(gaussianBlob, "Gaussian Blob",
                                               -1.5, 1.5, 30,
                                               -1.5, 1.5, 30,
                                               -1.5, 1.5, 30,
                                               "viz_scalar3d_gaussian.mml");

    // Example 2: Sinusoidal Wave
    std::cout << "2. Triply Periodic Wave: f = sin(pi*x) * sin(pi*y) * sin(pi*z)\n";
    std::cout << "   Creates a 3D checkerboard-like pattern\n";
    ScalarFunctionFromStdFunc<3> sinusoidalWave([](const VectorN<Real, 3>& v) -> Real {
        return static_cast<Real>(std::sin(Constants::PI * v[0]) *
               std::sin(Constants::PI * v[1]) *
               std::sin(Constants::PI * v[2]));
    });
    Visualizer::VisualizeScalarFunc3DCartesian(sinusoidalWave, "Triply Periodic Wave",
                                               -2.0, 2.0, 30,
                                               -2.0, 2.0, 30,
                                               -2.0, 2.0, 30,
                                               "viz_scalar3d_wave.mml");

    // Example 3: Gyroid - a triply periodic minimal surface (famous in materials science!)
    std::cout << "3. Gyroid: sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x)\n";
    std::cout << "   Famous minimal surface, used in 3D printing and material science\n";
    ScalarFunction<3> gyroid{[](const VectorN<Real, 3>& v) {
        Real scale = Constants::PI / 2.0;  // One full period
        Real x = v[0] * scale;
        Real y = v[1] * scale;
        Real z = v[2] * scale;
        return std::sin(x) * std::cos(y) + 
               std::sin(y) * std::cos(z) + 
               std::sin(z) * std::cos(x);
    }};
    Visualizer::VisualizeScalarFunc3DCartesian(gyroid, "Gyroid (Minimal Surface)",
                                               -2.0, 2.0, 35,
                                               -2.0, 2.0, 35,
                                               -2.0, 2.0, 35,
                                               "viz_scalar3d_gyroid.mml");

    // Example 4: Sphere Distance Function (SDF)
    std::cout << "4. Signed Distance to Sphere: sqrt(x^2+y^2+z^2) - r\n";
    std::cout << "   Negative inside, zero on surface, positive outside\n";
    ScalarFunction<3> sphereDistance{[](const VectorN<Real, 3>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return r - 1.0;  // Distance to unit sphere
    }};
    Visualizer::VisualizeScalarFunc3DCartesian(sphereDistance, "Sphere Distance Function",
                                               -1.5, 1.5, 30,
                                               -1.5, 1.5, 30,
                                               -1.5, 1.5, 30,
                                               "viz_scalar3d_sphere_sdf.mml");

    // Example 5: Torus Distance Function
    std::cout << "5. Signed Distance to Torus\n";
    std::cout << "   Donut-shaped distance field\n";
    ScalarFunction<3> torusDistance{[](const VectorN<Real, 3>& v) {
        Real R = 0.7;  // Major radius
        Real r = 0.3;  // Minor radius
        Real q = std::sqrt(v[0]*v[0] + v[1]*v[1]) - R;
        return std::sqrt(q*q + v[2]*v[2]) - r;
    }};
    Visualizer::VisualizeScalarFunc3DCartesian(torusDistance, "Torus Distance Function",
                                               -1.5, 1.5, 35,
                                               -1.5, 1.5, 35,
                                               -1.0, 1.0, 25,
                                               "viz_scalar3d_torus_sdf.mml");

    std::cout << "\n3D Scalar function examples complete!\n";
}
