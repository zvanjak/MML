/**
 * @file show_parametric_curve_3d.cpp
 * @brief Demonstrates visualization of 3D parametric curves
 * 
 * Uses the MML Visualizer infrastructure to display parametric curves in 3D.
 * Cross-platform: works on Windows, Linux, and macOS.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/VectorN.h"
#include "core/Curves.h"
#include "tools/Visualizer.h"
#endif

using namespace MML;

void Show_Parametric_Curve_3D_Examples()
{
    std::cout << "\n=== 3D Parametric Curve Visualization Examples ===\n\n";

    // Example 1: Helix
    std::cout << "1. Helix: (cos(t), sin(t), t/5)\n";
    ParametricCurve<3> helix{[](Real t) {
        return VectorN<Real, 3>{std::cos(t), std::sin(t), t / 5.0} * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(helix, "Helix",
                                      0.0, 10*Constants::PI, 500,
                                      "viz_curve3d_helix.txt");

    // Example 2: Using predefined Helix curve from library
    std::cout << "2. Library Helix (radius=5, pitch=1)\n";
    Curves::HelixCurve libHelix(5.0, 1.0);
    Visualizer::VisualizeParamCurve3D(libHelix, "Library Helix",
                                      -20.0, 20.0, 500,
                                      "viz_curve3d_lib_helix.txt");

    // Example 3: Toroidal spiral
    std::cout << "3. Toroidal Spiral\n";
    Curves::ToroidalSpiralCurve toroid(50.0);
    Visualizer::VisualizeParamCurve3D(toroid, "Toroidal Spiral",
                                      0.0, 2*Constants::PI, 2000,
                                      "viz_curve3d_toroidal.txt");

    // Example 4: Trefoil knot
    std::cout << "4. Trefoil Knot\n";
    ParametricCurve<3> trefoil{[](Real t) {
        return VectorN<Real, 3>{
            std::sin(t) + 2*std::sin(2*t),
            std::cos(t) - 2*std::cos(2*t),
            -std::sin(3*t)
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(trefoil, "Trefoil Knot",
                                      0.0, 2*Constants::PI, 300,
                                      "viz_curve3d_trefoil.txt");

    // Example 5: Viviani's curve (intersection of sphere and cylinder)
    std::cout << "5. Viviani's Curve\n";
    ParametricCurve<3> viviani{[](Real t) {
        Real a = 2.0;
        return VectorN<Real, 3>{
            a * (1 + std::cos(t)),
            a * std::sin(t),
            2 * a * std::sin(t / 2)
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(viviani, "Viviani's Curve",
                                      0.0, 4*Constants::PI, 300,
                                      "viz_curve3d_viviani.txt");

    // Example 6: Conical helix (helix on a cone)
    std::cout << "6. Conical Helix\n";
    ParametricCurve<3> conicalHelix{[](Real t) {
        Real r = t / 10.0;  // radius increases with t
        return VectorN<Real, 3>{
            r * std::cos(t),
            r * std::sin(t),
            t / 5.0
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(conicalHelix, "Conical Helix",
                                      0.1, 10*Constants::PI, 500,
                                      "viz_curve3d_conical_helix.txt");

    // Example 7: Figure-8 in 3D space
    std::cout << "7. 3D Figure-8\n";
    ParametricCurve<3> figure8_3d{[](Real t) {
        return VectorN<Real, 3>{
            std::sin(t),
            std::sin(t) * std::cos(t),
            std::sin(2*t) / 2
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(figure8_3d, "3D Figure-8",
                                      0.0, 2*Constants::PI, 200,
                                      "viz_curve3d_figure8.txt");

    // Example 8: Spherical spiral
    std::cout << "8. Spherical Spiral (loxodrome)\n";
    ParametricCurve<3> sphericalSpiral{[](Real t) {
        Real R = 3.0;  // sphere radius
        Real theta = t;  // polar angle
        Real phi = 4 * t;  // azimuthal angle
        return VectorN<Real, 3>{
            R * std::sin(theta) * std::cos(phi),
            R * std::sin(theta) * std::sin(phi),
            R * std::cos(theta)
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(sphericalSpiral, "Spherical Spiral",
                                      0.0, Constants::PI, 500,
                                      "viz_curve3d_spherical_spiral.txt");

    // Example 9: Cinquefoil knot (5_1 knot)
    std::cout << "9. Cinquefoil Knot\n";
    ParametricCurve<3> cinquefoil{[](Real t) {
        return VectorN<Real, 3>{
            std::cos(2*t) * (2 + std::cos(5*t)),
            std::sin(2*t) * (2 + std::cos(5*t)),
            std::sin(5*t)
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(cinquefoil, "Cinquefoil Knot",
                                      0.0, 2*Constants::PI, 500,
                                      "viz_curve3d_cinquefoil.txt");

    // Example 10: Slinky curve
    std::cout << "10. Slinky Curve\n";
    ParametricCurve<3> slinky{[](Real t) {
        Real R = 5.0;  // major radius
        Real r = 1.0;  // minor radius
        Real n = 20;   // number of coils
        return VectorN<Real, 3>{
            (R + r * std::cos(n * t)) * std::cos(t),
            (R + r * std::cos(n * t)) * std::sin(t),
            r * std::sin(n * t)
        } * 50.0;
    }};
    Visualizer::VisualizeParamCurve3D(slinky, "Slinky Curve",
                                      0.0, 2*Constants::PI, 1000,
                                      "viz_curve3d_slinky.txt");

    std::cout << "\n3D parametric curve visualization examples complete!\n";
}
