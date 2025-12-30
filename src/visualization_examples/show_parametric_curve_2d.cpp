/**
 * @file show_parametric_curve_2d.cpp
 * @brief Demonstrates visualization of 2D parametric curves
 * 
 * Uses the MML Visualizer infrastructure to display parametric curves in 2D.
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

void Show_Parametric_Curve_2D_Examples()
{
    std::cout << "\n=== 2D Parametric Curve Visualization Examples ===\n\n";

    // Example 1: Circle
    std::cout << "1. Circle: (cos(t), sin(t))\n";
    ParametricCurve<2> circle{[](Real t) {
        return VectorN<Real, 2>{std::cos(t), std::sin(t)};
    }};
    Visualizer::VisualizeParamCurve2D(circle, "Circle",
                                      0.0, 2*Constants::PI, 100,
                                      "viz_curve2d_circle.txt");

    // Example 2: Ellipse
    std::cout << "2. Ellipse: (3*cos(t), 2*sin(t))\n";
    ParametricCurve<2> ellipse{[](Real t) {
        return VectorN<Real, 2>{3.0*std::cos(t), 2.0*std::sin(t)};
    }};
    Visualizer::VisualizeParamCurve2D(ellipse, "Ellipse (a=3, b=2)",
                                      0.0, 2*Constants::PI, 100,
                                      "viz_curve2d_ellipse.txt");

    // Example 3: Lissajous curve
    std::cout << "3. Lissajous curve: (sin(3t), sin(4t))\n";
    ParametricCurve<2> lissajous{[](Real t) {
        return VectorN<Real, 2>{std::sin(3*t), std::sin(4*t)};
    }};
    Visualizer::VisualizeParamCurve2D(lissajous, "Lissajous (3:4)",
                                      0.0, 2*Constants::PI, 500,
                                      "viz_curve2d_lissajous.txt");

    // Example 4: Cardioid
    std::cout << "4. Cardioid: a*(1-cos(t)) in polar\n";
    ParametricCurve<2> cardioid{[](Real t) {
        Real r = 2.0 * (1.0 - std::cos(t));
        return VectorN<Real, 2>{r * std::cos(t), r * std::sin(t)};
    }};
    Visualizer::VisualizeParamCurve2D(cardioid, "Cardioid",
                                      0.0, 2*Constants::PI, 200,
                                      "viz_curve2d_cardioid.txt");

    // Example 5: Spiral (Archimedean)
    std::cout << "5. Archimedean spiral: r = a*t\n";
    ParametricCurve<2> spiral{[](Real t) {
        Real r = 0.2 * t;
        return VectorN<Real, 2>{r * std::cos(t), r * std::sin(t)};
    }};
    Visualizer::VisualizeParamCurve2D(spiral, "Archimedean Spiral",
                                      0.0, 6*Constants::PI, 500,
                                      "viz_curve2d_spiral.txt");

    // Example 6: Figure-8 (Lemniscate)
    std::cout << "6. Lemniscate (figure-8)\n";
    ParametricCurve<2> lemniscate{[](Real t) {
        Real denom = 1.0 + std::sin(t)*std::sin(t);
        Real a = 3.0;
        return VectorN<Real, 2>{
            a * std::cos(t) / denom,
            a * std::sin(t) * std::cos(t) / denom
        };
    }};
    Visualizer::VisualizeParamCurve2D(lemniscate, "Lemniscate",
                                      0.0, 2*Constants::PI, 200,
                                      "viz_curve2d_lemniscate.txt");

    // Example 7: Rose curve (4 petals)
    std::cout << "7. Rose curve: r = cos(2t)\n";
    ParametricCurve<2> rose{[](Real t) {
        Real r = 3.0 * std::cos(2*t);
        return VectorN<Real, 2>{r * std::cos(t), r * std::sin(t)};
    }};
    Visualizer::VisualizeParamCurve2D(rose, "4-Petal Rose",
                                      0.0, 2*Constants::PI, 200,
                                      "viz_curve2d_rose.txt");

    // Example 8: Epicycloid (3 cusps)
    std::cout << "8. Epicycloid (3 cusps)\n";
    ParametricCurve<2> epicycloid{[](Real t) {
        Real R = 3.0;  // outer circle radius
        Real r = 1.0;  // rolling circle radius
        return VectorN<Real, 2>{
            (R + r) * std::cos(t) - r * std::cos((R + r) * t / r),
            (R + r) * std::sin(t) - r * std::sin((R + r) * t / r)
        };
    }};
    Visualizer::VisualizeParamCurve2D(epicycloid, "Epicycloid (3 cusps)",
                                      0.0, 2*Constants::PI, 300,
                                      "viz_curve2d_epicycloid.txt");

    // Example 9: Butterfly curve
    std::cout << "9. Butterfly curve\n";
    ParametricCurve<2> butterfly{[](Real t) {
        Real factor = std::exp(std::cos(t)) - 2*std::cos(4*t) - std::pow(std::sin(t/12), 5);
        return VectorN<Real, 2>{
            std::sin(t) * factor,
            std::cos(t) * factor
        };
    }};
    Visualizer::VisualizeParamCurve2D(butterfly, "Butterfly Curve",
                                      0.0, 12*Constants::PI, 1000,
                                      "viz_curve2d_butterfly.txt");

    std::cout << "\n2D parametric curve visualization examples complete!\n";
}
