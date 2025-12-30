/**
 * @file show_scalar_function.cpp
 * @brief Demonstrates visualization of 2D scalar functions (surfaces)
 * 
 * Uses the MML Visualizer infrastructure to display 2D scalar functions
 * as 3D surfaces.
 * Cross-platform: works on Windows, Linux, and macOS.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/VectorN.h"
#include "core/Fields.h"
#include "tools/Visualizer.h"
#endif

using namespace MML;

void Show_Scalar_Function_Examples()
{
    std::cout << "\n=== Scalar Function (Surface) Visualization Examples ===\n\n";

    // Example 1: Simple paraboloid
    std::cout << "1. Paraboloid: f(x,y) = x^2 + y^2\n";
    ScalarFunction<2> paraboloid{[](const VectorN<Real, 2>& v) {
        return v[0]*v[0] + v[1]*v[1];
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(paraboloid, "Paraboloid z = x^2 + y^2",
                                               -5.0, 5.0, 30, -5.0, 5.0, 30,
                                               "viz_scalar_paraboloid.txt");

    // Example 2: Saddle surface (hyperbolic paraboloid)
    std::cout << "2. Saddle: f(x,y) = x^2 - y^2\n";
    ScalarFunction<2> saddle{[](const VectorN<Real, 2>& v) {
        return v[0]*v[0] - v[1]*v[1];
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(saddle, "Saddle z = x^2 - y^2",
                                               -5.0, 5.0, 30, -5.0, 5.0, 30,
                                               "viz_scalar_saddle.txt");

    // Example 3: Monkey saddle
    std::cout << "3. Monkey Saddle: f(x,y) = x^3 - 3xy^2\n";
    ScalarFunction<2> monkeySaddle{[](const VectorN<Real, 2>& v) {
        return v[0]*v[0]*v[0] - 3*v[0]*v[1]*v[1];
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(monkeySaddle, "Monkey Saddle z = x^3 - 3xy^2",
                                               -3.0, 3.0, 40, -3.0, 3.0, 40,
                                               "viz_scalar_monkey_saddle.txt");

    // Example 4: Ripple surface
    std::cout << "4. Ripple: f(x,y) = sin(sqrt(x^2 + y^2))\n";
    ScalarFunction<2> ripple{[](const VectorN<Real, 2>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        if (r < 1e-10) return Real(1.0);  // sinc(0) = 1
        return 50.0 * std::sin(r) / r;
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(ripple, "Sinc Ripple",
                                               -15.0, 15.0, 50, -15.0, 15.0, 50,
                                               "viz_scalar_ripple.txt");

    // Example 5: Gaussian bump
    std::cout << "5. Gaussian: f(x,y) = exp(-(x^2 + y^2)/10)\n";
    ScalarFunction<2> gaussian{[](const VectorN<Real, 2>& v) {
        return 50.0 * std::exp(-(v[0]*v[0] + v[1]*v[1]) / 10.0);
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(gaussian, "Gaussian Bump",
                                               -8.0, 8.0, 40, -8.0, 8.0, 40,
                                               "viz_scalar_gaussian.txt");

    // Example 6: Two-bump surface
    std::cout << "6. Two Bumps: sum of two Gaussians\n";
    ScalarFunction<2> twoBumps{[](const VectorN<Real, 2>& v) {
        Real bump1 = std::exp(-((v[0]-2)*(v[0]-2) + v[1]*v[1]) / 2.0);
        Real bump2 = std::exp(-((v[0]+2)*(v[0]+2) + v[1]*v[1]) / 2.0) * 0.8;
        return 50.0 * (bump1 + bump2);
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(twoBumps, "Two Gaussian Bumps",
                                               -6.0, 6.0, 50, -4.0, 4.0, 40,
                                               "viz_scalar_two_bumps.txt");

    // Example 7: Egg crate surface
    std::cout << "7. Egg Crate: f(x,y) = sin(x)*sin(y)\n";
    ScalarFunction<2> eggCrate{[](const VectorN<Real, 2>& v) {
        return 20.0 * std::sin(v[0]) * std::sin(v[1]);
    }};
    Visualizer::VisualizeScalarFunc2DCartesian(eggCrate, "Egg Crate sin(x)*sin(y)",
                                               -2*Constants::PI, 2*Constants::PI, 40,
                                               -2*Constants::PI, 2*Constants::PI, 40,
                                               "viz_scalar_egg_crate.txt");

    std::cout << "\nScalar function visualization examples complete!\n";
}
