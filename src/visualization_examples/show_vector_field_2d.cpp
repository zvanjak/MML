/**
 * @file show_vector_field_2d.cpp
 * @brief Demonstrates visualization of 2D vector fields
 * 
 * Uses the MML Visualizer infrastructure to display 2D vector fields.
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

void Show_Vector_Field_2D_Examples()
{
    std::cout << "\n=== 2D Vector Field Visualization Examples ===\n\n";

    // Example 1: Constant flow field
    std::cout << "1. Constant flow: F(x,y) = (1, 0.5)\n";
    VectorFunction<2> constantFlow{[](const VectorN<Real, 2>& v) {
        return VectorN<Real, 2>{1.0, 0.5};
    }};
    Visualizer::VisualizeVectorField2DCartesian(constantFlow, "Constant Flow",
                                                -5.0, 5.0, 15, -5.0, 5.0, 15,
                                                "viz_vector2d_constant.txt");

    // Example 2: Radial field (source/sink)
    std::cout << "2. Radial field (source): F(x,y) = (x, y) / r\n";
    VectorFunction<2> radialField{[](const VectorN<Real, 2>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        if (r < 0.1) r = 0.1;  // Avoid division by zero
        return VectorN<Real, 2>{v[0]/r, v[1]/r};
    }};
    Visualizer::VisualizeVectorField2DCartesian(radialField, "Radial Source Field",
                                                -5.0, 5.0, 15, -5.0, 5.0, 15,
                                                "viz_vector2d_radial.txt");

    // Example 3: Rotational field (vortex)
    std::cout << "3. Rotational field (vortex): F(x,y) = (-y, x) / r\n";
    VectorFunction<2> vortexField{[](const VectorN<Real, 2>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        if (r < 0.1) r = 0.1;
        return VectorN<Real, 2>{-v[1]/r, v[0]/r};
    }};
    Visualizer::VisualizeVectorField2DCartesian(vortexField, "Vortex Field",
                                                -5.0, 5.0, 15, -5.0, 5.0, 15,
                                                "viz_vector2d_vortex.txt");

    // Example 4: Saddle point field
    std::cout << "4. Saddle point: F(x,y) = (x, -y)\n";
    VectorFunction<2> saddleField{[](const VectorN<Real, 2>& v) {
        return VectorN<Real, 2>{v[0], -v[1]};
    }};
    Visualizer::VisualizeVectorField2DCartesian(saddleField, "Saddle Point Field",
                                                -5.0, 5.0, 15, -5.0, 5.0, 15,
                                                "viz_vector2d_saddle.txt");

    // Example 5: Dipole field (two opposite charges)
    std::cout << "5. Dipole field: two opposite charges at (-2,0) and (2,0)\n";
    VectorFunction<2> dipoleField{[](const VectorN<Real, 2>& v) {
        // Positive charge at (-2, 0)
        Real dx1 = v[0] + 2.0;
        Real dy1 = v[1];
        Real r1_sq = dx1*dx1 + dy1*dy1;
        if (r1_sq < 0.25) r1_sq = 0.25;
        
        // Negative charge at (2, 0)
        Real dx2 = v[0] - 2.0;
        Real dy2 = v[1];
        Real r2_sq = dx2*dx2 + dy2*dy2;
        if (r2_sq < 0.25) r2_sq = 0.25;
        
        // Field = sum of contributions (positive attracts, negative repels in our convention)
        Real Ex = dx1/r1_sq - dx2/r2_sq;
        Real Ey = dy1/r1_sq - dy2/r2_sq;
        
        return VectorN<Real, 2>{Ex, Ey};
    }};
    Visualizer::VisualizeVectorField2DCartesian(dipoleField, "Electric Dipole Field",
                                                -6.0, 6.0, 20, -4.0, 4.0, 15,
                                                "viz_vector2d_dipole.txt");

    // Example 6: Shear flow
    std::cout << "6. Shear flow: F(x,y) = (y, 0)\n";
    VectorFunction<2> shearFlow{[](const VectorN<Real, 2>& v) {
        return VectorN<Real, 2>{v[1], 0.0};
    }};
    Visualizer::VisualizeVectorField2DCartesian(shearFlow, "Shear Flow",
                                                -5.0, 5.0, 15, -5.0, 5.0, 15,
                                                "viz_vector2d_shear.txt");

    // Example 7: Double vortex
    std::cout << "7. Double vortex: two counter-rotating vortices\n";
    VectorFunction<2> doubleVortex{[](const VectorN<Real, 2>& v) {
        // Vortex 1 at (-2, 0), counter-clockwise
        Real dx1 = v[0] + 2.0;
        Real dy1 = v[1];
        Real r1 = std::sqrt(dx1*dx1 + dy1*dy1);
        if (r1 < 0.3) r1 = 0.3;
        
        // Vortex 2 at (2, 0), clockwise
        Real dx2 = v[0] - 2.0;
        Real dy2 = v[1];
        Real r2 = std::sqrt(dx2*dx2 + dy2*dy2);
        if (r2 < 0.3) r2 = 0.3;
        
        Real Vx = -dy1/r1 + dy2/r2;
        Real Vy = dx1/r1 - dx2/r2;
        
        return VectorN<Real, 2>{Vx, Vy};
    }};
    Visualizer::VisualizeVectorField2DCartesian(doubleVortex, "Double Vortex",
                                                -6.0, 6.0, 20, -4.0, 4.0, 15,
                                                "viz_vector2d_double_vortex.txt");

    std::cout << "\n2D vector field visualization examples complete!\n";
}
