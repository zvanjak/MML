/**
 * @file show_vector_field_3d.cpp
 * @brief Demonstrates visualization of 3D vector fields
 * 
 * Uses the MML Visualizer infrastructure to display 3D vector fields.
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

void Show_Vector_Field_3D_Examples()
{
    std::cout << "\n=== 3D Vector Field Visualization Examples ===\n\n";

    // Example 1: Constant flow field
    std::cout << "1. Constant flow: F(x,y,z) = (1, 0, 0)\n";
    VectorFunction<3> constantFlow{[](const VectorN<Real, 3>& v) {
        return VectorN<Real, 3>{1.0, 0.0, 0.0};
    }};
    Visualizer::VisualizeVectorField3DCartesian(constantFlow, "Constant Flow",
                                                -5.0, 5.0, 8, -5.0, 5.0, 8, -5.0, 5.0, 8,
                                                "viz_vector3d_constant.txt");

    // Example 2: Radial field (3D source)
    std::cout << "2. Radial field (3D source): F = r/|r|\n";
    VectorFunction<3> radialField{[](const VectorN<Real, 3>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (r < 0.1) r = 0.1;
        return VectorN<Real, 3>{v[0]/r, v[1]/r, v[2]/r};
    }};
    Visualizer::VisualizeVectorField3DCartesian(radialField, "3D Radial Source",
                                                -5.0, 5.0, 8, -5.0, 5.0, 8, -5.0, 5.0, 8,
                                                "viz_vector3d_radial.txt");

    // Example 3: Gravitational field of two masses
    std::cout << "3. Two-body gravity: masses at (3,0,0) and (-3,0,0)\n";
    VectorFunction<3> gravityField{[](const VectorN<Real, 3>& v) {
        const VectorN<Real, 3> m1_pos{3.0, 0.0, 0.0};
        const VectorN<Real, 3> m2_pos{-3.0, 0.0, 0.0};
        const Real G = 10.0;
        const Real M1 = 100.0;
        const Real M2 = 100.0;
        
        VectorN<Real, 3> r1 = v - m1_pos;
        VectorN<Real, 3> r2 = v - m2_pos;
        
        Real d1 = r1.NormL2();
        Real d2 = r2.NormL2();
        
        if (d1 < 0.5) d1 = 0.5;
        if (d2 < 0.5) d2 = 0.5;
        
        // F = -GM * r / |r|^3
        VectorN<Real, 3> F1 = r1 * (-G * M1 / (d1 * d1 * d1));
        VectorN<Real, 3> F2 = r2 * (-G * M2 / (d2 * d2 * d2));
        
        return F1 + F2;
    }};
    Visualizer::VisualizeVectorField3DCartesian(gravityField, "Two-Body Gravity",
                                                -8.0, 8.0, 10, -5.0, 5.0, 7, -5.0, 5.0, 7,
                                                "viz_vector3d_gravity.txt");

    // Example 4: Magnetic field of a current loop (simplified)
    std::cout << "4. Magnetic field pattern (simplified dipole)\n";
    VectorFunction<3> magneticField{[](const VectorN<Real, 3>& v) {
        // Simplified magnetic dipole field aligned with z-axis
        Real x = v[0], y = v[1], z = v[2];
        Real r_sq = x*x + y*y + z*z;
        if (r_sq < 0.25) r_sq = 0.25;
        Real r = std::sqrt(r_sq);
        Real r5 = r_sq * r_sq * r;
        
        // Dipole field: B = (3(m·r)r - m*r^2) / r^5, m along z
        Real m = 100.0;  // magnetic moment magnitude
        Real factor = m / r5;
        
        Real Bx = 3.0 * z * x * factor;
        Real By = 3.0 * z * y * factor;
        Real Bz = (3.0 * z * z - r_sq) * factor;
        
        return VectorN<Real, 3>{Bx, By, Bz};
    }};
    Visualizer::VisualizeVectorField3DCartesian(magneticField, "Magnetic Dipole",
                                                -4.0, 4.0, 8, -4.0, 4.0, 8, -4.0, 4.0, 8,
                                                "viz_vector3d_magnetic.txt");

    // Example 5: Swirling flow (helical)
    std::cout << "5. Helical flow: rotation + upward flow\n";
    VectorFunction<3> helicalFlow{[](const VectorN<Real, 3>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        if (r < 0.1) r = 0.1;
        
        // Tangential component (rotation around z-axis)
        Real Vx = -v[1] / r;
        Real Vy = v[0] / r;
        // Vertical component
        Real Vz = 0.5;
        
        return VectorN<Real, 3>{Vx, Vy, Vz};
    }};
    Visualizer::VisualizeVectorField3DCartesian(helicalFlow, "Helical Flow",
                                                -5.0, 5.0, 8, -5.0, 5.0, 8, -5.0, 5.0, 8,
                                                "viz_vector3d_helical.txt");

    // Example 6: Saddle flow in 3D
    std::cout << "6. 3D Saddle: F = (x, -y, 0)\n";
    VectorFunction<3> saddleFlow{[](const VectorN<Real, 3>& v) {
        return VectorN<Real, 3>{v[0], -v[1], 0.0};
    }};
    Visualizer::VisualizeVectorField3DCartesian(saddleFlow, "3D Saddle Flow",
                                                -5.0, 5.0, 8, -5.0, 5.0, 8, -3.0, 3.0, 5,
                                                "viz_vector3d_saddle.txt");

    std::cout << "\n3D vector field visualization examples complete!\n";
}
