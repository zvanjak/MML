///////////////////////////////////////////////////////////////////////////////////////////
/// @file vector_fields_visualization.cpp
/// @brief Chapter 02: Vector Field Visualization Demos
/// @details Demonstrates visualization of vector fields:
///          - 2D vector fields F: R² → R² displayed as arrow plots
///          - 3D vector fields F: R³ → R³ displayed as 3D arrow plots
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
/// @brief Demo: 2D Vector Fields
/// @details Vector fields F: R² → R² visualized as arrow plots.
///          Shows classic fields from physics: vortex, dipole, saddle point.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_VectorField2D()
{
    std::cout << "\n=== Demo: 2D Vector Fields ===\n";
    
    //-------------------------------------------------------------------------
    // Vortex Field - circular flow around the origin
    // Classic rotational field: F(x,y) = (-y, x) / r
    // Models whirlpools, hurricanes, tornadoes
    //-------------------------------------------------------------------------
    VectorFunction<2> vortex{ [](const VectorN<Real, 2>& v) {
        Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        if (r < 0.1) r = 0.1;  // Avoid singularity at origin
        return VectorN<Real, 2>{ -v[1]/r, v[0]/r };
    }};
    
    Visualizer::VisualizeVectorField2DCartesian(
        vortex, "Vortex Field: F = (-y, x)/r",
        -5.0, 5.0, 20, -5.0, 5.0, 20,
        "ch02_vortex_2d.mml");
    
    std::cout << "   Visualized: Vortex - circular flow, curl ≠ 0\n";
    
    //-------------------------------------------------------------------------
    // Electric Dipole Field - two opposite charges
    // Positive charge at (-2, 0), negative at (2, 0)
    // Fundamental in electrostatics, models molecules
    //-------------------------------------------------------------------------
    VectorFunction<2> dipole{ [](const VectorN<Real, 2>& v) {
        // Positive charge at (-2, 0) - field points away
        Real dx1 = v[0] + 2.0, dy1 = v[1];
        Real r1_sq = dx1*dx1 + dy1*dy1;
        if (r1_sq < 0.25) r1_sq = 0.25;
        
        // Negative charge at (2, 0) - field points toward
        Real dx2 = v[0] - 2.0, dy2 = v[1];
        Real r2_sq = dx2*dx2 + dy2*dy2;
        if (r2_sq < 0.25) r2_sq = 0.25;
        
        // Electric field: E = k*q*r/|r|³ (inverse square law)
        Real Ex = dx1/(r1_sq) - dx2/(r2_sq);
        Real Ey = dy1/(r1_sq) - dy2/(r2_sq);
        
        return VectorN<Real, 2>{ Ex, Ey };
    }};
    
    Visualizer::VisualizeVectorField2DCartesian(
        dipole, "Electric Dipole Field",
        -6.0, 6.0, 25, -4.0, 4.0, 17,
        "ch02_dipole_2d.mml");
    
    std::cout << "   Visualized: Electric dipole - field lines from + to -\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 3D Vector Fields
/// @details Vector fields F: R³ → R³ visualized as 3D arrow plots.
///          Shows gravitational and magnetic fields.
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_Demo_VectorField3D()
{
    std::cout << "\n=== Demo: 3D Vector Fields ===\n";
    
    //-------------------------------------------------------------------------
    // Two-Body Gravitational Field
    // Two massive objects at (±3, 0, 0), Newton's law of gravitation
    // F = -GMm*r/|r|³ (attractive, inverse square law)
    //-------------------------------------------------------------------------
    VectorFunction<3> gravity{ [](const VectorN<Real, 3>& v) {
        const VectorN<Real, 3> m1_pos{ 3.0, 0.0, 0.0 };
        const VectorN<Real, 3> m2_pos{ -3.0, 0.0, 0.0 };
        const Real G = 10.0;   // Gravitational constant (scaled)
        const Real M = 100.0;  // Mass of each body
        
        VectorN<Real, 3> r1 = v - m1_pos;
        VectorN<Real, 3> r2 = v - m2_pos;
        
        Real d1 = r1.NormL2();
        Real d2 = r2.NormL2();
        
        if (d1 < 0.5) d1 = 0.5;  // Avoid singularity near masses
        if (d2 < 0.5) d2 = 0.5;
        
        // F = -GM * r / |r|³ for each mass
        VectorN<Real, 3> F1 = r1 * (-G * M / (d1 * d1 * d1));
        VectorN<Real, 3> F2 = r2 * (-G * M / (d2 * d2 * d2));
        
        return F1 + F2;
    }};
    
    Visualizer::VisualizeVectorField3DCartesian(
        gravity, "Two-Body Gravitational Field",
        -8.0, 8.0, 12, -5.0, 5.0, 8, -5.0, 5.0, 8,
        "ch02_gravity_3d.mml");
    
    std::cout << "   Visualized: Two-body gravity - attractive inverse-square\n";
    
    //-------------------------------------------------------------------------
    // Magnetic Dipole Field
    // Simplified model of a bar magnet aligned with z-axis
    // B = (μ₀/4π) * (3(m·r)r/r⁵ - m/r³)
    //-------------------------------------------------------------------------
    VectorFunction<3> magnetic_dipole{ [](const VectorN<Real, 3>& v) {
        Real x = v[0], y = v[1], z = v[2];
        Real r_sq = x*x + y*y + z*z;
        if (r_sq < 0.25) r_sq = 0.25;
        Real r = std::sqrt(r_sq);
        Real r5 = r_sq * r_sq * r;
        Real r3 = r_sq * r;
        
        // Magnetic moment along z-axis: m = (0, 0, 1)
        // B = 3(m·r)r/r⁵ - m/r³
        Real m_dot_r = z;  // (0,0,1) · (x,y,z) = z
        
        Real Bx = 3.0 * m_dot_r * x / r5;
        Real By = 3.0 * m_dot_r * y / r5;
        Real Bz = 3.0 * m_dot_r * z / r5 - 1.0 / r3;
        
        return VectorN<Real, 3>{ Bx, By, Bz };
    }};
    
    Visualizer::VisualizeVectorField3DCartesian(
        magnetic_dipole, "Magnetic Dipole Field (Bar Magnet)",
        -3.0, 3.0, 10, -3.0, 3.0, 10, -3.0, 3.0, 10,
        "ch02_magnetic_dipole_3d.mml");
    
    std::cout << "   Visualized: Magnetic dipole - classic bar magnet field\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief Entry point for Vector Field visualization demos
///////////////////////////////////////////////////////////////////////////////////////////
void Chapter02_VectorFieldVisualization()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║          CHAPTER 02: Vector Field Visualization                   ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
    
    Chapter02_Demo_VectorField2D();
    Chapter02_Demo_VectorField3D();
}
