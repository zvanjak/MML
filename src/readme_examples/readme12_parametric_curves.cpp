///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme12_parametric_curves.cpp                                      ///
///  Description: README example - Parametric Curves                                  ///
///               Demonstrates curve properties, Frenet frame, curvature              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/VectorN.h"
#include "base/Function.h"
#include "core/Curves.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Readme_ParametricCurves()
{
    std::cout << std::endl;
    std::cout << "=== Parametric Curves ===" << std::endl;

    // Define a 3D helix: r(t) = (cos(t), sin(t), 0.2t)
    Curves::CurveCartesian3D helix([](Real t) -> VectorN<Real, 3> { 
        return {cos(t), sin(t), 0.2*t}; 
    });

    Real t = Constants::PI / 4;
    std::cout << "Helix: r(t) = (cos(t), sin(t), 0.2t)" << std::endl;
    std::cout << "Evaluated at t = π/4:" << std::endl;

    // Curve properties at parameter t
    auto pos = helix(t);                    // Position on curve
    auto tangent = helix.getTangent(t);     // Tangent vector dr/dt
    auto unit_tan = helix.getTangentUnit(t);// Unit tangent T
    auto normal = helix.getNormal(t);       // Normal vector (acceleration)
    auto binormal = helix.getBinormal(t);   // Binormal B = T × N

    std::cout << std::setprecision(6);
    std::cout << std::endl << "  Position:       " << pos << std::endl;
    std::cout << "  Tangent dr/dt:  " << tangent << std::endl;
    std::cout << "  Unit tangent:   " << unit_tan << std::endl;
    std::cout << "  Normal (accel): " << normal << std::endl;
    std::cout << "  Binormal:       " << binormal << std::endl;

    // Curvature κ (Frenet-Serret apparatus)
    Real curvature = helix.getCurvature(t);
    std::cout << std::endl << "  Curvature κ:    " << curvature << std::endl;

    // Verify orthogonality of Frenet frame
    auto normal_unit = helix.getNormalUnit(t);
    Real t_dot_n = ScalarProduct(unit_tan, normal_unit);
    std::cout << "  T·N (should be 0): " << t_dot_n << std::endl;

    // Predefined curves
    std::cout << std::endl << "Predefined curves:" << std::endl;
    
    Curves::LemniscateCurve lemniscate;
    auto lem_pos = lemniscate(0.5);
    std::cout << "  Lemniscate at t=0.5: " << lem_pos << std::endl;

    Curves::ToroidalSpiralCurve torus(5, 2);  // (n, R) - n wraps, major radius R
    auto torus_pos = torus(1.0);
    std::cout << "  Toroidal spiral at t=1: " << torus_pos << std::endl;

    Curves::Circle3DXZCurve circle(3.0);  // Circle in XZ plane, radius 3
    auto circle_pos = circle(Constants::PI/2);
    std::cout << "  Circle (r=3) in XZ at t=π/2: " << circle_pos << std::endl;

    // 2D curves
    Curves::ArchimedeanSpiralCurve spiral(0.5);  // r = 0.5*t
    auto spiral_pos = spiral(4*Constants::PI);
    std::cout << "  Archimedean spiral at t=4π: " << spiral_pos << std::endl;

    std::cout << std::endl;
}
