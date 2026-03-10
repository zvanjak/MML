///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_fields.cpp                                                ///
///  Description: Brief demonstration of Fields.h - scalar and vector fields          ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/Fields.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Docs_Demo_Fields()
{
    std::cout << "=== Fields Demo ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    // Inverse radial potential (like gravity: Φ = -GM/r)
    VectorN<Real, 3> point{1.0, 1.0, 1.0};
    Real r = point.NormL2();
    
    std::cout << "Point: (1, 1, 1), |r| = " << r << std::endl;
    
    // Scalar potential field
    Real potential = Fields::InverseRadialPotentialFieldCart(point);
    std::cout << "Inverse radial potential Φ = 1/|r| = " << potential << std::endl;
    
    // With gravitational constant (example: G*M = -10)
    Fields::InverseRadialFieldCart gravField(-10.0);
    std::cout << "Gravitational potential (G*M=-10): " << gravField(point) << std::endl;
    
    // Force field (F = -∇Φ)
    VectorN<Real, 3> force = Fields::InverseRadialPotentialForceFieldCart(point);
    std::cout << "Force field at point: (" << force[0] << ", " << force[1] << ", " << force[2] << ")" << std::endl;
    std::cout << "Force magnitude: " << force.NormL2() << std::endl;
    
    // Spherical coordinates version
    VectorN<Real, 3> spherPoint{2.0, 0.5, 0.3}; // r=2, theta, phi
    Real potentialSpher = Fields::InverseRadialPotentialFieldSpher(spherPoint);
    std::cout << "\nSpherical at r=2: potential = " << potentialSpher << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}
