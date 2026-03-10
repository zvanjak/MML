///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_vector_types.cpp                                          ///
///  Description: Brief demonstration of VectorTypes.h - Vec2, Vec3, Vec4 types       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Docs_Demo_VectorTypes()
{
    std::cout << "=== VectorTypes Demo ===" << std::endl;
    
    // 2D vectors
    Vector2Cartesian v2a(3.0, 4.0);
    Vector2Cartesian v2b(1.0, 2.0);
    std::cout << "v2a = (" << v2a.X() << ", " << v2a.Y() << ")" << std::endl;
    std::cout << "v2a.NormL2() = " << v2a.NormL2() << " (should be 5)" << std::endl;
    
    Vector2Cartesian v2sum = v2a + v2b;
    std::cout << "v2a + v2b = (" << v2sum.X() << ", " << v2sum.Y() << ")" << std::endl;
    
    // 3D vectors
    Vector3Cartesian v3a(1.0, 2.0, 3.0);
    Vector3Cartesian v3b(4.0, 5.0, 6.0);
    std::cout << "\nv3a = (" << v3a.X() << ", " << v3a.Y() << ", " << v3a.Z() << ")" << std::endl;
    
    // Dot product (operator*)
    Real dot = v3a * v3b;
    std::cout << "v3a · v3b = " << dot << std::endl;
    
    // Cross product (free function VectorProduct)
    Vector3Cartesian cross = VectorProduct(v3a, v3b);
    std::cout << "v3a × v3b = (" << cross.X() << ", " << cross.Y() << ", " << cross.Z() << ")" << std::endl;
    
    // Unit vector
    Vector3Cartesian unit = v3a.GetAsUnitVector();
    std::cout << "v3a normalized: norm = " << unit.NormL2() << std::endl;
    
    // 4D vectors (Vector4Minkowski - 4D Minkowski spacetime vector)
    Vector4Minkowski v4{1.0, 2.0, 3.0, 1.0};
    std::cout << "\nv4 = (" << v4[0] << ", " << v4[1] << ", " << v4[2] << ", " << v4[3] << ")" << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}
