///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_coord_system.cpp                                          ///
///  Description: Brief demonstration of CoordSystem.h - coordinate systems           ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "core/CoordSystem.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Docs_Demo_CoordSystem()
{
    std::cout << "=== CoordSystem Demo ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    
    // Create reference frames
    ReferenceFrame3D worldFrame;
    ReferenceFrame3D localFrame(&worldFrame);
    std::cout << "Created world frame and child local frame" << std::endl;
    
    // Coordinate transformations
    std::cout << "\nSpherical coordinate system:" << std::endl;
    
    // Point in Cartesian: (1, 1, 1)
    VectorN<Real, 3> cartPoint{1.0, 1.0, 1.0};
    std::cout << "Cartesian point: (" << cartPoint[0] << ", " << cartPoint[1] << ", " << cartPoint[2] << ")" << std::endl;
    
    // Convert to spherical (r, theta, phi) using transformation class
    Vector3Spherical spherPoint = CoordTransfCartToSpher.transf(Vector3Cartesian(cartPoint));
    std::cout << "Spherical (r, theta, phi): (" << spherPoint.R() << ", " << spherPoint.Theta() << ", " << spherPoint.Phi() << ")" << std::endl;
    
    // Convert back
    Vector3Cartesian backToCart = CoordTransfSpherToCart.transf(spherPoint);
    std::cout << "Back to Cartesian: (" << backToCart.X() << ", " << backToCart.Y() << ", " << backToCart.Z() << ")" << std::endl;
    
    // Cylindrical system
    std::cout << "\nCylindrical coordinate system:" << std::endl;
    Vector3Cylindrical cylPoint = CoordTransfCartToCyl.transf(Vector3Cartesian(cartPoint));
    std::cout << "Cylindrical (r, phi, z): (" << cylPoint.R() << ", " << cylPoint.Phi() << ", " << cylPoint.Z() << ")" << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}
