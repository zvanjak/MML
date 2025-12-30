#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Geometry.h"
#endif

using namespace MML;

void Docs_Demo_Geometry_2D_Points()
{
    std::cout << "=== 2D Points Demo ===" << std::endl << std::endl;
    
    // Point2Cartesian - construction and access
    std::cout << "Point2Cartesian:" << std::endl;
    Point2Cartesian p1;             // Default: (0, 0)
    Point2Cartesian p2(3.0, 4.0);   // (3, 4)
    
    std::cout << "  p1 (default): (" << p1.X() << ", " << p1.Y() << ")" << std::endl;
    std::cout << "  p2: (" << p2.X() << ", " << p2.Y() << ")" << std::endl;
    
    // Distance
    Real dist = p1.Dist(p2);
    std::cout << "  Distance p1 to p2: " << dist << std::endl;
    
    // Modification
    p2.X() = 5.0;
    p2.Y() = 12.0;
    std::cout << "  p2 after modification: (" << p2.X() << ", " << p2.Y() << ")" << std::endl;
    std::cout << "  New distance: " << p1.Dist(p2) << std::endl << std::endl;
    
    // Point arithmetic
    std::cout << "Point2Cartesian arithmetic:" << std::endl;
    Point2Cartesian a(1, 2), b(3, 4);
    Point2Cartesian sum = a + b;
    Point2Cartesian diff = b - a;
    Point2Cartesian scaled = a * 2.0;
    
    std::cout << "  a + b = (" << sum.X() << ", " << sum.Y() << ")" << std::endl;
    std::cout << "  b - a = (" << diff.X() << ", " << diff.Y() << ")" << std::endl;
    std::cout << "  a * 2 = (" << scaled.X() << ", " << scaled.Y() << ")" << std::endl << std::endl;
    
    // Equality comparison
    std::cout << "Equality comparison:" << std::endl;
    Point2Cartesian c(1.0, 2.0);
    Point2Cartesian d(1.0 + 1e-12, 2.0 + 1e-12);
    std::cout << "  a == c: " << (a == c ? "true" : "false") << std::endl;
    std::cout << "  a.IsEqual(d, 1e-10): " << (a.IsEqual(d, 1e-10) ? "true" : "false") << std::endl << std::endl;
    
    // Point2Polar - construction and conversion
    std::cout << "Point2Polar:" << std::endl;
    Point2Polar polar(5.0, Constants::PI / 4);  // r=5, angle=45°
    std::cout << "  Polar: r=" << polar.R() << ", phi=" << polar.Phi() << " rad" << std::endl;
    
    Point2Cartesian fromPolar = polar.TransfToCart();
    std::cout << "  Converted to Cartesian: (" << fromPolar.X() << ", " << fromPolar.Y() << ")" << std::endl << std::endl;
}

void Docs_Demo_Geometry_3D_Points()
{
    std::cout << "=== 3D Points Demo ===" << std::endl << std::endl;
    
    // Point3Cartesian - construction and access
    std::cout << "Point3Cartesian:" << std::endl;
    Point3Cartesian p1;                 // Default: (0, 0, 0)
    Point3Cartesian p2(1.0, 2.0, 3.0);  // (1, 2, 3)
    
    std::cout << "  p1 (default): (" << p1.X() << ", " << p1.Y() << ", " << p1.Z() << ")" << std::endl;
    std::cout << "  p2: (" << p2.X() << ", " << p2.Y() << ", " << p2.Z() << ")" << std::endl;
    
    // Distance
    Real dist = p1.Dist(p2);
    std::cout << "  Distance p1 to p2: " << dist << std::endl << std::endl;
    
    // 3D Point arithmetic
    std::cout << "Point3Cartesian arithmetic:" << std::endl;
    Point3Cartesian origin(0, 0, 0);
    Point3Cartesian offset(10, 20, 30);
    Point3Cartesian pt(5, 5, 5);
    
    Point3Cartesian translated = pt + offset;
    Point3Cartesian scaled_pt = pt * 2.0;
    Point3Cartesian midpoint = (pt + translated) / 2.0;
    
    std::cout << "  pt + offset = (" << translated.X() << ", " << translated.Y() << ", " << translated.Z() << ")" << std::endl;
    std::cout << "  pt * 2 = (" << scaled_pt.X() << ", " << scaled_pt.Y() << ", " << scaled_pt.Z() << ")" << std::endl;
    std::cout << "  Midpoint = (" << midpoint.X() << ", " << midpoint.Y() << ", " << midpoint.Z() << ")" << std::endl << std::endl;
    
    // Point3Spherical - construction and conversion
    std::cout << "Point3Spherical:" << std::endl;
    // Point on unit sphere at (θ=90°, φ=45°)
    Point3Spherical sph(1.0, Constants::PI/2, Constants::PI/4);
    std::cout << "  Spherical: r=" << sph.R() << ", theta=" << sph.Theta() << ", phi=" << sph.Phi() << std::endl;
    
    Point3Cartesian fromSph = sph.TransfToCart();
    std::cout << "  Converted to Cartesian: (" << fromSph.X() << ", " << fromSph.Y() << ", " << fromSph.Z() << ")" << std::endl << std::endl;
    
    // Point3Cylindrical - construction and conversion
    std::cout << "Point3Cylindrical:" << std::endl;
    // Point at radius 5 from z-axis, angle 30°, height 10
    Point3Cylindrical cyl(5.0, Constants::PI/6, 10.0);
    std::cout << "  Cylindrical: r=" << cyl.R() << ", phi=" << cyl.Phi() << ", z=" << cyl.Z() << std::endl;
    
    Point3Cartesian fromCyl = cyl.TransfToCart();
    std::cout << "  Converted to Cartesian: (" << fromCyl.X() << ", " << fromCyl.Y() << ", " << fromCyl.Z() << ")" << std::endl << std::endl;
}

void Docs_Demo_Geometry()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*****               Geometry - Basic Points           *****" << std::endl;
    std::cout << "***********************************************************" << std::endl << std::endl;
    
    Docs_Demo_Geometry_2D_Points();
    Docs_Demo_Geometry_3D_Points();
}
