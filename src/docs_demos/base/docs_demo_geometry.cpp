///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_geometry.cpp
///  Description: Demo for Geometry.h and related headers - Geometric primitives
///////////////////////////////////////////////////////////////////////////////////////////
#include "base/Geometry/Geometry.h"
#include <iostream>
#include <iomanip>

using namespace MML;

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 2D Point classes
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Points2D()
{
    std::cout << "=== 2D Points ===" << std::endl;
    
    // Cartesian 2D points
    Point2Cartesian p1;              // Default: origin
    Point2Cartesian p2(3.0, 4.0);    // From coordinates
    Point2Cartesian p3(1.0, 1.0);
    
    std::cout << "Point2Cartesian:" << std::endl;
    std::cout << "  p1 (origin): " << p1 << std::endl;
    std::cout << "  p2: " << p2 << std::endl;
    std::cout << "  p2.X() = " << p2.X() << ", p2.Y() = " << p2.Y() << std::endl;
    
    // Distance
    std::cout << "\nDistance p1 to p2: " << p1.Dist(p2) << std::endl;
    
    // Arithmetic
    std::cout << "\nPoint arithmetic:" << std::endl;
    std::cout << "  p2 + p3 = " << (p2 + p3) << std::endl;
    std::cout << "  p2 - p3 = " << (p2 - p3) << std::endl;
    std::cout << "  p2 * 2.0 = " << (p2 * 2.0) << std::endl;
    std::cout << "  p2 / 2.0 = " << (p2 / 2.0) << std::endl;
    
    // Polar 2D points
    Point2Polar pol(5.0, Constants::PI / 4);  // r=5, theta=45°
    std::cout << "\nPoint2Polar (r=5, θ=π/4):" << std::endl;
    std::cout << "  R() = " << pol.R() << ", Phi() = " << pol.Phi() << " rad" << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 3D Point classes
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Points3D()
{
    std::cout << "=== 3D Points ===" << std::endl;
    
    // Cartesian 3D
    Point3Cartesian p1;
    Point3Cartesian p2(1.0, 2.0, 3.0);
    Point3Cartesian p3(4.0, 5.0, 6.0);
    
    std::cout << "Point3Cartesian:" << std::endl;
    std::cout << "  p2: " << p2 << std::endl;
    std::cout << "  X=" << p2.X() << ", Y=" << p2.Y() << ", Z=" << p2.Z() << std::endl;
    
    // Distance
    std::cout << "  Distance p2 to p3: " << p2.Dist(p3) << std::endl;
    
    // Spherical 3D (r, theta, phi)
    Point3Spherical sph(1.0, Constants::PI / 2, Constants::PI / 4);
    std::cout << "\nPoint3Spherical (r=1, θ=π/2, φ=π/4):" << std::endl;
    std::cout << "  R=" << sph.R() << ", Theta=" << sph.Theta() << ", Phi=" << sph.Phi() << std::endl;
    
    // Cylindrical 3D (r, phi, z)
    Point3Cylindrical cyl(2.0, Constants::PI / 3, 5.0);
    std::cout << "\nPoint3Cylindrical (r=2, φ=π/3, z=5):" << std::endl;
    std::cout << "  R=" << cyl.R() << ", Phi=" << cyl.Phi() << ", Z=" << cyl.Z() << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: Coordinate transformations
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_CoordTransform()
{
    std::cout << "=== Coordinate Transformations ===" << std::endl;
    
    // Cartesian to Spherical
    Point3Cartesian cart(1.0, 1.0, 1.0);
    Point3Spherical sph;
    CoordTransfCart2Spher(cart, sph);
    
    std::cout << "Cartesian (1, 1, 1) → Spherical:" << std::endl;
    std::cout << "  R = " << sph.R() << std::endl;
    std::cout << "  Theta = " << sph.Theta() << " rad (" << sph.Theta() * 180/Constants::PI << "°)" << std::endl;
    std::cout << "  Phi = " << sph.Phi() << " rad (" << sph.Phi() * 180/Constants::PI << "°)" << std::endl;
    
    // Spherical back to Cartesian
    Point3Cartesian cart2;
    CoordTransfSpherToCart(sph, cart2);
    std::cout << "\nBack to Cartesian: " << cart2 << std::endl;
    
    // Cartesian to Cylindrical
    Point3Cylindrical cyl;
    CoordTransfCart2Cyl(cart, cyl);
    std::cout << "\nCartesian (1, 1, 1) → Cylindrical:" << std::endl;
    std::cout << "  R = " << cyl.R() << std::endl;
    std::cout << "  Phi = " << cyl.Phi() << " rad" << std::endl;
    std::cout << "  Z = " << cyl.Z() << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 2D shapes - Triangle
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Triangle()
{
    std::cout << "=== Triangle ===" << std::endl;
    
    // Create triangle from side lengths
    Triangle tri(3.0, 4.0, 5.0);  // 3-4-5 right triangle
    
    std::cout << "Right triangle with sides 3, 4, 5:" << std::endl;
    std::cout << "  Area: " << tri.Area() << std::endl;
    std::cout << "  Perimeter: " << tri.Perimeter() << std::endl;
    std::cout << "  Inradius: " << tri.Inradius() << std::endl;
    std::cout << "  Circumradius: " << tri.Circumradius() << std::endl;
    
    // Triangle types
    std::cout << "  IsRight: " << (tri.IsRight() ? "true" : "false") << std::endl;
    std::cout << "  IsEquilateral: " << (tri.IsEquilateral() ? "true" : "false") << std::endl;
    std::cout << "  IsIsosceles: " << (tri.IsIsosceles() ? "true" : "false") << std::endl;
    
    // Equilateral triangle
    Triangle equi(5.0, 5.0, 5.0);
    std::cout << "\nEquilateral triangle (5, 5, 5):" << std::endl;
    std::cout << "  Area: " << equi.Area() << std::endl;
    std::cout << "  IsEquilateral: " << (equi.IsEquilateral() ? "true" : "false") << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 2D shapes - Circle and related
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Circle()
{
    std::cout << "=== Circle and Related Shapes ===" << std::endl;
    
    // Circle
    Circle circ(5.0);  // radius = 5
    std::cout << "Circle (r=5):" << std::endl;
    std::cout << "  Radius: " << circ.Radius() << std::endl;
    std::cout << "  Diameter: " << circ.Diameter() << std::endl;
    std::cout << "  Circumference: " << circ.Circumference() << std::endl;
    std::cout << "  Area: " << circ.Area() << std::endl;
    
    // Ellipse
    Ellipse ellipse(5.0, 3.0);  // a=5 (semi-major), b=3 (semi-minor)
    std::cout << "\nEllipse (a=5, b=3):" << std::endl;
    std::cout << "  Semi-major axis: " << ellipse.SemiMajorAxis() << std::endl;
    std::cout << "  Semi-minor axis: " << ellipse.SemiMinorAxis() << std::endl;
    std::cout << "  Area: " << ellipse.Area() << std::endl;
    std::cout << "  Eccentricity: " << ellipse.Eccentricity() << std::endl;
    std::cout << "  Perimeter (approx): " << ellipse.Perimeter() << std::endl;
    
    // Circular sector
    CircularSector sector(5.0, Constants::PI / 4);  // r=5, angle=45°
    std::cout << "\nCircular Sector (r=5, θ=π/4):" << std::endl;
    std::cout << "  Arc length: " << sector.ArcLength() << std::endl;
    std::cout << "  Area: " << sector.Area() << std::endl;
    std::cout << "  Perimeter: " << sector.Perimeter() << std::endl;
    
    // Annulus (ring)
    Annulus ring(5.0, 3.0);  // outer=5, inner=3
    std::cout << "\nAnnulus (R=5, r=3):" << std::endl;
    std::cout << "  Area: " << ring.Area() << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 2D shapes - Polygons
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Polygons()
{
    std::cout << "=== 2D Polygons ===" << std::endl;
    
    // Regular polygon
    RegularPolygon hex(6, 5.0);  // hexagon, circumradius=5
    std::cout << "Regular Hexagon (n=6, R=5):" << std::endl;
    std::cout << "  Number of sides: " << hex.NumSides() << std::endl;
    std::cout << "  Side length: " << hex.SideLength() << std::endl;
    std::cout << "  Area: " << hex.Area() << std::endl;
    std::cout << "  Perimeter: " << hex.Perimeter() << std::endl;
    std::cout << "  Apothem: " << hex.Apothem() << std::endl;
    std::cout << "  Interior angle: " << hex.InteriorAngle() * 180/Constants::PI << "°" << std::endl;
    
    // Rectangle
    Rectangle rect(4.0, 3.0);
    std::cout << "\nRectangle (4 × 3):" << std::endl;
    std::cout << "  Area: " << rect.Area() << std::endl;
    std::cout << "  Perimeter: " << rect.Perimeter() << std::endl;
    std::cout << "  Diagonal: " << rect.Diagonal() << std::endl;
    
    // Parallelogram
    Parallelogram para(5.0, 3.0, Constants::PI / 3);  // base=5, side=3, angle=60°
    std::cout << "\nParallelogram (base=5, side=3, θ=60°):" << std::endl;
    std::cout << "  Area: " << para.Area() << std::endl;
    
    // Rhombus
    Rhombus rhomb(5.0, Constants::PI / 4);  // side=5, angle=45°
    std::cout << "\nRhombus (side=5, θ=45°):" << std::endl;
    std::cout << "  Area: " << rhomb.Area() << std::endl;
    std::cout << "  Perimeter: " << rhomb.Perimeter() << std::endl;
    
    // Trapezoid
    Trapezoid trap(6.0, 4.0, 3.0);  // base1=6, base2=4, height=3
    std::cout << "\nTrapezoid (b1=6, b2=4, h=3):" << std::endl;
    std::cout << "  Area: " << trap.Area() << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: 3D shapes - Sphere, Cylinder, Cone
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_3DShapes()
{
    std::cout << "=== 3D Geometric Shapes ===" << std::endl;
    
    // Sphere
    SphereGeom sphere(5.0);  // radius = 5
    std::cout << "Sphere (r=5):" << std::endl;
    std::cout << "  Radius: " << sphere.Radius() << std::endl;
    std::cout << "  Surface area: " << sphere.SurfaceArea() << std::endl;
    std::cout << "  Volume: " << sphere.Volume() << std::endl;
    
    // Cylinder
    CylinderGeom cyl(3.0, 10.0);  // radius=3, height=10
    std::cout << "\nCylinder (r=3, h=10):" << std::endl;
    std::cout << "  Radius: " << cyl.Radius() << std::endl;
    std::cout << "  Height: " << cyl.Height() << std::endl;
    std::cout << "  Lateral surface area: " << cyl.LateralSurfaceArea() << std::endl;
    std::cout << "  Total surface area: " << cyl.SurfaceArea() << std::endl;
    std::cout << "  Volume: " << cyl.Volume() << std::endl;
    
    // Cone
    ConeGeom cone(3.0, 4.0);  // radius=3, height=4
    std::cout << "\nCone (r=3, h=4):" << std::endl;
    std::cout << "  Slant height: " << cone.SlantHeight() << std::endl;
    std::cout << "  Lateral surface area: " << cone.LateralSurfaceArea() << std::endl;
    std::cout << "  Total surface area: " << cone.SurfaceArea() << std::endl;
    std::cout << "  Volume: " << cone.Volume() << std::endl;
    
    // Frustum (truncated cone)
    Frustum frust(4.0, 2.0, 3.0);  // R=4, r=2, h=3
    std::cout << "\nFrustum (R=4, r=2, h=3):" << std::endl;
    std::cout << "  Slant height: " << frust.SlantHeight() << std::endl;
    std::cout << "  Lateral surface area: " << frust.LateralSurfaceArea() << std::endl;
    std::cout << "  Volume: " << frust.Volume() << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: Equality comparisons with tolerance
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Comparisons()
{
    std::cout << "=== Point Comparisons ===" << std::endl;
    
    Point2Cartesian p1(1.0, 2.0);
    Point2Cartesian p2(1.0, 2.0);
    Point2Cartesian p3(1.0 + 1e-10, 2.0);
    Point2Cartesian p4(1.5, 2.5);
    
    std::cout << "p1 = " << p1 << std::endl;
    std::cout << "p3 = " << p3 << " (tiny offset)" << std::endl;
    
    std::cout << "\nExact comparison:" << std::endl;
    std::cout << "  p1 == p2: " << (p1 == p2 ? "true" : "false") << std::endl;
    std::cout << "  p1 == p3: " << (p1 == p3 ? "true" : "false") << std::endl;
    
    std::cout << "\nTolerant comparison:" << std::endl;
    std::cout << "  p1.IsEqualTo(p3, 1e-9): " << (p1.IsEqualTo(p3, 1e-9) ? "true" : "false") << std::endl;
    std::cout << "  p1.IsEqualTo(p4, 0.1): " << (p1.IsEqualTo(p4, 0.1) ? "true" : "false") << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Demo: Practical application - Distance calculations
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry_Practical()
{
    std::cout << "=== Practical Example: GPS Locations ===" << std::endl;
    
    // Using spherical coordinates for Earth locations
    // Note: Simplified, not geodetic
    Real earthRadius = 6371.0;  // km
    
    // New York (approx): 40.7128°N, 74.0060°W
    Real nyLat = 40.7128 * Constants::PI / 180;
    Real nyLon = -74.0060 * Constants::PI / 180;
    
    // London (approx): 51.5074°N, 0.1278°W
    Real lonLat = 51.5074 * Constants::PI / 180;
    Real lonLon = -0.1278 * Constants::PI / 180;
    
    // Convert to Cartesian for distance calculation
    Point3Cartesian ny(
        earthRadius * std::cos(nyLat) * std::cos(nyLon),
        earthRadius * std::cos(nyLat) * std::sin(nyLon),
        earthRadius * std::sin(nyLat)
    );
    
    Point3Cartesian london(
        earthRadius * std::cos(lonLat) * std::cos(lonLon),
        earthRadius * std::cos(lonLat) * std::sin(lonLon),
        earthRadius * std::sin(lonLat)
    );
    
    Real straightLineDist = ny.Dist(london);
    std::cout << "Straight-line distance NY to London: " << straightLineDist << " km" << std::endl;
    std::cout << "(Great-circle distance would be larger)" << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// @brief Main entry point for Geometry demo
//////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Geometry()
{
    std::cout << "////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "///     Geometry.h - Geometric Primitives Demo           ///" << std::endl;
    std::cout << "////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    Docs_Demo_Geometry_Points2D();
    Docs_Demo_Geometry_Points3D();
    Docs_Demo_Geometry_CoordTransform();
    Docs_Demo_Geometry_Triangle();
    Docs_Demo_Geometry_Circle();
    Docs_Demo_Geometry_Polygons();
    Docs_Demo_Geometry_3DShapes();
    Docs_Demo_Geometry_Comparisons();
    Docs_Demo_Geometry_Practical();
    
    std::cout << "=== Geometry Demo Complete ===" << std::endl;
}
