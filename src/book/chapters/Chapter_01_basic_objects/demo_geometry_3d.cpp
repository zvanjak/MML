///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_geometry_3d.cpp
///  Description: Comprehensive 3D Geometry demonstrations for Chapter 01
///               Points, vectors, lines, planes, and coordinate systems
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Geometry/Geometry.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/BaseUtils.h"
#include "mml/core/CoordTransf/CoordTransfSpherical.h"
#include "mml/core/CoordTransf/CoordTransfCylindrical.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        3D Points - Multiple Coordinate Systems                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_3D_Points()
{
    std::cout << "\n=== 3D Point Types ===\n\n";

    // Cartesian coordinates (x, y, z)
    std::cout << "--- Cartesian Points ---\n";
    Point3Cartesian origin;
    Point3Cartesian p1(1.0, 2.0, 3.0);
    Point3Cartesian p2(4.0, 5.0, 6.0);

    std::cout << "origin: (" << origin.X() << ", " << origin.Y() << ", " << origin.Z() << ")\n";
    std::cout << "p1: (" << p1.X() << ", " << p1.Y() << ", " << p1.Z() << ")\n";
    std::cout << "Distance p1 to p2: " << p1.Dist(p2) << "\n";

    // Spherical coordinates (r, θ, φ)
    std::cout << "\n--- Spherical Points (r, θ, φ) ---\n";
    Point3Spherical sph1(1.0, Constants::PI/2, 0.0);      // On X-axis
    Point3Spherical sph2(1.0, Constants::PI/2, Constants::PI/4);  // 45° in XY plane
    Point3Spherical sph3(1.0, 0.0, 0.0);                  // On Z-axis (north pole)

    std::cout << "sph1 (r=1, θ=π/2, φ=0): R=" << sph1.R() 
              << ", Theta=" << sph1.Theta() << ", Phi=" << sph1.Phi() << "\n";
    std::cout << "sph2 (r=1, θ=π/2, φ=π/4): 45° in XY plane\n";
    std::cout << "sph3 (r=1, θ=0, φ=0): North pole\n";

    // Cylindrical coordinates (r, φ, z)
    std::cout << "\n--- Cylindrical Points (r, φ, z) ---\n";
    Point3Cylindrical cyl1(2.0, 0.0, 5.0);                  // On X-axis, z=5
    Point3Cylindrical cyl2(2.0, Constants::PI/3, 5.0);      // 60° in XY, z=5

    std::cout << "cyl1 (r=2, φ=0, z=5): R=" << cyl1.R() 
              << ", Phi=" << cyl1.Phi() << ", Z=" << cyl1.Z() << "\n";
    std::cout << "cyl2 (r=2, φ=60°, z=5): 60° rotation in XY\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Coordinate Transformations                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_3D_CoordTransform()
{
    std::cout << "\n=== Coordinate Transformations ===\n\n";

    // Coordinate transformations work with Vectors (position vectors from origin)
    // Start with a Cartesian vector
    Vector3Cartesian cart(1.0, 1.0, 1.0);
    std::cout << "Original Cartesian vector: " << cart << "\n\n";

    // Cartesian → Spherical
    std::cout << "--- Cartesian → Spherical ---\n";
    Vector3Spherical sph = CoordTransfCartToSpher.transf(cart);
    
    std::cout << "  R     = " << sph[0] << " (distance from origin = √3 ≈ 1.732)\n";
    std::cout << "  Theta = " << sph[1] << " rad (" 
              << sph[1] * 180/Constants::PI << "°) (polar angle)\n";
    std::cout << "  Phi   = " << sph[2] << " rad (" 
              << sph[2] * 180/Constants::PI << "°) (azimuthal angle = 45°)\n";

    // Spherical → Cartesian (verify roundtrip)
    std::cout << "\n--- Spherical → Cartesian (roundtrip) ---\n";
    Vector3Cartesian cart2 = CoordTransfSpherToCart.transf(sph);
    std::cout << "  Back to Cartesian: " << cart2 << "\n";

    // Cartesian → Cylindrical
    std::cout << "\n--- Cartesian → Cylindrical ---\n";
    Vector3Cylindrical cyl = CoordTransfCartToCyl.transf(cart);
    
    std::cout << "  R   = " << cyl[0] << " (radial distance in XY = √2 ≈ 1.414)\n";
    std::cout << "  Phi = " << cyl[1] << " rad (" 
              << cyl[1] * 180/Constants::PI << "°)\n";
    std::cout << "  Z   = " << cyl[2] << "\n";

    // Practical example: Earth surface point
    std::cout << "\n--- Practical: Earth Coordinates ---\n";
    double R_earth = 6371.0;  // km
    double lat = 40.7 * Constants::PI / 180.0;   // NYC latitude in radians
    double lon = -74.0 * Constants::PI / 180.0;  // NYC longitude in radians (west is negative)
    
    // Spherical coords: (r, theta, phi) where theta is polar angle from Z-axis
    Vector3Spherical nyc_sph(R_earth, Constants::PI/2 - lat, lon);
    Vector3Cartesian nyc_cart = CoordTransfSpherToCart.transf(nyc_sph);
    
    std::cout << "NYC (approx): lat=40.7°N, lon=74°W\n";
    std::cout << "  Spherical: R=" << R_earth << "km, θ=" << (Constants::PI/2 - lat) 
              << " rad, φ=" << lon << " rad\n";
    std::cout << "  Cartesian: (" << nyc_cart[0] << ", " << nyc_cart[1] << ", " << nyc_cart[2] << ") km\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        3D Vectors                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_3D_Vectors()
{
    std::cout << "\n=== 3D Vectors ===\n\n";

    // Construction
    std::cout << "--- Vector Construction ---\n";
    Vector3Cartesian v1;                              // Default: (0, 0, 0)
    Vector3Cartesian v2(1.0, 2.0, 3.0);               // From components
    Vector3Cartesian v3(Point3Cartesian(0.0, 0.0, 0.0), 
                        Point3Cartesian(1.0, 1.0, 1.0));  // From two points

    std::cout << "v1 (default): " << v1 << "\n";
    std::cout << "v2 (1, 2, 3): " << v2 << "\n";
    std::cout << "v3 (origin to (1,1,1)): " << v3 << "\n";

    // Standard basis vectors
    std::cout << "\n--- Standard Basis ---\n";
    Vector3Cartesian i_hat(1.0, 0.0, 0.0);
    Vector3Cartesian j_hat(0.0, 1.0, 0.0);
    Vector3Cartesian k_hat(0.0, 0.0, 1.0);
    std::cout << "î = " << i_hat << "\n";
    std::cout << "ĵ = " << j_hat << "\n";
    std::cout << "k̂ = " << k_hat << "\n";

    // Operations
    std::cout << "\n--- Vector Operations ---\n";
    Vector3Cartesian a(1.0, 2.0, 3.0);
    Vector3Cartesian b(4.0, 5.0, 6.0);
    
    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "a + b = " << (a + b) << "\n";
    std::cout << "a - b = " << (a - b) << "\n";
    std::cout << "2.0 * a = " << (2.0 * a) << "\n";
    std::cout << "a.NormL2() = " << a.NormL2() << " (√14 ≈ 3.74)\n";

    // Dot product
    std::cout << "\n--- Dot Product ---\n";
    std::cout << "a · b = " << ScalarProduct(a, b) << " (1*4 + 2*5 + 3*6 = 32)\n";

    // Cross product
    std::cout << "\n--- Cross Product ---\n";
    auto cross = VectorProduct(i_hat, j_hat);
    std::cout << "î × ĵ = " << cross << " (should be k̂)\n";
    
    auto ab_cross = VectorProduct(a, b);
    std::cout << "a × b = " << ab_cross << "\n";

    // Parallel and perpendicular checks
    std::cout << "\n--- Geometric Relationships ---\n";
    Vector3Cartesian parallel(2.0, 4.0, 6.0);  // 2*a
    std::cout << "a = " << a << "\n";
    std::cout << "parallel = " << parallel << " (= 2a)\n";
    std::cout << "a.IsParallelTo(parallel): " << (a.IsParallelTo(parallel) ? "true" : "false") << "\n";
    std::cout << "î.IsPerpendicularTo(ĵ): " << (i_hat.IsPerpendicularTo(j_hat) ? "true" : "false") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Lines in 3D                                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_3D_Lines()
{
    std::cout << "\n=== Lines in 3D ===\n\n";

    // Coordinate axes
    std::cout << "--- Coordinate Axes ---\n";
    Line3D x_axis(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(1.0, 0.0, 0.0));
    Line3D y_axis(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(0.0, 1.0, 0.0));
    Line3D z_axis(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(0.0, 0.0, 1.0));

    std::cout << "X-axis: origin, direction (1, 0, 0)\n";
    std::cout << "Y-axis: origin, direction (0, 1, 0)\n";
    std::cout << "Z-axis: origin, direction (0, 0, 1)\n";

    // General lines
    std::cout << "\n--- General Lines ---\n";
    Line3D line1(Point3Cartesian(1.0, 2.0, 3.0), Vector3Cartesian(1.0, 1.0, 1.0));
    Line3D line2(Point3Cartesian(0.0, 0.0, 1.0), Point3Cartesian(1.0, 1.0, 2.0));

    std::cout << "line1: through (1,2,3), direction (1,1,1)\n";
    std::cout << "line2: through (0,0,1) and (1,1,2)\n";

    // Line segments
    std::cout << "\n--- Line Segments ---\n";
    SegmentLine3D seg1(Point3Cartesian(0.0, 0.0, 0.0), Point3Cartesian(1.0, 1.0, 1.0));
    SegmentLine3D seg2(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(0.0, 1.0, 0.0), 5.0);

    std::cout << "seg1: from origin to (1,1,1), length = " << seg1.Length() << " (√3 ≈ 1.73)\n";
    std::cout << "seg2: from origin, direction (0,1,0), length = " << seg2.Length() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Planes in 3D                                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_3D_Planes()
{
    std::cout << "\n=== Planes in 3D ===\n\n";

    // Different ways to define planes
    std::cout << "--- Plane Construction ---\n";
    
    // From point and normal
    Plane3D xy_plane(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(0.0, 0.0, 1.0));
    std::cout << "XY plane: through origin, normal (0, 0, 1)\n";
    
    // From three points
    Plane3D plane2(Point3Cartesian(1.0, 0.0, 0.0), 
                   Point3Cartesian(0.0, 1.0, 0.0), 
                   Point3Cartesian(0.0, 0.0, 1.0));
    std::cout << "plane2: through (1,0,0), (0,1,0), (0,0,1)\n";
    
    // From equation ax + by + cz + d = 0
    Plane3D plane3(1.0, 2.0, 3.0, 4.0);  // x + 2y + 3z + 4 = 0
    std::cout << "plane3: x + 2y + 3z + 4 = 0\n";

    // Coordinate planes
    Plane3D xz_plane(1.0, 3.0, 2.0);  // Using normal direction shorthand
    std::cout << "xz_plane: defined with direction (1, 3, 2)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Plane Operations                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Plane_Operations()
{
    std::cout << "\n=== Plane Operations ===\n\n";

    Plane3D plane(1.0, 3.0, -2.0, 4.0);  // x + 3y - 2z + 4 = 0
    Point3Cartesian test_point(1.0, 2.0, 3.0);

    std::cout << "Plane: x + 3y - 2z + 4 = 0\n";
    std::cout << "Test point: (" << test_point.X() << ", " << test_point.Y() << ", " << test_point.Z() << ")\n\n";

    // Point-plane relationships
    std::cout << "--- Point-Plane Relationships ---\n";
    std::cout << "Distance to plane: " << plane.DistToPoint(test_point) << "\n";
    std::cout << "Is point on plane: " << (plane.IsPointOnPlane(test_point) ? "yes" : "no") << "\n";
    
    auto projection = plane.ProjectionToPlane(test_point);
    std::cout << "Projection onto plane: (" << projection.X() << ", " << projection.Y() << ", " << projection.Z() << ")\n";

    // Line-plane relationships
    std::cout << "\n--- Line-Plane Relationships ---\n";
    Line3D test_line(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(1.0, 1.0, 1.0));
    
    std::cout << "Test line: through origin, direction (1,1,1)\n";
    std::cout << "Is line on plane: " << (plane.IsLineOnPlane(test_line) ? "yes" : "no") << "\n";
    std::cout << "Angle to line: " << plane.AngleToLine(test_line) << " rad\n";

    Point3Cartesian intersection;
    bool intersects = plane.IntersectionWithLine(test_line, intersection);
    if (intersects) {
        std::cout << "Intersection point: (" << intersection.X() << ", " 
                  << intersection.Y() << ", " << intersection.Z() << ")\n";
    }

    // Plane-plane relationships
    std::cout << "\n--- Plane-Plane Relationships ---\n";
    Plane3D plane2(2.0, 6.0, -4.0, 8.0);  // 2x + 6y - 4z + 8 = 0 (parallel to first)
    Plane3D plane3(1.0, 0.0, 0.0, -1.0);  // x - 1 = 0 (YZ plane shifted)

    std::cout << "plane2: 2x + 6y - 4z + 8 = 0 (parallel to first?)\n";
    std::cout << "IsParallelToPlane: " << (plane.IsParallelToPlane(plane2) ? "yes" : "no") << "\n";
    
    std::cout << "\nplane3: x = 1 (YZ plane shifted)\n";
    std::cout << "IsPerpendicularToPlane: " << (plane.IsPerpendicularToPlane(plane3) ? "yes" : "no") << "\n";
    std::cout << "Angle between planes: " << plane.AngleToPlane(plane3) << " rad\n";

    Line3D intersection_line;
    if (plane.IntersectionWithPlane(plane3, intersection_line)) {
        std::cout << "Planes intersect along a line\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_3D_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // Computer Graphics: Camera lookAt direction
    std::cout << "--- Graphics: Camera Direction ---\n";
    Point3Cartesian camera_pos(0.0, 5.0, 10.0);
    Point3Cartesian target(0.0, 0.0, 0.0);
    
    Vector3Cartesian look_dir(camera_pos, target);
    double dist = camera_pos.Dist(target);
    
    std::cout << "Camera at: (" << camera_pos.X() << ", " << camera_pos.Y() << ", " << camera_pos.Z() << ")\n";
    std::cout << "Looking at: (" << target.X() << ", " << target.Y() << ", " << target.Z() << ")\n";
    std::cout << "Look direction: " << look_dir << "\n";
    std::cout << "Distance to target: " << dist << "\n";

    // Physics: Magnetic field around wire (cross product)
    std::cout << "\n--- Physics: Magnetic Field (Biot-Savart) ---\n";
    Vector3Cartesian current_dir(0.0, 0.0, 1.0);  // Current in Z direction
    Vector3Cartesian r(1.0, 0.0, 0.0);            // Point at distance r on X-axis
    
    auto B_dir = VectorProduct(current_dir, r);   // B ∝ I × r̂
    std::cout << "Current direction: " << current_dir << " (along Z)\n";
    std::cout << "Observation point: " << r << " (on X-axis)\n";
    std::cout << "B field direction: " << B_dir << " (along -Y, circles the wire)\n";

    // Engineering: Distance from point to plane (clearance check)
    std::cout << "\n--- Engineering: Clearance Check ---\n";
    Plane3D floor(0.0, 0.0, 1.0, 0.0);  // z = 0 (floor plane)
    Point3Cartesian cable_point(5.0, 3.0, 2.5);  // Cable hanging point
    
    double clearance = floor.DistToPoint(cable_point);
    std::cout << "Floor plane: z = 0\n";
    std::cout << "Cable at: (" << cable_point.X() << ", " << cable_point.Y() << ", " << cable_point.Z() << ")\n";
    std::cout << "Clearance above floor: " << clearance << " m\n";

    // Robotics: Check if points are coplanar
    std::cout << "\n--- Robotics: Coplanarity Check ---\n";
    Point3Cartesian p1(0.0, 0.0, 0.0);
    Point3Cartesian p2(1.0, 0.0, 0.0);
    Point3Cartesian p3(0.0, 1.0, 0.0);
    Point3Cartesian p4(1.0, 1.0, 0.0);  // All in Z=0 plane
    
    Plane3D work_surface(p1, p2, p3);
    bool p4_on_surface = work_surface.IsPointOnPlane(p4);
    std::cout << "Points p1, p2, p3, p4 - all in z=0 plane?\n";
    std::cout << "p4 on work surface: " << (p4_on_surface ? "yes (coplanar)" : "no") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Geometry_3D()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                     3D GEOMETRY IN MML                        ****\n";
    std::cout << "****    Points, Vectors, Lines, Planes, and Transformations        ****\n";
    std::cout << "***********************************************************************\n";

    Demo_3D_Points();
    Demo_3D_CoordTransform();
    Demo_3D_Vectors();
    Demo_3D_Lines();
    Demo_3D_Planes();
    Demo_Plane_Operations();
    Demo_3D_Applications();
}