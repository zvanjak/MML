///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_geometry_2d.cpp
///  Description: Comprehensive 2D Geometry demonstrations for Chapter 01
///               Points, vectors, lines, shapes, and coordinate systems
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Geometry/Geometry.h"
#include "mml/base/Geometry/Geometry2D.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/BaseUtils.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        2D Points - Cartesian and Polar                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Points()
{
    std::cout << "\n=== 2D Point Types ===\n\n";

    // Cartesian coordinates
    std::cout << "--- Cartesian Points ---\n";
    Point2Cartesian origin;              // Default: (0, 0)
    Point2Cartesian p1(3.0, 4.0);        // Explicit coordinates
    Point2Cartesian p2(1.0, 1.0);

    std::cout << "origin: (" << origin.X() << ", " << origin.Y() << ")\n";
    std::cout << "p1: (" << p1.X() << ", " << p1.Y() << ")\n";
    std::cout << "p2: (" << p2.X() << ", " << p2.Y() << ")\n";

    // Distance calculation
    std::cout << "\nDistance from origin to p1: " << origin.Dist(p1) << " (should be 5)\n";
    std::cout << "Distance p1 to p2: " << p1.Dist(p2) << "\n";

    // Point arithmetic (translations)
    std::cout << "\n--- Point Arithmetic ---\n";
    auto sum = p1 + p2;
    auto diff = p1 - p2;
    auto scaled = p1 * 2.0;
    auto divided = p1 / 2.0;
    std::cout << "p1 + p2 = (" << sum.X() << ", " << sum.Y() << ")\n";
    std::cout << "p1 - p2 = (" << diff.X() << ", " << diff.Y() << ")\n";
    std::cout << "p1 * 2.0 = (" << scaled.X() << ", " << scaled.Y() << ")\n";
    std::cout << "p1 / 2.0 = (" << divided.X() << ", " << divided.Y() << ")\n";

    // Polar coordinates
    std::cout << "\n--- Polar Points ---\n";
    Point2Polar pol1(5.0, 0.0);                        // (r=5, θ=0) → (5, 0)
    Point2Polar pol2(5.0, Constants::PI / 4);          // (r=5, θ=45°)
    Point2Polar pol3(5.0, Constants::PI / 2);          // (r=5, θ=90°) → (0, 5)

    std::cout << "pol1 (r=5, θ=0°):   R=" << pol1.R() << ", Phi=" << pol1.Phi() << " rad\n";
    std::cout << "pol2 (r=5, θ=45°):  R=" << pol2.R() << ", Phi=" << pol2.Phi() << " rad\n";
    std::cout << "pol3 (r=5, θ=90°):  R=" << pol3.R() << ", Phi=" << pol3.Phi() << " rad\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        2D Vectors                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Vectors()
{
    std::cout << "\n=== 2D Vectors ===\n\n";

    // Construction
    std::cout << "--- Vector Construction ---\n";
    Vector2Cartesian v1;                          // Default: (0, 0)
    Vector2Cartesian v2(3.0, 4.0);                // From components
    Vector2Cartesian v3(Point2Cartesian(1.0, 1.0), 
                        Point2Cartesian(4.0, 5.0));  // From two points

    std::cout << "v1 (default): " << v1 << "\n";
    std::cout << "v2 (3, 4):    " << v2 << "\n";
    std::cout << "v3 (from points): " << v3 << "\n";

    // Operations
    std::cout << "\n--- Vector Operations ---\n";
    Vector2Cartesian a(1.0, 0.0);
    Vector2Cartesian b(0.0, 1.0);

    std::cout << "a = " << a << ", b = " << b << "\n";
    std::cout << "a + b = " << (a + b) << "\n";
    std::cout << "a - b = " << (a - b) << "\n";
    std::cout << "3.0 * a = " << (3.0 * a) << "\n";
    std::cout << "b / 2.0 = " << (b / 2.0) << "\n";

    // Magnitude and normalization
    std::cout << "\n--- Magnitude and Direction ---\n";
    std::cout << "v2.NormL2() = " << v2.NormL2() << " (should be 5)\n";

    // Dot product
    std::cout << "\n--- Dot Product ---\n";
    Vector2Cartesian v_a(1.0, 2.0);
    Vector2Cartesian v_b(3.0, -3.0);
    std::cout << "v_a = " << v_a << ", v_b = " << v_b << "\n";
    std::cout << "ScalarProduct(v_a, v_b) = " << ScalarProduct(v_a, v_b) << " (1*3 + 2*(-3) = -3)\n";

    // Perpendicular check
    Vector2Cartesian perp1(1.0, 0.0);
    Vector2Cartesian perp2(0.0, 1.0);
    std::cout << "perp1 · perp2 = " << ScalarProduct(perp1, perp2) << " (perpendicular → 0)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Lines and Segments                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Lines()
{
    std::cout << "\n=== Lines and Line Segments ===\n\n";

    // Line2D: infinite line
    std::cout << "--- Infinite Lines ---\n";
    Line2D x_axis(Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0));  // X-axis
    Line2D y_axis(Point2Cartesian(0.0, 0.0), Vector2Cartesian(0.0, 1.0)); // Y-axis
    Line2D diagonal(Point2Cartesian(0.0, 0.0), Vector2Cartesian(1.0, 1.0)); // y = x

    std::cout << "X-axis: through origin, direction (1, 0)\n";
    std::cout << "Y-axis: through origin, direction (0, 1)\n";
    std::cout << "Diagonal: through origin, direction (1, 1)\n";

    // SegmentLine2D: finite segment
    std::cout << "\n--- Line Segments ---\n";
    SegmentLine2D seg1(Point2Cartesian(0.0, 0.0), Point2Cartesian(3.0, 4.0));
    SegmentLine2D seg2(Point2Cartesian(1.0, 1.0), Vector2Cartesian(1.0, 0.0), 5.0);  // Length 5

    std::cout << "seg1: from (0,0) to (3,4), length = " << seg1.Length() << " (should be 5)\n";
    std::cout << "seg2: from (1,1), direction (1,0), length = " << seg2.Length() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        2D Shapes - Triangles                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Triangles()
{
    std::cout << "\n=== Triangles ===\n\n";

    // Create triangles from points
    Triangle2D tri(Point2Cartesian(0.0, 0.0), 
                   Point2Cartesian(4.0, 0.0), 
                   Point2Cartesian(0.0, 3.0));

    std::cout << "Right triangle with vertices at (0,0), (4,0), (0,3):\n";
    std::cout << "  (This is a 3-4-5 right triangle)\n\n";

    // Triangle from side lengths (abstract)
    std::cout << "--- Abstract Triangle (from side lengths) ---\n";
    Triangle abstract_tri(3.0, 4.0, 5.0);  // 3-4-5 right triangle
    
    std::cout << "Triangle with sides 3, 4, 5:\n";
    std::cout << "  Area:          " << abstract_tri.Area() << " (3*4/2 = 6)\n";
    std::cout << "  Perimeter:     " << abstract_tri.Perimeter() << " (3+4+5 = 12)\n";
    std::cout << "  Inradius:      " << abstract_tri.Inradius() << "\n";
    std::cout << "  Circumradius:  " << abstract_tri.Circumradius() << " (hypotenuse/2 = 2.5)\n";

    // Triangle properties
    std::cout << "\n--- Triangle Type Checks ---\n";
    std::cout << "  IsRight:       " << (abstract_tri.IsRight() ? "true" : "false") << "\n";
    std::cout << "  IsEquilateral: " << (abstract_tri.IsEquilateral() ? "true" : "false") << "\n";
    std::cout << "  IsIsosceles:   " << (abstract_tri.IsIsosceles() ? "true" : "false") << "\n";

    // Equilateral triangle
    std::cout << "\n--- Equilateral Triangle ---\n";
    Triangle equi(5.0, 5.0, 5.0);
    std::cout << "Triangle (5, 5, 5):\n";
    std::cout << "  Area: " << equi.Area() << "\n";
    std::cout << "  IsEquilateral: " << (equi.IsEquilateral() ? "true" : "false") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        2D Shapes - Circles and Ellipses                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Circles()
{
    std::cout << "\n=== Circles and Ellipses ===\n\n";

    // Circle
    std::cout << "--- Circle ---\n";
    Circle circ(5.0);  // radius = 5
    
    std::cout << "Circle (r = 5):\n";
    std::cout << "  Radius:        " << circ.Radius() << "\n";
    std::cout << "  Diameter:      " << circ.Diameter() << "\n";
    std::cout << "  Circumference: " << circ.Circumference() << " (2πr ≈ 31.42)\n";
    std::cout << "  Area:          " << circ.Area() << " (πr² ≈ 78.54)\n";

    // Ellipse
    std::cout << "\n--- Ellipse ---\n";
    Ellipse ellipse(5.0, 3.0);  // semi-major a=5, semi-minor b=3
    
    std::cout << "Ellipse (a=5, b=3):\n";
    std::cout << "  Semi-major axis: " << ellipse.SemiMajor() << "\n";
    std::cout << "  Semi-minor axis: " << ellipse.SemiMinor() << "\n";
    std::cout << "  Area:            " << ellipse.Area() << " (πab ≈ 47.12)\n";
    std::cout << "  Eccentricity:    " << ellipse.Eccentricity() << " (0 = circle, 1 = parabola)\n";
    std::cout << "  Perimeter:       " << ellipse.Perimeter() << " (Ramanujan approx)\n";

    // Circular sector (pizza slice)
    std::cout << "\n--- Circular Sector ---\n";
    CircularSector sector(5.0, Constants::PI / 4);  // r=5, angle=45°
    
    std::cout << "Sector (r=5, θ=45°):\n";
    std::cout << "  Arc length: " << sector.ArcLength() << "\n";
    std::cout << "  Area:       " << sector.Area() << "\n";
    std::cout << "  Perimeter:  " << sector.Perimeter() << " (arc + 2r)\n";

    // Annulus (ring/washer)
    std::cout << "\n--- Annulus (Ring) ---\n";
    Annulus ring(5.0, 3.0);  // outer R=5, inner r=3
    
    std::cout << "Annulus (R=5, r=3):\n";
    std::cout << "  Area: " << ring.Area() << " (π(R²-r²) = π(25-9) ≈ 50.27)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        2D Shapes - Polygons                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Polygons()
{
    std::cout << "\n=== Polygons ===\n\n";

    // Regular polygons
    std::cout << "--- Regular Polygons ---\n";
    
    RegularPolygon hex(6, 5.0);  // hexagon, circumradius=5
    std::cout << "Regular Hexagon (n=6, R=5):\n";
    std::cout << "  Number of sides:   " << hex.NumSides() << "\n";
    std::cout << "  Side length:       " << hex.SideLength() << "\n";
    std::cout << "  Area:              " << hex.Area() << "\n";
    std::cout << "  Perimeter:         " << hex.Perimeter() << "\n";
    std::cout << "  Apothem:           " << hex.Apothem() << " (inscribed radius)\n";
    std::cout << "  Interior angle:    " << hex.InteriorAngle() * 180/Constants::PI << "°\n";

    RegularPolygon pent(5, 5.0);  // pentagon
    std::cout << "\nRegular Pentagon (n=5, R=5):\n";
    std::cout << "  Area:           " << pent.Area() << "\n";
    std::cout << "  Interior angle: " << pent.InteriorAngle() * 180/Constants::PI << "°\n";

    // Rectangle
    std::cout << "\n--- Rectangle ---\n";
    Rectangle rect(4.0, 3.0);
    std::cout << "Rectangle (4 × 3):\n";
    std::cout << "  Area:      " << rect.Area() << " (12)\n";
    std::cout << "  Perimeter: " << rect.Perimeter() << " (14)\n";
    std::cout << "  Diagonal:  " << rect.Diagonal() << " (5)\n";

    // General polygon from vertices
    std::cout << "\n--- General Polygon ---\n";
    Polygon2D poly({Point2Cartesian(0.0, 0.0), 
                    Point2Cartesian(4.0, 0.0), 
                    Point2Cartesian(4.0, 3.0),
                    Point2Cartesian(0.0, 3.0)});  // Rectangle as polygon
    std::cout << "Polygon with 4 vertices (rectangle shape)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_2D_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // Game development: Collision detection setup
    std::cout << "--- Game Dev: Bounding Circle ---\n";
    Point2Cartesian player_pos(100.0, 200.0);
    double player_radius = 25.0;
    Point2Cartesian enemy_pos(130.0, 220.0);
    double enemy_radius = 20.0;
    
    double dist = player_pos.Dist(enemy_pos);
    double min_dist = player_radius + enemy_radius;
    bool collision = dist < min_dist;
    
    std::cout << "Player at (" << player_pos.X() << ", " << player_pos.Y() << "), radius " << player_radius << "\n";
    std::cout << "Enemy at (" << enemy_pos.X() << ", " << enemy_pos.Y() << "), radius " << enemy_radius << "\n";
    std::cout << "Distance: " << dist << ", Min safe distance: " << min_dist << "\n";
    std::cout << "Collision: " << (collision ? "YES!" : "No") << "\n";

    // Architecture: Room layout
    std::cout << "\n--- Architecture: Room Dimensions ---\n";
    Rectangle room(5.0, 4.0);  // 5m × 4m
    std::cout << "Room: 5m × 4m\n";
    std::cout << "  Floor area:      " << room.Area() << " m²\n";
    std::cout << "  Wall perimeter:  " << room.Perimeter() << " m\n";
    std::cout << "  Diagonal:        " << room.Diagonal() << " m\n";

    // Engineering: Pipe cross-section
    std::cout << "\n--- Engineering: Pipe Cross-section ---\n";
    double outer_diameter = 10.0;  // cm
    double wall_thickness = 0.5;   // cm
    Annulus pipe(outer_diameter/2, outer_diameter/2 - wall_thickness);
    std::cout << "Pipe: outer Ø=" << outer_diameter << "cm, wall=" << wall_thickness << "cm\n";
    std::cout << "  Cross-sectional area of material: " << pipe.Area() << " cm²\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Geometry_2D()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                     2D GEOMETRY IN MML                        ****\n";
    std::cout << "****           Points, Vectors, Lines, and Shapes                  ****\n";
    std::cout << "***********************************************************************\n";

    Demo_2D_Points();
    Demo_2D_Vectors();
    Demo_2D_Lines();
    Demo_2D_Triangles();
    Demo_2D_Circles();
    Demo_2D_Polygons();
    Demo_2D_Applications();
}