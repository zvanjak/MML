#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Geometry2D.h"
#include "base/Geometry3D.h"
#endif

using namespace MML;

void Docs_Demo_2D_Lines()
{
    std::cout << "=== 2D Lines and Segments ===" << std::endl << std::endl;
    
    // Line2D - construction and basic operations
    std::cout << "Line2D:" << std::endl;
    Point2Cartesian p(0, 0);
    Vector2Cartesian dir(1, 1);
    Line2D line1(p, dir);  // Point + direction
    
    Point2Cartesian a(0, 0), b(2, 2);
    Line2D line2(a, b);    // Two points
    
    std::cout << "  line1 start: (" << line1.StartPoint().X() << ", " << line1.StartPoint().Y() << ")" << std::endl;
    std::cout << "  line1 direction: (" << line1.Direction().X() << ", " << line1.Direction().Y() << ")" << std::endl;
    
    // Parametric evaluation
    Point2Cartesian p_at_2 = line1(2.0);
    std::cout << "  Point at t=2: (" << p_at_2.X() << ", " << p_at_2.Y() << ")" << std::endl << std::endl;
    
    // SegmentLine2D
    std::cout << "SegmentLine2D:" << std::endl;
    Point2Cartesian seg_a(0, 0), seg_b(4, 3);
    SegmentLine2D seg1(seg_a, seg_b);
    
    std::cout << "  Start: (" << seg1.StartPoint().X() << ", " << seg1.StartPoint().Y() << ")" << std::endl;
    std::cout << "  End: (" << seg1.EndPoint().X() << ", " << seg1.EndPoint().Y() << ")" << std::endl;
    std::cout << "  Length: " << seg1.Length() << std::endl;
    
    // Midpoint
    Point2Cartesian mid = seg1.PointOnSegment(0.5);
    std::cout << "  Midpoint: (" << mid.X() << ", " << mid.Y() << ")" << std::endl << std::endl;
}

void Docs_Demo_2D_Triangle()
{
    std::cout << "=== 2D Triangle ===" << std::endl << std::endl;
    
    // Construction
    Point2Cartesian p1(0, 0), p2(3, 0), p3(0, 4);
    Triangle2D tri(p1, p2, p3);
    
    std::cout << "Right triangle (3-4-5):" << std::endl;
    std::cout << "  Side A: " << tri.A() << std::endl;
    std::cout << "  Side B: " << tri.B() << std::endl;
    std::cout << "  Side C: " << tri.C() << std::endl;
    std::cout << "  Area: " << tri.Area() << std::endl;
    std::cout << "  Is right: " << (tri.IsRight() ? "yes" : "no") << std::endl;
    std::cout << "  Is isosceles: " << (tri.IsIsosceles() ? "yes" : "no") << std::endl << std::endl;
}

void Docs_Demo_2D_Polygon()
{
    std::cout << "=== 2D Polygon ===" << std::endl << std::endl;
    
    // Square with side length 2
    Polygon2D square({
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    });
    
    std::cout << "Square (side=2):" << std::endl;
    std::cout << "  Vertices: " << square.NumVertices() << std::endl;
    std::cout << "  Area: " << square.Area() << std::endl;
    std::cout << "  Perimeter: " << square.Perimeter() << std::endl << std::endl;
}

void Docs_Demo_3D_Lines()
{
    std::cout << "=== 3D Lines and Segments ===" << std::endl << std::endl;
    
    // Line3D - construction
    Point3Cartesian p(0, 0, 0);
    Vector3Cartesian dir(1, 0, 0);
    Line3D line1(p, dir);  // x-axis
    
    std::cout << "Line3D (x-axis):" << std::endl;
    std::cout << "  Start: (" << line1.StartPoint().X() << ", " << line1.StartPoint().Y() << ", " << line1.StartPoint().Z() << ")" << std::endl;
    std::cout << "  Direction: (" << line1.Direction().X() << ", " << line1.Direction().Y() << ", " << line1.Direction().Z() << ")" << std::endl;
    
    // Point at t=5
    Point3Cartesian p_at_5 = line1(5.0);
    std::cout << "  Point at t=5: (" << p_at_5.X() << ", " << p_at_5.Y() << ", " << p_at_5.Z() << ")" << std::endl;
    
    // Distance from point to line
    Point3Cartesian test_pt(0, 3, 4);
    Real dist = line1.Dist(test_pt);
    std::cout << "  Distance from (0,3,4) to line: " << dist << std::endl;
    
    // Nearest point on line
    Point3Cartesian nearest = line1.NearestPointOnLine(test_pt);
    std::cout << "  Nearest point: (" << nearest.X() << ", " << nearest.Y() << ", " << nearest.Z() << ")" << std::endl << std::endl;
    
    // SegmentLine3D
    std::cout << "SegmentLine3D:" << std::endl;
    Point3Cartesian seg_a(0, 0, 0), seg_b(1, 1, 1);
    SegmentLine3D seg(seg_a, seg_b);
    
    std::cout << "  Length: " << seg.Length() << std::endl;
    Point3Cartesian seg_pt = seg(0.5);
    std::cout << "  Point at t=0.5: (" << seg_pt.X() << ", " << seg_pt.Y() << ", " << seg_pt.Z() << ")" << std::endl << std::endl;
}

void Docs_Demo_3D_Plane()
{
    std::cout << "=== 3D Plane ===" << std::endl << std::endl;
    
    // XY plane (z = 0)
    Plane3D plane = Plane3D::GetXYPlane();
    
    std::cout << "XY Plane (z=0):" << std::endl;
    std::cout << "  Equation: " << plane.A() << "x + " << plane.B() << "y + " << plane.C() << "z + " << plane.D() << " = 0" << std::endl;
    
    // Point on plane check
    Point3Cartesian on_plane(1, 2, 0);
    Point3Cartesian off_plane(1, 2, 3);
    std::cout << "  (1,2,0) on plane: " << (plane.IsPointOnPlane(on_plane) ? "yes" : "no") << std::endl;
    std::cout << "  (1,2,3) on plane: " << (plane.IsPointOnPlane(off_plane) ? "yes" : "no") << std::endl;
    std::cout << "  Distance to (1,2,3): " << plane.DistToPoint(off_plane) << std::endl;
    
    // Line-plane intersection
    std::cout << std::endl << "Line-Plane Intersection:" << std::endl;
    Line3D line(Point3Cartesian(0, 0, 5), Vector3Cartesian(0, 0, -1));
    Point3Cartesian intersection;
    if (plane.IntersectionWithLine(line, intersection)) {
        std::cout << "  Intersection at: (" << intersection.X() << ", " << intersection.Y() << ", " << intersection.Z() << ")" << std::endl;
    }
    
    // Plane from three points
    std::cout << std::endl << "Plane from 3 points:" << std::endl;
    Point3Cartesian p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1);
    Plane3D plane2(p1, p2, p3);
    std::cout << "  Equation: " << plane2.A() << "x + " << plane2.B() << "y + " << plane2.C() << "z + " << plane2.D() << " = 0" << std::endl;
    
    // Plane-plane intersection
    std::cout << std::endl << "Plane-Plane Intersection:" << std::endl;
    Plane3D yz_plane = Plane3D::GetYZPlane();
    Line3D inter_line;
    if (plane.IntersectionWithPlane(yz_plane, inter_line)) {
        std::cout << "  XY and YZ planes intersect along y-axis" << std::endl;
        std::cout << "  Direction: (" << inter_line.Direction().X() << ", " << inter_line.Direction().Y() << ", " << inter_line.Direction().Z() << ")" << std::endl;
    }
    std::cout << std::endl;
}

void Docs_Demo_3D_Triangle()
{
    std::cout << "=== 3D Triangle ===" << std::endl << std::endl;
    
    // Right triangle in XY plane
    Point3Cartesian p1(0, 0, 0), p2(3, 0, 0), p3(0, 4, 0);
    Triangle3D tri(p1, p2, p3);
    
    std::cout << "Right triangle (3-4-5) in XY plane:" << std::endl;
    std::cout << "  Area: " << tri.Area() << std::endl;
    std::cout << "  Is right: " << (tri.IsRight() ? "yes" : "no") << std::endl;
    
    Point3Cartesian centroid = tri.Centroid();
    std::cout << "  Centroid: (" << centroid.X() << ", " << centroid.Y() << ", " << centroid.Z() << ")" << std::endl << std::endl;
}

void Docs_Demo_Geometry_2D_3D()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*****           Geometry 2D & 3D Primitives           *****" << std::endl;
    std::cout << "***********************************************************" << std::endl << std::endl;
    
    Docs_Demo_2D_Lines();
    Docs_Demo_2D_Triangle();
    Docs_Demo_2D_Polygon();
    Docs_Demo_3D_Lines();
    Docs_Demo_3D_Plane();
    Docs_Demo_3D_Triangle();
}