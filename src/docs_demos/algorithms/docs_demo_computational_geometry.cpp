///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_computational_geometry.cpp                                ///
///  Description: Brief demonstration of ComputationalGeometry.h                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "algorithms/ComputationalGeometry.h"
#endif

#include <iostream>
#include <vector>

using namespace MML;

void Docs_Demo_ComputationalGeometry()
{
    std::cout << "=== Computational Geometry Demo ===" << std::endl;
    
    // Convex Hull
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0), Point2Cartesian(1, 1), Point2Cartesian(2, 0),
        Point2Cartesian(1, -1), Point2Cartesian(0.5, 0.5), Point2Cartesian(1.5, 0.5)
    };
    
    std::cout << "Computing convex hull of 6 points..." << std::endl;
    Polygon2D hull = CompGeometry::ConvexHull2D::Compute(points);
    std::cout << "Hull has " << hull.NumVertices() << " vertices" << std::endl;
    
    // Segment intersection
    Point2Cartesian p1(0, 0), q1(2, 2);
    Point2Cartesian p2(0, 2), q2(2, 0);
    
    SegmentIntersection inter = CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    std::cout << "\nSegment (0,0)-(2,2) with (0,2)-(2,0):" << std::endl;
    std::cout << "  Intersects: " << (inter.type == SegmentIntersectionType::Point ? "Yes" : "No") << std::endl;
    if (inter.type == SegmentIntersectionType::Point)
        std::cout << "  At: (" << inter.point.X() << ", " << inter.point.Y() << ")" << std::endl;
    
    // Point in polygon (using Polygon2D::Contains method)
    Polygon2D square({Point2Cartesian(0,0), Point2Cartesian(2,0), 
                      Point2Cartesian(2,2), Point2Cartesian(0,2)});
    Point2Cartesian inside(1, 1);
    Point2Cartesian outside(3, 3);
    
    std::cout << "\nPoint in polygon (square 0-2 x 0-2):" << std::endl;
    std::cout << "  (1,1) inside: " << (square.Contains(inside) ? "Yes" : "No") << std::endl;
    std::cout << "  (3,3) inside: " << (square.Contains(outside) ? "Yes" : "No") << std::endl;
    
    // Polygon area (using Polygon2D::Area method)
    Real area = square.Area();
    std::cout << "  Square area: " << area << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}

