///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_computational_geometry.cpp                                ///
///  Description: Brief demonstration of ComputationalGeometry.h                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/algorithms/ComputationalGeometry.h>
#endif

#include <iostream>
#include <limits>
#include <vector>

using namespace MML;

void Docs_Demo_CircleVoronoiContracts()
{
    using namespace MML::CompGeometry;
    const Real scale = Real{1e-6};
    const auto circle = Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{
        {0, 0}, {2 * scale, 0}, {scale, std::sqrt(Real{3}) * scale}});
    const Real offset = Real{1e9};
    const auto translated = Voronoi::ComputeVoronoi({
        {offset, offset}, {offset + 2, offset}, {offset, offset + 2}});
    const auto collinear = Voronoi::ComputeVoronoi({{2, 0}, {0, 0}, {1, 0}});
    if (!IsWithinAbs(circle.Radius() / scale, Real{2} / std::sqrt(Real{3}), Precision::EqualityTolerance) ||
        translated.NumVertices() != 1 || translated.vertices[0].X() != offset + 1 ||
        translated.vertices[0].Y() != offset + 1 || collinear.NumEdges() != 2)
        throw GeometryError("Circle and Voronoi construction demo failed");
    std::cout << "Tiny enclosing radius: " << circle.Radius()
              << "; translated Voronoi center offset: " << translated.vertices[0].X() - offset
              << "; collinear bisectors: " << collinear.NumEdges() << '\n';

    const int exponent = 3 * std::numeric_limits<Real>::max_exponent / 4;
    const Real queryScale = std::scalbn(Real{1}, exponent);
    const auto farSites = Voronoi::ComputeVoronoi({{0, 0}, {queryScale, 0}});
    const int nearest = farSites.FindNearestSite({2 * queryScale, 0});
    const Real spacing = std::scalbn(Real{1}, -exponent);
    const auto tinySites = Voronoi::ComputeVoronoi({{0, 0}, {spacing, 0}, {2 * spacing, 0}});
    const auto middle = tinySites.GetClippedCell(1, -spacing, -spacing, 3 * spacing, spacing);
    if (nearest != 1 || middle.NumVertices() != 4)
        throw GeometryError("Voronoi scaled query and clipping demo failed");
    int left = 0, right = 0;
    for (const auto& vertex : middle.Vertices()) {
        const Real x = vertex.X() / spacing, y = vertex.Y() / spacing;
        if ((x != Real{0.5} && x != Real{1.5}) || std::abs(y) != 1)
            throw GeometryError("Voronoi clipped middle cell violates its independent strip bounds");
        if (x == Real{0.5}) ++left;
        else ++right;
    }
    if (left != 2 || right != 2)
        throw GeometryError("Voronoi clipped middle cell lost a strip boundary");
    std::cout << "Huge-scale nearest site: " << nearest
              << "; tiny clipped middle cell: " << middle.NumVertices() << " vertices\n";
}

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
