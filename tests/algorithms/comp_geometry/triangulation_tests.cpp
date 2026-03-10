/// @file triangulation_tests.cpp
/// @brief Tests for triangulation algorithms
/// @details Tests EarClipTriangulation (simple polygon) and DelaunayTriangulation (Bowyer-Watson)

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include "algorithms/ComputationalGeometry.h"

#include <set>
#include <random>

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::CompGeometry::TriangulationTests {

// ============================================================================
// EAR CLIPPING TRIANGULATION TESTS
// ============================================================================

TEST_CASE("EarClipTriangulation - Triangle", "[ComputationalGeometry][Triangulation]")
{
    Polygon2D triangle({
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    });

    auto triangles = MML::CompGeometry::Triangulation::EarClipTriangulation(triangle);
    
    REQUIRE(triangles.size() == 1);
    REQUIRE_THAT(triangles[0].Area(), WithinAbs(REAL(2.0), REAL(1e-10)));
}

TEST_CASE("EarClipTriangulation - Convex quadrilateral", "[ComputationalGeometry][Triangulation]")
{
    Polygon2D square({
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    });

    auto triangles = MML::CompGeometry::Triangulation::EarClipTriangulation(square);
    
    REQUIRE(triangles.size() == 2);
    
    // Total area should equal original polygon area
    Real totalArea = 0;
    for (const auto& tri : triangles)
        totalArea += tri.Area();
    
    REQUIRE_THAT(totalArea, WithinAbs(REAL(4.0), REAL(1e-10)));
}

TEST_CASE("EarClipTriangulation - Concave polygon (L-shape)", "[ComputationalGeometry][Triangulation]")
{
    // L-shaped polygon (non-convex)
    Polygon2D lShape({
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(0, 2)
    });

    REQUIRE_FALSE(lShape.IsConvex());
    
    auto triangles = MML::CompGeometry::Triangulation::EarClipTriangulation(lShape);
    
    // n-2 triangles for n vertices
    REQUIRE(triangles.size() == 4);
    
    // Total area should equal original polygon area
    Real totalArea = 0;
    for (const auto& tri : triangles)
        totalArea += tri.Area();
    
    Real expectedArea = lShape.Area();
    REQUIRE_THAT(totalArea, WithinAbs(expectedArea, REAL(1e-10)));
}

TEST_CASE("EarClipTriangulation - Arrow/chevron shape", "[ComputationalGeometry][Triangulation]")
{
    // Arrow pointing right (concave)
    Polygon2D arrow({
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 1),
        Point2Cartesian(0, 2),
        Point2Cartesian(1, 1)  // concave vertex
    });

    auto triangles = MML::CompGeometry::Triangulation::EarClipTriangulation(arrow);
    
    REQUIRE(triangles.size() == 2);
    
    Real totalArea = 0;
    for (const auto& tri : triangles)
        totalArea += tri.Area();
    
    REQUIRE_THAT(totalArea, WithinAbs(arrow.Area(), REAL(1e-10)));
}

TEST_CASE("EarClipTriangulation - Pentagon", "[ComputationalGeometry][Triangulation]")
{
    // Regular pentagon
    std::vector<Point2Cartesian> vertices;
    for (int i = 0; i < 5; i++)
    {
        Real angle = 2.0 * Constants::PI * i / 5.0 - Constants::PI / 2.0;
        vertices.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
    }
    Polygon2D pentagon(vertices);

    auto triangles = MML::CompGeometry::Triangulation::EarClipTriangulation(pentagon);
    
    REQUIRE(triangles.size() == 3);
    
    Real totalArea = 0;
    for (const auto& tri : triangles)
        totalArea += tri.Area();
    
    REQUIRE_THAT(totalArea, WithinAbs(pentagon.Area(), REAL(1e-8)));
}

TEST_CASE("EarClipTriangulation - Star shape (highly non-convex)", "[ComputationalGeometry][Triangulation]")
{
    // 5-pointed star
    std::vector<Point2Cartesian> vertices;
    for (int i = 0; i < 5; i++)
    {
        // Outer point
        Real outerAngle = 2.0 * Constants::PI * i / 5.0 - Constants::PI / 2.0;
        vertices.push_back(Point2Cartesian(2.0 * std::cos(outerAngle), 2.0 * std::sin(outerAngle)));
        
        // Inner point
        Real innerAngle = outerAngle + Constants::PI / 5.0;
        vertices.push_back(Point2Cartesian(0.8 * std::cos(innerAngle), 0.8 * std::sin(innerAngle)));
    }
    Polygon2D star(vertices);

    auto triangles = MML::CompGeometry::Triangulation::EarClipTriangulation(star);
    
    // 10-vertex polygon -> 8 triangles
    REQUIRE(triangles.size() == 8);
    
    Real totalArea = 0;
    for (const auto& tri : triangles)
        totalArea += tri.Area();
    
    REQUIRE_THAT(totalArea, WithinAbs(star.Area(), REAL(1e-8)));
}

// ============================================================================
// DELAUNAY TRIANGULATION TESTS (Bowyer-Watson Algorithm)
// ============================================================================

TEST_CASE("DelaunayTriangulation - Simple triangle", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 3);
    REQUIRE(dt.NumTriangles() == 1);
    
    // Single triangle should use all 3 vertices
    auto tri = dt.triangles[0];
    std::set<int> indices(tri.begin(), tri.end());
    REQUIRE(indices.size() == 3);
    REQUIRE(indices.count(0) == 1);
    REQUIRE(indices.count(1) == 1);
    REQUIRE(indices.count(2) == 1);
}

TEST_CASE("DelaunayTriangulation - Square (4 points)", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 4);
    REQUIRE(dt.NumTriangles() == 2);  // 4 points -> 2 triangles
    
    // Total vertex references should be 6 (2 triangles * 3 vertices)
    std::multiset<int> allVertices;
    for (const auto& tri : dt.triangles) {
        for (int v : tri) {
            allVertices.insert(v);
        }
    }
    REQUIRE(allVertices.size() == 6);
}

TEST_CASE("DelaunayTriangulation - Pentagon (5 points)", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(3, 1.5),
        Point2Cartesian(1, 3),
        Point2Cartesian(-1, 1.5)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 5);
    REQUIRE(dt.NumTriangles() == 3);  // 5 points -> 3 triangles
}

TEST_CASE("DelaunayTriangulation - Regular grid", "[ComputationalGeometry][Delaunay]")
{
    // 3x3 grid = 9 points
    std::vector<Point2Cartesian> points;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            points.push_back(Point2Cartesian(i, j));
        }
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 9);
    // For n points in general position: 2n - 2 - k triangles (k = hull vertices)
    // 9 points with 4 hull vertices -> expect ~12 triangles (may vary slightly)
    REQUIRE(dt.NumTriangles() >= 8);
    REQUIRE(dt.NumTriangles() <= 14);
}

TEST_CASE("DelaunayTriangulation - Delaunay property", "[ComputationalGeometry][Delaunay]")
{
    // The key property: no point should be inside any triangle's circumcircle
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(4, 3),
        Point2Cartesian(0, 3),
        Point2Cartesian(2, 1.5)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // Check Delaunay property for each triangle
    for (int t = 0; t < dt.NumTriangles(); t++) {
        const auto& tri = dt.triangles[t];
        const Point2Cartesian& a = dt.points[tri[0]];
        const Point2Cartesian& b = dt.points[tri[1]];
        const Point2Cartesian& c = dt.points[tri[2]];
        
        // Compute circumcircle center and radius
        Real ax = a.X(), ay = a.Y();
        Real bx = b.X(), by = b.Y();
        Real cx = c.X(), cy = c.Y();
        
        Real d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        if (std::abs(d) < 1e-10) continue;  // Degenerate triangle
        
        Real ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d;
        Real uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d;
        Real r2 = (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy);
        
        // Check all other points are outside (or on) the circumcircle
        for (int i = 0; i < dt.NumPoints(); i++) {
            if (i == tri[0] || i == tri[1] || i == tri[2]) continue;
            
            const Point2Cartesian& p = dt.points[i];
            Real dist2 = (p.X() - ux) * (p.X() - ux) + (p.Y() - uy) * (p.Y() - uy);
            
            // Allow small tolerance
            REQUIRE(dist2 >= r2 - 1e-6);
        }
    }
}

TEST_CASE("DelaunayTriangulation - LocatePoint", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(4, 4),
        Point2Cartesian(0, 4),
        Point2Cartesian(2, 2)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // Point at center should be found
    Point2Cartesian center(2, 2);
    int triIdx = dt.LocatePoint(center);
    REQUIRE(triIdx >= 0);
    
    // Point inside convex hull
    Point2Cartesian inside(1, 1);
    triIdx = dt.LocatePoint(inside);
    REQUIRE(triIdx >= 0);
    
    // Point outside convex hull
    Point2Cartesian outside(10, 10);
    triIdx = dt.LocatePoint(outside);
    REQUIRE(triIdx == -1);
}

TEST_CASE("DelaunayTriangulation - Interpolation", "[ComputationalGeometry][Delaunay]")
{
    // Triangle with known function values
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };
    
    // f(x,y) = x + y at vertices: f(0,0)=0, f(2,0)=2, f(1,2)=3
    std::vector<Real> values = {0.0, 2.0, 3.0};
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // Interpolate at centroid (1, 2/3)
    Point2Cartesian centroid(1.0, 2.0/3.0);
    Real interp = dt.Interpolate(centroid, values);
    Real expected = centroid.X() + centroid.Y();  // ~1.667
    REQUIRE_THAT(interp, WithinAbs(expected, 0.1));
    
    // At a vertex
    Real atV0 = dt.Interpolate(points[0], values);
    REQUIRE_THAT(atV0, WithinAbs(0.0, 1e-6));
}

TEST_CASE("DelaunayTriangulation - ConvexHullEdges", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(4, 4),
        Point2Cartesian(0, 4),
        Point2Cartesian(2, 2)  // Interior point
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    auto hullEdges = dt.GetConvexHullEdges();
    
    // Square has 4 hull edges
    REQUIRE(hullEdges.size() == 4);
    
    // All edges should connect adjacent hull vertices (0,1,2,3)
    std::set<int> hullVertices;
    for (const auto& edge : hullEdges) {
        hullVertices.insert(edge.first);
        hullVertices.insert(edge.second);
    }
    // Hull vertices are 0,1,2,3 (vertex 4 is interior)
    REQUIRE(hullVertices.count(4) == 0);  // Interior point not on hull
}

TEST_CASE("DelaunayTriangulation - Neighbor info", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // 2 triangles sharing an edge should be neighbors
    REQUIRE(dt.NumTriangles() == 2);
    
    // Each triangle should have at least one neighbor
    bool hasNeighbor0 = (dt.neighbors[0][0] >= 0 || dt.neighbors[0][1] >= 0 || dt.neighbors[0][2] >= 0);
    bool hasNeighbor1 = (dt.neighbors[1][0] >= 0 || dt.neighbors[1][1] >= 0 || dt.neighbors[1][2] >= 0);
    REQUIRE(hasNeighbor0);
    REQUIRE(hasNeighbor1);
}

TEST_CASE("DelaunayTriangulate - Returns Triangle2D vector", "[ComputationalGeometry][Delaunay]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };
    
    auto triangles = MML::CompGeometry::Triangulation::DelaunayTriangulate(points);
    
    REQUIRE(triangles.size() == 1);
    
    // Check triangle area is positive (CCW orientation)
    Real area = triangles[0].Area();
    REQUIRE(area > 0);
}

TEST_CASE("DelaunayTriangulation - Random points", "[ComputationalGeometry][Delaunay]")
{
    // Generate random points
    std::mt19937 rng(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<Real> dist(0, 10);
    
    std::vector<Point2Cartesian> points;
    for (int i = 0; i < 20; i++) {
        points.push_back(Point2Cartesian(dist(rng), dist(rng)));
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 20);
    REQUIRE(dt.NumTriangles() > 0);
    
    // Basic sanity: triangles should cover the convex hull
    // For n points, we expect roughly 2n - 2 - k triangles
    REQUIRE(dt.NumTriangles() >= 10);
}

// ============================================================================
// DELAUNAY TRIANGULATION - EDGE CASES
// ============================================================================

TEST_CASE("DelaunayTriangulation - Empty and minimal inputs", "[ComputationalGeometry][Delaunay][EdgeCase]")
{
    SECTION("Empty point set")
    {
        std::vector<Point2Cartesian> empty;
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(empty);
        REQUIRE(dt.NumPoints() == 0);
        REQUIRE(dt.NumTriangles() == 0);
    }
    
    SECTION("Single point")
    {
        std::vector<Point2Cartesian> single = {Point2Cartesian(5, 3)};
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(single);
        REQUIRE(dt.NumPoints() == 1);
        REQUIRE(dt.NumTriangles() == 0);
    }
    
    SECTION("Two points")
    {
        std::vector<Point2Cartesian> two = {Point2Cartesian(0, 0), Point2Cartesian(1, 1)};
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(two);
        REQUIRE(dt.NumPoints() == 2);
        REQUIRE(dt.NumTriangles() == 0);
    }
}

TEST_CASE("DelaunayTriangulation - Collinear points", "[ComputationalGeometry][Delaunay][EdgeCase]")
{
    SECTION("Three collinear points")
    {
        std::vector<Point2Cartesian> collinear = {
            Point2Cartesian(0, 0),
            Point2Cartesian(1, 1),
            Point2Cartesian(2, 2)
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(collinear);
        // Collinear points form a degenerate triangle (area = 0)
        // Algorithm may produce 0 or 1 degenerate triangle
        REQUIRE(dt.NumPoints() == 3);
    }
    
    SECTION("Four collinear points")
    {
        std::vector<Point2Cartesian> collinear = {
            Point2Cartesian(0, 0),
            Point2Cartesian(1, 0),
            Point2Cartesian(2, 0),
            Point2Cartesian(3, 0)
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(collinear);
        REQUIRE(dt.NumPoints() == 4);
    }
    
    SECTION("Mix of collinear and non-collinear")
    {
        std::vector<Point2Cartesian> points = {
            Point2Cartesian(0, 0),
            Point2Cartesian(1, 0),
            Point2Cartesian(2, 0),  // Collinear with first two
            Point2Cartesian(1, 1)   // Not collinear
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
        REQUIRE(dt.NumPoints() == 4);
        REQUIRE(dt.NumTriangles() >= 2);  // Should form valid triangulation
    }
}

TEST_CASE("DelaunayTriangulation - Duplicate points", "[ComputationalGeometry][Delaunay][EdgeCase]")
{
    SECTION("Exact duplicates")
    {
        std::vector<Point2Cartesian> points = {
            Point2Cartesian(0, 0),
            Point2Cartesian(0, 0),  // Duplicate
            Point2Cartesian(1, 0),
            Point2Cartesian(0.5, 1)
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
        // Should handle gracefully - may produce triangles with duplicate vertices
        REQUIRE(dt.NumPoints() == 4);
    }
    
    SECTION("Near-duplicate points")
    {
        std::vector<Point2Cartesian> points = {
            Point2Cartesian(0, 0),
            Point2Cartesian(1e-12, 1e-12),  // Very close
            Point2Cartesian(1, 0),
            Point2Cartesian(0.5, 1)
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
        REQUIRE(dt.NumPoints() == 4);
    }
}

TEST_CASE("DelaunayTriangulation - Co-circular points", "[ComputationalGeometry][Delaunay][EdgeCase]")
{
    // Points on a circle - degenerate case where multiple triangulations are valid
    std::vector<Point2Cartesian> points;
    int numPoints = 8;
    for (int i = 0; i < numPoints; i++) {
        Real angle = 2 * Constants::PI * i / numPoints;
        points.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == numPoints);
    REQUIRE(dt.NumTriangles() >= numPoints - 2);  // At least n-2 triangles
    
    // All points should be used
    std::set<int> usedVertices;
    for (const auto& tri : dt.triangles) {
        for (int v : tri) {
            usedVertices.insert(v);
        }
    }
    REQUIRE(usedVertices.size() == static_cast<size_t>(numPoints));
}

// ============================================================================
// DELAUNAY TRIANGULATION - NUMERICAL ROBUSTNESS
// ============================================================================

TEST_CASE("DelaunayTriangulation - Large coordinates", "[ComputationalGeometry][Delaunay][Numerical]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(1e6, 1e6),
        Point2Cartesian(1e6 + 1, 1e6),
        Point2Cartesian(1e6 + 0.5, 1e6 + 1)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 3);
    REQUIRE(dt.NumTriangles() == 1);
}

TEST_CASE("DelaunayTriangulation - Moderately small coordinates", "[ComputationalGeometry][Delaunay][Numerical]")
{
    // Note: The Bowyer-Watson algorithm involves determinant calculations with squared terms.
    // For very small coordinates (< 1e-3), numerical precision can become problematic.
    // For extremely small coordinates, users should normalize their data first.
    // This test uses coordinates that should work reliably with double precision.
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0.001, 0.001),
        Point2Cartesian(0.002, 0.001),
        Point2Cartesian(0.0015, 0.002)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 3);
    REQUIRE(dt.NumTriangles() == 1);
}

TEST_CASE("DelaunayTriangulation - Mixed scale coordinates", "[ComputationalGeometry][Delaunay][Numerical]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1000, 0),
        Point2Cartesian(500, 0.001),  // Very small Y relative to X span
        Point2Cartesian(500, 1000)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 4);
    REQUIRE(dt.NumTriangles() >= 2);
}

TEST_CASE("DelaunayTriangulation - Nearly collinear points", "[ComputationalGeometry][Delaunay][Numerical]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 1e-10),  // Nearly on the line y=0
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 1)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 4);
    REQUIRE(dt.NumTriangles() >= 2);
}

// ============================================================================
// DELAUNAY TRIANGULATION - PROPERTY VERIFICATION
// ============================================================================

TEST_CASE("DelaunayTriangulation - Triangle count formula", "[ComputationalGeometry][Delaunay][Property]")
{
    // For n points with h on convex hull: triangles = 2n - 2 - h
    SECTION("Square with center point: n=5, h=4, expect 4 triangles")
    {
        std::vector<Point2Cartesian> points = {
            Point2Cartesian(0, 0),
            Point2Cartesian(2, 0),
            Point2Cartesian(2, 2),
            Point2Cartesian(0, 2),
            Point2Cartesian(1, 1)  // Center
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
        // Expected: 2*5 - 2 - 4 = 4 triangles
        REQUIRE(dt.NumTriangles() == 4);
    }
    
    SECTION("Triangle with center: n=4, h=3, expect 3 triangles")
    {
        std::vector<Point2Cartesian> points = {
            Point2Cartesian(0, 0),
            Point2Cartesian(3, 0),
            Point2Cartesian(1.5, 3),
            Point2Cartesian(1.5, 1)  // Interior
        };
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
        // Expected: 2*4 - 2 - 3 = 3 triangles
        REQUIRE(dt.NumTriangles() == 3);
    }
    
    SECTION("Hexagon: n=6, h=6, expect 4 triangles")
    {
        std::vector<Point2Cartesian> points;
        for (int i = 0; i < 6; i++) {
            Real angle = 2 * Constants::PI * i / 6;
            points.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
        }
        auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
        // Expected: 2*6 - 2 - 6 = 4 triangles
        REQUIRE(dt.NumTriangles() == 4);
    }
}

TEST_CASE("DelaunayTriangulation - All points used", "[ComputationalGeometry][Delaunay][Property]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(4, 3),
        Point2Cartesian(0, 3),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(3, 2)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // All 7 points should appear in at least one triangle
    std::set<int> usedVertices;
    for (const auto& tri : dt.triangles) {
        for (int v : tri) {
            usedVertices.insert(v);
        }
    }
    REQUIRE(usedVertices.size() == 7);
}

TEST_CASE("DelaunayTriangulation - No overlapping triangles", "[ComputationalGeometry][Delaunay][Property]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(3, 3),
        Point2Cartesian(0, 3),
        Point2Cartesian(1.5, 1.5)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // Sum of triangle areas should equal convex hull area
    Real totalTriArea = 0;
    for (int i = 0; i < dt.NumTriangles(); i++) {
        Triangle2D tri = dt.GetTriangle(i);
        totalTriArea += std::abs(tri.Area());
    }
    
    // Convex hull is the 3x3 square = 9
    REQUIRE_THAT(totalTriArea, WithinAbs(9.0, 1e-6));
}

TEST_CASE("DelaunayTriangulation - Consistent orientation", "[ComputationalGeometry][Delaunay][Property]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2),
        Point2Cartesian(1, 1)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // All triangles should have consistent orientation (all positive or all negative area)
    // Check signed area
    int positiveCount = 0;
    int negativeCount = 0;
    
    for (int i = 0; i < dt.NumTriangles(); i++) {
        const auto& tri = dt.triangles[i];
        const Point2Cartesian& a = dt.points[tri[0]];
        const Point2Cartesian& b = dt.points[tri[1]];
        const Point2Cartesian& c = dt.points[tri[2]];
        
        Real signedArea = (b.X() - a.X()) * (c.Y() - a.Y()) - 
                          (c.X() - a.X()) * (b.Y() - a.Y());
        
        if (signedArea > 0) positiveCount++;
        else if (signedArea < 0) negativeCount++;
    }
    
    // All should be same orientation (allowing for numerical issues with zero-area)
    bool consistent = (positiveCount == 0 || negativeCount == 0);
    REQUIRE(consistent);
}

TEST_CASE("DelaunayTriangulation - Delaunay property exhaustive", "[ComputationalGeometry][Delaunay][Property]")
{
    // More thorough test of the empty circumcircle property
    std::mt19937 rng(12345);
    std::uniform_real_distribution<Real> dist(-10, 10);
    
    std::vector<Point2Cartesian> points;
    for (int i = 0; i < 30; i++) {
        points.push_back(Point2Cartesian(dist(rng), dist(rng)));
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // Check Delaunay property for ALL triangles
    int violations = 0;
    for (int t = 0; t < dt.NumTriangles(); t++) {
        const auto& tri = dt.triangles[t];
        const Point2Cartesian& a = dt.points[tri[0]];
        const Point2Cartesian& b = dt.points[tri[1]];
        const Point2Cartesian& c = dt.points[tri[2]];
        
        // Compute circumcircle
        Real ax = a.X(), ay = a.Y();
        Real bx = b.X(), by = b.Y();
        Real cx = c.X(), cy = c.Y();
        
        Real d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        if (std::abs(d) < 1e-10) continue;  // Degenerate
        
        Real ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + 
                   (cx*cx + cy*cy) * (ay - by)) / d;
        Real uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + 
                   (cx*cx + cy*cy) * (bx - ax)) / d;
        Real r2 = (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy);
        
        // Check ALL other points
        for (int i = 0; i < dt.NumPoints(); i++) {
            if (i == tri[0] || i == tri[1] || i == tri[2]) continue;
            
            const Point2Cartesian& p = dt.points[i];
            Real dist2 = (p.X() - ux) * (p.X() - ux) + (p.Y() - uy) * (p.Y() - uy);
            
            if (dist2 < r2 - 1e-6) {
                violations++;
            }
        }
    }
    
    REQUIRE(violations == 0);
}

// ============================================================================
// DELAUNAY TRIANGULATION - INTERPOLATION TESTS
// ============================================================================

TEST_CASE("DelaunayTriangulation - Linear function interpolation exact", "[ComputationalGeometry][Delaunay][Interpolation]")
{
    // For linear function f(x,y) = ax + by + c, Delaunay interpolation should be exact
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(3, 3),
        Point2Cartesian(0, 3),
        Point2Cartesian(1.5, 1.5)
    };
    
    // f(x,y) = 2x + 3y + 1
    auto f = [](Real x, Real y) { return 2*x + 3*y + 1; };
    
    std::vector<Real> values;
    for (const auto& p : points) {
        values.push_back(f(p.X(), p.Y()));
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // Test interpolation at several interior points
    std::vector<Point2Cartesian> testPoints = {
        Point2Cartesian(1, 1),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(2, 2),
        Point2Cartesian(0.5, 0.5)
    };
    
    for (const auto& p : testPoints) {
        Real interp = dt.Interpolate(p, values);
        Real expected = f(p.X(), p.Y());
        REQUIRE_THAT(interp, WithinAbs(expected, 1e-6));
    }
}

TEST_CASE("DelaunayTriangulation - Interpolation at vertices exact", "[ComputationalGeometry][Delaunay][Interpolation]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };
    
    std::vector<Real> values = {10.0, 20.0, 30.0};
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    // At each vertex, interpolation should return exact value
    for (size_t i = 0; i < points.size(); i++) {
        Real interp = dt.Interpolate(points[i], values);
        REQUIRE_THAT(interp, WithinAbs(values[i], 1e-6));
    }
}

TEST_CASE("DelaunayTriangulation - Interpolation outside returns default", "[ComputationalGeometry][Delaunay][Interpolation]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(0.5, 1)
    };
    
    std::vector<Real> values = {1.0, 2.0, 3.0};
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    Point2Cartesian outside(10, 10);
    Real defaultVal = -999.0;
    Real result = dt.Interpolate(outside, values, defaultVal);
    
    REQUIRE(result == defaultVal);
}

// ============================================================================
// DELAUNAY TRIANGULATION - SPECIAL CONFIGURATIONS
// ============================================================================

TEST_CASE("DelaunayTriangulation - L-shaped point set", "[ComputationalGeometry][Delaunay][Config]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(0, 2)
    };
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 6);
    REQUIRE(dt.NumTriangles() >= 4);  // L-shape should have at least 4 triangles
}

TEST_CASE("DelaunayTriangulation - Star pattern", "[ComputationalGeometry][Delaunay][Config]")
{
    std::vector<Point2Cartesian> points;
    
    // Center point
    points.push_back(Point2Cartesian(0, 0));
    
    // Outer points (star tips)
    for (int i = 0; i < 5; i++) {
        Real angle = 2 * Constants::PI * i / 5;
        points.push_back(Point2Cartesian(2 * std::cos(angle), 2 * std::sin(angle)));
    }
    
    // Inner points (between tips)
    for (int i = 0; i < 5; i++) {
        Real angle = 2 * Constants::PI * i / 5 + Constants::PI / 5;
        points.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 11);
    REQUIRE(dt.NumTriangles() >= 10);
}

TEST_CASE("DelaunayTriangulation - Clustered points", "[ComputationalGeometry][Delaunay][Config]")
{
    std::vector<Point2Cartesian> points;
    
    // Cluster 1 around (0,0)
    points.push_back(Point2Cartesian(0, 0));
    points.push_back(Point2Cartesian(0.1, 0));
    points.push_back(Point2Cartesian(0, 0.1));
    
    // Cluster 2 around (10,10)
    points.push_back(Point2Cartesian(10, 10));
    points.push_back(Point2Cartesian(10.1, 10));
    points.push_back(Point2Cartesian(10, 10.1));
    
    // Single point far away
    points.push_back(Point2Cartesian(5, 20));
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 7);
    REQUIRE(dt.NumTriangles() >= 5);
    
    // All points should be used
    std::set<int> used;
    for (const auto& tri : dt.triangles) {
        for (int v : tri) used.insert(v);
    }
    REQUIRE(used.size() == 7);
}

TEST_CASE("DelaunayTriangulation - Large point set", "[ComputationalGeometry][Delaunay][Performance]")
{
    std::mt19937 rng(54321);
    std::uniform_real_distribution<Real> dist(0, 100);
    
    std::vector<Point2Cartesian> points;
    for (int i = 0; i < 100; i++) {
        points.push_back(Point2Cartesian(dist(rng), dist(rng)));
    }
    
    auto dt = MML::CompGeometry::Triangulation::ComputeDelaunay(points);
    
    REQUIRE(dt.NumPoints() == 100);
    // For 100 points, expect roughly 180-190 triangles (2n - 2 - h where h ~ 10-20)
    REQUIRE(dt.NumTriangles() >= 150);
    REQUIRE(dt.NumTriangles() <= 200);
    
    // Verify Delaunay property on random sample
    std::uniform_int_distribution<int> triDist(0, dt.NumTriangles() - 1);
    for (int test = 0; test < 20; test++) {
        int t = triDist(rng);
        const auto& tri = dt.triangles[t];
        const Point2Cartesian& a = dt.points[tri[0]];
        const Point2Cartesian& b = dt.points[tri[1]];
        const Point2Cartesian& c = dt.points[tri[2]];
        
        Real ax = a.X(), ay = a.Y();
        Real bx = b.X(), by = b.Y();
        Real cx = c.X(), cy = c.Y();
        
        Real d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        if (std::abs(d) < 1e-10) continue;
        
        Real ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + 
                   (cx*cx + cy*cy) * (ay - by)) / d;
        Real uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + 
                   (cx*cx + cy*cy) * (bx - ax)) / d;
        Real r2 = (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy);
        
        for (int i = 0; i < dt.NumPoints(); i++) {
            if (i == tri[0] || i == tri[1] || i == tri[2]) continue;
            const Point2Cartesian& p = dt.points[i];
            Real dist2 = (p.X() - ux) * (p.X() - ux) + (p.Y() - uy) * (p.Y() - uy);
            REQUIRE(dist2 >= r2 - 1e-5);
        }
    }
}

} // namespace MML::Tests::Algorithms::CompGeometry::TriangulationTests


