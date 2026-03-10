/// @file convex_hull_tests.cpp
/// @brief Tests for 2D Convex Hull algorithms (Andrew's monotone chain)
/// @details Tests ConvexHull computation from point sets and polygons

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include "algorithms/ComputationalGeometry.h"

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::CompGeometry::ConvexHullTests
{

// ============================================================================
// CONVEX HULL TESTS
// ============================================================================

TEST_CASE("ConvexHull - Empty and degenerate cases", "[ComputationalGeometry][ConvexHull]")
{
    SECTION("Empty point set")
    {
        std::vector<Point2Cartesian> empty;
        Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(empty);
        REQUIRE(hull.NumVertices() == 0);
    }

    SECTION("Single point")
    {
        std::vector<Point2Cartesian> single = { Point2Cartesian(5, 3) };
        Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(single);
        REQUIRE(hull.NumVertices() == 1);
        REQUIRE_THAT(hull[0].X(), WithinAbs(REAL(5.0), REAL(1e-10)));
    }

    SECTION("Two points")
    {
        std::vector<Point2Cartesian> two = { 
            Point2Cartesian(0, 0), 
            Point2Cartesian(1, 1) 
        };
        Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(two);
        REQUIRE(hull.NumVertices() == 2);
    }

    SECTION("Collinear points")
    {
        std::vector<Point2Cartesian> collinear = {
            Point2Cartesian(0, 0),
            Point2Cartesian(1, 1),
            Point2Cartesian(2, 2),
            Point2Cartesian(3, 3)
        };
        Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(collinear);
        // Hull of collinear points should have 2 vertices (endpoints)
        REQUIRE(hull.NumVertices() == 2);
    }
}

TEST_CASE("ConvexHull - Square", "[ComputationalGeometry][ConvexHull]")
{
    // Points forming a square with an interior point
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1),
        Point2Cartesian(0.5, 0.5)  // Interior point
    };

    Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(points);
    
    REQUIRE(hull.NumVertices() == 4);
    REQUIRE(hull.IsConvex());
    REQUIRE_THAT(hull.Area(), WithinAbs(REAL(1.0), REAL(1e-10)));
}

TEST_CASE("ConvexHull - Triangle", "[ComputationalGeometry][ConvexHull]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };

    Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(points);
    
    REQUIRE(hull.NumVertices() == 3);
    REQUIRE(hull.IsConvex());
    REQUIRE_THAT(hull.Area(), WithinAbs(REAL(2.0), REAL(1e-10)));
}

TEST_CASE("ConvexHull - Random point cloud", "[ComputationalGeometry][ConvexHull]")
{
    // Create a Circle2D of points with some interior points
    std::vector<Point2Cartesian> points;
    
    // Add Circle2D boundary points
    for (int i = 0; i < 8; i++)
    {
        Real angle = i * Constants::PI / 4.0;
        points.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
    }
    
    // Add interior points
    points.push_back(Point2Cartesian(0, 0));
    points.push_back(Point2Cartesian(0.5, 0));
    points.push_back(Point2Cartesian(0, 0.5));

    Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(points);
    
    REQUIRE(hull.NumVertices() == 8);
    REQUIRE(hull.IsConvex());
    
    // All original boundary points should be on the hull
    // (interior points excluded)
}

TEST_CASE("ConvexHull - Duplicate points", "[ComputationalGeometry][ConvexHull]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(0, 0),  // duplicate
        Point2Cartesian(1, 0),
        Point2Cartesian(1, 0),  // duplicate
        Point2Cartesian(1, 1),
        Point2Cartesian(0, 1)
    };

    Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(points);
    
    REQUIRE(hull.NumVertices() == 4);
    REQUIRE(hull.IsConvex());
}

TEST_CASE("ConvexHull - From Polygon2D", "[ComputationalGeometry][ConvexHull]")
{
    // Non-convex polygon (L-shape)
    // Shape:    (0,2)---(1,2)
    //             |       |
    //             |   (1,1)---(2,1)
    //             |             |
    //           (0,0)---------(2,0)
    //
    // Convex hull: (0,0) -> (2,0) -> (2,1) -> (1,2) -> (0,2) = 5 vertices
    Polygon2D lShape({
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 1),
        Point2Cartesian(1, 1),
        Point2Cartesian(1, 2),
        Point2Cartesian(0, 2)
    });

    REQUIRE_FALSE(lShape.IsConvex());
    
    Polygon2D hull = MML::CompGeometry::ConvexHull2D::Compute(lShape);
    
    REQUIRE(hull.IsConvex());
    REQUIRE(hull.NumVertices() == 5);  // (0,0), (2,0), (2,1), (1,2), (0,2)
}

} // namespace MML::Tests::Algorithms::CompGeometry::ConvexHullTests

