/// @file circles_tests.cpp
/// @brief Tests for circle-related algorithms
/// @details Tests MinimumEnclosingCircle (Welzl) and ClosestPair (divide-and-conquer)

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include "algorithms/ComputationalGeometry.h"

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::CompGeometry::CirclesTests
{

// ============================================================================
// MINIMUM ENCLOSING CIRCLE TESTS
// ============================================================================

TEST_CASE("MinimumEnclosingCircle - Empty and single point", "[ComputationalGeometry][MEC]")
{
    SECTION("Empty point set")
    {
        std::vector<Point2Cartesian> empty;
        Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(empty);
        REQUIRE_THAT(mec.Radius(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
    }

    SECTION("Single point")
    {
        std::vector<Point2Cartesian> single = { Point2Cartesian(3, 4) };
        Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(single);
        REQUIRE_THAT(mec.Radius(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(mec.Center().X(), WithinAbs(REAL(3.0), TOL(1e-10, 1e-5)));
        REQUIRE_THAT(mec.Center().Y(), WithinAbs(REAL(4.0), TOL(1e-10, 1e-5)));
    }
}

TEST_CASE("MinimumEnclosingCircle - Two points", "[ComputationalGeometry][MEC]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0)
    };

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(points);
    
    // Circle2D should have diameter = distance between points
    REQUIRE_THAT(mec.Radius(), WithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(mec.Center().X(), WithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(mec.Center().Y(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("MinimumEnclosingCircle - Three points (non-collinear)", "[ComputationalGeometry][MEC]")
{
    // Right triangle: circumcircle has hypotenuse as diameter
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(0, 3)
    };

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(points);
    
    // Hypotenuse length = 5, so radius = 2.5
    REQUIRE_THAT(mec.Radius(), WithinAbs(REAL(2.5), TOL(1e-10, 1e-5)));
    
    // Center should be midpoint of hypotenuse
    REQUIRE_THAT(mec.Center().X(), WithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(mec.Center().Y(), WithinAbs(REAL(1.5), TOL(1e-10, 1e-5)));
}

TEST_CASE("MinimumEnclosingCircle - Square with interior points", "[ComputationalGeometry][MEC]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2),
        Point2Cartesian(1, 1),  // center
        Point2Cartesian(0.5, 0.5),  // interior
        Point2Cartesian(1.5, 1.5)   // interior
    };

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(points);
    
    // Circle2D should pass through diagonal corners
    // Diagonal = 2√2, radius = √2
    REQUIRE_THAT(mec.Radius(), WithinAbs(std::sqrt(REAL(2.0)), TOL(1e-8, 1e-4)));
    REQUIRE_THAT(mec.Center().X(), WithinAbs(REAL(1.0), TOL(1e-8, 1e-4)));
    REQUIRE_THAT(mec.Center().Y(), WithinAbs(REAL(1.0), TOL(1e-8, 1e-4)));
}

TEST_CASE("MinimumEnclosingCircle - Points on Circle2D boundary", "[ComputationalGeometry][MEC]")
{
    // Unit Circle2D points
    std::vector<Point2Cartesian> points;
    for (int i = 0; i < 8; i++)
    {
        Real angle = i * Constants::PI / 4.0;
        points.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
    }

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(points);
    
    REQUIRE_THAT(mec.Radius(), WithinAbs(REAL(1.0), TOL(1e-8, 1e-4)));
    REQUIRE_THAT(mec.Center().X(), WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
    REQUIRE_THAT(mec.Center().Y(), WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
}

TEST_CASE("MinimumEnclosingCircle - Collinear points", "[ComputationalGeometry][MEC]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(4, 0)
    };

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(points);
    
    // Should enclose all points
    REQUIRE_THAT(mec.Radius(), WithinAbs(REAL(2.0), TOL(1e-8, 1e-4)));
    REQUIRE_THAT(mec.Center().X(), WithinAbs(REAL(2.0), TOL(1e-8, 1e-4)));
    
    // Verify all points are inside
    for (const auto& p : points)
    {
        REQUIRE(mec.Contains(p));
    }
}

TEST_CASE("MinimumEnclosingCircle - From Polygon2D", "[ComputationalGeometry][MEC]")
{
    Polygon2D triangle({
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(2, 3)
    });

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(triangle);
    
    // All vertices should be inside or on the Circle2D
    for (int i = 0; i < triangle.NumVertices(); i++)
    {
        REQUIRE(mec.Contains(triangle[i]));
    }
}

TEST_CASE("MinimumEnclosingCircle - All points inside result", "[ComputationalGeometry][MEC]")
{
    // Random-ish point cloud
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(1, 2),
        Point2Cartesian(3, 1),
        Point2Cartesian(2, 4),
        Point2Cartesian(4, 3),
        Point2Cartesian(0, 3),
        Point2Cartesian(2.5, 2.5)
    };

    Circle2D mec = MML::CompGeometry::Circles::MinimumEnclosingCircle(points);
    
    // All points must be inside or on the Circle2D
    for (const auto& p : points)
    {
        Real dist = p.Dist(mec.Center());
        REQUIRE(dist <= mec.Radius() + TOL(1e-8, 1e-4));
    }
}

// ============================================================================
// CLOSEST PAIR OF POINTS TESTS
// ============================================================================

TEST_CASE("ClosestPair - Two points", "[ComputationalGeometry][ClosestPair]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 4)
    };

    auto result = MML::CompGeometry::Circles::ClosestPair(points);
    
    REQUIRE_THAT(result.distance, WithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ClosestPair - Three points", "[ComputationalGeometry][ClosestPair]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),  // Closest to (0,0)
        Point2Cartesian(5, 5)
    };

    auto result = MML::CompGeometry::Circles::ClosestPair(points);
    
    REQUIRE_THAT(result.distance, WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ClosestPair - Collinear points", "[ComputationalGeometry][ClosestPair]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(5, 0),
        Point2Cartesian(6, 0)  // Closest pair: (5,0) and (6,0)
    };

    auto result = MML::CompGeometry::Circles::ClosestPair(points);
    
    REQUIRE_THAT(result.distance, WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ClosestPair - Grid of points", "[ComputationalGeometry][ClosestPair]")
{
    // 3x3 grid with spacing 1
    std::vector<Point2Cartesian> points;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            points.push_back(Point2Cartesian(i, j));
        }
    }

    auto result = MML::CompGeometry::Circles::ClosestPair(points);
    
    // Closest pair should have distance 1 (adjacent points)
    REQUIRE_THAT(result.distance, WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("ClosestPair - Random cloud", "[ComputationalGeometry][ClosestPair]")
{
    std::vector<Point2Cartesian> points = {
        Point2Cartesian(1, 2),
        Point2Cartesian(5, 8),
        Point2Cartesian(3, 1),
        Point2Cartesian(7, 4),
        Point2Cartesian(3.1, 1.1),  // Very close to (3,1)
        Point2Cartesian(9, 2)
    };

    auto result = MML::CompGeometry::Circles::ClosestPair(points);
    
    // (3,1) and (3.1,1.1) are closest: distance = sqrt(0.01+0.01) = sqrt(0.02)
    Real expected = std::sqrt(REAL(0.02));
    REQUIRE_THAT(result.distance, WithinAbs(expected, TOL(1e-8, 1e-4)));
}

} // namespace MML::Tests::Algorithms::CompGeometry::CirclesTests

