/// @file intersections_tests.cpp
/// @brief Tests for line/segment/ray intersection algorithms
/// @details Tests SegmentIntersection, LineIntersection2D, LineIntersection3D, RayTriangleIntersection

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include "algorithms/ComputationalGeometry.h"

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::CompGeometry::IntersectionsTests {

// ============================================================================
// SEGMENT INTERSECTION TESTS
// ============================================================================

TEST_CASE("SegmentIntersection - Crossing segments", "[ComputationalGeometry][SegmentIntersection]")
{
    // X pattern
    Point2Cartesian p1(0, 0), q1(2, 2);  // diagonal /
    Point2Cartesian p2(0, 2), q2(2, 0);  // diagonal \

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::Point);
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("SegmentIntersection - Parallel non-intersecting", "[ComputationalGeometry][SegmentIntersection]")
{
    Point2Cartesian p1(0, 0), q1(2, 0);
    Point2Cartesian p2(0, 1), q2(2, 1);

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::None);
}

TEST_CASE("SegmentIntersection - T intersection (endpoint touch)", "[ComputationalGeometry][SegmentIntersection]")
{
    Point2Cartesian p1(0, 0), q1(2, 0);  // horizontal
    Point2Cartesian p2(1, 0), q2(1, 2);  // vertical, touching midpoint

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::Point);
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("SegmentIntersection - Collinear overlapping", "[ComputationalGeometry][SegmentIntersection]")
{
    Point2Cartesian p1(0, 0), q1(3, 0);
    Point2Cartesian p2(1, 0), q2(4, 0);  // overlaps from 1 to 3

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::Overlap);
    
    // Overlap should be from (1,0) to (3,0)
    Real minX = std::min(result.overlapStart.X(), result.overlapEnd.X());
    Real maxX = std::max(result.overlapStart.X(), result.overlapEnd.X());
    REQUIRE_THAT(minX, WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(maxX, WithinAbs(REAL(3.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("SegmentIntersection - Collinear non-overlapping", "[ComputationalGeometry][SegmentIntersection]")
{
    Point2Cartesian p1(0, 0), q1(1, 0);
    Point2Cartesian p2(2, 0), q2(3, 0);

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::Collinear);
}

TEST_CASE("SegmentIntersection - Endpoint to endpoint", "[ComputationalGeometry][SegmentIntersection]")
{
    Point2Cartesian p1(0, 0), q1(1, 1);
    Point2Cartesian p2(1, 1), q2(2, 0);

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::Point);
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("SegmentIntersection - No intersection (skew)", "[ComputationalGeometry][SegmentIntersection]")
{
    Point2Cartesian p1(0, 0), q1(1, 1);
    Point2Cartesian p2(2, 0), q2(3, 1);

    auto result = MML::CompGeometry::Intersections::IntersectSegments(p1, q1, p2, q2);
    
    REQUIRE(result.type == SegmentIntersectionType::None);
}

TEST_CASE("SegmentIntersection - Using SegmentLine2D", "[ComputationalGeometry][SegmentIntersection]")
{
    SegmentLine2D seg1(Point2Cartesian(0, 0), Point2Cartesian(2, 2));
    SegmentLine2D seg2(Point2Cartesian(0, 2), Point2Cartesian(2, 0));

    auto result = MML::CompGeometry::Intersections::IntersectSegments(seg1, seg2);
    
    REQUIRE(result.IsPointIntersection());
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

// ============================================================================
// LINE-LINE INTERSECTION 2D TESTS
// ============================================================================

TEST_CASE("LineIntersection2D - Crossing lines", "[ComputationalGeometry][LineIntersection]")
{
    // Horizontal line y = 1
    Point2Cartesian p1(0, 1);
    Vector2Cartesian d1(1, 0);
    
    // Vertical line x = 2
    Point2Cartesian p2(2, 0);
    Vector2Cartesian d2(0, 1);

    auto result = MML::CompGeometry::Intersections::IntersectLines2D(p1, d1, p2, d2);
    
    REQUIRE(result.IsPoint());
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("LineIntersection2D - Parallel lines", "[ComputationalGeometry][LineIntersection]")
{
    Point2Cartesian p1(0, 0);
    Vector2Cartesian d1(1, 1);
    
    Point2Cartesian p2(0, 1);  // Parallel line, offset by 1 in y
    Vector2Cartesian d2(1, 1);

    auto result = MML::CompGeometry::Intersections::IntersectLines2D(p1, d1, p2, d2);
    
    REQUIRE(result.IsParallel());
}

TEST_CASE("LineIntersection2D - Coincident lines", "[ComputationalGeometry][LineIntersection]")
{
    Point2Cartesian p1(0, 0);
    Vector2Cartesian d1(1, 1);
    
    Point2Cartesian p2(2, 2);  // Same line, different point
    Vector2Cartesian d2(2, 2);  // Same direction (scaled)

    auto result = MML::CompGeometry::Intersections::IntersectLines2D(p1, d1, p2, d2);
    
    REQUIRE(result.IsCoincident());
}

TEST_CASE("LineIntersection2D - Diagonal crossing", "[ComputationalGeometry][LineIntersection]")
{
    // Line y = x
    Point2Cartesian p1(0, 0);
    Vector2Cartesian d1(1, 1);
    
    // Line y = -x + 4
    Point2Cartesian p2(0, 4);
    Vector2Cartesian d2(1, -1);

    auto result = MML::CompGeometry::Intersections::IntersectLines2D(p1, d1, p2, d2);
    
    REQUIRE(result.IsPoint());
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("LineIntersection2D - Using Line2D objects", "[ComputationalGeometry][LineIntersection]")
{
    Line2D line1(Point2Cartesian(0, 0), Vector2Cartesian(1, 0));  // y = 0
    Line2D line2(Point2Cartesian(0, 0), Vector2Cartesian(0, 1));  // x = 0

    auto result = MML::CompGeometry::Intersections::IntersectLines2D(line1, line2);
    
    REQUIRE(result.IsPoint());
    REQUIRE_THAT(result.point.X(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point.Y(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
}

// ============================================================================
// LINE-LINE INTERSECTION 3D TESTS
// ============================================================================

TEST_CASE("LineIntersection3D - Intersecting lines", "[ComputationalGeometry][LineIntersection3D]")
{
    // Line along X axis
    Point3Cartesian p1(0, 0, 0);
    Vector3Cartesian d1(1, 0, 0);
    
    // Line along Y axis
    Point3Cartesian p2(0, 0, 0);
    Vector3Cartesian d2(0, 1, 0);

    auto result = MML::CompGeometry::Intersections::IntersectLines3D(p1, d1, p2, d2);
    
    REQUIRE(result.IsIntersecting());
    REQUIRE_THAT(result.point1.X(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point1.Y(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
    REQUIRE_THAT(result.point1.Z(), WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("LineIntersection3D - Parallel lines", "[ComputationalGeometry][LineIntersection3D]")
{
    Point3Cartesian p1(0, 0, 0);
    Vector3Cartesian d1(1, 0, 0);
    
    Point3Cartesian p2(0, 1, 0);  // Parallel, offset in Y
    Vector3Cartesian d2(1, 0, 0);

    auto result = MML::CompGeometry::Intersections::IntersectLines3D(p1, d1, p2, d2);
    
    REQUIRE(result.IsParallel());
}

TEST_CASE("LineIntersection3D - Skew lines", "[ComputationalGeometry][LineIntersection3D]")
{
    // Line along X axis at z=0
    Point3Cartesian p1(0, 0, 0);
    Vector3Cartesian d1(1, 0, 0);
    
    // Line along Y axis at z=1 (skew - don't intersect)
    Point3Cartesian p2(0, 0, 1);
    Vector3Cartesian d2(0, 1, 0);

    auto result = MML::CompGeometry::Intersections::IntersectLines3D(p1, d1, p2, d2);
    
    REQUIRE(result.IsSkew());
    REQUIRE_THAT(result.distance, WithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
}

TEST_CASE("LineIntersection3D - Coincident lines", "[ComputationalGeometry][LineIntersection3D]")
{
    Point3Cartesian p1(0, 0, 0);
    Vector3Cartesian d1(1, 1, 1);
    
    Point3Cartesian p2(2, 2, 2);  // Same line
    Vector3Cartesian d2(1, 1, 1);

    auto result = MML::CompGeometry::Intersections::IntersectLines3D(p1, d1, p2, d2);
    
    REQUIRE(result.IsCoincident());
}

TEST_CASE("LineIntersection3D - Using Line3D objects", "[ComputationalGeometry][LineIntersection3D]")
{
    Line3D line1(Point3Cartesian(0, 0, 0), Vector3Cartesian(1, 0, 0));
    Line3D line2(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 0, 1));

    auto result = MML::CompGeometry::Intersections::IntersectLines3D(line1, line2);
    
    REQUIRE(result.IsIntersecting());
}

// ============================================================================
// RAY-TRIANGLE INTERSECTION TESTS (Moller-Trumbore Algorithm)
// ============================================================================

TEST_CASE("RayTriangleIntersection - Direct hit through center", "[ComputationalGeometry][RayTriangle]")
{
    // Triangle in XY plane at z=0
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(2, 0, 0);
    Point3Cartesian v2(1, 2, 0);
    
    // Ray from above, pointing down through centroid
    Point3Cartesian origin(1, 0.666667, 5);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == true);
    REQUIRE_THAT(hit.t, WithinAbs(5.0, 1e-6));
    REQUIRE_THAT(hit.point.Z(), WithinAbs(0.0, 1e-6));
    
    // Barycentric coords should sum to ~1
    REQUIRE_THAT(hit.u + hit.v + hit.BarycentricW(), WithinAbs(1.0, TOL(1e-10, 1e-5)));
}

TEST_CASE("RayTriangleIntersection - Miss (parallel ray)", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 0);
    Point3Cartesian v2(0, 1, 0);
    
    // Ray parallel to triangle plane
    Point3Cartesian origin(0.5, 0.5, 1);
    Vector3Cartesian dir(1, 0, 0);  // Parallel to XY plane
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == false);
}

TEST_CASE("RayTriangleIntersection - Miss (outside triangle)", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 0);
    Point3Cartesian v2(0, 1, 0);
    
    // Ray pointing at plane but outside triangle
    Point3Cartesian origin(2, 2, 1);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == false);
}

TEST_CASE("RayTriangleIntersection - Hit at vertex", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 0);
    Point3Cartesian v2(0, 1, 0);
    
    // Ray hitting exactly at v0
    Point3Cartesian origin(0, 0, 1);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == true);
    REQUIRE_THAT(hit.t, WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(hit.point.X(), WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(hit.point.Y(), WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(hit.point.Z(), WithinAbs(0.0, 1e-6));
}

TEST_CASE("RayTriangleIntersection - Hit on edge", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(2, 0, 0);
    Point3Cartesian v2(0, 2, 0);
    
    // Ray hitting the v0-v1 edge at midpoint
    Point3Cartesian origin(1, 0, 5);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == true);
    REQUIRE_THAT(hit.point.X(), WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(hit.point.Y(), WithinAbs(0.0, 1e-6));
    // v = 0 on edge v0-v1
    REQUIRE_THAT(hit.v, WithinAbs(0.0, 1e-6));
}

TEST_CASE("RayTriangleIntersection - Triangle behind ray", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 0);
    Point3Cartesian v2(0, 1, 0);
    
    // Ray starting at z=1, pointing UP (away from triangle at z=0)
    Point3Cartesian origin(0.25, 0.25, 1);
    Vector3Cartesian dir(0, 0, 1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == false);  // Triangle is behind ray origin
}

TEST_CASE("RayTriangleIntersection - Backface culling", "[ComputationalGeometry][RayTriangle]")
{
    // Triangle with normal pointing in +Z direction (CCW when viewed from +Z)
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 0);
    Point3Cartesian v2(0, 1, 0);
    
    // Ray from below pointing up (hitting backface)
    Point3Cartesian origin(0.25, 0.25, -1);
    Vector3Cartesian dir(0, 0, 1);
    
    // Without culling - should hit
    auto hit_no_cull = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2, false);
    REQUIRE(hit_no_cull.hit == true);
    
    // With culling - should miss (backface)
    auto hit_cull = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2, true);
    REQUIRE(hit_cull.hit == false);
}

TEST_CASE("RayTriangleIntersection - Non-unit direction", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 0);
    Point3Cartesian v2(0, 1, 0);
    
    // Ray with non-unit direction vector
    Point3Cartesian origin(0.25, 0.25, 10);
    Vector3Cartesian dir(0, 0, -5);  // Length 5, not normalized
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == true);
    // t should be 2 (10 units / 5 per unit)
    REQUIRE_THAT(hit.t, WithinAbs(2.0, 1e-6));
    REQUIRE_THAT(hit.point.Z(), WithinAbs(0.0, 1e-6));
}

TEST_CASE("RayTriangleIntersection - Using Triangle3D object", "[ComputationalGeometry][RayTriangle]")
{
    Triangle3D triangle(
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(2, 0, 0),
        Point3Cartesian(1, 2, 0)
    );
    
    Point3Cartesian origin(1, 0.666667, 3);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, triangle);
    
    REQUIRE(hit.hit == true);
    REQUIRE_THAT(hit.t, WithinAbs(3.0, 1e-6));
}

TEST_CASE("RayTriangleIntersection - Using Line3D as ray", "[ComputationalGeometry][RayTriangle]")
{
    Triangle3D triangle(
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(1, 0, 0),
        Point3Cartesian(0, 1, 0)
    );
    
    Line3D ray(Point3Cartesian(0.25, 0.25, 2), Vector3Cartesian(0, 0, -1));
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(ray, triangle);
    
    REQUIRE(hit.hit == true);
    REQUIRE_THAT(hit.t, WithinAbs(2.0, 1e-6));
}

TEST_CASE("RayTriangleIntersection - Tilted triangle", "[ComputationalGeometry][RayTriangle]")
{
    // Triangle tilted in 3D space
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(1, 0, 1);  // Tilted up
    Point3Cartesian v2(0, 1, 0);
    
    // Ray from above
    Point3Cartesian origin(0.25, 0.25, 5);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == true);
    // Should hit somewhere between z=0 and z=1
    REQUIRE(hit.point.Z() >= 0.0);
    REQUIRE(hit.point.Z() <= 1.0);
}

TEST_CASE("RayTriangleIntersection - Barycentric coordinates verify", "[ComputationalGeometry][RayTriangle]")
{
    Point3Cartesian v0(0, 0, 0);
    Point3Cartesian v1(3, 0, 0);
    Point3Cartesian v2(0, 3, 0);
    
    // Hit point at (1, 1, 0)
    Point3Cartesian origin(1, 1, 1);
    Vector3Cartesian dir(0, 0, -1);
    
    auto hit = MML::CompGeometry::Intersections::IntersectRayTriangle(origin, dir, v0, v1, v2);
    
    REQUIRE(hit.hit == true);
    
    // Verify: hitPoint = w*v0 + u*v1 + v*v2
    Real w = hit.BarycentricW();
    Point3Cartesian reconstructed = w * v0 + hit.u * v1 + hit.v * v2;
    
    REQUIRE_THAT(reconstructed.X(), WithinAbs(hit.point.X(), 1e-6));
    REQUIRE_THAT(reconstructed.Y(), WithinAbs(hit.point.Y(), 1e-6));
    REQUIRE_THAT(reconstructed.Z(), WithinAbs(hit.point.Z(), 1e-6));
}

} // namespace MML::Tests::Algorithms::CompGeometry::IntersectionsTests



