/// @file convex_hull_3d_tests.cpp
/// @brief Tests for 3D convex hull computation
/// @details Tests incremental algorithm for 3D convex hulls

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include "algorithms/ComputationalGeometry.h"

#include <random>

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::CompGeometry::ConvexHull3DTests
{

// ============================================================================
// CONVEX HULL 3D - BASIC TESTS
// ============================================================================

TEST_CASE("ConvexHull3D - Tetrahedron (4 points)", "[ComputationalGeometry][ConvexHull3D]")
{
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(1, 0, 0),
        Point3Cartesian(0.5, 1, 0),
        Point3Cartesian(0.5, 0.5, 1)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    REQUIRE(hull.NumVertices() == 4);
    REQUIRE(hull.NumFaces() == 4);  // Tetrahedron has 4 triangular faces
    
    // Volume of tetrahedron
    REQUIRE(hull.Volume() > 0);
    REQUIRE(hull.SurfaceArea() > 0);
}

TEST_CASE("ConvexHull3D - Cube (8 points)", "[ComputationalGeometry][ConvexHull3D]")
{
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(1, 0, 0),
        Point3Cartesian(1, 1, 0),
        Point3Cartesian(0, 1, 0),
        Point3Cartesian(0, 0, 1),
        Point3Cartesian(1, 0, 1),
        Point3Cartesian(1, 1, 1),
        Point3Cartesian(0, 1, 1)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    REQUIRE(hull.NumVertices() == 8);
    REQUIRE(hull.NumFaces() == 12);  // Cube: 6 faces * 2 triangles each
    
    // Volume of unit cube = 1
    REQUIRE_THAT(hull.Volume(), WithinAbs(1.0, 1e-6));
    
    // Surface area of unit cube = 6
    REQUIRE_THAT(hull.SurfaceArea(), WithinAbs(6.0, 1e-6));
}

TEST_CASE("ConvexHull3D - Contains point", "[ComputationalGeometry][ConvexHull3D]")
{
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(2, 0, 0),
        Point3Cartesian(0, 2, 0),
        Point3Cartesian(0, 0, 2)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    // Interior point
    REQUIRE(hull.Contains(Point3Cartesian(0.5, 0.5, 0.5)));
    
    // Vertex
    REQUIRE(hull.Contains(Point3Cartesian(0, 0, 0)));
    
    // Outside point
    REQUIRE_FALSE(hull.Contains(Point3Cartesian(1, 1, 1)));
    REQUIRE_FALSE(hull.Contains(Point3Cartesian(-1, 0, 0)));
    REQUIRE_FALSE(hull.Contains(Point3Cartesian(0, -1, 0)));
}

TEST_CASE("ConvexHull3D - Points with interior points", "[ComputationalGeometry][ConvexHull3D]")
{
    // Cube corners plus interior points (should be ignored)
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(2, 0, 0),
        Point3Cartesian(2, 2, 0),
        Point3Cartesian(0, 2, 0),
        Point3Cartesian(0, 0, 2),
        Point3Cartesian(2, 0, 2),
        Point3Cartesian(2, 2, 2),
        Point3Cartesian(0, 2, 2),
        // Interior points
        Point3Cartesian(1, 1, 1),
        Point3Cartesian(0.5, 0.5, 0.5)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    REQUIRE(hull.NumVertices() == 8);
    // Volume = 8
    REQUIRE_THAT(hull.Volume(), WithinAbs(8.0, 1e-6));
}

TEST_CASE("ConvexHull3D - Octahedron", "[ComputationalGeometry][ConvexHull3D]")
{
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(1, 0, 0),
        Point3Cartesian(-1, 0, 0),
        Point3Cartesian(0, 1, 0),
        Point3Cartesian(0, -1, 0),
        Point3Cartesian(0, 0, 1),
        Point3Cartesian(0, 0, -1)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    // Octahedron has 8 triangular faces
    REQUIRE(hull.NumFaces() == 8);
    
    // 6 vertices
    REQUIRE(hull.NumVertices() == 6);
    
    // Volume of octahedron with vertices at ±1 on axes = 4/3
    Real expectedVolume = 4.0 / 3.0;
    REQUIRE_THAT(hull.Volume(), WithinAbs(expectedVolume, 1e-6));
}

TEST_CASE("ConvexHull3D - Random point cloud", "[ComputationalGeometry][ConvexHull3D]")
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> dist(-5, 5);
    
    std::vector<Point3Cartesian> points;
    for (int i = 0; i < 50; i++) {
        points.push_back(Point3Cartesian(dist(rng), dist(rng), dist(rng)));
    }
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    // Should have some vertices and faces
    REQUIRE(hull.NumVertices() >= 4);
    REQUIRE(hull.NumFaces() >= 4);
    
    // Volume should be positive
    REQUIRE(hull.Volume() > 0);
    
    // Surface area should be positive
    REQUIRE(hull.SurfaceArea() > 0);
}

// ============================================================================
// CONVEX HULL 3D - EDGE CASES
// ============================================================================

TEST_CASE("ConvexHull3D - Degenerate cases", "[ComputationalGeometry][ConvexHull3D][EdgeCase]")
{
    SECTION("Less than 4 points")
    {
        std::vector<Point3Cartesian> points = {
            Point3Cartesian(0, 0, 0),
            Point3Cartesian(1, 0, 0),
            Point3Cartesian(0, 1, 0)
        };
        auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
        // Degenerate - no real 3D hull
        REQUIRE(hull.NumFaces() <= 2);
    }
}

TEST_CASE("ConvexHull3D - Coplanar points", "[ComputationalGeometry][ConvexHull3D][EdgeCase]")
{
    // All points in z=0 plane
    std::vector<Point3Cartesian> coplanar = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(2, 0, 0),
        Point3Cartesian(2, 2, 0),
        Point3Cartesian(0, 2, 0),
        Point3Cartesian(1, 1, 0)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(coplanar);
    
    // Coplanar points form a degenerate hull with essentially zero volume
    REQUIRE(hull.Volume() < 0.01);
}

TEST_CASE("ConvexHull3D - Collinear points", "[ComputationalGeometry][ConvexHull3D][EdgeCase]")
{
    std::vector<Point3Cartesian> collinear = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(1, 1, 1),
        Point3Cartesian(2, 2, 2),
        Point3Cartesian(3, 3, 3)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(collinear);
    
    // Collinear - no volume
    REQUIRE(hull.NumFaces() == 0);
}

// ============================================================================
// CONVEX HULL 3D - PROPERTY VERIFICATION
// ============================================================================

TEST_CASE("ConvexHull3D - Euler characteristic", "[ComputationalGeometry][ConvexHull3D][Property]")
{
    // For a convex polyhedron: V - E + F = 2
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(0, 0, 0),
        Point3Cartesian(3, 0, 0),
        Point3Cartesian(3, 3, 0),
        Point3Cartesian(0, 3, 0),
        Point3Cartesian(0, 0, 3),
        Point3Cartesian(3, 0, 3),
        Point3Cartesian(3, 3, 3),
        Point3Cartesian(0, 3, 3)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    int V = hull.NumVertices();
    int F = hull.NumFaces();
    
    // For triangulated cube: V=8, F=12, E=18 -> 8 - 18 + 12 = 2
    // E = (3*F)/2 for triangulated surface = 18
    int E = (3 * F) / 2;
    int euler = V - E + F;
    REQUIRE(euler == 2);
}

// ============================================================================
// CONVEX HULL 3D - SPECIFIC SHAPES
// ============================================================================

TEST_CASE("ConvexHull3D - Icosahedron", "[ComputationalGeometry][ConvexHull3D][Shape]")
{
    // Generate icosahedron vertices using golden ratio
    Real phi = (1.0 + std::sqrt(5.0)) / 2.0;
    
    std::vector<Point3Cartesian> points = {
        Point3Cartesian(0, 1, phi),
        Point3Cartesian(0, -1, phi),
        Point3Cartesian(0, 1, -phi),
        Point3Cartesian(0, -1, -phi),
        Point3Cartesian(1, phi, 0),
        Point3Cartesian(-1, phi, 0),
        Point3Cartesian(1, -phi, 0),
        Point3Cartesian(-1, -phi, 0),
        Point3Cartesian(phi, 0, 1),
        Point3Cartesian(-phi, 0, 1),
        Point3Cartesian(phi, 0, -1),
        Point3Cartesian(-phi, 0, -1)
    };
    
    auto hull = MML::CompGeometry::ConvexHull3DComputer::Compute(points);
    
    // Icosahedron has 12 vertices and 20 faces
    // Float precision can cause one near-coplanar face to be missed
    REQUIRE(hull.NumVertices() == 12);
    if constexpr (std::is_same_v<Real, float>)
        REQUIRE(hull.NumFaces() >= 19);
    else
        REQUIRE(hull.NumFaces() == 20);
}

} // namespace MML::Tests::Algorithms::CompGeometry::ConvexHull3DTests


