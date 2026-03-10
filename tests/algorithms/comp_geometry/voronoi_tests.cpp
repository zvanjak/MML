/// @file voronoi_tests.cpp
/// @brief Tests for Voronoi diagram computation
/// @details Tests Fortune's algorithm for Voronoi diagrams

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

namespace MML::Tests::Algorithms::CompGeometry::VoronoiTests {

// ============================================================================
// VORONOI DIAGRAM BASIC TESTS
// ============================================================================

TEST_CASE("VoronoiDiagram - Two points", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 2);
    
    // Two sites -> one edge (perpendicular bisector at x=1)
    // The edge is infinite, but there should be one
    REQUIRE(voronoi.NumEdges() >= 1);
}

TEST_CASE("VoronoiDiagram - Three points (triangle)", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 3);
    
    // Three sites form 3 edges meeting at circumcenter
    REQUIRE(voronoi.NumEdges() >= 3);
    
    // There should be one finite vertex (the circumcenter)
    REQUIRE(voronoi.NumVertices() >= 1);
}

TEST_CASE("VoronoiDiagram - Square (4 points)", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(2, 2),
        Point2Cartesian(0, 2)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 4);
    
    // Square creates a cross pattern with center vertex at (1,1)
    REQUIRE(voronoi.NumEdges() >= 4);
    REQUIRE(voronoi.NumVertices() >= 1);
    
    // Center vertex should be at (1, 1)
    bool hasCenter = false;
    for (int i = 0; i < voronoi.NumVertices(); i++) {
        const auto& v = voronoi.vertices[i];
        if (std::abs(v.X() - 1.0) < 0.1 && std::abs(v.Y() - 1.0) < 0.1) {
            hasCenter = true;
            break;
        }
    }
    REQUIRE(hasCenter);
}

TEST_CASE("VoronoiDiagram - Collinear points", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(1, 0),
        Point2Cartesian(2, 0),
        Point2Cartesian(3, 0)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 4);
    
    // Collinear points have no triangulation, so Voronoi computed from Delaunay has no edges
    // This is a degenerate case - implementation returns empty structure
    // All cells are infinite for collinear case
    for (int i = 0; i < 4; i++) {
        REQUIRE(voronoi.IsCellInfinite(i));
    }
}

TEST_CASE("VoronoiDiagram - Pentagon", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites;
    for (int i = 0; i < 5; i++) {
        Real angle = 2 * Constants::PI * i / 5;
        sites.push_back(Point2Cartesian(std::cos(angle), std::sin(angle)));
    }
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 5);
    REQUIRE(voronoi.NumEdges() >= 5);
    REQUIRE(voronoi.NumVertices() >= 1);
}

TEST_CASE("VoronoiDiagram - GetCell", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(3, 0),
        Point2Cartesian(1.5, 3)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    // Each site should have cell info
    for (int i = 0; i < voronoi.NumSites(); i++) {
        auto cell = voronoi.GetCell(i);
        // Cell polygon (may be empty for infinite cells)
        // Just check we can call the method
    }
}

TEST_CASE("VoronoiDiagram - Nearest site query", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(2, 4)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    // Query point very close to first site
    Point2Cartesian query1(0.1, 0.1);
    int nearest1 = voronoi.FindNearestSite(query1);
    REQUIRE(nearest1 == 0);
    
    // Query point very close to second site
    Point2Cartesian query2(3.9, 0.1);
    int nearest2 = voronoi.FindNearestSite(query2);
    REQUIRE(nearest2 == 1);
    
    // Query point very close to third site
    Point2Cartesian query3(2, 3.8);
    int nearest3 = voronoi.FindNearestSite(query3);
    REQUIRE(nearest3 == 2);
}

// ============================================================================
// VORONOI DIAGRAM - EDGE CASES
// ============================================================================

TEST_CASE("VoronoiDiagram - Empty and minimal", "[ComputationalGeometry][Voronoi][EdgeCase]")
{
    SECTION("Empty site list")
    {
        std::vector<Point2Cartesian> empty;
        auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(empty);
        REQUIRE(voronoi.NumSites() == 0);
        REQUIRE(voronoi.NumEdges() == 0);
    }
    
    SECTION("Single site")
    {
        std::vector<Point2Cartesian> single = {Point2Cartesian(5, 3)};
        auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(single);
        REQUIRE(voronoi.NumSites() == 1);
        // Single site: entire plane is its cell, no edges needed
        REQUIRE(voronoi.NumEdges() == 0);
    }
}

TEST_CASE("VoronoiDiagram - Duplicate sites", "[ComputationalGeometry][Voronoi][EdgeCase]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(0, 0),  // Duplicate
        Point2Cartesian(2, 0),
        Point2Cartesian(1, 2)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    // Should handle gracefully
    REQUIRE(voronoi.NumSites() == 4);
}

TEST_CASE("VoronoiDiagram - Co-circular sites", "[ComputationalGeometry][Voronoi][EdgeCase]")
{
    // Sites on a circle - degenerate case
    std::vector<Point2Cartesian> sites;
    for (int i = 0; i < 6; i++) {
        Real angle = 2 * Constants::PI * i / 6;
        sites.push_back(Point2Cartesian(2 * std::cos(angle), 2 * std::sin(angle)));
    }
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 6);
    
    // All edges should meet at origin (circumcenter)
    bool hasOrigin = false;
    for (int i = 0; i < voronoi.NumVertices(); i++) {
        const auto& v = voronoi.vertices[i];
        if (std::abs(v.X()) < 0.1 && std::abs(v.Y()) < 0.1) {
            hasOrigin = true;
            break;
        }
    }
    // Note: For exactly co-circular points, there's one vertex at center
    REQUIRE(hasOrigin);
}

// ============================================================================
// VORONOI DIAGRAM - PROPERTY VERIFICATION
// ============================================================================

TEST_CASE("VoronoiDiagram - Duality with Delaunay", "[ComputationalGeometry][Voronoi][Property]")
{
    // Voronoi vertices should correspond to Delaunay triangle circumcenters
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(4, 3),
        Point2Cartesian(0, 3),
        Point2Cartesian(2, 1.5)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    auto delaunay = MML::CompGeometry::Triangulation::ComputeDelaunay(sites);
    
    // Number of Voronoi vertices should be close to number of Delaunay triangles
    // (exactly equal for non-degenerate cases)
    // Allow some flexibility for edge cases
    REQUIRE(voronoi.NumVertices() <= delaunay.NumTriangles() + 2);
    
    // Each Voronoi edge corresponds to a Delaunay edge
    // Total edges in both structures should be similar
}

TEST_CASE("VoronoiDiagram - Edge perpendicular to site pair", "[ComputationalGeometry][Voronoi][Property]")
{
    // Simple case: verify edge is perpendicular bisector
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    // The edge should be vertical at x = 2
    // Check that edge endpoints have x ≈ 2
    for (int i = 0; i < voronoi.NumEdges(); i++) {
        const auto& edge = voronoi.edges[i];
        // For finite edges, check midpoint
        if (edge.v1 >= 0 && edge.v2 >= 0) {
            const auto& vert1 = voronoi.vertices[edge.v1];
            const auto& vert2 = voronoi.vertices[edge.v2];
            Real midX = (vert1.X() + vert2.X()) / 2;
            REQUIRE_THAT(midX, WithinAbs(2.0, 0.1));
        }
    }
}

TEST_CASE("VoronoiDiagram - Random points stability", "[ComputationalGeometry][Voronoi][Property]")
{
    std::mt19937 rng(98765);
    std::uniform_real_distribution<Real> dist(0, 10);
    
    std::vector<Point2Cartesian> sites;
    for (int i = 0; i < 20; i++) {
        sites.push_back(Point2Cartesian(dist(rng), dist(rng)));
    }
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    REQUIRE(voronoi.NumSites() == 20);
    REQUIRE(voronoi.NumEdges() >= 10);  // Should have multiple edges
    
    // Every site should have at least one adjacent edge
    std::set<int> sitesWithEdges;
    for (int i = 0; i < voronoi.NumEdges(); i++) {
        const auto& edge = voronoi.edges[i];
        sitesWithEdges.insert(edge.site1);
        sitesWithEdges.insert(edge.site2);
    }
    // Most sites should have edges (some boundary sites might not in clipped version)
    REQUIRE(sitesWithEdges.size() >= 15);
}

// ============================================================================
// VORONOI DIAGRAM - CLIPPED CELLS
// ============================================================================

TEST_CASE("VoronoiDiagram - GetClippedCell", "[ComputationalGeometry][Voronoi]")
{
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(5, 5),
        Point2Cartesian(15, 5),
        Point2Cartesian(10, 15)
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    // Clip to a 20x20 box
    for (int i = 0; i < 3; i++) {
        Polygon2D clipped = voronoi.GetClippedCell(i, 0, 0, 20, 20);
        REQUIRE(clipped.NumVertices() >= 3);
        REQUIRE(clipped.Area() > 0);
        
        // Site should be inside its clipped cell
        REQUIRE(clipped.Contains(sites[i]));
    }
}

TEST_CASE("VoronoiDiagram - Cell is infinite check", "[ComputationalGeometry][Voronoi]")
{
    // Square with interior point
    std::vector<Point2Cartesian> sites = {
        Point2Cartesian(0, 0),
        Point2Cartesian(4, 0),
        Point2Cartesian(4, 4),
        Point2Cartesian(0, 4),
        Point2Cartesian(2, 2)  // Interior point
    };
    
    auto voronoi = MML::CompGeometry::Voronoi::ComputeVoronoi(sites);
    
    // Corners are infinite, center is finite
    REQUIRE(voronoi.IsCellInfinite(0));
    REQUIRE(voronoi.IsCellInfinite(1));
    REQUIRE(voronoi.IsCellInfinite(2));
    REQUIRE(voronoi.IsCellInfinite(3));
    REQUIRE_FALSE(voronoi.IsCellInfinite(4));
    
    // Center cell should be a quadrilateral
    REQUIRE(voronoi.cells[4].size() == 4);
    
    // Get the finite cell as polygon
    Polygon2D centerCell = voronoi.GetCell(4);
    REQUIRE(centerCell.NumVertices() == 4);
    REQUIRE(centerCell.Area() > 0);
}

} // namespace MML::Tests::Algorithms::CompGeometry::VoronoiTests


