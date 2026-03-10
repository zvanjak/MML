///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        kdtree_tests.cpp                                                    ///
///  Description: Unit tests for KD-Tree spatial data structure                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "MMLBase.h"
#include "mml/algorithms/CompGeometry/KDTree.h"

#include <random>
#include <chrono>

using namespace MML;

namespace MML::Tests::Algorithms::CompGeometry::KDTreeTests {
using Catch::Approx;

/******************************************************************************/
/*****                    Test Helpers                                    *****/
/******************************************************************************/

// Generate random 2D points in [0, range) x [0, range)
std::vector<VectorN<Real, 2>> generateRandom2DPoints(int n, Real range = 100.0, int seed = 42) {
	std::mt19937 rng(seed);
	std::uniform_real_distribution<Real> dist(0.0, range);
	
	std::vector<VectorN<Real, 2>> points;
	points.reserve(n);
	for (int i = 0; i < n; ++i) {
		points.push_back({dist(rng), dist(rng)});
	}
	return points;
}

// Generate random 3D points
std::vector<VectorN<Real, 3>> generateRandom3DPoints(int n, Real range = 100.0, int seed = 42) {
	std::mt19937 rng(seed);
	std::uniform_real_distribution<Real> dist(0.0, range);
	
	std::vector<VectorN<Real, 3>> points;
	points.reserve(n);
	for (int i = 0; i < n; ++i) {
		points.push_back({dist(rng), dist(rng), dist(rng)});
	}
	return points;
}

// Generate random Point2Cartesian points
std::vector<Point2Cartesian> generateRandom2DPointsCartesian(int n, Real range = 100.0, int seed = 42) {
	std::mt19937 rng(seed);
	std::uniform_real_distribution<Real> dist(0.0, range);
	
	std::vector<Point2Cartesian> points;
	points.reserve(n);
	for (int i = 0; i < n; ++i) {
		points.emplace_back(dist(rng), dist(rng));
	}
	return points;
}

// Generate random Point3Cartesian points
std::vector<Point3Cartesian> generateRandom3DPointsCartesian(int n, Real range = 100.0, int seed = 42) {
	std::mt19937 rng(seed);
	std::uniform_real_distribution<Real> dist(0.0, range);
	
	std::vector<Point3Cartesian> points;
	points.reserve(n);
	for (int i = 0; i < n; ++i) {
		points.emplace_back(dist(rng), dist(rng), dist(rng));
	}
	return points;
}

// Brute force nearest neighbor for VectorN (for verification)
template<int Dim>
std::pair<size_t, Real> bruteForceNearest(const std::vector<VectorN<Real, Dim>>& points,
                                          const VectorN<Real, Dim>& query) {
	size_t bestIdx = 0;
	Real bestDistSq = std::numeric_limits<Real>::max();
	
	for (size_t i = 0; i < points.size(); ++i) {
		Real distSq = 0;
		for (int d = 0; d < Dim; ++d) {
			Real diff = points[i][d] - query[d];
			distSq += diff * diff;
		}
		if (distSq < bestDistSq) {
			bestDistSq = distSq;
			bestIdx = i;
		}
	}
	
	return {bestIdx, std::sqrt(bestDistSq)};
}

/******************************************************************************/
/*****                    KDTree<2> Tests                                 *****/
/******************************************************************************/

TEST_CASE("KDTree<2>_EmptyTree", "[algorithms][kdtree]") {
	KDTree<2> tree;
	
	REQUIRE(tree.empty());
	REQUIRE(tree.size() == 0);
	
	// Query on empty tree should return not found
	VectorN<Real, 2> query({5.0, 5.0});
	auto result = tree.findNearest(query);
	REQUIRE_FALSE(result.found);
	REQUIRE(result.status == AlgorithmStatus::InvalidInput);
}

TEST_CASE("KDTree<2>_SinglePoint", "[algorithms][kdtree]") {
	KDTree<2> tree;
	std::vector<VectorN<Real, 2>> points = {{5.0, 5.0}};
	
	tree.build(points);
	
	REQUIRE_FALSE(tree.empty());
	REQUIRE(tree.size() == 1);
	
	auto result = tree.findNearest({3.0, 3.0});
	REQUIRE(result.found);
	REQUIRE(result.index == 0);
	REQUIRE(result.distance == Approx(std::sqrt(8.0)).margin(1e-10));
}

TEST_CASE("KDTree<2>_FourCorners", "[algorithms][kdtree]") {
	KDTree<2> tree;
	std::vector<VectorN<Real, 2>> points = {
		{0.0, 0.0},   // index 0: bottom-left
		{10.0, 0.0},  // index 1: bottom-right
		{0.0, 10.0},  // index 2: top-left
		{10.0, 10.0}  // index 3: top-right
	};
	
	tree.build(points);
	
	// Query near each corner
	SECTION("Near bottom-left") {
		auto result = tree.findNearest({1.0, 1.0});
		REQUIRE(result.found);
		REQUIRE(result.index == 0);
	}
	
	SECTION("Near bottom-right") {
		auto result = tree.findNearest({9.0, 1.0});
		REQUIRE(result.found);
		REQUIRE(result.index == 1);
	}
	
	SECTION("Near top-left") {
		auto result = tree.findNearest({1.0, 9.0});
		REQUIRE(result.found);
		REQUIRE(result.index == 2);
	}
	
	SECTION("Near top-right") {
		auto result = tree.findNearest({9.0, 9.0});
		REQUIRE(result.found);
		REQUIRE(result.index == 3);
	}
	
	SECTION("Center - should pick one") {
		auto result = tree.findNearest({5.0, 5.0});
		REQUIRE(result.found);
		// All corners are equidistant, so any is acceptable
		Real expectedDist = std::sqrt(50.0);
		REQUIRE(result.distance == Approx(expectedDist).margin(1e-10));
	}
}

TEST_CASE("KDTree<2>_NearestNeighborCorrectness", "[algorithms][kdtree]") {
	const int N = 1000;
	const int numQueries = 100;
	
	auto points = generateRandom2DPoints(N);
	KDTree<2> tree;
	tree.build(points);
	
	// Generate random queries
	std::mt19937 rng(123);
	std::uniform_real_distribution<Real> dist(0.0, 100.0);
	
	for (int q = 0; q < numQueries; ++q) {
		VectorN<Real, 2> query({dist(rng), dist(rng)});
		
		auto kdResult = tree.findNearest(query);
		auto [bfIndex, bfDist] = bruteForceNearest(points, query);
		
		REQUIRE(kdResult.found);
		REQUIRE(kdResult.index == bfIndex);
		REQUIRE(kdResult.distance == Approx(bfDist).margin(1e-10));
	}
}

TEST_CASE("KDTree<2>_KNearestNeighbors", "[algorithms][kdtree]") {
	std::vector<VectorN<Real, 2>> points = {
		{0.0, 0.0},
		{1.0, 0.0},
		{2.0, 0.0},
		{3.0, 0.0},
		{4.0, 0.0}
	};
	
	KDTree<2> tree;
	tree.build(points);
	
	// Query at 1.5 - nearest should be index 1 (at 1.0) and 2 (at 2.0)
	VectorN<Real, 2> query({1.5, 0.0});
	
	SECTION("k=1") {
		auto result = tree.findKNearest(query, 1);
		REQUIRE(result.indices.size() == 1);
		// Both index 1 and 2 are equidistant
		REQUIRE((result.indices[0] == 1 || result.indices[0] == 2));
		REQUIRE(result.distances[0] == Approx(0.5).margin(1e-10));
	}
	
	SECTION("k=2") {
		auto result = tree.findKNearest(query, 2);
		REQUIRE(result.indices.size() == 2);
		// Should contain indices 1 and 2
		bool hasOne = std::find(result.indices.begin(), result.indices.end(), 1) != result.indices.end();
		bool hasTwo = std::find(result.indices.begin(), result.indices.end(), 2) != result.indices.end();
		REQUIRE(hasOne);
		REQUIRE(hasTwo);
	}
	
	SECTION("k=3") {
		auto result = tree.findKNearest(query, 3);
		REQUIRE(result.indices.size() == 3);
		// Verify sorted by distance
		for (size_t i = 1; i < result.distances.size(); ++i) {
			REQUIRE(result.distances[i] >= result.distances[i-1]);
		}
	}
	
	SECTION("k > n should return all points") {
		auto result = tree.findKNearest(query, 10);
		REQUIRE(result.indices.size() == 5);
	}
}

TEST_CASE("KDTree<2>_RadiusSearch", "[algorithms][kdtree]") {
	std::vector<VectorN<Real, 2>> points = {
		{0.0, 0.0},   // 0
		{1.0, 0.0},   // 1
		{0.0, 1.0},   // 2
		{10.0, 10.0}, // 3
		{10.5, 10.0}  // 4
	};
	
	KDTree<2> tree;
	tree.build(points);
	
	SECTION("Radius captures cluster") {
		VectorN<Real, 2> query({0.0, 0.0});
		auto result = tree.findInRadius(query, 1.5);
		
		// Should find indices 0, 1, 2 (all within distance 1.5 from origin)
		REQUIRE(result.indices.size() == 3);
	}
	
	SECTION("Small radius") {
		VectorN<Real, 2> query({0.0, 0.0});
		auto result = tree.findInRadius(query, 0.5);
		
		// Only index 0 is within 0.5
		REQUIRE(result.indices.size() == 1);
		REQUIRE(result.indices[0] == 0);
	}
	
	SECTION("No points in radius") {
		VectorN<Real, 2> query({50.0, 50.0});
		auto result = tree.findInRadius(query, 1.0);
		
		REQUIRE(result.indices.empty());
	}
	
	SECTION("Negative radius is invalid") {
		VectorN<Real, 2> query({0.0, 0.0});
		auto result = tree.findInRadius(query, -1.0);
		
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
	}
}

TEST_CASE("KDTree<2>_BoxSearch", "[algorithms][kdtree]") {
	std::vector<VectorN<Real, 2>> points = {
		{0.0, 0.0},
		{1.0, 1.0},
		{2.0, 2.0},
		{5.0, 5.0},
		{8.0, 8.0}
	};
	
	KDTree<2> tree;
	tree.build(points);
	
	SECTION("Box captures subset") {
		VectorN<Real, 2> minCorner({0.5, 0.5});
		VectorN<Real, 2> maxCorner({5.5, 5.5});
		
		auto result = tree.findInBox(minCorner, maxCorner);
		
		// Should find indices 1, 2, 3
		REQUIRE(result.size() == 3);
	}
	
	SECTION("Tiny box") {
		VectorN<Real, 2> minCorner({0.9, 0.9});
		VectorN<Real, 2> maxCorner({1.1, 1.1});
		
		auto result = tree.findInBox(minCorner, maxCorner);
		
		REQUIRE(result.size() == 1);
		REQUIRE(result[0] == 1);
	}
	
	SECTION("Empty box") {
		VectorN<Real, 2> minCorner({100.0, 100.0});
		VectorN<Real, 2> maxCorner({110.0, 110.0});
		
		auto result = tree.findInBox(minCorner, maxCorner);
		
		REQUIRE(result.empty());
	}
}

/******************************************************************************/
/*****                    KDTree<3> Tests                                 *****/
/******************************************************************************/

TEST_CASE("KDTree<3>_BasicOperations", "[algorithms][kdtree]") {
	std::vector<VectorN<Real, 3>> points = {
		{0.0, 0.0, 0.0},
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{1.0, 1.0, 1.0}
	};
	
	KDTree<3> tree;
	tree.build(points);
	
	REQUIRE(tree.size() == 5);
	
	// Nearest to origin should be origin
	auto result = tree.findNearest({0.1, 0.1, 0.1});
	REQUIRE(result.found);
	REQUIRE(result.index == 0);
	
	// Nearest to (1,1,1) should be (1,1,1)
	result = tree.findNearest({0.9, 0.9, 0.9});
	REQUIRE(result.index == 4);
}

TEST_CASE("KDTree<3>_NearestNeighborCorrectness", "[algorithms][kdtree]") {
	const int N = 500;
	const int numQueries = 50;
	
	auto points = generateRandom3DPoints(N);
	KDTree<3> tree;
	tree.build(points);
	
	std::mt19937 rng(456);
	std::uniform_real_distribution<Real> dist(0.0, 100.0);
	
	for (int q = 0; q < numQueries; ++q) {
		VectorN<Real, 3> query({dist(rng), dist(rng), dist(rng)});
		
		auto kdResult = tree.findNearest(query);
		auto [bfIndex, bfDist] = bruteForceNearest(points, query);
		
		REQUIRE(kdResult.found);
		REQUIRE(kdResult.index == bfIndex);
		REQUIRE(kdResult.distance == Approx(bfDist).margin(1e-10));
	}
}

/******************************************************************************/
/*****                    KDTree2D Specialization Tests                   *****/
/******************************************************************************/

TEST_CASE("KDTree2D_Point2CartesianAPI", "[algorithms][kdtree]") {
	std::vector<Point2Cartesian> points = {
		{0.0, 0.0},
		{1.0, 0.0},
		{0.0, 1.0},
		{1.0, 1.0}
	};
	
	KDTree2D tree;
	tree.build(points);
	
	REQUIRE(tree.size() == 4);
	
	SECTION("findNearest") {
		auto result = tree.findNearest({0.1, 0.1});
		REQUIRE(result.found);
		REQUIRE(result.index == 0);
	}
	
	SECTION("findKNearest") {
		auto result = tree.findKNearest({0.5, 0.5}, 2);
		REQUIRE(result.indices.size() == 2);
	}
	
	SECTION("findInRadius") {
		auto result = tree.findInRadius({0.0, 0.0}, 1.5);
		// All 4 points are within 1.5 of origin
		REQUIRE(result.indices.size() == 4);
	}
	
	SECTION("findInBox") {
		auto result = tree.findInBox({0.0, 0.0}, {0.5, 0.5});
		REQUIRE(result.size() == 1);
		REQUIRE(result[0] == 0);
	}
	
	SECTION("point accessor") {
		const Point2Cartesian& p = tree.point(2);
		REQUIRE(p.X() == Approx(0.0));
		REQUIRE(p.Y() == Approx(1.0));
	}
}

TEST_CASE("KDTree2D_VerifyAgainstBruteForce", "[algorithms][kdtree]") {
	const int N = 500;
	auto points = generateRandom2DPointsCartesian(N);
	
	KDTree2D tree;
	tree.build(points);
	
	std::mt19937 rng(789);
	std::uniform_real_distribution<Real> dist(0.0, 100.0);
	
	for (int q = 0; q < 50; ++q) {
		Point2Cartesian query(dist(rng), dist(rng));
		
		auto kdResult = tree.findNearest(query);
		auto bfResult = FindNearestBruteForce2D(points, query);
		
		REQUIRE(kdResult.index == bfResult.index);
		REQUIRE(kdResult.distance == Approx(bfResult.distance).margin(1e-10));
	}
}

/******************************************************************************/
/*****                    KDTree3D Specialization Tests                   *****/
/******************************************************************************/

TEST_CASE("KDTree3D_Point3CartesianAPI", "[algorithms][kdtree]") {
	std::vector<Point3Cartesian> points = {
		{0.0, 0.0, 0.0},
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{1.0, 1.0, 1.0}
	};
	
	KDTree3D tree;
	tree.build(points);
	
	REQUIRE(tree.size() == 5);
	
	SECTION("findNearest") {
		auto result = tree.findNearest({0.1, 0.1, 0.1});
		REQUIRE(result.found);
		REQUIRE(result.index == 0);
	}
	
	SECTION("findKNearest") {
		auto result = tree.findKNearest({0.5, 0.5, 0.5}, 3);
		REQUIRE(result.indices.size() == 3);
	}
	
	SECTION("findInRadius") {
		auto result = tree.findInRadius({0.0, 0.0, 0.0}, 1.1);
		// Origin, (1,0,0), (0,1,0), (0,0,1) are within 1.1
		REQUIRE(result.indices.size() == 4);
	}
	
	SECTION("point accessor") {
		const Point3Cartesian& p = tree.point(4);
		REQUIRE(p.X() == Approx(1.0));
		REQUIRE(p.Y() == Approx(1.0));
		REQUIRE(p.Z() == Approx(1.0));
	}
}

TEST_CASE("KDTree3D_VerifyAgainstBruteForce", "[algorithms][kdtree]") {
	const int N = 500;
	auto points = generateRandom3DPointsCartesian(N);
	
	KDTree3D tree;
	tree.build(points);
	
	std::mt19937 rng(321);
	std::uniform_real_distribution<Real> dist(0.0, 100.0);
	
	for (int q = 0; q < 50; ++q) {
		Point3Cartesian query(dist(rng), dist(rng), dist(rng));
		
		auto kdResult = tree.findNearest(query);
		auto bfResult = FindNearestBruteForce3D(points, query);
		
		REQUIRE(kdResult.index == bfResult.index);
		REQUIRE(kdResult.distance == Approx(bfResult.distance).margin(1e-10));
	}
}

/******************************************************************************/
/*****                    Statistics Tests                                *****/
/******************************************************************************/

TEST_CASE("KDTree_Statistics", "[algorithms][kdtree]") {
	const int N = 1000;
	auto points = generateRandom2DPoints(N);
	
	KDTree<2> tree;
	tree.build(points);
	
	const auto& stats = tree.stats();
	
	REQUIRE(stats.num_points == N);
	REQUIRE(stats.num_nodes > 0);
	REQUIRE(stats.num_leaf_nodes > 0);
	REQUIRE(stats.tree_depth > 0);
	REQUIRE(stats.build_time_ms >= 0.0);
	
	// For a balanced tree, depth should be roughly log2(n/leaf_size)
	// With default max_leaf_size=10, depth should be around 7-10 for 1000 points
	REQUIRE(stats.tree_depth <= 20);  // Very generous upper bound
}

/******************************************************************************/
/*****                    Configuration Tests                             *****/
/******************************************************************************/

TEST_CASE("KDTree_Configurations", "[algorithms][kdtree]") {
	const int N = 500;
	auto points = generateRandom2DPoints(N);
	
	SECTION("Default config") {
		KDTree<2> tree;
		tree.build(points, KDTreeConfig::Default());
		REQUIRE(tree.size() == N);
	}
	
	SECTION("High performance config") {
		KDTree<2> tree;
		tree.build(points, KDTreeConfig::HighPerformance());
		// High performance uses larger leaves, so fewer nodes
		const auto& stats = tree.stats();
		REQUIRE(stats.num_nodes > 0);
	}
	
	SECTION("Low memory config") {
		KDTree<2> tree;
		tree.build(points, KDTreeConfig::LowMemory());
		// Low memory uses smaller leaves, so more nodes
		const auto& stats = tree.stats();
		REQUIRE(stats.num_nodes > 0);
	}
}

/******************************************************************************/
/*****                    Edge Cases                                      *****/
/******************************************************************************/

TEST_CASE("KDTree_EdgeCases", "[algorithms][kdtree]") {
	SECTION("All points same location") {
		std::vector<VectorN<Real, 2>> points(100, {5.0, 5.0});
		KDTree<2> tree;
		tree.build(points);
		
		auto result = tree.findNearest({5.0, 5.0});
		REQUIRE(result.found);
		REQUIRE(result.distance == Approx(0.0).margin(1e-10));
	}
	
	SECTION("Points on a line") {
		std::vector<VectorN<Real, 2>> points;
		for (int i = 0; i < 100; ++i) {
			points.push_back({static_cast<Real>(i), 0.0});
		}
		
		KDTree<2> tree;
		tree.build(points);
		
		auto result = tree.findNearest({50.5, 0.0});
		REQUIRE(result.found);
		REQUIRE((result.index == 50 || result.index == 51));
	}
	
	SECTION("k-NN with k=0") {
		std::vector<VectorN<Real, 2>> points = {{1.0, 1.0}};
		KDTree<2> tree;
		tree.build(points);
		
		auto result = tree.findKNearest({0.0, 0.0}, 0);
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
	}
	
	SECTION("Clear and rebuild") {
		std::vector<VectorN<Real, 2>> points = {{1.0, 1.0}, {2.0, 2.0}};
		KDTree<2> tree;
		tree.build(points);
		
		REQUIRE(tree.size() == 2);
		
		tree.clear();
		REQUIRE(tree.empty());
		REQUIRE(tree.size() == 0);
		
		// Rebuild with different points
		std::vector<VectorN<Real, 2>> newPoints = {{5.0, 5.0}};
		tree.build(newPoints);
		REQUIRE(tree.size() == 1);
	}
}

/******************************************************************************/
/*****                    Performance Tests (Smoke Tests)                 *****/
/******************************************************************************/

TEST_CASE("KDTree_Performance_Build", "[algorithms][kdtree][.performance]") {
	// Build large tree - marked with .performance tag, skipped in regular runs
	const int N = 100000;
	auto points = generateRandom3DPoints(N);
	
	KDTree<3> tree;
	
	auto start = std::chrono::high_resolution_clock::now();
	tree.build(points);
	auto end = std::chrono::high_resolution_clock::now();
	
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	
	INFO("Build time for " << N << " 3D points: " << duration << " ms");
	REQUIRE(tree.size() == N);
	
	// Queries should still be fast
	std::mt19937 rng(999);
	std::uniform_real_distribution<Real> dist(0.0, 100.0);
	
	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 10000; ++i) {
		VectorN<Real, 3> query({dist(rng), dist(rng), dist(rng)});
		auto result = tree.findNearest(query);
		REQUIRE(result.found);
	}
	end = std::chrono::high_resolution_clock::now();
	
	auto queryTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	INFO("10000 queries: " << queryTime << " ms");
}

} // namespace MML::Tests::Algorithms::CompGeometry::KDTreeTests
