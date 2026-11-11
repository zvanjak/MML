#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/CompGeometry/KDTree.h>
#endif
#include <catch2/catch_all.hpp>
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <numeric>

using namespace MML;

namespace {
	Real LargeKDScale() { return std::scalbn(Real(1), std::numeric_limits<Real>::max_exponent * 2 / 3); }
	Real SmallKDScale() { return std::scalbn(Real(1), std::numeric_limits<Real>::min_exponent * 2 / 3); }

	template<int Dim>
	void CheckKDOracle(Real scale, bool balanced) {
		std::vector<VectorN<Real, Dim>> points;
		for (int i = 0; i < 37; ++i) {
			VectorN<Real, Dim> point;
			for (int d = 0; d < Dim; ++d)
				point[d] = Real(d == 0 ? 4 * (i - 18) : (i * (d + 1)) % 13 - 6) * scale;
			points.push_back(point);
		}
		KDTreeConfig config;
		config.max_leaf_size = 1;
		config.balanced = balanced;
		KDTree<Dim> tree;
		tree.build(points, config);
		for (int q = -3; q <= 3; ++q) {
			VectorN<Real, Dim> query;
			for (int d = 0; d < Dim; ++d) query[d] = (Real(q) + Real(d + 1) / 4) * scale;
			std::vector<Real> squared;
			for (const auto& point : points) {
				Real sum = 0;
				for (int d = 0; d < Dim; ++d) {
					const Real difference = point[d] / scale - query[d] / scale;
					sum += difference * difference;
				}
				squared.push_back(sum);
			}
			auto ordered = squared;
			std::sort(ordered.begin(), ordered.end());
			const auto nearest = tree.findNearest(query);
			REQUIRE(nearest.status == AlgorithmStatus::Success);
			REQUIRE(nearest.found);
			REQUIRE(nearest.index < points.size());
			CHECK(squared[nearest.index] == ordered.front());
			CHECK(IsNearlyEqual(nearest.distance / scale, std::sqrt(ordered.front()), Precision::EqualityTolerance));

			const auto neighbors = tree.findKNearest(query, 7);
			REQUIRE(neighbors.status == AlgorithmStatus::Success);
			REQUIRE(neighbors.indices.size() == 7);
			REQUIRE(neighbors.distances.size() == neighbors.indices.size());
			CHECK(std::is_sorted(neighbors.distances.begin(), neighbors.distances.end()));
			auto unique = neighbors.indices;
			std::sort(unique.begin(), unique.end());
			CHECK(std::adjacent_find(unique.begin(), unique.end()) == unique.end());
			for (size_t i = 0; i < neighbors.indices.size(); ++i) {
				REQUIRE(neighbors.indices[i] < points.size());
				CHECK(squared[neighbors.indices[i]] == ordered[i]);
				CHECK(IsNearlyEqual(neighbors.distances[i] / scale, std::sqrt(ordered[i]), Precision::EqualityTolerance));
			}
			const Real radius = 9;
			const auto inside = tree.findInRadius(query, radius * scale);
			REQUIRE(inside.status == AlgorithmStatus::Success);
			REQUIRE(inside.indices.size() == inside.distances.size());
			std::vector<size_t> expected;
			for (size_t i = 0; i < squared.size(); ++i)
				if (squared[i] <= radius * radius) expected.push_back(i);
			auto actual = inside.indices;
			std::sort(actual.begin(), actual.end());
			CHECK(actual == expected);
			for (size_t i = 0; i < inside.indices.size(); ++i) {
				REQUIRE(inside.indices[i] < points.size());
				CHECK(IsNearlyEqual(inside.distances[i] / scale, std::sqrt(squared[inside.indices[i]]), Precision::EqualityTolerance));
			}
			VectorN<Real, Dim> lower, upper;
			for (int d = 0; d < Dim; ++d) {
				lower[d] = -5 * scale;
				upper[d] = 5 * scale;
			}
			expected.clear();
			for (size_t i = 0; i < points.size(); ++i) {
				bool include = true;
				for (int d = 0; d < Dim; ++d)
					include = include && points[i][d] >= lower[d] && points[i][d] <= upper[d];
				if (include) expected.push_back(i);
			}
			actual = tree.findInBox(lower, upper);
			std::sort(actual.begin(), actual.end());
			CHECK(actual == expected);
		}
	}
}

TEST_CASE("KDTree contracts reject invalid leaf limits even for empty rebuilds",
	"[kdtree][kdtree-contract][kdtree-red]") {
	for (const int leafSize : {0, -1, INT_MIN}) {
		KDTreeConfig config;
		config.max_leaf_size = leafSize;
		KDTree<2> tree;
		tree.build({{7, 2}});
		CHECK_THROWS_AS(tree.build({}, config), ArgumentError);
		CHECK(tree.size() == 1);
		CHECK(tree.findNearest({7, 2}).found);
		KDTree2D wrapper;
		wrapper.build({{7, 2}});
		CHECK_THROWS_AS(wrapper.build({}, config), ArgumentError);
		CHECK(wrapper.size() == 1);
		if (wrapper.size() == 1) CHECK(wrapper.point(0).X() == 7);
	}
}

TEST_CASE("KDTree contracts large and tiny distances keep actual candidates and radius membership",
	"[kdtree][kdtree-contract][kdtree-red]") {
	const Real reportedScale = std::pow(Real(10), Real(std::min(200, std::numeric_limits<Real>::max_exponent10 * 2 / 3)));
	for (const Real scale : {LargeKDScale(), SmallKDScale(), reportedScale})
		for (const bool balanced : {true, false}) {
			KDTreeConfig config;
			config.max_leaf_size = 1;
			config.balanced = balanced;
			KDTree<2> tree;
			tree.build({{2 * scale, 0}, {scale, 0}}, config);
			const auto nearest = tree.findNearest({0, 0});
			CHECK(nearest.status == AlgorithmStatus::Success);
			CHECK(nearest.found);
			CHECK(nearest.index == 1);
			CHECK(nearest.distance == scale);
			const auto neighbors = tree.findKNearest({0, 0}, 2);
			CHECK(neighbors.indices == std::vector<size_t>{1, 0});
			CHECK(neighbors.distances == std::vector<Real>{scale, 2 * scale});
			const auto radius = tree.findInRadius({0, 0}, scale);
			CHECK(radius.indices == std::vector<size_t>{1});
			CHECK(radius.distances == std::vector<Real>{scale});
			const auto brute2 = FindNearestBruteForce2D({{2 * scale, 0}, {scale, 0}}, {0, 0});
			CHECK(brute2.found);
			CHECK(brute2.index == 1);
			CHECK(brute2.distance == scale);
			const auto brute3 = FindNearestBruteForce3D({{2 * scale, 0, 0}, {scale, 0, 0}}, {0, 0, 0});
			CHECK(brute3.found);
			CHECK(brute3.index == 1);
			CHECK(brute3.distance == scale);
		}
}

TEST_CASE("KDTree contracts independent scaled oracles cover every query and split mode",
	"[kdtree][kdtree-contract][kdtree-red]") {
	for (const Real scale : {Real(1), LargeKDScale(), SmallKDScale()})
		for (const bool balanced : {true, false}) {
			CheckKDOracle<2>(scale, balanced);
			CheckKDOracle<3>(scale, balanced);
			CheckKDOracle<7>(scale, balanced);
		}
}

TEST_CASE("KDTree contracts unrepresentable requested distances have no usable payload",
	"[kdtree][kdtree-contract][kdtree-red]") {
	const Real maximum = std::numeric_limits<Real>::max();
	KDTree<2> tree;
	tree.build({{maximum, 0}});
	const auto nearest = tree.findNearest({-maximum, 0});
	CHECK(nearest.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(nearest.found);
	CHECK(nearest.index == SIZE_MAX);
	CHECK(nearest.nodes_visited > 0);
	const auto neighbors = tree.findKNearest({-maximum, 0}, 1);
	CHECK(neighbors.status == AlgorithmStatus::NumericalInstability);
	CHECK(neighbors.indices.empty());
	CHECK(neighbors.distances.empty());
	CHECK(neighbors.nodes_visited > 0);
	const auto radius = tree.findInRadius({-maximum, 0}, maximum);
	CHECK(radius.status == AlgorithmStatus::Success);
	CHECK(radius.indices.empty());
	CHECK(radius.distances.empty());
	const auto brute2 = FindNearestBruteForce2D({{maximum, 0}}, {-maximum, 0});
	CHECK(brute2.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(brute2.found);
	CHECK(brute2.index == SIZE_MAX);
	const auto brute3 = FindNearestBruteForce3D({{maximum, 0, 0}}, {-maximum, 0, 0});
	CHECK(brute3.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(brute3.found);
	CHECK(brute3.index == SIZE_MAX);
}

TEST_CASE("KDTree contracts unused overflowing distances do not poison finite selected results",
	"[kdtree][kdtree-contract][kdtree-red]") {
	const Real coordinate = std::numeric_limits<Real>::max() * (Real(3) / 4);
	KDTree<2> tree;
	tree.build({{coordinate, coordinate}, {1, 0}});
	CHECK(tree.findNearest({0, 0}).index == 1);
	const auto one = tree.findKNearest({0, 0}, 1);
	CHECK(one.status == AlgorithmStatus::Success);
	CHECK(one.indices == std::vector<size_t>{1});
	CHECK(one.distances == std::vector<Real>{1});
	const auto both = tree.findKNearest({0, 0}, 2);
	CHECK(both.status == AlgorithmStatus::NumericalInstability);
	CHECK(both.indices.empty());
	CHECK(both.distances.empty());
	const auto inside = tree.findInRadius({0, 0}, 2);
	CHECK(inside.status == AlgorithmStatus::Success);
	CHECK(inside.indices == std::vector<size_t>{1});
	CHECK(inside.distances == std::vector<Real>{1});
}

TEST_CASE("KDTree contracts finite maximum and subnormal distances remain usable",
	"[kdtree][kdtree-contract][kdtree-red]") {
	for (const Real distance : {std::numeric_limits<Real>::max(), std::numeric_limits<Real>::denorm_min(), Real(0)}) {
		KDTree<3> tree;
		tree.build({{distance, 0, 0}});
		const auto nearest = tree.findNearest({0, 0, 0});
		CHECK(nearest.status == AlgorithmStatus::Success);
		CHECK(nearest.found);
		CHECK(nearest.index == 0);
		CHECK(nearest.distance == distance);
		const auto neighbors = tree.findKNearest({0, 0, 0}, 1);
		CHECK(neighbors.indices == std::vector<size_t>{0});
		CHECK(neighbors.distances == std::vector<Real>{distance});
		const auto boundary = tree.findInRadius({0, 0, 0}, distance);
		CHECK(boundary.indices == std::vector<size_t>{0});
		CHECK(boundary.distances == std::vector<Real>{distance});
		const auto zero = tree.findInRadius({0, 0, 0}, 0);
		CHECK(zero.indices.empty() == (distance != 0));
	}
}

TEST_CASE("KDTree contracts midpoint splitting is honored and collapsed partitions terminate",
	"[kdtree][kdtree-contract][kdtree-red]") {
	const std::vector<VectorN<Real, 1>> skewed{{0}, {1}, {2}, {4}, {8}, {16}, {32}, {64}, {128}};
	KDTreeConfig config;
	config.max_leaf_size = 1;
	KDTree<1> median, midpoint;
	median.build(skewed, config);
	config.balanced = false;
	midpoint.build(skewed, config);
	CHECK(midpoint.stats().tree_depth > median.stats().tree_depth);
	for (const Real query : {Real(0), Real(3), Real(127)}) {
		const auto a = median.findNearest({query}), b = midpoint.findNearest({query});
		CHECK(a.found);
		CHECK(b.found);
		CHECK(a.distance == b.distance);
	}
	midpoint.build(std::vector<VectorN<Real, 1>>(31, {1}), config);
	CHECK(midpoint.stats().tree_depth == 5);
	CHECK(midpoint.findKNearest({1}, 40).indices.size() == 31);
	const Real adjacent = std::nextafter(Real(1), Real(2));
	midpoint.build({{1}, {adjacent}}, config);
	CHECK(midpoint.stats().tree_depth == 1);
	CHECK(midpoint.findNearest({adjacent}).distance == 0);
	const Real maximum = std::numeric_limits<Real>::max();
	midpoint.build({{-maximum}, {maximum}}, config);
	CHECK(midpoint.stats().tree_depth == 1);
	CHECK(midpoint.findNearest({maximum}).distance == 0);
	CHECK(midpoint.findInRadius({0}, maximum).indices.size() == 2);
	midpoint.build({{maximum / 2}, {maximum}}, config);
	CHECK(midpoint.stats().tree_depth == 1);
	CHECK(midpoint.findNearest({maximum}).index == 1);
}

TEST_CASE("KDTree contracts invalid coordinates and query controls fail explicitly",
	"[kdtree][kdtree-contract][kdtree-red]") {
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		KDTree<2> tree;
		tree.build({{7, 2}});
		CHECK_THROWS_AS(tree.build({{invalid, 0}}), NumericInputError);
		CHECK(tree.size() == 1);
		CHECK(tree.point(0)[0] == 7);
		const auto nearest = tree.findNearest({invalid, 0});
		CHECK(nearest.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(nearest.found);
		CHECK(nearest.index == SIZE_MAX);
		const auto neighbors = tree.findKNearest({invalid, 0}, 1);
		CHECK(neighbors.status == AlgorithmStatus::InvalidInput);
		CHECK(neighbors.indices.empty());
		CHECK(neighbors.distances.empty());
		CHECK(tree.findInRadius({invalid, 0}, 1).status == AlgorithmStatus::InvalidInput);
		CHECK(tree.findInRadius({0, 0}, invalid).status == AlgorithmStatus::InvalidInput);
		CHECK_THROWS_AS(tree.findInBox({invalid, 0}, {1, 1}), NumericInputError);
		CHECK_THROWS_AS(tree.findInBox({2, 0}, {1, 1}), ArgumentError);
		const auto brute2 = FindNearestBruteForce2D({{invalid, 0}}, {0, 0});
		CHECK(brute2.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(brute2.found);
		CHECK(brute2.index == SIZE_MAX);
		const auto brute3 = FindNearestBruteForce3D({{1, 0, 0}}, {invalid, 0, 0});
		CHECK(brute3.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(brute3.found);
		CHECK(brute3.index == SIZE_MAX);
	}
}

TEST_CASE("KDTree contracts singleton invalid leaf limits are rejected without recursion",
	"[kdtree][kdtree-contract]") {
	for (const int leafSize : {0, -1, INT_MIN}) {
		KDTreeConfig config;
		config.max_leaf_size = leafSize;
		KDTree<2> tree;
		CHECK_THROWS_AS(tree.build({{1, 2}}, config), ArgumentError);
		CHECK(tree.empty());
	}
}

TEST_CASE("KDTree contracts wrappers retain ownership across rejected builds and moves",
	"[kdtree][kdtree-contract]") {
	const Real scale = LargeKDScale();
	KDTreeConfig config;
	config.max_leaf_size = 1;
	config.balanced = false;
	std::vector<Point2Cartesian> points2{{2 * scale, 0}, {scale, 0}};
	KDTree2D tree2;
	tree2.build(points2, config);
	points2[1] = Point2Cartesian(0, 0);
	CHECK(tree2.point(1).X() == scale);
	CHECK(tree2.findNearest({0, 0}).index == 1);
	CHECK(tree2.findKNearest({0, 0}, 1).distances == std::vector<Real>{scale});
	CHECK(tree2.findInRadius({0, 0}, scale).indices == std::vector<size_t>{1});
	CHECK_THROWS_AS(tree2.build({{std::numeric_limits<Real>::infinity(), 0}}), NumericInputError);
	CHECK(tree2.point(1).X() == scale);
	std::vector<Point3Cartesian> points3{{2 * scale, 0, 0}, {scale, 0, 0}};
	KDTree3D tree3;
	tree3.build(points3, config);
	auto moved = std::move(tree3);
	points3.clear();
	CHECK(moved.point(1).X() == scale);
	CHECK(moved.findNearest({0, 0, 0}).index == 1);
	CHECK(moved.findKNearest({0, 0, 0}, 1).distances == std::vector<Real>{scale});
	CHECK(moved.findInRadius({0, 0, 0}, scale).indices == std::vector<size_t>{1});
	config.max_leaf_size = 0;
	CHECK_THROWS_AS(moved.build({{0, 0, 0}}, config), ArgumentError);
	CHECK(moved.size() == 2);
	CHECK(moved.point(1).X() == scale);
	moved.clear();
	CHECK(moved.empty());
	CHECK(moved.size() == 0);
}

TEST_CASE("KDTree contracts empty and duplicate query controls retain their semantics",
	"[kdtree][kdtree-contract]") {
	KDTree<2> tree;
	CHECK_FALSE(tree.findNearest({0, 0}).found);
	CHECK(tree.findNearest({0, 0}).status == AlgorithmStatus::InvalidInput);
	CHECK(tree.findKNearest({0, 0}, 3).indices.empty());
	CHECK(tree.findKNearest({0, 0}, 3).status == AlgorithmStatus::Success);
	CHECK(tree.findKNearest({0, 0}, 0).status == AlgorithmStatus::InvalidInput);
	CHECK(tree.findInRadius({0, 0}, 1).indices.empty());
	CHECK(tree.findInRadius({0, 0}, 1).status == AlgorithmStatus::Success);
	CHECK(tree.findInRadius({0, 0}, -1).status == AlgorithmStatus::InvalidInput);
	CHECK(tree.findInBox({0, 0}, {1, 1}).empty());
	CHECK_FALSE(FindNearestBruteForce2D({}, {0, 0}).found);
	CHECK_FALSE(FindNearestBruteForce3D({}, {0, 0, 0}).found);
	KDTreeConfig config;
	config.max_leaf_size = INT_MAX;
	tree.build({{0, 0}, {0, 0}, {1, 0}}, config);
	CHECK(tree.stats().num_nodes == 1);
	const auto duplicates = tree.findKNearest({0, 0}, 20);
	CHECK(duplicates.indices.size() == 3);
	CHECK(duplicates.distances == std::vector<Real>{0, 0, 1});
	CHECK(tree.findInRadius({0, 0}, 0).indices.size() == 2);
	tree.build({});
	CHECK(tree.empty());
	CHECK(tree.size() == 0);
	CHECK(tree.stats().num_nodes == 0);
}
