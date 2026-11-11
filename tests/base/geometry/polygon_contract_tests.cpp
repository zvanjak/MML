#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/Geometry2D/Geometry2DPolygon.h>
#include <mml/algorithms/CompGeometry/RobustPredicates.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using namespace MML;

namespace {
	void CheckPolygonRing(std::vector<Point2Cartesian> vertices, bool simple) {
		for (int direction = 0; direction < 2; ++direction) {
			for (std::size_t rotation = 0; rotation < vertices.size(); ++rotation) {
				INFO("direction " << direction << ", rotation " << rotation);
				CHECK(Polygon2D(vertices).IsSimple() == simple);
				std::rotate(vertices.begin(), vertices.begin() + 1, vertices.end());
			}
			std::reverse(vertices.begin(), vertices.end());
		}
	}
}

TEST_CASE("Polygon topology rejects crossings touches and retraced edges", "[Polygon2D][IsSimple][polygon-contract]")
{
	const std::vector<std::vector<Point2Cartesian>> invalid{
		{{0, 0}, {2, 0}, {1, 0}, {1, 1}},
		{{0, 0}, {2, 2}, {2, 0}, {0, 2}},
		{{0, 0}, {4, 0}, {4, 4}, {2, 4}, {2, 0}, {1, 4}, {0, 4}},
		{{0, 0}, {2, 0}, {1, 1}, {2, 2}, {0, 2}, {1, 1}},
		{{0, 0}, {4, 0}, {4, 2}, {1, 2}, {1, 0}, {3, 0}, {3, 1}, {0, 1}},
		{{0, 0}, {1, 0}, {2, 0}},
		{{0, 0}, {1, 0}, {2, 0}, {3, 0}}
	};
	for (const auto& vertices : invalid)
		CheckPolygonRing(vertices, false);
}

TEST_CASE("Polygon topology preserves concavity and collinear subdivisions", "[Polygon2D][IsSimple][polygon-contract]")
{
	const std::vector<std::vector<Point2Cartesian>> valid{
		{{0, 0}, {4, 0}, {4, 4}, {0, 4}},
		{{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 4}, {0, 4}},
		{{0, 0}, {1, 0}, {4, 0}, {4, 2}, {4, 4}, {0, 4}, {0, 1}},
		{{0, 0}, {4, 0}, {4, 4}, {3, 4}, {3, 1}, {1, 1}, {1, 4}, {0, 4}}
	};
	for (const auto& vertices : valid)
		CheckPolygonRing(vertices, true);
}

TEST_CASE("Polygon topology requires implicit closure and nonzero edges", "[Polygon2D][IsSimple][polygon-contract]")
{
	CHECK_FALSE(Polygon2D{}.IsSimple());
	CHECK_FALSE((Polygon2D{{0, 0}}.IsSimple()));
	CHECK_FALSE((Polygon2D{{0, 0}, {1, 0}}.IsSimple()));
	CheckPolygonRing({{0, 0}, {2, 0}, {2, 2}, {0, 2}, {0, 0}}, false);
	CheckPolygonRing({{0, 0}, {2, 0}, {2, 0}, {2, 2}, {0, 2}}, false);
	CheckPolygonRing({{0, 0}, {-Real(0), Real(0)}, {1, 1}}, false);

	Polygon2D polygon{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
	auto& retainedVertex = polygon[1];
	CHECK(polygon.IsSimple());
	retainedVertex = polygon[2];
	CHECK_FALSE(polygon.IsSimple());
	retainedVertex = Point2Cartesian(2, 0);
	CHECK(polygon.IsSimple());
	CHECK(polygon.NumVertices() == 4);
	CHECK(polygon[1].X() == 2);
	CHECK(polygon[1].Y() == 0);
}

TEST_CASE("Polygon topology does not expand nearly touching boundaries", "[Polygon2D][IsSimple][polygon-contract]")
{
	for (Real gap : {std::nextafter(Real(1), Real(2)) - 1, std::numeric_limits<Real>::denorm_min()}) {
		REQUIRE(gap > 0);
		CheckPolygonRing({{0, 0}, {4, 0}, {4, 4}, {2, 4}, {2, gap}, {1, gap}, {1, 4}, {0, 4}}, true);
	}
	CheckPolygonRing({{0, 0}, {4, 0}, {4, 4}, {2, 4}, {2, 0}, {1, 0}, {1, 4}, {0, 4}}, false);
}

TEST_CASE("Polygon topology is invariant under exact power of two scaling", "[Polygon2D][IsSimple][polygon-contract]")
{
	const int largeExponent = std::numeric_limits<Real>::max_exponent - 12;
	const int smallExponent = std::numeric_limits<Real>::min_exponent + std::numeric_limits<Real>::digits + 12;
	const std::array<std::array<int, 2>, 5> scales{{
		{0, 0}, {largeExponent, largeExponent}, {smallExponent, smallExponent},
		{largeExponent, smallExponent}, {smallExponent, largeExponent}
	}};
	for (const auto& scale : scales) {
		INFO("x exponent " << scale[0] << ", y exponent " << scale[1]);
		for (bool simple : {false, true}) {
			std::vector<Point2Cartesian> vertices = simple
				? std::vector<Point2Cartesian>{{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 4}, {0, 4}}
				: std::vector<Point2Cartesian>{{0, 0}, {2, 0}, {1, 0}, {1, 1}};
			for (auto& vertex : vertices)
				vertex = Point2Cartesian(std::ldexp(vertex.X() + 32, scale[0]),
					std::ldexp(vertex.Y() - 32, scale[1]));
			CheckPolygonRing(vertices, simple);
		}
	}
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	CheckPolygonRing({{0, 0}, {tiny, 0}, {0, tiny}}, true);
}

TEST_CASE("Shared orientation handles cancellation and extreme exponent spreads", "[Polygon2D][Predicates][polygon-contract]")
{
	using Predicate = CompGeometry::RobustPredicates;
	const Real k = std::ldexp(Real(1), (std::numeric_limits<Real>::digits + 1) / 2);
	const Point2Cartesian a(0, 0), b(k, k - 1), c(k + 1, k);
	CHECK(Predicate::Orientation2D(a, b, c) == 1);
	CHECK(Predicate::Orientation2D(a, c, b) == -1);
	CheckPolygonRing({a, b, c}, true);

	const Real huge = std::numeric_limits<Real>::max();
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	const Point2Cartesian first(huge, huge), second(-huge, -huge), third(tiny, 0);
	CHECK(Predicate::Orientation2D(first, second, third) == 1);
	CHECK(Predicate::Orientation2D(first, third, second) == -1);
	CHECK(Predicate::Orientation2D(first, second, Point2Cartesian(0, 0)) == 0);
	CheckPolygonRing({first, second, third}, true);
}

TEST_CASE("Shared orientation agrees with an integer determinant oracle", "[Polygon2D][Predicates][polygon-contract]")
{
	const int exponent = std::numeric_limits<Real>::max_exponent / 2 + 8;
	const std::array<std::array<int, 2>, 5> scales{{
		{0, 0}, {exponent, exponent}, {-exponent, -exponent},
		{exponent, -exponent}, {-exponent, exponent}
	}};
	std::uint32_t state = 0x361f0b95;
	for (int sample = 0; sample < 128; ++sample) {
		std::array<std::int64_t, 6> coordinates{};
		for (auto& coordinate : coordinates) {
			state = state * 1664525u + 1013904223u;
			coordinate = static_cast<std::int64_t>(state % 17) - 8;
		}
		const auto determinant = (coordinates[2] - coordinates[0]) * (coordinates[5] - coordinates[1])
			- (coordinates[3] - coordinates[1]) * (coordinates[4] - coordinates[0]);
		const int expected = (determinant > 0) - (determinant < 0);
		for (const auto& scale : scales) {
			std::array<Point2Cartesian, 3> points;
			for (std::size_t i = 0; i < points.size(); ++i)
				points[i] = Point2Cartesian(std::ldexp(Real(coordinates[2 * i]), scale[0]),
					std::ldexp(Real(coordinates[2 * i + 1]), scale[1]));
			INFO("sample " << sample << ", exponents " << scale[0] << ", " << scale[1]);
			CHECK(CompGeometry::RobustPredicates::Orientation2D(points[0], points[1], points[2]) == expected);
		}
	}
}

TEST_CASE("Polygon topology rejects nonfinite coordinates explicitly", "[Polygon2D][IsSimple][polygon-contract]")
{
	for (Real invalid : {std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity(),
		std::numeric_limits<Real>::quiet_NaN()}) {
		for (int vertex = 0; vertex < 4; ++vertex) {
			Polygon2D polygon{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
			const auto original = polygon[vertex];
			polygon[vertex].X() = invalid;
			CHECK_THROWS_AS(polygon.IsSimple(), GeometryError);
			polygon[vertex] = original;
			polygon[vertex].Y() = invalid;
			CHECK_THROWS_AS(polygon.IsSimple(), GeometryError);
		}
		CHECK_THROWS_AS(CompGeometry::RobustPredicates::Orientation2D(
			Point2Cartesian(0, 0), Point2Cartesian(1, invalid), Point2Cartesian(0, 1)), GeometryError);
	}
}
