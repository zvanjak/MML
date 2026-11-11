#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/CompGeometry/Triangulation.h>
#endif
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <vector>

namespace {
using namespace MML;
using namespace MML::CompGeometry;
using Edge = std::pair<int, int>;

long double Cross(const Point2Cartesian& a, const Point2Cartesian& b, const Point2Cartesian& c) {
	return (static_cast<long double>(b.X()) - a.X()) * (static_cast<long double>(c.Y()) - a.Y())
		- (static_cast<long double>(b.Y()) - a.Y()) * (static_cast<long double>(c.X()) - a.X());
}

std::vector<Point2Cartesian> Transform(const std::vector<Point2Cartesian>& points, int exponent, Real shift) {
	std::vector<Point2Cartesian> result;
	for (const auto& p : points)
		result.emplace_back(std::ldexp(p.X() + shift, exponent), std::ldexp(p.Y() + shift, exponent));
	return result;
}

void CheckMesh(const DelaunayTriangulation& mesh, const std::vector<Point2Cartesian>& reference) {
	REQUIRE(mesh.points.size() == reference.size());
	REQUIRE(mesh.neighbors.size() == mesh.triangles.size());
	REQUIRE_FALSE(mesh.triangles.empty());
	std::map<Edge, std::vector<std::pair<int, int>>> incidence;
	std::set<int> used;
	long double totalArea = 0;
	for (int t = 0; t < mesh.NumTriangles(); ++t) {
		const auto& triangle = mesh.triangles[t];
		for (int v : triangle) {
			REQUIRE(v >= 0);
			REQUIRE(v < mesh.NumPoints());
			used.insert(v);
		}
		const auto& a = reference[triangle[0]];
		const auto& b = reference[triangle[1]];
		const auto& c = reference[triangle[2]];
		const auto area = Cross(a, b, c);
		CHECK(area > 0);
		totalArea += area;
		for (int e = 0; e < 3; ++e)
			incidence[std::minmax(triangle[e], triangle[(e + 1) % 3])].push_back({t, e});

		// Independent, normalized determinant oracle: never form extreme-scale powers.
		for (int v = 0; v < mesh.NumPoints(); ++v) {
			const auto& d = reference[v];
			const long double ax = static_cast<long double>(a.X()) - d.X(), ay = static_cast<long double>(a.Y()) - d.Y();
			const long double bx = static_cast<long double>(b.X()) - d.X(), by = static_cast<long double>(b.Y()) - d.Y();
			const long double cx = static_cast<long double>(c.X()) - d.X(), cy = static_cast<long double>(c.Y()) - d.Y();
			const auto first = (ax * ax + ay * ay) * (bx * cy - by * cx);
			const auto second = (bx * bx + by * by) * (cx * ay - cy * ax);
			const auto third = (cx * cx + cy * cy) * (ax * by - ay * bx);
			const auto budget = 128 * std::numeric_limits<long double>::epsilon()
				* (std::abs(first) + std::abs(second) + std::abs(third));
			CHECK(first + second + third <= budget);
		}
	}
	REQUIRE(used.size() == reference.size());
	std::size_t boundary = 0;
	long double boundaryArea = 0;
	for (const auto& [edge, faces] : incidence) {
		REQUIRE(faces.size() <= 2);
		for (const auto& [t, e] : faces) {
			if (faces.size() == 1) {
				CHECK(mesh.neighbors[t][e] == -1);
				++boundary;
				const auto& tri = mesh.triangles[t];
				const auto& a = reference[tri[e]];
				const auto& b = reference[tri[(e + 1) % 3]];
				boundaryArea += static_cast<long double>(a.X()) * b.Y() - static_cast<long double>(a.Y()) * b.X();
				for (const auto& p : reference)
					CHECK(Cross(a, b, p) >= 0);
			} else {
				const auto other = faces[0].first == t ? faces[1] : faces[0];
				CHECK(mesh.neighbors[t][e] == other.first);
				CHECK(mesh.neighbors[other.first][other.second] == t);
				const auto& tri = mesh.triangles[t];
				const auto& adjacent = mesh.triangles[other.first];
				CHECK(tri[e] == adjacent[(other.second + 1) % 3]);
				CHECK(tri[(e + 1) % 3] == adjacent[other.second]);
			}
		}
	}
	CHECK(mesh.GetConvexHullEdges().size() == boundary);
	CHECK(mesh.triangles.size() == 2 * reference.size() - 2 - boundary);
	CHECK(std::abs(totalArea - boundaryArea) <= 128 * std::numeric_limits<long double>::epsilon() * std::abs(boundaryArea));
}
}

TEST_CASE("Delaunay retains skinny three-point topology in every input order", "[triangulation-contract][Delaunay]") {
	const std::vector<Point2Cartesian> points{{0, 0}, {1, 0}, {Real(0.5), Real(1e-6)}};
	std::array<int, 3> order{0, 1, 2};
	do {
		std::vector<Point2Cartesian> permuted{points[order[0]], points[order[1]], points[order[2]]};
		const auto mesh = Triangulation::ComputeDelaunay(permuted);
		REQUIRE(mesh.NumTriangles() == 1);
		CheckMesh(mesh, permuted);
	} while (std::next_permutation(order.begin(), order.end()));
	CheckMesh(Triangulation::ComputeDelaunay({{0, 0}, {1, 0}, {Real(0.5), 1}}), {{0, 0}, {1, 0}, {Real(0.5), 1}});
}

TEST_CASE("Delaunay hull coverage and neighbors survive scale translation and insertion order", "[triangulation-contract][Delaunay]") {
	const Real thin = std::ldexp(Real(1), -std::numeric_limits<Real>::digits / 3);
	const std::vector<std::vector<Point2Cartesian>> clouds{
		{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {3, thin}, {2, thin}, {1, thin}, {0, thin}, {Real(1.5), thin / 2}},
		{{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}},
		{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {1, 1}, {3, 1}, {2, 3}},
		{{0, 0}, {1, thin}, {2, 0}, {3, 2 * thin}, {4, thin}, {2, 3 * thin}}
	};
	std::mt19937 rng(697);
	for (const auto& cloud : clouds) {
		for (int exponent : {0, std::numeric_limits<Real>::max_exponent / 2, std::numeric_limits<Real>::min_exponent / 2}) {
			for (Real shift : {Real(0), std::ldexp(Real(1), std::numeric_limits<Real>::digits / 4)}) {
				auto reference = cloud;
				for (int permutation = 0; permutation < 3; ++permutation) {
					CAPTURE(exponent, shift, permutation, cloud.size());
					const auto transformed = Transform(reference, exponent, shift);
					CheckMesh(Triangulation::ComputeDelaunay(transformed), reference);
					std::shuffle(reference.begin(), reference.end(), rng);
				}
			}
		}
	}
}

TEST_CASE("Delaunay duplicates retain first indices and collinear input has no faces", "[triangulation-contract][Delaunay]") {
	const std::vector<Point2Cartesian> points{{0, 0}, {0, 0}, {1, 0}, {Real(0.5), 1}, {1, 0}};
	const auto mesh = Triangulation::ComputeDelaunay(points);
	REQUIRE(mesh.NumPoints() == 5);
	REQUIRE(mesh.NumTriangles() == 1);
	CHECK(std::set<int>(mesh.triangles[0].begin(), mesh.triangles[0].end()) == std::set<int>{0, 2, 3});
	for (int i = 0; i < mesh.NumPoints(); ++i) {
		CHECK(mesh.points[i].X() == points[i].X());
		CHECK(mesh.points[i].Y() == points[i].Y());
	}
	for (const auto& line : std::vector<std::vector<Point2Cartesian>>{
		{}, {{1, 1}}, {{0, 0}, {1, 1}}, {{1, 1}, {1, 1}, {1, 1}},
		{{0, 0}, {1, 1}, {2, 2}, {-1, -1}}}) {
		const auto empty = Triangulation::ComputeDelaunay(line);
		CHECK(empty.points.size() == line.size());
		CHECK(empty.triangles.empty());
		CHECK(empty.neighbors.empty());
		CHECK(empty.LocatePoint({0, 0}) == -1);
	}
}

TEST_CASE("Ear clipping uses local winding at translated and extreme-scale geometry", "[triangulation-contract][Triangulation]") {
	const std::vector<std::vector<Point2Cartesian>> polygons{
		{{0, 0}, {0, 1}, {1, 1}, {1, 0}},
		{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}},
		{{0, 0}, {0, 1}, {1, 0}}
	};
	for (const auto& polygon : polygons) {
		for (int exponent : {0, std::numeric_limits<Real>::max_exponent / 2, std::numeric_limits<Real>::min_exponent / 2}) {
			const Real shift = std::numeric_limits<Real>::digits >= 53 ? Real(1e9) : Real(1024);
			auto vertices = Transform(polygon, exponent, shift);
			for (bool reverse : {false, true}) {
				CAPTURE(polygon.size(), exponent, reverse);
				if (reverse) std::reverse(vertices.begin(), vertices.end());
				const auto clipped = Triangulation::EarClipTriangulationDetailed(Polygon2D(vertices));
				REQUIRE(clipped.IsSuccess());
				CHECK_FALSE(clipped.used_fan_fallback);
				CHECK(clipped.unclipped_vertices == 0);
				REQUIRE(clipped.triangles.size() == polygon.size() - 2);
				long double total = 0;
				for (const auto& tri : clipped.triangles) {
					auto normalized = [&](const Point2Cartesian& p) {
						return Point2Cartesian(std::ldexp(p.X(), -exponent) - shift, std::ldexp(p.Y(), -exponent) - shift);
					};
					const auto area = Cross(normalized(tri.Pnt1()), normalized(tri.Pnt2()), normalized(tri.Pnt3()));
					CHECK(area > 0);
					total += area;
				}
				const long double expected = polygon.size() == 6 ? 6 : polygon.size() == 4 ? 2 : 1;
				CHECK(total == expected);
			}
		}
	}
	const Polygon2D degenerate({{0, 0}, {1, 0}, {2, 0}, {3, 0}});
	const auto fallback = Triangulation::EarClipTriangulationDetailed(degenerate);
	CHECK(fallback.status == AlgorithmStatus::InvalidInput);
	CHECK(fallback.used_fan_fallback);
	CHECK(fallback.unclipped_vertices == 4);
	CHECK_FALSE(fallback.error_message.empty());
	CHECK_THROWS_AS(Triangulation::EarClipTriangulation(degenerate), GeometryError);
}

TEST_CASE("Barycentric location and interpolation do not use absolute area cutoffs", "[triangulation-contract][Delaunay]") {
	for (int exponent : {0, 3 * std::numeric_limits<Real>::max_exponent / 4, 3 * std::numeric_limits<Real>::min_exponent / 4}) {
		const Real side = exponent == 0 ? Real(1e-6) : std::ldexp(Real(1), exponent);
		CAPTURE(exponent, side);
		DelaunayTriangulation mesh;
		mesh.points = {{0, 0}, {side, 0}, {0, side}};
		mesh.triangles = {{0, 1, 2}};
		mesh.neighbors = {{-1, -1, -1}};
		const Point2Cartesian inside(side / 4, side / 4), outside(side * 2, side * 2);
		CHECK(mesh.LocatePoint(inside) == 0);
		CHECK(mesh.LocatePoint(outside) == -1);
		CHECK(mesh.Interpolate(inside, {1, 1, 1}, -77) == 1);
		CHECK(mesh.Interpolate(outside, {1, 1, 1}, -77) == -77);
		CHECK(mesh.Interpolate(inside, {2, 6, 10}, -77) == Catch::Approx(5));
		CHECK(mesh.Interpolate({side / 2, 0}, {2, 6, 10}, -77) == Catch::Approx(4));
		for (int i = 0; i < 3; ++i)
			CHECK(mesh.Interpolate(mesh.points[i], {2, 6, 10}, -77) == 2 + 4 * i);
		const Real maximum = std::numeric_limits<Real>::max();
		CHECK(mesh.Interpolate(inside, {maximum, maximum, maximum}, -77) == maximum);
	}
	const auto translated = Triangulation::ComputeDelaunay({{8192, 8192}, {8193, 8192}, {8192, 8193}});
	CHECK(translated.Interpolate({Real(8192.25), Real(8192.25)}, {2, 6, 10}, -77) == 5);
	const Real maximum = std::numeric_limits<Real>::max(), tiny = std::numeric_limits<Real>::denorm_min();
	const auto mixed = Triangulation::ComputeDelaunay({{-maximum, 0}, {maximum, 0}, {0, tiny}});
	REQUIRE(mixed.NumTriangles() == 1);
	CHECK(mixed.Interpolate({0, 0}, {1, 3, 99}, -77) == 2);
	const Real epsilon = std::numeric_limits<Real>::epsilon();
	const auto cancellation = Triangulation::ComputeDelaunay({{0, 0}, {1, 1 + epsilon}, {1 + epsilon, 1 + 2 * epsilon}});
	REQUIRE(cancellation.NumTriangles() == 1);
	CHECK(cancellation.Interpolate({Real(0.5), (1 + epsilon) / 2}, {1, 3, 99}, -77) == 2);
}

TEST_CASE("Planar predicates retain incircle signs across the finite exponent range", "[triangulation-contract][Predicates]") {
	for (int exponent : {0, std::numeric_limits<Real>::max_exponent / 2, std::numeric_limits<Real>::min_exponent / 2}) {
		const Real scale = std::ldexp(Real(1), exponent);
		const Point2Cartesian a(scale, 0), b(0, scale), c(-scale, 0), on(0, -scale);
		const Point2Cartesian inside(0, std::nextafter(-scale, Real(0)));
		const Point2Cartesian outside(0, std::nextafter(-scale, -std::numeric_limits<Real>::infinity()));
		CHECK(RobustPredicates::InCircle2D(a, b, c, on) == 0);
		CHECK(RobustPredicates::InCircle2D(a, b, c, inside) == 1);
		CHECK(RobustPredicates::InCircle2D(a, b, c, outside) == -1);
		CHECK(RobustPredicates::InCircle2D(a, c, b, inside) == -1);
	}
	const Real large = std::numeric_limits<Real>::max(), tiny = std::numeric_limits<Real>::denorm_min();
	CHECK(RobustPredicates::InCircle2D({large, 0}, {0, large}, {-large, 0}, {0, tiny}) == 1);
	CHECK(RobustPredicates::InCircle2D({large, 0}, {0, large}, {-large, 0}, {tiny, -large}) == -1);
	CHECK(RobustPredicates::InCircle2D({tiny, 0}, {0, tiny}, {-tiny, 0}, {0, 0}) == 1);
}

TEST_CASE("Triangulation invalid coordinates and interpolation sizes fail explicitly", "[triangulation-contract][Delaunay]") {
	const auto mesh = Triangulation::ComputeDelaunay({{0, 0}, {1, 0}, {0, 1}});
	CHECK_THROWS_AS(mesh.Interpolate({0, 0}, {1, 2}, -77), ArgumentError);
	for (Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(Triangulation::ComputeDelaunay({{invalid, 0}}), GeometryError);
		CHECK_THROWS_AS(mesh.LocatePoint({invalid, 0}), GeometryError);
		CHECK_THROWS_AS(mesh.Interpolate({0, 0}, {invalid, 1, 2}), GeometryError);
		CHECK_THROWS_AS(Triangulation::EarClipTriangulationDetailed(Polygon2D({{0, 0}, {1, 0}, {0, invalid}})), GeometryError);
		CHECK_THROWS_AS(RobustPredicates::InCircle2D({0, 0}, {1, 0}, {0, 1}, {invalid, 0}), GeometryError);
	}
	DelaunayTriangulation degenerate;
	degenerate.points = {{0, 0}, {1, 0}, {2, 0}};
	degenerate.triangles = {{0, 1, 2}};
	CHECK(degenerate.LocatePoint({3, 0}) == -1);
}
