#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/CompGeometry/Circles.h>
#include <mml/algorithms/CompGeometry/VoronoiDiagram.h>
#endif
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>

using namespace MML;
using namespace MML::CompGeometry;

namespace {
	// Normal inputs whose squares under/overflow Real, but whose coordinates and
	// skinny-triangle circumcenters remain finite in every supported precision.
	constexpr int scaleExponent = 3 * std::numeric_limits<Real>::max_exponent / 4;
	const Real largeOffset = std::scalbn(Real{1}, std::numeric_limits<Real>::digits - 4);
	struct OracleCircle {
		long double x = 0, y = 0, squaredRadius = std::numeric_limits<long double>::infinity();
	};

	OracleCircle BruteForceUnitCircle(const std::vector<Point2Cartesian>& points) {
		OracleCircle best;
		const auto consider = [&](long double x, long double y, long double radiusSquared) {
			for (const auto& point : points) {
				const long double dx = point.X() - x, dy = point.Y() - y;
				if (dx * dx + dy * dy > radiusSquared + Precision::RoundoffZeroThreshold)
					return;
			}
			if (radiusSquared < best.squaredRadius)
				best = {x, y, radiusSquared};
		};
		for (std::size_t i = 0; i < points.size(); ++i) {
			consider(points[i].X(), points[i].Y(), 0);
			for (std::size_t j = i + 1; j < points.size(); ++j) {
				const long double x = (points[i].X() + points[j].X()) / 2.L;
				const long double y = (points[i].Y() + points[j].Y()) / 2.L;
				consider(x, y, (points[i].X() - x) * (points[i].X() - x) + (points[i].Y() - y) * (points[i].Y() - y));
				for (std::size_t k = j + 1; k < points.size(); ++k) {
					const long double ax = points[i].X(), ay = points[i].Y();
					const long double bx = points[j].X(), by = points[j].Y();
					const long double cx = points[k].X(), cy = points[k].Y();
					const long double determinant = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
					if (determinant == 0) continue;
					const long double a2 = ax * ax + ay * ay, b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
					const long double ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / determinant;
					const long double uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / determinant;
					consider(ux, uy, (ux - ax) * (ux - ax) + (uy - ay) * (uy - ay));
				}
			}
		}
		return best;
	}

	std::vector<Point2Cartesian> OracleClippedCell(const std::vector<Point2Cartesian>& sites,
		int site, Real minX, Real minY, Real maxX, Real maxY) {
		struct Line { long double a, b, c; };
		std::vector<Line> lines{{-1, 0, -minX}, {1, 0, maxX}, {0, -1, -minY}, {0, 1, maxY}};
		const long double sx = sites[site].X(), sy = sites[site].Y();
		for (const auto& other : sites) {
			const long double ox = other.X(), oy = other.Y();
			if (ox == sx && oy == sy) continue;
			lines.push_back({2 * (ox - sx), 2 * (oy - sy), ox * ox + oy * oy - sx * sx - sy * sy});
		}
		std::vector<Point2Cartesian> result;
		// Independent oracle: intersect every pair of unit-scale halfplane boundaries,
		// then retain feasible intersections; no iterative polygon clipping.
		for (std::size_t i = 0; i < lines.size(); ++i) {
			for (std::size_t j = i + 1; j < lines.size(); ++j) {
				const auto& p = lines[i];
				const auto& q = lines[j];
				const long double determinant = p.a * q.b - q.a * p.b;
				if (determinant == 0) continue;
				const long double x = (p.c * q.b - q.c * p.b) / determinant;
				const long double y = (p.a * q.c - q.a * p.c) / determinant;
				const bool feasible = std::all_of(lines.begin(), lines.end(), [&](const Line& line) {
					const long double scale = std::abs(line.a * x) + std::abs(line.b * y) + std::abs(line.c);
					return line.a * x + line.b * y - line.c <= Precision::RoundoffZeroThreshold * scale;
				});
				if (!feasible) continue;
				const Point2Cartesian candidate{Real(x), Real(y)};
				if (std::none_of(result.begin(), result.end(), [&](const auto& point) {
					return std::hypot(point.X() - candidate.X(), point.Y() - candidate.Y()) <= Precision::EqualityTolerance;
				}))
					result.push_back(candidate);
			}
		}
		return result;
	}

	void CheckClippedCell(const VoronoiDiagram& diagram, const std::vector<Point2Cartesian>& unitSites,
		int site, Real scale, Real offset) {
		const auto expected = OracleClippedCell(unitSites, site, -2, -2, 4, 4);
		const auto actual = diagram.GetClippedCell(site, offset - 2 * scale, offset - 2 * scale,
			offset + 4 * scale, offset + 4 * scale);
		CAPTURE(site, scale, offset);
		REQUIRE(actual.NumVertices() == static_cast<int>(expected.size()));
		long double twiceArea = 0;
		for (int i = 0; i < actual.NumVertices(); ++i) {
			const Real x = (actual[i].X() - offset) / scale;
			const Real y = (actual[i].Y() - offset) / scale;
			REQUIRE(std::isfinite(x));
			REQUIRE(std::isfinite(y));
			CHECK(x >= -2);
			CHECK(x <= 4);
			CHECK(y >= -2);
			CHECK(y <= 4);
			CHECK(std::any_of(expected.begin(), expected.end(), [&](const auto& point) {
				return std::hypot(x - point.X(), y - point.Y()) <= Precision::EqualityTolerance;
			}));
			for (const auto& other : unitSites) {
				const long double ownX = static_cast<long double>(x) - unitSites[site].X();
				const long double ownY = static_cast<long double>(y) - unitSites[site].Y();
				const long double otherX = static_cast<long double>(x) - other.X();
				const long double otherY = static_cast<long double>(y) - other.Y();
				const long double ownSquared = ownX * ownX + ownY * ownY;
				const long double otherSquared = otherX * otherX + otherY * otherY;
				CHECK(ownSquared <= otherSquared + Precision::RoundoffZeroThreshold * (ownSquared + otherSquared));
			}
			const auto& next = actual[(i + 1) % actual.NumVertices()];
			CHECK_FALSE((actual[i].X() == next.X() && actual[i].Y() == next.Y()));
			twiceArea += static_cast<long double>(x) * ((next.Y() - offset) / scale)
				- static_cast<long double>(y) * ((next.X() - offset) / scale);
		}
		CHECK(twiceArea > 0);
	}
}

TEST_CASE("Voronoi nearest sites retain finite ordering across exponent extremes",
	"[CircleVoronoiContract][Voronoi][VoronoiQuery]")
{
	for (const Real scale : {std::scalbn(Real{1}, -scaleExponent), Real{1},
		std::scalbn(Real{1}, scaleExponent), std::numeric_limits<Real>::denorm_min() * 8}) {
		const auto diagram = Voronoi::ComputeVoronoi({{0, 0}, {scale, 0}, {scale, 0}});
		CHECK(diagram.FindNearestSite({2 * scale, 0}) == 1);
		CHECK(diagram.FindNearestSite({scale, 0}) == 1);
		CHECK(diagram.FindNearestSite({scale / 2, 0}) == 0);
	}
	const Real huge = std::numeric_limits<Real>::max();
	const auto wide = Voronoi::ComputeVoronoi({{-huge, -huge}, {-huge, 0}});
	CHECK(wide.FindNearestSite({huge, huge}) == 1);
	for (const Real offset : {largeOffset, -largeOffset}) {
		const auto translated = Voronoi::ComputeVoronoi({{offset, offset}, {offset + 2, offset}});
		CHECK(translated.FindNearestSite({offset + 3, offset + 1}) == 1);
		CHECK(translated.FindNearestSite({offset + 1, offset}) == 0);
	}
}

TEST_CASE("Voronoi clipped cells match independent halfplanes at tiny huge and translated scales",
	"[CircleVoronoiContract][Voronoi][VoronoiQuery]")
{
	for (const auto& unitSites : {std::vector<Point2Cartesian>{{2, 0}, {0, 0}, {1, 0}, {1, 0}},
		std::vector<Point2Cartesian>{{0, 0}, {2, 0}, {0, 2}},
		std::vector<Point2Cartesian>{{0, 0}, {-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 0}},
		std::vector<Point2Cartesian>{{1, 1}}}) {
		for (const auto& transform : std::array<std::array<Real, 2>, 7>{{
			{std::scalbn(Real{1}, -scaleExponent), 0}, {1, 0},
			{std::scalbn(Real{1}, scaleExponent), 0}, {1, largeOffset}, {1, -largeOffset},
			{std::scalbn(Real{1}, std::numeric_limits<Real>::max_exponent - 5),
				std::scalbn(Real{1}, std::numeric_limits<Real>::max_exponent - 1)},
			{std::numeric_limits<Real>::denorm_min() * 8, 0}}}) {
			const Real scale = transform[0], offset = transform[1];
			std::vector<Point2Cartesian> sites;
			for (const auto& point : unitSites)
				sites.emplace_back(offset + point.X() * scale, offset + point.Y() * scale);
			const auto diagram = Voronoi::ComputeVoronoi(sites);
			for (int site = 0; site < diagram.NumSites(); ++site)
				CheckClippedCell(diagram, unitSites, site, scale, offset);
		}
	}
}

TEST_CASE("Voronoi clipping does not need representable spans or padded boxes",
	"[CircleVoronoiContract][Voronoi][VoronoiQuery]")
{
	const Real huge = std::numeric_limits<Real>::max();
	const auto single = Voronoi::ComputeVoronoi({{huge, huge}});
	const auto box = single.GetClippedCell(0, -huge, -huge, huge, huge);
	REQUIRE(box.NumVertices() == 4);
	for (const auto& vertex : box.Vertices()) {
		CHECK(std::abs(vertex.X()) == huge);
		CHECK(std::abs(vertex.Y()) == huge);
	}
	const auto pair = Voronoi::ComputeVoronoi({{-huge, 0}, {huge, 0}});
	const auto left = pair.GetClippedCell(0, -huge, -huge, huge, huge);
	REQUIRE(left.NumVertices() == 4);
	int bisectorVertices = 0;
	for (const auto& vertex : left.Vertices()) {
		CHECK(std::isfinite(vertex.X()));
		CHECK(std::isfinite(vertex.Y()));
		CHECK(vertex.X() <= 0);
		CHECK(std::abs(vertex.Y()) == huge);
		if (vertex.X() == 0) ++bisectorVertices;
	}
	CHECK(bisectorVertices == 2);

	const Real tiny = std::numeric_limits<Real>::denorm_min() * 8;
	const auto mixed = Voronoi::ComputeVoronoi({{0, 0}, {2 * tiny, 0}});
	const auto narrow = mixed.GetClippedCell(0, 0, -1, huge, 1);
	REQUIRE(narrow.NumVertices() == 4);
	int tinyBoundary = 0;
	for (const auto& vertex : narrow.Vertices()) {
		CHECK((vertex.X() == 0 || vertex.X() == tiny));
		CHECK(std::abs(vertex.Y()) == 1);
		if (vertex.X() == tiny) ++tinyBoundary;
	}
	CHECK(tinyBoundary == 2);
}

TEST_CASE("Voronoi oblique clipped cells match exhaustive unit halfplane intersections",
	"[CircleVoronoiContract][Voronoi][VoronoiQuery]")
{
	const std::array<Point2Cartesian, 6> lattice{{{-1, -1}, {2, 0}, {0, 3}, {3, 2}, {1, 1}, {0, 1}}};
	for (int mask = 1; mask < 64; ++mask) {
		std::vector<Point2Cartesian> unitSites;
		for (int i = 0; i < 6; ++i)
			if ((mask & (1 << i)) != 0) unitSites.push_back(lattice[i]);
		for (const Real scale : {std::scalbn(Real{1}, -scaleExponent), Real{1}, std::scalbn(Real{1}, scaleExponent)}) {
			std::vector<Point2Cartesian> sites;
			for (const auto& point : unitSites)
				sites.emplace_back(point.X() * scale, point.Y() * scale);
			const auto diagram = Voronoi::ComputeVoronoi(sites);
			for (int site = 0; site < diagram.NumSites(); ++site)
				CheckClippedCell(diagram, unitSites, site, scale, 0);
		}
	}
}

TEST_CASE("Voronoi queries validate inputs and preserve empty duplicate and boundary contracts",
	"[CircleVoronoiContract][Voronoi][VoronoiQuery]")
{
	const auto empty = Voronoi::ComputeVoronoi({});
	CHECK(empty.FindNearestSite({0, 0}) == -1);
	CHECK(empty.GetClippedCell(0, -1, -1, 1, 1).NumVertices() == 0);
	const auto diagram = Voronoi::ComputeVoronoi({{0, 0}, {2, 0}, {0, 0}});
	CHECK(diagram.GetClippedCell(-1, -1, -1, 1, 1).NumVertices() == 0);
	CHECK(diagram.GetClippedCell(3, -1, -1, 1, 1).NumVertices() == 0);
	CHECK(diagram.GetClippedCell(0, 2, -1, 3, 1).NumVertices() == 0);
	CHECK(diagram.GetClippedCell(0, 1, -1, 2, 1).NumVertices() == 0);
	CHECK(diagram.GetClippedCell(0, 0, -1, 0, 1).NumVertices() == 0);
	CHECK_THROWS_AS(diagram.GetClippedCell(0, 1, -1, -1, 1), GeometryError);
	CHECK_THROWS_AS(diagram.GetClippedCell(0, -1, 1, 1, -1), GeometryError);
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(empty.FindNearestSite({invalid, 0}), GeometryError);
		CHECK_THROWS_AS(diagram.FindNearestSite({0, invalid}), GeometryError);
		CHECK_THROWS_AS(diagram.GetClippedCell(0, invalid, -1, 1, 1), GeometryError);
		CHECK_THROWS_AS(diagram.GetClippedCell(0, -1, invalid, 1, 1), GeometryError);
		CHECK_THROWS_AS(diagram.GetClippedCell(0, -1, -1, invalid, 1), GeometryError);
		CHECK_THROWS_AS(diagram.GetClippedCell(0, -1, -1, 1, invalid), GeometryError);
		CHECK_THROWS_AS(empty.GetClippedCell(0, invalid, -1, 1, 1), GeometryError);
		auto corrupt = diagram;
		corrupt.sites[0] = {invalid, 0};
		CHECK_THROWS_AS(corrupt.FindNearestSite({0, 0}), GeometryError);
		CHECK_THROWS_AS(corrupt.GetClippedCell(0, -1, -1, 1, 1), GeometryError);
	}
	const auto first = diagram.GetClippedCell(0, -1, -1, 1, 1);
	const auto duplicate = diagram.GetClippedCell(2, -1, -1, 1, 1);
	REQUIRE(first.NumVertices() == 4);
	REQUIRE(duplicate.NumVertices() == first.NumVertices());
	for (int i = 0; i < first.NumVertices(); ++i) {
		CHECK(first[i].X() == duplicate[i].X());
		CHECK(first[i].Y() == duplicate[i].Y());
	}
}

TEST_CASE("Minimum circles match an independent exhaustive support oracle on scaled lattice subsets",
	"[CircleVoronoiContract][MEC]")
{
	for (int mask = 1; mask < 512; ++mask) {
		const Real scale = std::scalbn(Real{1}, (mask % 3 - 1) * scaleExponent);
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale > 0);
		std::vector<Point2Cartesian> unit, scaled;
		for (int slot = 0; slot < 9; ++slot)
			if ((mask & (1 << slot)) != 0) {
				unit.emplace_back(Real(slot % 3 - 1), Real(slot / 3 - 1));
				scaled.emplace_back(unit.back().X() * scale, unit.back().Y() * scale);
			}
		const auto expected = BruteForceUnitCircle(unit);
		const auto actual = Circles::MinimumEnclosingCircle(scaled);
		CHECK(IsWithinAbs(actual.Center().X() / scale, Real(expected.x), Precision::EqualityTolerance));
		CHECK(IsWithinAbs(actual.Center().Y() / scale, Real(expected.y), Precision::EqualityTolerance));
		CHECK(IsWithinAbs(actual.Radius() / scale, Real(std::sqrt(expected.squaredRadius)), Precision::EqualityTolerance));
		for (const auto& point : scaled)
			CHECK(std::hypot((point.X() - actual.Center().X()) / scale, (point.Y() - actual.Center().Y()) / scale) <=
				(actual.Radius() / scale) * (Real{1} + Precision::RoundoffZeroThreshold));
	}
}

TEST_CASE("Minimum enclosing circles retain small and enormous noncollinear triangles",
	"[CircleVoronoiContract][MEC]")
{
	for (const Real scale : {Real{1e-6}, std::scalbn(Real{1}, -scaleExponent), Real{1},
		std::scalbn(Real{1}, scaleExponent)}) {
		const std::array<Point2Cartesian, 3> input{{{0, 0}, {2 * scale, 0}, {scale, std::sqrt(Real{3}) * scale}}};
		std::array<int, 3> order{0, 1, 2};
		do {
			const std::vector<Point2Cartesian> points{input[order[0]], input[order[1]], input[order[2]]};
			const auto circle = Circles::MinimumEnclosingCircle(points);
			CHECK(IsWithinAbs(circle.Radius() / scale, Real{2} / std::sqrt(Real{3}), Precision::EqualityTolerance));
			CHECK(IsWithinAbs(circle.Center().X() / scale, Real{1}, Precision::EqualityTolerance));
			CHECK(IsWithinAbs(circle.Center().Y() / scale, Real{1} / std::sqrt(Real{3}), Precision::EqualityTolerance));
			for (const auto& point : points)
				CHECK(std::hypot((point.X() - circle.Center().X()) / scale,
					(point.Y() - circle.Center().Y()) / scale) <=
					(circle.Radius() / scale) * (Real{1} + Precision::RoundoffZeroThreshold));
		} while (std::next_permutation(order.begin(), order.end()));
	}
}

TEST_CASE("Circle construction preserves translated and diametral controls", "[CircleVoronoiContract][MEC]")
{
	// Keep unit separations representable even at a large translated origin.
	for (const Real offset : {Real{0}, largeOffset, -largeOffset}) {
		const auto circle = Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{
			{offset, offset}, {offset + 2, offset}, {offset, offset + 2}});
		CHECK(circle.Center().X() == offset + 1);
		CHECK(circle.Center().Y() == offset + 1);
		CHECK(IsWithinAbs(circle.Radius(), std::sqrt(Real{2}), Precision::EqualityTolerance));
	}
	const Real huge = std::numeric_limits<Real>::max();
	const auto diameter = Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{{-huge, 0}, {huge, 0}});
	CHECK(diameter.Center().X() == 0);
	CHECK(diameter.Radius() == huge);
	const auto same = Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{{huge, huge}, {huge, huge}});
	CHECK(same.Center().X() == huge);
	CHECK(same.Center().Y() == huge);
	CHECK(same.Radius() == 0);
	CHECK_THROWS_AS(Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{{-huge, -huge}, {huge, huge}}), GeometryError);
}

TEST_CASE("Voronoi circumcenters and boundary directions use local scale safe arithmetic",
	"[CircleVoronoiContract][Voronoi]")
{
	for (const Real scale : {std::scalbn(Real{1}, -scaleExponent), Real{1}, std::scalbn(Real{1}, scaleExponent)}) {
		const auto diagram = Voronoi::ComputeVoronoi({{0, 0}, {2 * scale, 0}, {0, 2 * scale}});
		REQUIRE(diagram.NumVertices() == 1);
		REQUIRE(diagram.NumEdges() == 3);
		CHECK(IsWithinAbs(diagram.vertices[0].X() / scale, Real{1}, Precision::EqualityTolerance));
		CHECK(IsWithinAbs(diagram.vertices[0].Y() / scale, Real{1}, Precision::EqualityTolerance));
		for (const auto& edge : diagram.edges) {
			CHECK(std::isfinite(edge.direction.X()));
			CHECK(std::isfinite(edge.direction.Y()));
			CHECK(IsWithinAbs(std::hypot(edge.direction.X(), edge.direction.Y()), Real{1}, Precision::EqualityTolerance));
			const auto& first = diagram.sites[edge.site1];
			const auto& second = diagram.sites[edge.site2];
			CHECK(std::abs(((second.X() - first.X()) / scale) * edge.direction.X() +
				((second.Y() - first.Y()) / scale) * edge.direction.Y()) <= Precision::EqualityTolerance);
		}
	}
	for (const Real offset : {largeOffset, -largeOffset}) {
		const auto diagram = Voronoi::ComputeVoronoi({{offset, offset}, {offset + 2, offset}, {offset, offset + 2}});
		REQUIRE(diagram.NumVertices() == 1);
		CHECK(diagram.vertices[0].X() == offset + 1);
		CHECK(diagram.vertices[0].Y() == offset + 1);
	}
}

TEST_CASE("Collinear Voronoi topology contains every adjacent distinct site bisector",
	"[CircleVoronoiContract][Voronoi]")
{
	for (const Real scale : {std::scalbn(Real{1}, -scaleExponent), Real{1}, std::scalbn(Real{1}, scaleExponent)}) {
		const auto diagram = Voronoi::ComputeVoronoi({{2 * scale, 4 * scale}, {0, 0}, {scale, 2 * scale}, {scale, 2 * scale}});
		REQUIRE(diagram.NumSites() == 4);
		REQUIRE(diagram.NumVertices() == 0);
		REQUIRE(diagram.NumEdges() == 2);
		std::set<std::pair<int, int>> adjacent;
		for (const auto& edge : diagram.edges) {
			CHECK(edge.v1 == -1);
			CHECK(edge.v2 == -1);
			CHECK(edge.IsInfinite());
			CHECK(IsWithinAbs(std::hypot(edge.direction.X(), edge.direction.Y()), Real{1}, Precision::EqualityTolerance));
			CHECK(IsWithinAbs(edge.direction.X() + 2 * edge.direction.Y(), Real{0}, Precision::EqualityTolerance));
			adjacent.emplace(std::min(edge.site1, edge.site2), std::max(edge.site1, edge.site2));
		}
		CHECK(adjacent == std::set<std::pair<int, int>>{{0, 2}, {1, 2}});
		for (int site = 0; site < diagram.NumSites(); ++site)
			CHECK(diagram.IsCellInfinite(site));
	}
	const auto same = Voronoi::ComputeVoronoi({{1, 1}, {1, 1}});
	CHECK(same.NumEdges() == 0);
	CHECK(same.IsCellInfinite(0));
	CHECK(same.IsCellInfinite(1));
	const auto pair = Voronoi::ComputeVoronoi({{0, 0}, {2, 0}});
	REQUIRE(pair.NumEdges() == 1);
	CHECK(IsWithinAbs(std::hypot(pair.edges[0].direction.X(), pair.edges[0].direction.Y()), Real{1}, Precision::EqualityTolerance));
}

TEST_CASE("Circle and Voronoi entry points reject nonfinite coordinates before shortcuts",
	"[CircleVoronoiContract][MEC][Voronoi]")
{
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		for (const auto& points : {std::vector<Point2Cartesian>{{invalid, 0}},
			std::vector<Point2Cartesian>{{0, 0}, {0, invalid}},
			std::vector<Point2Cartesian>{{0, 0}, {1, 0}, {invalid, 1}}}) {
			CHECK_THROWS_AS(Circles::MinimumEnclosingCircle(points), GeometryError);
			CHECK_THROWS_AS(Voronoi::ComputeVoronoi(points), GeometryError);
		}
	}
	CHECK(Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{}).Radius() == 0);
	CHECK(Voronoi::ComputeVoronoi({}).NumEdges() == 0);
}

TEST_CASE("Voronoi degenerate topology and unrepresentable measurements are explicit",
	"[CircleVoronoiContract][MEC][Voronoi]")
{
	const std::vector<Point2Cartesian> sites{{0, 0}, {-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 0}};
	const auto duplicate = Voronoi::ComputeVoronoi(sites);
	REQUIRE_FALSE(duplicate.IsCellInfinite(0));
	CHECK_FALSE(duplicate.IsCellInfinite(5));
	CHECK_FALSE(duplicate.cells[0].empty());
	CHECK(duplicate.cells[0] == duplicate.cells[5]);
	for (const auto& edge : duplicate.edges)
		CHECK(edge.site1 != 5);
	for (const auto& edge : duplicate.edges)
		CHECK(edge.site2 != 5);
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	const std::vector<Point2Cartesian> skinny{{-1, 0}, {1, 0}, {0, tiny}};
	CHECK_THROWS_AS(Voronoi::ComputeVoronoi(skinny), GeometryError);
	const auto enclosing = Circles::MinimumEnclosingCircle(skinny);
	CHECK(enclosing.Center().X() == 0);
	CHECK(enclosing.Center().Y() == 0);
	CHECK(enclosing.Radius() == 1);
	for (const Real scale : {std::scalbn(Real{1}, -scaleExponent), Real{1}, std::scalbn(Real{1}, scaleExponent)}) {
		const Real height = Real{1e-8};
		const auto diagram = Voronoi::ComputeVoronoi({{0, 0}, {2 * scale, 0}, {scale, height * scale}});
		REQUIRE(diagram.NumVertices() == 1);
		const Real expectedY = (height * height - Real{1}) / (Real{2} * height);
		CHECK(IsWithinAbs(diagram.vertices[0].X() / scale, Real{1}, Precision::EqualityTolerance));
		CHECK(IsWithinAbs(diagram.vertices[0].Y() / scale / expectedY, Real{1}, Precision::EqualityTolerance));
	}
}
