#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/Geometry2D/Geometry2DLines.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DLines.h>
#endif

#include <array>
#include <cmath>
#include <limits>

using namespace MML;

namespace {
	constexpr Real LineRoundoff = 128 * Precision::Eps;

	void CheckScaled(Real actual, Real expected, Real scale = 1) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual / scale, Catch::Matchers::WithinAbs(expected / scale, LineRoundoff));
	}

	void CheckPoint(const Point2Cartesian& actual, const Point2Cartesian& expected, Real scale = 1) {
		CheckScaled(actual.X(), expected.X(), scale);
		CheckScaled(actual.Y(), expected.Y(), scale);
	}

	void CheckPoint(const Pnt3Cart& actual, const Pnt3Cart& expected, Real scale = 1) {
		CheckScaled(actual.X(), expected.X(), scale);
		CheckScaled(actual.Y(), expected.Y(), scale);
		CheckScaled(actual.Z(), expected.Z(), scale);
	}

	Pnt3Cart CyclePoint(Real x, Real y, Real z, int cycle) {
		if (cycle == 0) return Pnt3Cart(x, y, z);
		if (cycle == 1) return Pnt3Cart(z, x, y);
		return Pnt3Cart(y, z, x);
	}

	Vec3Cart CycleDirection(Real x, Real y, int cycle) {
		const auto point = CyclePoint(x, y, 0, cycle);
		return Vec3Cart(point.X(), point.Y(), point.Z());
	}
}

TEST_CASE("Near-parallel 3D lines have finite analytic closest pairs",
	"[line-segment-contract][Line3D]") {
	for (Real slope : {TOL(1e-9, 1e-3), -TOL(1e-9, 1e-3), TOL(1e-8, 1e-2), TOL(1e-6, 1e-1)}) {
		for (Real scale : {Real(1e-3), Real(1), Real(1e3)}) {
			for (Real gap : {Real(0), 2 * scale}) {
				for (int cycle = 0; cycle < 3; ++cycle) {
					for (Real sign : {Real(-1), Real(1)}) {
						CAPTURE(slope, scale, gap, cycle, sign);
						const Line3D first(CyclePoint(0, 0, 0, cycle), CycleDirection(1, 0, cycle));
						const Line3D second(CyclePoint(0, scale, gap, cycle),
							CycleDirection(sign, sign * slope, cycle));
						const Real x = -scale / slope;
						const Real t2 = x * std::hypot(Real(1), slope) / sign;
						const auto result = first.Intersection(second);
						CHECK(result.type == (gap == 0 ? LineIntersectionType3D::Point : LineIntersectionType3D::Skew));
						CheckScaled(result.t1, x, std::abs(x));
						CheckScaled(result.t2, t2, std::abs(t2));
						CheckPoint(result.point1, CyclePoint(x, 0, 0, cycle), std::abs(x));
						CheckPoint(result.point2, CyclePoint(x, 0, gap, cycle), std::abs(x));
						CheckScaled(result.distance, gap, scale);
						CheckScaled(first.Dist(second), gap, scale);

						Real distance = -1;
						Pnt3Cart point1, point2;
						CHECK(first.Dist(second, distance, point1, point2));
						CheckScaled(distance, gap, scale);
						CheckPoint(point1, result.point1, std::abs(x));
						CheckPoint(point2, result.point2, std::abs(x));
						const Vec3Cart connector(point1, point2);
						CheckScaled(connector * first.Direction(), 0, std::abs(x));
						CheckScaled(connector * second.Direction(), 0, std::abs(x));

						const auto reversed = second.Intersection(first);
						CHECK(reversed.type == result.type);
						CheckScaled(reversed.distance, gap, scale);
						CheckScaled(reversed.t1, result.t2, std::abs(t2));
						CheckScaled(reversed.t2, result.t1, std::abs(x));
						CheckPoint(reversed.point1, result.point2, std::abs(x));
						CheckPoint(reversed.point2, result.point1, std::abs(x));
						Pnt3Cart intersection(11, 12, 13);
						CHECK(first.Intersects(second, intersection) == (gap == 0));
						if (gap == 0) CheckPoint(intersection, result.point1, std::abs(x));
						else CheckPoint(intersection, Pnt3Cart(11, 12, 13));
					}
				}
			}
		}
	}
}

TEST_CASE("Parallel and coincident 3D results contain a representative closest pair",
	"[line-segment-contract][Line3D]") {
	for (Real sign : {Real(-1), Real(1)}) {
		for (Real separation : {Real(0), Real(1)}) {
			const Line3D first(Pnt3Cart(2, 3, 4), Vec3Cart(1, 0, 0));
			const Line3D second(Pnt3Cart(7, 3 + 3 * separation, 4 + 4 * separation), Vec3Cart(sign, 0, 0));
			const auto result = first.Intersection(second);
			CHECK(result.type == (separation == 0 ? LineIntersectionType3D::Coincident : LineIntersectionType3D::Parallel));
			CheckPoint(result.point1, Pnt3Cart(7, 3, 4));
			CheckPoint(result.point2, second.StartPoint());
			CheckPoint(first(result.t1), result.point1);
			CheckPoint(second(result.t2), result.point2);
			CheckScaled(result.t1, 5);
			CheckScaled(result.t2, 0);
			CheckScaled(result.distance, 5 * separation);
			CheckScaled(first.Dist(second), result.distance);
			CheckScaled(second.Dist(first), result.distance);

			Real distance = -1;
			Pnt3Cart point1(11, 12, 13), point2(14, 15, 16);
			CHECK_FALSE(first.Dist(second, distance, point1, point2));
			CheckScaled(distance, result.distance);
			CheckPoint(point1, Pnt3Cart(11, 12, 13));
			CheckPoint(point2, Pnt3Cart(14, 15, 16));
			Pnt3Cart intersection(11, 12, 13);
			CHECK(first.Intersects(second, intersection) == (separation == 0));
			CheckPoint(intersection, separation == 0 ? first.StartPoint() : Pnt3Cart(11, 12, 13));
		}
	}
}

TEST_CASE("All 3D line result surfaces follow the public parallel tolerance",
	"[line-segment-contract][Line3D]") {
	const Line3D first(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
	for (Real factor : {Real(0.5), Real(2)}) {
		for (Real sign : {Real(-1), Real(1)}) {
			const Real slope = factor * Precision::Line3DIsParallelTolerance;
			const Line3D second(Pnt3Cart(7, 1, 0), Vec3Cart(sign, sign * slope, 0));
			const bool parallel = factor < 1;
			CHECK(first.IsParallel(second) == parallel);
			const auto result = first.Intersection(second);
			CHECK(result.type == (parallel ? LineIntersectionType3D::Parallel : LineIntersectionType3D::Point));
			CheckScaled(result.distance, parallel ? Real(1) : Real(0));
			CheckScaled(first.Dist(second), result.distance);
			Real distance = -1;
			Pnt3Cart point1(11, 12, 13), point2(14, 15, 16);
			CHECK(first.Dist(second, distance, point1, point2) == !parallel);
			CheckScaled(distance, result.distance);
			if (parallel) {
				CheckPoint(point1, Pnt3Cart(11, 12, 13));
				CheckPoint(point2, Pnt3Cart(14, 15, 16));
				CheckPoint(result.point1, Pnt3Cart(7, 0, 0));
				CheckPoint(result.point2, second.StartPoint());
			}
		}
	}
}

TEST_CASE("3D line solves explicitly reject invalid or unrepresentable pair results",
	"[line-segment-contract][Line3D]") {
	const Line3D first(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
	for (Real invalid : {Real(0), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity()}) {
		Line3D second(Pnt3Cart(0, 1, 0), Vec3Cart(0, 1, 0));
		second.Direction() = Vec3Cart(invalid, 0, 0);
		CHECK_THROWS_AS(first.Intersection(second), GeometryError);
		CHECK_THROWS_AS(first.Dist(second), GeometryError);
		Real distance = -1;
		Pnt3Cart point1, point2;
		CHECK_THROWS_AS(first.Dist(second, distance, point1, point2), GeometryError);
	}
	const Line3D farLine(Pnt3Cart(0, std::numeric_limits<Real>::max() / 2, 0), Vec3Cart(1, TOL(1e-9, 1e-3), 0));
	CheckScaled(first.Dist(farLine), 0);
	CHECK_THROWS_AS(first.Intersection(farLine), GeometryError);
	const Real coordinate = std::numeric_limits<Real>::max() * Real(0.75);
	const Line3D diagonal(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 1));
	const Line3D farParallel(Pnt3Cart(coordinate, coordinate, coordinate), Vec3Cart(1, 1, 1));
	Real distance = -1;
	Pnt3Cart point1(11, 12, 13), point2(14, 15, 16);
	CHECK_FALSE(diagonal.Dist(farParallel, distance, point1, point2));
	CheckScaled(distance, 0, coordinate);
	CheckPoint(point1, Pnt3Cart(11, 12, 13));
	CheckPoint(point2, Pnt3Cart(14, 15, 16));
	CHECK_THROWS_AS(diagonal.Intersection(farParallel), GeometryError);
}

TEST_CASE("Point-segment intersections test actual incidence in either operand order",
	"[line-segment-contract][SegmentLine2D]") {
	const std::array<Point2Cartesian, 7> points{{
		{0, 0}, {1, 1}, {2, 2}, {1, 0}, {0, 1}, {3, 3}, {-1, -1}
	}};
	for (bool reverse : {false, true}) {
		const SegmentLine2D diagonal(reverse ? Point2Cartesian(2, 2) : Point2Cartesian(0, 0),
			reverse ? Point2Cartesian(0, 0) : Point2Cartesian(2, 2));
		for (std::size_t index = 0; index < points.size(); ++index) {
			const auto& point = points[index];
			const SegmentLine2D singleton(point, point);
			const bool intersects = index < 3;
			for (bool swap : {false, true}) {
				const auto& first = swap ? singleton : diagonal;
				const auto& second = swap ? diagonal : singleton;
				const auto result = first.Intersection(second);
				CHECK(result.type == (intersects ? SegmentIntersectionType::Point : SegmentIntersectionType::None));
				Point2Cartesian output(11, 12);
				CHECK(first.Intersects(second, &output) == intersects);
				CheckPoint(output, intersects ? point : Point2Cartesian(11, 12));
				if (intersects) CheckPoint(result.point, point);
			}
		}
	}
}

TEST_CASE("Two point-segments use geometric incidence tolerance",
	"[line-segment-contract][SegmentLine2D]") {
	const Point2Cartesian point(1, 2);
	const SegmentLine2D singleton(point, point);
	for (Real offset : {Real(0), Precision::GeometryEpsilon / 2, 2 * Precision::GeometryEpsilon, Real(1)}) {
		const Point2Cartesian otherPoint(1 + offset, 2);
		const SegmentLine2D other(otherPoint, otherPoint);
		const bool intersects = offset <= Precision::GeometryEpsilon;
		CHECK(singleton.Intersection(other).HasIntersection() == intersects);
		CHECK(other.Intersection(singleton).HasIntersection() == intersects);
	}
}

TEST_CASE("Nonzero 2D segments retain projection and distance across scales and reversal",
	"[line-segment-contract][SegmentLine2D]") {
	for (Real scale : {Real(1e-6), Real(1e-12), TOL(1e-100, 1e-15), TOL(1e-200, 1e-30), Real(1), TOL(1e100, 1e15), TOL(1e200, 1e30)}) {
		const Point2Cartesian start(0, 0), end(3 * scale, 4 * scale);
		for (bool reverse : {false, true}) {
			const SegmentLine2D segment(reverse ? end : start, reverse ? start : end);
			const Point2Cartesian query(-Real(1.25) * scale, Real(2.5) * scale);
			CheckPoint(segment.ClosestPoint(query), Point2Cartesian(Real(0.75) * scale, scale), scale);
			CheckScaled(segment.DistanceToPoint(query), Real(2.5) * scale, scale);
			CheckPoint(segment.ClosestPoint(Point2Cartesian(-3 * scale, -4 * scale)), start, scale);
			CheckScaled(segment.DistanceToPoint(Point2Cartesian(-3 * scale, -4 * scale)), 5 * scale, scale);
			CheckPoint(segment.ClosestPoint(Point2Cartesian(6 * scale, 8 * scale)), end, scale);
			CheckScaled(segment.DistanceToPoint(Point2Cartesian(6 * scale, 8 * scale)), 5 * scale, scale);
			CheckPoint(segment.ClosestPoint(start), start, scale);
			CheckPoint(segment.ClosestPoint(end), end, scale);
		}
	}
	const Real tiny = std::numeric_limits<Real>::min();
	const SegmentLine2D horizontal(Point2Cartesian(0, 0), Point2Cartesian(tiny, 0));
	CheckPoint(horizontal.ClosestPoint(Point2Cartesian(tiny / 2, 1)), Point2Cartesian(tiny / 2, 0), tiny);
	CheckScaled(horizontal.DistanceToPoint(Point2Cartesian(tiny / 2, 1)), 1);
}

TEST_CASE("2D segment degeneracy and ordinary intersection classifications remain usable",
	"[line-segment-contract][SegmentLine2D]") {
	const SegmentLine2D singleton(Point2Cartesian(2, 3), Point2Cartesian(2, 3));
	CheckPoint(singleton.ClosestPoint(Point2Cartesian(5, 7)), Point2Cartesian(2, 3));
	CheckScaled(singleton.DistanceToPoint(Point2Cartesian(5, 7)), 5);
	const SegmentLine2D horizontal(Point2Cartesian(0, 0), Point2Cartesian(4, 0));
	CHECK(horizontal.Intersection(SegmentLine2D(Point2Cartesian(2, -1), Point2Cartesian(2, 1))).IsPointIntersection());
	CHECK(horizontal.Intersection(SegmentLine2D(Point2Cartesian(4, 0), Point2Cartesian(5, 1))).IsPointIntersection());
	const auto overlap = horizontal.Intersection(SegmentLine2D(Point2Cartesian(3, 0), Point2Cartesian(1, 0)));
	CHECK(overlap.type == SegmentIntersectionType::Overlap);
	CheckPoint(overlap.overlapStart, Point2Cartesian(1, 0));
	CheckPoint(overlap.overlapEnd, Point2Cartesian(3, 0));
	CHECK(horizontal.Intersection(SegmentLine2D(Point2Cartesian(5, 0), Point2Cartesian(6, 0))).type
		== SegmentIntersectionType::Collinear);
	CHECK_FALSE(horizontal.Intersects(SegmentLine2D(Point2Cartesian(0, 1), Point2Cartesian(4, 1))));
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(horizontal.ClosestPoint(Point2Cartesian(invalid, 0)), GeometryError);
		CHECK_THROWS_AS(horizontal.DistanceToPoint(Point2Cartesian(0, invalid)), GeometryError);
	}
}
