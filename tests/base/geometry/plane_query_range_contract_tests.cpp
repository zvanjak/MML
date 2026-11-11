#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/Geometry3D/Geometry3DPlane.h>
#endif
#include <array>
#include <limits>

using namespace MML;
namespace {
	const Real Max = std::numeric_limits<Real>::max();
	const Real Min = std::numeric_limits<Real>::denorm_min();
	const Real LargeCoefficient = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 16);
	Plane3D Coefficients(Real a, Real b, Real c, Real d) {
		auto plane = Plane3D::XYPlane();
		plane.A() = a; plane.B() = b; plane.C() = c; plane.D() = d;
		return plane;
	}
	void Relative(Real value, Real expected) {
		REQUIRE(std::isfinite(value));
		CHECK(std::abs(value / expected - 1) <= 64 * Constants::Eps);
	}
	void CheckQueryPoint(const Pnt3Cart& value, const Pnt3Cart& expected, Real scale = 1) {
		const std::array<Real, 3> actual{value.X(), value.Y(), value.Z()};
		const std::array<Real, 3> reference{expected.X(), expected.Y(), expected.Z()};
		for (int i = 0; i < 3; ++i) {
			REQUIRE(std::isfinite(actual[i]));
			CHECK(std::abs(actual[i] / scale - reference[i] / scale) <= 128 * Constants::Eps);
		}
	}
}

TEST_CASE("Plane queries retain finite results after overflowing dot sums", "[plane-query-range][Plane3D]") {
	const Plane3D plane(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 1));
	const Pnt3Cart p(Max, Max, -Max / 2);
	Relative(plane.DistToPoint(p), (Max / 2) * std::sqrt(Real(3)));
	CheckQueryPoint(plane.ProjectionToPlane(p), Pnt3Cart(Max / 2, Max / 2, -Max), Max);
	CHECK_THROWS_AS(plane.DistToPoint(Pnt3Cart(Max, Max, Max)), GeometryError);
	CheckQueryPoint(plane.ProjectionToPlane(Pnt3Cart(Max, Max, Max)), Pnt3Cart(0, 0, 0), Max);
}

TEST_CASE("Plane queries support finite nonunit and mutable coefficient scales", "[plane-query-range][Plane3D]") {
	for (Real scale : {Min, std::numeric_limits<Real>::min() * 16, Real(1e-8), Real(1), LargeCoefficient, -LargeCoefficient}) {
		CAPTURE(scale);
		const auto plane = Coefficients(scale, scale, scale, -scale);
		Relative(plane.DistToPoint(Pnt3Cart(4, Real(0.5), -1)), Real(2.5) / std::sqrt(Real(3)));
		CheckQueryPoint(plane.ProjectionToPlane(Pnt3Cart(4, Real(0.5), -1)),
			Pnt3Cart(4 - Real(2.5) / 3, Real(0.5) - Real(2.5) / 3, -1 - Real(2.5) / 3));
		CheckQueryPoint(plane.getPointOnPlane(), Pnt3Cart(Real(1) / 3, Real(1) / 3, Real(1) / 3));
		const Real tiny = Min * 16;
		CheckQueryPoint(Coefficients(scale, 2 * scale, -scale, 0).ProjectionToPlane(Pnt3Cart(4 * tiny, tiny / 2, -tiny)),
			Pnt3Cart(3 * tiny, -Real(1.5) * tiny, 0), tiny);
	}
}

TEST_CASE("Plane residuals retain cancellation and the absolute membership contract", "[plane-query-range][Plane3D]") {
	const Real huge = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 32);
	const auto plane = Coefficients(1, 1, 1, -1);
	for (const auto p : {Pnt3Cart(huge, 1, -huge), Pnt3Cart(1, huge, -huge), Pnt3Cart(huge, -huge, 1)}) {
		CHECK(plane.IsPointOnPlane(p));
		CHECK(plane.DistToPoint(p) == 0);
		CheckQueryPoint(plane.ProjectionToPlane(p), p, huge);
	}
	const Real u = Constants::Eps;
	const auto products = Coefficients(1 + u, 1, 1, 0);
	Relative(products.DistToPoint(Pnt3Cart(1 + u, -(1 + 2 * u), 0)),
		u * u / std::hypot(1 + u, Real(1), Real(1)));
	const Real tolerance = Precision::Plane3DIsPointOnPlaneTolerance;
	const auto nonunit = Coefficients(2, 0, 0, 0);
	CHECK(nonunit.IsPointOnPlane(Pnt3Cart(tolerance / 4, 0, 0)));
	CHECK_FALSE(nonunit.IsPointOnPlane(Pnt3Cart(tolerance, 0, 0)));
	CHECK_FALSE(nonunit.IsPointOnPlane(Pnt3Cart(tolerance / 2, 0, 0)));
	CHECK_FALSE(nonunit.IsPointOnPlane(Pnt3Cart(0, 0, 0), 0));
}

TEST_CASE("Plane line intersections keep extended parameters until final coordinates", "[plane-query-range][Plane3D]") {
	const auto plane = Coefficients(1, 0, 0, -Max);
	const Line3D line(Pnt3Cart(-Max, 2, 3), Vec3Cart(1, 0, 0));
	Pnt3Cart out(7, 8, 9);
	REQUIRE(plane.IntersectionWithLine(line, out));
	CheckQueryPoint(out, Pnt3Cart(Max, 2, 3), Max);
	CheckQueryPoint(plane.ProjectionToPlane(line.StartPoint()), Pnt3Cart(Max, 2, 3), Max);
	for (Real scale : {Min, Real(1), Max}) {
		const auto scaled = Coefficients(scale, 0, 0, -scale);
		Line3D ray(Pnt3Cart(3, 2, 1), Vec3Cart(1, 0, 0));
		ray.Direction() = Vec3Cart(scale, 0, 0);
		REQUIRE(scaled.IntersectionWithLine(ray, out));
		CheckQueryPoint(out, Pnt3Cart(1, 2, 1));
	}
	out = Pnt3Cart(7, 8, 9);
	CHECK_FALSE(plane.IntersectionWithLine(Line3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0)), out));
	CheckQueryPoint(out, Pnt3Cart(7, 8, 9));
	const auto xy = Plane3D::XYPlane();
	Line3D on(Pnt3Cart(Max, Max, 0), Vec3Cart(1, 0, 0));
	on.Direction() = Vec3Cart(Max, Max, 0);
	CHECK(xy.IsLineOnPlane(on));
}

TEST_CASE("Plane plane intersections and distances ignore coefficient magnitude", "[plane-query-range][Plane3D]") {
	for (Real scale : {Min, Real(1), Max}) {
		const auto first = Coefficients(scale, 0, 0, -scale);
		const auto second = Coefficients(0, scale, 0, -scale);
		Line3D out(Pnt3Cart(7, 8, 9), Vec3Cart(0, 0, 1));
		REQUIRE(first.IntersectionWithPlane(second, out));
		CheckQueryPoint(out.StartPoint(), Pnt3Cart(1, 1, 0));
		Relative(std::abs(out.Direction().Z()), Real(1));
		CHECK(first.DistToPlane(second) == 0);
		CHECK(first.IsPerpendicularToPlane(second));
		Relative(first.AngleToPlane(second), Constants::PI / 2);
		Relative(first.AngleToLine(Line3D(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0))), Constants::PI / 2);
		CHECK(first.Normal().X() == scale);
		CHECK(first.IsParallelToPlane(Coefficients(-scale, 0, 0, 0)));
		Relative(first.DistToPlane(Coefficients(-scale, 0, 0, 0)), Real(1));
	}
	const auto first = Coefficients(1, 0, 0, -Max);
	const auto second = Coefficients(1, 1, 1, Max);
	Line3D line;
	REQUIRE(first.IntersectionWithPlane(second, line));
	CheckQueryPoint(line.StartPoint(), Pnt3Cart(Max, -Max, -Max), Max);
	const auto remote = Coefficients(Min, 0, 0, Max);
	CHECK(remote.DistToPlane(remote) == 0);
	CHECK_THROWS_AS(remote.DistToPlane(Coefficients(Min, 0, 0, -Max)), GeometryError);
}

TEST_CASE("Plane query outputs explicitly reject genuine range loss", "[plane-query-range][Plane3D]") {
	CHECK_THROWS_AS(Coefficients(1, 3, 0, 0).DistToPoint(Pnt3Cart(Min, 0, 0)), GeometryError);
	CHECK_THROWS_AS(Coefficients(1, 1, 0, 0).ProjectionToPlane(Pnt3Cart(Min, 0, 0)), GeometryError);
	CHECK_THROWS_AS(Coefficients(1, 2, 2, 0).ProjectionToPlane(Pnt3Cart(-Max, Max, Max)), GeometryError);
	Real x = 7, y = 8, z = 9;
	CHECK_THROWS_AS(Coefficients(Min, 0, 0, Max).getCoordAxisSegments(x, y, z), GeometryError);
	CHECK(x == 7); CHECK(y == 8); CHECK(z == 9);
	Plane3D::XYPlane().getCoordAxisSegments(x, y, z);
	CHECK(std::isinf(x)); CHECK(std::isinf(y)); CHECK(z == 0);
	const auto remote = Coefficients(1, 0, 0, -Max);
	Pnt3Cart point(7, 8, 9);
	CHECK_THROWS_AS(remote.IntersectionWithLine(Line3D(Pnt3Cart(0, 0, 0), Vec3Cart(1, 2, 0)), point), GeometryError);
	CheckQueryPoint(point, Pnt3Cart(7, 8, 9));
	Line3D line(Pnt3Cart(7, 8, 9), Vec3Cart(0, 0, 1));
	CHECK_THROWS_AS(remote.IntersectionWithPlane(Coefficients(1, 1, 0, Max), line), GeometryError);
	CheckQueryPoint(line.StartPoint(), Pnt3Cart(7, 8, 9));
}

TEST_CASE("Plane queries reject invalid coefficients points directions and tolerances", "[plane-query-range][Plane3D]") {
	for (Real bad : {Real(0), Constants::PosInf, std::numeric_limits<Real>::quiet_NaN()}) {
		const auto plane = Coefficients(bad, 0, 0, 0);
		const Pnt3Cart p(1, 2, 3);
		Pnt3Cart out(7, 8, 9);
		Line3D line(p, Vec3Cart(1, 0, 0)), result(p, Vec3Cart(0, 1, 0));
		CHECK_THROWS_AS(plane.DistToPoint(p), GeometryError);
		CHECK_THROWS_AS(plane.IsPointOnPlane(p), GeometryError);
		CHECK_THROWS_AS(plane.ProjectionToPlane(p), GeometryError);
		CHECK_THROWS_AS(plane.getPointOnPlane(), GeometryError);
		CHECK_THROWS_AS(plane.IntersectionWithLine(line, out), GeometryError);
		CHECK_THROWS_AS(plane.IntersectionWithPlane(Plane3D::XYPlane(), result), GeometryError);
		CHECK_THROWS_AS(plane.DistToPlane(Plane3D::XYPlane()), GeometryError);
		CheckQueryPoint(out, Pnt3Cart(7, 8, 9));
	}
	auto plane = Plane3D::XYPlane();
	plane.D() = Constants::PosInf;
	CHECK_THROWS_AS(plane.DistToPoint(Pnt3Cart(0, 0, 0)), GeometryError);
	plane = Plane3D::XYPlane();
	CHECK_THROWS_AS(plane.ProjectionToPlane(Pnt3Cart(Constants::PosInf, 0, 0)), GeometryError);
	Line3D line(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
	line.Direction() = Vec3Cart(0, 0, 0);
	Pnt3Cart out;
	CHECK_THROWS_AS(plane.IntersectionWithLine(line, out), GeometryError);
	for (Real bad : {Real(-1), Constants::PosInf, std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(plane.IsPointOnPlane(Pnt3Cart(0, 0, 0), bad), GeometryError);
		CHECK_THROWS_AS(plane.IntersectionWithLine(line, out, bad), GeometryError);
	}
}
