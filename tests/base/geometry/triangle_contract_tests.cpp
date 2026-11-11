#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/Geometry2D/Geometry2DTriangle.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DSurfaces.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

using namespace MML;

static_assert(std::is_convertible_v<TriangleSurface3D*, IParametricSurface<3>*>);
static_assert(std::is_convertible_v<const TriangleSurface3D*, const IParametricSurface<3>*>);
static_assert(std::is_convertible_v<TriangleSurface3D*, IVectorFunctionNM<2, 3>*>);

namespace {
	constexpr Real TriangleRoundoff = 128 * Precision::Eps;
	constexpr std::array<std::array<int, 3>, 6> TriangleVertexOrders{{
		{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}
	}};

	void CheckTriangleValue(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		if (expected == 0)
			CHECK(IsWithinAbs(actual, 0, TriangleRoundoff));
		else
			CHECK(IsWithinAbs(actual / expected, 1, TriangleRoundoff));
	}

	Real VertexAngle(const Point2Cartesian& origin, const Point2Cartesian& first, const Point2Cartesian& second) {
		const Real x1 = first.X() - origin.X(), y1 = first.Y() - origin.Y();
		const Real x2 = second.X() - origin.X(), y2 = second.Y() - origin.Y();
		return std::atan2(std::abs(x1 * y2 - y1 * x2), x1 * x2 + y1 * y2);
	}

	void CheckCurrentTriangle(const Triangle2D& triangle, bool right, bool isosceles, bool equilateral) {
		const auto p1 = triangle.Pnt1(), p2 = triangle.Pnt2(), p3 = triangle.Pnt3();
		const Real a = std::hypot(p2.X() - p1.X(), p2.Y() - p1.Y());
		const Real b = std::hypot(p3.X() - p2.X(), p3.Y() - p2.Y());
		const Real c = std::hypot(p1.X() - p3.X(), p1.Y() - p3.Y());
		const Real perimeter = a + b + c;
		const Real angleA = VertexAngle(p3, p1, p2);
		const Real area = a * c * std::sin(VertexAngle(p1, p2, p3)) / 2;
		CheckTriangleValue(triangle.A(), a);
		CheckTriangleValue(triangle.B(), b);
		CheckTriangleValue(triangle.C(), c);
		CheckTriangleValue(triangle.Perimeter(), perimeter);
		CheckTriangleValue(triangle.Area(), area);
		CheckTriangleValue(triangle.AngleA(), angleA);
		CheckTriangleValue(triangle.AngleB(), VertexAngle(p1, p2, p3));
		CheckTriangleValue(triangle.AngleC(), VertexAngle(p2, p1, p3));
		CheckTriangleValue(triangle.InRadius(), 2 * area / perimeter);
		CheckTriangleValue(triangle.CircumRadius(), a / (2 * std::sin(angleA)));
		const auto center = triangle.Incenter();
		CheckTriangleValue(center.X(), (a * p3.X() + b * p1.X() + c * p2.X()) / perimeter);
		CheckTriangleValue(center.Y(), (a * p3.Y() + b * p1.Y() + c * p2.Y()) / perimeter);
		CHECK(triangle.IsRight() == right);
		CHECK(triangle.IsIsosceles() == isosceles);
		CHECK(triangle.IsEquilateral() == equilateral);
	}

	void CheckSurfaceConsumer(const TriangleSurface3D& surface) {
		const IParametricSurface<3>* pointer = &surface;
		const IParametricSurface<3>& reference = surface;
		const IVectorFunctionNM<2, 3>& vectorFunction = surface;
		std::unique_ptr<IParametricSurface<3>> owned = std::make_unique<TriangleSurface3D>(surface);
		CheckTriangleValue(pointer->getMinU(), 0);
		CheckTriangleValue(pointer->getMaxU(), 4);
		CheckTriangleValue(pointer->getMinW(2), 0);
		CheckTriangleValue(pointer->getMaxW(0), 0);
		CheckTriangleValue(pointer->getMaxW(2), 1);
		CheckTriangleValue(pointer->getMaxW(4), 0);
		const std::array<std::array<Real, 2>, 5> parameters{{
			{0, 0}, {4, 0}, {2, 1}, {1, Real(0.5)}, {3, Real(0.25)}
		}};
		for (const auto& parameter : parameters) {
			const VectorN<Real, 2> coordinates{parameter[0], parameter[1]};
			const VectorN<Real, 3> expected{5, parameter[0] - 2, parameter[1]};
			for (const auto& actual : {reference(parameter[0], parameter[1]),
				reference(coordinates), vectorFunction(coordinates), (*owned)(coordinates)})
				for (int component = 0; component < 3; ++component)
					CheckTriangleValue(actual[component], expected[component]);
		}
	}
}

TEST_CASE("Retained triangle vertex and coordinate references update every side-derived query",
	"[triangle-contract][Triangle2D]") {
	Triangle2D triangle(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
	auto& first = triangle.Pnt1();
	auto& second = triangle.Pnt2();
	auto& third = triangle.Pnt3();
	auto& firstX = first.X();
	CheckCurrentTriangle(triangle, true, false, false);
	firstX = -3;
	CheckCurrentTriangle(triangle, false, true, false);
	second = Point2Cartesian(6, 0);
	CheckCurrentTriangle(triangle, false, false, false);
	third.Y() = 8;
	CheckCurrentTriangle(triangle, false, false, false);
	first = Point2Cartesian(0, 0);
	second = Point2Cartesian(2, 0);
	third = Point2Cartesian(1, std::sqrt(Real(3)));
	CheckCurrentTriangle(triangle, false, true, true);
}

TEST_CASE("Triangle copies and moves cannot preserve side data predating a retained-reference write",
	"[triangle-contract][Triangle2D]") {
	Triangle2D triangle(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
	auto& second = triangle.Pnt2();
	CheckTriangleValue(triangle.A(), 3);
	second.X() = 6;
	Triangle2D copied(triangle);
	Triangle2D moved(std::move(triangle));
	Triangle2D assigned(Point2Cartesian(0, 0), Point2Cartesian(1, 0), Point2Cartesian(0, 1));
	assigned = copied;
	Triangle2D moveAssigned = assigned;
	moveAssigned = std::move(copied);
	for (const auto* current : {&moved, &assigned, &moveAssigned})
		CheckCurrentTriangle(*current, true, false, false);
}

TEST_CASE("Coordinate triangle areas retain skinny heights in every vertex order",
	"[triangle-contract][Triangle2D]") {
	struct Example {
		std::array<Point2Cartesian, 3> points;
		Real area;
	};
	const Real base = std::ldexp(Real(1), 27), height = std::ldexp(Real(1), -23);
	// Resolve the translated height by at least eight coordinate ULPs; keep the original double offset.
	const int shiftExponent = std::min(10, std::numeric_limits<Real>::digits - 27);
	const Real shiftX = std::ldexp(Real(1), 20), shiftY = std::ldexp(Real(1), shiftExponent);
	const Real spacing = std::nextafter(shiftY, std::numeric_limits<Real>::infinity()) - shiftY;
	REQUIRE(height >= 8 * spacing);
	REQUIRE((shiftY + height) - shiftY == height);
	REQUIRE(base * height / 2 == 8);
	const std::array<Example, 3> examples{{
		{{Point2Cartesian(0, 0), Point2Cartesian(Real(1e8), 0), Point2Cartesian(Real(5e7), Real(1e-7L))}, 5},
		{{Point2Cartesian(shiftX, shiftY), Point2Cartesian(shiftX + base, shiftY),
			Point2Cartesian(shiftX + base / 2, shiftY + height)}, 8},
		{{Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4)}, 6}
	}};
	for (const auto& example : examples) {
		for (std::size_t index = 0; index < TriangleVertexOrders.size(); ++index) {
			const auto& order = TriangleVertexOrders[index];
			const Triangle2D triangle(example.points[order[0]], example.points[order[1]], example.points[order[2]]);
			CheckTriangleValue(triangle.Area(), example.area);
			CheckTriangleValue(triangle.SignedArea(), index < 3 ? example.area : -example.area);
			CHECK(triangle.IsCounterClockwise() == (index < 3));
			CHECK(triangle.IsClockwise() == (index >= 3));
		}
	}
}

TEST_CASE("Coordinate triangles do not invent heights erased from their stored vertices",
	"[triangle-contract][Triangle2D]") {
	const Real base = std::ldexp(Real(1), 27), height = std::ldexp(Real(1), -23);
	// Match the original float fixture's coordinate spacing in every Real type.
	const Real shiftY = std::ldexp(Real(1), std::numeric_limits<Real>::digits - 14);
	const Real spacing = std::nextafter(shiftY, std::numeric_limits<Real>::infinity()) - shiftY;
	REQUIRE(height > 0);
	REQUIRE(height < spacing / 2);
	REQUIRE(shiftY + height == shiftY);
	REQUIRE(base * height / 2 == 8);
	const std::array<Point2Cartesian, 3> points{{
		Point2Cartesian(0, shiftY), Point2Cartesian(base, shiftY), Point2Cartesian(base / 2, shiftY + height)
	}};
	for (const auto& order : TriangleVertexOrders) {
		CAPTURE(order[0], order[1], order[2]);
		const Triangle2D triangle(points[order[0]], points[order[1]], points[order[2]]);
		CHECK(triangle.Area() == 0);
		CHECK(triangle.SignedArea() == 0);
		CHECK_FALSE(triangle.IsClockwise());
		CHECK_FALSE(triangle.IsCounterClockwise());
		CHECK(triangle.InRadius() == 0);
		CHECK_THROWS_AS(triangle.CircumRadius(), GeometryError);
	}
}

TEST_CASE("Coordinate triangle areas are zero for repeated and collinear vertices",
	"[triangle-contract][Triangle2D]") {
	const std::array<std::array<Point2Cartesian, 3>, 4> examples{{
		{Point2Cartesian(0, 0), Point2Cartesian(1, 1), Point2Cartesian(2, 2)},
		{Point2Cartesian(4, 7), Point2Cartesian(4, 7), Point2Cartesian(-2, 5)},
		{Point2Cartesian(3, 4), Point2Cartesian(3, 4), Point2Cartesian(3, 4)},
		{Point2Cartesian(Real(-1e8), 0), Point2Cartesian(0, 0), Point2Cartesian(Real(1e8), 0)}
	}};
	for (const auto& points : examples) {
		const Triangle2D triangle(points[0], points[1], points[2]);
		CheckTriangleValue(triangle.Area(), 0);
		CheckTriangleValue(triangle.SignedArea(), 0);
		CHECK_FALSE(triangle.IsClockwise());
		CHECK_FALSE(triangle.IsCounterClockwise());
	}
}

TEST_CASE("Skinny coordinate triangle area propagates into radii and incenter",
	"[triangle-contract][Triangle2D]") {
	const Triangle2D triangle(Point2Cartesian(0, 0), Point2Cartesian(Real(1e8), 0),
		Point2Cartesian(Real(5e7), Real(1e-7L)));
	CheckTriangleValue(triangle.InRadius(), Real(5e-8L));
	CHECK_NOTHROW(CheckTriangleValue(triangle.CircumRadius(), Real(1.25e22L)));
	CheckTriangleValue(triangle.Incenter().X(), Real(5e7));
	CheckTriangleValue(triangle.Incenter().Y(), Real(5e-8L));
}

TEST_CASE("Triangle surfaces are usable through parametric and vector-function consumers",
	"[triangle-contract][TriangleSurface3D]") {
	CHECK((std::is_convertible_v<TriangleSurface3D*, IParametricSurface<3>*>));
	CHECK((std::has_virtual_destructor_v<IParametricSurface<3>>));
	const std::array<Pnt3Cart, 3> points{{{5, -2, 0}, {5, 2, 0}, {5, 0, 1}}};
	for (int first = 0; first < 3; ++first) {
		const TriangleSurface3D surface(points[first], points[(first + 1) % 3], points[(first + 2) % 3]);
		CheckSurfaceConsumer(surface);
	}
}
