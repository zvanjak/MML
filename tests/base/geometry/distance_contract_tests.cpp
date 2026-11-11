#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/GeometrySpherical.h>
#endif

#include <array>
#include <cmath>
#include <limits>

using namespace MML;

namespace {
	constexpr Real DistanceRoundoff = 128 * Precision::Eps;
	constexpr Real TinyScale = TOL(1e-200, 1e-25);
	constexpr Real LargeScale = TOL(1e200, 1e25);
	constexpr Real FarRadius = TOL(1e300, 1e35);

	void CheckDistance(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		if (expected == 0) CHECK(actual == 0);
		else CHECK_THAT(actual / expected, Catch::Matchers::WithinAbs(1, DistanceRoundoff));
	}

	template<class Point>
	void CheckPair(const Point& first, const Point& second, Real expected) {
		REQUIRE(std::isfinite(expected));
		REQUIRE(expected >= 0);
		CheckDistance(first.Dist(second), expected);
		CheckDistance(second.Dist(first), expected);
		CheckDistance(first.Dist(first), 0);
		if (expected > 0) {
			CHECK_FALSE(first.IsEqualTo(second, expected / 2));
			CHECK(first.IsEqualTo(second, 2 * expected));
		}
	}
}

TEST_CASE("Cartesian point distances preserve representable large and tiny norms",
	"[distance-contract][Point2Cartesian][Point3Cartesian]") {
	for (Real scale : {TinyScale, Real(1), LargeScale}) {
		REQUIRE(scale > 0);
		CheckPair(Point2Cartesian(0, 0), Point2Cartesian(3 * scale, 4 * scale), 5 * scale);
		CheckPair(Pnt3Cart(0, 0, 0), Pnt3Cart(scale, 2 * scale, 2 * scale), 3 * scale);
	}
	const Real large = LargeScale, adjacent = std::nextafter(large, std::numeric_limits<Real>::infinity());
	CheckPair(Point2Cartesian(large, large), Point2Cartesian(adjacent, large), adjacent - large);
	CheckPair(Pnt3Cart(large, large, large), Pnt3Cart(large, adjacent, large), adjacent - large);
}

TEST_CASE("Equal-angle curvilinear point distances retain small radial differences",
	"[distance-contract][Point2Polar][Point3Spherical][Point3Cylindrical]") {
	const Real radius = Real(1e8);
	const Real increment = TOL(Real(1), std::nextafter(radius, std::numeric_limits<Real>::infinity()) - radius);
	REQUIRE(radius + increment > radius);
	for (Real angle : {Real(0), Real(0.7), Constants::PI}) {
		CheckPair(Point2Polar(radius, angle), Point2Polar(radius + increment, angle), increment);
		CheckPair(Pnt3Sph(radius, Real(0.7), angle), Pnt3Sph(radius + increment, Real(0.7), angle), increment);
		CheckPair(Pnt3Cyl(radius, angle, 2), Pnt3Cyl(radius + increment, angle, 2), increment);
	}
}

TEST_CASE("Distance contracts distinguish rounded coincidence from adjacent radial coordinates",
	"[distance-contract][Point2Polar][Point3Spherical][Point3Cylindrical]") {
	const Real radius = std::ldexp(Real(1), std::numeric_limits<Real>::digits);
	const Real rounded = radius + 1;
	const Real adjacent = std::nextafter(radius, std::numeric_limits<Real>::infinity());
	REQUIRE(rounded == radius);
	REQUIRE(adjacent - radius == 2);
	CheckPair(Point2Polar(radius, Real(0.7)), Point2Polar(rounded, Real(0.7)), 0);
	CheckPair(Pnt3Sph(radius, Real(0.7), Real(0.3)), Pnt3Sph(rounded, Real(0.7), Real(0.3)), 0);
	CheckPair(Pnt3Cyl(radius, Real(0.7), 2), Pnt3Cyl(rounded, Real(0.7), 2), 0);
	CheckPair(Point2Polar(radius, Real(0.7)), Point2Polar(adjacent, Real(0.7)), 2);
	CheckPair(Pnt3Sph(radius, Real(0.7), Real(0.3)), Pnt3Sph(adjacent, Real(0.7), Real(0.3)), 2);
	CheckPair(Pnt3Cyl(radius, Real(0.7), 2), Pnt3Cyl(adjacent, Real(0.7), 2), 2);
}

TEST_CASE("Cartesian distance stress scales cross the active squared-arithmetic range",
	"[distance-contract][Point2Cartesian][Point3Cartesian]") {
	const Real small = std::ldexp(Real(1), std::numeric_limits<Real>::min_exponent);
	const Real large = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 4);
	REQUIRE(small > 0);
	REQUIRE(small < std::sqrt(std::numeric_limits<Real>::denorm_min()));
	REQUIRE(large > std::sqrt(std::numeric_limits<Real>::max()));
	REQUIRE(std::isfinite(5 * large));
	for (Real scale : {small, large}) {
		CheckPair(Point2Cartesian(0, 0), Point2Cartesian(3 * scale, 4 * scale), 5 * scale);
		CheckPair(Pnt3Cart(0, 0, 0), Pnt3Cart(scale, 2 * scale, 2 * scale), 3 * scale);
	}
}

TEST_CASE("Curvilinear angular distances avoid cancelling cosines and premature underflow",
	"[distance-contract][Point2Polar][Point3Spherical][Point3Cylindrical]") {
	for (Real radius : {TinyScale, Real(1), LargeScale}) {
		const Real angle = Real(1e-9);
		const Real expected = radius * angle;
		CheckPair(Point2Polar(radius, 0), Point2Polar(radius, angle), expected);
		CheckPair(Pnt3Sph(radius, 0, 0), Pnt3Sph(radius, angle, 0), expected);
		CheckPair(Pnt3Sph(radius, Constants::PI / 2, 0), Pnt3Sph(radius, Constants::PI / 2, angle), expected);
		CheckPair(Pnt3Cyl(radius, 0, 0), Pnt3Cyl(radius, angle, 0), expected);
	}
	CheckPair(Point2Polar(LargeScale, 0), Point2Polar(LargeScale, TinyScale), 1);
	CheckPair(Point2Polar(LargeScale, TinyScale), Point2Polar(LargeScale, 2 * Constants::PI), 1);
	CheckPair(Point2Polar(-LargeScale, TinyScale), Point2Polar(LargeScale, Constants::PI), 1);
	CheckPair(Pnt3Sph(-LargeScale, TinyScale, 0), Pnt3Sph(LargeScale, Constants::PI, 0), 1);
	CheckPair(Pnt3Sph(FarRadius, TinyScale, 0),
		Pnt3Sph(FarRadius, TinyScale, TinyScale), (FarRadius * TinyScale) * TinyScale);
	CheckPair(Pnt3Cyl(4 * LargeScale, 0, 0), Pnt3Cyl(4 * LargeScale, TinyScale, 3), 5);
	CheckPair(Pnt3Cyl(TinyScale, 0, 0), Pnt3Cyl(TinyScale, TinyScale, 3), 3);
	const Real smallestAngle = std::numeric_limits<Real>::denorm_min();
	CheckPair(Point2Polar(FarRadius, 0), Point2Polar(FarRadius, smallestAngle), FarRadius * smallestAngle);
}

TEST_CASE("Curvilinear distances preserve signed radii and noncanonical coordinate conventions",
	"[distance-contract][Point2Polar][Point3Spherical][Point3Cylindrical]") {
	for (Real scale : {TinyScale, Real(1), LargeScale}) {
		const Point2Polar polar1(-2 * scale, Real(0.3)), polar2(3 * scale, Real(7.2));
		const Real x = -2 * std::cos(Real(0.3)) - 3 * std::cos(Real(7.2));
		const Real y = -2 * std::sin(Real(0.3)) - 3 * std::sin(Real(7.2));
		CheckPair(polar1, polar2, scale * std::hypot(x, y));
		const Real theta1 = Real(-0.4), theta2 = Real(3.7), phi1 = Real(0.2), phi2 = Real(-1.3);
		const Real dx = -2 * std::sin(theta1) * std::cos(phi1) - 3 * std::sin(theta2) * std::cos(phi2);
		const Real dy = -2 * std::sin(theta1) * std::sin(phi1) - 3 * std::sin(theta2) * std::sin(phi2);
		const Real dz = -2 * std::cos(theta1) - 3 * std::cos(theta2);
		CheckPair(Pnt3Sph(-2 * scale, theta1, phi1), Pnt3Sph(3 * scale, theta2, phi2),
			scale * std::hypot(dx, dy, dz));
	}
	CheckPair(Point2Polar(-3, 0), Point2Polar(3, Constants::PI), 0);
	CheckPair(Pnt3Sph(-3, 0, 1), Pnt3Sph(3, Constants::PI, 2), 0);
	CheckPair(Pnt3Sph(3, Constants::PI, 1), Pnt3Sph(3, Constants::PI, 2), 0);
}

TEST_CASE("Cartesian-to-curvilinear conversions share stable radius and polar-angle norms",
	"[distance-contract][Point2Polar][Point3Spherical][Point3Cylindrical]") {
	for (Real scale : {TinyScale, Real(1), LargeScale}) {
		const Point2Polar polar(Point2Cartesian(3 * scale, 4 * scale));
		const Pnt3Cart cart(3 * scale, 4 * scale, 12 * scale);
		const Pnt3Sph spherical(cart);
		const Pnt3Cyl cylindrical(cart);
		CheckDistance(polar.R(), 5 * scale);
		CheckDistance(spherical.R(), 13 * scale);
		CheckDistance(cylindrical.R(), 5 * scale);
		CheckDistance(spherical.Theta(), std::atan2(Real(5), Real(12)));
		CheckDistance(spherical.Phi(), std::atan2(Real(4), Real(3)));
		CheckDistance(SphericalGeometryCalculator::RadiusFromCartesian(cart), 13 * scale);
	}
}

TEST_CASE("Great-circle distances resolve tiny equatorial and meridian separations",
	"[distance-contract][GeometrySpherical]") {
	for (Real degrees : {Real(1e-9), TOL(1e-100, 1e-15), TinyScale}) {
		REQUIRE(degrees > 0);
		const Real expected = degrees * (Constants::PI / 180);
		CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 0, degrees), expected);
		CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, degrees, 0), expected);
		CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, degrees, 0, 0), expected);
		CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(LargeScale, 0, 0, 0, degrees),
			LargeScale * expected);
	}
	const Real smallest = std::numeric_limits<Real>::denorm_min();
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, TinyScale, 0, 360),
		TinyScale * (Constants::PI / 180));
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(LargeScale, 0, TinyScale, 0, 360),
		Constants::PI / 180);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(FarRadius, 0, 0, 0, smallest),
		(FarRadius * smallest) * (Constants::PI / 180));
}

TEST_CASE("Physical great-circle distance can survive an unrepresentable angular result",
	"[distance-contract][GeometrySpherical]") {
	const Real degrees = std::numeric_limits<Real>::denorm_min();
	REQUIRE(degrees > 0);
	REQUIRE(degrees * (Constants::PI / 180) == 0);
	CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 0, degrees), GeometryError);
	const Real expected = (FarRadius * degrees) * (Constants::PI / 180);
	REQUIRE(std::isnormal(expected));
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(FarRadius, 0, 0, 0, degrees), expected);
}

TEST_CASE("Great-circle distances preserve antipodes poles wrapping and physical units",
	"[distance-contract][GeometrySpherical]") {
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 0, 180), Constants::PI);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(90, 20, -90, 87), Constants::PI);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(90, 0, 90, 123), 0);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(-90, 0, -90, 123), 0);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(20, 30, 20, 390), 0);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 179, 0, -179), Constants::PI / 90);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 90, 0), Constants::PI / 2);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(2, 0, 0, 0, 90), Constants::PI);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 0, 0, 90), 0);
	const Real longitude = TOL(180 - Real(1e-7), std::nextafter(Real(180), Real(0)));
	REQUIRE(longitude < 180);
	const Real angle = SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 0, longitude);
	CHECK_THAT(angle, Catch::Matchers::WithinAbs(longitude * (Constants::PI / 180), DistanceRoundoff));
	CHECK(angle < Constants::PI);
	const Real hugeLongitude = std::numeric_limits<Real>::max();
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, hugeLongitude, 0, -hugeLongitude),
		std::abs(std::remainder(-2 * std::remainder(hugeLongitude, Real(360)), Real(360))) * (Constants::PI / 180));
}

TEST_CASE("Distance operations reject nonfinite inputs and overflowing required results",
	"[distance-contract][GeometrySpherical]") {
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(Point2Cartesian().Dist(Point2Cartesian(invalid, 0)), GeometryError);
		CHECK_THROWS_AS(Pnt3Cart().Dist(Pnt3Cart(0, 0, invalid)), GeometryError);
		CHECK_THROWS_AS(Point2Polar(0, invalid).Dist(Point2Polar()), GeometryError);
		CHECK_THROWS_AS(Pnt3Sph(0, invalid, 0).Dist(Pnt3Sph()), GeometryError);
		CHECK_THROWS_AS(Pnt3Cyl(0, 0, invalid).Dist(Pnt3Cyl()), GeometryError);
		CHECK_THROWS_AS(Point2Polar(Point2Cartesian(invalid, 0)), GeometryError);
		CHECK_THROWS_AS(Pnt3Sph(Pnt3Cart(0, invalid, 0)), GeometryError);
		CHECK_THROWS_AS(Pnt3Cyl(Pnt3Cart(0, 0, invalid)), GeometryError);
		for (int index = 0; index < 4; ++index) {
			std::array<Real, 4> arguments{0, 0, 1, 1};
			arguments[index] = invalid;
			CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(
				arguments[0], arguments[1], arguments[2], arguments[3]), GeometryError);
		}
		CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(invalid, 0, 0, 1, 1), GeometryError);
	}
	const Real largest = std::numeric_limits<Real>::max();
	CHECK_THROWS_AS(Point2Cartesian(largest, 0).Dist(Point2Cartesian(-largest, 0)), GeometryError);
	CHECK_THROWS_AS(Point2Polar(largest, 0).Dist(Point2Polar(largest, Constants::PI)), GeometryError);
	const Real tiny = std::numeric_limits<Real>::min();
	CHECK_THROWS_AS(Point2Polar(tiny, 0).Dist(Point2Polar(tiny, tiny)), GeometryError);
	CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(largest, 0, 0, 0, 180), GeometryError);
}

TEST_CASE("Great-circle geographic validation accepts periodic longitude and zero radius",
	"[distance-contract][GeometrySpherical]") {
	for (Real latitude : {Real(-91), Real(91)}) {
		CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(latitude, 0, 0, 0), GeometryError);
		CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, latitude, 0), GeometryError);
		CHECK_THROWS_AS(SphericalGeometryCalculator::SphericalFromLatLong(latitude, 0), GeometryError);
	}
	CHECK_THROWS_AS(SphericalGeometryCalculator::DistanceBetweenLatLong(-1, 0, 0, 1, 1), GeometryError);
	CheckDistance(SphericalGeometryCalculator::SphericalFromLatLong(0, 270).Phi(), 3 * Constants::PI / 2);
	CheckDistance(SphericalGeometryCalculator::DistanceBetweenLatLong(0, 0, 0, 1, 1), 0);
}
