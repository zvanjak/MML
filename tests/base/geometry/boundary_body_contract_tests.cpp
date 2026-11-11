#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Geometry/Geometry3DBodies/Geometry3DBodyBase.h>
#endif

#include <cmath>
#include <limits>

using namespace MML;

namespace {
	Real UnitLower(Real x) {
		return x >= 0 && x <= 1 ? Real(0) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real UnitUpper(Real x) {
		return x >= 0 && x <= 1 ? Real(1) : std::numeric_limits<Real>::quiet_NaN();
	}
	bool InUnitSquare(Real x, Real y) { return x >= 0 && x <= 1 && y >= 0 && y <= 1; }
	Real UnitBottom(Real x, Real y) {
		return InUnitSquare(x, y) ? Real(0) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real UnitRoof(Real x, Real y) {
		return InUnitSquare(x, y) ? 1 + x + y : std::numeric_limits<Real>::quiet_NaN();
	}
	Real UnitTop(Real x, Real y) {
		return InUnitSquare(x, y) ? Real(1) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real AliasedRoof(Real x, Real y) {
		if (!InUnitSquare(x, y)) return std::numeric_limits<Real>::quiet_NaN();
		const Real wave = std::sin(64 * Constants::PI * x);
		return 1 + wave * wave;
	}
	Real CurvedLower(Real x) { return UnitLower(x) + x * x; }
	Real CurvedUpper(Real x) { return std::fma(x, x, UnitUpper(x)); }
	bool InCurvedDomain(Real x, Real y) {
		return x >= 0 && x <= 1 && y >= CurvedLower(x) && y <= CurvedUpper(x);
	}
	Real CurvedBottom(Real x, Real y) {
		return InCurvedDomain(x, y) ? Real(0) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real CurvedRoof(Real x, Real y) {
		return InCurvedDomain(x, y) ? Real(1) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real ShearedLower(Real x) { return UnitLower(x) + x; }
	Real ShearedUpper(Real x) { return UnitUpper(x) + x; }
	bool InShearedDomain(Real x, Real y) { return x >= 0 && x <= 1 && y >= x && y <= x + 1; }
	Real ShearedBottom(Real x, Real y) {
		return InShearedDomain(x, y) ? Real(0) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real ShearedRoof(Real x, Real y) {
		return InShearedDomain(x, y) ? 1 + y : std::numeric_limits<Real>::quiet_NaN();
	}
	Real PinchedBottom(Real x, Real y) {
		return x >= 0 && x <= 1 && y >= 0 && y <= x ? Real(0) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real PinchedRoof(Real x, Real y) { return PinchedBottom(x, y) + 1; }
	Real ThrowY(Real) { throw GeometryError("Boundary callback must not run for guaranteed bounds"); }
	Real ThrowZ(Real, Real) { throw GeometryError("Boundary callback must not run for guaranteed bounds"); }
	Real InvalidZ(Real, Real) { return std::numeric_limits<Real>::quiet_NaN(); }
	Real ReversedLower(Real) { return 2; }
	Real TinyUpper(Real x) { return UnitLower(x) + std::ldexp(Real(1), -20); }
	Real NarrowWidth() { return std::ldexp(Real(1), -std::numeric_limits<Real>::digits / 3 - 3); }
	bool InNarrowX(Real x) { return x >= 4 && x <= 4 + NarrowWidth(); }
	Real NarrowLower(Real x) {
		return InNarrowX(x) ? Real(8) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real NarrowUpper(Real x) { return NarrowLower(x) + NarrowWidth(); }
	bool InNarrowDomain(Real x, Real y) {
		return InNarrowX(x) && y >= 8 && y <= 8 + NarrowWidth();
	}
	Real NarrowBottom(Real x, Real y) {
		return InNarrowDomain(x, y) ? Real(0) : std::numeric_limits<Real>::quiet_NaN();
	}
	Real NarrowRoof(Real x, Real y) {
		return InNarrowDomain(x, y) ? 1 + (x - 4) / NarrowWidth() : std::numeric_limits<Real>::quiet_NaN();
	}
	void CheckBoundaryValue(Real actual, Real expected, Real tolerance = TOL(2e-8, 2e-3)) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual, Catch::Matchers::WithinAbs(expected, tolerance));
	}
}

TEST_CASE("Function bounded solids do not advertise sampled extrema as guarantees", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, UnitLower, UnitUpper, UnitBottom, AliasedRoof, 1);
	CHECK_THROWS_AS(body.getBoundingBox(), GeometryError);
	CHECK_THROWS_AS(body.getBoundingSphere(), GeometryError);
	const auto sampled = body.EstimateBoundingBox();
	CHECK(sampled.Max().Z() < Real(1) + Real(1) / 16);
	CheckBoundaryValue(body.Volume(), 1);
	body.SetBoundingBox(Box3D(0, 0, 0, 1, 1, 2));
	const Pnt3Cart peak(Real(1) / 128, Real(1) / 2, 2);
	CHECK(body.getBoundingBox().Contains(peak));
	CHECK(body.getBoundingSphere().Contains(peak));
	CHECK(body.EstimateBoundingSphere().Radius() < body.getBoundingSphere().Radius());

	BoundaryGeometrySettings refined;
	refined.integrationIntervals = 256;
	CheckBoundaryValue(body.EstimateVolume(refined), Real(3) / 2);
	CheckBoundaryValue(body.EstimateBoundingBox(refined).Max().Z(), 2);
}

TEST_CASE("Function bounded area stays inside closed rectangular callback domains", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, UnitLower, UnitUpper, UnitBottom, UnitRoof, 1);
	CheckBoundaryValue(body.Volume(), 2);
	CheckBoundaryValue(body.SurfaceArea(), 9 + std::sqrt(Real(3)));
	const auto center = body.getCenter();
	CheckBoundaryValue(center.X(), Real(13) / 24);
	CheckBoundaryValue(center.Y(), Real(13) / 24);
	CheckBoundaryValue(center.Z(), Real(25) / 24);
	BoundaryGeometrySettings coarse;
	coarse.integrationIntervals = 2;
	coarse.derivativeStep = Precision::BoundaryDerivativeMaximumStep;
	CheckBoundaryValue(body.EstimateSurfaceArea(coarse), 9 + std::sqrt(Real(3)));
}

TEST_CASE("Function bounded area respects moving curved domains at every stencil node", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, CurvedLower, CurvedUpper, CurvedBottom, CurvedRoof, 1);
	CheckBoundaryValue(body.Volume(), 1);
	CheckBoundaryValue(body.SurfaceArea(), 4 + std::sqrt(Real(5)) + std::asinh(Real(2)) / 2);
}

TEST_CASE("Curved boundary callbacks and guards share represented endpoints", "[boundary-body][boundary-contract][regression]")
{
	for (int panel = 1; panel < Precision::BoundaryIntegrationIntervals; ++panel) {
		const Real center = Real(panel) / Precision::BoundaryIntegrationIntervals;
		for (Real x : {center - Precision::CbrtEps, center, center + Precision::CbrtEps}) {
			const Real lower = CurvedLower(x), upper = CurvedUpper(x);
			CAPTURE(x, lower, upper);
			CHECK(InCurvedDomain(x, lower));
			CHECK(InCurvedDomain(x, upper));
			CHECK_FALSE(InCurvedDomain(x, std::nextafter(lower, Constants::NegInf)));
			CHECK_FALSE(InCurvedDomain(x, std::nextafter(upper, Constants::PosInf)));
		}
	}
}

TEST_CASE("Function bounded area includes the parameterization chain rule", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, ShearedLower, ShearedUpper, ShearedBottom, ShearedRoof, 1);
	CheckBoundaryValue(body.Volume(), 2);
	CheckBoundaryValue(body.SurfaceArea(), 5 + 5 * std::sqrt(Real(2)));
	const auto center = body.getCenter();
	CheckBoundaryValue(center.X(), Real(13) / 24);
	CheckBoundaryValue(center.Y(), Real(13) / 12);
	CheckBoundaryValue(center.Z(), Real(25) / 24);
}

TEST_CASE("Function bounded area handles domains pinched to zero width at an endpoint", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, UnitLower, ShearedLower, PinchedBottom, PinchedRoof, 1);
	CheckBoundaryValue(body.Volume(), Real(1) / 2);
	CheckBoundaryValue(body.SurfaceArea(), 3 + std::sqrt(Real(2)));
}

TEST_CASE("Guaranteed function bounds are independent of numerical measurements", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, ThrowY, ThrowY, ThrowZ, ThrowZ, 1);
	Box3D supplied(0, -2, -3, 1, 4, 5);
	body.SetBoundingBox(supplied);
	supplied.Max().Z() = 0;
	const IBody& polymorphic = body;
	CHECK(polymorphic.getBoundingBox().Max().Z() == 5);
	const auto sphere = polymorphic.getBoundingSphere();
	for (Real x : {Real(0), Real(1)})
		for (Real y : {Real(-2), Real(4)})
			for (Real z : {Real(-3), Real(5)})
				CHECK(sphere.Center().Dist(Pnt3Cart(x, y, z)) <= sphere.Radius());
	CHECK_THROWS_AS(body.Volume(), GeometryError);
	auto copy = body;
	CHECK(copy.getBoundingBox().Max().Z() == 5);
	body.ClearBoundingBox();
	CHECK_THROWS_AS(body.getBoundingSphere(), GeometryError);
	CHECK(copy.getBoundingSphere().Radius() == sphere.Radius());
}

TEST_CASE("Guaranteed function bounds reject stale public boundary definitions", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity body(0, 1, UnitLower, UnitUpper, UnitBottom, UnitRoof, 1);
	body.SetBoundingBox(Box3D(0, 0, 0, 1, 1, 3));
	body._x2 = Real(1) / 2;
	CHECK_THROWS_AS(body.getBoundingBox(), GeometryError);
	body.SetBoundingBox(Box3D(0, 0, 0, 1, 1, 3));
	CHECK(body.getBoundingSphere().Radius() > 0);
	body._z2 = UnitTop;
	CHECK_THROWS_AS(body.getBoundingSphere(), GeometryError);
	body.SetBoundingBox(Box3D(0, 0, 0, 1, 1, 1));
	CHECK(body.getBoundingBox().Max().Z() == 1);
}

TEST_CASE("Guaranteed function spheres handle point boxes and large translated boxes", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity point(0, 0, ThrowY, ThrowY, ThrowZ, ThrowZ, 1);
	point.SetBoundingBox(Box3D(0, 2, 3, 0, 2, 3));
	const auto pointSphere = point.getBoundingSphere();
	CHECK(pointSphere.Radius() == 0);
	CHECK(pointSphere.Center().Y() == 2);
	CHECK(pointSphere.Center().Z() == 3);

	SolidBodyWithBoundaryConstDensity body(0, 1, ThrowY, ThrowY, ThrowZ, ThrowZ, 1);
	const Real large = std::numeric_limits<Real>::max() / 2;
	const Real upper = large + large / 2;
	body.SetBoundingBox(Box3D(0, large, 0, 1, upper, 0));
	const auto sphere = body.getBoundingSphere();
	CHECK(std::isfinite(sphere.Center().Y()));
	CHECK(std::isfinite(sphere.Radius()));
	for (Real x : {Real(0), Real(1)})
		for (Real y : {large, upper})
			CHECK(sphere.Center().Dist(Pnt3Cart(x, y, 0)) <= sphere.Radius());
}

TEST_CASE("Function geometry exposes usable resolution controls on narrow domains", "[boundary-body][boundary-contract]")
{
	const Real width = NarrowWidth();
	SolidBodyWithBoundaryConstDensity body(4, 4 + width, NarrowLower, NarrowUpper, NarrowBottom, NarrowRoof, 1);
	CHECK(std::isfinite(body.SurfaceArea()));
	BoundaryGeometrySettings settings;
	settings.derivativeStep = Real(1) / 8;
	CheckBoundaryValue(body.EstimateSurfaceArea(settings) / width, width + std::hypot(width, Real(1)) + 6);
	CheckBoundaryValue(body.EstimateVolume(settings) / (width * width), Real(3) / 2);
}

TEST_CASE("Function geometry validates definition samples settings and supplied boxes", "[boundary-body][boundary-contract]")
{
	CHECK_THROWS_AS((SolidBodyWithBoundaryConstDensity(1, 0, UnitLower, UnitUpper, UnitBottom, UnitTop, 1)), GeometryError);
	CHECK_THROWS_AS((SolidBodyWithBoundaryConstDensity(0, 1, nullptr, UnitUpper, UnitBottom, UnitTop, 1)), GeometryError);
	SolidBodyWithBoundaryConstDensity body(0, 1, UnitLower, UnitUpper, UnitBottom, UnitTop, 1);
	for (int intervals : {-2, 0, 1, 3}) {
		BoundaryGeometrySettings settings;
		settings.integrationIntervals = intervals;
		CHECK_THROWS_AS(body.EstimateVolume(settings), GeometryError);
		CHECK_THROWS_AS(body.EstimateSurfaceArea(settings), GeometryError);
		CHECK_THROWS_AS(body.EstimateCenter(settings), GeometryError);
		CHECK_THROWS_AS(body.EstimateBoundingBox(settings), GeometryError);
	}
	for (Real step : {Real(0), Real(-1), Real(1), std::numeric_limits<Real>::infinity(),
		std::numeric_limits<Real>::quiet_NaN()}) {
		BoundaryGeometrySettings settings;
		settings.derivativeStep = step;
		CHECK_THROWS_AS(body.EstimateSurfaceArea(settings), GeometryError);
	}
	BoundaryGeometrySettings unresolved;
	unresolved.derivativeStep = std::numeric_limits<Real>::denorm_min();
	CHECK_THROWS_AS(body.EstimateSurfaceArea(unresolved), GeometryError);
	CHECK_THROWS_AS(body.SetBoundingBox(Box3D(0, 0, 0, Real(1) / 2, 1, 1)), GeometryError);
	CHECK_THROWS_AS(body.SetBoundingBox(Box3D(0, 0, 0, 1, 1, std::numeric_limits<Real>::infinity())), GeometryError);
	Box3D unordered;
	unordered.Min().Y() = 2;
	CHECK_THROWS_AS(body.SetBoundingBox(unordered), GeometryError);
	body._y1 = ReversedLower;
	CHECK_THROWS_AS(body.Volume(), GeometryError);
	body._y1 = UnitLower;
	body._z2 = UnitBottom;
	body._z1 = UnitTop;
	CHECK_THROWS_AS(body.EstimateBoundingBox(), GeometryError);
	body._z1 = UnitBottom;
	body._z2 = InvalidZ;
	CHECK_THROWS_AS(body.Volume(), GeometryError);
	CHECK_THROWS_AS(body.SurfaceArea(), GeometryError);
	CHECK_THROWS_AS(body.getCenter(), GeometryError);
	CHECK_THROWS_AS(body.EstimateBoundingBox(), GeometryError);
	body._z2 = nullptr;
	CHECK_THROWS_AS(body.IsInside(Pnt3Cart(Real(1) / 2, Real(1) / 2, Real(1) / 2)), GeometryError);
}

TEST_CASE("Function geometry preserves zero and small volume distinctions", "[boundary-body][boundary-contract]")
{
	SolidBodyWithBoundaryConstDensity flat(0, 1, UnitLower, UnitUpper, UnitBottom, UnitBottom, 1);
	CHECK(flat.Volume() == 0);
	CHECK_THROWS_AS(flat.getCenter(), GeometryError);
	flat.SetBoundingBox(Box3D(0, 0, 0, 1, 1, 0));
	CHECK(std::isfinite(flat.getBoundingSphere().Radius()));
	const Real width = std::ldexp(Real(1), -20);
	SolidBodyWithBoundaryConstDensity tiny(0, width, UnitLower, TinyUpper, UnitBottom, UnitTop, 1);
	CheckBoundaryValue(tiny.Volume() / (width * width), 1);
	const auto center = tiny.getCenter();
	CheckBoundaryValue(center.X() / width, Real(1) / 2);
	CheckBoundaryValue(center.Y() / width, Real(1) / 2);
	CheckBoundaryValue(center.Z(), Real(1) / 2);
}
