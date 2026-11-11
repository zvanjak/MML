#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Integration/SurfaceIntegration.h>
#endif

#include <array>
#include <limits>

using namespace MML;

namespace {
	const Real TinyDirection = std::numeric_limits<Real>::min() * Real(16);
	const Real SmallSide = std::sqrt(std::numeric_limits<Real>::min()) * Real(16);
	const Real LargeSide = 1 / SmallSide;

	void CheckRelative(Real actual, Real expected) {
		REQUIRE(std::isfinite(actual));
		CHECK_THAT(actual / expected, Catch::Matchers::WithinAbs(Real(1),
			Real(64) * std::numeric_limits<Real>::epsilon()));
	}

	struct QuadraticFlux : IVectorFunction<3> {
		Real side;
		bool translated = false;
		mutable int calls = 0;
		explicit QuadraticFlux(Real length) : side(length) {}
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override {
			++calls;
			const Real coordinate = translated ? p[1] : p[0] / side;
			const Real density = (1 + coordinate * coordinate) / side / side;
			return translated ? VectorN<Real, 3>{density, 0, 0} : VectorN<Real, 3>{0, 0, density};
		}
	};
}

TEST_CASE("Plane directions ignore physical magnitude", "[plane-rectangle-scale][Plane3D]") {
	for (const Real magnitude : {std::numeric_limits<Real>::denorm_min(), TinyDirection,
		Real(1e-8), Real(1), 1 / TinyDirection, std::numeric_limits<Real>::max()}) {
		CAPTURE(magnitude);
		const Plane3D plane(Pnt3Cart(2, 3, 4), Vec3Cart(magnitude, -magnitude, magnitude));
		const Real component = 1 / std::sqrt(Real(3));
		CheckRelative(plane.Normal().X(), component);
		CheckRelative(plane.Normal().Y(), -component);
		CheckRelative(plane.Normal().Z(), component);
		CHECK_THAT(plane.DistToPoint(Pnt3Cart(2, 3, 4)), Catch::Matchers::WithinAbs(Real(0),
			Real(64) * std::numeric_limits<Real>::epsilon()));
	}
}

TEST_CASE("Three-point planes retain tiny large and cancelling cross products", "[plane-rectangle-scale][Plane3D]") {
	for (const Real side : {TinyDirection, Real(1e-8), Real(1), 1 / TinyDirection}) {
		CAPTURE(side);
		const Plane3D plane(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0), Pnt3Cart(side, side, 0));
		CHECK(plane.Normal().Z() == Real(1));
		const Plane3D reversed(Pnt3Cart(0, 0, 0), Pnt3Cart(side, side, 0), Pnt3Cart(side, 0, 0));
		CHECK(reversed.Normal().Z() == Real(-1));
	}
	const Real huge = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 1);
	const Real tiny = 1 / huge;
	const Plane3D spanning(Pnt3Cart(-huge, 0, 0), Pnt3Cart(0, 0, 0), Pnt3Cart(huge, tiny, 0));
	CHECK(spanning.Normal().Z() == Real(1));
	const Real adjacent = std::nextafter(Real(1), Real(2));
	const Plane3D cancellation(Pnt3Cart(0, 0, 0), Pnt3Cart(adjacent, 1, 0), Pnt3Cart(1, 2 - adjacent, 0));
	CHECK(cancellation.Normal().Z() == Real(-1)); // Determinant is -(adjacent - 1)^2.
}

TEST_CASE("Plane offsets do not overflow before cancellation", "[plane-rectangle-scale][Plane3D]") {
	const Real large = std::numeric_limits<Real>::max();
	const Plane3D plane(Pnt3Cart(large, large, -large), Vec3Cart(1, 1, 1));
	CheckRelative(plane.D(), -large / std::sqrt(Real(3)));
	const Plane3D residual(Pnt3Cart(1 / TinyDirection, 1 / TinyDirection, TinyDirection), Vec3Cart(1, -1, 1));
	CheckRelative(residual.D(), -TinyDirection / std::sqrt(Real(3)));
}

TEST_CASE("Original microscopic square constructs before integration", "[plane-rectangle-scale][RectSurface3D]") {
	const Real side = Real(1e-8);
	const RectSurface3D square(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0),
		Pnt3Cart(side, side, 0), Pnt3Cart(0, side, 0));
	CheckRelative(square.getArea(), side * side);
	CHECK(square.getNormal().Z() == Real(1));
}

TEST_CASE("Quadrilateral areas and frames survive physical scale", "[plane-rectangle-scale][RectSurface3D]") {
	for (const Real side : {SmallSide, Real(1e-8), Real(1), LargeSide,
		std::sqrt(std::numeric_limits<Real>::max()) * Real(0.9)}) {
		CAPTURE(side);
		const RectSurface3D square(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0),
			Pnt3Cart(side, side, 0), Pnt3Cart(0, side, 0));
		CheckRelative(square.getArea(), side * side);
		CHECK(square.getNormal().Z() == Real(1));
		CHECK(square._localX.X() == Real(1));
		CHECK(square._localY.Y() == Real(1));
		CheckRelative(square.getCenter().X(), side / 2);
		CheckRelative(square(0, 0)[1], side / 2);
		const RectSurface3D reversed(square._pnt1, square._pnt4, square._pnt3, square._pnt2);
		CheckRelative(reversed.getArea(), square.getArea());
		CHECK(reversed.getNormal().Z() == Real(-1));
	}
	const Real huge = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 1);
	const Real tiny = 1 / huge;
	const RectSurface3D slender(Pnt3Cart(-huge, 0, 0), Pnt3Cart(0, 0, 0),
		Pnt3Cart(huge, tiny, 0), Pnt3Cart(0, tiny, 0));
	CheckRelative(slender.getArea(), Real(1));
	CHECK(slender.getNormal().Z() == Real(1));
	const RectSurface3D anisotropic(Pnt3Cart(0, 0, 0), Pnt3Cart(1 / TinyDirection, 0, 0),
		Pnt3Cart(1 / TinyDirection, TinyDirection, 0), Pnt3Cart(0, TinyDirection, 0));
	CheckRelative(anisotropic.getArea(), Real(1));
}

TEST_CASE("Quadrilateral coplanarity is relative but plane membership remains absolute", "[plane-rectangle-scale][RectSurface3D]") {
	for (const Real side : {SmallSide, Real(1), LargeSide}) {
		CAPTURE(side);
		const Real tolerance = Precision::Plane3DIsPointOnPlaneTolerance;
		CHECK_NOTHROW(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0),
			Pnt3Cart(side, side, 0), Pnt3Cart(0, side, side * tolerance / 4)));
		CHECK_THROWS_AS(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0),
			Pnt3Cart(side, side, 0), Pnt3Cart(0, side, side / 4)), GeometryError);
	}
	const Plane3D plane = Plane3D::XYPlane();
	CHECK(plane.IsPointOnPlane(Pnt3Cart(0, 0, Precision::Plane3DIsPointOnPlaneTolerance / 2)));
	CHECK_FALSE(plane.IsPointOnPlane(Pnt3Cart(0, 0, Precision::Plane3DIsPointOnPlaneTolerance * 2)));
}

TEST_CASE("Geometry rejects genuine degeneracy nonfinite inputs and unrepresentable results", "[plane-rectangle-scale]") {
	const Pnt3Cart origin(0, 0, 0);
	CHECK_THROWS_AS(Plane3D(origin, Vec3Cart(0, 0, 0)), GeometryError);
	CHECK_THROWS_AS(Plane3D(origin, Pnt3Cart(1, 1, 1), Pnt3Cart(2, 2, 2)), GeometryError);
	CHECK_THROWS_AS(RectSurface3D(origin, Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 0), origin), GeometryError);
	for (const Real bad : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(Plane3D(origin, Vec3Cart(bad, 0, 1)), GeometryError);
		CHECK_THROWS_AS(Plane3D(Pnt3Cart(bad, 0, 0), Vec3Cart(0, 0, 1)), GeometryError);
		CHECK_THROWS_AS(RectSurface3D(origin, Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 0), Pnt3Cart(0, 1, bad)), GeometryError);
	}
	const Real max = std::numeric_limits<Real>::max();
	CHECK_THROWS_AS(Plane3D(Pnt3Cart(max, max, max), Vec3Cart(1, 1, 1)), GeometryError);
	CHECK_THROWS_AS(Plane3D(Pnt3Cart(std::numeric_limits<Real>::denorm_min(), 0, 0), Vec3Cart(1, 3, 0)), GeometryError);
	for (const Real side : {std::numeric_limits<Real>::min(), std::numeric_limits<Real>::max() / 4}) {
		const RectSurface3D square(origin, Pnt3Cart(side, 0, 0), Pnt3Cart(side, side, 0), Pnt3Cart(0, side, 0));
		CHECK_THROWS_AS(square.getArea(), GeometryError);
	}
}

TEST_CASE("Adaptive rectangular refinement retains scaled flux and large-coordinate centers", "[plane-rectangle-scale][geometric-integration]") {
	for (const Real side : {SmallSide, Real(1e-8), Real(1), LargeSide}) {
		CAPTURE(side);
		const RectSurface3D square(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0),
			Pnt3Cart(side, side, 0), Pnt3Cart(0, side, 0));
		QuadraticFlux field(side);
		const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, square, Real(1e-4), 6);
		CHECK(result.converged);
		CHECK(result.status == AlgorithmStatus::Success);
		CHECK(result.function_evaluations == field.calls);
		CHECK(field.calls > 5);
		CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(4) / 3, Real(1e-4)));
	}
	const Real offset = std::numeric_limits<Real>::max() * Real(0.75);
	const RectSurface3D translated(Pnt3Cart(offset, 0, 0), Pnt3Cart(offset, 1, 0),
		Pnt3Cart(offset, 1, 1), Pnt3Cart(offset, 0, 1));
	CHECK(translated.getCenter().X() == offset);
	CHECK(translated(0, 0)[0] == offset);
	QuadraticFlux field(1);
	field.translated = true;
	const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, translated, Real(1e-4), 6);
	CHECK(result.converged);
	CHECK(result.function_evaluations == field.calls);
	CHECK(field.calls > 5);
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(4) / 3, Real(1e-4)));
}
