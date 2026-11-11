#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Geodesic.h>
#include <mml/core/DifferentialGeometry/InducedMetric.h>
#endif

using namespace MML;
using namespace MML::DifferentialGeometry;
using namespace MML::Surfaces;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::DifferentialGeometryTests
{
	TEST_CASE("InducedMetric2D reproduces first fundamental form coefficients", "[DifferentialGeometry][InducedMetric]")
	{
		Sphere sphere(static_cast<MML::Real>(2.0));
		InducedMetric2D metric(sphere);

		Real u = Constants::PI / static_cast<MML::Real>(3.0);
		Real w = Constants::PI / static_cast<MML::Real>(5.0);
		VectorN<Real, 2> pos{u, w};

		Real E, F, G;
		sphere.GetFirstNormalFormCoefficients(u, w, E, F, G);

		REQUIRE_THAT(metric.component(0, 0, pos), WithinAbs(E, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(metric.component(0, 1, pos), WithinAbs(F, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(metric.component(1, 0, pos), WithinAbs(F, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(metric.component(1, 1, pos), WithinAbs(G, TOL(1e-10, 1e-5)));

		auto matrix = FirstFundamentalFormMatrix(sphere, u, w);
		REQUIRE_THAT(matrix(0, 0), WithinAbs(E, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(matrix(0, 1), WithinAbs(F, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(matrix(1, 0), WithinAbs(F, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(matrix(1, 1), WithinAbs(G, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("InducedMetric2D exposes geodesic-compatible plane metric", "[DifferentialGeometry][InducedMetric][Geodesic]")
	{
		PlaneSurface plane(Vec3Cart{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, Vec3Cart{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
		InducedMetric2D metric(plane);

		VectorN<Real, 2> pos{static_cast<MML::Real>(0.25), static_cast<MML::Real>(-0.5)};
		REQUIRE_THAT(metric.component(0, 0, pos), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(metric.component(0, 1, pos), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(metric.component(1, 1, pos), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-8, 1e-4)));

		auto solution = IntegrateSurfaceGeodesicFixedStep(
			plane,
			VectorN<Real, 2>{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)},
			VectorN<Real, 2>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5)},
			static_cast<MML::Real>(0.0),
			static_cast<MML::Real>(1.0),
			40);

		Vector<Real> finalState = solution.getXValuesAtEnd();
		auto finalPosition = GeodesicPositionFromState<2>(finalState);
		auto finalVelocity = GeodesicVelocityFromState<2>(finalState);

		REQUIRE_THAT(finalPosition[0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(finalPosition[1], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(finalVelocity[0], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(finalVelocity[1], WithinAbs(static_cast<MML::Real>(0.5), TOL(1e-5, 1e-3)));
	}

	TEST_CASE("Intrinsic Gaussian curvature agrees for plane and unit sphere", "[DifferentialGeometry][InducedMetric][Curvature]")
	{
		PlaneSurface plane(Vec3Cart{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)}, Vec3Cart{static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)});
		REQUIRE_THAT(GaussianCurvatureIntrinsic(plane, static_cast<MML::Real>(0.25), static_cast<MML::Real>(0.5)), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-3)));

		Sphere sphere(static_cast<MML::Real>(1.0));
		Real u = Constants::PI / static_cast<MML::Real>(3.0);
		Real w = Constants::PI / static_cast<MML::Real>(4.0);
		Real intrinsic = GaussianCurvatureIntrinsic(sphere, u, w);
		Real extrinsic = sphere.GaussianCurvature(u, w);

		REQUIRE_THAT(intrinsic, WithinAbs(static_cast<MML::Real>(1.0), TOL(5e-3, 5e-2)));
		REQUIRE_THAT(intrinsic, WithinAbs(extrinsic, TOL(5e-3, 5e-2)));
	}
}
