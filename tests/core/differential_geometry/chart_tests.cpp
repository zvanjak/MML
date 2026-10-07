#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/DifferentialGeometry/Chart.h>
#include <mml/core/Surfaces.h>
#endif

using namespace MML;
using namespace MML::DifferentialGeometry;
using namespace MML::Surfaces;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::DifferentialGeometryTests
{
	TEST_CASE("FunctionChart maps points and tangent vectors for a plane patch", "[DifferentialGeometry][Chart]")
	{
		PlaneSurface plane(Vec3Cart{REAL(0.0), REAL(0.0), REAL(0.0)}, Vec3Cart{REAL(0.0), REAL(0.0), REAL(1.0)});
		FunctionChart<2, 3, Parameter2, Cartesian3> chart(plane);
		Point<Real, 2, Parameter2> parameterPoint{ REAL(0.25), REAL(-0.5) };

		auto ambientPoint = map_point(chart, parameterPoint);
		REQUIRE_THAT(ambientPoint[2], WithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));

		TangentVector<2, Parameter2> parameterVector{ REAL(2.0), REAL(-3.0) };
		auto ambientVector = push_forward(chart, parameterVector, parameterPoint);

		Real ambientNorm = std::sqrt(ambientVector[0] * ambientVector[0] + ambientVector[1] * ambientVector[1] + ambientVector[2] * ambientVector[2]);
		Real parameterNorm = std::sqrt(parameterVector[0] * parameterVector[0] + parameterVector[1] * parameterVector[1]);
		REQUIRE_THAT(ambientNorm, WithinAbs(parameterNorm, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(ambientVector[2], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("FunctionChart pulls back covectors and forms", "[DifferentialGeometry][Chart]")
	{
		PlaneSurface plane(Vec3Cart{REAL(0.0), REAL(0.0), REAL(0.0)}, Vec3Cart{REAL(0.0), REAL(0.0), REAL(1.0)});
		FunctionChart<2, 3, Parameter2, Cartesian3> chart(plane);
		Point<Real, 2, Parameter2> parameterPoint{ REAL(0.0), REAL(0.0) };

		Covector<3, Cartesian3> dz{ REAL(0.0), REAL(0.0), REAL(1.0) };
		auto dzPull = pull_back(chart, dz, parameterPoint);
		REQUIRE_THAT(dzPull[0], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(dzPull[1], WithinAbs(REAL(0.0), TOL(1e-8, 1e-4)));

		Form2<3, Cartesian3> areaXY;
		areaXY.SetAlternatingComponent(REAL(1.0), 0, 1);
		auto pulledArea = pull_back(chart, areaXY, parameterPoint);
		TangentVector<2, Parameter2> eu{ REAL(1.0), REAL(0.0) };
		TangentVector<2, Parameter2> ev{ REAL(0.0), REAL(1.0) };

		REQUIRE(pulledArea.IsAlternating());
		REQUIRE_THAT(std::abs(pulledArea(eu, ev)), WithinAbs(REAL(1.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("FunctionChart induced metric matches unit sphere chart", "[DifferentialGeometry][Chart]")
	{
		Sphere sphere(REAL(1.0));
		FunctionChart<2, 3, Parameter2, Cartesian3> chart(sphere);
		Point<Real, 2, Parameter2> parameterPoint{ Constants::PI / REAL(3.0), Constants::PI / REAL(4.0) };

		Metric<2, Parameter2> metric = chart.induced_metric(parameterPoint);

		REQUIRE_THAT(metric(0, 0), WithinAbs(REAL(1.0), TOL(1e-6, 1e-3)));
		REQUIRE_THAT(metric(0, 1), WithinAbs(REAL(0.0), TOL(1e-6, 1e-3)));
		REQUIRE_THAT(metric(1, 0), WithinAbs(REAL(0.0), TOL(1e-6, 1e-3)));
		REQUIRE_THAT(metric(1, 1), WithinAbs(std::sin(parameterPoint[0]) * std::sin(parameterPoint[0]), TOL(1e-6, 1e-3)));
	}
}