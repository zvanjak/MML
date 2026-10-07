#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/DifferentialGeometry/Atlas.h>
#endif

using namespace MML;
using namespace MML::DifferentialGeometry;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::DifferentialGeometryTests
{
	TEST_CASE("Unit sphere stereographic atlas maps chart origins to poles", "[DifferentialGeometry][Atlas]")
	{
		auto atlas = MakeUnitSphereStereographicAtlas();
		REQUIRE(atlas.chartCount() == 2);
		REQUIRE(atlas.chartName(0) == "north");
		REQUIRE(atlas.chartName(1) == "south");

		Point<Real, 2, Parameter2> origin{ REAL(0.0), REAL(0.0) };
		auto northPole = atlas.chart(0).map_point(origin);
		auto southPole = atlas.chart(1).map_point(origin);

		REQUIRE_THAT(northPole[0], WithinAbs(REAL(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(northPole[1], WithinAbs(REAL(0.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(northPole[2], WithinAbs(REAL(1.0), TOL(1e-12, 1e-6)));
		REQUIRE_THAT(southPole[2], WithinAbs(-REAL(1.0), TOL(1e-12, 1e-6)));
		REQUIRE(atlas.selectChartForAmbientPoint(northPole) == 0);
		REQUIRE(atlas.selectChartForAmbientPoint(southPole) == 1);
	}

	TEST_CASE("Stereographic atlas transitions round trip on chart overlap", "[DifferentialGeometry][Atlas]")
	{
		auto atlas = MakeUnitSphereStereographicAtlas();
		Point<Real, 2, Parameter2> northCoords{ REAL(1.5), REAL(0.5) };

		auto southCoords = atlas.transition(0, 1, northCoords);
		auto northRoundTrip = atlas.transition(1, 0, southCoords);
		auto northPoint = atlas.chart(0).map_point(northCoords);
		auto southPoint = atlas.chart(1).map_point(southCoords);

		for (int i = 0; i < 2; i++)
			REQUIRE_THAT(northRoundTrip[i], WithinAbs(northCoords[i], TOL(1e-12, 1e-6)));
		for (int i = 0; i < 3; i++)
			REQUIRE_THAT(southPoint[i], WithinAbs(northPoint[i], TOL(1e-12, 1e-6)));
	}

	TEST_CASE("Stereographic atlas overlap metrics agree at self-transition points", "[DifferentialGeometry][Atlas]")
	{
		auto atlas = MakeUnitSphereStereographicAtlas();
		Point<Real, 2, Parameter2> northCoords{ REAL(1.0), REAL(0.0) };
		auto southCoords = atlas.transition(0, 1, northCoords);

		Metric<2, Parameter2> northMetric = atlas.chart(0).induced_metric(northCoords);
		Metric<2, Parameter2> southMetric = atlas.chart(1).induced_metric(southCoords);

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				REQUIRE_THAT(southMetric(i, j), WithinAbs(northMetric(i, j), TOL(1e-10, 1e-6)));
	}

	TEST_CASE("Stereographic transition rejects pole points outside overlap", "[DifferentialGeometry][Atlas]")
	{
		auto atlas = MakeUnitSphereStereographicAtlas();
		Point<Real, 2, Parameter2> pole{ REAL(0.0), REAL(0.0) };
		REQUIRE_THROWS_AS(atlas.transition(0, 1, pole), std::invalid_argument);
	}
}