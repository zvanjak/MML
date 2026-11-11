#define MML_USE_SINGLE_HEADER
#include <MML.h>
#include <test_beds/parametric_curves_test_bed.h>
#include <test_beds/parametric_surfaces_test_bed.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Shared showcase fixtures explicitly opt in on top of the standalone distribution",
          "[single-header][showcase][testbeds]")
{
	using namespace MML;
	const auto& planar = TestBeds::ParametricCurvesTestBed::getShowcaseCurves2D();
	const auto& spatial = TestBeds::ParametricCurvesTestBed::getShowcaseCurves3D();
	const auto& surfaces = TestBeds::ParametricSurfaceTestBed::getShowcaseSurfaces();
	REQUIRE(planar.size() == 5);
	REQUIRE(spatial.size() == 2);
	REQUIRE(surfaces.size() == 6);
	REQUIRE(planar.front().curve(Real(0))[0] == Real(1));
	REQUIRE(spatial.front().curve(Real(2))[2] == Real(8));
	REQUIRE(surfaces.front().surface(Real(1), Real(2))[2] == Real(-11));
}
