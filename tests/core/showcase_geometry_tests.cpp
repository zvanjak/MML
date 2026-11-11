#include <catch2/catch_all.hpp>
#include <test_beds/parametric_curves_test_bed.h>
#include <test_beds/parametric_surfaces_test_bed.h>

#include <cmath>
#include <set>
#include <string>
#include <type_traits>

using namespace MML;
namespace CurveCases = MML::TestBeds::ShowcaseCurves;
namespace SurfaceCases = MML::TestBeds::ShowcaseSurfaces;

namespace {
	template<int Dimension>
	void RequirePoint(const VectorN<Real, Dimension>& actual, const VectorN<Real, Dimension>& expected) {
		for (int i = 0; i < Dimension; ++i)
			REQUIRE(IsNearlyEqual(actual[i], expected[i], Precision::EqualityTolerance, Precision::EqualityTolerance));
	}
}

TEST_CASE("Showcase curve definitions preserve formulas, parameters and domains", "[showcase][curves][testbeds]")
{
	STATIC_REQUIRE(std::is_base_of_v<Curves::ICurveCartesian2D, CurveCases::ButterflyCurve>);
	STATIC_REQUIRE(std::is_base_of_v<Curves::ICurveCartesian3D, CurveCases::TwistedCubicCurve>);
	CurveCases::LemniscateCurve lemniscate;
	RequirePoint(lemniscate(Real(0)), VectorN<Real, 2>{1, 0});
	REQUIRE(lemniscate.getMinT() == Constants::NegInf);
	REQUIRE(lemniscate.getMaxT() == Constants::PosInf);
	CurveCases::DeltoidCurve deltoid(2);
	RequirePoint(deltoid(Real(0)), VectorN<Real, 2>{8, 0});
	CurveCases::AstroidCurve astroid(3);
	RequirePoint(astroid(Real(0)), VectorN<Real, 2>{3, 0});
	RequirePoint(astroid(Constants::PI / 2), VectorN<Real, 2>{0, 3});
	REQUIRE_THROWS_AS(CurveCases::AstroidCurve(Real(0)), ArgumentError);
	REQUIRE_THROWS_AS(CurveCases::AstroidCurve(Real(-1)), ArgumentError);
	CurveCases::EpitrochoidCurve epitrochoid(1, Real(1) / 2, 3);
	RequirePoint(epitrochoid(Real(0)), VectorN<Real, 2>{Real(1) / 2, 0});
	CurveCases::ButterflyCurve butterfly(2);
	RequirePoint(butterfly(Real(0)), VectorN<Real, 2>{0, 2 * (std::exp(Real(1)) - 2)});
	REQUIRE(butterfly.getMinT() == 0);
	REQUIRE(butterfly.getMaxT() == 12 * Constants::PI);
	CurveCases::TwistedCubicCurve cubic;
	RequirePoint(cubic(Real(2)), VectorN<Real, 3>{2, 4, 8});
	for (Real t : {Real(-2), Real(0), Real(1) / 3, Real(2)})
		RequirePoint(cubic(t), TestBeds::TwistedCubic_pos(t));
	CurveCases::ToroidalSpiralCurve toroidal(3, Real(2));
	RequirePoint(toroidal(Real(0)), VectorN<Real, 3>{8, 0, 2});
	RequirePoint(toroidal(Constants::PI / 2), VectorN<Real, 3>{0, 6, 0});
	REQUIRE(toroidal.getMinT() == Constants::NegInf);
	REQUIRE(toroidal.getMaxT() == Constants::PosInf);
	RequirePoint(CurveCases::ToroidalSpiralCurve(Real(3))(0), VectorN<Real, 3>{12, 0, 3});
	RequirePoint(CurveCases::ToroidalSpiralCurve(3)(0), VectorN<Real, 3>{4, 0, 1});
}

TEST_CASE("Showcase surface definitions preserve defaults, scale and branch behavior", "[showcase][surfaces][testbeds]")
{
	STATIC_REQUIRE(std::is_base_of_v<Surfaces::ISurfaceCartesian, SurfaceCases::DiniSurface>);
	SurfaceCases::MonkeySaddle saddle;
	RequirePoint(saddle(1, 2), VectorN<Real, 3>{1, 2, -11});
	REQUIRE(saddle.getMinU() == -10);
	REQUIRE(saddle.getMaxW() == 10);
	SurfaceCases::MobiusStrip mobius(2);
	RequirePoint(mobius(0, Real(1) / 2), VectorN<Real, 3>{3, 0, 0});
	RequirePoint(mobius(0, Real(1) / 3), mobius(2 * Constants::PI, -Real(1) / 3));
	SurfaceCases::Helicoid helicoid(2);
	RequirePoint(helicoid(3, Constants::PI / 2), VectorN<Real, 3>{0, 3, Constants::PI});
	REQUIRE(helicoid.getMinU() == -30);
	REQUIRE(helicoid.getMaxW() == 4 * Constants::PI);
	SurfaceCases::EnneperSurface enneper(2);
	RequirePoint(enneper(1, 2), VectorN<Real, 3>{Real(28) / 3, Real(8) / 3, -6});
	for (Real u : {Real(-1), Real(0), Real(1)})
		for (Real v : {Real(-1), Real(0), Real(1)})
			RequirePoint(enneper(u, v), 2 * TestBeds::EnneperSurface(u, v));
	SurfaceCases::KleinBottle klein(2);
	RequirePoint(klein(0, 0), VectorN<Real, 3>{14, 0, 0});
	// The branch at pi is intentionally the second half of the original immersion.
	for (Real u : {Constants::PI / 2, Constants::PI, 3 * Constants::PI / 2}) {
		const Real v = Constants::PI / 3;
		const Real r = 2 - std::cos(v);
		const Real x = u < Constants::PI
			? 6 * std::cos(u) * (1 + std::sin(u)) + r * std::cos(u / 2) * std::cos(v)
			: 6 * std::cos(u) * (1 + std::sin(u)) - r * std::cos(u / 2) * std::cos(v);
		const Real y = u < Constants::PI
			? 16 * std::sin(u) + r * std::sin(u / 2) * std::cos(v)
			: 16 * std::sin(u);
		RequirePoint(klein(u, v), VectorN<Real, 3>{2 * x, 2 * y, 2 * r * std::sin(v)});
	}
	SurfaceCases::DiniSurface dini(2, Real(1) / 3);
	const Real u = Constants::PI / 2, v = Constants::PI / 3;
	RequirePoint(dini(u, v), VectorN<Real, 3>{
		2 * std::cos(u) * std::sin(v), 2 * std::sin(u) * std::sin(v),
		2 * (std::cos(v) + std::log(std::tan(v / 2))) + u / 3});
	REQUIRE(dini.getMaxW() == 2);
	REQUIRE(dini.getMinU() == SurfaceCases::DiniSurface().getMinU());
}

TEST_CASE("Existing surface fixture parameterizations remain distinct and usable", "[showcase][surfaces][testbeds]")
{
	const Real u = Real(1) / 2, v = Real(1) / 3;
	RequirePoint(SurfaceCases::MobiusStrip()(u, v), TestBeds::MobiusStrip(u, 2 * v));
	RequirePoint(SurfaceCases::Helicoid()(u, v), TestBeds::Helicoid(v, u));
	REQUIRE_FALSE(IsNearlyEqual(SurfaceCases::KleinBottle()(u, v)[0], TestBeds::KleinBottle(u, v)[0]));
	REQUIRE_FALSE(IsNearlyEqual(SurfaceCases::DiniSurface()(u, v)[2], TestBeds::DiniSurface(u, v)[2]));
	REQUIRE(TestBeds::ParametricCurvesTestBed::getNumTestCurves() == 13);
	REQUIRE(TestBeds::ParametricSurfaceTestBed::getNumTestSurfaces3() == 18);
}

TEST_CASE("Shared showcase catalogs have stable lifetimes and finite samples", "[showcase][testbeds]")
{
	const auto& curves2 = TestBeds::ParametricCurvesTestBed::getShowcaseCurves2D();
	const auto& curves3 = TestBeds::ParametricCurvesTestBed::getShowcaseCurves3D();
	const auto& surfaces = TestBeds::ParametricSurfaceTestBed::getShowcaseSurfaces();
	REQUIRE(curves2.size() == 5);
	REQUIRE(curves3.size() == 2);
	REQUIRE(surfaces.size() == 6);
	REQUIRE(&curves2 == &TestBeds::ParametricCurvesTestBed::getShowcaseCurves2D());
	REQUIRE(&curves3 == &TestBeds::ParametricCurvesTestBed::getShowcaseCurves3D());
	REQUIRE(&surfaces == &TestBeds::ParametricSurfaceTestBed::getShowcaseSurfaces());
	std::set<std::string> names;
	auto sampleCurves = [&](const auto& curves) {
		for (const auto& entry : curves) {
			REQUIRE(names.insert(entry.name).second);
			REQUIRE(entry.start < entry.end);
			for (int i = 0; i <= 8; ++i) {
				const auto point = entry.curve(entry.start + (entry.end - entry.start) * Real(i) / 8);
				for (int coordinate = 0; coordinate < point.size(); ++coordinate)
					REQUIRE(std::isfinite(point[coordinate]));
			}
		}
	};
	sampleCurves(curves2);
	sampleCurves(curves3);
	for (const auto& entry : surfaces) {
		REQUIRE(names.insert(entry.name).second);
		const auto& surface = entry.surface;
		for (int i = 0; i <= 4; ++i)
			for (int j = 0; j <= 4; ++j) {
				const auto point = surface(
					surface.getMinU() + (surface.getMaxU() - surface.getMinU()) * Real(i) / 4,
					surface.getMinW() + (surface.getMaxW() - surface.getMinW()) * Real(j) / 4);
				for (int coordinate = 0; coordinate < 3; ++coordinate)
					REQUIRE(std::isfinite(point[coordinate]));
			}
	}
}

TEST_CASE("Fixture curves and surfaces retain core numerical geometry behavior", "[showcase][curves][surfaces]")
{
	CurveCases::AstroidCurve astroid(3);
	const auto tangent = astroid.getTangent(Constants::PI / 4);
	REQUIRE(tangent.NormL2() > 0);
	REQUIRE(std::isfinite(Curves::getCurvature2D(astroid, Constants::PI / 4)));
	SurfaceCases::Helicoid helicoid;
	REQUIRE(IsWithinAbs(helicoid.Normal(1, 1).NormL2(), Real(1), Precision::LooseTolerance));
	SurfaceCases::EnneperSurface enneper;
	REQUIRE(IsWithinAbs(enneper.Normal(0, 0).NormL2(), Real(1), Precision::LooseTolerance));
	REQUIRE(std::isfinite(enneper.GaussianCurvature(0, 0)));
}
