#include <test_beds/parametric_curves_test_bed.h>
#include <test_beds/parametric_surfaces_test_bed.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace MML;

void Docs_Demo_ShowcaseGeometry()
{
	TestBeds::ShowcaseCurves::ButterflyCurve butterfly(2);
	TestBeds::ShowcaseSurfaces::EnneperSurface enneper(2);
	std::cout << "Shared fixture butterfly: " << butterfly(Real(0)) << '\n'
	          << "Shared fixture Enneper: " << enneper(Real(1), Real(1)) << '\n';

	const auto& planar = TestBeds::ParametricCurvesTestBed::getShowcaseCurves2D();
	const auto& spatial = TestBeds::ParametricCurvesTestBed::getShowcaseCurves3D();
	const auto& surfaces = TestBeds::ParametricSurfaceTestBed::getShowcaseSurfaces();
	auto sampleCurves = [](const auto& curves) {
		for (const auto& entry : curves) {
			const Real t = (entry.start + entry.end) / 2;
			const auto point = entry.curve(t);
			for (int i = 0; i < point.size(); ++i)
				if (!std::isfinite(point[i]))
					throw std::runtime_error(std::string("Nonfinite showcase curve sample: ") + entry.name);
			std::cout << entry.name << ": " << point << '\n';
		}
	};
	sampleCurves(planar);
	sampleCurves(spatial);
	for (const auto& entry : surfaces) {
		const auto& surface = entry.surface;
		const Real u = (surface.getMinU() + surface.getMaxU()) / 2;
		const Real v = (surface.getMinW() + surface.getMaxW()) / 2;
		const auto point = surface(u, v);
		for (int i = 0; i < 3; ++i)
			if (!std::isfinite(point[i]))
				throw std::runtime_error(std::string("Nonfinite showcase surface sample: ") + entry.name);
		std::cout << entry.name << ": " << point << '\n';
	}

	TestBeds::ShowcaseCurves::DeltoidCurve deltoid(2);
	TestBeds::ShowcaseCurves::AstroidCurve astroid(3);
	TestBeds::ShowcaseCurves::EpitrochoidCurve epitrochoid(1, Real(1) / 2, 3);
	TestBeds::ShowcaseCurves::LemniscateCurve lemniscate;
	TestBeds::ShowcaseCurves::ToroidalSpiralCurve toroidal(5, 2);
	TestBeds::ShowcaseSurfaces::Helicoid helicoid(Real(2) / 5);
	std::cout << "Parameterized fixtures: " << deltoid(0) << ' ' << astroid(Constants::PI / 4)
	          << ' ' << epitrochoid(0) << ' ' << lemniscate(0) << ' ' << toroidal(0)
	          << ' ' << helicoid(1, 1) << '\n';
	const Real u = Real(1) / 2, v = Real(1) / 3;
	const auto fullWidth = TestBeds::ShowcaseSurfaces::MobiusStrip()(u, v);
	const auto halfWidth = TestBeds::MobiusStrip(u, 2 * v);
	if (!IsWithinAbs((fullWidth - halfWidth).NormL2(), Real(0), Precision::EqualityTolerance))
		throw std::runtime_error("Showcase and legacy Mobius width conventions were confused");
	std::cout << "Fixture catalogs: " << planar.size() << " planar curves, " << spatial.size()
	          << " spatial curves, " << surfaces.size() << " surfaces; no GUI launched.\n";
}
