#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/DifferentialGeometry/Frames.h>
#include <mml/core/Surfaces.h>
#endif

using namespace MML;
using namespace MML::DifferentialGeometry;
using namespace MML::Surfaces;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::DifferentialGeometryTests
{
	class UnitCircleCurve2D : public IParametricCurve<2>
	{
	public:
		Real getMinT() const override { return static_cast<MML::Real>(0.0); }
		Real getMaxT() const override { return static_cast<MML::Real>(2.0) * Constants::PI; }
		VectorN<Real, 2> operator()(Real t) const override { return {std::cos(t), std::sin(t)}; }
	};

	class SphereEquatorParameterCurve : public IParametricCurve<2>
	{
	public:
		Real getMinT() const override { return static_cast<MML::Real>(0.0); }
		Real getMaxT() const override { return static_cast<MML::Real>(2.0) * Constants::PI; }
		VectorN<Real, 2> operator()(Real t) const override { return {Constants::PI / static_cast<MML::Real>(2.0), t}; }
	};

	TEST_CASE("ComputeFrenetFrame gives unit circle tangent normal and curvature", "[DifferentialGeometry][Frames]")
	{
		UnitCircleCurve2D circle;
		auto frame = ComputeFrenetFrame(circle, static_cast<MML::Real>(0.0));

		REQUIRE_THAT(frame.speed, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-7, 1e-4)));
		REQUIRE_THAT(frame.curvature, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(frame.tangent[0], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-7, 1e-4)));
		REQUIRE_THAT(frame.tangent[1], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-7, 1e-4)));
		REQUIRE_THAT(frame.normal[0], WithinAbs(-static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(frame.normal[1], WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-3)));
	}

	TEST_CASE("ComputeDarbouxFrame classifies sphere equator as a geodesic", "[DifferentialGeometry][Frames][Darboux]")
	{
		Sphere sphere(static_cast<MML::Real>(1.0));
		SphereEquatorParameterCurve equator;
		auto frame = ComputeDarbouxFrame(sphere, equator, static_cast<MML::Real>(0.0));

		REQUIRE_THAT(frame.curvature, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(frame.geodesicCurvature, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(std::abs(frame.normalCurvature), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(frame.tangent[1], WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-7, 1e-4)));
	}
}
