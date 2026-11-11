#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/DifferentialGeometry/FormIntegration.h>
#endif

using namespace MML;
using namespace MML::DifferentialGeometry;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Core::DifferentialGeometryTests
{
	class RotationOneForm2D : public IFormField<2, 1, Cartesian2>
	{
	public:
		Form1<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& point) const override
		{
			Form1<2, Cartesian2> form;
			form.Component(0) = -static_cast<MML::Real>(0.5) * point[1];
			form.Component(1) = static_cast<MML::Real>(0.5) * point[0];
			return form;
		}
	};

	class RotationOneForm3D : public IFormField<3, 1, Cartesian3>
	{
	public:
		Form1<3, Cartesian3> operator()(const Point<Real, 3, Cartesian3>& point) const override
		{
			Form1<3, Cartesian3> form;
			form.Component(0) = -static_cast<MML::Real>(0.5) * point[1];
			form.Component(1) = static_cast<MML::Real>(0.5) * point[0];
			form.Component(2) = static_cast<MML::Real>(0.0);
			return form;
		}
	};

	class UnitSquareSurface : public IParametricSurfaceRect<3>
	{
	public:
		VectorN<Real, 3> operator()(Real u, Real w) const override { return {u, w, static_cast<MML::Real>(0.0)}; }
		Real getMinU() const override { return static_cast<MML::Real>(0.0); }
		Real getMaxU() const override { return static_cast<MML::Real>(1.0); }
		Real getMinW() const override { return static_cast<MML::Real>(0.0); }
		Real getMaxW() const override { return static_cast<MML::Real>(1.0); }
	};

	TEST_CASE("IntegrateOneForm integrates a rotation form around the unit circle", "[DifferentialGeometry][FormIntegration]")
	{
		RotationOneForm2D form;
		ParametricCurveFromStdFunc<2> circle(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI,
			[](Real t) { return VectorN<Real, 2>{ std::cos(t), std::sin(t) }; });

		IntegrationResult result = IntegrateOneFormDetailed(form, circle, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0) * Constants::PI);

		REQUIRE(result.converged);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI, TOL(1e-5, 1e-3)));
	}

	TEST_CASE("GreenTheoremResidual vanishes for a square and constant exterior derivative", "[DifferentialGeometry][FormIntegration][Green]")
	{
		RotationOneForm2D form;
		Real residual = GreenTheoremResidual(form, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 60, 60);

		REQUIRE_THAT(residual, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-3)));
	}

	TEST_CASE("IntegrateTwoForm and StokesResidual agree on a flat square patch", "[DifferentialGeometry][FormIntegration][Stokes]")
	{
		RotationOneForm3D form;
		UnitSquareSurface surface;
		auto exteriorDerivative = ExteriorDerivative(form);

		IntegrationResult interior = IntegrateTwoFormDetailed(exteriorDerivative, surface, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 40, 40);
		Real residual = StokesResidual(form, surface, 40, 40);

		REQUIRE(interior.converged);
		REQUIRE_THAT(interior.value, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(residual, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-3)));
	}
}
