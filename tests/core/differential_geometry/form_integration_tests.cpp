#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/DifferentialGeometry/FormIntegration.h>
#include <mml/core/Derivation/DerivativeProviders.h>
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
		VectorN<Real, 3> operator()(Real u, Real w) const override {
			if (u < 0 || u > 1 || w < 0 || w > 1)
				throw DomainError("UnitSquareSurface: outside closed domain");
			return {u, w, 0};
		}
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

		REQUIRE(interior.status == AlgorithmStatus::Completed);
		REQUIRE_FALSE(interior.converged);
		REQUIRE_FALSE(interior.error_estimated);
		REQUIRE_THAT(interior.value, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(residual, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-5, 1e-3)));
	}

	template<int N, class Frame>
	struct QuadraticAreaForm : IFormField<N, 2, Frame> {
		mutable int calls = 0;
		DifferentialForm<Real, N, 2, Frame> operator()(const Point<Real, N, Frame>& p) const override {
			++calls;
			DifferentialForm<Real, N, 2, Frame> result;
			result.SetAlternatingComponent(p[0] * p[0], 0, 1);
			return result;
		}
	};

	struct ExponentialOneForm : IFormField<2, 1, Cartesian2> {
		mutable int calls = 0;
		Form1<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& p) const override {
			++calls;
			return {std::exp(p[0]), 0};
		}
	};

	struct ClosedCurve : IParametricCurve<2> {
		mutable int calls = 0;
		Real getMinT() const override { return 0; }
		Real getMaxT() const override { return 1; }
		VectorN<Real, 2> operator()(Real t) const override {
			if (t < 0 || t > 1) throw DomainError("ClosedCurve: outside closed domain");
			++calls;
			return {t * t, 0};
		}
	};

	struct ProviderCurve : ClosedCurve, ICurveDerivativeProvider<2> {
		mutable int derivativeCalls = 0;
		VectorN<Real, 2> EvaluateDerivative(Real t, int order) const override {
			if (t < 0 || t > 1 || order != 1) throw DomainError("ProviderCurve: invalid derivative request");
			++derivativeCalls;
			return {2 * t, 0};
		}
	};

	struct ContextRestore {
		AlgorithmContext saved = AlgorithmContext::Get();
		~ContextRestore() { AlgorithmContext::Get() = saved; }
	};

	TEST_CASE("Form midpoint rules report unchecked accuracy and real sample work", "[FormIntegration][form-contract]")
	{
		QuadraticAreaForm<2, Cartesian2> rectangleForm;
		const auto coarse = IntegrateTwoFormOverRectangleDetailed(rectangleForm, 0, 1, 0, 1, 1, 1);
		const auto fine = IntegrateTwoFormOverRectangleDetailed(rectangleForm, 0, 1, 0, 1, 8, 4);
		CHECK(coarse.value == Real(0.25));
		CHECK_THAT(Real(1) / 3 - fine.value, WithinAbs((Real(1) / 3 - coarse.value) / 64, Precision::EqualityTolerance));
		for (const auto* result : {&coarse, &fine}) {
			CHECK(result->status == AlgorithmStatus::Completed);
			CHECK(result->IsUsable());
			CHECK_FALSE(result->converged);
			CHECK_FALSE(result->error_estimated);
			CHECK(std::isnan(result->error_estimate));
			CHECK(result->iterations == 1);
		}
		CHECK(coarse.function_evaluations == 1);
		CHECK(fine.function_evaluations == 32);
		CHECK(rectangleForm.calls == 33);

		QuadraticAreaForm<3, Cartesian3> surfaceForm;
		UnitSquareSurface surface;
		const auto patch = IntegrateTwoFormDetailed(surfaceForm, surface, 0, 1, 0, 1, 1, 1, 2);
		CHECK_THAT(patch.value, WithinAbs(coarse.value, Precision::EqualityTolerance));
		CHECK(patch.status == AlgorithmStatus::Completed);
		CHECK_FALSE(patch.error_estimated);
		CHECK(patch.function_evaluations == surfaceForm.calls);
	}

	TEST_CASE("Form curve integration respects closed domains and analytic providers", "[FormIntegration][form-contract]")
	{
		ExponentialOneForm form;
		ClosedCurve numerical;
		const auto result = IntegrateOneFormDetailed(form, numerical, 0, 1, Precision::DefaultTolerance);
		CHECK(result.IsUsable());
		CHECK_THAT(result.value, WithinAbs(std::exp(Real(1)) - 1, Precision::DefaultTolerance));
		CHECK(result.function_evaluations == form.calls);

		ProviderCurve provider;
		form.calls = 0;
		const auto analytic = IntegrateOneFormDetailed(form, provider, 0, 1, Precision::DefaultTolerance);
		CHECK(analytic.IsUsable());
		CHECK_THAT(analytic.value, WithinAbs(result.value, Precision::DefaultTolerance));
		CHECK(provider.derivativeCalls == analytic.function_evaluations);
		CHECK(provider.calls == analytic.function_evaluations);
		CHECK(form.calls == analytic.function_evaluations);
		const auto reversed = IntegrateOneFormDetailed(form, provider, 1, 0, Precision::DefaultTolerance, 2);
		CHECK(reversed.IsUsable());
		CHECK_THAT(reversed.value, WithinAbs(-analytic.value, Precision::DefaultTolerance));
	}

	TEST_CASE("Form pairing and reparameterization preserve oriented integrals", "[FormIntegration][form-contract]")
	{
		struct CustomFrame {};
		struct ExactDifferential : IFormField<2, 1, CustomFrame> {
			Real scale = 1;
			Form1<2, CustomFrame> operator()(const Point<Real, 2, CustomFrame>& p) const override {
				return {scale * p[1], scale * p[0]};
			}
		} form;
		ParametricCurveFromStdFunc<2> parabola(0, 1, [](Real t) { return VectorN<Real, 2>{t, t*t}; });
		ParametricCurveFromStdFunc<2> reparameterized(0, 1, [](Real t) { return VectorN<Real, 2>{t*t, t*t*t*t}; });
		// Endpoint stencils are second order with h=O(eps^(1/3)), hence
		// O(eps^(2/3)) derivative error in addition to quadrature error.
		const Real derivativeBound = std::max(Precision::DefaultTolerance,
			Real(2) * std::pow(Precision::Eps, Real(2) / 3));
		CHECK_THAT(IntegrateOneForm(form, parabola, 0, 1, Precision::DefaultTolerance), WithinAbs(Real(1), derivativeBound));
		CHECK_THAT(IntegrateOneForm(form, reparameterized, 0, 1, Precision::DefaultTolerance), WithinAbs(Real(1), derivativeBound));
		form.scale = 2;
		CHECK_THAT(IntegrateOneForm(form, parabola, 1, 0, Precision::DefaultTolerance), WithinAbs(Real(-2), 2 * derivativeBound));
	}

	TEST_CASE("Rectangle form boundaries retain failed edge diagnostics and work", "[FormIntegration][form-contract]")
	{
		ContextRestore restore;
		AlgorithmContext::Get().trapezoidIntegrationMaxSteps = 7;
		ExponentialOneForm form;
		const auto result = IntegrateOneFormAroundRectangleDetailed(form, 0, 1, 0, 1, Precision::EqualityTolerance);
		CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK_FALSE(result.IsUsable());
		CHECK_FALSE(result.converged);
		CHECK(result.error_estimated);
		CHECK(result.error_estimate > 0);
		CHECK(result.iterations == 28);
		CHECK(result.function_evaluations == 260);
		CHECK(result.function_evaluations == form.calls);
		CHECK_FALSE(result.error_message.empty());
		CHECK_THROWS_AS(IntegrateOneFormAroundRectangle(form, 0, 1, 0, 1, Precision::EqualityTolerance), ConvergenceError);
		CHECK_THROWS_AS(GreenTheoremResidual(form, 0, 1, 0, 1, 2, 2, Precision::EqualityTolerance), ConvergenceError);
	}

	TEST_CASE("Boundary accuracy is checked after all four oriented edges are accumulated", "[FormIntegration][form-contract]")
	{
		struct Form : IFormField<2, 1, Cartesian2> {
			Form1<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& p) const override {
				return {(-1 + 100 * p[1]) * std::exp(p[0]), std::exp(p[1])};
			}
		} form;
		const auto result = IntegrateOneFormAroundRectangleDetailed(form, 0, 1, 0, 1, Precision::DefaultTolerance);
		CAPTURE(result.value, result.error_estimate, result.error_message);
		CHECK(result.IsUsable());
		CHECK(result.converged);
		CHECK(result.error_estimated);
		// Four-edge aggregation permits the existing 50*eps roundoff floor
		// on the sum of absolute edge integrals, not only on their difference.
		const Real absoluteEdges = 102 * (std::exp(Real(1)) - 1);
		const Real comparisonBound = std::max(100 * Precision::DefaultTolerance,
			Precision::IntegrationRoundoffFloorScale * absoluteEdges);
		CHECK_THAT(result.value, WithinAbs(-100 * (std::exp(Real(1)) - 1), comparisonBound));
		CHECK(result.error_estimate > 0);
	}

	TEST_CASE("Form controls reject unsafe geometry and grid budgets before callbacks", "[FormIntegration][form-contract][form-budget]")
	{
		QuadraticAreaForm<2, Cartesian2> form;
		const int huge = std::numeric_limits<int>::max();
		CHECK_THROWS_AS(IntegrateTwoFormOverRectangleDetailed(form, 0, 1, 0, 1, huge, 2), ResourceLimitError);
		CHECK(form.calls == 0);
		CHECK_THROWS(IntegrateTwoFormOverRectangleDetailed(form, 0, Constants::PosInf, 0, 1, 1, 1));
		CHECK_THROWS(IntegrateTwoFormOverRectangleDetailed(form, 0, 1, 0, 1, 0, 1));
		CHECK_THROWS_AS(GreenTheoremResidualDetailed(RotationOneForm2D{}, 0, 1, 0, 1, huge, 2), ResourceLimitError);
		CHECK(form.calls == 0);
		QuadraticAreaForm<3, Cartesian3> surfaceForm;
		UnitSquareSurface surface;
		CHECK_THROWS_AS(IntegrateTwoFormDetailed(surfaceForm, surface, 0, 1, 0, 1, huge, 2), ResourceLimitError);
		CHECK_THROWS(IntegrateTwoFormDetailed(surfaceForm, surface, 0, 1, 0, 1, 1, 1, -1));
		CHECK_THROWS(IntegrateTwoFormDetailed(surfaceForm, surface, -1, 1, 0, 1, 1, 1));
		CHECK_THROWS_AS(StokesResidualDetailed(RotationOneForm3D{}, surface, huge, 2), ResourceLimitError);
		CHECK(surfaceForm.calls == 0);
		ExponentialOneForm oneForm;
		ClosedCurve curve;
		CHECK_THROWS(IntegrateOneFormDetailed(oneForm, curve, -1, 1));
		CHECK_THROWS(IntegrateOneFormDetailed(oneForm, curve, 0, 1, Defaults::LineIntegralPrecision, std::numeric_limits<Real>::quiet_NaN()));
		CHECK(curve.calls == 0);
		CHECK(oneForm.calls == 0);
		struct TinyCurve : ClosedCurve {
			Real getMaxT() const override { return std::nextafter(Real(0), Real(1)); }
		} tiny;
		CHECK_THROWS_AS(IntegrateOneFormDetailed(oneForm, tiny, 0, tiny.getMaxT()), NumericalMethodError);
		struct StopAtFirstSample : IFormField<2, 2, Cartesian2> {
			mutable int calls = 0;
			Form2<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>&) const override {
				++calls;
				throw NumericalMethodError("Accepted exact counter boundary; stop before huge loop");
			}
		} stop;
		CHECK_THROWS_AS(IntegrateTwoFormOverRectangleDetailed(stop, 0, 1, 0, 1, huge, 1), NumericalMethodError);
		CHECK(stop.calls == 1);
	}

	TEST_CASE("Two-form pullbacks preserve orientation in a non-Cartesian four-dimensional frame", "[FormIntegration][form-contract]")
	{
		struct PatchFrame {};
		struct Patch : IParametricSurfaceRect<4> {
			bool swapped = false;
			Real getMinU() const override { return 0; }
			Real getMaxU() const override { return 1; }
			Real getMinW() const override { return 0; }
			Real getMaxW() const override { return 1; }
			VectorN<Real, 4> operator()(Real u, Real w) const override {
				if (u < 0 || u > 1 || w < 0 || w > 1) throw DomainError("Patch: outside domain");
				return swapped ? VectorN<Real, 4>{2*w, 3*u, w, u} : VectorN<Real, 4>{2*u, 3*w, u, w};
			}
		} patch;
		QuadraticAreaForm<4, PatchFrame> form;
		const auto forward = IntegrateTwoFormDetailed(form, patch, 0, 1, 0, 1, 4, 4);
		patch.swapped = true;
		const auto reverse = IntegrateTwoFormDetailed(form, patch, 0, 1, 0, 1, 4, 4);
		const Real midpointOracle = Real(63) / 8;
		// The linear patch has zero truncation error, but its two numerical
		// tangent differences contribute O(eps/h) roundoff to the pullback.
		const Real tangentBound = std::max(Precision::DefaultTolerance,
			Real(2) * Precision::Eps / PrecisionValues<Real>::DerivativeStepSize * midpointOracle);
		CHECK_THAT(forward.value, WithinAbs(midpointOracle, tangentBound));
		CHECK_THAT(reverse.value, WithinAbs(-midpointOracle, tangentBound));
		CHECK(forward.status == AlgorithmStatus::Completed);
		CHECK(reverse.function_evaluations == 16);
		CHECK(form.calls == 32);
	}

	TEST_CASE("Form integration rejects nonfinite samples and weighted totals", "[FormIntegration][form-contract]")
	{
		struct ConstantAreaForm : IFormField<2, 2, Cartesian2> {
			Real value = Constants::PosInf;
			Form2<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>&) const override {
				Form2<2, Cartesian2> result;
				result.SetAlternatingComponent(value, 0, 1);
				return result;
			}
		} form;
		CHECK_THROWS_AS(IntegrateTwoFormOverRectangleDetailed(form, 0, 1, 0, 1, 1, 1), NumericalMethodError);
		form.value = std::numeric_limits<Real>::max();
		CHECK_THROWS_AS(IntegrateTwoFormOverRectangleDetailed(form, 0, 2, 0, 2, 1, 1), NumericalMethodError);
		ParametricCurveFromStdFunc<2> badCurve(0, 1, [](Real) { return VectorN<Real, 2>{std::numeric_limits<Real>::quiet_NaN(), 0}; });
		ExponentialOneForm oneForm;
		CHECK_THROWS_AS(IntegrateOneFormDetailed(oneForm, badCurve, 0, 1), NumericalMethodError);
		struct BadProvider : ProviderCurve {
			VectorN<Real, 2> EvaluateDerivative(Real, int) const override { return {Constants::PosInf, 0}; }
		} badProvider;
		CHECK_THROWS_AS(IntegrateOneFormDetailed(oneForm, badProvider, 0, 1), NumericalMethodError);
		form.value = 1 / std::numeric_limits<Real>::max();
		const auto wide = IntegrateTwoFormOverRectangleDetailed(form, -std::numeric_limits<Real>::max(),
			std::numeric_limits<Real>::max(), 0, 1, 1, 1);
		CHECK_THAT(wide.value, WithinAbs(Real(2), Precision::EqualityTolerance));
	}

	TEST_CASE("Theorem residual results own component diagnostics without claiming an accuracy check", "[FormIntegration][form-contract]")
	{
		const auto green = [] {
			RotationOneForm2D form;
			return GreenTheoremResidualDetailed(form, 0, 1, 0, 1, 4, 3);
		}();
		CHECK(green.status == AlgorithmStatus::Completed);
		CHECK(green.IsUsable());
		CHECK_FALSE(green.converged);
		CHECK_FALSE(green.error_estimated);
		CHECK(std::isnan(green.error_estimate));
		CHECK(green.boundary.converged);
		CHECK(green.boundary.error_estimated);
		CHECK(green.interior.status == AlgorithmStatus::Completed);
		CHECK(green.interior.function_evaluations == 12);
		CHECK(green.function_evaluations == green.boundary.function_evaluations + 12);
		CHECK(green.iterations == green.boundary.iterations + 1);
		CHECK_THAT(green.value, WithinAbs(Real(0), Precision::DefaultTolerance));

		RotationOneForm3D form;
		UnitSquareSurface surface;
		const auto stokes = StokesResidualDetailed(form, surface, 4, 3);
		CHECK(stokes.status == AlgorithmStatus::Completed);
		CHECK(stokes.boundary.IsUsable());
		CHECK_FALSE(stokes.interior.error_estimated);
		CHECK(stokes.function_evaluations == stokes.boundary.function_evaluations + 12);
		CHECK_THAT(stokes.value, WithinAbs(Real(0), Precision::DefaultTolerance));
	}

	TEST_CASE("Theorem residuals retain failed boundary uncertainty even with an unchecked interior", "[FormIntegration][form-contract]")
	{
		ContextRestore restore;
		AlgorithmContext::Get().trapezoidIntegrationMaxSteps = 7;
		ExponentialOneForm form;
		const auto green = GreenTheoremResidualDetailed(form, 0, 1, 0, 1, 2, 2, Precision::EqualityTolerance);
		CHECK(green.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(green.boundary.status == green.status);
		CHECK(green.boundary.error_estimated);
		CHECK(green.boundary.error_estimate > 0);
		CHECK(green.interior.status == AlgorithmStatus::Completed);
		CHECK_FALSE(green.error_estimated);
		CHECK(std::isnan(green.error_estimate));
		CHECK(green.function_evaluations == 264);
		CHECK(green.iterations == 29);
		CHECK_THROWS_AS(green.ValueOrThrow("Green"), ConvergenceError);

		struct Exponential3D : IFormField<3, 1, Cartesian3> {
			Form1<3, Cartesian3> operator()(const Point<Real, 3, Cartesian3>& p) const override {
				return {std::exp(p[0]), 0, 0};
			}
		} form3;
		UnitSquareSurface surface;
		const auto stokes = StokesResidualDetailed(form3, surface, 2, 2, Precision::EqualityTolerance);
		CHECK(stokes.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(stokes.boundary.error_estimate > 0);
		CHECK_FALSE(stokes.error_estimated);
		CHECK(stokes.function_evaluations == 264);
		CHECK_THROWS_AS(StokesResidual(form3, surface, 2, 2, Precision::EqualityTolerance), ConvergenceError);
	}
}
