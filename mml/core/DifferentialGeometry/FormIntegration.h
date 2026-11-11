///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/FormIntegration.h                              ///
///  Description: Integration helpers for differential form fields                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_FORM_INTEGRATION_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_FORM_INTEGRATION_H

#include <mml/MMLBase.h>
#include <array>
#include <utility>

#include <mml/base/Function.h>
#include <mml/base/DifferentialGeometry/TypedFields.h>
#include <mml/core/Derivation.h>
#include <mml/core/Derivation/DerivativeProviders.h>
#include <mml/core/Integration.h>
#include <mml/core/DifferentialGeometry/FieldOperations.h>

namespace MML::DifferentialGeometry
{
	namespace Detail
	{
		template<int N>
		void RequireFiniteVector(const VectorN<Real, N>& value)
		{
			for (int i = 0; i < N; ++i)
				IntegrationDetail::RequireFinite(value[i], "Form integration geometry");
		}

		template<class Form>
		Form CheckedForm(Form value)
		{
			for (int i = 0; i < Form::ComponentCount; ++i)
				IntegrationDetail::RequireFinite(value[i], "Form integration component");
			return value;
		}

		inline void ValidateRectangle(Real x0, Real x1, Real y0, Real y1)
		{
			ValidateBounds(x0, x1, "Form integration first interval");
			ValidateBounds(y0, y1, "Form integration second interval");
			if (x1 <= x0 || y1 <= y0)
				throw ArgumentError("Form integration: rectangle bounds must be increasing");
		}

		inline int ValidateGrid(Real x0, Real x1, Real y0, Real y1, int numX, int numY)
		{
			ValidateRectangle(x0, x1, y0, y1);
			if (numX <= 0 || numY <= 0)
				throw ArgumentError("Form integration: subdivision counts must be positive");
			if (numY > std::numeric_limits<int>::max() / numX)
				throw ResourceLimitError("Form integration: midpoint grid exceeds evaluation counter");
			return numX * numY;
		}

		template<int N, class FrameTag>
		Real WeightedTwoForm(const IFormField<N, 2, FrameTag>& form, const VectorN<Real, N>& point,
			const VectorN<Real, N>& first, const VectorN<Real, N>& second, Real halfU, Real halfW, int cells)
		{
			RequireFiniteVector(first);
			RequireFiniteVector(second);
			const auto value = CheckedForm(form(Point<Real, N, FrameTag>(point)));
			IntegrationDetail::Sum pairing;
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
					pairing.Add(IntegrationDetail::Product(
						{value.Component(i, j), first[i], second[j], halfU, halfW, Real(4)}, Real(cells)));
			IntegrationDetail::RequireFinite(pairing.Value(), "Form integration weighted pairing");
			return pairing.Value();
		}

		inline IntegrationResult SumBoundary(const std::array<IntegrationResult, 4>& edges, Real eps)
		{
			auto result = edges[0];
			for (size_t i = 1; i < edges.size(); ++i)
				result = IntegrationDetail::CombineIntegrals(result, edges[i], eps, i + 1 == edges.size());
			return result;
		}

		template<int N, class FrameTag>
		class OneFormLineIntegrand : public IRealFunction
		{
			const IFormField<N, 1, FrameTag>& _form;
			const IParametricCurve<N>& _curve;
			Real _derivativeStep;

		public:
			OneFormLineIntegrand(const IFormField<N, 1, FrameTag>& form,
				const IParametricCurve<N>& curve,
				Real derivativeStep)
				: _form(form), _curve(curve), _derivativeStep(derivativeStep)
			{
			}

			Real operator()(Real t) const override
			{
				auto sample = [&](Real parameter) {
					auto point = _curve(parameter);
					RequireFiniteVector(point);
					return point;
				};
				const auto point = sample(t);
				VectorN<Real, N> tangent;
				if (const auto* provider = dynamic_cast<const ICurveDerivativeProvider<N>*>(&_curve))
					tangent = provider->EvaluateDerivative(t, 1);
				else {
					const Real h = _derivativeStep > 0 ? _derivativeStep : Derivation::ScaleStep(Derivation::NDer4_h, t);
					if (std::isfinite(h) && std::isfinite(t - 2 * h) && std::isfinite(t + 2 * h) &&
						t - 2 * h >= _curve.getMinT() && t + 2 * h <= _curve.getMaxT() &&
						t - h < t && t < t + h)
						tangent = Derivation::NDer4(_curve, t, h);
					else {
						const Real boundedStep = _derivativeStep > 0 ? _derivativeStep
							: Derivation::ScaleStep(Derivation::NDer2_h, t);
						tangent = Derivation::Detail::BoundedFirstDerivative(
							sample, point, t, _curve.getMinT(), _curve.getMaxT(), boundedStep);
					}
				}
				RequireFiniteVector(tangent);
				const auto value = CheckedForm(_form(Point<Real, N, FrameTag>(point)));
				IntegrationDetail::Sum pairing;
				for (int i = 0; i < N; ++i)
					pairing.Add(IntegrationDetail::Product({value.Component(i), tangent[i]}));
				IntegrationDetail::RequireFinite(pairing.Value(), "One-form pairing");
				return pairing.Value();
			}
		};

		template<int N, class FrameTag>
		Real TwoFormSurfaceMidpointSum(const IFormField<N, 2, FrameTag>& form,
			const IParametricSurfaceRect<N>& surface,
			Real u0,
			Real u1,
			Real w0,
			Real w1,
			int numU,
			int numW,
			Real derivativeStep)
		{
			const int cells = ValidateGrid(u0, u1, w0, w1, numU, numW);
			ValidateNonNegativeTolerance(derivativeStep, "Form integration derivative step");
			if (!(surface.getMinU() <= u0 && u1 <= surface.getMaxU() &&
				surface.getMinW() <= w0 && w1 <= surface.getMaxW()))
				throw ArgumentError("Form integration: requested range lies outside surface domain");
			const Real halfU = IntegrationDetail::HalfWidth(u0, u1);
			const Real halfW = IntegrationDetail::HalfWidth(w0, w1);
			IntegrationDetail::Sum total;
			auto sample = [&](Real u, Real w) {
				auto point = surface(u, w);
				RequireFiniteVector(point);
				return point;
			};

			for (int i = 0; i < numU; i++) {
				for (int j = 0; j < numW; j++) {
					const Real u = std::lerp(u0, u1, (Real(i) + Real(0.5)) / numU);
					const Real w = std::lerp(w0, w1, (Real(j) + Real(0.5)) / numW);
					const auto point = sample(u, w);
					const Real stepU = derivativeStep > 0 ? derivativeStep
						: Derivation::ScaleStep(PrecisionValues<Real>::DerivativeStepSize, u);
					const Real stepW = derivativeStep > 0 ? derivativeStep
						: Derivation::ScaleStep(PrecisionValues<Real>::DerivativeStepSize, w);
					const auto tangentU = Derivation::Detail::BoundedFirstDerivative(
						[&](Real t) { return sample(t, w); }, point, u, u0, u1, stepU);
					const auto tangentW = Derivation::Detail::BoundedFirstDerivative(
						[&](Real t) { return sample(u, t); }, point, w, w0, w1, stepW);
					total.Add(WeightedTwoForm(form, point, tangentU, tangentW, halfU, halfW, cells));
				}
			}

			IntegrationDetail::RequireFinite(total.Value(), "Two-form surface integral");
			return total.Value();
		}
	} // namespace Detail

	template<int N, class FrameTag>
	IntegrationResult IntegrateOneFormDetailed(const IFormField<N, 1, FrameTag>& form,
		const IParametricCurve<N>& curve,
		Real t0,
		Real t1,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		ValidateBounds(t0, t1, "One-form integration interval");
		ValidateNonNegativeTolerance(derivativeStep, "One-form derivative step");
		const Real lower = curve.getMinT(), upper = curve.getMaxT();
		if (!(lower <= upper && lower <= std::min(t0, t1) && std::max(t0, t1) <= upper))
			throw ArgumentError("One-form integration: interval lies outside curve domain");
		Detail::OneFormLineIntegrand<N, FrameTag> helper(form, curve, derivativeStep);
		return IntegrateTrapCore(helper, t0, t1, eps);
	}

	template<int N, class FrameTag>
	Real IntegrateOneForm(const IFormField<N, 1, FrameTag>& form,
		const IParametricCurve<N>& curve,
		Real t0,
		Real t1,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		return IntegrateOneFormDetailed<N>(form, curve, t0, t1, eps, derivativeStep).ValueOrThrow("IntegrateOneForm");
	}

	template<int N, class FrameTag>
	IntegrationResult IntegrateTwoFormDetailed(const IFormField<N, 2, FrameTag>& form,
		const IParametricSurfaceRect<N>& surface,
		Real u0,
		Real u1,
		Real w0,
		Real w1,
		int numU = 40,
		int numW = 40,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		Real value = Detail::TwoFormSurfaceMidpointSum<N, FrameTag>(form, surface, u0, u1, w0, w1, numU, numW, derivativeStep);
		return IntegrationResult::FixedRule(value, numU * numW);
	}

	template<int N, class FrameTag>
	Real IntegrateTwoForm(const IFormField<N, 2, FrameTag>& form,
		const IParametricSurfaceRect<N>& surface,
		Real u0,
		Real u1,
		Real w0,
		Real w1,
		int numU = 40,
		int numW = 40,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		return IntegrateTwoFormDetailed<N>(form, surface, u0, u1, w0, w1, numU, numW, derivativeStep).ValueOrThrow("IntegrateTwoForm");
	}

	template<class FrameTag>
	IntegrationResult IntegrateTwoFormOverRectangleDetailed(const IFormField<2, 2, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		int numX = 80,
		int numY = 80)
	{
		const int cells = Detail::ValidateGrid(x0, x1, y0, y1, numX, numY);
		const Real halfX = IntegrationDetail::HalfWidth(x0, x1);
		const Real halfY = IntegrationDetail::HalfWidth(y0, y1);
		const VectorN<Real, 2> ex{1, 0}, ey{0, 1};
		IntegrationDetail::Sum total;

		for (int i = 0; i < numX; i++) {
			for (int j = 0; j < numY; j++) {
				const VectorN<Real, 2> point{
					std::lerp(x0, x1, (Real(i) + Real(0.5)) / numX),
					std::lerp(y0, y1, (Real(j) + Real(0.5)) / numY)};
				total.Add(Detail::WeightedTwoForm(form, point, ex, ey, halfX, halfY, cells));
			}
		}

		IntegrationDetail::RequireFinite(total.Value(), "Two-form rectangle integral");
		return IntegrationResult::FixedRule(total.Value(), cells);
	}

	template<class FrameTag>
	Real IntegrateTwoFormOverRectangle(const IFormField<2, 2, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		int numX = 80,
		int numY = 80)
	{
		return IntegrateTwoFormOverRectangleDetailed(form, x0, x1, y0, y1, numX, numY).ValueOrThrow("IntegrateTwoFormOverRectangle");
	}

	template<class FrameTag>
	IntegrationResult IntegrateOneFormAroundRectangleDetailed(const IFormField<2, 1, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		Detail::ValidateRectangle(x0, x1, y0, y1);
		ParametricCurveFromStdFunc<2> bottom(x0, x1, [=](Real t) { return VectorN<Real, 2>{ t, y0 }; });
		ParametricCurveFromStdFunc<2> right(y0, y1, [=](Real t) { return VectorN<Real, 2>{ x1, t }; });
		ParametricCurveFromStdFunc<2> top(x0, x1, [=](Real t) { return VectorN<Real, 2>{ t, y1 }; });
		ParametricCurveFromStdFunc<2> left(y0, y1, [=](Real t) { return VectorN<Real, 2>{ x0, t }; });

		return Detail::SumBoundary({
			IntegrateOneFormDetailed(form, bottom, x0, x1, eps, derivativeStep),
			IntegrateOneFormDetailed(form, right, y0, y1, eps, derivativeStep),
			IntegrateOneFormDetailed(form, top, x1, x0, eps, derivativeStep),
			IntegrateOneFormDetailed(form, left, y1, y0, eps, derivativeStep)}, eps);
	}

	template<class FrameTag>
	Real IntegrateOneFormAroundRectangle(const IFormField<2, 1, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		return IntegrateOneFormAroundRectangleDetailed(form, x0, x1, y0, y1, eps, derivativeStep).ValueOrThrow("IntegrateOneFormAroundRectangle");
	}

	struct FormTheoremResult : IntegrationResult
	{
		IntegrationResult boundary;
		IntegrationResult interior;

		FormTheoremResult(IntegrationResult boundaryResult, IntegrationResult interiorResult, Real eps)
			: boundary(std::move(boundaryResult)), interior(std::move(interiorResult))
		{
			auto negativeInterior = interior;
			negativeInterior.value = -negativeInterior.value;
			static_cast<IntegrationResult&>(*this) = IntegrationDetail::CombineIntegrals(boundary, negativeInterior, eps);
		}
	};

	template<class FrameTag>
	FormTheoremResult GreenTheoremResidualDetailed(const IFormField<2, 1, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		int numX = 80,
		int numY = 80,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		Detail::ValidateGrid(x0, x1, y0, y1, numX, numY);
		auto exteriorDerivative = ExteriorDerivative(form, derivativeStep);
		auto boundary = IntegrateOneFormAroundRectangleDetailed(form, x0, x1, y0, y1, eps, derivativeStep);
		auto interior = IntegrateTwoFormOverRectangleDetailed(exteriorDerivative, x0, x1, y0, y1, numX, numY);
		return FormTheoremResult(std::move(boundary), std::move(interior), eps);
	}

	template<class FrameTag>
	Real GreenTheoremResidual(const IFormField<2, 1, FrameTag>& form,
		Real x0, Real x1, Real y0, Real y1, int numX = 80, int numY = 80,
		Real eps = Defaults::LineIntegralPrecision, Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		return GreenTheoremResidualDetailed(form, x0, x1, y0, y1, numX, numY, eps, derivativeStep)
			.ValueOrThrow("GreenTheoremResidual");
	}

	template<class FrameTag>
	FormTheoremResult StokesResidualDetailed(const IFormField<3, 1, FrameTag>& form,
		const IParametricSurfaceRect<3>& surface,
		int numU = 40,
		int numW = 40,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		Real u0 = surface.getMinU();
		Real u1 = surface.getMaxU();
		Real w0 = surface.getMinW();
		Real w1 = surface.getMaxW();
		Detail::ValidateGrid(u0, u1, w0, w1, numU, numW);

		ParametricCurveFromStdFunc<3> bottom(u0, u1, [&](Real t) { return surface(t, w0); });
		ParametricCurveFromStdFunc<3> right(w0, w1, [&](Real t) { return surface(u1, t); });
		ParametricCurveFromStdFunc<3> top(u0, u1, [&](Real t) { return surface(t, w1); });
		ParametricCurveFromStdFunc<3> left(w0, w1, [&](Real t) { return surface(u0, t); });

		auto boundary = Detail::SumBoundary({
			IntegrateOneFormDetailed(form, bottom, u0, u1, eps, derivativeStep),
			IntegrateOneFormDetailed(form, right, w0, w1, eps, derivativeStep),
			IntegrateOneFormDetailed(form, top, u1, u0, eps, derivativeStep),
			IntegrateOneFormDetailed(form, left, w1, w0, eps, derivativeStep)}, eps);

		auto exteriorDerivative = ExteriorDerivative(form, derivativeStep);
		auto interior = IntegrateTwoFormDetailed(exteriorDerivative, surface, u0, u1, w0, w1, numU, numW, derivativeStep);
		return FormTheoremResult(std::move(boundary), std::move(interior), eps);
	}

	template<class FrameTag>
	Real StokesResidual(const IFormField<3, 1, FrameTag>& form,
		const IParametricSurfaceRect<3>& surface, int numU = 40, int numW = 40,
		Real eps = Defaults::LineIntegralPrecision, Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		return StokesResidualDetailed(form, surface, numU, numW, eps, derivativeStep).ValueOrThrow("StokesResidual");
	}
} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_FORM_INTEGRATION_H
