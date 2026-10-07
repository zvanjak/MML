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

#include <mml/base/Function.h>
#include <mml/base/DifferentialGeometry/TypedFields.h>
#include <mml/core/Derivation.h>
#include <mml/core/Integration.h>
#include <mml/core/DifferentialGeometry/FieldOperations.h>

namespace MML::DifferentialGeometry
{
	namespace Detail
	{
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
				VectorN<Real, N> point = _curve(t);
				VectorN<Real, N> tangent = _derivativeStep > REAL(0.0)
					? Derivation::NDer4(_curve, t, _derivativeStep)
					: Derivation::NDer4(_curve, t);

				return _form(Point<Real, N, FrameTag>(point))(TangentVector<N, FrameTag>(tangent));
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
			if (numU <= 0 || numW <= 0)
				throw ArgumentError("TwoFormSurfaceMidpointSum - subdivision counts must be positive");
			if (u1 <= u0 || w1 <= w0)
				throw ArgumentError("TwoFormSurfaceMidpointSum - invalid parameter interval");

			Real du = (u1 - u0) / numU;
			Real dw = (w1 - w0) / numW;
			Real h = derivativeStep > REAL(0.0) ? derivativeStep : PrecisionValues<Real>::DerivativeStepSize;
			Real total = REAL(0.0);

			for (int i = 0; i < numU; i++) {
				for (int j = 0; j < numW; j++) {
					Real u = u0 + (i + REAL(0.5)) * du;
					Real w = w0 + (j + REAL(0.5)) * dw;
					VectorN<Real, N> point = surface(u, w);
					VectorN<Real, N> tangentU = (surface(u + h, w) - surface(u - h, w)) / (REAL(2.0) * h);
					VectorN<Real, N> tangentW = (surface(u, w + h) - surface(u, w - h)) / (REAL(2.0) * h);

					total += form(Point<Real, N, FrameTag>(point))(
						TangentVector<N, FrameTag>(tangentU),
						TangentVector<N, FrameTag>(tangentW)) * du * dw;
				}
			}

			return total;
		}
	} // namespace Detail

	template<int N, class FrameTag>
	IntegrationResult IntegrateOneForm(const IFormField<N, 1, FrameTag>& form,
		const IParametricCurve<N>& curve,
		Real t0,
		Real t1,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = REAL(0.0))
	{
		Detail::OneFormLineIntegrand<N, FrameTag> helper(form, curve, derivativeStep);
		return IntegrateTrap(helper, t0, t1, eps);
	}

	template<int N, class FrameTag>
	IntegrationResult IntegrateTwoForm(const IFormField<N, 2, FrameTag>& form,
		const IParametricSurfaceRect<N>& surface,
		Real u0,
		Real u1,
		Real w0,
		Real w1,
		int numU = 40,
		int numW = 40,
		Real derivativeStep = REAL(0.0))
	{
		Real value = Detail::TwoFormSurfaceMidpointSum<N, FrameTag>(form, surface, u0, u1, w0, w1, numU, numW, derivativeStep);
		return IntegrationResult(value, REAL(0.0), numU * numW, true);
	}

	template<class FrameTag>
	IntegrationResult IntegrateTwoFormOverRectangle(const IFormField<2, 2, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		int numX = 80,
		int numY = 80)
	{
		if (numX <= 0 || numY <= 0)
			throw ArgumentError("IntegrateTwoFormOverRectangle - subdivision counts must be positive");
		if (x1 <= x0 || y1 <= y0)
			throw ArgumentError("IntegrateTwoFormOverRectangle - invalid rectangle bounds");

		Real dx = (x1 - x0) / numX;
		Real dy = (y1 - y0) / numY;
		TangentVector<2, FrameTag> ex{ REAL(1.0), REAL(0.0) };
		TangentVector<2, FrameTag> ey{ REAL(0.0), REAL(1.0) };
		Real total = REAL(0.0);

		for (int i = 0; i < numX; i++) {
			for (int j = 0; j < numY; j++) {
				Point<Real, 2, FrameTag> point{ x0 + (i + REAL(0.5)) * dx, y0 + (j + REAL(0.5)) * dy };
				total += form(point)(ex, ey) * dx * dy;
			}
		}

		return IntegrationResult(total, REAL(0.0), numX * numY, true);
	}

	template<class FrameTag>
	IntegrationResult IntegrateOneFormAroundRectangle(const IFormField<2, 1, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = REAL(0.0))
	{
		ParametricCurveFromStdFunc<2> bottom(x0, x1, [=](Real t) { return VectorN<Real, 2>{ t, y0 }; });
		ParametricCurveFromStdFunc<2> right(y0, y1, [=](Real t) { return VectorN<Real, 2>{ x1, t }; });
		ParametricCurveFromStdFunc<2> top(x1, x0, [=](Real t) { return VectorN<Real, 2>{ t, y1 }; });
		ParametricCurveFromStdFunc<2> left(y1, y0, [=](Real t) { return VectorN<Real, 2>{ x0, t }; });

		Real value = IntegrateOneForm(form, bottom, x0, x1, eps, derivativeStep).value
			+ IntegrateOneForm(form, right, y0, y1, eps, derivativeStep).value
			+ IntegrateOneForm(form, top, x1, x0, eps, derivativeStep).value
			+ IntegrateOneForm(form, left, y1, y0, eps, derivativeStep).value;

		return IntegrationResult(value, REAL(0.0), 4, true);
	}

	template<class FrameTag>
	Real GreenTheoremResidual(const IFormField<2, 1, FrameTag>& form,
		Real x0,
		Real x1,
		Real y0,
		Real y1,
		int numX = 80,
		int numY = 80,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = REAL(0.0))
	{
		auto exteriorDerivative = ExteriorDerivative(form, derivativeStep);
		Real boundary = IntegrateOneFormAroundRectangle(form, x0, x1, y0, y1, eps, derivativeStep).value;
		Real interior = IntegrateTwoFormOverRectangle(exteriorDerivative, x0, x1, y0, y1, numX, numY).value;
		return boundary - interior;
	}

	template<class FrameTag>
	Real StokesResidual(const IFormField<3, 1, FrameTag>& form,
		const IParametricSurfaceRect<3>& surface,
		int numU = 40,
		int numW = 40,
		Real eps = Defaults::LineIntegralPrecision,
		Real derivativeStep = REAL(0.0))
	{
		Real u0 = surface.getMinU();
		Real u1 = surface.getMaxU();
		Real w0 = surface.getMinW();
		Real w1 = surface.getMaxW();

		ParametricCurveFromStdFunc<3> bottom(u0, u1, [&](Real t) { return surface(t, w0); });
		ParametricCurveFromStdFunc<3> right(w0, w1, [&](Real t) { return surface(u1, t); });
		ParametricCurveFromStdFunc<3> top(u1, u0, [&](Real t) { return surface(t, w1); });
		ParametricCurveFromStdFunc<3> left(w1, w0, [&](Real t) { return surface(u0, t); });

		Real boundary = IntegrateOneForm(form, bottom, u0, u1, eps, derivativeStep).value
			+ IntegrateOneForm(form, right, w0, w1, eps, derivativeStep).value
			+ IntegrateOneForm(form, top, u1, u0, eps, derivativeStep).value
			+ IntegrateOneForm(form, left, w1, w0, eps, derivativeStep).value;

		auto exteriorDerivative = ExteriorDerivative(form, derivativeStep);
		Real interior = IntegrateTwoForm(exteriorDerivative, surface, u0, u1, w0, w1, numU, numW, derivativeStep).value;
		return boundary - interior;
	}
} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_FORM_INTEGRATION_H