///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationParametricCurve.h                                         ///
///  Description: Derivatives of parametric curves                                    ///
///               Tangent, normal, binormal, curvature, torsion calculations          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_PARAMETRIC_CURVE_H
#define MML_DERIVATION_PARAMETRIC_CURVE_H

#include <mml/MMLBase.h>

#include "DerivationBase.h"
#include "FirstDerivativeStencil.h"

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::One>(
				[&](int offset) { return f(t + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}

		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer1(f, t, ScaleStep(NDer1_h, t), error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>;
			auto result = Detail::EvaluateFirstDerivativeStencilWithOffsets<Detail::FirstDerivativeOrder::Two>(
				[&](int offset) { return f(t + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr,
				Stencil::value_offsets, Stencil::error_offsets,
				[](const auto& at, const auto& norm, Real step) {
					auto diff = at(1) - at(-1);
					return Constants::Eps * norm((at(1) + at(-1)) / (REAL(2.0) * step))
					     + norm((at(2) - at(-2)) / REAL(2.0) - diff) / (REAL(6.0) * step);
				});
			if (error) *error = result.error;
			return result.value;
		}

		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer2(f, t, ScaleStep(NDer2_h, t), error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Four>;
			auto result = Detail::EvaluateFirstDerivativeStencilWithOffsets<Detail::FirstDerivativeOrder::Four>(
				[&](int offset) { return f(t + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr,
				Stencil::value_offsets, Stencil::error_offsets,
				[](const auto& at, const auto& norm, Real step) {
					Real truncation = norm(at(3) - at(-3)) / REAL(2.0) + REAL(2.0) * norm(at(-2) - at(2))
					                + REAL(5.0) * norm(at(1) - at(-1)) / REAL(2.0);
					return std::abs(truncation) / (REAL(30.0) * step)
					     + Constants::Eps * (norm(at(2)) + norm(at(-2)) + REAL(8.0) * (norm(at(-1)) + norm(at(1)))) / (REAL(12.0) * step);
				});
			if (error) *error = result.error;
			return result.value;
		}

		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer4(f, t, ScaleStep(NDer4_h, t), error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::Six>(
				[&](int offset) { return f(t + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}

		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer6(f, t, ScaleStep(NDer6_h, t), error);
		}	
		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::Eight>(
				[&](int offset) { return f(t + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}

		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer8(f, t, ScaleStep(NDer8_h, t), error);
		}

		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********  NOTE: Direct finite difference formulas (not derivatives of derivatives!)                        ********/
		/********        NSecDer2: 3 function evals (O(h²) accuracy)                                                ********/
		/********        NSecDer4: 5 function evals (O(h⁴) accuracy)                                                ********/
		/********************************************************************************************************************/
		
		// f''(t) ≈ [f(t-h) - 2f(t) + f(t+h)] / h²
		// Second-order accurate (O(h²)), 3 function evaluations
		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> y0 = f(t);
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			
			Real h2 = h * h;
			VectorN<Real, N> result = (ymh - 2.0 * y0 + yh) / h2;
			
			if (error)
			{
				// Error estimate using 4th derivative approximation
				VectorN<Real, N> y2h = f(t + 2*h);
				VectorN<Real, N> ym2h = f(t - 2*h);
				Real f4_approx = (ym2h - 4.0*ymh + 6.0*y0 - 4.0*yh + y2h).NormL2() / h2;
				
				*error = f4_approx * h2 / 12.0 + 
				         Constants::Eps * (ymh.NormL2() + 2*y0.NormL2() + yh.NormL2()) / h2;
			}
			
			return result;
		}
		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer2(f, t, ScaleStep(NDer2_h, t), error);
		}

		// f''(t) ≈ [-f(t-2h) + 16f(t-h) - 30f(t) + 16f(t+h) - f(t+2h)] / (12h²)
		// Fourth-order accurate (O(h⁴)), 5 function evaluations
		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> y0 = f(t);
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y2h = f(t + 2*h);
			VectorN<Real, N> ym2h = f(t - 2*h);
			
			Real h2 = h * h;
			VectorN<Real, N> result = (-ym2h + 16.0*ymh - 30.0*y0 + 16.0*yh - y2h) / (12.0 * h2);
			
			if (error)
			{
				// Error estimate using 6th derivative approximation
				VectorN<Real, N> y3h = f(t + 3*h);
				VectorN<Real, N> ym3h = f(t - 3*h);
				Real f6_approx = (ym3h - 6.0*ym2h + 15.0*ymh - 20.0*y0 + 15.0*yh - 6.0*y2h + y3h).NormL2() / h2;
				
				*error = f6_approx * h2 * h2 / 90.0 + 
				         Constants::Eps * (ym2h.NormL2() + 16*ymh.NormL2() + 30*y0.NormL2() + 
				                           16*yh.NormL2() + y2h.NormL2()) / (12.0 * h2);
			}
			
			return result;
		}
		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer4(f, t, ScaleStep(NDer4_h, t), error);
		}

		/********************************************************************************************************************/
		/********                                       THIRD DERIVATIVES                                            ********/
		/********  NOTE: Direct finite difference formulas (not derivatives of derivatives!)                        ********/
		/********        NThirdDer2: 4 function evals (O(h²) accuracy)                                              ********/
		/********        NThirdDer4: 6 function evals (O(h⁴) accuracy)                                              ********/
		/********************************************************************************************************************/
		
		// f'''(t) ≈ [-f(t-2h) + 2f(t-h) - 2f(t+h) + f(t+2h)] / (2h³)
		// Second-order accurate (O(h²)), 4 function evaluations
		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y2h = f(t + 2*h);
			VectorN<Real, N> ym2h = f(t - 2*h);
			
			Real h3 = h * h * h;
			VectorN<Real, N> result = (-ym2h + 2.0*ymh - 2.0*yh + y2h) / (2.0 * h3);
			
			if (error)
			{
				// Error estimate using 5th derivative approximation
				VectorN<Real, N> y3h = f(t + 3*h);
				VectorN<Real, N> ym3h = f(t - 3*h);
				Real f5_approx = (ym3h - 3.0*ym2h + 5.0*ymh - 5.0*yh + 3.0*y2h - y3h).NormL2() / h3;
				
				*error = f5_approx * h * h / 4.0 + 
				         Constants::Eps * (ym2h.NormL2() + 2*ymh.NormL2() + 2*yh.NormL2() + y2h.NormL2()) / (2.0 * h3);
			}
			
			return result;
		}
		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			// Use larger step size for third derivatives (h³ in denominator needs bigger h)
			return NThirdDer2(f, t, ScaleStep(NDer4_h, t), error);
		}

		// f'''(t) ≈ [f(t-3h) - 8f(t-2h) + 13f(t-h) - 13f(t+h) + 8f(t+2h) - f(t+3h)] / (8h³)
		// Fourth-order accurate (O(h⁴)), 6 function evaluations
		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y2h = f(t + 2*h);
			VectorN<Real, N> ym2h = f(t - 2*h);
			VectorN<Real, N> y3h = f(t + 3*h);
			VectorN<Real, N> ym3h = f(t - 3*h);
			
			Real h3 = h * h * h;
			VectorN<Real, N> result = (ym3h - 8.0*ym2h + 13.0*ymh - 13.0*yh + 8.0*y2h - y3h) / (8.0 * h3);
			
			if (error)
			{
				// Error estimate using 7th derivative approximation
				VectorN<Real, N> y4h = f(t + 4*h);
				VectorN<Real, N> ym4h = f(t - 4*h);
				Real f7_approx = (ym4h - 4.0*ym3h + 9.0*ym2h - 13.0*ymh + 13.0*yh - 9.0*y2h + 4.0*y3h - y4h).NormL2() / h3;
				
				*error = f7_approx * h * h * h * h / 120.0 + 
				         Constants::Eps * (ym3h.NormL2() + 8*ym2h.NormL2() + 13*ymh.NormL2() + 
				                           13*yh.NormL2() + 8*y2h.NormL2() + y3h.NormL2()) / (8.0 * h3);
			}
			
			return result;
		}
		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer4(f, t, ScaleStep(NDer4_h, t), error);
		}
		
		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		template<int N>
		static inline VectorN<Real, N>(*DeriveCurve)(const IParametricCurve<N>& f, 
																								 Real x, Real* error) = Derivation::NDer4;
		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveSec)(const IParametricCurve<N>& f, 
																										Real x, Real* error) = Derivation::NSecDer4;
		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveThird)(const IParametricCurve<N>& f, 
																											Real x, Real* error) = Derivation::NThirdDer4;

	}
}

#endif // MML_DERIVATION_PARAMETRIC_CURVE_H