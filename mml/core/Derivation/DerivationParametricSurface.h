///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationParametricSurface.h                                       ///
///  Description: Derivatives of parametric surfaces                                  ///
///               Tangent planes, normal vectors, curvature                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_PARAMETRIC_SURFACE_H
#define MML_DERIVATION_PARAMETRIC_SURFACE_H

#include <mml/MMLBase.h>

#include "DerivationBase.h"
#include "FirstDerivativeStencil.h"

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>

namespace MML
{
	namespace Derivation
	{
		///////////////////////             FIRST DERIVATIONS              //////////////////////////
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer1_u(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::One>(
				[&](int offset) { return f(u + offset * h, w); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static VectorN<Real, N> NDer1_u(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer1_u(f, u, w, ScaleStep(NDer1_h, u), error);
		}
		template <int N>
		static VectorN<Real, N> NDer1_w(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::One>(
				[&](int offset) { return f(u, w + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static VectorN<Real, N> NDer1_w(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer1_w(f, u, w, ScaleStep(NDer1_h, w), error);
		}
	
		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer2_u(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>;
			auto result = Detail::EvaluateFirstDerivativeStencilWithOffsets<Detail::FirstDerivativeOrder::Two>(
				[&](int offset) { return f(u + offset * h, w); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr,
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
		static VectorN<Real, N> NDer2_u(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer2_u(f, u, w, ScaleStep(NDer2_h, u), error);
		}
		template <int N>
		static VectorN<Real, N> NDer2_w(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>;
			auto result = Detail::EvaluateFirstDerivativeStencilWithOffsets<Detail::FirstDerivativeOrder::Two>(
				[&](int offset) { return f(u, w + offset * h); }, [](const auto& value) { return value.NormL2(); }, h, error != nullptr,
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
		static VectorN<Real, N> NDer2_w(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer2_w(f, u, w, ScaleStep(NDer2_h, w), error);
		}

		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer2_uu(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(u + h, w);
			VectorN<Real, N> ym = f(u - h, w);
			VectorN<Real, N> y0 = f(u, w);
			VectorN<Real, N> diff = yh - 2 * y0 + ym;
			if (error)
			{
				Real ypph = diff.NormL2() / (h * h);
				*error = ypph / 2 + (yh.NormL2() + ym.NormL2() + y0.NormL2()) * Constants::Eps / (h * h);
			}
			return diff / (h * h);
		}
		template <int N>
		static VectorN<Real, N> NDer2_uu(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer2_uu(f, u, w, ScaleStep(NDer2_h, u), error);
		}

		template <int N>
		static VectorN<Real, N> NDer2_uw(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h_u, Real h_w, Real* error = nullptr)
		{
			// Mixed partial derivative: ∂²f/∂u∂w
			// Central difference: [f(u+hu,w+hw) - f(u+hu,w-hw) - f(u-hu,w+hw) + f(u-hu,w-hw)] / (4*hu*hw)
			VectorN<Real, N> fpp = f(u + h_u, w + h_w);
			VectorN<Real, N> fpm = f(u + h_u, w - h_w);
			VectorN<Real, N> fmp = f(u - h_u, w + h_w);
			VectorN<Real, N> fmm = f(u - h_u, w - h_w);
			
			VectorN<Real, N> diff = fpp - fpm - fmp + fmm;
			
			if (error)
			{
				Real norm_sum = fpp.NormL2() + fpm.NormL2() + fmp.NormL2() + fmm.NormL2();
				*error = norm_sum * Constants::Eps / (4 * h_u * h_w);
			}
			return diff / (4 * h_u * h_w);
		}
		template <int N>
		static VectorN<Real, N> NDer2_uw(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			return NDer2_uw<N>(f, u, w, h, h, error);
		}
		template <int N>
		static VectorN<Real, N> NDer2_uw(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			// Scale step by max parameter magnitude for stability at large u, w
			Real scale = std::max(std::abs(u), std::abs(w));
			return NDer2_uw(f, u, w, ScaleStep(NDer2_h, scale), error);
		}

		template <int N>
		static VectorN<Real, N> NDer2_ww(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(u, w + h);
			VectorN<Real, N> ym = f(u, w - h);
			VectorN<Real, N> y0 = f(u, w);
			VectorN<Real, N> diff = yh - 2 * y0 + ym;
			if (error)
			{
				Real ypph = diff.NormL2() / (h * h);
				*error = ypph / 2 + (yh.NormL2() + ym.NormL2() + y0.NormL2()) * Constants::Eps / (h * h);
			}
			return diff / (h * h);
		}
		template <int N>
		static VectorN<Real, N> NDer2_ww(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer2_ww(f, u, w, ScaleStep(NDer2_h, w), error);
		}
	}
}

#endif // MML_DERIVATION_PARAMETRIC_SURFACE_H