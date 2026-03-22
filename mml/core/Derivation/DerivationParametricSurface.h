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

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/Vector/VectorN.h"
#include "base/Matrix/MatrixNM.h"

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
			VectorN<Real, N> yh = f(u + h, w);
			VectorN<Real, N> y0 = f(u, w);
			VectorN<Real, N> diff = yh - y0;

			if (error)
			{
				VectorN<Real, N> ym = f(u - h, w);
				VectorN<Real, N> ypph_vec = yh - 2 * y0 + ym;

				Real ypph = ypph_vec.NormL2() / h;

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Eps / h;
			}
			return diff / h;
		}
		template <int N>
		static VectorN<Real, N> NDer1_u(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer1_u(f, u, w, ScaleStep(NDer1_h, u), error);
		}
		template <int N>
		static VectorN<Real, N> NDer1_w(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(u, w + h);
			VectorN<Real, N> y0 = f(u, w);
			VectorN<Real, N> diff = yh - y0;

			if (error)
			{
				VectorN<Real, N> ym = f(u, w - h);
				VectorN<Real, N> ypph_vec = yh - 2 * y0 + ym;

				Real ypph = ypph_vec.NormL2() / h;

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Eps / h;
			}
			return diff / h;
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
			VectorN<Real, N> yh = f(u + h, w);
			VectorN<Real, N> ymh = f(u - h, w);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = f(u + 2 * h, w);
				VectorN<Real, N> ymth = f(u - 2 * h, w);

				*error = Constants::Eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}
		template <int N>
		static VectorN<Real, N> NDer2_u(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer2_u(f, u, w, ScaleStep(NDer2_h, u), error);
		}
		template <int N>
		static VectorN<Real, N> NDer2_w(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(u, w + h);
			VectorN<Real, N> ymh = f(u, w - h);
			VectorN<Real, N> diff = yh - ymh;
			if (error)
			{
				VectorN<Real, N> yth = f(u, w + 2 * h);
				VectorN<Real, N> ymth = f(u, w - 2 * h);
				*error = Constants::Eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
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