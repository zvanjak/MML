///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationParametricSurface.h                                       ///
///  Description: Derivatives of parametric surfaces                                  ///
///               Tangent planes, normal vectors, curvature                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_PARAMETRIC_SURFACE_H
#define MML_DERIVATION_PARAMETRIC_SURFACE_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

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
			return NDer1_u(f, u, w, NDer1_h, error);
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
			return NDer1_w(f, u, w, NDer1_h, error);
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
			return NDer2_u(f, u, w, NDer2_h, error);
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
			return NDer2_w(f, u, w, NDer2_h, error);
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
			return NDer2_uu(f, u, w, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer2_uw(const IParametricSurfaceRect<N>& f, Real u, Real w, Real h, Real* error = nullptr)
		{
			// Correct formula for mixed partial derivative: ∂²f/∂u∂w
			// Use central difference: [f(u+h,w+h) - f(u+h,w-h) - f(u-h,w+h) + f(u-h,w-h)] / (4h²)
			VectorN<Real, N> fpp = f(u + h, w + h);  // f(u+h, w+h)
			VectorN<Real, N> fpm = f(u + h, w - h);  // f(u+h, w-h)
			VectorN<Real, N> fmp = f(u - h, w + h);  // f(u-h, w+h)
			VectorN<Real, N> fmm = f(u - h, w - h);  // f(u-h, w-h)
			
			VectorN<Real, N> diff = fpp - fpm - fmp + fmm;
			
			if (error)
			{
				// Error estimate for mixed partial
				Real norm_sum = fpp.NormL2() + fpm.NormL2() + fmp.NormL2() + fmm.NormL2();
				*error = norm_sum * Constants::Eps / (4 * h * h);
			}
			return diff / (4 * h * h);
		}
		template <int N>
		static VectorN<Real, N> NDer2_uw(const IParametricSurfaceRect<N>& f, Real u, Real w, Real* error = nullptr)
		{
			return NDer2_uw(f, u, w, NDer2_h, error);
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
			return NDer2_ww(f, u, w, NDer2_h, error);
		}
	}
}

#endif // MML_DERIVATION_PARAMETRIC_SURFACE_H