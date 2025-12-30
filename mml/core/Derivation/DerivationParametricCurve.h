///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationParametricCurve.h                                         ///
///  Description: Derivatives of parametric curves                                    ///
///               Tangent, normal, binormal, curvature, torsion calculations          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_PARAMETRIC_CURVE_H
#define MML_DERIVATION_PARAMETRIC_CURVE_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"

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
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> y0 = f(t);
			VectorN<Real, N> diff = yh - y0;

			if (error)
			{
				VectorN<Real, N> ym = f(t - h);
				VectorN<Real, N> ypph_vec = yh - 2 * y0 + ym;

				Real ypph = ypph_vec.NormL2() / h;

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Eps / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer1(f, t, NDer1_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = f(t + 2 * h);
				VectorN<Real, N> ymth = f(t - 2 * h);

				*error = Constants::Eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer2(f, t, NDer2_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y2h = f(t + 2 * h);
			VectorN<Real, N> ym2h = f(t - 2 * h);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = f(t + 3 * h);
				VectorN<Real, N> ym3h = f(t - 3 * h);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer4(f, t, NDer4_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = f(t - 2 * h) - f(t + 2 * h);
			VectorN<Real, N> y3 = f(t + 3 * h) - f(t - 3 * h);

			if (error)
			{
				VectorN<Real, N> y7 = (f(t + 4 * h) - f(t - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer6(f, t, NDer6_h, error);
		}	
		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = f(t - 2 * h) - f(t + 2 * h);
			VectorN<Real, N> y3 = f(t + 3 * h) - f(t - 3 * h);
			VectorN<Real, N> y4 = f(t - 4 * h) - f(t + 4 * h);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (f(t + 5 * h) - f(t - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Eps / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer8(f, t, NDer8_h, error);
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
			return NSecDer2(f, t, NDer2_h, error);
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
			return NSecDer4(f, t, NDer4_h, error);
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
			return NThirdDer2(f, t, NDer4_h, error);
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
			return NThirdDer4(f, t, NDer4_h, error);
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