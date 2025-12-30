///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationScalarFunction.h                                          ///
///  Description: Partial derivatives of scalar functions f:R^n->R                    ///
///               Gradient, Hessian, directional derivatives                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_SCALAR_FUNCTION_H
#define MML_DERIVATION_SCALAR_FUNCTION_H

#include "MMLBase.h"

#include "DerivationBase.h"

#include "base/VectorN.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x = point;
			Real y0 = f(x);

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = orig_x - h;
				Real ym = f(x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer1Partial(f, deriv_index, point, NDer1_h, error);
		}
		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, point, NDer1_h, error);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			auto    x = point;
			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f(x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f(x);

				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, deriv_index, point, NDer2_h, error);
		}
		
		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, point, NDer2_h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f(x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f(x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 
																				8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, point, NDer4_h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x);

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x[deriv_index] = orig_x + 4 * h;
				Real y4h = f(x);

				x[deriv_index] = orig_x - 4 * h;
				Real ym4h = f(x);

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, deriv_index, point, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, point, NDer6_h, error);
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x);

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x);

			x[deriv_index] = orig_x + 4 * h;
			Real y4h = f(x);

			x[deriv_index] = orig_x - 4 * h;
			Real ym4h = f(x);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x[deriv_index] = orig_x + 5 * h;
				Real y5h = f(x);

				x[deriv_index] = orig_x - 5 * h;
				Real ym5h = f(x);

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, deriv_index, point, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, point, NDer8_h, error);
		}
		
		/********************************************************************************************************************/
		/********                                      SECOND PARTIAL DERIVATIVES                                    ********/
		/********  NOTE: Direct finite difference formulas optimized for minimal function evaluations                ********/
		/********        For mixed partials (der_ind1 != der_ind2): ∂²f/∂x∂y                                        ********/
		/********        For pure partials (der_ind1 == der_ind2): ∂²f/∂x²                                          ********/
		/********                                                                                                    ********/
		/********        NSecDer2Partial: O(h²) accuracy, 4-9 function evals depending on pure/mixed                ********/
		/********        NSecDer4Partial: O(h⁴) accuracy, 5-13 function evals depending on pure/mixed               ********/
		/********************************************************************************************************************/
		
		// Second-order accurate (O(h²)) second partial derivative
		// For pure second partial (∂²f/∂xᵢ²): [f(x-h) - 2f(x) + f(x+h)] / h² - 3 function evaluations
		// For mixed partial (∂²f/∂xᵢ∂xⱼ): [f(x+h_i+h_j) - f(x+h_i-h_j) - f(x-h_i+h_j) + f(x-h_i-h_j)] / (4h²) - 4 evaluations
		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real h2 = h * h;
			
			if (der_ind1 == der_ind2)
			{
				// Pure second partial ∂²f/∂xᵢ² - use standard 3-point formula
				auto x_plus = point;
				auto x_minus = point;
				
				x_plus[der_ind1] = point[der_ind1] + h;
				x_minus[der_ind1] = point[der_ind1] - h;
				
				Real f0 = f(point);
				Real fh = f(x_plus);
				Real fmh = f(x_minus);
				
				Real result = (fmh - 2.0*f0 + fh) / h2;
				
				if (error)
				{
					// Error estimate using 4th derivative
					auto x_2plus = point;
					auto x_2minus = point;
					x_2plus[der_ind1] = point[der_ind1] + 2*h;
					x_2minus[der_ind1] = point[der_ind1] - 2*h;
					
					Real f4_approx = std::abs(f(x_2minus) - 4.0*fmh + 6.0*f0 - 4.0*fh + f(x_2plus)) / h2;
					*error = f4_approx * h2 / 12.0 + 
					         Constants::Eps * (std::abs(fmh) + 2*std::abs(f0) + std::abs(fh)) / h2;
				}
				
				return result;
			}
			else
			{
				// Mixed partial ∂²f/∂xᵢ∂xⱼ - use 4-point cross-difference formula
				auto x_pp = point;  // +h in both directions
				auto x_pm = point;  // +h in i, -h in j
				auto x_mp = point;  // -h in i, +h in j
				auto x_mm = point;  // -h in both directions
				
				x_pp[der_ind1] = point[der_ind1] + h;
				x_pp[der_ind2] = point[der_ind2] + h;
				
				x_pm[der_ind1] = point[der_ind1] + h;
				x_pm[der_ind2] = point[der_ind2] - h;
				
				x_mp[der_ind1] = point[der_ind1] - h;
				x_mp[der_ind2] = point[der_ind2] + h;
				
				x_mm[der_ind1] = point[der_ind1] - h;
				x_mm[der_ind2] = point[der_ind2] - h;
				
				Real fpp = f(x_pp);
				Real fpm = f(x_pm);
				Real fmp = f(x_mp);
				Real fmm = f(x_mm);
				
				Real result = (fpp - fpm - fmp + fmm) / (4.0 * h2);
				
				if (error)
				{
					// Approximate error using higher-order differences
					*error = Constants::Eps * (std::abs(fpp) + std::abs(fpm) + 
					                           std::abs(fmp) + std::abs(fmm)) / (4.0 * h2);
				}
				
				return result;
			}
		}
		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer2Partial(f, der_ind1, der_ind2, point, NDer2_h, error);
		}

		// Fourth-order accurate (O(h⁴)) second partial derivative
		// For pure second partial (∂²f/∂xᵢ²): [-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)] / (12h²) - 5 evaluations
		// For mixed partial (∂²f/∂xᵢ∂xⱼ): Higher-order cross-difference stencil - 13 evaluations
		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real h2 = h * h;
			
			if (der_ind1 == der_ind2)
			{
				// Pure second partial ∂²f/∂xᵢ² - use 5-point formula
				auto x = point;
				Real x_orig = point[der_ind1];
				
				Real f0 = f(point);
				
				x[der_ind1] = x_orig + h;
				Real fh = f(x);
				
				x[der_ind1] = x_orig - h;
				Real fmh = f(x);
				
				x[der_ind1] = x_orig + 2*h;
				Real f2h = f(x);
				
				x[der_ind1] = x_orig - 2*h;
				Real fm2h = f(x);
				
				Real result = (-fm2h + 16.0*fmh - 30.0*f0 + 16.0*fh - f2h) / (12.0 * h2);
				
				if (error)
				{
					// Error estimate using 6th derivative
					x[der_ind1] = x_orig + 3*h;
					Real f3h = f(x);
					x[der_ind1] = x_orig - 3*h;
					Real fm3h = f(x);
					
					Real f6_approx = std::abs(fm3h - 6.0*fm2h + 15.0*fmh - 20.0*f0 + 15.0*fh - 6.0*f2h + f3h) / h2;
					*error = f6_approx * h2 * h2 / 90.0 + 
					         Constants::Eps * (std::abs(fm2h) + 16*std::abs(fmh) + 30*std::abs(f0) + 
					                           16*std::abs(fh) + std::abs(f2h)) / (12.0 * h2);
				}
				
				return result;
			}
			else
			{
				// Mixed partial ∂²f/∂xᵢ∂xⱼ using Richardson extrapolation for O(h⁴) accuracy
				// 
				// Standard 4-point cross-difference formula (O(h²)):
				//   D_h = [f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)] / (4h²)
				//
				// Richardson extrapolation with step h and 2h:
				//   D_h  = exact + c·h² + O(h⁴)
				//   D_2h = exact + 4c·h² + O(h⁴)
				//   (4·D_h - D_2h) / 3 = exact + O(h⁴)
				//
				// Simplifying: (16·d1 - d2) / (48·h²) where:
				//   d1 = f(±h,±h) cross-difference
				//   d2 = f(±2h,±2h) cross-difference
				
				Real x_i = point[der_ind1];
				Real x_j = point[der_ind2];
				
				auto eval = [&](Real di, Real dj) {
					auto x = point;
					x[der_ind1] = x_i + di * h;
					x[der_ind2] = x_j + dj * h;
					return f(x);
				};
				
				// Cross-difference at step h
				Real f_p1_p1 = eval( 1,  1);
				Real f_p1_m1 = eval( 1, -1);
				Real f_m1_p1 = eval(-1,  1);
				Real f_m1_m1 = eval(-1, -1);
				Real d1 = (f_p1_p1 - f_p1_m1 - f_m1_p1 + f_m1_m1);
				
				// Cross-difference at step 2h  
				Real f_p2_p2 = eval( 2,  2);
				Real f_p2_m2 = eval( 2, -2);
				Real f_m2_p2 = eval(-2,  2);
				Real f_m2_m2 = eval(-2, -2);
				Real d2 = (f_p2_p2 - f_p2_m2 - f_m2_p2 + f_m2_m2);
				
				// Richardson extrapolation: (4·D_h - D_2h) / 3 = (16·d1 - d2) / (48·h²)
				Real result = (16.0 * d1 - d2) / (48.0 * h2);
				
				if (error)
				{
					*error = Constants::Eps * (std::abs(f_p1_p1) + std::abs(f_p1_m1) + 
					                           std::abs(f_m1_p1) + std::abs(f_m1_m1)) / (4.0 * h2);
				}
				
				return result;
			}
		}
		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer4Partial(f, der_ind1, der_ind2, point, NDer4_h, error);
		}

		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		template<int N>
		static inline Real(*DerivePartial)(const IScalarFunction<N>& f, int deriv_index, 
																			 const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;
		template<int N>
		static inline Real(*DeriveSecPartial)(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
																					const VectorN<Real, N>& point, Real* error) = Derivation::NSecDer4Partial;
		template<int N>
		static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																											VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;
	}
}

#endif // MML_DERIVATION_SCALAR_FUNCTION_H
