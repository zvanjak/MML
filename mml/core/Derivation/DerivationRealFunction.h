///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationRealFunction.h                                            ///
///  Description: Numerical derivatives of real-valued functions f:R->R               ///
///               Forward, backward, central differences, higher order                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_REAL_FUNCTION_H
#define MML_DERIVATION_REAL_FUNCTION_H

#include "MMLBase.h"

#include "DerivationBase.h"

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		static Real NDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real y0 = f(x);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = f(x - h);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
			}
			return diff / h;
		}
		static Real NDer1(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^1/2
			// Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
			// Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
			// This approximation will get better as we move to higher orders of accuracy.
			return NDer1(f, x, NDer1_h, error);
		}
		static Real NDer1(const IRealFunction& f, Real x)
		{
			return NDer1(f, x, NDer1_h, nullptr);
		}
		
		static Real NDer1Left(const IRealFunction& f, Real x, Real* error = nullptr) 
		{ return NDer1(f, x - 2 * NDer1_h, NDer1_h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real* error = nullptr) 
		{ return NDer1(f, x + 2 * NDer1_h, NDer1_h, error); }
		static Real NDer1Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) 
		{ return NDer1(f, x - 2 * h, h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) 
		{ return NDer1(f, x + 2 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		static Real NDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			
			// Check for non-finite function values
			if (std::isnan(yh) || std::isinf(yh) || std::isnan(ymh) || std::isinf(ymh))
				throw NumericalMethodError("Non-finite function values in numerical derivative calculation");
			
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = f(x + 2 * h);
				Real ym2h = f(x - 2 * h);
				*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + 
								 std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}
			
			Real result = diff / (2 * h);
			
			// Check for non-finite result
			if (std::isnan(result) || std::isinf(result))
				throw NumericalMethodError("Non-finite result in numerical derivative calculation");
			
			return result;
		}
		static Real NDer2(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^2/3
			// See the previous discussion to understand determination of h and the error bound.
			// Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]

			return NDer2(f, x, NDer2_h, error);
		}
		static Real NDer2(const IRealFunction& f, Real x)
		{
			return NDer2(f, x, NDer2_h, nullptr);
		}
		
		static Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x - 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x + 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x - 3 * h, h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x + 3 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		static Real NDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y2h = f(x + 2 * h);
			Real ym2h = f(x - 2 * h);

			// Check for non-finite function values
			if (!std::isfinite(yh) || !std::isfinite(ymh) || !std::isfinite(y2h) || !std::isfinite(ym2h))
				throw NumericalMethodError("Non-finite function values in numerical derivative calculation");

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				// Mathematica code to extract the remainder:
				// Series[(f[x-2*h]+ 8*f[x+h] - 8*f[x-h] - f[x+2*h])/(12*h), {h, 0, 7}]
				Real y3h = f(x + 3 * h);
				Real ym3h = f(x - 3 * h);

				// Error from fifth derivative:
				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				// Error from function evaluation:
				*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) + 
																				8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			
			Real result = (y2 + 8 * y1) / (12 * h);
			
			// Check for non-finite result
			if (!std::isfinite(result))
				throw NumericalMethodError("Non-finite result in numerical derivative calculation");
			
			return result;
		}
		static Real NDer4(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^4/5
			return NDer4(f, x, NDer4_h, error);
		}
		static Real NDer4(const IRealFunction& f, Real x)
		{
			return NDer4(f, x, NDer4_h, nullptr);
		}

		static Real NDer4Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x - 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x + 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x - 4 * h, h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x + 4 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		static Real NDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			const Real eps = (std::numeric_limits<Real>::epsilon)();

			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y1 = yh - ymh;
			Real y2 = f(x - 2 * h) - f(x + 2 * h);
			Real y3 = f(x + 3 * h) - f(x - 3 * h);

			if (error)
			{
				// Mathematica code to generate fd scheme for 7th derivative:
				// Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
				// Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
				// Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
				Real y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		static Real NDer6(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, error);
		}
		static Real NDer6(const IRealFunction& f, Real x)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, nullptr);
		}

		static Real NDer6Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x - 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x + 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x - 5 * h, h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x + 5 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		static Real NDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y1 = yh - ymh;
			Real y2 = f(x - 2 * h) - f(x + 2 * h);
			Real y3 = f(x + 3 * h) - f(x - 3 * h);
			Real y4 = f(x - 4 * h) - f(x + 4 * h);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				// Mathematica code to generate fd scheme for 7th derivative:
				// Sum[(-1)^i*Binomial[9, i]*(f[x+(4-i)*h] + f[x+(5-i)*h])/2, {i, 0, 9}]
				// Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
				// Series[(f[x+5*h]-f[x- 5*h])/2 + 4*(f[x-4*h] - f[x+4*h]) + 27*(f[x+3*h] - f[x-3*h])/2 + 24*(f[x-2*h]  - f[x+2*h]) + 21*(f[x+h] - f[x-h]), {h, 0, 15}]
				Real f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}
		static Real NDer8(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^8/9.
			// In Real precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
			// Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
			// Mathematica code to get the error:
			// Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
			// If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.

			return NDer8(f, x, NDer8_h, error);
		}
		static Real NDer8(const IRealFunction& f, Real x)
		{
			return NDer8(f, x, NDer8_h, nullptr);
		}

		static Real NDer8Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x - 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x + 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x - 6 * h, h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x + 6 * h, h, error); }

		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********  NOTE: Direct finite difference formulas (not derivatives of derivatives!)                        ********/
		/********        NSecDer2: 3 function evals (O(h²) accuracy)                                                ********/
		/********        NSecDer4: 5 function evals (O(h⁴) accuracy)                                                ********/
		/********************************************************************************************************************/
		
		// f''(x) ≈ [f(x-h) - 2f(x) + f(x+h)] / h²
		// Second-order accurate (O(h²)), 3 function evaluations
		static Real NSecDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real y0 = f(x);
			Real yh = f(x + h);
			Real ymh = f(x - h);
			
			Real h2 = h * h;
			Real result = (ymh - 2.0 * y0 + yh) / h2;
			
			if (error)
			{
				// Error estimate using 4th derivative approximation
				Real y2h = f(x + 2*h);
				Real ym2h = f(x - 2*h);
				Real f4_approx = std::abs(ym2h - 4*ymh + 6*y0 - 4*yh + y2h) / h2;
				
				*error = f4_approx * h2 / 12.0 + 
				         Constants::Eps * (std::abs(ymh) + 2*std::abs(y0) + std::abs(yh)) / h2;
			}
			
			return result;
		}
		static Real NSecDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer2(f, x, NDer2_h, error);
		}

		// f''(x) ≈ [-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)] / (12h²)
		// Fourth-order accurate (O(h⁴)), 5 function evaluations
		static Real NSecDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real y0 = f(x);
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y2h = f(x + 2*h);
			Real ym2h = f(x - 2*h);
			
			Real h2 = h * h;
			Real result = (-ym2h + 16.0*ymh - 30.0*y0 + 16.0*yh - y2h) / (12.0 * h2);
			
			if (error)
			{
				// Error estimate using 6th derivative approximation
				Real y3h = f(x + 3*h);
				Real ym3h = f(x - 3*h);
				Real f6_approx = std::abs(ym3h - 6*ym2h + 15*ymh - 20*y0 + 15*yh - 6*y2h + y3h) / h2;
				
				*error = f6_approx * h2 * h2 / 90.0 + 
				         Constants::Eps * (std::abs(ym2h) + 16*std::abs(ymh) + 30*std::abs(y0) + 
				                           16*std::abs(yh) + std::abs(y2h)) / (12.0 * h2);
			}
			
			return result;
		}
		static Real NSecDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer4(f, x, NDer4_h, error);
		}

		/********************************************************************************************************************/
		/********                                       THIRD DERIVATIVES                                            ********/
		/********  NOTE: Direct finite difference formulas (not derivatives of derivatives!)                        ********/
		/********        NThirdDer2: 4 function evals (O(h²) accuracy)                                              ********/
		/********        NThirdDer4: 6 function evals (O(h⁴) accuracy)                                              ********/
		/********************************************************************************************************************/
		
		// f'''(x) ≈ [-f(x-2h) + 2f(x-h) - 2f(x+h) + f(x+2h)] / (2h³)
		// Second-order accurate (O(h²)), 4 function evaluations
		static Real NThirdDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y2h = f(x + 2*h);
			Real ym2h = f(x - 2*h);
			
			Real h3 = h * h * h;
			Real result = (-ym2h + 2.0*ymh - 2.0*yh + y2h) / (2.0 * h3);
			
			if (error)
			{
				// Error estimate using 5th derivative approximation
				Real y3h = f(x + 3*h);
				Real ym3h = f(x - 3*h);
				Real f5_approx = std::abs(ym3h - 3*ym2h + 5*ymh - 5*yh + 3*y2h - y3h) / h3;
				
				*error = f5_approx * h * h / 4.0 + 
				         Constants::Eps * (std::abs(ym2h) + 2*std::abs(ymh) + 2*std::abs(yh) + std::abs(y2h)) / (2.0 * h3);
			}
			
			return result;
		}
		static Real NThirdDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			// Use larger step size for third derivatives (h³ in denominator needs bigger h)
			return NThirdDer2(f, x, NDer4_h, error);
		}

		// f'''(x) ≈ [f(x-3h) - 8f(x-2h) + 13f(x-h) - 13f(x+h) + 8f(x+2h) - f(x+3h)] / (8h³)
		// Fourth-order accurate (O(h⁴)), 6 function evaluations
		static Real NThirdDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y2h = f(x + 2*h);
			Real ym2h = f(x - 2*h);
			Real y3h = f(x + 3*h);
			Real ym3h = f(x - 3*h);
			
			Real h3 = h * h * h;
			Real result = (ym3h - 8.0*ym2h + 13.0*ymh - 13.0*yh + 8.0*y2h - y3h) / (8.0 * h3);
			
			if (error)
			{
				// Error estimate using 7th derivative approximation
				Real y4h = f(x + 4*h);
				Real ym4h = f(x - 4*h);
				Real f7_approx = std::abs(ym4h - 4*ym3h + 9*ym2h - 13*ymh + 13*yh - 9*y2h + 4*y3h - y4h) / h3;
				
				*error = f7_approx * h * h * h * h / 120.0 + 
				         Constants::Eps * (std::abs(ym3h) + 8*std::abs(ym2h) + 13*std::abs(ymh) + 
				                           13*std::abs(yh) + 8*std::abs(y2h) + std::abs(y3h)) / (8.0 * h3);
			}
			
			return result;
		}
		static Real NThirdDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer4(f, x, NDer4_h, error);
		}

		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		static inline Real(*Derive)(const IRealFunction& f, Real x) = Derivation::NDer4;
		static inline Real(*DeriveErr)(const IRealFunction& f, 
																	 Real x, Real* error) = Derivation::NDer4;
		static inline Real(*DeriveSec)(const IRealFunction& f, 
																	 Real x, Real* error) = Derivation::NSecDer4;
		static inline Real(*DeriveThird)(const IRealFunction& f, 
																		 Real x, Real* error) = Derivation::NThirdDer2;
	}
}

#endif // MML_DERIVATION_REAL_FUNCTION_H