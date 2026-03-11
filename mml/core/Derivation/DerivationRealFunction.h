///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationRealFunction.h                                            ///
///  Description: Numerical derivatives of real-valued functions f:R->R               ///
///               Forward, backward, central differences, higher order                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
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
			return NDer1(f, x, ScaleStep(NDer1_h, x), error);
		}
		static Real NDer1(const IRealFunction& f, Real x)
		{
			return NDer1(f, x, ScaleStep(NDer1_h, x), nullptr);
		}
		
		static Real NDer1Left(const IRealFunction& f, Real x, Real* error = nullptr) 
		{ Real h = ScaleStep(NDer1_h, x); return NDer1(f, x - 2 * h, h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real* error = nullptr) 
		{ Real h = ScaleStep(NDer1_h, x); return NDer1(f, x + 2 * h, h, error); }
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
			return NDer2(f, x, ScaleStep(NDer2_h, x), error);
		}
		static Real NDer2(const IRealFunction& f, Real x)
		{
			return NDer2(f, x, ScaleStep(NDer2_h, x), nullptr);
		}
		
		static Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer2_h, x); return NDer2(f, x - 2 * h, h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer2_h, x); return NDer2(f, x + 2 * h, h, error); }
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
			return NDer4(f, x, ScaleStep(NDer4_h, x), error);
		}
		static Real NDer4(const IRealFunction& f, Real x)
		{
			return NDer4(f, x, ScaleStep(NDer4_h, x), nullptr);
		}

		static Real NDer4Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer4_h, x); return NDer4(f, x - 4 * h, h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer4_h, x); return NDer4(f, x + 4 * h, h, error); }
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
			return NDer6(f, x, ScaleStep(NDer6_h, x), error);
		}
		static Real NDer6(const IRealFunction& f, Real x)
		{
			return NDer6(f, x, ScaleStep(NDer6_h, x), nullptr);
		}

		static Real NDer6Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer6_h, x); return NDer6(f, x - 5 * h, h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer6_h, x); return NDer6(f, x + 5 * h, h, error); }
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
			return NDer8(f, x, ScaleStep(NDer8_h, x), error);
		}
		static Real NDer8(const IRealFunction& f, Real x)
		{
			return NDer8(f, x, ScaleStep(NDer8_h, x), nullptr);
		}

		static Real NDer8Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer8_h, x); return NDer8(f, x - 6 * h, h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer8_h, x); return NDer8(f, x + 6 * h, h, error); }
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
			return NSecDer2(f, x, ScaleStep(NDer2_h, x), error);
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
			return NSecDer4(f, x, ScaleStep(NDer4_h, x), error);
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
			return NThirdDer2(f, x, ScaleStep(NDer4_h, x), error);
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
			return NThirdDer4(f, x, ScaleStep(NDer4_h, x), error);
		}

		/********************************************************************************************************************/
		/********              RICHARDSON EXTRAPOLATION DERIVATIVES (dfridr-style adaptive)                          ********/
		/********  Uses tableau-based Richardson extrapolation to achieve high accuracy automatically               ********/
		/********  Based on Numerical Recipes "dfridr" algorithm with MML Richardson utility                        ********/
		/********************************************************************************************************************/
		
		/// @brief First derivative using Richardson extrapolation (dfridr algorithm).
		/// @details Builds a Richardson tableau using successively smaller step sizes,
		///          automatically finding the optimal step size. Uses central differences.
		///          Typically achieves ~10-12 digits of accuracy for smooth functions.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate derivative
		/// @param h Initial step size (will be refined automatically). Default: 0.1
		/// @param error Output parameter for error estimate (optional)
		/// @param max_iter Maximum tableau depth (default 10)
		/// @param con Step shrinking factor (default 1.4)
		/// @return Derivative estimate with error estimate if requested
		static Real NDerRichardson(const IRealFunction& f, Real x, Real h = 0.1, 
		                           Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			const Real con2 = con * con;
			const Real big = std::numeric_limits<Real>::max();
			const Real safe = 2.0;

			if (h == 0.0)
			{
				if (error) *error = big;
				return 0.0;
			}

			// Allocate tableau (triangular matrix stored as 2D)
			std::vector<std::vector<Real>> a(max_iter, std::vector<Real>(max_iter));

			Real hh = h;
			a[0][0] = (f(x + hh) - f(x - hh)) / (2.0 * hh);
			Real err = big;
			Real ans = a[0][0];

			for (int i = 1; i < max_iter; i++)
			{
				// Shrink step size
				hh /= con;
				a[0][i] = (f(x + hh) - f(x - hh)) / (2.0 * hh);

				// Richardson extrapolation: eliminate error terms
				Real fac = con2;
				for (int j = 1; j <= i; j++)
				{
					// A_j(h) = (con^(2j) * A_{j-1}(h/con) - A_{j-1}(h)) / (con^(2j) - 1)
					a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
					fac *= con2;

					// Error estimate: max difference from neighboring tableau entries
					Real errt = std::max(std::abs(a[j][i] - a[j - 1][i]),
					                     std::abs(a[j][i] - a[j - 1][i - 1]));

					// Track best answer (lowest error)
					if (errt <= err)
					{
						err = errt;
						ans = a[j][i];
					}
				}

				// Convergence check: error growing means we've passed optimal h
				if (std::abs(a[i][i] - a[i - 1][i - 1]) >= safe * err)
					break;
			}

			if (error) *error = err;
			return ans;
		}

		/// @brief First derivative using Richardson extrapolation with default parameters.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Derivative estimate
		static Real NDerRichardson(const IRealFunction& f, Real x, Real* error)
		{
			return NDerRichardson(f, x, 0.1, error, 10, 1.4);
		}

		/// @brief Second derivative using Richardson extrapolation.
		/// @details Uses central difference formula for second derivative with Richardson 
		///          extrapolation to achieve high accuracy automatically.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate second derivative
		/// @param h Initial step size (will be refined automatically). Default: 0.1
		/// @param error Output parameter for error estimate (optional)
		/// @param max_iter Maximum tableau depth (default 10)
		/// @param con Step shrinking factor (default 1.4)
		/// @return Second derivative estimate
		static Real NSecDerRichardson(const IRealFunction& f, Real x, Real h = 0.1,
		                              Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			const Real con2 = con * con;
			const Real big = std::numeric_limits<Real>::max();
			const Real safe = 2.0;

			if (h == 0.0)
			{
				if (error) *error = big;
				return 0.0;
			}

			std::vector<std::vector<Real>> a(max_iter, std::vector<Real>(max_iter));

			Real hh = h;
			Real fx = f(x);
			a[0][0] = (f(x + hh) - 2.0 * fx + f(x - hh)) / (hh * hh);
			Real err = big;
			Real ans = a[0][0];

			for (int i = 1; i < max_iter; i++)
			{
				hh /= con;
				a[0][i] = (f(x + hh) - 2.0 * fx + f(x - hh)) / (hh * hh);

				Real fac = con2;
				for (int j = 1; j <= i; j++)
				{
					a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
					fac *= con2;

					Real errt = std::max(std::abs(a[j][i] - a[j - 1][i]),
					                     std::abs(a[j][i] - a[j - 1][i - 1]));

					if (errt <= err)
					{
						err = errt;
						ans = a[j][i];
					}
				}

				if (std::abs(a[i][i] - a[i - 1][i - 1]) >= safe * err)
					break;
			}

			if (error) *error = err;
			return ans;
		}

		/// @brief Second derivative using Richardson extrapolation with default parameters.
		static Real NSecDerRichardson(const IRealFunction& f, Real x, Real* error)
		{
			return NSecDerRichardson(f, x, 0.1, error, 10, 1.4);
		}

		/// @brief Third derivative using Richardson extrapolation.
		/// @details Uses central difference formula for third derivative with Richardson
		///          extrapolation. Requires more iterations for convergence.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate third derivative
		/// @param h Initial step size (default 0.2, larger than for lower derivatives)
		/// @param error Output parameter for error estimate (optional)
		/// @param max_iter Maximum tableau depth (default 10)
		/// @param con Step shrinking factor (default 1.4)
		/// @return Third derivative estimate
		static Real NThirdDerRichardson(const IRealFunction& f, Real x, Real h = 0.2,
		                                Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			const Real con2 = con * con;
			const Real big = std::numeric_limits<Real>::max();
			const Real safe = 2.0;

			if (h == 0.0)
			{
				if (error) *error = big;
				return 0.0;
			}

			std::vector<std::vector<Real>> a(max_iter, std::vector<Real>(max_iter));

			Real hh = h;
			// f'''(x) ≈ [f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)] / (2h³)
			Real h3 = hh * hh * hh;
			a[0][0] = (f(x + 2*hh) - 2.0*f(x + hh) + 2.0*f(x - hh) - f(x - 2*hh)) / (2.0 * h3);
			Real err = big;
			Real ans = a[0][0];

			for (int i = 1; i < max_iter; i++)
			{
				hh /= con;
				h3 = hh * hh * hh;
				a[0][i] = (f(x + 2*hh) - 2.0*f(x + hh) + 2.0*f(x - hh) - f(x - 2*hh)) / (2.0 * h3);

				Real fac = con2;
				for (int j = 1; j <= i; j++)
				{
					a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
					fac *= con2;

					Real errt = std::max(std::abs(a[j][i] - a[j - 1][i]),
					                     std::abs(a[j][i] - a[j - 1][i - 1]));

					if (errt <= err)
					{
						err = errt;
						ans = a[j][i];
					}
				}

				if (std::abs(a[i][i] - a[i - 1][i - 1]) >= safe * err)
					break;
			}

			if (error) *error = err;
			return ans;
		}

		/// @brief Third derivative using Richardson extrapolation with default parameters.
		static Real NThirdDerRichardson(const IRealFunction& f, Real x, Real* error)
		{
			return NThirdDerRichardson(f, x, 0.2, error, 10, 1.4);
		}

		/********************************************************************************************************************/
		/********                          ADAPTIVE STEP SIZE DERIVATIVES (Aliases)                                  ********/
		/********  These are aliases for Richardson extrapolation methods, which automatically                       ********/
		/********  find the optimal step size through tableau-based refinement.                                     ********/
		/********************************************************************************************************************/

		/// @brief Adaptive first derivative (alias for NDerRichardson).
		/// @details Automatically finds optimal step size using Richardson extrapolation.
		///          This solves the fundamental problem of choosing h: too small causes
		///          cancellation error, too large causes truncation error.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Derivative estimate
		static Real NDerAdaptive(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NDerRichardson(f, x, 0.1, error, 10, 1.4);
		}

		/// @brief Adaptive second derivative (alias for NSecDerRichardson).
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate second derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Second derivative estimate
		static Real NSecDerAdaptive(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDerRichardson(f, x, 0.1, error, 10, 1.4);
		}

		/// @brief Adaptive third derivative (alias for NThirdDerRichardson).
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate third derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Third derivative estimate
		static Real NThirdDerAdaptive(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDerRichardson(f, x, 0.2, error, 10, 1.4);
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