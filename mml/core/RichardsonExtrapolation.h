///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RichardsonExtrapolation.h                                           ///
///  Description: Richardson extrapolation algorithms for improved accuracy           ///
///               Used by Romberg integration, dfridr derivatives, and more           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_RICHARDSON_EXTRAPOLATION_H
#define MML_RICHARDSON_EXTRAPOLATION_H

#include <vector>
#include <cmath>
#include <limits>

#include "MMLBase.h"

namespace MML
{
	/// @brief Result of Richardson extrapolation
	struct RichardsonResult
	{
		Real value;           ///< Extrapolated value
		Real error_estimate;  ///< Estimated error from tableau differences
		int  iterations;      ///< Number of refinement steps used
		bool converged;       ///< True if converged within tolerance

		RichardsonResult(Real v = 0.0, Real e = 0.0, int it = 0, bool conv = false)
			: value(v), error_estimate(e), iterations(it), converged(conv) {}
	};

	namespace Richardson
	{
		/*************************************************************************/
		/*****         NEVILLE'S ALGORITHM (Polynomial Extrapolation)        *****/
		/*****         Used by Romberg integration - extrapolates to h=0     *****/
		/*************************************************************************/

		/// @brief Extrapolate sequence of values to limit using Neville's algorithm.
		/// @details Given values y[i] computed at "positions" x[i] (typically h² values),
		///          extrapolates to x=0 using polynomial interpolation.
		/// @param y Vector of computed values (at least 2 elements)
		/// @param x Vector of positions/step sizes (same size as y), must approach 0
		/// @return RichardsonResult with extrapolated value and error estimate
		static RichardsonResult NevilleExtrapolate(const std::vector<Real>& y, 
		                                            const std::vector<Real>& x)
		{
			int n = static_cast<int>(y.size());
			if (n < 2 || x.size() != y.size())
				return RichardsonResult(y.empty() ? 0.0 : y[0], 0.0, 0, false);

			std::vector<Real> c(n), d(n);
			
			// Find index of x closest to 0
			int ns = 0;
			Real dif = std::abs(x[0]);
			for (int i = 0; i < n; i++)
			{
				Real dift = std::abs(x[i]);
				if (dift < dif)
				{
					ns = i;
					dif = dift;
				}
				c[i] = y[i];
				d[i] = y[i];
			}

			// Initial best guess
			Real result = y[ns];
			Real error = 0.0;
			ns--;

			// Build Neville tableau
			for (int m = 1; m < n; m++)
			{
				for (int i = 0; i < n - m; i++)
				{
					Real ho = x[i];
					Real hp = x[i + m];
					Real w = c[i + 1] - d[i];
					Real den = ho - hp;

					if (std::abs(den) < Precision::DivisionSafetyThreshold)
						return RichardsonResult(result, std::abs(error), m, false);

					den = w / den;
					d[i] = hp * den;
					c[i] = ho * den;
				}

				// Choose correction from c or d based on position
				if (2 * (ns + 1) < (n - m))
					error = c[ns + 1];
				else
					error = d[ns--];

				result += error;
			}

			return RichardsonResult(result, std::abs(error), n, true);
		}

		/*************************************************************************/
		/*****         RICHARDSON TABLEAU (Ridders' Method Style)            *****/
		/*****         Used by dfridr - geometric step reduction              *****/
		/*************************************************************************/

		/// @brief Richardson extrapolation with geometric step reduction (dfridr-style).
		/// @details Builds a tableau of estimates computed at h, h/con, h/con², ...
		///          and extrapolates to h→0. Automatically finds optimal estimate.
		/// @tparam Evaluator Callable that computes estimate at given step size: Real(Real h)
		/// @param eval Function object: eval(h) returns approximation computed with step h
		/// @param h0 Initial step size
		/// @param con Step reduction factor (typically 1.4)
		/// @param max_iter Maximum tableau size (typically 10)
		/// @param safe Safety factor for convergence check (typically 2.0)
		/// @return RichardsonResult with best estimate and error
		template<typename Evaluator>
		static RichardsonResult Extrapolate(Evaluator eval, Real h0, 
		                                    Real con = 1.4, int max_iter = 10, Real safe = 2.0)
		{
			const Real con2 = con * con;
			const Real big = std::numeric_limits<Real>::max();

			if (std::abs(h0) < Precision::DivisionSafetyThreshold)
				return RichardsonResult(0.0, big, 0, false);

			// Allocate tableau (triangular, but use full matrix for clarity)
			std::vector<std::vector<Real>> a(max_iter, std::vector<Real>(max_iter));

			Real hh = h0;
			a[0][0] = eval(hh);
			Real err = big;
			Real ans = a[0][0];

			int final_iter = 1;
			for (int i = 1; i < max_iter; i++)
			{
				hh /= con;
				a[0][i] = eval(hh);

				// Richardson extrapolation: eliminate error terms
				Real fac = con2;
				for (int j = 1; j <= i; j++)
				{
					// A_j(h) = (con^(2j) * A_{j-1}(h/con) - A_{j-1}(h)) / (con^(2j) - 1)
					a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
					fac *= con2;

					// Error estimate: max difference from previous column
					Real errt = std::max(std::abs(a[j][i] - a[j - 1][i]), 
					                     std::abs(a[j][i] - a[j - 1][i - 1]));

					// Track best answer (lowest error)
					if (errt <= err)
					{
						err = errt;
						ans = a[j][i];
					}
				}

				final_iter = i + 1;

				// Convergence check: error growing means we've passed optimal h
				if (std::abs(a[i][i] - a[i - 1][i - 1]) >= safe * err)
					break;
			}

			return RichardsonResult(ans, err, final_iter, err < big);
		}

		/*************************************************************************/
		/*****         SIMPLE RICHARDSON (Known Error Order)                 *****/
		/*****         When error order p is known: A_true ≈ A(h) + c*h^p    *****/
		/*************************************************************************/

		/// @brief Simple Richardson extrapolation for known error order.
		/// @details Given A(h) and A(h/r) where error is O(h^p), computes improved estimate.
		///          Formula: A_improved = (r^p * A(h/r) - A(h)) / (r^p - 1)
		/// @param a_h Value computed at step h
		/// @param a_h_over_r Value computed at step h/r
		/// @param r Step ratio (typically 2)
		/// @param p Error order (e.g., 2 for trapezoidal rule)
		/// @return Improved estimate
		static Real SimpleExtrapolate(Real a_h, Real a_h_over_r, Real r, int p)
		{
			Real rp = std::pow(r, p);
			return (rp * a_h_over_r - a_h) / (rp - 1.0);
		}

		/// @brief Estimate error from simple Richardson extrapolation.
		/// @param a_h Value at step h
		/// @param a_h_over_r Value at step h/r
		/// @param r Step ratio
		/// @param p Error order
		/// @return Estimated error in a_h_over_r
		static Real SimpleErrorEstimate(Real a_h, Real a_h_over_r, Real r, int p)
		{
			Real rp = std::pow(r, p);
			return std::abs(a_h_over_r - a_h) / (rp - 1.0);
		}

	} // namespace Richardson

} // namespace MML

#endif // MML_RICHARDSON_EXTRAPOLATION_H
