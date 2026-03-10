///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Integration1D.h                                                     ///
///  Description: 1D numerical integration (Trapezoidal, Simpson, Romberg, Gauss)     ///
///               Adaptive quadrature with error control and convergence diagnostics  ///
///                                                                                   ///
///  Methods:     TRAP     - Extended trapezoidal rule (adaptive)                     ///
///               SIMPSON  - Simpson's rule with Richardson extrapolation             ///
///               ROMBERG  - Romberg integration (exponentially fast for smooth f)    ///
///               GAUSS10  - 10-point Gauss-Legendre (exact for poly deg ≤ 19)        ///
///                                                                                   ///
///  Usage:       IntegrationResult r = IntegrateSimpson(f, a, b, eps);               ///
///               if (r.converged) std::cout << r.value << " ± " << r.error_estimate; ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_1D_H
#define MML_INTEGRATION_1D_H

#include <vector>

#include "MMLBase.h"
#include "IntegrationBase.h"
#include "GaussKronrod.h"

#include "interfaces/IFunction.h"


namespace MML
{
	/// Integration method selection for 1D/2D/3D integration
	enum IntegrationMethod { TRAP, SIMPSON, ROMBERG, GAUSS10, GAUSS10KRONROD21 };

	/// @brief Abstract base for quadrature refinement
	class IQuadrature {
	protected:
		int _currStep = 0;

	public:
		// Returns the nth stage of refinement of the extended trapezoidal rule.
		virtual Real next() = 0;
	};

	/// @brief Trapezoidal rule integrator with progressive refinement
	/// @details Each call to next() doubles the number of evaluation points
	class TrapIntegrator : IQuadrature
	{
	public:
		Real _a, _b, _currSum;
		const IRealFunction& _func;

		TrapIntegrator(const IRealFunction& func, Real a, Real b) 
			:	_func(func), _a(a), _b(b), _currSum(0.0)	{	}

		Real next() {
			Real x, sum, del;
			int subDivNum, j;

			_currStep++;
			if (_currStep == 1) {
				return (_currSum = 0.5 * (_b - _a) * (_func(_a) + _func(_b)));
			}
			else {
				for (subDivNum = 1, j = 1; j < _currStep - 1; j++)
					subDivNum *= 2;

				del = (_b - _a) / subDivNum;
				x = _a + 0.5 * del;

				for (sum = 0.0, j = 0; j < subDivNum; j++, x += del)
					sum += _func(x);

				_currSum = 0.5 * (_currSum + (_b - _a) * sum / subDivNum);

				return _currSum;
			}
		}
	};

	/// @brief Extended trapezoidal rule with adaptive refinement
	/// @param func Function to integrate
	/// @param a,b Integration bounds
	/// @param eps Desired relative accuracy (default: 1e-6)
	/// @return IntegrationResult with value, error_estimate, iterations, converged flag
	static IntegrationResult IntegrateTrap(const IRealFunction& func, Real a, Real b,
																		 const Real eps = Defaults::TrapezoidIntegrationEPS)	{
		int		j;
		Real	currSum, oldSum = 0.0;

		TrapIntegrator t(func, a, b);

		for (j = 0; j < Defaults::TrapezoidIntegrationMaxSteps; j++)
		{
			currSum = t.next();

			if (j > 5) {
				Real error = std::abs(currSum - oldSum);
				if ( error < eps * std::abs(oldSum) || (currSum == 0.0 && oldSum == 0.0) )
				{
					return IntegrationResult(currSum, error, j, true);
				}
			}

			oldSum = currSum;
		}
		// Did not converge within max iterations
		Real final_error = std::abs(currSum - oldSum);
		return IntegrationResult(currSum, final_error, j, false);
	}

	/// @brief Simpson's rule with adaptive refinement
	/// @details Uses 4/3 Richardson extrapolation on trapezoidal estimates
	/// @param func Function to integrate
	/// @param a,b Integration bounds
	/// @param eps Desired relative accuracy (default: 1e-8)
	/// @return IntegrationResult with value, error_estimate, iterations, converged flag
	static IntegrationResult IntegrateSimpson(const IRealFunction& func, Real a, Real b,
																			 const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		int j;
		Real currSum, st, ost = 0.0, oldSum = 0.0;

		TrapIntegrator t(func, a, b);

		for (j = 0; j < Defaults::SimpsonIntegrationMaxSteps; j++)
		{
			st = t.next();

			currSum = (4.0 * st - ost) / 3.0;

			if (j > 5) {
				Real error = std::abs(currSum - oldSum);
				if (error < eps * std::abs(oldSum) || (currSum == 0.0 && oldSum == 0.0))
				{
					return IntegrationResult(currSum, error, j, true);
				}
			}

			oldSum = currSum;
			ost = st;
		}
		// Did not converge within max iterations
		Real final_error = std::abs(currSum - oldSum);
		return IntegrationResult(currSum, final_error, j, false);
	}


	/**
	 * @brief IntegrateRomberg - Romberg integration using Richardson extrapolation.
	 *
	 * Romberg integration is an adaptive method that achieves very high accuracy by
	 * combining the trapezoidal rule with Richardson extrapolation. It builds a tableau
	 * of estimates and extrapolates to the limit h→0, eliminating error terms successively.
	 *
	 * ALGORITHM:
	 * 1. Compute successive trapezoidal approximations T(h), T(h/2), T(h/4), ...
	 * 2. Apply Richardson extrapolation (polynomial interpolation) to extrapolate to h=0
	 * 3. The extrapolation uses Neville's algorithm on the sequence h², (h/2)², (h/4)², ...
	 * 4. Convergence is detected when the extrapolated value stabilizes
	 *
	 * CONVERGENCE:
	 * For analytic functions, convergence is exponentially fast: O(h^(2K)) where K is
	 * the number of extrapolation steps (default K=5). This makes Romberg one of the
	 * fastest methods for smooth, analytic integrands.
	 *
	 * ADVANTAGES:
	 * - Very fast convergence for smooth functions (faster than Simpson)
	 * - Achieves high accuracy with relatively few function evaluations
	 * - Built-in error estimation from extrapolation
	 *
	 * DISADVANTAGES:
	 * - Requires smooth (analytic) integrand for optimal performance
	 * - Poor performance on functions with discontinuous derivatives
	 * - More complex than simple quadrature rules
	 *
	 * USE CASES:
	 * - High-precision integration of smooth functions
	 * - Scientific computing where accuracy is paramount
	 * - When function evaluations are expensive but function is analytic
	 *
	 * COMPLEXITY: O(2^n) function evaluations for n refinement steps
	 *
	 * @param func The function to integrate
	 * @param a Lower bound of integration
	 * @param b Upper bound of integration
	 * @param eps Desired relative accuracy (default: Defaults::RombergIntegrationEPS)
	 * @return IntegrationResult with value, error estimate, iterations, convergence status
	 *
	 * @note Reference: Numerical Recipes §4.3 "Romberg Integration"
	 */
	static IntegrationResult IntegrateRomberg(const IRealFunction& func, Real a, Real b,
																						const Real eps = Defaults::RombergIntegrationEPS)
	{
		const int JMAX = Defaults::RombergIntegrationMaxSteps;
		const int K = Defaults::RombergIntegrationUsedPnts;  // Number of points for polynomial extrapolation
		
		std::vector<Real> s(JMAX);      // Successive trapezoidal approximations
		std::vector<Real> h(JMAX + 1);  // Step sizes squared: h[j] = (initial_h / 2^j)^2
		
		h[0] = 1.0;
		TrapIntegrator trap(func, a, b);
		
		Real ss = 0.0;  // Extrapolated result
		Real dss = 0.0; // Error estimate from extrapolation
		
		for (int j = 0; j < JMAX; j++)
		{
			s[j] = trap.next();
			
			if (j >= K - 1)
			{
				// Perform polynomial extrapolation using Neville's algorithm
				// Extrapolate to h=0 using the last K points
				std::vector<Real> c(K), d(K);
				
				int ns = 0;
				Real dif = std::abs(h[j - K + 1]);
				
				// Initialize c and d arrays with the last K values of s
				for (int i = 0; i < K; i++)
				{
					Real dift = std::abs(h[j - K + 1 + i]);
					if (dift < dif)
					{
						ns = i;
						dif = dift;
					}
					c[i] = s[j - K + 1 + i];
					d[i] = s[j - K + 1 + i];
				}
				
				// Initial best guess
				ss = s[j - K + 1 + ns];
				ns--;
				
				// Neville's algorithm: build up the extrapolation tableau
				for (int m = 1; m < K; m++)
				{
					for (int i = 0; i < K - m; i++)
					{
						Real ho = h[j - K + 1 + i];
						Real hp = h[j - K + 1 + i + m];
						Real w = c[i + 1] - d[i];
						Real den = ho - hp;
						
						if (den == 0.0)
						{
							// This should not happen with proper h values
							return IntegrationResult(ss, std::abs(dss), j + 1, false);
						}
						
						den = w / den;
						d[i] = hp * den;
						c[i] = ho * den;
					}
					
					// Decide which correction to add (from c or d)
					if (2 * (ns + 1) < (K - m))
						dss = c[ns + 1];
					else
						dss = d[ns--];
					
					ss += dss;
				}
				
				// Check for convergence
				if (std::abs(dss) <= eps * std::abs(ss) || (ss == 0.0 && dss == 0.0))
				{
					return IntegrationResult(ss, std::abs(dss), j + 1, true);
				}
			}
			
			// Prepare for next iteration: h_new = h_old / 4 (because we're extrapolating in h²)
			h[j + 1] = 0.25 * h[j];
		}
		
		// Did not converge within max iterations
		return IntegrationResult(ss, std::abs(dss), JMAX, false);
	}


	/**
	 * @brief IntegrateGauss10 - Ten-point Gauss-Legendre quadrature.
	 *
	 * Gauss-Legendre quadrature achieves the highest algebraic accuracy possible
	 * for a given number of function evaluations. With 10 points, it exactly
	 * integrates polynomials up to degree 19 (2n-1 rule).
	 *
	 * ALGORITHM:
	 * Uses pre-computed abscissas (zeros of Legendre polynomial P_10) and weights
	 * to evaluate: ∫f(x)dx ≈ Σ w_i * f(x_i)
	 *
	 * The nodes and weights are transformed from [-1,1] to [a,b] via linear mapping.
	 *
	 * CHARACTERISTICS:
	 * - Fixed-order method: exactly 10 function evaluations per call
	 * - No adaptive refinement or error estimation
	 * - Extremely fast for smooth functions
	 * - Exact for polynomials of degree ≤ 19
	 *
	 * ADVANTAGES:
	 * - Very fast (only 10 function evaluations)
	 * - High accuracy for smooth functions
	 * - No iteration overhead
	 * - Deterministic behavior
	 *
	 * DISADVANTAGES:
	 * - No built-in error estimation (returned error_estimate is 0)
	 * - Cannot adapt to function complexity
	 * - May be inaccurate for highly oscillatory or singular functions
	 * - Fixed precision - cannot be improved without using more points
	 *
	 * USE CASES:
	 * - Quick integration of smooth functions
	 * - Inner loops where speed is critical
	 * - Polynomial or near-polynomial integrands
	 * - When approximate error bounds are known a priori
	 *
	 * @param func The function to integrate
	 * @param a Lower bound of integration
	 * @param b Upper bound of integration
	 * @return IntegrationResult with value, 0 error estimate (not available), 1 iteration, converged=true
	 *
	 * @note Error estimate is 0 because this is a non-adaptive method.
	 *       Use adaptive methods (Romberg, Simpson) when error bounds are needed.
	 * @note Reference: Abramowitz & Stegun, Table 25.4 for nodes and weights
	 */
	static IntegrationResult IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)
	{
		// Pre-computed Gauss-Legendre nodes (positive only, symmetric about 0)
		static const Real x[] = { 0.1488743389816312, 0.4333953941292472,
		                          0.6794095682990244, 0.8650633666889845, 0.9739065285171717 };
		// Pre-computed Gauss-Legendre weights (for the corresponding nodes)
		static const Real w[] = { 0.2955242247147529, 0.2692667193099963,
		                          0.2190863625159821, 0.1494513491505806, 0.0666713443086881 };
		
		// Transform from [-1,1] to [a,b]: x_mapped = xm + xr * x_i
		Real xm = 0.5 * (b + a);  // Midpoint
		Real xr = 0.5 * (b - a);  // Half-width (scale factor)
		
		Real s = 0;
		for (int j = 0; j < 5; j++) 
		{
			Real dx = xr * x[j];
			s += w[j] * (func(xm + dx) + func(xm - dx));
		}
		Real result = s * xr;
		
		// Return IntegrationResult for API consistency
		// Note: Gauss quadrature has no built-in error estimate (non-adaptive)
		// iterations=1 indicates single-pass evaluation (10 function calls)
		return IntegrationResult(result, 0.0, 1, true);
	}

	/**
	 * @brief 21-point Gauss-Kronrod integration (G10K21 rule)
	 * @details Uses 21-point Kronrod rule with embedded 10-point Gauss for error estimation.
	 *          Provides built-in error estimate from difference between Gauss and Kronrod results.
	 *          More accurate than GAUSS10 with only 21 function evaluations.
	 * @param func Function to integrate
	 * @param a Lower bound of integration
	 * @param b Upper bound of integration
	 * @return IntegrationResult with value, error_estimate, iterations=21, converged=true
	 */
	static IntegrationResult IntegrateGK21(const IRealFunction& func, Real a, Real b)
	{
		// Wrap in lambda to avoid passing abstract class by value to template
		auto gkResult = Integration::IntegrateGK21([&func](Real x) { return func(x); }, a, b);
		return IntegrationResult(gkResult.value, gkResult.error_estimate, gkResult.function_evals, gkResult.converged);
	}

	/**
	 * @brief Unified integration function with explicit method selection.
	 *
	 * This is the recommended way to perform 1D integration - explicitly specify
	 * the integration method rather than relying on global state.
	 *
	 * @param func The function to integrate
	 * @param a Lower bound of integration
	 * @param b Upper bound of integration
	 * @param method Integration method to use (default: TRAP)
	 * @param eps Desired relative accuracy (not used for GAUSS10)
	 * @return IntegrationResult with value, error estimate, iterations, convergence status
	 *
	 * Example:
	 * @code
	 *   auto result = Integrate(f, 0.0, 1.0, SIMPSON, 1e-8);
	 *   if (!result.converged) { ... }
	 * @endcode
	 */
	static IntegrationResult Integrate(const IRealFunction& func, Real a, Real b,
	                                   IntegrationMethod method = TRAP,
	                                   Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		switch (method)
		{
		case TRAP:
			return IntegrateTrap(func, a, b, eps);
		case SIMPSON:
			return IntegrateSimpson(func, a, b, eps);
		case ROMBERG:
			return IntegrateRomberg(func, a, b, eps);
		case GAUSS10:
			return IntegrateGauss10(func, a, b);
		case GAUSS10KRONROD21:
			return IntegrateGK21(func, a, b);
		default:
			return IntegrateTrap(func, a, b, eps);
		}
	}

	/**
	 * @brief Template-based integration with compile-time method selection.
	 *
	 * Zero-overhead method selection using compile-time dispatch.
	 * Use this when the method is known at compile time for best performance.
	 *
	 * @tparam Method Integration method (compile-time constant)
	 * @param func The function to integrate
	 * @param a Lower bound of integration
	 * @param b Upper bound of integration
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, iterations, convergence status
	 *
	 * Example:
	 * @code
	 *   auto result = Integrate<SIMPSON>(f, 0.0, 1.0, 1e-10);
	 * @endcode
	 */
	template<IntegrationMethod Method = TRAP>
	static IntegrationResult Integrate(const IRealFunction& func, Real a, Real b,
	                                   Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		if constexpr (Method == TRAP)
			return IntegrateTrap(func, a, b, eps);
		else if constexpr (Method == SIMPSON)
			return IntegrateSimpson(func, a, b, eps);
		else if constexpr (Method == ROMBERG)
			return IntegrateRomberg(func, a, b, eps);
		else if constexpr (Method == GAUSS10)
			return IntegrateGauss10(func, a, b);
		else if constexpr (Method == GAUSS10KRONROD21)
			return IntegrateGK21(func, a, b);
		else
			return IntegrateTrap(func, a, b, eps);
	}
}

#endif // MML_INTEGRATION_1D_H