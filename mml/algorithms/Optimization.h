///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization.h                                                      ///
///  Description: 1D optimization algorithms (Golden Section, Brent's method)         ///
///               Line search and bracketing utilities                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_H
#define MML_OPTIMIZATION_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"
#include "base/Function.h"

namespace MML
{
	/**
	 * Custom exception for optimization algorithm failures.
	 */
	class OptimizationError : public std::runtime_error
	{
		int _iterations;
	public:
		OptimizationError(std::string inMessage, int iterations = 0) 
			: std::runtime_error(inMessage), _iterations(iterations) { }
		
		int iterations() const noexcept { return _iterations; }
	};

	/**
	 * @brief Struct to hold results from minimization algorithms.
	 * 
	 * Contains the location of the minimum, the function value at that point,
	 * the number of iterations used, and convergence status.
	 */
	struct MinimizationResult
	{
		Real xmin;              ///< Location of the minimum
		Real fmin;              ///< Function value at the minimum
		int iterations;         ///< Number of iterations performed
		bool converged;         ///< Whether the algorithm converged
		
		MinimizationResult() : xmin(0), fmin(0), iterations(0), converged(false) {}
		MinimizationResult(Real x, Real f, int iter, bool conv = true) 
			: xmin(x), fmin(f), iterations(iter), converged(conv) {}
	};

	/**
	 * @brief Struct to hold a bracket around a minimum.
	 * 
	 * Contains three points (ax, bx, cx) and their function values (fa, fb, fc)
	 * such that bx is between ax and cx, and fb < min(fa, fc).
	 * This guarantees a minimum exists in the interval [ax, cx].
	 */
	struct MinimumBracket
	{
		Real ax, bx, cx;        ///< Three points bracketing the minimum
		Real fa, fb, fc;        ///< Function values at those points
		bool valid;             ///< Whether a valid bracket was found
		
		MinimumBracket() : ax(0), bx(0), cx(0), fa(0), fb(0), fc(0), valid(false) {}
	};

	/**
	 * @brief Namespace containing function minimization algorithms.
	 * 
	 * Implements algorithms from Numerical Recipes Chapter 10:
	 * - BracketMinimum: Initially bracketing a minimum (§10.1)
	 * - GoldenSectionSearch: Golden section search in one dimension (§10.2)
	 * - BrentMinimize: Parabolic interpolation and Brent's method (§10.3)
	 * - BrentMinimizeWithDeriv: One-dimensional search with first derivatives (§10.4)
	 * 
	 * USAGE PATTERN:
	 * @code
	 * // 1. Define or create a function to minimize
	 * RealFunctionFromStdFunc func([](Real x) { return (x-2)*(x-2) + 1; });
	 * 
	 * // 2. Find initial bracket
	 * MinimumBracket bracket = Minimization::BracketMinimum(func, 0.0, 1.0);
	 * 
	 * // 3. Refine using Brent's method
	 * MinimizationResult result = Minimization::BrentMinimize(func, bracket);
	 * 
	 * // Result: xmin ≈ 2.0, fmin ≈ 1.0
	 * @endcode
	 */
	namespace Minimization
	{
		/*********************************************************************/
		/*****                  10.1 Bracket Minimum                     *****/
		/*********************************************************************/
		/**
		 * @brief Find an initial bracket containing a minimum of a function.
		 * 
		 * Given an initial interval [a, b], this function searches for a triplet
		 * (ax, bx, cx) such that fb < fa and fb < fc, guaranteeing that a minimum
		 * exists between ax and cx (for a continuous function).
		 * 
		 * ALGORITHM:
		 * The algorithm uses the golden ratio (≈1.618) for geometric expansion and
		 * parabolic interpolation to efficiently find a bracketing triplet:
		 * 
		 * 1. Start with two points a, b and evaluate f(a) and f(b)
		 * 2. If f(b) > f(a), swap so we're going "downhill"
		 * 3. Use golden ratio to generate a third point c beyond b
		 * 4. While f(b) > f(c) (still going downhill):
		 *    - Try parabolic extrapolation for next point
		 *    - Apply various safeguards to prevent runaway expansion
		 *    - Update triplet and continue
		 * 5. Return when f(b) < f(c), ensuring bracket is valid
		 * 
		 * GOLDEN RATIO:
		 * The golden ratio φ = (1 + √5)/2 ≈ 1.618034 is used because:
		 * - It provides optimal interval expansion rate
		 * - Parabolic interpolation through three points spaced by φ tends to
		 *   converge well for smooth functions
		 * 
		 * PARABOLIC INTERPOLATION:
		 * When applicable, the algorithm fits a parabola through the three points
		 * and uses its minimum as the next guess. This can dramatically accelerate
		 * bracketing for smooth, unimodal functions.
		 * 
		 * SAFEGUARDS:
		 * - GLIMIT (100x interval) prevents unbounded expansion
		 * - TINY prevents division by zero in parabolic step
		 * - Various geometric checks ensure monotonic progress
		 * 
		 * COMPLEXITY:
		 * - Typically O(log(distance_to_minimum)) function evaluations
		 * - Worst case: function may evaluate many times if minimum is very far
		 * 
		 * @param func   Real function to minimize
		 * @param a      Initial left point of search interval
		 * @param b      Initial right point of search interval  
		 * @return MinimumBracket containing the bracket (ax,bx,cx) with function values
		 * 
		 * POSTCONDITION (if valid):
		 * - bx is between ax and cx
		 * - fb <= fa and fb <= fc
		 * - A minimum exists in [ax, cx]
		 * 
		 * @note Based on Numerical Recipes §10.1 (mnbrak)
		 * @see GoldenSectionSearch, BrentMinimize for methods to refine the bracket
		 */
		static MinimumBracket BracketMinimum(const IRealFunction& func, Real a, Real b)
		{
			const Real GOLD = 1.618034;    // Golden ratio
			const Real GLIMIT = 100.0;     // Maximum magnification for parabolic step
			const Real TINY = 1.0e-20;     // Small number to prevent division by zero

			MinimumBracket result;
			result.ax = a;
			result.bx = b;
			
			result.fa = func(result.ax);
			result.fb = func(result.bx);
			
			// Ensure we're going downhill from ax to bx
			if (result.fb > result.fa)
			{
				std::swap(result.ax, result.bx);
				std::swap(result.fa, result.fb);
			}
			
			// First guess for cx using golden ratio
			result.cx = result.bx + GOLD * (result.bx - result.ax);
			result.fc = func(result.cx);
			
			Real fu;
			
			// Keep bracketing until we have fb < fc (minimum bracketed)
			while (result.fb > result.fc)
			{
				// Parabolic extrapolation to find potential minimum
				Real r = (result.bx - result.ax) * (result.fb - result.fc);
				Real q = (result.bx - result.cx) * (result.fb - result.fa);
				Real u = result.bx - ((result.bx - result.cx) * q - (result.bx - result.ax) * r) /
				        (2.0 * std::copysign(std::max(std::abs(q - r), TINY), q - r));
				Real ulim = result.bx + GLIMIT * (result.cx - result.bx);
				
				// Test various possibilities for u
				if ((result.bx - u) * (u - result.cx) > 0.0)
				{
					// u is between bx and cx
					fu = func(u);
					if (fu < result.fc)
					{
						// Minimum between bx and cx
						result.ax = result.bx;
						result.bx = u;
						result.fa = result.fb;
						result.fb = fu;
						result.valid = true;
						return result;
					}
					else if (fu > result.fb)
					{
						// Minimum between ax and u
						result.cx = u;
						result.fc = fu;
						result.valid = true;
						return result;
					}
					// Parabolic fit was useless, use golden section
					u = result.cx + GOLD * (result.cx - result.bx);
					fu = func(u);
				}
				else if ((result.cx - u) * (u - ulim) > 0.0)
				{
					// u is between cx and limit
					fu = func(u);
					if (fu < result.fc)
					{
						// Shift triplet
						result.bx = result.cx;
						result.cx = u;
						u = result.cx + GOLD * (result.cx - result.bx);
						result.fb = result.fc;
						result.fc = fu;
						fu = func(u);
					}
				}
				else if ((u - ulim) * (ulim - result.cx) >= 0.0)
				{
					// u is beyond the limit
					u = ulim;
					fu = func(u);
				}
				else
				{
					// Reject parabolic u, use golden section
					u = result.cx + GOLD * (result.cx - result.bx);
					fu = func(u);
				}
				
				// Shift triplet (ax, bx, cx) <- (bx, cx, u)
				result.ax = result.bx;
				result.bx = result.cx;
				result.cx = u;
				result.fa = result.fb;
				result.fb = result.fc;
				result.fc = fu;
			}
			
			result.valid = true;
			return result;
		}

		/*********************************************************************/
		/*****          10.2 Golden Section Search                       *****/
		/*********************************************************************/
		/**
		 * @brief Find minimum using golden section search.
		 * 
		 * Given a bracketing triplet (ax, bx, cx) where fb < fa and fb < fc,
		 * this algorithm successively narrows the bracket until the minimum
		 * is found to the desired tolerance.
		 * 
		 * ALGORITHM:
		 * Golden section search is based on dividing intervals using the golden ratio
		 * φ = (√5-1)/2 ≈ 0.618034. At each step:
		 * 
		 * 1. The largest interval is divided by the golden ratio
		 * 2. A new point is evaluated, creating two overlapping subintervals
		 * 3. The interval containing the minimum is kept
		 * 4. Process repeats until interval width < tolerance
		 * 
		 * WHY GOLDEN RATIO?
		 * The golden ratio has a unique property: if you divide an interval by φ,
		 * you can reuse one of the evaluation points in the next iteration.
		 * This means only ONE new function evaluation per iteration!
		 * 
		 * Alternative ratio (like bisection) would require TWO evaluations.
		 * Golden section achieves the same O(log N) convergence with half the work.
		 * 
		 * CONVERGENCE:
		 * - Linear convergence with ratio 0.618
		 * - After k iterations, interval reduced by factor φ^k
		 * - For 10^-8 tolerance: ~40 iterations needed
		 * - Each iteration requires exactly 1 function evaluation
		 * 
		 * COMPLEXITY:
		 * - O(log(1/tol)) function evaluations
		 * - Guaranteed convergence for any unimodal function
		 * 
		 * WHEN TO USE:
		 * - Function is expensive to evaluate
		 * - Derivatives not available
		 * - Function may not be smooth enough for Brent's method
		 * - Simplicity is more important than raw speed
		 * 
		 * @param func    Real function to minimize
		 * @param bracket Pre-computed bracket from BracketMinimum
		 * @param tol     Fractional tolerance (default: 3e-8, ≈ √(machine epsilon))
		 * @return MinimizationResult with xmin, fmin, and iteration count
		 * 
		 * @note Based on Numerical Recipes §10.2
		 * @see BracketMinimum to obtain initial bracket
		 * @see BrentMinimize for typically faster convergence
		 */
		static MinimizationResult GoldenSectionSearch(const IRealFunction& func, 
		                                              const MinimumBracket& bracket,
		                                              Real tol = 3.0e-8)
		{
			const Real R = 0.61803399;  // Golden ratio conjugate (1 - φ) = (√5-1)/2
			const Real C = 1.0 - R;     // = 1/φ = 2/(√5+1)
			
			MinimizationResult result;
			result.converged = false;
			result.iterations = 0;
			
			Real x0 = bracket.ax;
			Real x3 = bracket.cx;
			Real x1, x2;
			
			// Make x0 to x3 the smaller segment
			if (std::abs(bracket.cx - bracket.bx) > std::abs(bracket.bx - bracket.ax))
			{
				x1 = bracket.bx;
				x2 = bracket.bx + C * (bracket.cx - bracket.bx);
			}
			else
			{
				x2 = bracket.bx;
				x1 = bracket.bx - C * (bracket.bx - bracket.ax);
			}
			
			Real f1 = func(x1);
			Real f2 = func(x2);
			
			// Main iteration loop
			while (std::abs(x3 - x0) > tol * (std::abs(x1) + std::abs(x2)))
			{
				result.iterations++;
				
				if (f2 < f1)
				{
					// Minimum is in [x1, x3], narrow from left
					x0 = x1;
					x1 = x2;
					x2 = R * x2 + C * x3;
					f1 = f2;
					f2 = func(x2);
				}
				else
				{
					// Minimum is in [x0, x2], narrow from right
					x3 = x2;
					x2 = x1;
					x1 = R * x1 + C * x0;
					f2 = f1;
					f1 = func(x1);
				}
			}
			
			// Return best point found
			if (f1 < f2)
			{
				result.xmin = x1;
				result.fmin = f1;
			}
			else
			{
				result.xmin = x2;
				result.fmin = f2;
			}
			
			result.converged = true;
			return result;
		}

		/**
		 * @brief Convenience function combining bracketing and golden section search.
		 * 
		 * @param func  Real function to minimize
		 * @param a     Initial left search bound
		 * @param b     Initial right search bound
		 * @param tol   Fractional tolerance
		 * @return MinimizationResult with xmin, fmin, and iteration count
		 */
		static MinimizationResult GoldenSectionSearch(const IRealFunction& func,
		                                              Real a, Real b,
		                                              Real tol = 3.0e-8)
		{
			MinimumBracket bracket = BracketMinimum(func, a, b);
			if (!bracket.valid)
				throw OptimizationError("Failed to bracket minimum in GoldenSectionSearch");
			return GoldenSectionSearch(func, bracket, tol);
		}

		/*********************************************************************/
		/*****      10.3 Brent's Method (Parabolic Interpolation)        *****/
		/*********************************************************************/
		/**
		 * @brief Find minimum using Brent's method with parabolic interpolation.
		 * 
		 * Brent's method combines the reliability of golden section search with
		 * the speed of parabolic interpolation. It's the standard choice for
		 * one-dimensional minimization when derivatives are not available.
		 * 
		 * ALGORITHM:
		 * Brent's method maintains a bracketing interval [a,b] and tracks three
		 * points: x (best point so far), w (second best), v (previous w).
		 * 
		 * At each iteration:
		 * 1. Fit a parabola through x, w, v
		 * 2. If parabolic step is:
		 *    - Within the bracket AND
		 *    - Less than half the previous step AND
		 *    - Not too close to endpoints
		 *    → Accept parabolic step (fast quadratic convergence)
		 * 3. Otherwise:
		 *    → Use golden section step (guaranteed progress)
		 * 
		 * CONVERGENCE:
		 * - Superlinear convergence when parabolic steps accepted
		 * - Near-quadratic for smooth functions
		 * - Falls back to linear (golden section) for difficult functions
		 * - Never worse than golden section
		 * 
		 * SAFEGUARDS:
		 * - Parabolic step must stay within bracket
		 * - Step must be decreasing (not just bouncing)
		 * - Minimum step size prevents infinite loops near tolerance
		 * - Track e (extent of movement) to detect stalling
		 * 
		 * WHY BRENT'S METHOD?
		 * - Best of both worlds: fast when possible, robust always
		 * - No function derivatives needed
		 * - Guaranteed convergence for continuous functions
		 * - Industry standard for 1D minimization
		 * 
		 * COMPLEXITY:
		 * - Best case: O(log log(1/tol)) for smooth unimodal functions
		 * - Worst case: O(log(1/tol)) - same as golden section
		 * - Typically 2-3x faster than pure golden section
		 * 
		 * @param func    Real function to minimize
		 * @param bracket Pre-computed bracket from BracketMinimum
		 * @param tol     Fractional tolerance (default: 3e-8)
		 * @param maxIter Maximum iterations (default: 100)
		 * @return MinimizationResult with xmin, fmin, and iteration count
		 * 
		 * @throws OptimizationError if maximum iterations exceeded
		 * 
		 * @note Based on Numerical Recipes §10.3
		 * @see BracketMinimum to obtain initial bracket
		 * @see BrentMinimizeWithDeriv for version using derivatives
		 */
		static MinimizationResult BrentMinimize(const IRealFunction& func,
		                                        const MinimumBracket& bracket,
		                                        Real tol = 3.0e-8,
		                                        int maxIter = 100)
		{
			const Real CGOLD = 0.3819660;  // Golden ratio for section (1 - 1/φ)
			const Real ZEPS = std::numeric_limits<Real>::epsilon() * 1.0e-3;
			
			MinimizationResult result;
			result.converged = false;
			result.iterations = 0;
			
			Real a, b, d = 0.0, etemp, fu, fv, fw, fx;
			Real p, q, r, tol1, tol2, u, v, w, x, xm;
			Real e = 0.0;  // Distance moved on the step before last
			
			// a and b must be in ascending order
			a = (bracket.ax < bracket.cx) ? bracket.ax : bracket.cx;
			b = (bracket.ax > bracket.cx) ? bracket.ax : bracket.cx;
			
			// Initialize x, w, v to the best point (bx from bracket)
			x = w = v = bracket.bx;
			fw = fv = fx = func(x);
			
			for (int iter = 0; iter < maxIter; iter++)
			{
				result.iterations++;
				xm = 0.5 * (a + b);  // Midpoint
				tol1 = tol * std::abs(x) + ZEPS;
				tol2 = 2.0 * tol1;
				
				// Test for convergence
				if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a)))
				{
					result.xmin = x;
					result.fmin = fx;
					result.converged = true;
					return result;
				}
				
				// Try parabolic interpolation
				if (std::abs(e) > tol1)
				{
					// Fit parabola
					r = (x - w) * (fx - fv);
					q = (x - v) * (fx - fw);
					p = (x - v) * q - (x - w) * r;
					q = 2.0 * (q - r);
					if (q > 0.0) p = -p;
					q = std::abs(q);
					etemp = e;
					e = d;
					
					// Accept parabolic step?
					if (std::abs(p) >= std::abs(0.5 * q * etemp) || 
					    p <= q * (a - x) || p >= q * (b - x))
					{
						// No, use golden section instead
						e = (x >= xm) ? a - x : b - x;
						d = CGOLD * e;
					}
					else
					{
						// Yes, take parabolic step
						d = p / q;
						u = x + d;
						// Don't step too close to boundaries
						if (u - a < tol2 || b - u < tol2)
							d = std::copysign(tol1, xm - x);
					}
				}
				else
				{
					// Golden section step
					e = (x >= xm) ? a - x : b - x;
					d = CGOLD * e;
				}
				
				// Function evaluation (never closer than tol1)
				u = (std::abs(d) >= tol1) ? x + d : x + std::copysign(tol1, d);
				fu = func(u);
				
				// Update bracket and best points
				if (fu <= fx)
				{
					if (u >= x) a = x; else b = x;
					// Shift (v, w, x) <- (w, x, u)
					v = w; w = x; x = u;
					fv = fw; fw = fx; fx = fu;
				}
				else
				{
					if (u < x) a = u; else b = u;
					if (fu <= fw || w == x)
					{
						v = w; w = u;
						fv = fw; fw = fu;
					}
					else if (fu <= fv || v == x || v == w)
					{
						v = u;
						fv = fu;
					}
				}
			}
			
			// Failed to converge
			result.xmin = x;
			result.fmin = fx;
			result.converged = false;
			throw OptimizationError("Maximum iterations exceeded in BrentMinimize", maxIter);
		}

		/**
		 * @brief Convenience function combining bracketing and Brent's method.
		 * 
		 * @param func    Real function to minimize
		 * @param a       Initial left search bound
		 * @param b       Initial right search bound
		 * @param tol     Fractional tolerance
		 * @param maxIter Maximum iterations
		 * @return MinimizationResult with xmin, fmin, and iteration count
		 */
		static MinimizationResult BrentMinimize(const IRealFunction& func,
		                                        Real a, Real b,
		                                        Real tol = 3.0e-8,
		                                        int maxIter = 100)
		{
			MinimumBracket bracket = BracketMinimum(func, a, b);
			if (!bracket.valid)
				throw OptimizationError("Failed to bracket minimum in BrentMinimize");
			return BrentMinimize(func, bracket, tol, maxIter);
		}

		/*********************************************************************/
		/*****   10.4 One-Dimensional Search with First Derivatives      *****/
		/*********************************************************************/
		/**
		 * @brief Find minimum using Brent's method with derivative information.
		 * 
		 * When derivatives are available, this method can converge faster than
		 * standard Brent's method by using derivative information to guide
		 * the search direction and perform secant-like updates.
		 * 
		 * ALGORITHM:
		 * Similar to BrentMinimize but additionally tracks derivatives:
		 * 
		 * 1. Maintains bracket [a,b] and three points (x,w,v) with function
		 *    values AND derivatives
		 * 2. Uses secant method with derivatives for extrapolation:
		 *    - d1 = derivative-based step from (x,dx) and (w,dw)
		 *    - d2 = derivative-based step from (x,dx) and (v,dv)
		 * 3. Takes smaller step if it's acceptable, else bisection
		 * 4. Uses sign of derivative to determine which half of bracket to keep
		 * 
		 * DERIVATIVE USAGE:
		 * - Derivative sign indicates direction to minimum
		 * - If f'(x) < 0: minimum is to the right of x
		 * - If f'(x) > 0: minimum is to the left of x
		 * - Secant update: x_new = x - f'(x) * (x - x_prev) / (f'(x) - f'(x_prev))
		 * 
		 * CONVERGENCE:
		 * - Superlinear (order ≈ 1.618) when secant steps accepted
		 * - Faster than derivative-free Brent for smooth functions
		 * - Falls back to bisection when needed
		 * 
		 * WHEN TO USE:
		 * - Derivative is cheap to compute (automatic differentiation, closed form)
		 * - Function is smooth with well-behaved derivative
		 * - Maximum speed is important
		 * 
		 * WHEN NOT TO USE:
		 * - Derivative is expensive (numerical differentiation)
		 * - Function has discontinuous derivatives
		 * - Simplicity is preferred over speed
		 * 
		 * @param func    Real function to minimize (must have operator())
		 * @param dfunc   Derivative function (df/dx)
		 * @param bracket Pre-computed bracket from BracketMinimum
		 * @param tol     Fractional tolerance (default: 3e-8)
		 * @param maxIter Maximum iterations (default: 100)
		 * @return MinimizationResult with xmin, fmin, and iteration count
		 * 
		 * @throws OptimizationError if maximum iterations exceeded
		 * 
		 * @note Based on Numerical Recipes §10.4 (dbrent)
		 * @see BrentMinimize for version without derivatives
		 */
		static MinimizationResult BrentMinimizeWithDeriv(const IRealFunction& func,
		                                                  const IRealFunction& dfunc,
		                                                  const MinimumBracket& bracket,
		                                                  Real tol = 3.0e-8,
		                                                  int maxIter = 100)
		{
			const Real ZEPS = std::numeric_limits<Real>::epsilon() * 1.0e-3;
			
			MinimizationResult result;
			result.converged = false;
			result.iterations = 0;
			
			bool ok1, ok2;
			Real a, b, d = 0.0, d1, d2, du, dv, dw, dx, e = 0.0;
			Real fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;
			
			// a and b must be in ascending order
			a = (bracket.ax < bracket.cx) ? bracket.ax : bracket.cx;
			b = (bracket.ax > bracket.cx) ? bracket.ax : bracket.cx;
			
			// Initialize x, w, v to the best point
			x = w = v = bracket.bx;
			fw = fv = fx = func(x);
			dw = dv = dx = dfunc(x);
			
			for (int iter = 0; iter < maxIter; iter++)
			{
				result.iterations++;
				xm = 0.5 * (a + b);
				tol1 = tol * std::abs(x) + ZEPS;
				tol2 = 2.0 * tol1;
				
				// Test for convergence
				if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a)))
				{
					result.xmin = x;
					result.fmin = fx;
					result.converged = true;
					return result;
				}
				
				if (std::abs(e) > tol1)
				{
					// Use derivatives for secant extrapolation
					d1 = 2.0 * (b - a);
					d2 = d1;
					
					if (dw != dx) d1 = (w - x) * dx / (dx - dw);
					if (dv != dx) d2 = (v - x) * dx / (dx - dv);
					
					u1 = x + d1;
					u2 = x + d2;
					
					// Check if either step is acceptable
					ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
					ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
					
					olde = e;
					e = d;
					
					if (ok1 || ok2)
					{
						// Choose the smaller step if both are okay
						if (ok1 && ok2)
							d = (std::abs(d1) < std::abs(d2)) ? d1 : d2;
						else if (ok1)
							d = d1;
						else
							d = d2;
						
						if (std::abs(d) <= std::abs(0.5 * olde))
						{
							u = x + d;
							if (u - a < tol2 || b - u < tol2)
								d = std::copysign(tol1, xm - x);
						}
						else
						{
							// Bisection
							e = (dx >= 0.0) ? a - x : b - x;
							d = 0.5 * e;
						}
					}
					else
					{
						// Bisection
						e = (dx >= 0.0) ? a - x : b - x;
						d = 0.5 * e;
					}
				}
				else
				{
					// Bisection
					e = (dx >= 0.0) ? a - x : b - x;
					d = 0.5 * e;
				}
				
				// Function evaluation
				if (std::abs(d) >= tol1)
				{
					u = x + d;
					fu = func(u);
				}
				else
				{
					u = x + std::copysign(tol1, d);
					fu = func(u);
					if (fu > fx)
					{
						// Minimum bracketed, return
						result.xmin = x;
						result.fmin = fx;
						result.converged = true;
						return result;
					}
				}
				
				du = dfunc(u);
				
				// Update bracket and best points
				if (fu <= fx)
				{
					if (u >= x) a = x; else b = x;
					// Shift (v,w,x) <- (w,x,u)
					v = w; fv = fw; dv = dw;
					w = x; fw = fx; dw = dx;
					x = u; fx = fu; dx = du;
				}
				else
				{
					if (u < x) a = u; else b = u;
					if (fu <= fw || w == x)
					{
						v = w; fv = fw; dv = dw;
						w = u; fw = fu; dw = du;
					}
					else if (fu < fv || v == x || v == w)
					{
						v = u; fv = fu; dv = du;
					}
				}
			}
			
			// Failed to converge
			result.xmin = x;
			result.fmin = fx;
			result.converged = false;
			throw OptimizationError("Maximum iterations exceeded in BrentMinimizeWithDeriv", maxIter);
		}

		/**
		 * @brief Convenience function combining bracketing and derivative-based Brent.
		 * 
		 * @param func    Real function to minimize
		 * @param dfunc   Derivative function
		 * @param a       Initial left search bound
		 * @param b       Initial right search bound
		 * @param tol     Fractional tolerance
		 * @param maxIter Maximum iterations
		 * @return MinimizationResult with xmin, fmin, and iteration count
		 */
		static MinimizationResult BrentMinimizeWithDeriv(const IRealFunction& func,
		                                                  const IRealFunction& dfunc,
		                                                  Real a, Real b,
		                                                  Real tol = 3.0e-8,
		                                                  int maxIter = 100)
		{
			MinimumBracket bracket = BracketMinimum(func, a, b);
			if (!bracket.valid)
				throw OptimizationError("Failed to bracket minimum in BrentMinimizeWithDeriv");
			return BrentMinimizeWithDeriv(func, dfunc, bracket, tol, maxIter);
		}

		/*********************************************************************/
		/*****                 Maximization wrappers                     *****/
		/*********************************************************************/
		/**
		 * @brief Find maximum using golden section search (negates function).
		 * 
		 * @param func  Real function to maximize
		 * @param a     Initial left search bound
		 * @param b     Initial right search bound
		 * @param tol   Fractional tolerance
		 * @return MinimizationResult with xmax, fmax (note: fmax is actual max value)
		 */
		static MinimizationResult GoldenSectionMaximize(const IRealFunction& func,
		                                                 Real a, Real b,
		                                                 Real tol = 3.0e-8)
		{
			// Negate function to turn maximum into minimum
			RealFunctionFromStdFunc negFunc([&func](Real x) { return -func(x); });
			MinimizationResult result = GoldenSectionSearch(negFunc, a, b, tol);
			result.fmin = -result.fmin;  // Restore original function value
			return result;
		}

		/**
		 * @brief Find maximum using Brent's method (negates function).
		 * 
		 * @param func    Real function to maximize
		 * @param a       Initial left search bound
		 * @param b       Initial right search bound
		 * @param tol     Fractional tolerance
		 * @param maxIter Maximum iterations
		 * @return MinimizationResult with xmax, fmax
		 */
		static MinimizationResult BrentMaximize(const IRealFunction& func,
		                                         Real a, Real b,
		                                         Real tol = 3.0e-8,
		                                         int maxIter = 100)
		{
			RealFunctionFromStdFunc negFunc([&func](Real x) { return -func(x); });
			MinimizationResult result = BrentMinimize(negFunc, a, b, tol, maxIter);
			result.fmin = -result.fmin;
			return result;
		}

	}  // namespace Minimization
}  // namespace MML

#endif // MML_OPTIMIZATION_H
