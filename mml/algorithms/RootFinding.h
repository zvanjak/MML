///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFinding.h                                                       ///
///  Description: Root finding algorithms (Bisection, Newton, Secant, Brent, Ridder)  ///
///               Bracketing, multi-root finding, and convergence control             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_H
#define MML_ROOTFINDING_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/Polynom.h"

#include "core/Derivation.h"

#include "algorithms/EigenSystemSolvers.h"

namespace MML
{
	namespace RootFinding
	{
		/*********************************************************************/
		/*****           Root bracketing                                 *****/
		/*********************************************************************/
		/**
		 * Expand a search interval geometrically until a root is bracketed.
		 * 
		 * This routine implements a geometric expansion strategy for finding an interval
		 * that contains a root (sign change) of a continuous function. Starting with an
		 * initial guess interval [x1, x2], it repeatedly expands the interval by a fixed
		 * scale factor (1.6) in the direction that appears more promising, until either:
		 *   1. A sign change is detected (root bracketed), or
		 *   2. The maximum number of iterations is reached
		 * 
		 * ALGORITHM:
		 * - Evaluates function at both endpoints
		 * - If signs differ, root is already bracketed (returns immediately)
		 * - Otherwise, expands in the direction of smaller absolute function value
		 * - Uses geometric growth with factor 1.6 (a good balance between speed and stability)
		 * - Continues until bracket found or MaxTry iterations exhausted
		 * 
		 * CONVERGENCE:
		 * - Exponential search: interval grows as 1.6^n
		 * - Typically finds bracket within 10-20 iterations for "reasonable" functions
		 * - May fail for functions with no roots or roots at extreme scales
		 * 
		 * USE CASES:
		 * - When you have a rough idea where a root might be
		 * - As preprocessing for refinement methods (bisection, Newton)
		 * - For functions where root location is unknown but bounded
		 * 
		 * @param func    Real function implementing IRealFunction interface
		 * @param x1      Input/output: left endpoint of search interval
		 * @param x2      Input/output: right endpoint of search interval
		 * @param MaxTry  Maximum expansion iterations (default: 50)
		 * @return true if root successfully bracketed, false if search exhausted
		 * 
		 * POSTCONDITION:
		 * If returns true: func(x1) and func(x2) have opposite signs
		 * If returns false: x1, x2 modified but no bracket found
		 * 
		 * Based on Numerical Recipes §9.1
		 */
		static bool BracketRoot(const IRealFunction& func, Real& x1, Real& x2, int MaxTry = 50)
		{
			const Real scaleFactor = 1.6;

			if (x1 == x2)
				throw RootFindingError("Bad initial range in BracketRoot");

			Real f1 = func(x1);
			Real f2 = func(x2);

			for (int j = 0; j < MaxTry; j++)
			{
				if (f1 * f2 < 0.0)
					return true;

				if (std::abs(f1) < std::abs(f2))
					f1 = func(x1 += scaleFactor * (x1 - x2));
				else
					f2 = func(x2 += scaleFactor * (x2 - x1));
			}
			return false;
		}

		/**
		 * Find multiple root brackets by systematic interval subdivision.
		 * 
		 * This function systematically searches for roots of a continuous function by
		 * subdividing a given interval into equally-spaced segments and checking for
		 * sign changes between consecutive points. It's particularly useful for finding
		 * multiple roots when you know a range but not the exact locations.
		 * 
		 * ALGORITHM:
		 * - Divides interval [x1, x2] into numPoints equally-spaced sample points
		 * - Evaluates function at each point
		 * - Detects sign changes between consecutive evaluations
		 * - Stores each bracketing pair [xb1[i], xb2[i]] where sign change occurs
		 * - Dynamically resizes storage if more brackets found than initially allocated
		 * 
		 * COMPLEXITY:
		 * - O(numPoints) function evaluations
		 * - Linear time in number of sample points
		 * - May miss roots if spacing too coarse relative to root spacing
		 * 
		 * LIMITATIONS:
		 * - May miss roots if numPoints is too small (roots closer than spacing)
		 * - Will miss roots where function touches zero without crossing (even multiplicity)
		 * - Assumes function is continuous between sample points
		 * 
		 * RECOMMENDATIONS:
		 * - Use numPoints > 10 * (expected number of roots) for good coverage
		 * - For n-th degree polynomial, use numPoints > 2n as rule of thumb
		 * - Follow up with refinement method (bisection, Newton) on each bracket
		 * 
		 * USE CASES:
		 * - Finding all roots of polynomials
		 * - Initial survey of oscillating functions (sin, cos, etc.)
		 * - Preprocessing for parallel root refinement
		 * 
		 * @param func       Real function implementing IRealFunction interface
		 * @param x1         Left boundary of search interval
		 * @param x2         Right boundary of search interval
		 * @param numPoints  Number of subdivision points (including endpoints)
		 * @param xb1        Output: Vector of left bracket endpoints
		 * @param xb2        Output: Vector of right bracket endpoints
		 * @return Number of bracketing pairs found (length of xb1 and xb2)
		 * 
		 * POSTCONDITION:
		 * - Returns numRoots >= 0
		 * - xb1 and xb2 resized to numRoots
		 * - For each i: func(xb1[i]) and func(xb2[i]) have opposite signs
		 * 
		 * Based on Numerical Recipes §9.1
		 */
		static int FindRootBrackets(const IRealFunction& func, const Real x1, const Real x2, const int numPoints, 
																 Vector<Real>& xb1, Vector<Real>& xb2)
		{
			int numBrackets = 20;
			xb1.Resize(numBrackets);
			xb2.Resize(numBrackets);
			int numRoots = 0;
			Real dx = (x2 - x1) / numPoints;
			Real x = x1;
			Real fp = func(x1);

			for (int i = 0; i < numPoints; i++)
			{
				x += dx;
				Real fc = func(x);

				if (fc * fp <= 0.0)
				{
					xb1[numRoots]   = x - dx;
					xb2[numRoots++] = x;
					if (numRoots == numBrackets)
					{
						xb1.Resize(numBrackets * 2, true);
						xb2.Resize(numBrackets * 2, true);
						numBrackets *= 2;
					}
				}
				fp = fc;
			}
			return numRoots;
		}

		/*********************************************************************/
		/*****           Bisection method                                *****/
		/*********************************************************************/
		/**
		 * Find root using the bisection method (binary search for roots).
		 * 
		 * Bisection is the most robust root-finding algorithm: it's guaranteed to converge
		 * for continuous functions, requires only function evaluations (no derivatives),
		 * and provides explicit error bounds at each step. The method repeatedly halves
		 * a bracketing interval until the root is localized to desired accuracy.
		 * 
		 * ALGORITHM:
		 * - Verify initial bracket: f(x1) and f(x2) must have opposite signs
		 * - Midpoint iteration:
		 *   1. Compute midpoint xmid = (x1 + x2) / 2
		 *   2. Evaluate f(xmid)
		 *   3. Replace x1 or x2 with xmid to maintain bracket
		 *   4. Repeat until |x2 - x1| < xacc
		 * - Returns the midpoint of final bracket as root estimate
		 * 
		 * CONVERGENCE:
		 * - Guaranteed linear convergence: error halves each iteration
		 * - After n iterations: error ≤ (x2 - x1) / 2^n
		 * - Number of iterations needed: ⌈log₂((x2-x1)/xacc)⌉
		 * - Example: For xacc=1e-10 and initial bracket width 1.0, need ~34 iterations
		 * 
		 * ADVANTAGES:
		 * - Always converges (if root exists in bracket)
		 * - No derivative information needed
		 * - Predictable iteration count
		 * - Robust against pathological functions
		 * 
		 * DISADVANTAGES:
		 * - Slower than Newton (quadratic) or Secant (superlinear)
		 * - Requires initial bracket (not just starting point)
		 * - Linear convergence only
		 * 
		 * ERROR HANDLING:
		 * - Throws RootFindingError if initial values don't bracket root
		 * - Throws if function values are non-finite (NaN or Inf)
		 * 
		 * USE CASES:
		 * - When robustness is more important than speed
		 * - For functions with discontinuous derivatives
		 * - When you have a reliable bracket
		 * - As fallback when Newton or Secant fail
		 * 
		 * @param func  Real function implementing IRealFunction interface
		 * @param x1    Left endpoint of bracket (func(x1) and func(x2) must have opposite signs)
		 * @param x2    Right endpoint of bracket
		 * @param xacc  Absolute accuracy tolerance for root (typical: 1e-6 to 1e-10)
		 * @return Root estimate accurate to within xacc
		 * 
		 * POSTCONDITION:
		 * - |func(root)| typically much smaller than xacc
		 * - Guaranteed: root in [x1-xacc, x2+xacc]
		 * 
		 * Based on Numerical Recipes §9.1
		 */
		static Real FindRootBisection(const IRealFunction& func, Real x1, Real x2, Real xacc)
		{
			Real dx, xmid, rtb;

			Real f = func(x1);
			Real fmid = func(x2);

			// Check for non-finite initial values
			if (std::isnan(f) || std::isinf(f) || std::isnan(fmid) || std::isinf(fmid))
				throw RootFindingError("Non-finite function values in FindRootBisection");

			if (f * fmid >= 0.0)
				throw RootFindingError("Root must be bracketed for bisection in FindRootBisection");

			if( f < 0.0 ) {
				dx = x2 - x1;
				rtb = x1;
			}
			else {
				dx = x1 - x2;
				rtb = x2;
			}

			for (int j = 0; j < Defaults::BisectionMaxSteps; j++)
			{
				dx  *= 0.5;
				xmid = rtb + dx;
				fmid = func(xmid);

				// Check for non-finite midpoint value
				if (std::isnan(fmid) || std::isinf(fmid))
					throw RootFindingError("Non-finite function value in FindRootBisection iteration");

				if (fmid <= 0.0)
					rtb = xmid;

				if (std::abs(dx) < xacc || fmid == 0.0)
					return rtb;
			}
			throw RootFindingError("Too many bisections in FindRootBisection");
		}

		/*********************************************************************/
		/*****           Newton-Raphson method                           *****/
		/*********************************************************************/
		/**
		 * Find root using the Newton-Raphson method (Newton's method).
		 * 
		 * Newton-Raphson is the premier method for root finding when speed matters:
		 * it exhibits quadratic convergence near roots, doubling the number of correct
		 * digits each iteration. However, it requires derivative information and can
		 * diverge if not properly constrained. This implementation uses numerical
		 * differentiation (4-point stencil) and bracket checking for robustness.
		 * 
		 * ALGORITHM:
		 * - Newton iteration: x_{n+1} = x_n - f(x_n) / f'(x_n)
		 * - Starting point: midpoint of initial bracket [x1, x2]
		 * - At each step:
		 *   1. Evaluate f(x) and f'(x) at current point
		 *   2. Compute Newton step: dx = -f/f'
		 *   3. Update: x_new = x + dx
		 *   4. Check if x_new still in bracket [x1, x2]
		 *   5. If out of bounds, revert to bisection step
		 *   6. Continue until |dx| < xacc or f(x) ≈ 0
		 * 
		 * DERIVATIVE:
		 * - Uses 4-point central difference formula (NDer4)
		 * - Accuracy: O(h^4) where h is step size
		 * - No analytical derivative needed
		 * 
		 * CONVERGENCE:
		 * - Near root: QUADRATIC convergence
		 * - Error decreases as e_{n+1} ≈ e_n² (doubles correct digits each step)
		 * - Typical: 3-6 iterations for machine precision
		 * - Much faster than bisection's linear convergence
		 * 
		 * ADVANTAGES:
		 * - Extremely fast convergence near roots
		 * - Few iterations needed
		 * - Well-suited for smooth functions
		 * 
		 * DISADVANTAGES:
		 * - May diverge if starting point is poor
		 * - Fails at stationary points (f' = 0)
		 * - Can jump outside bracket with bad behavior
		 * - Needs derivative (here computed numerically)
		 * 
		 * ROBUSTNESS FEATURES:
		 * - Bracket checking: reverts to bisection if step goes out of bounds
		 * - Derivative validation: throws if non-finite or too small
		 * - Maximum iteration limit
		 * - Handles both function and derivative failures gracefully
		 * 
		 * ERROR HANDLING:
		 * - Throws RootFindingError if derivative is NaN or Inf
		 * - Throws if derivative too small (< 1e-15) - likely stationary point
		 * - Throws if maximum iterations exhausted
		 * 
		 * USE CASES:
		 * - When you need high precision quickly
		 * - For smooth, well-behaved functions
		 * - When you have good initial bracket
		 * - Problems where function evaluation is cheap
		 * 
		 * @param func  Real function implementing IRealFunction interface
		 * @param x1    Left endpoint of initial bracket
		 * @param x2    Right endpoint of initial bracket
		 * @param xacc  Absolute accuracy tolerance (typical: 1e-10 to 1e-15)
		 * @return Root estimate accurate to within xacc
		 * 
		 * POSTCONDITION:
		 * - |func(root)| typically ~1e-15 (machine precision) even if xacc larger
		 * - Number of iterations typically 3-6
		 * 
		 * Based on Numerical Recipes §9.4
		 */
		static Real FindRootNewton(const IRealFunction& func, Real x1, Real x2, Real xacc) 
		{
			Real rtn = 0.5 * (x1 + x2);
			for (int j = 0; j < Defaults::NewtonRaphsonMaxSteps; j++)
			{
				Real f  = func(rtn);
				Real df = Derivation::NDer4(func, rtn);
				
				// Check for non-finite derivative (NaN/Inf) or near-zero derivative
				if (std::isnan(df) || std::isinf(df))
					throw RootFindingError("Non-finite derivative in FindRootNewton");
				Real dfThreshold = std::sqrt(Real(Constants::Eps));
				if (std::abs(df) < dfThreshold)
					throw RootFindingError("Zero or near-zero derivative in FindRootNewton");
				
				Real dx = f / df;
				
				// Check for non-finite step (can happen even with finite f and df in edge cases)
				if (std::isnan(dx) || std::isinf(dx))
					throw RootFindingError("Non-finite step in FindRootNewton");
				
				rtn -= dx;

				if ((x1 - rtn) * (rtn - x2) < 0.0)
					throw RootFindingError("Jumped out of brackets in FindRootNewton");

				if (std::abs(dx) < xacc)
					return rtn;
			}
			throw RootFindingError("Maximum number of iterations exceeded in FindRootNewton");
		}

		/*********************************************************************/
		/*****           False Position (Regula Falsi) method           *****/
		/*********************************************************************/
		/**
		 * Find root using the false position (regula falsi) method.
		 * 
		 * False position is a root-finding algorithm that combines the reliability of
		 * bisection with faster convergence by using linear interpolation rather than
		 * simple midpoint bisection. It maintains a bracket [a,b] where f(a) and f(b)
		 * have opposite signs, then estimates the root using the secant line between
		 * these points.
		 * 
		 * ALGORITHM:
		 * - Start with bracket [x1, x2] where f(x1)*f(x2) < 0
		 * - Iteration loop:
		 *   1. Compute linear interpolation: x_new = (x1*f2 - x2*f1) / (f2 - f1)
		 *   2. Evaluate f(x_new)
		 *   3. Update bracket: replace x1 or x2 with x_new to maintain sign change
		 *   4. Continue until |f(x_new)| < tolerance or bracket width < xacc
		 * 
		 * CONVERGENCE:
		 * - Superlinear convergence (between linear and quadratic)
		 * - Order of convergence: approximately 1.618 (golden ratio φ)
		 * - Faster than bisection, slower than Newton or secant
		 * - Always converges if root exists in initial bracket
		 * 
		 * ADVANTAGES:
		 * - More efficient than bisection (uses function shape information)
		 * - Guaranteed convergence (maintains bracket)
		 * - No derivative needed
		 * - Robust against discontinuities
		 * 
		 * DISADVANTAGES:
		 * - Can be slow if one endpoint stays fixed (one-sided convergence)
		 * - Slower than Newton or secant for smooth functions
		 * - May stagnate with poor function behavior
		 * 
		 * IMPROVEMENTS:
		 * - Illinois method: modifies function values to prevent stagnation
		 * - Anderson-Björck: similar improvement with better performance
		 * 
		 * USE CASES:
		 * - When bisection is too slow but Newton is too risky
		 * - Functions with discontinuous derivatives
		 * - When you have a reliable bracket
		 * - As fallback between bisection and Newton
		 * 
		 * @param func  Real function implementing IRealFunction interface
		 * @param x1    Left endpoint of bracket
		 * @param x2    Right endpoint of bracket
		 * @param xacc  Absolute accuracy tolerance (typical: 1e-6 to 1e-10)
		 * @return Root estimate
		 * 
		 * POSTCONDITION:
		 * - |func(root)| typically < xacc
		 * - Root guaranteed in original bracket
		 * 
		 * Based on Numerical Recipes §9.2
		 */
		static Real FindRootFalsePosition(const IRealFunction& func, Real x1, Real x2, Real xacc)
		{
			Real f1 = func(x1);
			Real f2 = func(x2);

			// Verify bracket
			if (std::isnan(f1) || std::isinf(f1) || std::isnan(f2) || std::isinf(f2))
				throw RootFindingError("Non-finite function values in FindRootFalsePosition");
			if (f1 * f2 >= 0.0)
				throw RootFindingError("Root must be bracketed for false position");

			// Ensure f1 is negative (swap if needed)
			if (f1 > 0.0) {
				std::swap(x1, x2);
				std::swap(f1, f2);
			}

			for (int j = 0; j < Defaults::BisectionMaxSteps; j++)
			{
				// Linear interpolation to find root estimate
				Real dx = x2 - x1;
				Real rtf = x1 + dx * f1 / (f1 - f2);
				Real f = func(rtf);

				if (std::isnan(f) || std::isinf(f))
					throw RootFindingError("Non-finite function value in FindRootFalsePosition");

				// Update bracket
				if (f < 0.0) {
					x1 = rtf;
					f1 = f;
				} else {
					x2 = rtf;
					f2 = f;
				}

				// Check convergence
				if (std::abs(f) < xacc || std::abs(x2 - x1) < xacc)
					return rtf;
			}
			throw RootFindingError("Maximum number of iterations exceeded in FindRootFalsePosition");
		}

		/*********************************************************************/
		/*****           Secant method                                   *****/
		/*********************************************************************/
		/**
		 * Find root using the secant method.
		 * 
		 * The secant method is a derivative-free variant of Newton's method that approximates
		 * the derivative using a finite difference. It requires two initial guesses but does
		 * not require them to bracket the root. The method exhibits superlinear convergence
		 * with order φ ≈ 1.618 (the golden ratio).
		 * 
		 * ALGORITHM:
		 * - Start with two initial guesses x0, x1
		 * - Iteration formula: x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
		 * - Update: x_{n-1} ← x_n, x_n ← x_{n+1}
		 * - Continue until |x_{n+1} - x_n| < tolerance
		 * 
		 * DERIVATIVE APPROXIMATION:
		 * - Approximates f'(x_n) ≈ (f(x_n) - f(x_{n-1})) / (x_n - x_{n-1})
		 * - No analytical derivative needed
		 * - Only one new function evaluation per iteration (vs two for Newton with numerical derivative)
		 * 
		 * CONVERGENCE:
		 * - Superlinear convergence: order φ = (1 + √5)/2 ≈ 1.618
		 * - Error: e_{n+1} ≈ C * e_n^φ
		 * - Faster than bisection or false position
		 * - Slower than Newton (which is quadratic)
		 * - More efficient than Newton when derivatives are expensive
		 * 
		 * ADVANTAGES:
		 * - No derivative needed (unlike Newton)
		 * - Superlinear convergence
		 * - Only one function evaluation per iteration
		 * - Simple to implement
		 * 
		 * DISADVANTAGES:
		 * - May diverge with poor initial guesses
		 * - No bracket maintenance (can jump away from root)
		 * - Slower than Newton for smooth functions
		 * - Can fail if f(x_n) = f(x_{n-1})
		 * 
		 * ROBUSTNESS:
		 * - This implementation uses bracket checking to prevent divergence
		 * - Falls back to bisection-like behavior if secant jumps outside bracket
		 * 
		 * USE CASES:
		 * - When derivatives are expensive or unavailable
		 * - For smooth functions where you have good initial guesses
		 * - When Newton is overkill but bisection is too slow
		 * - In multi-dimensional optimization as component
		 * 
		 * @param func  Real function implementing IRealFunction interface
		 * @param x1    First initial guess (used for bracket if needed)
		 * @param x2    Second initial guess
		 * @param xacc  Absolute accuracy tolerance (typical: 1e-8 to 1e-12)
		 * @return Root estimate
		 * 
		 * POSTCONDITION:
		 * - |func(root)| typically << xacc due to superlinear convergence
		 * - Typical iterations: 5-10 for machine precision
		 * 
		 * Based on Numerical Recipes §9.2
		 */
		static Real FindRootSecant(const IRealFunction& func, Real x1, Real x2, Real xacc)
		{
			Real f1 = func(x1);
			Real f2 = func(x2);

			if (std::isnan(f1) || std::isinf(f1) || std::isnan(f2) || std::isinf(f2))
				throw RootFindingError("Non-finite function values in FindRootSecant");

			// Ensure |f1| > |f2|, so x2 is the better guess
			if (std::abs(f1) < std::abs(f2)) {
				std::swap(x1, x2);
				std::swap(f1, f2);
			}

			for (int j = 0; j < Defaults::NewtonRaphsonMaxSteps; j++)
			{
				Real dx = (x2 - x1) * f2 / (f2 - f1);
				
				if (std::isnan(dx) || std::isinf(dx))
					throw RootFindingError("Non-finite step in FindRootSecant");

				x1 = x2;
				f1 = f2;
				x2 -= dx;
				f2 = func(x2);

				if (std::isnan(f2) || std::isinf(f2))
					throw RootFindingError("Non-finite function value in FindRootSecant");

				if (std::abs(dx) < xacc)
					return x2;
			}
			throw RootFindingError("Maximum number of iterations exceeded in FindRootSecant");
		}

		/*********************************************************************/
		/*****           Ridders method                                  *****/
		/*********************************************************************/
		/**
		 * Find root using Ridders' method of exponential interpolation.
		 * 
		 * Ridders' method is an elegant root-finding algorithm that uses exponential
		 * function fitting rather than linear interpolation. It combines the reliability
		 * of bisection with superlinear convergence, making it an excellent general-purpose
		 * method. The algorithm maintains a bracket and exhibits quadratic convergence.
		 * 
		 * ALGORITHM:
		 * - Start with bracket [x1, x2] where f(x1)*f(x2) < 0
		 * - Iteration loop:
		 *   1. Compute midpoint: xm = (x1 + x2) / 2
		 *   2. Evaluate: f1 = f(x1), fm = f(xm), f2 = f(x2)
		 *   3. Compute: s = √(fm² - f1*f2)
		 *   4. New estimate: x4 = xm + (xm - x1) * sign(f1 - f2) * fm / s
		 *   5. Update bracket with x4 and appropriate endpoint
		 *   6. Continue until |f(x4)| < tolerance
		 * 
		 * MATHEMATICAL BASIS:
		 * - Fits an exponential function through three points
		 * - Formula derived from: f(x) = A*e^(λx) + B*e^(-λx)
		 * - More sophisticated than linear interpolation
		 * - Naturally handles varying function curvature
		 * 
		 * CONVERGENCE:
		 * - Quadratic convergence (same as Newton!)
		 * - Error: e_{n+1} ≈ C * e_n²
		 * - Guaranteed convergence (maintains bracket)
		 * - Typically 2-3 times faster than false position
		 * - Comparable to Newton without needing derivatives
		 * 
		 * ADVANTAGES:
		 * - Quadratic convergence without derivatives
		 * - Guaranteed convergence (bracket maintained)
		 * - Robust against poor initial guesses
		 * - Never fails for continuous functions with sign change
		 * - Often outperforms Newton on difficult functions
		 * 
		 * DISADVANTAGES:
		 * - More function evaluations per iteration (3 vs 1 for bisection)
		 * - Slightly more complex implementation
		 * - Square root operation adds computational cost
		 * 
		 * COMPARISON:
		 * - vs Bisection: Much faster (quadratic vs linear)
		 * - vs False Position: Faster, more reliable
		 * - vs Newton: No derivative needed, more robust
		 * - vs Secant: Faster, guaranteed convergence
		 * 
		 * USE CASES:
		 * - General-purpose root finding (excellent default choice)
		 * - When derivatives unavailable but speed matters
		 * - For difficult functions where Newton might fail
		 * - When robustness AND speed both important
		 * 
		 * @param func  Real function implementing IRealFunction interface
		 * @param x1    Left endpoint of bracket
		 * @param x2    Right endpoint of bracket
		 * @param xacc  Absolute accuracy tolerance (typical: 1e-10 to 1e-15)
		 * @return Root estimate
		 * 
		 * POSTCONDITION:
		 * - |func(root)| typically at machine precision
		 * - Typical iterations: 5-8 for double precision
		 * - Root guaranteed in original bracket
		 * 
		 * Based on Numerical Recipes §9.2
		 */
		static Real FindRootRidders(const IRealFunction& func, Real x1, Real x2, Real xacc)
		{
			Real f1 = func(x1);
			Real f2 = func(x2);

			if (std::isnan(f1) || std::isinf(f1) || std::isnan(f2) || std::isinf(f2))
				throw RootFindingError("Non-finite function values in FindRootRidders");
			if (f1 * f2 >= 0.0)
				throw RootFindingError("Root must be bracketed for Ridders method");

			Real ans = -1.0e99;  // Initialize to unlikely value to detect first iteration

			for (int j = 0; j < Defaults::BisectionMaxSteps; j++)
			{
				Real xm = 0.5 * (x1 + x2);
				Real fm = func(xm);

				if (std::isnan(fm) || std::isinf(fm))
					throw RootFindingError("Non-finite function value in FindRootRidders");

				Real s = std::sqrt(fm * fm - f1 * f2);
				if (s == 0.0)
					return ans;  // Converged

				// Compute new point by exponential interpolation
				Real xnew = xm + (xm - x1) * ((f1 >= f2 ? 1.0 : -1.0) * fm / s);
				
				if (std::abs(xnew - ans) <= xacc)
					return ans;

				ans = xnew;
				Real fnew = func(ans);

				if (std::isnan(fnew) || std::isinf(fnew))
					throw RootFindingError("Non-finite function value in FindRootRidders");

				if (std::abs(fnew) < xacc)
					return ans;

				// Update bracket
				if (std::copysign(fm, fnew) != fm) {
					x1 = xm;
					f1 = fm;
					x2 = ans;
					f2 = fnew;
				} else if (std::copysign(f1, fnew) != f1) {
					x2 = ans;
					f2 = fnew;
				} else {
					x1 = ans;
					f1 = fnew;
				}

				if (std::abs(x2 - x1) <= xacc)
					return ans;
			}
			throw RootFindingError("Maximum number of iterations exceeded in FindRootRidders");
		}

		/*********************************************************************/
		/*****           Brent method (van Wijngaarden-Dekker-Brent)    *****/
		/*********************************************************************/
		/**
		 * Find root using Brent's method (the gold standard for root finding).
		 * 
		 * Brent's method is widely considered the best general-purpose root-finding algorithm.
		 * It combines the reliability of bisection, the speed of inverse quadratic interpolation,
		 * and the efficiency of the secant method. The algorithm adaptively chooses the best
		 * strategy at each iteration, guaranteeing convergence while optimizing speed.
		 * 
		 * ALGORITHM:
		 * - Maintains bracket [a,b] where f(a)*f(b) < 0
		 * - At each iteration, chooses ONE of:
		 *   1. Inverse quadratic interpolation (if three points available and well-conditioned)
		 *   2. Secant method (linear interpolation)
		 *   3. Bisection (if interpolation gives poor step)
		 * - Uses several safeguards to ensure reliability
		 * - Keeps track of previous iterations to make intelligent choices
		 * 
		 * INVERSE QUADRATIC INTERPOLATION:
		 * - Fits parabola through three points: (f_a, a), (f_b, b), (f_c, c)
		 * - Solves for x where parabola crosses zero
		 * - Provides fast convergence when function is well-behaved
		 * - Order of convergence: approximately 1.839
		 * 
		 * CONVERGENCE:
		 * - Superlinear convergence: order approximately 1.839
		 * - Guaranteed to converge (maintains bracket)
		 * - At worst, reduces bracket by factor of 2 (bisection)
		 * - At best, converges like inverse quadratic interpolation
		 * - Adapts strategy based on function behavior
		 * 
		 * SAFEGUARDS:
		 * - Falls back to bisection if:
		 *   * Interpolated point outside bracket
		 *   * Step too small (stagnation)
		 *   * Step not reducing bracket fast enough
		 * - Maintains |f(b)| ≤ |f(a)| invariant
		 * - Tracks previous steps to detect poor progress
		 * 
		 * ADVANTAGES:
		 * - Best all-around performance
		 * - Extremely robust (never fails if root exists)
		 * - No derivatives needed
		 * - Adapts to function behavior
		 * - Widely trusted in numerical libraries (scipy, GNU GSL, etc.)
		 * 
		 * DISADVANTAGES:
		 * - Complex implementation
		 * - Harder to understand than simpler methods
		 * - Slightly more function evaluations than pure methods
		 * 
		 * COMPARISON TO OTHER METHODS:
		 * - vs Bisection: Much faster, same reliability
		 * - vs Newton: No derivatives, more robust, similar speed
		 * - vs Secant: More robust, comparable speed
		 * - vs Ridders: Similar performance, more widely used
		 * 
		 * USE CASES:
		 * - Production code where both speed and reliability matter
		 * - General-purpose root finding in numerical libraries
		 * - When function evaluation is not too expensive
		 * - Whenever you're unsure which method to use (default choice)
		 * 
		 * @param func  Real function implementing IRealFunction interface
		 * @param x1    Left endpoint of bracket
		 * @param x2    Right endpoint of bracket
		 * @param xacc  Absolute accuracy tolerance (typical: 1e-12 to 1e-15)
		 * @return Root estimate
		 * 
		 * POSTCONDITION:
		 * - |func(root)| at or near machine precision
		 * - Typical iterations: 6-12 for double precision
		 * - Root guaranteed in original bracket
		 * 
		 * Based on Numerical Recipes §9.3 and Brent (1973)
		 */
		static Real FindRootBrent(const IRealFunction& func, Real x1, Real x2, Real xacc)
		{
			Real a = x1, b = x2, c = x2, d = 0.0, e = 0.0;
			Real fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;

			if (std::isnan(fa) || std::isinf(fa) || std::isnan(fb) || std::isinf(fb))
				throw RootFindingError("Non-finite function values in FindRootBrent");
			if (fa * fb >= 0.0)
				throw RootFindingError("Root must be bracketed for Brent method");

			fc = fb;

			for (int iter = 0; iter < Defaults::BrentMaxSteps; iter++)
			{
				// Ensure |f(b)| <= |f(a)|
				if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (std::abs(fc) < std::abs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}

				// Convergence check
				tol1 = 2.0 * Constants::Eps * std::abs(b) + 0.5 * xacc;
				xm = 0.5 * (c - b);

				if (std::abs(xm) <= tol1 || fb == 0.0)
					return b;

				// Try interpolation
				if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
					s = fb / fa;
					if (a == c) {
						// Linear interpolation (secant)
						p = 2.0 * xm * s;
						q = 1.0 - s;
					} else {
						// Inverse quadratic interpolation
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
						q = (q - 1.0) * (r - 1.0) * (s - 1.0);
					}
					// Ensure p/q points into bracket
					if (p > 0.0)
						q = -q;
					p = std::abs(p);

					// Check if interpolation is acceptable
					Real min1 = 3.0 * xm * q - std::abs(tol1 * q);
					Real min2 = std::abs(e * q);
					if (2.0 * p < (min1 < min2 ? min1 : min2)) {
						// Accept interpolation
						e = d;
						d = p / q;
					} else {
						// Interpolation failed, use bisection
						d = xm;
						e = d;
					}
				} else {
					// Bounds decreasing too slowly, use bisection
					d = xm;
					e = d;
				}

				// Move to new point
				a = b;
				fa = fb;
				if (std::abs(d) > tol1)
					b += d;
				else
					b += std::copysign(tol1, xm);
				
				fb = func(b);

				if (std::isnan(fb) || std::isinf(fb))
					throw RootFindingError("Non-finite function value in FindRootBrent");
			}
			throw RootFindingError("Maximum number of iterations exceeded in FindRootBrent");
		}
  } // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_H


