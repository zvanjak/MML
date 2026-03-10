///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingBracketing.h                                             ///
///  Description: Bracket-based root finding algorithms                               ///
///               Bracketing utilities, Bisection, and False Position methods         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_BRACKETING_H
#define MML_ROOTFINDING_BRACKETING_H

#include "mml/MMLBase.h"
#include "mml/core/NumericValidation.h"

#include "mml/interfaces/IFunction.h"
#include "mml/base/Vector/Vector.h"

#include "mml/algorithms/RootFinding/RootFindingBase.h"

namespace MML {
	namespace RootFinding {

		/*********************************************************************/
		/*****           Forward Declarations                            *****/
		/*********************************************************************/
		// Config-based overloads declared here so simple versions can call them
		static RootFindingResult FindRootBisection(const IRealFunction& func, Real x1, Real x2, 
												   const RootFindingConfig& config);
		static RootFindingResult FindRootFalsePosition(const IRealFunction& func, Real x1, Real x2,
													   const RootFindingConfig& config);

		/*********************************************************************/
		/*****           Root bracketing                                 *****/
		/*********************************************************************/
		/// Expand a search interval geometrically until a root is bracketed.
		///
		/// This routine implements a geometric expansion strategy for finding an interval
		/// that contains a root (sign change) of a continuous function. Starting with an
		/// initial guess interval [x1, x2], it repeatedly expands the interval by a fixed
		/// scale factor (1.6) in the direction that appears more promising, until either:
		/// 1. A sign change is detected (root bracketed), or
		/// 2. The maximum number of iterations is reached
		///
		/// ALGORITHM:
		/// - Evaluates function at both endpoints
		/// - If signs differ, root is already bracketed (returns immediately)
		/// - Otherwise, expands in the direction of smaller absolute function value
		/// - Uses geometric growth with factor 1.6 (a good balance between speed and stability)
		/// - Continues until bracket found or MaxTry iterations exhausted
		///
		/// CONVERGENCE:
		/// - Exponential search: interval grows as 1.6^n
		/// - Typically finds bracket within 10-20 iterations for "reasonable" functions
		/// - May fail for functions with no roots or roots at extreme scales
		///
		/// USE CASES:
		/// - When you have a rough idea where a root might be
		/// - As preprocessing for refinement methods (bisection, Newton)
		/// - For functions where root location is unknown but bounded
		///
		/// @param func    Real function implementing IRealFunction interface
		/// @param x1      Input/output: left endpoint of search interval
		/// @param x2      Input/output: right endpoint of search interval
		/// @param MaxTry  Maximum expansion iterations (default: 50)
		/// @return true if root successfully bracketed, false if search exhausted
		///
		/// POSTCONDITION:
		/// If returns true: func(x1) and func(x2) have opposite signs
		/// If returns false: x1, x2 modified but no bracket found
		///
		/// Based on Numerical Recipes §9.1
		static bool BracketRoot(const IRealFunction& func, Real& x1, Real& x2, int MaxTry = 50) {
			const Real scaleFactor = 1.6;

			if (x1 == x2)
				throw RootFindingError("Bad initial range in BracketRoot");

			Real f1 = func(x1);
			Real f2 = func(x2);

			for (int j = 0; j < MaxTry; j++) {
				if (f1 * f2 < 0.0)
					return true;

				if (std::abs(f1) < std::abs(f2))
					f1 = func(x1 += scaleFactor * (x1 - x2));
				else
					f2 = func(x2 += scaleFactor * (x2 - x1));
			}
			return false;
		}

		/// Find multiple root brackets by systematic interval subdivision.
		///
		/// This function systematically searches for roots of a continuous function by
		/// subdividing a given interval into equally-spaced segments and checking for
		/// sign changes between consecutive points. It's particularly useful for finding
		/// multiple roots when you know a range but not the exact locations.
		///
		/// ALGORITHM:
		/// - Divides interval [x1, x2] into numPoints equally-spaced sample points
		/// - Evaluates function at each point
		/// - Detects sign changes between consecutive evaluations
		/// - Stores each bracketing pair [xb1[i], xb2[i]] where sign change occurs
		/// - Dynamically resizes storage if more brackets found than initially allocated
		///
		/// COMPLEXITY:
		/// - O(numPoints) function evaluations
		/// - Linear time in number of sample points
		/// - May miss roots if spacing too coarse relative to root spacing
		///
		/// LIMITATIONS:
		/// - May miss roots if numPoints is too small (roots closer than spacing)
		/// - Will miss roots where function touches zero without crossing (even multiplicity)
		/// - Assumes function is continuous between sample points
		///
		/// RECOMMENDATIONS:
		/// - Use numPoints > 10 * (expected number of roots) for good coverage
		/// - For n-th degree polynomial, use numPoints > 2n as rule of thumb
		/// - Follow up with refinement method (bisection, Newton) on each bracket
		///
		/// USE CASES:
		/// - Finding all roots of polynomials
		/// - Initial survey of oscillating functions (sin, cos, etc.)
		/// - Preprocessing for parallel root refinement
		///
		/// @param func       Real function implementing IRealFunction interface
		/// @param x1         Left boundary of search interval
		/// @param x2         Right boundary of search interval
		/// @param numPoints  Number of subdivision points (including endpoints)
		/// @param xb1        Output: Vector of left bracket endpoints
		/// @param xb2        Output: Vector of right bracket endpoints
		/// @return Number of bracketing pairs found (length of xb1 and xb2)
		///
		/// POSTCONDITION:
		/// - Returns numRoots >= 0
		/// - xb1 and xb2 resized to numRoots
		/// - For each i: func(xb1[i]) and func(xb2[i]) have opposite signs
		///
		/// Based on Numerical Recipes §9.1
		static int FindRootBrackets(const IRealFunction& func, const Real x1, const Real x2, const int numPoints, Vector<Real>& xb1,
									Vector<Real>& xb2) {
			int numBrackets = 20;
			xb1.Resize(numBrackets);
			xb2.Resize(numBrackets);
			int numRoots = 0;
			Real dx = (x2 - x1) / numPoints;
			Real x = x1;
			Real fp = func(x1);

			for (int i = 0; i < numPoints; i++) {
				x += dx;
				Real fc = func(x);

				if (fc * fp <= 0.0) {
					xb1[numRoots] = x - dx;
					xb2[numRoots++] = x;
					if (numRoots == numBrackets) {
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
		/// Find root using the bisection method (binary search for roots).
		///
		/// Bisection is the most robust root-finding algorithm: it's guaranteed to converge
		/// for continuous functions, requires only function evaluations (no derivatives),
		/// and provides explicit error bounds at each step. The method repeatedly halves
		/// a bracketing interval until the root is localized to desired accuracy.
		///
		/// ALGORITHM:
		/// - Verify initial bracket: f(x1) and f(x2) must have opposite signs
		/// - Midpoint iteration:
		/// 1. Compute midpoint xmid = (x1 + x2) / 2
		/// 2. Evaluate f(xmid)
		/// 3. Replace x1 or x2 with xmid to maintain bracket
		/// 4. Repeat until |x2 - x1| < xacc
		/// - Returns the midpoint of final bracket as root estimate
		///
		/// CONVERGENCE:
		/// - Guaranteed linear convergence: error halves each iteration
		/// - After n iterations: error ≤ (x2 - x1) / 2^n
		/// - Number of iterations needed: ⌈log₂((x2-x1)/xacc)⌉
		/// - Example: For xacc=1e-10 and initial bracket width 1.0, need ~34 iterations
		///
		/// ADVANTAGES:
		/// - Always converges (if root exists in bracket)
		/// - No derivative information needed
		/// - Predictable iteration count
		/// - Robust against pathological functions
		///
		/// DISADVANTAGES:
		/// - Slower than Newton (quadratic) or Secant (superlinear)
		/// - Requires initial bracket (not just starting point)
		/// - Linear convergence only
		///
		/// ERROR HANDLING:
		/// - Throws RootFindingError if initial values don't bracket root
		/// - Throws if function values are non-finite (NaN or Inf)
		///
		/// USE CASES:
		/// - When robustness is more important than speed
		/// - For functions with discontinuous derivatives
		/// - When you have a reliable bracket
		/// - As fallback when Newton or Secant fail
		///
		/// @param func  Real function implementing IRealFunction interface
		/// @param x1    Left endpoint of bracket (func(x1) and func(x2) must have opposite signs)
		/// @param x2    Right endpoint of bracket
		/// @param xacc  Absolute accuracy tolerance for root (typical: 1e-6 to 1e-10)
		/// @return Root estimate accurate to within xacc
		///
		/// POSTCONDITION:
		/// - |func(root)| typically much smaller than xacc
		/// - Guaranteed: root in [x1-xacc, x2+xacc]
		///
		/// Based on Numerical Recipes §9.1
		/// @see FindRootBisection(const IRealFunction&, Real, Real, const RootFindingConfig&)
		static Real FindRootBisection(const IRealFunction& func, Real x1, Real x2, Real xacc) {
			RootFindingConfig config;
			config.tolerance = xacc;
			config.max_iterations = 0;  // Use default from Defaults namespace
			
			auto result = FindRootBisection(func, x1, x2, config);
			if (!result.converged)
				throw RootFindingError(result.error_message);
			return result.root;
		}

		/// Find root using the bisection method with configurable parameters.
		///
		/// This overload provides full control over convergence parameters and
		/// returns detailed diagnostic information via RootFindingResult.
		///
		/// @param func   Real function implementing IRealFunction interface
		/// @param x1     Left endpoint of bracket
		/// @param x2     Right endpoint of bracket
		/// @param config Configuration parameters (tolerance, max_iterations, etc.)
		/// @return RootFindingResult with root value and convergence diagnostics
		///
		/// @example
		/// RootFindingConfig config;
		/// config.tolerance = 1e-12;
		/// config.max_iterations = 100;
		/// auto result = FindRootBisection(f, 0, 1, config);
		/// if (result.converged) {
		///     std::cout << "Root: " << result.root << " (iterations: " << result.iterations_used << ")\n";
		/// }
		static RootFindingResult FindRootBisection(const IRealFunction& func, Real x1, Real x2, 
												   const RootFindingConfig& config) {
			RootFindingResult result;
			Real dx, xmid, rtb;

			Real f = func(x1);
			Real fmid = func(x2);

			// Check for non-finite initial values using shared validation helper
			if (!IsFunctionValueValid(f) || !IsFunctionValueValid(fmid)) {
				result.error_message = "Non-finite function values in FindRootBisection";
				return result;
			}

			if (f * fmid >= 0.0) {
				result.error_message = "Root must be bracketed for bisection in FindRootBisection";
				return result;
			}

			if (f < 0.0) {
				dx = x2 - x1;
				rtb = x1;
			} else {
				dx = x1 - x2;
				rtb = x2;
			}

			int maxIter = (config.max_iterations > 0) ? config.max_iterations : Defaults::BisectionMaxSteps;

			for (int j = 0; j < maxIter; j++) {
				dx *= 0.5;
				xmid = rtb + dx;
				fmid = func(xmid);

				// Check for non-finite midpoint value using shared validation helper
				if (!IsFunctionValueValid(fmid)) {
					result.root = rtb;
					result.function_value = fmid;
					result.iterations_used = j + 1;
					result.achieved_tolerance = std::abs(dx);
					result.error_message = "Non-finite function value in FindRootBisection iteration";
					return result;
				}

				if (fmid <= 0.0)
					rtb = xmid;

				if (config.verbose) {
					std::cout << "Bisection iter " << j << ": x=" << rtb 
							  << ", f(x)=" << fmid << ", dx=" << dx << std::endl;
				}

				if (std::abs(dx) < config.tolerance || fmid == 0.0) {
					result.root = rtb;
					result.function_value = func(rtb);
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(dx);
					return result;
				}
			}
			
			result.root = rtb;
			result.function_value = func(rtb);
			result.iterations_used = maxIter;
			result.achieved_tolerance = std::abs(dx);
			result.error_message = "Too many bisections in FindRootBisection";
			return result;
		}

		/*********************************************************************/
		/*****           False Position (Regula Falsi) method           *****/
		/*********************************************************************/
		/// Find root using the false position (regula falsi) method.
		///
		/// False position is a root-finding algorithm that combines the reliability of
		/// bisection with faster convergence by using linear interpolation rather than
		/// simple midpoint bisection. It maintains a bracket [a,b] where f(a) and f(b)
		/// have opposite signs, then estimates the root using the secant line between
		/// these points.
		///
		/// ALGORITHM:
		/// - Start with bracket [x1, x2] where f(x1)*f(x2) < 0
		/// - Iteration loop:
		/// 1. Compute linear interpolation: x_new = (x1*f2 - x2*f1) / (f2 - f1)
		/// 2. Evaluate f(x_new)
		/// 3. Update bracket: replace x1 or x2 with x_new to maintain sign change
		/// 4. Continue until |f(x_new)| < tolerance or bracket width < xacc
		///
		/// CONVERGENCE:
		/// - Superlinear convergence (between linear and quadratic)
		/// - Order of convergence: approximately 1.618 (golden ratio φ)
		/// - Faster than bisection, slower than Newton or secant
		/// - Always converges if root exists in initial bracket
		///
		/// ADVANTAGES:
		/// - More efficient than bisection (uses function shape information)
		/// - Guaranteed convergence (maintains bracket)
		/// - No derivative needed
		/// - Robust against discontinuities
		///
		/// DISADVANTAGES:
		/// - Can be slow if one endpoint stays fixed (one-sided convergence)
		/// - Slower than Newton or secant for smooth functions
		/// - May stagnate with poor function behavior
		///
		/// IMPROVEMENTS:
		/// - Illinois method: modifies function values to prevent stagnation
		/// - Anderson-Björck: similar improvement with better performance
		///
		/// USE CASES:
		/// - When bisection is too slow but Newton is too risky
		/// - Functions with discontinuous derivatives
		/// - When you have a reliable bracket
		/// - As fallback between bisection and Newton
		///
		/// @param func  Real function implementing IRealFunction interface
		/// @param x1    Left endpoint of bracket
		/// @param x2    Right endpoint of bracket
		/// @param xacc  Absolute accuracy tolerance (typical: 1e-6 to 1e-10)
		/// @return Root estimate
		///
		/// POSTCONDITION:
		/// - |func(root)| typically < xacc
		/// - Root guaranteed in original bracket
		///
		/// Based on Numerical Recipes §9.2
		/// @see FindRootFalsePosition(const IRealFunction&, Real, Real, const RootFindingConfig&)
		static Real FindRootFalsePosition(const IRealFunction& func, Real x1, Real x2, Real xacc) {
			RootFindingConfig config;
			config.tolerance = xacc;
			config.max_iterations = 0;
			
			auto result = FindRootFalsePosition(func, x1, x2, config);
			if (!result.converged)
				throw RootFindingError(result.error_message);
			return result.root;
		}

		/// Find root using false position method with configurable parameters.
		static RootFindingResult FindRootFalsePosition(const IRealFunction& func, Real x1, Real x2,
													   const RootFindingConfig& config) {
			RootFindingResult result;
			Real f1 = func(x1);
			Real f2 = func(x2);

			if (!IsFunctionValueValid(f1) || !IsFunctionValueValid(f2)) {
				result.error_message = "Non-finite function values in FindRootFalsePosition";
				return result;
			}
			if (f1 * f2 >= 0.0) {
				result.error_message = "Root must be bracketed for false position";
				return result;
			}

			if (f1 > 0.0) {
				std::swap(x1, x2);
				std::swap(f1, f2);
			}

			int maxIter = (config.max_iterations > 0) ? config.max_iterations : Defaults::BisectionMaxSteps;
			Real rtf = x1;

			for (int j = 0; j < maxIter; j++) {
				Real dx = x2 - x1;
				rtf = x1 + dx * f1 / (f1 - f2);
				Real f = func(rtf);

				if (!IsFunctionValueValid(f)) {
					result.root = rtf;
					result.function_value = f;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite function value in FindRootFalsePosition";
					return result;
				}

				if (f < 0.0) {
					x1 = rtf;
					f1 = f;
				} else {
					x2 = rtf;
					f2 = f;
				}

				if (config.verbose) {
					std::cout << "FalsePosition iter " << j << ": x=" << rtf 
							  << ", f(x)=" << f << std::endl;
				}

				if (std::abs(f) < config.tolerance || std::abs(x2 - x1) < config.tolerance) {
					result.root = rtf;
					result.function_value = f;
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(x2 - x1);
					return result;
				}
			}
			
			result.root = rtf;
			result.function_value = func(rtf);
			result.iterations_used = maxIter;
			result.achieved_tolerance = std::abs(x2 - x1);
			result.error_message = "Maximum number of iterations exceeded in FindRootFalsePosition";
			return result;
		}

	} // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_BRACKETING_H
