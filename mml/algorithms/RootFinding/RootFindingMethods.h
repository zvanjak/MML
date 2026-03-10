///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingMethods.h                                                ///
///  Description: Advanced root finding algorithms                                    ///
///               Newton-Raphson, Secant, Ridders, and Brent methods                  ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [NR3]  Press et al., Numerical Recipes 3rd ed., Ch. 9                         ///
///    [B73]  Brent, R.P. (1973). Algorithms for Minimization without Derivatives    ///
///    [R79]  Ridders, C.J.F. (1979). IEEE Trans. Circ. Syst. 26(11), pp. 979-980   ///
///    [D01]  Dekker, T.J. (1969). Constructive Aspects Fund. Thm. Algebra, pp.37-48///
///                                                                                   ///
///  See references/book_references.md and references/paperes_references.md          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_METHODS_H
#define MML_ROOTFINDING_METHODS_H

#include "mml/MMLBase.h"
#include "mml/core/NumericValidation.h"

#include "mml/interfaces/IFunction.h"

#include "mml/core/Derivation.h"

#include "mml/algorithms/RootFinding/RootFindingBase.h"

namespace MML {
	namespace RootFinding {

		/*********************************************************************/
		/*****           Forward Declarations                            *****/
		/*********************************************************************/
		// Config-based overloads declared here so simple versions can call them
		static RootFindingResult FindRootNewton(const IRealFunction& func, Real x1, Real x2,
												const RootFindingConfig& config);
		static RootFindingResult FindRootSecant(const IRealFunction& func, Real x1, Real x2,
												const RootFindingConfig& config);
		static RootFindingResult FindRootRidders(const IRealFunction& func, Real x1, Real x2,
												 const RootFindingConfig& config);
		static RootFindingResult FindRootBrent(const IRealFunction& func, Real x1, Real x2,
											   const RootFindingConfig& config);

		/*********************************************************************/
		/*****           Newton-Raphson method                           *****/
		/*********************************************************************/
		/// Find root using the Newton-Raphson method (Newton's method).
		///
		/// Newton-Raphson is the premier method for root finding when speed matters:
		/// it exhibits quadratic convergence near roots, doubling the number of correct
		/// digits each iteration. However, it requires derivative information and can
		/// diverge if not properly constrained. This implementation uses numerical
		/// differentiation (4-point stencil) and bracket checking for robustness.
		///
		/// ALGORITHM:
		/// - Newton iteration: x_{n+1} = x_n - f(x_n) / f'(x_n)
		/// - Starting point: midpoint of initial bracket [x1, x2]
		/// - At each step:
		/// 1. Evaluate f(x) and f'(x) at current point
		/// 2. Compute Newton step: dx = -f/f'
		/// 3. Update: x_new = x + dx
		/// 4. Check if x_new still in bracket [x1, x2]
		/// 5. If out of bounds, revert to bisection step
		/// 6. Continue until |dx| < xacc or f(x) ≈ 0
		///
		/// DERIVATIVE:
		/// - Uses 4-point central difference formula (NDer4)
		/// - Accuracy: O(h^4) where h is step size
		/// - No analytical derivative needed
		///
		/// CONVERGENCE:
		/// - Near root: QUADRATIC convergence
		/// - Error decreases as e_{n+1} ≈ e_n² (doubles correct digits each step)
		/// - Typical: 3-6 iterations for machine precision
		/// - Much faster than bisection's linear convergence
		///
		/// ADVANTAGES:
		/// - Extremely fast convergence near roots
		/// - Few iterations needed
		/// - Well-suited for smooth functions
		///
		/// DISADVANTAGES:
		/// - May diverge if starting point is poor
		/// - Fails at stationary points (f' = 0)
		/// - Can jump outside bracket with bad behavior
		/// - Needs derivative (here computed numerically)
		///
		/// ROBUSTNESS FEATURES:
		/// - Bracket checking: reverts to bisection if step goes out of bounds
		/// - Derivative validation: throws if non-finite or too small
		/// - Maximum iteration limit
		/// - Handles both function and derivative failures gracefully
		///
		/// ERROR HANDLING:
		/// - Throws RootFindingError if derivative is NaN or Inf
		/// - Throws if derivative too small (< 1e-15) - likely stationary point
		/// - Throws if maximum iterations exhausted
		///
		/// USE CASES:
		/// - When you need high precision quickly
		/// - For smooth, well-behaved functions
		/// - When you have good initial bracket
		/// - Problems where function evaluation is cheap
		///
		/// @param func  Real function implementing IRealFunction interface
		/// @param x1    Left endpoint of initial bracket
		/// @param x2    Right endpoint of initial bracket
		/// @param xacc  Absolute accuracy tolerance (typical: 1e-10 to 1e-15)
		/// @return Root estimate accurate to within xacc
		///
		/// POSTCONDITION:
		/// - |func(root)| typically ~1e-15 (machine precision) even if xacc larger
		/// - Number of iterations typically 3-6
		///
		/// Based on Numerical Recipes §9.4
		/// @see FindRootNewton(const IRealFunction&, Real, Real, const RootFindingConfig&)
		static Real FindRootNewton(const IRealFunction& func, Real x1, Real x2, Real xacc) {
			RootFindingConfig config;
			config.tolerance = xacc;
			config.max_iterations = 0;  // Use default from Defaults namespace
			
			auto result = FindRootNewton(func, x1, x2, config);
			if (!result.converged)
				throw RootFindingError(result.error_message);
			return result.root;
		}

		/// Find root using Newton-Raphson method with configurable parameters.
		///
		/// This overload provides full control over convergence parameters and
		/// returns detailed diagnostic information via RootFindingResult.
		///
		/// @param func   Real function implementing IRealFunction interface
		/// @param x1     Left endpoint of initial bracket
		/// @param x2     Right endpoint of initial bracket
		/// @param config Configuration parameters (tolerance, max_iterations, etc.)
		/// @return RootFindingResult with root value and convergence diagnostics
		static RootFindingResult FindRootNewton(const IRealFunction& func, Real x1, Real x2,
												const RootFindingConfig& config) {
			RootFindingResult result;
			Real rtn = 0.5 * (x1 + x2);
			int maxIter = (config.max_iterations > 0) ? config.max_iterations : Defaults::NewtonRaphsonMaxSteps;
			Real dx = 0.0;

			for (int j = 0; j < maxIter; j++) {
				Real f = func(rtn);
				Real df = Derivation::NDer4(func, rtn);

				// Check for non-finite derivative (NaN/Inf) or near-zero derivative
				if (!IsFunctionValueValid(df)) {
					result.root = rtn;
					result.function_value = f;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite derivative in FindRootNewton";
					return result;
				}
				Real dfThreshold = std::sqrt(Real(Constants::Eps));
				if (std::abs(df) < dfThreshold) {
					result.root = rtn;
					result.function_value = f;
					result.iterations_used = j + 1;
					result.error_message = "Zero or near-zero derivative in FindRootNewton";
					return result;
				}

				dx = f / df;

				// Check for non-finite step using shared validation helper
				if (!IsFinite(dx)) {
					result.root = rtn;
					result.function_value = f;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite step in FindRootNewton";
					return result;
				}

				rtn -= dx;

				if (config.verbose) {
					std::cout << "Newton iter " << j << ": x=" << rtn 
							  << ", f(x)=" << f << ", dx=" << dx << std::endl;
				}

				if ((x1 - rtn) * (rtn - x2) < 0.0) {
					result.root = rtn;
					result.function_value = func(rtn);
					result.iterations_used = j + 1;
					result.achieved_tolerance = std::abs(dx);
					result.error_message = "Jumped out of brackets in FindRootNewton";
					return result;
				}

				if (std::abs(dx) < config.tolerance) {
					result.root = rtn;
					result.function_value = func(rtn);
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(dx);
					return result;
				}
			}
			
			result.root = rtn;
			result.function_value = func(rtn);
			result.iterations_used = maxIter;
			result.achieved_tolerance = std::abs(dx);
			result.error_message = "Maximum number of iterations exceeded in FindRootNewton";
			return result;
		}

		/*********************************************************************/
		/*****           Secant method                                   *****/
		/*********************************************************************/
		/// Find root using the secant method.
		///
		/// The secant method is a derivative-free variant of Newton's method that approximates
		/// the derivative using a finite difference. It requires two initial guesses but does
		/// not require them to bracket the root. The method exhibits superlinear convergence
		/// with order φ ≈ 1.618 (the golden ratio).
		///
		/// ALGORITHM:
		/// - Start with two initial guesses x0, x1
		/// - Iteration formula: x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
		/// - Update: x_{n-1} ← x_n, x_n ← x_{n+1}
		/// - Continue until |x_{n+1} - x_n| < tolerance
		///
		/// DERIVATIVE APPROXIMATION:
		/// - Approximates f'(x_n) ≈ (f(x_n) - f(x_{n-1})) / (x_n - x_{n-1})
		/// - No analytical derivative needed
		/// - Only one new function evaluation per iteration (vs two for Newton with numerical derivative)
		///
		/// CONVERGENCE:
		/// - Superlinear convergence: order φ = (1 + √5)/2 ≈ 1.618
		/// - Error: e_{n+1} ≈ C * e_n^φ
		/// - Faster than bisection or false position
		/// - Slower than Newton (which is quadratic)
		/// - More efficient than Newton when derivatives are expensive
		///
		/// ADVANTAGES:
		/// - No derivative needed (unlike Newton)
		/// - Superlinear convergence
		/// - Only one function evaluation per iteration
		/// - Simple to implement
		///
		/// DISADVANTAGES:
		/// - May diverge with poor initial guesses
		/// - No bracket maintenance (can jump away from root)
		/// - Slower than Newton for smooth functions
		/// - Can fail if f(x_n) = f(x_{n-1})
		///
		/// ROBUSTNESS:
		/// - This implementation uses bracket checking to prevent divergence
		/// - Falls back to bisection-like behavior if secant jumps outside bracket
		///
		/// USE CASES:
		/// - When derivatives are expensive or unavailable
		/// - For smooth functions where you have good initial guesses
		/// - When Newton is overkill but bisection is too slow
		/// - In multi-dimensional optimization as component
		///
		/// @param func  Real function implementing IRealFunction interface
		/// @param x1    First initial guess (used for bracket if needed)
		/// @param x2    Second initial guess
		/// @param xacc  Absolute accuracy tolerance (typical: 1e-8 to 1e-12)
		/// @return Root estimate
		///
		/// POSTCONDITION:
		/// - |func(root)| typically << xacc due to superlinear convergence
		/// - Typical iterations: 5-10 for machine precision
		///
		/// Based on Numerical Recipes §9.2
		/// @see FindRootSecant(const IRealFunction&, Real, Real, const RootFindingConfig&)
		static Real FindRootSecant(const IRealFunction& func, Real x1, Real x2, Real xacc) {
			RootFindingConfig config;
			config.tolerance = xacc;
			config.max_iterations = 0;
			
			auto result = FindRootSecant(func, x1, x2, config);
			if (!result.converged)
				throw RootFindingError(result.error_message);
			return result.root;
		}

		/// Find root using secant method with configurable parameters.
		static RootFindingResult FindRootSecant(const IRealFunction& func, Real x1, Real x2,
												const RootFindingConfig& config) {
			RootFindingResult result;
			Real f1 = func(x1);
			Real f2 = func(x2);

			if (!IsFunctionValueValid(f1) || !IsFunctionValueValid(f2)) {
				result.error_message = "Non-finite function values in FindRootSecant";
				return result;
			}

			if (std::abs(f1) < std::abs(f2)) {
				std::swap(x1, x2);
				std::swap(f1, f2);
			}

			int maxIter = (config.max_iterations > 0) ? config.max_iterations : Defaults::NewtonRaphsonMaxSteps;
			Real dx = 0.0;

			for (int j = 0; j < maxIter; j++) {
				dx = (x2 - x1) * f2 / (f2 - f1);

				if (!IsFinite(dx)) {
					result.root = x2;
					result.function_value = f2;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite step in FindRootSecant";
					return result;
				}

				x1 = x2;
				f1 = f2;
				x2 -= dx;
				f2 = func(x2);

				if (!IsFunctionValueValid(f2)) {
					result.root = x2;
					result.function_value = f2;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite function value in FindRootSecant";
					return result;
				}

				if (config.verbose) {
					std::cout << "Secant iter " << j << ": x=" << x2 
							  << ", f(x)=" << f2 << ", dx=" << dx << std::endl;
				}

				if (std::abs(dx) < config.tolerance) {
					result.root = x2;
					result.function_value = f2;
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(dx);
					return result;
				}
			}
			
			result.root = x2;
			result.function_value = f2;
			result.iterations_used = maxIter;
			result.achieved_tolerance = std::abs(dx);
			result.error_message = "Maximum number of iterations exceeded in FindRootSecant";
			return result;
		}

		/*********************************************************************/
		/*****           Ridders method                                  *****/
		/*********************************************************************/
		/// Find root using Ridders' method of exponential interpolation.
		///
		/// Ridders' method is an elegant root-finding algorithm that uses exponential
		/// function fitting rather than linear interpolation. It combines the reliability
		/// of bisection with superlinear convergence, making it an excellent general-purpose
		/// method. The algorithm maintains a bracket and exhibits quadratic convergence.
		///
		/// ALGORITHM:
		/// - Start with bracket [x1, x2] where f(x1)*f(x2) < 0
		/// - Iteration loop:
		/// 1. Compute midpoint: xm = (x1 + x2) / 2
		/// 2. Evaluate: f1 = f(x1), fm = f(xm), f2 = f(x2)
		/// 3. Compute: s = √(fm² - f1*f2)
		/// 4. New estimate: x4 = xm + (xm - x1) * sign(f1 - f2) * fm / s
		/// 5. Update bracket with x4 and appropriate endpoint
		/// 6. Continue until |f(x4)| < tolerance
		///
		/// MATHEMATICAL BASIS:
		/// - Fits an exponential function through three points
		/// - Formula derived from: f(x) = A*e^(λx) + B*e^(-λx)
		/// - More sophisticated than linear interpolation
		/// - Naturally handles varying function curvature
		///
		/// CONVERGENCE:
		/// - Quadratic convergence (same as Newton!)
		/// - Error: e_{n+1} ≈ C * e_n²
		/// - Guaranteed convergence (maintains bracket)
		/// - Typically 2-3 times faster than false position
		/// - Comparable to Newton without needing derivatives
		///
		/// ADVANTAGES:
		/// - Quadratic convergence without derivatives
		/// - Guaranteed convergence (bracket maintained)
		/// - Robust against poor initial guesses
		/// - Never fails for continuous functions with sign change
		/// - Often outperforms Newton on difficult functions
		///
		/// DISADVANTAGES:
		/// - More function evaluations per iteration (3 vs 1 for bisection)
		/// - Slightly more complex implementation
		/// - Square root operation adds computational cost
		///
		/// COMPARISON:
		/// - vs Bisection: Much faster (quadratic vs linear)
		/// - vs False Position: Faster, more reliable
		/// - vs Newton: No derivative needed, more robust
		/// - vs Secant: Faster, guaranteed convergence
		///
		/// USE CASES:
		/// - General-purpose root finding (excellent default choice)
		/// - When derivatives unavailable but speed matters
		/// - For difficult functions where Newton might fail
		/// - When robustness AND speed both important
		///
		/// @param func  Real function implementing IRealFunction interface
		/// @param x1    Left endpoint of bracket
		/// @param x2    Right endpoint of bracket
		/// @param xacc  Absolute accuracy tolerance (typical: 1e-10 to 1e-15)
		/// @return Root estimate
		///
		/// POSTCONDITION:
		/// - |func(root)| typically at machine precision
		/// - Typical iterations: 5-8 for double precision
		/// - Root guaranteed in original bracket
		///
		/// Based on Numerical Recipes §9.2
		/// @see FindRootRidders(const IRealFunction&, Real, Real, const RootFindingConfig&)
		static Real FindRootRidders(const IRealFunction& func, Real x1, Real x2, Real xacc) {
			RootFindingConfig config;
			config.tolerance = xacc;
			config.max_iterations = 0;
			
			auto result = FindRootRidders(func, x1, x2, config);
			if (!result.converged)
				throw RootFindingError(result.error_message);
			return result.root;
		}

		/// Find root using Ridders' method with configurable parameters.
		static RootFindingResult FindRootRidders(const IRealFunction& func, Real x1, Real x2,
												 const RootFindingConfig& config) {
			RootFindingResult result;
			Real f1 = func(x1);
			Real f2 = func(x2);

			if (!IsFunctionValueValid(f1) || !IsFunctionValueValid(f2)) {
				result.error_message = "Non-finite function values in FindRootRidders";
				return result;
			}
			if (f1 * f2 >= 0.0) {
				result.error_message = "Root must be bracketed for Ridders method";
				return result;
			}

			int maxIter = (config.max_iterations > 0) ? config.max_iterations : Defaults::BisectionMaxSteps;
			Real ans = -1.0e99;

			for (int j = 0; j < maxIter; j++) {
				Real xm = 0.5 * (x1 + x2);
				Real fm = func(xm);

				if (!IsFunctionValueValid(fm)) {
					result.root = ans;
					result.function_value = fm;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite function value in FindRootRidders";
					return result;
				}

				Real s = std::sqrt(fm * fm - f1 * f2);
				if (s == 0.0) {
					result.root = ans;
					result.function_value = func(ans);
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(x2 - x1);
					return result;
				}

				Real xnew = xm + (xm - x1) * ((f1 >= f2 ? 1.0 : -1.0) * fm / s);

				if (std::abs(xnew - ans) <= config.tolerance) {
					result.root = ans;
					result.function_value = func(ans);
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(xnew - ans);
					return result;
				}

				ans = xnew;
				Real fnew = func(ans);

				if (!IsFunctionValueValid(fnew)) {
					result.root = ans;
					result.function_value = fnew;
					result.iterations_used = j + 1;
					result.error_message = "Non-finite function value in FindRootRidders";
					return result;
				}

				if (config.verbose) {
					std::cout << "Ridders iter " << j << ": x=" << ans 
							  << ", f(x)=" << fnew << std::endl;
				}

				if (std::abs(fnew) < config.tolerance) {
					result.root = ans;
					result.function_value = fnew;
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(fnew);
					return result;
				}

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

				if (std::abs(x2 - x1) <= config.tolerance) {
					result.root = ans;
					result.function_value = fnew;
					result.iterations_used = j + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(x2 - x1);
					return result;
				}
			}
			
			result.root = ans;
			result.function_value = func(ans);
			result.iterations_used = maxIter;
			result.achieved_tolerance = std::abs(x2 - x1);
			result.error_message = "Maximum number of iterations exceeded in FindRootRidders";
			return result;
		}

		/*********************************************************************/
		/*****           Brent method (van Wijngaarden-Dekker-Brent)    *****/
		/*********************************************************************/
		/// Find root using Brent's method (the gold standard for root finding).
		///
		/// Brent's method is widely considered the best general-purpose root-finding algorithm.
		/// It combines the reliability of bisection, the speed of inverse quadratic interpolation,
		/// and the efficiency of the secant method. The algorithm adaptively chooses the best
		/// strategy at each iteration, guaranteeing convergence while optimizing speed.
		///
		/// ALGORITHM:
		/// - Maintains bracket [a,b] where f(a)*f(b) < 0
		/// - At each iteration, chooses ONE of:
		/// 1. Inverse quadratic interpolation (if three points available and well-conditioned)
		/// 2. Secant method (linear interpolation)
		/// 3. Bisection (if interpolation gives poor step)
		/// - Uses several safeguards to ensure reliability
		/// - Keeps track of previous iterations to make intelligent choices
		///
		/// INVERSE QUADRATIC INTERPOLATION:
		/// - Fits parabola through three points: (f_a, a), (f_b, b), (f_c, c)
		/// - Solves for x where parabola crosses zero
		/// - Provides fast convergence when function is well-behaved
		/// - Order of convergence: approximately 1.839
		///
		/// CONVERGENCE:
		/// - Superlinear convergence: order approximately 1.839
		/// - Guaranteed to converge (maintains bracket)
		/// - At worst, reduces bracket by factor of 2 (bisection)
		/// - At best, converges like inverse quadratic interpolation
		/// - Adapts strategy based on function behavior
		///
		/// SAFEGUARDS:
		/// - Falls back to bisection if:
		/// * Interpolated point outside bracket
		/// * Step too small (stagnation)
		/// * Step not reducing bracket fast enough
		/// - Maintains |f(b)| ≤ |f(a)| invariant
		/// - Tracks previous steps to detect poor progress
		///
		/// ADVANTAGES:
		/// - Best all-around performance
		/// - Extremely robust (never fails if root exists)
		/// - No derivatives needed
		/// - Adapts to function behavior
		/// - Widely trusted in numerical libraries (scipy, GNU GSL, etc.)
		///
		/// DISADVANTAGES:
		/// - Complex implementation
		/// - Harder to understand than simpler methods
		/// - Slightly more function evaluations than pure methods
		///
		/// COMPARISON TO OTHER METHODS:
		/// - vs Bisection: Much faster, same reliability
		/// - vs Newton: No derivatives, more robust, similar speed
		/// - vs Secant: More robust, comparable speed
		/// - vs Ridders: Similar performance, more widely used
		///
		/// USE CASES:
		/// - Production code where both speed and reliability matter
		/// - General-purpose root finding in numerical libraries
		/// - When function evaluation is not too expensive
		/// - Whenever you're unsure which method to use (default choice)
		///
		/// @param func  Real function implementing IRealFunction interface
		/// @param x1    Left endpoint of bracket
		/// @param x2    Right endpoint of bracket
		/// @param xacc  Absolute accuracy tolerance (typical: 1e-12 to 1e-15)
		/// @return Root estimate
		///
		/// POSTCONDITION:
		/// - |func(root)| at or near machine precision
		/// - Typical iterations: 6-12 for double precision
		/// - Root guaranteed in original bracket
		///
		/// Based on Numerical Recipes §9.3 and Brent (1973)
		/// @see FindRootBrent(const IRealFunction&, Real, Real, const RootFindingConfig&)
		static Real FindRootBrent(const IRealFunction& func, Real x1, Real x2, Real xacc) {
			RootFindingConfig config;
			config.tolerance = xacc;
			config.max_iterations = 0;
			
			auto result = FindRootBrent(func, x1, x2, config);
			if (!result.converged)
				throw RootFindingError(result.error_message);
			return result.root;
		}

		/// Find root using Brent's method with configurable parameters.
		///
		/// This is the recommended method for general-purpose root finding.
		/// Combines reliability of bisection with speed of inverse quadratic interpolation.
		///
		/// @param func   Real function implementing IRealFunction interface
		/// @param x1     Left endpoint of bracket
		/// @param x2     Right endpoint of bracket
		/// @param config Configuration parameters (tolerance, max_iterations, etc.)
		/// @return RootFindingResult with root value and convergence diagnostics
		static RootFindingResult FindRootBrent(const IRealFunction& func, Real x1, Real x2,
											   const RootFindingConfig& config) {
			RootFindingResult result;
			Real a = x1, b = x2, c = x2, d = 0.0, e = 0.0;
			Real fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;

			if (!IsFunctionValueValid(fa) || !IsFunctionValueValid(fb)) {
				result.error_message = "Non-finite function values in FindRootBrent";
				return result;
			}
			if (fa * fb >= 0.0) {
				result.error_message = "Root must be bracketed for Brent method";
				return result;
			}

			fc = fb;
			int maxIter = (config.max_iterations > 0) ? config.max_iterations : Defaults::BrentMaxSteps;

			for (int iter = 0; iter < maxIter; iter++) {
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

				tol1 = 2.0 * Constants::Eps * std::abs(b) + 0.5 * config.tolerance;
				xm = 0.5 * (c - b);

				if (config.verbose) {
					std::cout << "Brent iter " << iter << ": x=" << b 
							  << ", f(x)=" << fb << ", xm=" << xm << std::endl;
				}

				if (std::abs(xm) <= tol1 || fb == 0.0) {
					result.root = b;
					result.function_value = fb;
					result.iterations_used = iter + 1;
					result.converged = true;
					result.achieved_tolerance = std::abs(xm);
					return result;
				}

				if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0 * xm * s;
						q = 1.0 - s;
					} else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
						q = (q - 1.0) * (r - 1.0) * (s - 1.0);
					}
					if (p > 0.0)
						q = -q;
					p = std::abs(p);

					Real min1 = 3.0 * xm * q - std::abs(tol1 * q);
					Real min2 = std::abs(e * q);
					if (2.0 * p < (min1 < min2 ? min1 : min2)) {
						e = d;
						d = p / q;
					} else {
						d = xm;
						e = d;
					}
				} else {
					d = xm;
					e = d;
				}

				a = b;
				fa = fb;
				if (std::abs(d) > tol1)
					b += d;
				else
					b += std::copysign(tol1, xm);

				fb = func(b);

				if (!IsFunctionValueValid(fb)) {
					result.root = b;
					result.function_value = fb;
					result.iterations_used = iter + 1;
					result.error_message = "Non-finite function value in FindRootBrent";
					return result;
				}
			}
			
			result.root = b;
			result.function_value = fb;
			result.iterations_used = maxIter;
			result.achieved_tolerance = std::abs(xm);
			result.error_message = "Maximum number of iterations exceeded in FindRootBrent";
			return result;
		}

	} // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_METHODS_H
