///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Multidim/QuasiNewton.h                                                ///
///  Description: Conjugate gradient, BFGS, and L-BFGS minimization                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_QUASI_NEWTON_H
#define MML_OPTIMIZATION_QUASI_NEWTON_H

#include <mml/algorithms/Optimization/Multidim/LineSearch.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace MML::Optimization {

	///            Conjugate Gradient Method (Fletcher-Reeves-Polak-Ribiere)///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Conjugate gradient method for minimization with derivatives
     * 
     * The conjugate gradient method is one of the most effective methods for
     * minimizing smooth functions when gradients are available. It generates
     * search directions that are conjugate with respect to the Hessian.
     * 
     * Reference: Numerical Recipes Chapter 10.8
     * 
     * Two variants are implemented:
     * - Fletcher-Reeves: gamma = (g_new · g_new) / (g_old · g_old)
     * - Polak-Ribiere:   gamma = (g_new · (g_new - g_old)) / (g_old · g_old)
     * 
     * Polak-Ribiere is generally preferred as it automatically resets to
     * steepest descent when progress stalls.
     * 
     * The algorithm:
     * 1. Compute gradient at starting point
     * 2. Set initial search direction to negative gradient (steepest descent)
     * 3. Minimize along search direction
     * 4. Compute new gradient
     * 5. Update search direction using conjugate gradient formula
     * 6. Repeat until convergence
     */
	class ConjugateGradient {
	public:
		enum class Method {
			FletcherReeves, // Original CG formula
			PolakRibiere	// Generally preferred variant
		};

	private:
		Real _ftol;		// Function tolerance
		Real _gtol;		// Gradient tolerance
		int _maxIter;	// Maximum iterations
		int _iter;		// Iteration counter
		Real _fret;		// Current function value
		Method _method; // CG variant

	public:
		ConjugateGradient(Real ftol = PrecisionValues<Real>::OptimizationTolerance, Real gtol = PrecisionValues<Real>::OptimizationGradientTolerance, int maxIter = 200, Method method = Method::PolakRibiere)
			: _ftol(ftol)
			, _gtol(gtol)
			, _maxIter(maxIter)
			, _iter(0)
			, _fret(0.0)
			, _method(method) {}

		Real getFtol() const { return _ftol; }
		void setFtol(Real ftol) { _ftol = ftol; }

		Real getGtol() const { return _gtol; }
		void setGtol(Real gtol) { _gtol = gtol; }

		int getMaxIter() const { return _maxIter; }
		void setMaxIter(int maxIter) { _maxIter = maxIter; }

		Method getMethod() const { return _method; }
		void setMethod(Method method) { _method = method; }

		int getIterations() const { return _iter; }
		Real getCurrentFValue() const { return _fret; }

		/**
         * @brief Minimize a differentiable function using conjugate gradient
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimOptimizationResult Minimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start) {
			// Validate inputs
			ValidateVectorFinite<N>(start, "ConjugateGradient::Minimize starting point");
			ValidateMultidimTolerance(_ftol, "ConjugateGradient ftol");
			ValidateMultidimTolerance(_gtol, "ConjugateGradient gtol");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("ConjugateGradient: maxIter must be positive");
			
			const Real EPS = Precision::Eps;   // NR dfpmin: machine epsilon

			VectorN<Real, N> p = start;
			VectorN<Real, N> g, h, xi;

			// Evaluate function and gradient at starting point
			Real fp = func(p);
			ValidateMultidimFunctionValue(fp, "ConjugateGradient initial evaluation");
			func.Gradient(p, xi);
			ValidateVectorFinite<N>(xi, "ConjugateGradient initial gradient");

			// Initialize: g = -gradient, h = xi = g (steepest descent)
			for (int j = 0; j < N; ++j) {
				g[j] = -xi[j];
				xi[j] = h[j] = g[j];
			}

			for (_iter = 0; _iter < _maxIter; ++_iter) {
				// Line minimization along xi
				_fret = LineMinimizer::Minimize(func, p, xi, _ftol);

				// Check function convergence
				if (2.0 * std::abs(_fret - fp) <= _ftol * (std::abs(_fret) + std::abs(fp) + EPS)) {
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					return MultidimOptimizationResult(result, _fret, _iter + 1, true);
				}

				fp = _fret;

				// Compute new gradient
				func.Gradient(p, xi);

				// Check gradient convergence
				Real test = 0.0;
				Real den = std::max(std::abs(fp), static_cast<MML::Real>(1.0));
				for (int j = 0; j < N; ++j) {
					Real temp = std::abs(xi[j]) * std::max(std::abs(p[j]), static_cast<MML::Real>(1.0)) / den;
					if (temp > test)
						test = temp;
				}
				if (test < _gtol) {
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					return MultidimOptimizationResult(result, _fret, _iter + 1, true);
				}

				// Compute gamma (CG update coefficient)
				Real gg = 0.0, dgg = 0.0;
				for (int j = 0; j < N; ++j) {
					gg += g[j] * g[j];

					if (_method == Method::FletcherReeves) {
						dgg += xi[j] * xi[j]; // Fletcher-Reeves
					} else {
						dgg += (xi[j] + g[j]) * xi[j]; // Polak-Ribiere
					}
				}

				if (gg == 0.0) {
					// Gradient is zero - we're done
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					return MultidimOptimizationResult(result, _fret, _iter + 1, true);
				}

				Real gam = dgg / gg;

				// Update search direction
				for (int j = 0; j < N; ++j) {
					g[j] = -xi[j];
					xi[j] = h[j] = g[j] + gam * h[j];
				}
			}

			// Max iterations reached
			Vector<Real> result(N);
			for (int j = 0; j < N; ++j)
				result[j] = p[j];
			return MultidimOptimizationResult(result, _fret, _iter, false);
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///               Shared quasi-Newton line-search helpers               ///
	///////////////////////////////////////////////////////////////////////////
	namespace detail {

		/// Outcome of an inexact line search. `accepted == false` means the point was left
		/// unchanged and the caller must not interpret the unchanged objective as convergence.
		/// `stepBelowResolution` marks the NR "TOLX" case: every finite trial along a genuine
		/// descent direction failed to decrease f down to a step that no longer changes x by a
		/// relative ZeroTolerance; this is x-convergence evidence, unlike a nonfinite valley.
		struct QuasiNewtonLineSearchResult {
			Real f = std::numeric_limits<Real>::quiet_NaN();
			Real alpha = 0;
			bool accepted = false;
			bool stepBelowResolution = false;
			const char* reason = "";
		};

		/// NR dfpmin gradient test: max_j |g_j| * max(|x_j|,1) / max(|f|,1) < gtol
		template<int N>
		inline bool QuasiNewtonGradientConverged(const VectorN<Real, N>& g, const VectorN<Real, N>& x, Real f, Real gtol) {
			Real test = 0.0;
			Real den = std::max(std::abs(f), static_cast<MML::Real>(1.0));
			for (int j = 0; j < N; ++j) {
				Real temp = std::abs(g[j]) * std::max(std::abs(x[j]), static_cast<MML::Real>(1.0)) / den;
				if (temp > test)
					test = temp;
			}
			return test < gtol;
		}

		/// Scale a search direction down to the trust length stpmax; returns the directional
		/// derivative, resetting to steepest descent if the direction is not a descent direction.
		template<int N>
		inline Real PrepareDescentDirection(const VectorN<Real, N>& g, VectorN<Real, N>& d, Real stpmax) {
			Real dnorm = 0.0;
			for (int j = 0; j < N; ++j)
				dnorm += d[j] * d[j];
			dnorm = std::sqrt(dnorm);
			if (dnorm > stpmax) {
				Real scale = stpmax / dnorm;
				for (int j = 0; j < N; ++j)
					d[j] *= scale;
			}

			Real slope = 0.0;
			for (int j = 0; j < N; ++j)
				slope += g[j] * d[j];
			if (slope >= 0.0) {
				for (int j = 0; j < N; ++j)
					d[j] = -g[j];
				slope = 0.0;
				for (int j = 0; j < N; ++j)
					slope += g[j] * d[j];
			}
			return slope;
		}

		/// Smallest step that still changes x by a relative ZeroTolerance along d.
		template<int N>
		inline Real MinimumMeaningfulStep(const VectorN<Real, N>& x, const VectorN<Real, N>& d) {
			Real test = 0.0;
			for (int j = 0; j < N; ++j) {
				Real temp = std::abs(d[j]) / std::max(std::abs(x[j]), static_cast<MML::Real>(1.0));
				if (temp > test)
					test = temp;
			}
			return test > 0 ? Precision::ZeroTolerance / test : std::numeric_limits<Real>::infinity();
		}

		/**
		 * @brief Bounded Armijo backtracking line search (NR lnsrch with safeguards).
		 *
		 * Differences from the textbook version: nonfinite trial values contract the step
		 * before any interpolation, the loop is bounded, and falling below the minimum
		 * meaningful step is reported as a failure instead of accepting a no-op move.
		 */
		template<int N>
		inline QuasiNewtonLineSearchResult BacktrackingLineSearch(const IDifferentiableScalarFunction<N>& func,
			VectorN<Real, N>& x, Real f, const VectorN<Real, N>& g, VectorN<Real, N>& d, Real stpmax) {
			const Real ALF = Precision::ArmijoSufficientDecrease;
			const int maxBacktracks = 60;

			QuasiNewtonLineSearchResult result;
			Real slope = PrepareDescentDirection<N>(g, d, stpmax);
			if (slope >= 0.0) {
				result.reason = "zero gradient: no descent direction";
				return result;
			}
			const Real alamin = MinimumMeaningfulStep<N>(x, d);

			Real alam = 1.0;
			Real alam2 = 0.0, f2 = 0.0;
			bool havePrevious = false;
			bool sawNonfinite = false;

			for (int iter = 0; iter < maxBacktracks; ++iter) {
				if (alam < alamin) {
					result.stepBelowResolution = !sawNonfinite;
					result.reason = sawNonfinite
						? "step fell below the minimum meaningful step while contracting away from nonfinite trials"
						: "step fell below the minimum meaningful step without sufficient decrease";
					return result;
				}

				VectorN<Real, N> xnew;
				for (int j = 0; j < N; ++j)
					xnew[j] = x[j] + alam * d[j];
				Real fnew = func(xnew);

				if (!std::isfinite(fnew)) {
					// Invalid trial: contract geometrically, do not feed it to the interpolant
					sawNonfinite = true;
					alam *= static_cast<MML::Real>(0.5);
					havePrevious = false;
					continue;
				}

				if (fnew <= f + ALF * alam * slope) {
					x = xnew;
					result.f = fnew;
					result.alpha = alam;
					result.accepted = true;
					return result;
				}

				Real tmplam;
				if (!havePrevious) {
					tmplam = -slope / (static_cast<MML::Real>(2.0) * (fnew - f - alam * slope));
				} else {
					Real rhs1 = fnew - f - alam * slope;
					Real rhs2 = f2 - f - alam2 * slope;
					Real a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
					Real b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
					if (a == 0.0) {
						tmplam = -slope / (static_cast<MML::Real>(2.0) * b);
					} else {
						Real disc = b * b - static_cast<MML::Real>(3.0) * a * slope;
						if (disc < 0.0)
							tmplam = static_cast<MML::Real>(0.5) * alam;
						else if (b <= 0.0)
							tmplam = (-b + std::sqrt(disc)) / (static_cast<MML::Real>(3.0) * a);
						else
							tmplam = -slope / (b + std::sqrt(disc));
					}
					if (tmplam > static_cast<MML::Real>(0.5) * alam)
						tmplam = static_cast<MML::Real>(0.5) * alam;
				}
				if (!std::isfinite(tmplam))
					tmplam = static_cast<MML::Real>(0.5) * alam;

				alam2 = alam;
				f2 = fnew;
				havePrevious = true;
				alam = std::max(tmplam, static_cast<MML::Real>(0.1) * alam);
			}

			result.reason = "backtracking budget exhausted without sufficient decrease";
			return result;
		}

	} // namespace detail

	///////////////////////////////////////////////////////////////////////////
	///               BFGS Quasi-Newton Method                              ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief BFGS (Broyden-Fletcher-Goldfarb-Shanno) quasi-Newton method
     * 
     * BFGS is one of the most effective methods for smooth unconstrained
     * optimization. It builds an approximation to the inverse Hessian matrix
     * using only gradient information.
     * 
     * The algorithm:
     * 1. Start with identity approximation to inverse Hessian
     * 2. Compute search direction p = -H * gradient
     * 3. Line search to find step size
     * 4. Update inverse Hessian approximation using BFGS formula
     * 5. Repeat until convergence
     * 
     * Advantages:
     * - Superlinear convergence rate
     * - Only requires gradient (not Hessian)
     * - Self-correcting: recovers from poor initial H estimate
     * 
     * Disadvantages:
     * - O(N²) storage for inverse Hessian approximation
     * - May be expensive for very large N
     */
	class BFGS {
	private:
		Real _ftol;	  // Function tolerance
		Real _gtol;	  // Gradient tolerance
		int _maxIter; // Maximum iterations
		int _iter;	  // Iteration counter
		Real _fret;	  // Current function value

	public:
		BFGS(Real ftol = PrecisionValues<Real>::OptimizationTolerance, Real gtol = PrecisionValues<Real>::OptimizationGradientTolerance, int maxIter = 200)
			: _ftol(ftol)
			, _gtol(gtol)
			, _maxIter(maxIter)
			, _iter(0)
			, _fret(0.0) {}

		Real getFtol() const { return _ftol; }
		void setFtol(Real ftol) { _ftol = ftol; }

		Real getGtol() const { return _gtol; }
		void setGtol(Real gtol) { _gtol = gtol; }

		int getMaxIter() const { return _maxIter; }
		void setMaxIter(int maxIter) { _maxIter = maxIter; }

		int getIterations() const { return _iter; }
		Real getCurrentFValue() const { return _fret; }

		/**
         * @brief Minimize using BFGS quasi-Newton method
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimOptimizationResult Minimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start) {
			// Validate inputs
			ValidateVectorFinite<N>(start, "BFGS::Minimize starting point");
			ValidateMultidimTolerance(_ftol, "BFGS ftol");
			ValidateMultidimTolerance(_gtol, "BFGS gtol");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("BFGS: maxIter must be positive");
			
			const Real EPS = Precision::Eps;   // NR dfpmin: machine epsilon
			const Real STPMX = 100.0; // Maximum step size

			VectorN<Real, N> p = start;
			VectorN<Real, N> g, dg, hdg, pnew, xi;
			Matrix<Real> hessin(N, N, 0.0);

			// Initialize inverse Hessian to identity
			for (int i = 0; i < N; ++i)
				hessin(i, i) = 1.0;

			// Evaluate function and gradient
			_fret = func(p);
			ValidateMultidimFunctionValue(_fret, "BFGS initial evaluation");
			func.Gradient(p, g);
			ValidateVectorFinite<N>(g, "BFGS initial gradient");

			// Initial search direction is steepest descent
			for (int j = 0; j < N; ++j)
				xi[j] = -g[j];

			// Compute maximum step size
			Real sum = 0.0;
			for (int j = 0; j < N; ++j)
				sum += p[j] * p[j];
			Real stpmax = STPMX * std::max(std::sqrt(sum), static_cast<Real>(N));

			for (_iter = 0; _iter < _maxIter; ++_iter) {
				// Line search along xi (strong Wolfe conditions)
				Real fp = _fret;
				VectorN<Real, N> pold = p;

				// Perform line minimization with Wolfe conditions
				detail::QuasiNewtonLineSearchResult ls = LineSearchWolfe(func, p, fp, g, xi, stpmax);

				if (!ls.accepted) {
					// No step was taken, so an unchanged objective is not convergence evidence.
					// A small gradient, or a step that fell below the x-resolution along a finite
					// descent direction (NR TOLX), may still certify a minimum.
					_fret = fp;
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					if (detail::QuasiNewtonGradientConverged<N>(g, p, _fret, _gtol) || ls.stepBelowResolution)
						return MultidimOptimizationResult(result, _fret, _iter + 1, true);
					MultidimOptimizationResult failed(result, _fret, _iter + 1, false);
					failed.status = AlgorithmStatus::Stalled;
					failed.error_message = std::string("BFGS: line search failed (") + ls.reason + ")";
					return failed;
				}
				_fret = ls.f;

				// Check for convergence on function value
				if (2.0 * std::abs(_fret - fp) <= _ftol * (std::abs(_fret) + std::abs(fp) + EPS)) {
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					return MultidimOptimizationResult(result, _fret, _iter + 1, true);
				}

				// Update xi to be the actual step taken
				for (int j = 0; j < N; ++j)
					xi[j] = p[j] - pold[j];

				// Save old gradient and compute new gradient
				for (int j = 0; j < N; ++j)
					dg[j] = g[j];
				func.Gradient(p, g);

				// Check gradient convergence
				if (detail::QuasiNewtonGradientConverged<N>(g, p, _fret, _gtol)) {
					Vector<Real> result(N);
					for (int j = 0; j < N; ++j)
						result[j] = p[j];
					return MultidimOptimizationResult(result, _fret, _iter + 1, true);
				}

				// Compute gradient difference
				for (int j = 0; j < N; ++j)
					dg[j] = g[j] - dg[j];

				// Compute H * dg
				for (int i = 0; i < N; ++i) {
					hdg[i] = 0.0;
					for (int j = 0; j < N; ++j)
						hdg[i] += hessin(i, j) * dg[j];
				}

				// BFGS update of inverse Hessian
				Real fac = 0.0, fae = 0.0, sumdg = 0.0, sumxi = 0.0;
				for (int j = 0; j < N; ++j) {
					fac += dg[j] * xi[j];
					fae += dg[j] * hdg[j];
					sumdg += dg[j] * dg[j];
					sumxi += xi[j] * xi[j];
				}

				if (fac > std::sqrt(EPS * sumdg * sumxi)) {
					fac = MML::Real(1.0) / fac;
					Real fad = MML::Real(1.0) / fae;

					// Vector that makes BFGS different from DFP
					for (int j = 0; j < N; ++j)
						dg[j] = fac * xi[j] - fad * hdg[j];

					// BFGS formula for updating inverse Hessian
					for (int i = 0; i < N; ++i) {
						for (int j = i; j < N; ++j) {
							hessin(i, j) += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
							hessin(j, i) = hessin(i, j);
						}
					}
				}

				// Compute new search direction
				for (int i = 0; i < N; ++i) {
					xi[i] = 0.0;
					for (int j = 0; j < N; ++j)
						xi[i] -= hessin(i, j) * g[j];
				}
			}

			// Max iterations reached
			Vector<Real> result(N);
			for (int j = 0; j < N; ++j)
				result[j] = p[j];
			return MultidimOptimizationResult(result, _fret, _iter, false);
		}

	private:
		/**
		 * @brief Line search satisfying strong Wolfe conditions
		 * 
		 * Implements Nocedal & Wright Algorithms 3.5 (bracket phase) and 3.6 (zoom phase).
		 * Strong Wolfe conditions guarantee that the BFGS update produces a positive-definite
		 * Hessian approximation, which is essential for convergence on non-convex problems.
		 * 
		 * Conditions enforced:
		 *   Sufficient decrease (Armijo): f(x + α·d) ≤ f(x) + c₁·α·∇f(x)·d
		 *   Curvature condition:          |∇f(x + α·d)·d| ≤ c₂·|∇f(x)·d|
		 * 
		 * Nonfinite trial values are treated as "too far" and contract the bracket. If no
		 * step with sufficient decrease exists above the minimum meaningful step, the search
		 * reports failure and leaves x unchanged instead of returning a zero step.
		 * 
		 * @reference Nocedal & Wright, "Numerical Optimization", 2nd ed., Chapter 3
		 */
		template<int N>
		detail::QuasiNewtonLineSearchResult LineSearchWolfe(const IDifferentiableScalarFunction<N>& func, VectorN<Real, N>& x,
							 Real phi0, const VectorN<Real, N>& g, VectorN<Real, N>& d, Real stpmax) {
			const Real c1 = Precision::ArmijoSufficientDecrease;    // Sufficient decrease parameter
			const Real c2 = MML::Real(0.9);       // Curvature condition parameter (0.9 for quasi-Newton)
			const Real alpha_max = 50.0;
			const int max_ls_iter = 25;

			detail::QuasiNewtonLineSearchResult result;
			Real dphi0 = detail::PrepareDescentDirection<N>(g, d, stpmax);
			if (dphi0 >= 0.0) {
				result.reason = "zero gradient: no descent direction";
				return result;
			}
			const Real alpha_min = detail::MinimumMeaningfulStep<N>(x, d);

			Real alpha_prev = 0.0;
			Real phi_prev = phi0;
			Real dphi_prev = dphi0;
			Real alpha_cur = 1.0; // Full Newton step

			auto evaluate = [&](Real alpha, Real& phi) {
				VectorN<Real, N> xtrial;
				for (int j = 0; j < N; ++j)
					xtrial[j] = x[j] + alpha * d[j];
				phi = func(xtrial);
				return xtrial;
			};
			auto finish = [&](Real alpha, Real phi) {
				for (int j = 0; j < N; ++j)
					x[j] += alpha * d[j];
				result.f = phi;
				result.alpha = alpha;
				result.accepted = true;
				return result;
			};

			for (int i = 0; i < max_ls_iter; ++i) {
				Real phi_cur;
				VectorN<Real, N> xtrial = evaluate(alpha_cur, phi_cur);

				// Armijo violated, nonfinite trial, or no decrease since the previous trial:
				// a Wolfe point lies in (alpha_prev, alpha_cur)
				if (!std::isfinite(phi_cur) || phi_cur > phi0 + c1 * alpha_cur * dphi0 || (i > 0 && phi_cur >= phi_prev)) {
					return WolfeZoom(func, x, d, phi0, dphi0, c1, c2, alpha_min,
						alpha_prev, alpha_cur, phi_prev, dphi_prev, phi_cur, result);
				}

				// Evaluate phi'(alpha_cur)
				VectorN<Real, N> g_trial;
				func.Gradient(xtrial, g_trial);
				Real dphi_cur = 0.0;
				for (int j = 0; j < N; ++j)
					dphi_cur += g_trial[j] * d[j];
				if (!std::isfinite(dphi_cur)) {
					return WolfeZoom(func, x, d, phi0, dphi0, c1, c2, alpha_min,
						alpha_prev, alpha_cur, phi_prev, dphi_prev, phi_cur, result);
				}

				// Check curvature condition (strong Wolfe satisfied)
				if (std::abs(dphi_cur) <= -c2 * dphi0)
					return finish(alpha_cur, phi_cur);

				// Slope is positive — minimum is between alpha_prev and alpha_cur
				if (dphi_cur >= 0.0) {
					return WolfeZoom(func, x, d, phi0, dphi0, c1, c2, alpha_min,
						alpha_cur, alpha_prev, phi_cur, dphi_cur, phi_prev, result);
				}

				// Advance bracket
				if (alpha_cur >= alpha_max)
					return finish(alpha_cur, phi_cur); // Armijo holds at the trust boundary
				alpha_prev = alpha_cur;
				phi_prev = phi_cur;
				dphi_prev = dphi_cur;
				alpha_cur = std::min(Real(2.0) * alpha_cur, alpha_max);
			}

			// Fallback: the last trial satisfied Armijo (otherwise we would have zoomed)
			return finish(alpha_prev, phi_prev);
		}

		/**
		 * @brief Zoom phase of Wolfe line search (Nocedal & Wright Algorithm 3.6)
		 * 
		 * Finds a step length satisfying strong Wolfe conditions within [alpha_lo, alpha_hi].
		 * alpha_lo always satisfies Armijo (or is 0). Trial points come from safeguarded
		 * quadratic interpolation of (phi_lo, dphi_lo, phi_hi), clamped to the inner 80% of the
		 * bracket, falling back to bisection when the model is unusable; this lets badly scaled
		 * problems contract by orders of magnitude per iteration while staying bounded.
		 * Nonfinite trial values contract alpha_hi. When the bracket collapses below alpha_min
		 * with alpha_lo == 0, the search fails explicitly.
		 */
		template<int N>
		detail::QuasiNewtonLineSearchResult WolfeZoom(const IDifferentiableScalarFunction<N>& func, VectorN<Real, N>& x, const VectorN<Real, N>& d,
					   Real phi0, Real dphi0, Real c1, Real c2, Real alpha_min,
					   Real alpha_lo, Real alpha_hi, Real phi_lo, Real dphi_lo, Real phi_hi,
					   detail::QuasiNewtonLineSearchResult& result) {
			const int max_zoom_iter = 40;
			const Real tau = MML::Real(0.1); // safeguard fraction of the bracket

			auto finish = [&](Real alpha, Real phi) {
				for (int j = 0; j < N; ++j)
					x[j] += alpha * d[j];
				result.f = phi;
				result.alpha = alpha;
				result.accepted = true;
				return result;
			};
			auto fail = [&](const char* reason) {
				result.reason = reason;
				return result;
			};
			bool sawNonfinite = !std::isfinite(phi_hi);

			for (int i = 0; i < max_zoom_iter; ++i) {
				const Real lo = std::min(alpha_lo, alpha_hi);
				const Real hi = std::max(alpha_lo, alpha_hi);
				const Real width = hi - lo;
				if (width < alpha_min) {
					if (alpha_lo > 0)
						return finish(alpha_lo, phi_lo);
					result.stepBelowResolution = !sawNonfinite;
					return fail(sawNonfinite
						? "bracket collapsed below the minimum meaningful step while contracting away from nonfinite trials"
						: "bracket collapsed below the minimum meaningful step without sufficient decrease");
				}

				// Quadratic model through (alpha_lo, phi_lo, dphi_lo) and (alpha_hi, phi_hi)
				Real alpha_j = MML::Real(0.5) * (alpha_lo + alpha_hi);
				if (std::isfinite(phi_hi) && std::isfinite(dphi_lo)) {
					const Real h = alpha_hi - alpha_lo;
					const Real denom = phi_hi - phi_lo - dphi_lo * h;
					if (denom > 0) {
						Real candidate = alpha_lo - dphi_lo * h * h / (MML::Real(2.0) * denom);
						if (std::isfinite(candidate))
							alpha_j = std::clamp(candidate, lo + tau * width, hi - tau * width);
					}
				}

				VectorN<Real, N> xtrial;
				for (int j = 0; j < N; ++j)
					xtrial[j] = x[j] + alpha_j * d[j];
				Real phi_j = func(xtrial);

				if (!std::isfinite(phi_j) || phi_j > phi0 + c1 * alpha_j * dphi0 || phi_j >= phi_lo) {
					if (!std::isfinite(phi_j))
						sawNonfinite = true;
					alpha_hi = alpha_j;
					phi_hi = phi_j;
					continue;
				}

				VectorN<Real, N> g_trial;
				func.Gradient(xtrial, g_trial);
				Real dphi_j = 0.0;
				for (int j = 0; j < N; ++j)
					dphi_j += g_trial[j] * d[j];
				if (!std::isfinite(dphi_j)) {
					sawNonfinite = true;
					alpha_hi = alpha_j;
					phi_hi = phi_j;
					continue;
				}

				// Strong Wolfe satisfied
				if (std::abs(dphi_j) <= -c2 * dphi0)
					return finish(alpha_j, phi_j);

				if (dphi_j * (alpha_hi - alpha_lo) >= 0.0) {
					alpha_hi = alpha_lo;
					phi_hi = phi_lo;
				}
				alpha_lo = alpha_j;
				phi_lo = phi_j;
				dphi_lo = dphi_j;
			}

			// Budget exhausted: alpha_lo satisfies Armijo whenever it is positive
			if (alpha_lo > 0)
				return finish(alpha_lo, phi_lo);
			return fail("zoom budget exhausted without sufficient decrease");
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///               L-BFGS Limited-Memory Quasi-Newton Method             ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief L-BFGS (Limited-memory BFGS) for large-scale optimization
     * 
     * L-BFGS is a quasi-Newton method optimized for problems with many variables
     * where storing the full N×N inverse Hessian is impractical. Instead of
     * storing the dense matrix, it maintains only the last m correction pairs
     * and uses them to implicitly represent the inverse Hessian.
     * 
     * Memory Comparison:
     * - BFGS:   O(N²) for inverse Hessian matrix
     * - L-BFGS: O(m×N) where m is typically 3-20
     * 
     * For N=10,000 variables:
     * - BFGS:   ~800 MB (100M doubles)
     * - L-BFGS: ~1.6 MB with m=10 (200K doubles)
     * 
     * The algorithm uses the two-loop recursion of Nocedal to efficiently
     * compute the search direction H_k * g without forming H_k explicitly.
     * 
     * @see Jorge Nocedal, "Updating Quasi-Newton Matrices with Limited Storage"
     *      Mathematics of Computation, Vol. 35, No. 151 (1980), pp. 773-782
     */
	class LBFGS {
	private:
		Real _ftol;           ///< Function value tolerance
		Real _gtol;           ///< Gradient tolerance
		int _maxIter;         ///< Maximum iterations
		int _memorySize;      ///< Number of correction pairs to store (m)
		int _iter;            ///< Iteration counter
		Real _fret;           ///< Current function value

	public:
		/**
         * @brief Construct L-BFGS optimizer
         * @param ftol Function value tolerance for convergence
         * @param gtol Gradient tolerance for convergence  
         * @param maxIter Maximum number of iterations
         * @param memorySize Number of correction pairs to store (typical: 3-20)
         */
		LBFGS(Real ftol = PrecisionValues<Real>::OptimizationTolerance, Real gtol = PrecisionValues<Real>::OptimizationGradientTolerance, int maxIter = 1000, int memorySize = 10)
			: _ftol(ftol)
			, _gtol(gtol)
			, _maxIter(maxIter)
			, _memorySize(memorySize)
			, _iter(0)
			, _fret(0.0) {}

		// Getters and setters
		Real getFtol() const { return _ftol; }
		void setFtol(Real ftol) { _ftol = ftol; }

		Real getGtol() const { return _gtol; }
		void setGtol(Real gtol) { _gtol = gtol; }

		int getMaxIter() const { return _maxIter; }
		void setMaxIter(int maxIter) { _maxIter = maxIter; }

		int getMemorySize() const { return _memorySize; }
		void setMemorySize(int m) { _memorySize = m; }

		int getIterations() const { return _iter; }
		Real getCurrentFValue() const { return _fret; }

		/**
         * @brief Minimize using L-BFGS limited-memory quasi-Newton method
         * @tparam N Dimension of the problem
         * @param func Differentiable scalar function to minimize
         * @param start Starting point
         * @return Optimization result with minimum location and value
         * @throws MultidimOptimizationInputError if inputs are invalid
         */
		template<int N>
		MultidimOptimizationResult Minimize(const IDifferentiableScalarFunction<N>& func, const VectorN<Real, N>& start) {
			// Validate inputs
			ValidateVectorFinite<N>(start, "LBFGS::Minimize starting point");
			ValidateMultidimTolerance(_ftol, "LBFGS ftol");
			ValidateMultidimTolerance(_gtol, "LBFGS gtol");
			if (_maxIter <= 0)
				throw MultidimOptimizationInputError("LBFGS: maxIter must be positive");
			if (_memorySize <= 0)
				throw MultidimOptimizationInputError("LBFGS: memorySize must be positive");

			const Real EPS = Precision::Eps;   // NR dfpmin: machine epsilon
			const Real STPMX = 100.0;

			// Current position and gradient
			VectorN<Real, N> x = start;
			VectorN<Real, N> g, g_old, q, r;

			// Circular buffers for correction pairs
			// s_k = x_{k+1} - x_k (position difference)
			// y_k = g_{k+1} - g_k (gradient difference)
			std::vector<VectorN<Real, N>> s_history(_memorySize);
			std::vector<VectorN<Real, N>> y_history(_memorySize);
			std::vector<Real> rho(_memorySize);  // 1 / (y_k^T s_k)
			std::vector<Real> alpha(_memorySize);

			int history_count = 0;  // Number of stored pairs
			int history_start = 0;  // Circular buffer start index

			// Initial function and gradient evaluation
			_fret = func(x);
			ValidateMultidimFunctionValue(_fret, "LBFGS initial evaluation");
			func.Gradient(x, g);
			ValidateVectorFinite<N>(g, "LBFGS initial gradient");

			// Compute maximum step size
			Real sum = 0.0;
			for (int j = 0; j < N; ++j)
				sum += x[j] * x[j];
			Real stpmax = STPMX * std::max(std::sqrt(sum), static_cast<Real>(N));

			for (_iter = 0; _iter < _maxIter; ++_iter) {
				Real fp = _fret;
				VectorN<Real, N> x_old = x;
				g_old = g;

				// =========================================================
				// Two-loop recursion to compute search direction r = H * g
				// Based on Nocedal's algorithm (Algorithm 7.4 in Nocedal & Wright)
				// =========================================================
				q = g;

				// First loop: backward through history
				for (int i = history_count - 1; i >= 0; --i) {
					int idx = (history_start + i) % _memorySize;
					alpha[idx] = rho[idx] * DotProduct(s_history[idx], q);
					for (int j = 0; j < N; ++j)
						q[j] -= alpha[idx] * y_history[idx][j];
				}

				// Initial Hessian approximation: H_0 = gamma * I
				// gamma = (s_{k-1}^T y_{k-1}) / (y_{k-1}^T y_{k-1})
				Real gamma = 1.0;
				if (history_count > 0) {
					int last_idx = (history_start + history_count - 1) % _memorySize;
					Real yTy = DotProduct(y_history[last_idx], y_history[last_idx]);
					Real sTy = DotProduct(s_history[last_idx], y_history[last_idx]);
					if (yTy > EPS)
						gamma = sTy / yTy;
				}

				// r = H_0 * q = gamma * q
				for (int j = 0; j < N; ++j)
					r[j] = gamma * q[j];

				// Second loop: forward through history
				for (int i = 0; i < history_count; ++i) {
					int idx = (history_start + i) % _memorySize;
					Real beta = rho[idx] * DotProduct(y_history[idx], r);
					for (int j = 0; j < N; ++j)
						r[j] += s_history[idx][j] * (alpha[idx] - beta);
				}

				// Search direction is -H*g
				VectorN<Real, N> p;
				for (int j = 0; j < N; ++j)
					p[j] = -r[j];

				// =========================================================
				// Line search along direction p
				// =========================================================
				detail::QuasiNewtonLineSearchResult ls = detail::BacktrackingLineSearch<N>(func, x, fp, g, p, stpmax);

				if (!ls.accepted) {
					// No step was taken: unchanged f is not convergence evidence, but a small
					// gradient or a step below the x-resolution along a finite descent direction is
					_fret = fp;
					if (detail::QuasiNewtonGradientConverged<N>(g, x, _fret, _gtol) || ls.stepBelowResolution)
						return MakeResult(x, true);
					MultidimOptimizationResult failed = MakeResult(x, false);
					failed.status = AlgorithmStatus::Stalled;
					failed.error_message = std::string("LBFGS: line search failed (") + ls.reason + ")";
					return failed;
				}
				_fret = ls.f;

				// Check for convergence on function value
				if (2.0 * std::abs(_fret - fp) <= _ftol * (std::abs(_fret) + std::abs(fp) + EPS)) {
					return MakeResult(x, true);
				}

				// Compute new gradient
				func.Gradient(x, g);

				// Check gradient convergence
				if (detail::QuasiNewtonGradientConverged<N>(g, x, _fret, _gtol)) {
					return MakeResult(x, true);
				}

				// =========================================================
				// Update history with new correction pair
				// =========================================================
				VectorN<Real, N> s_new, y_new;
				for (int j = 0; j < N; ++j) {
					s_new[j] = x[j] - x_old[j];
					y_new[j] = g[j] - g_old[j];
				}

				Real sTy = DotProduct(s_new, y_new);
				if (sTy > EPS) {  // Curvature condition - only update if positive
					int store_idx;
					if (history_count < _memorySize) {
						store_idx = history_count;
						history_count++;
					} else {
						// Circular buffer: overwrite oldest
						store_idx = history_start;
						history_start = (history_start + 1) % _memorySize;
					}

					s_history[store_idx] = s_new;
					y_history[store_idx] = y_new;
					rho[store_idx] = MML::Real(1.0) / sTy;
				}
			}

			// Max iterations reached
			return MakeResult(x, false);
		}

	private:
		/// @brief Create result struct from current state
		template<int N>
		MultidimOptimizationResult MakeResult(const VectorN<Real, N>& x, bool converged) {
			Vector<Real> result(N);
			for (int j = 0; j < N; ++j)
				result[j] = x[j];
			return MultidimOptimizationResult(result, _fret, _iter + 1, converged);
		}

		/// @brief Compute dot product of two VectorN
		template<int N>
		static Real DotProduct(const VectorN<Real, N>& a, const VectorN<Real, N>& b) {
			Real sum = 0.0;
			for (int j = 0; j < N; ++j)
				sum += a[j] * b[j];
			return sum;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                   Convenience Wrapper Functions                     ///

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_QUASI_NEWTON_H
