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
///  Usage:       IntegrationResult r = IntegrateSimpsonCore(f, a, b, eps);               ///
///               if (r.converged) std::cout << r.value << " ± " << r.error_estimate; ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_INTEGRATION_1D_H
#define MML_INTEGRATION_1D_H

#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <mml/MMLBase.h>
#include "IntegrationBase.h"
#include "GaussKronrod.h"

#include <mml/interfaces/IFunction.h>
#include <mml/base/Function.h>

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

	namespace IntegrationDetail
	{
		template<RealFunctionCallable Function>
		class TrapIntegratorCallable : IQuadrature
		{
		public:
			Real _a, _b, _currSum;
			Real _currAbsSum;   ///< Sampled absolute scale: |f|, or an integrand-provided bound for internal cancellation.
			Function& _func;
			int _evaluations = 0;
			IntegrationResult* _progress;
			bool _checkFinite;
			Real _sampleAbsoluteScale = 0;

			TrapIntegratorCallable(Function& func, Real a, Real b, IntegrationResult* progress = nullptr, bool checkFinite = true)
				: _a(a), _b(b), _currSum(0.0), _currAbsSum(0.0), _func(func), _progress(progress), _checkFinite(checkFinite) {
				ValidateBounds(a, b, "Trapezoidal integration");
			}

			Real Evaluate(Real x) {
				auto counted = [&](Real point) {
					if (_progress) _progress->function_evaluations = _evaluations;
					if constexpr (requires {
						{ _func.EvaluateWithAbsoluteScale(point).first } -> std::convertible_to<Real>;
						{ _func.EvaluateWithAbsoluteScale(point).second } -> std::convertible_to<Real>;
					}) {
						const auto evaluated = _func.EvaluateWithAbsoluteScale(point);
						_sampleAbsoluteScale = evaluated.second;
						if (_sampleAbsoluteScale < 0)
							throw NumericalMethodError("Integration sample has an invalid absolute scale");
						return evaluated.first;
					} else {
						const Real value = std::invoke(_func, point);
						_sampleAbsoluteScale = std::abs(value);
						return value;
					}
				};
				return Sample(counted, _evaluations, _checkFinite, x);
			}

			/// Absolute tolerance floor: 50 ulp of the integral of |f| (as in QUADPACK/GSL).
			/// Negligible for integrals of normal size; for integrals that are (near) zero it is
			/// the level below which a difference of partial sums is pure round-off, so such
			/// integrals can converge at all (analysis v4, CORE C5: the old floor eps*machine_eps
			/// made every zero-valued integral run to the iteration limit and report converged=false).
			Real absoluteFloor() const {
				return Precision::IntegrationRoundoffFloorScale * _currAbsSum;
			}

			/// Convergence tolerance for a requested relative precision eps: eps is taken relative
			/// to the integral of |f| (QUADPACK's resabs scale). For integrands without sign change
			/// this is the usual eps*|I|; for cancelling integrals (|I| << int|f|) it is the absolute
			/// accuracy the method can honestly deliver - a pure eps*|I| criterion can never be met
			/// there and used to run every oscillatory integral to the iteration limit.
			Real tolerance(Real eps) const {
				return eps * _currAbsSum + absoluteFloor();
			}

			Real next() {
				int subDivNum, j;
				_currStep++;
				const Real halfWidth = HalfWidth(_a, _b);
				if (_currStep == 1) {
					const Real fa = Evaluate(_a);
					const Real faScale = _sampleAbsoluteScale;
					const Real fb = Evaluate(_b);
					_currAbsSum = Product({std::abs(halfWidth), faScale}) +
						Product({std::abs(halfWidth), _sampleAbsoluteScale});
					_currSum = Product({halfWidth, fa}) + Product({halfWidth, fb});
				}
				else {
					for (subDivNum = 1, j = 1; j < _currStep - 1; j++)
						subDivNum *= 2;

					Sum sum, absSum;
					for (j = 0; j < subDivNum; ++j) {
						const Real x = std::lerp(_a, _b, (Real(j) + Real(.5)) / subDivNum);
						const Real fx = Evaluate(x);
						const Real contribution = Product({halfWidth, Real(1) / subDivNum, fx});
						sum.Add(contribution);
						absSum.Add(Product({std::abs(halfWidth), Real(1) / subDivNum, _sampleAbsoluteScale}));
					}
					_currSum = MML::Real(.5) * _currSum + sum.Value();
					_currAbsSum = MML::Real(.5) * _currAbsSum + absSum.Value();
				}
				RequireFinite(_currSum, "Trapezoidal integration");
				RequireFinite(_currAbsSum, "Trapezoidal absolute scale");
				return _currSum;
			}
		};

		inline void ValidateRefinement(Real a, Real b, Real eps, int steps) {
			ValidateBounds(a, b, "Integration");
			ValidateTolerance(eps, "Integration");
			if (steps <= 0 || steps >= std::numeric_limits<int>::digits)
				throw ArgumentError("Integration refinement count exceeds the positive supported counter range");
		}

		template<class Quadrature>
		IntegrationResult RefinementResult(Quadrature& quadrature, Real value, Real error, int iteration, bool converged) {
			RequireFinite(value, "Integration refinement");
			if (!std::isnan(error)) RequireFinite(error, "Integration error");
			IntegrationResult result(value, error, iteration, converged, !std::isnan(error));
			result.function_evaluations = quadrature._evaluations;
			if (!converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Integration refinement budget exhausted before tolerance was met";
			}
			if (quadrature._progress) Adopt(*quadrature._progress, result, "Integration");
			return result;
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateTrapImpl(Function& func, Real a, Real b,
			const Real eps = Defaults::TrapezoidIntegrationEPS, IntegrationResult* progress = nullptr, bool checkFinite = true)
		{
			ValidateRefinement(a, b, eps, Defaults::TrapezoidIntegrationMaxSteps);
			int j;
			Real currSum = 0, oldSum = 0.0;
			Real last_error = NoEstimate;

			TrapIntegratorCallable<Function> t(func, a, b, progress, checkFinite);
			for (j = 0; j < Defaults::TrapezoidIntegrationMaxSteps; j++)
			{
				currSum = t.next();

				if (j > 5) {
					Real error = std::abs(currSum - oldSum);
					if (error <= t.tolerance(eps))   // <= so an exactly-zero integral converges
					{
						return RefinementResult(t, currSum, error, j + 1, true);
					}
				}

				last_error = j > 0 ? std::abs(currSum - oldSum) : NoEstimate;
				static_cast<void>(RefinementResult(t, currSum, last_error, j + 1, false));
				oldSum = currSum;
			}
			return RefinementResult(t, currSum, last_error, j, false);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateSimpsonImpl(Function& func, Real a, Real b,
			const Real eps = Defaults::SimpsonIntegrationEPS, IntegrationResult* progress = nullptr, bool checkFinite = true)
		{
			ValidateRefinement(a, b, eps, Defaults::SimpsonIntegrationMaxSteps);
			int j;
			Real currSum = 0, st, ost = 0.0, oldSum = 0.0;
			Real last_error = NoEstimate;

			TrapIntegratorCallable<Function> t(func, a, b, progress, checkFinite);
			for (j = 0; j < Defaults::SimpsonIntegrationMaxSteps; j++)
			{
				st = t.next();

				currSum = st + (st / 3 - ost / 3);

				if (j > 5) {
					Real error = std::abs(currSum - oldSum);
					if (error <= t.tolerance(eps))   // <= so an exactly-zero integral converges
					{
						return RefinementResult(t, currSum, error, j + 1, true);
					}
				}

				last_error = j > 0 ? std::abs(currSum - oldSum) : NoEstimate;
				static_cast<void>(RefinementResult(t, currSum, last_error, j + 1, false));
				oldSum = currSum;
				ost = st;
			}
			return RefinementResult(t, currSum, last_error, j, false);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateRombergImpl(Function& func, Real a, Real b,
			const Real eps = Defaults::RombergIntegrationEPS, IntegrationResult* progress = nullptr, bool checkFinite = true)
		{
			const int JMAX = Defaults::RombergIntegrationMaxSteps;
			const int K = Defaults::RombergIntegrationUsedPnts;
			ValidateRefinement(a, b, eps, JMAX);
			if (K < 2)
				throw ArgumentError("Romberg interpolation order must be at least two");

			std::vector<Real> s(JMAX);
			std::vector<Real> h(JMAX + 1);

			h[0] = 1.0;
			TrapIntegratorCallable<Function> trap(func, a, b, progress, checkFinite);

			Real ss = 0.0;
			Real dss = 0.0;

			for (int j = 0; j < JMAX; j++)
			{
				s[j] = trap.next();
				static_cast<void>(RefinementResult(trap, s[j], NoEstimate, j + 1, false));

				if (j >= K - 1)
				{
					std::vector<Real> c(K), d(K);

					int ns = 0;
					Real dif = std::abs(h[j - K + 1]);

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

					ss = s[j - K + 1 + ns];
					ns--;

					for (int m = 1; m < K; m++)
					{
						for (int i = 0; i < K - m; i++)
						{
							Real ho = h[j - K + 1 + i];
							Real hp = h[j - K + 1 + i + m];
							Real w = c[i + 1] - d[i];
							Real den = ho - hp;

							if (std::abs(den) < Precision::DivisionSafetyThreshold)
							{
								return RefinementResult(trap, ss, std::abs(dss), j + 1, false);
							}

							den = w / den;
							d[i] = hp * den;
							c[i] = ho * den;
						}

						if (2 * (ns + 1) < (K - m))
							dss = c[ns + 1];
						else
							dss = d[ns--];

						ss += dss;
					}
					if (std::abs(dss) <= trap.tolerance(eps))
					{
						return RefinementResult(trap, ss, std::abs(dss), j + 1, true);
					}
					static_cast<void>(RefinementResult(trap, ss, std::abs(dss), j + 1, false));
				}

				h[j + 1] = MML::Real(0.25) * h[j];
			}

			return JMAX < K ? RefinementResult(trap, s.back(), NoEstimate, JMAX, false)
				: RefinementResult(trap, ss, std::abs(dss), JMAX, false);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateGauss10Impl(Function& func, const Real a, const Real b,
			IntegrationResult* progress = nullptr, bool checkFinite = true)
		{
			static const Real x[] = { MML::Real(0.1488743389816312), MML::Real(0.4333953941292472),
															MML::Real(0.6794095682990244), MML::Real(0.8650633666889845), MML::Real(0.9739065285171717) };
			static const Real w[] = { MML::Real(0.2955242247147529), MML::Real(0.2692667193099963),
															MML::Real(0.2190863625159821), MML::Real(0.1494513491505806), MML::Real(0.0666713443086881) };

			const Real xr = HalfWidth(a, b);
			const Real xm = std::midpoint(a, b);

			Sum sum;
			int evaluations = 0;
			auto counted = [&](Real point) {
				if (progress) progress->function_evaluations = evaluations;
				return std::invoke(func, point);
			};
			for (int j = 0; j < 5; j++)
			{
				Real dx = xr * x[j];
				sum.Add(Product({xr, w[j], Sample(counted, evaluations, checkFinite, xm + dx)}));
				sum.Add(Product({xr, w[j], Sample(counted, evaluations, checkFinite, xm - dx)}));
			}
			RequireFinite(sum.Value(), "Gauss10 integration");
			return IntegrationResult::FixedRule(sum.Value(), evaluations);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateGK21Impl(Function& func, Real a, Real b)
		{
			auto gkResult = Integration::IntegrateGK21Core([&func](Real x) { return std::invoke(func, x); }, a, b);
			return gkResult;
		}
	}

	/// @brief Trapezoidal rule integrator with progressive refinement
	/// @details Each call to next() doubles the number of evaluation points
	class TrapIntegrator : public IntegrationDetail::TrapIntegratorCallable<const IRealFunction>
	{
	public:
		TrapIntegrator(const IRealFunction& func, Real a, Real b)
			: IntegrationDetail::TrapIntegratorCallable<const IRealFunction>(func, a, b) { }
	};

	/// @brief Extended trapezoidal rule with adaptive refinement
	/// @param func Function to integrate
	/// @param a,b Integration bounds
	/// @param eps Desired relative accuracy (default: 1e-6)
	/// @return IntegrationResult with value, error_estimate, iterations, converged flag
	///
	/// @par Convergence model
	/// Uses successive-estimate differences against eps times the sampled integral
	/// of |f|, plus its named arithmetic roundoff floor. This remains meaningful
	/// for cancelling integrals; it is not a certified black-box error bound.
	static IntegrationResult IntegrateTrapCore(const IRealFunction& func, Real a, Real b,
																		 const Real eps = Defaults::TrapezoidIntegrationEPS)	{
		return IntegrationDetail::IntegrateTrapImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static inline auto IntegrateTrapDetailed(const IRealFunction& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateTrapCore(func, a, b, eps);
		r.algorithm_name = "IntegrateTrap";
		return r;
	}

	static inline Real IntegrateTrap(const IRealFunction& func, Real a, Real b,
																		 const Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateTrapCore(func, a, b, eps).ValueOrThrow("IntegrateTrap");
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateTrapCore(Function&& func, Real a, Real b,
		const Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrationDetail::IntegrateTrapImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateTrapDetailed(Function&& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateTrapCore(std::forward<Function>(func), a, b, eps);
		r.algorithm_name = "IntegrateTrap";
		return r;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real IntegrateTrap(Function&& func, Real a, Real b,
		const Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateTrapCore(std::forward<Function>(func), a, b, eps).ValueOrThrow("IntegrateTrap");
	}

	/// @brief Simpson's rule with adaptive refinement
	/// @details Uses 4/3 Richardson extrapolation on trapezoidal estimates
	/// @param func Function to integrate
	/// @param a,b Integration bounds
	/// @param eps Desired relative accuracy (default: 1e-8)
	/// @return IntegrationResult with value, error_estimate, iterations, converged flag
	///
	/// @par Convergence model
	/// Same mixed absolute+relative test as IntegrateTrap:
	/// Uses the sampled absolute-integral scale and arithmetic floor described above.
	static IntegrationResult IntegrateSimpsonCore(const IRealFunction& func, Real a, Real b,
																			 const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrationDetail::IntegrateSimpsonImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static inline auto IntegrateSimpsonDetailed(const IRealFunction& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateSimpsonCore(func, a, b, eps);
		r.algorithm_name = "IntegrateSimpson";
		return r;
	}

	static inline Real IntegrateSimpson(const IRealFunction& func, Real a, Real b,
																			 const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrateSimpsonCore(func, a, b, eps).ValueOrThrow("IntegrateSimpson");
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateSimpsonCore(Function&& func, Real a, Real b,
		const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrationDetail::IntegrateSimpsonImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateSimpsonDetailed(Function&& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateSimpsonCore(std::forward<Function>(func), a, b, eps);
		r.algorithm_name = "IntegrateSimpson";
		return r;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real IntegrateSimpson(Function&& func, Real a, Real b,
		const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrateSimpsonCore(std::forward<Function>(func), a, b, eps).ValueOrThrow("IntegrateSimpson");
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
	 * @par Convergence model
	 * Same mixed absolute+relative test as Trap/Simpson:
	 *   |dss| <= eps * integral(|f|) + arithmetic roundoff floor
	 *
	 * @note Reference: Numerical Recipes §4.3 "Romberg Integration"
	 */
	static IntegrationResult IntegrateRombergCore(const IRealFunction& func, Real a, Real b,
																						const Real eps = Defaults::RombergIntegrationEPS)
	{
		return IntegrationDetail::IntegrateRombergImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static inline auto IntegrateRombergDetailed(const IRealFunction& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateRombergCore(func, a, b, eps);
		r.algorithm_name = "IntegrateRomberg";
		return r;
	}

	static inline Real IntegrateRomberg(const IRealFunction& func, Real a, Real b,
																						const Real eps = Defaults::RombergIntegrationEPS)
	{
		return IntegrateRombergCore(func, a, b, eps).ValueOrThrow("IntegrateRomberg");
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateRombergCore(Function&& func, Real a, Real b,
		const Real eps = Defaults::RombergIntegrationEPS)
	{
		return IntegrationDetail::IntegrateRombergImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateRombergDetailed(Function&& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateRombergCore(std::forward<Function>(func), a, b, eps);
		r.algorithm_name = "IntegrateRomberg";
		return r;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real IntegrateRomberg(Function&& func, Real a, Real b,
		const Real eps = Defaults::RombergIntegrationEPS)
	{
		return IntegrateRombergCore(std::forward<Function>(func), a, b, eps).ValueOrThrow("IntegrateRomberg");
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
	 * - No built-in error estimation (returned error_estimate is NoEstimate)
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
	 * @return IntegrationResult with value, NoEstimate, error_estimated = false, 1 iteration
	 *
	 * @note The error estimate is unavailable because this is a fixed rule.
	 *       Use adaptive methods (Romberg, Simpson) when error bounds are needed.
	 * @note Reference: Abramowitz & Stegun, Table 25.4 for nodes and weights
	 */
	static IntegrationResult IntegrateGauss10Core(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrationDetail::IntegrateGauss10Impl(func, a, b);
	}

	static inline Real IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrateGauss10Core(func, a, b).ValueOrThrow("IntegrateGauss10");
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateGauss10Core(Function&& func, const Real a, const Real b)
	{
		return IntegrationDetail::IntegrateGauss10Impl(func, a, b);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateGauss10Detailed(Function&& func, const Real a, const Real b)
	{
		auto r = IntegrateGauss10Core(std::forward<Function>(func), a, b);
		r.algorithm_name = "IntegrateGauss10";
		return r;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real IntegrateGauss10(Function&& func, const Real a, const Real b)
	{
		return IntegrateGauss10Core(std::forward<Function>(func), a, b).ValueOrThrow("IntegrateGauss10");
	}

	/**
	 * @brief 21-point Gauss-Kronrod integration (G10K21 rule)
	 * @details Uses 21-point Kronrod rule with embedded 10-point Gauss for error estimation.
	 *          Provides built-in error estimate from difference between Gauss and Kronrod results.
	 *          More accurate than GAUSS10 with only 21 function evaluations.
	 * @param func Function to integrate
	 * @param a Lower bound of integration
	 * @param b Upper bound of integration
	 * @return Fixed-rule IntegrationResult with embedded error estimate and Completed status.
	 */
	static IntegrationResult IntegrateGK21Core(const IRealFunction& func, Real a, Real b)
	{
		return IntegrationDetail::IntegrateGK21Impl(func, a, b);
	}

	static inline Real IntegrateGK21(const IRealFunction& func, Real a, Real b)
	{
		return IntegrateGK21Core(func, a, b).ValueOrThrow("IntegrateGK21");
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateGK21Core(Function&& func, Real a, Real b)
	{
		return IntegrationDetail::IntegrateGK21Impl(func, a, b);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateGK21Detailed(Function&& func, Real a, Real b)
	{
		auto r = IntegrateGK21Core(std::forward<Function>(func), a, b);
		r.algorithm_name = "IntegrateGK21";
		return r;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real IntegrateGK21(Function&& func, Real a, Real b)
	{
		return IntegrateGK21Core(std::forward<Function>(func), a, b).ValueOrThrow("IntegrateGK21");
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
	 *   auto result = IntegrateCore(f, 0.0, 1.0, SIMPSON, 1e-8);
	 *   if (!result.converged) { ... }
	 * @endcode
	 */
	static IntegrationResult IntegrateCore(const IRealFunction& func, Real a, Real b,
	                                   IntegrationMethod method = TRAP,
	                                   Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		switch (method)
		{
		case TRAP:
			return IntegrateTrapCore(func, a, b, eps);
		case SIMPSON:
			return IntegrateSimpsonCore(func, a, b, eps);
		case ROMBERG:
			return IntegrateRombergCore(func, a, b, eps);
		case GAUSS10:
			return IntegrateGauss10Core(func, a, b);
		case GAUSS10KRONROD21:
			return IntegrateGK21Core(func, a, b);
		default:
			return IntegrateTrapCore(func, a, b, eps);
		}
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static inline auto IntegrateDetailed(const IRealFunction& func, Real a, Real b, IntegrationMethod method, Real eps)
	{
		auto r = IntegrateCore(func, a, b, method, eps);
		r.algorithm_name = "Integrate";
		return r;
	}

	static inline Real Integrate(const IRealFunction& func, Real a, Real b,
	                                   IntegrationMethod method = TRAP,
	                                   Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateCore(func, a, b, method, eps).ValueOrThrow("Integrate");
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateCore(Function&& func, Real a, Real b,
		IntegrationMethod method = TRAP,
		Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		switch (method)
		{
		case TRAP:
			return IntegrateTrapCore(std::forward<Function>(func), a, b, eps);
		case SIMPSON:
			return IntegrateSimpsonCore(std::forward<Function>(func), a, b, eps);
		case ROMBERG:
			return IntegrateRombergCore(std::forward<Function>(func), a, b, eps);
		case GAUSS10:
			return IntegrateGauss10Core(std::forward<Function>(func), a, b);
		case GAUSS10KRONROD21:
			return IntegrateGK21Core(std::forward<Function>(func), a, b);
		default:
			return IntegrateTrapCore(std::forward<Function>(func), a, b, eps);
		}
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateDetailed(Function&& func, Real a, Real b, IntegrationMethod method, Real eps)
	{
		auto r = IntegrateCore(std::forward<Function>(func), a, b, method, eps);
		r.algorithm_name = "Integrate";
		return r;
	}

	template<RealFunctionCallable Function>
		requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real Integrate(Function&& func, Real a, Real b,
		IntegrationMethod method = TRAP,
		Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateCore(std::forward<Function>(func), a, b, method, eps).ValueOrThrow("Integrate");
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
	static IntegrationResult IntegrateCore(const IRealFunction& func, Real a, Real b,
	                                   Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		if constexpr (Method == TRAP)
			return IntegrateTrapCore(func, a, b, eps);
		else if constexpr (Method == SIMPSON)
			return IntegrateSimpsonCore(func, a, b, eps);
		else if constexpr (Method == ROMBERG)
			return IntegrateRombergCore(func, a, b, eps);
		else if constexpr (Method == GAUSS10)
			return IntegrateGauss10Core(func, a, b);
		else if constexpr (Method == GAUSS10KRONROD21)
			return IntegrateGK21Core(func, a, b);
		else
			return IntegrateTrapCore(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<IntegrationMethod Method = TRAP>
	static auto IntegrateDetailed(const IRealFunction& func, Real a, Real b, Real eps)
	{
		auto r = IntegrateCore<Method>(func, a, b, eps);
		r.algorithm_name = "Integrate";
		return r;
	}

	template<IntegrationMethod Method = TRAP>
	static inline Real Integrate(const IRealFunction& func, Real a, Real b,
	                                   Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateCore<Method>(func, a, b, eps).ValueOrThrow("Integrate");
	}

	template<IntegrationMethod Method = TRAP, class Function>
		requires RealFunctionCallable<Function> && (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static IntegrationResult IntegrateCore(Function&& func, Real a, Real b,
		Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		if constexpr (Method == TRAP)
			return IntegrateTrapCore(std::forward<Function>(func), a, b, eps);
		else if constexpr (Method == SIMPSON)
			return IntegrateSimpsonCore(std::forward<Function>(func), a, b, eps);
		else if constexpr (Method == ROMBERG)
			return IntegrateRombergCore(std::forward<Function>(func), a, b, eps);
		else if constexpr (Method == GAUSS10)
			return IntegrateGauss10Core(std::forward<Function>(func), a, b);
		else if constexpr (Method == GAUSS10KRONROD21)
			return IntegrateGK21Core(std::forward<Function>(func), a, b);
		else
			return IntegrateTrapCore(std::forward<Function>(func), a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	template<IntegrationMethod Method = TRAP, class Function>
		requires RealFunctionCallable<Function> && (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static auto IntegrateDetailed(Function&& func, Real a, Real b, Real eps)
	{
		auto r = IntegrateCore<Method>(std::forward<Function>(func), a, b, eps);
		r.algorithm_name = "Integrate";
		return r;
	}

	template<IntegrationMethod Method = TRAP, class Function>
		requires RealFunctionCallable<Function> && (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
	static Real Integrate(Function&& func, Real a, Real b,
		Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrateCore<Method>(std::forward<Function>(func), a, b, eps).ValueOrThrow("Integrate");
	}

	/******************************************************************************/
	/*****              Detailed API - 1D Integration                         *****/
	/******************************************************************************/

	/// Helper to populate IntegrationDetailedResult from IntegrationResult and set
	/// failure status when the underlying method did not converge.
	namespace IntegrationDetail
	{
		/// @deprecated Thin alias of Adopt (kept for the single-header manifest order).
		inline void PopulateFromSimple(IntegrationResult& out, const IntegrationResult& r, const char* algorithm_name)
		{
			Adopt(out, r, algorithm_name);
		}
	} // namespace IntegrationDetail

	/// Trapezoidal integration with full diagnostics
	static IntegrationDetailedResult IntegrateTrapDetailed(
		const IRealFunction& func, Real a, Real b,
		const IntegrationConfig& config = {},
		const Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"IntegrateTrap", config,
			[&](IntegrationDetailedResult& result) {
				auto r = IntegrationDetail::IntegrateTrapImpl(func, a, b, eps, &result, config.check_finite);
				IntegrationDetail::PopulateFromSimple(result, r, "IntegrateTrap");
			});
	}

	/// Simpson integration with full diagnostics
	static IntegrationDetailedResult IntegrateSimpsonDetailed(
		const IRealFunction& func, Real a, Real b,
		const IntegrationConfig& config = {},
		const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"IntegrateSimpson", config,
			[&](IntegrationDetailedResult& result) {
				auto r = IntegrationDetail::IntegrateSimpsonImpl(func, a, b, eps, &result, config.check_finite);
				IntegrationDetail::PopulateFromSimple(result, r, "IntegrateSimpson");
			});
	}

	/// Romberg integration with full diagnostics
	static IntegrationDetailedResult IntegrateRombergDetailed(
		const IRealFunction& func, Real a, Real b,
		const IntegrationConfig& config = {},
		const Real eps = Defaults::RombergIntegrationEPS)
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"IntegrateRomberg", config,
			[&](IntegrationDetailedResult& result) {
				auto r = IntegrationDetail::IntegrateRombergImpl(func, a, b, eps, &result, config.check_finite);
				IntegrationDetail::PopulateFromSimple(result, r, "IntegrateRomberg");
			});
	}

	/// Gauss10 quadrature with full diagnostics
	static IntegrationDetailedResult IntegrateGauss10Detailed(
		const IRealFunction& func, Real a, Real b,
		const IntegrationConfig& config = {})
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"IntegrateGauss10", config,
			[&](IntegrationDetailedResult& result) {
				auto r = IntegrationDetail::IntegrateGauss10Impl(func, a, b, &result, config.check_finite);
				IntegrationDetail::PopulateFromSimple(result, r, "IntegrateGauss10");
				result.function_evaluations = 10;
			});
	}

	/// GK21 integration with full diagnostics
	static IntegrationDetailedResult IntegrateGK21Detailed(
		const IRealFunction& func, Real a, Real b,
		const IntegrationConfig& config = {})
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"IntegrateGK21", config,
			[&](IntegrationDetailedResult& result) {
				auto r = Integration::IntegrateGK21Impl([&](Real x) { return func(x); }, a, b, &result, config.check_finite);
				IntegrationDetail::PopulateFromSimple(result, r, "IntegrateGK21");
				result.function_evaluations = 21;
			});
	}

	/// Unified integration with runtime method selection and full diagnostics
	static inline IntegrationDetailedResult IntegrateDetailed(
		const IRealFunction& func, Real a, Real b,
		IntegrationMethod method = TRAP,
		const IntegrationConfig& config = {},
		Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		switch (method)
		{
		case TRAP:    return IntegrateTrapDetailed(func, a, b, config, eps);
		case SIMPSON: return IntegrateSimpsonDetailed(func, a, b, config, eps);
		case ROMBERG: return IntegrateRombergDetailed(func, a, b, config, eps);
		case GAUSS10: return IntegrateGauss10Detailed(func, a, b, config);
		case GAUSS10KRONROD21: return IntegrateGK21Detailed(func, a, b, config);
		default:      return IntegrateTrapDetailed(func, a, b, config, eps);
		}
	}

	/// Template-based integration with compile-time method selection and full diagnostics
	template<IntegrationMethod Method = TRAP>
	static inline IntegrationDetailedResult IntegrateDetailed(
		const IRealFunction& func, Real a, Real b,
		const IntegrationConfig& config = {},
		Real eps = Defaults::TrapezoidIntegrationEPS)
	{
		if constexpr (Method == TRAP)
			return IntegrateTrapDetailed(func, a, b, config, eps);
		else if constexpr (Method == SIMPSON)
			return IntegrateSimpsonDetailed(func, a, b, config, eps);
		else if constexpr (Method == ROMBERG)
			return IntegrateRombergDetailed(func, a, b, config, eps);
		else if constexpr (Method == GAUSS10)
			return IntegrateGauss10Detailed(func, a, b, config);
		else if constexpr (Method == GAUSS10KRONROD21)
			return IntegrateGK21Detailed(func, a, b, config);
		else
			return IntegrateTrapDetailed(func, a, b, config, eps);
	}
	namespace IntegrationDetail {
		template<class Compute>
		IntegrationResult CountedDetailed(const char* name, const IRealFunction& function,
			const IntegrationConfig& config, Compute&& compute) {
			int evaluations = 0;
			auto result = ExecuteIntegrationDetailed<IntegrationResult>(name, config, [&](IntegrationResult& out) {
				auto sample = [&](Real x) { return Sample(function, evaluations, config.check_finite, x); };
				::MML::Detail::RealFunctionCallableAdapter<decltype(sample)> adapter(sample);
				compute(adapter, out);
			});
			result.function_evaluations = evaluations;
			return result;
		}

		inline Real MethodTolerance(IntegrationMethod method) {
			switch (method) {
			case TRAP: return Defaults::TrapezoidIntegrationEPS;
			case SIMPSON: return Defaults::SimpsonIntegrationEPS;
			case ROMBERG: return Defaults::RombergIntegrationEPS;
			case GAUSS10:
			case GAUSS10KRONROD21: return Precision::IntegrationTolerance;
			}
			throw ArgumentError("Integration: unsupported method");
		}

		template<class Function>
		IntegrationResult ConfiguredAxis(Function& function, Real a, Real b, IntegrationMethod method,
			Real tolerance, const IntegrationConfig& config, IntegrationResult* progress = nullptr) {
			switch (method) {
			case TRAP: return IntegrateTrapImpl(function, a, b, tolerance, progress, config.check_finite);
			case SIMPSON: return IntegrateSimpsonImpl(function, a, b, tolerance, progress, config.check_finite);
			case ROMBERG: return IntegrateRombergImpl(function, a, b, tolerance, progress, config.check_finite);
			case GAUSS10: return IntegrateGauss10Impl(function, a, b, progress, config.check_finite);
			case GAUSS10KRONROD21:
				return Integration::IntegrateGK21Impl([&](Real x) { return function(x); }, a, b, progress, config.check_finite);
			}
			throw ArgumentError("Integration: unsupported method");
		}

		inline IntegrationResult NestedAxis(const IRealFunction& function, Real a, Real b, IntegrationMethod method,
			const IntegrationConfig& config, IntegrationResult* progress = nullptr) {
			return ConfiguredAxis(function, a, b, method, MethodTolerance(method), config, progress);
		}

		struct NestedErrors {
			Real maximum = 0;
			bool allEstimated = true;
			AlgorithmStatus failure = AlgorithmStatus::Completed;
			std::string message;

			void Observe(const IntegrationResult& inner) {
				if (!inner.IsUsable() && MML::IsUsable(failure)) {
					failure = inner.status;
					message = inner.error_message.empty() ? "An inner integral did not converge" : inner.error_message;
				}
				allEstimated = allEstimated && inner.error_estimated;
				if (inner.error_estimated) maximum = std::max(maximum, inner.error_estimate);
			}
			void Combine(IntegrationResult& result, Real a, Real b, IntegrationMethod method) const {
				if (result.error_estimated && allEstimated) {
					result.error_estimate += Product({2, std::abs(HalfWidth(a, b)), maximum});
					RequireFinite(result.error_estimate, "Nested integration error");
				} else {
					result.error_estimated = false;
					result.error_estimate = NoEstimate;
				}
				if (!MML::IsUsable(failure)) {
					result.status = failure;
					result.converged = false;
					result.error_message = "Nested integration: " + message;
				} else if (result.converged && (!result.error_estimated ||
					result.error_estimate > MethodTolerance(method) * std::abs(result.value) +
						Precision::IntegrationRoundoffFloorScale * std::abs(result.value))) {
					result.status = AlgorithmStatus::MaxIterationsExceeded;
					result.converged = false;
					result.error_message = "Nested integration aggregate uncertainty exceeds the requested accuracy";
				}
			}
		};

		struct EvaluationCountGuard {
			IntegrationResult& result;
			const int& count;
			~EvaluationCountGuard() { result.function_evaluations = count; }
		};
	}
}

#endif // MML_INTEGRATION_1D_H