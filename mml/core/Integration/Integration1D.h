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
			Real _currAbsSum;   ///< Same trapezoid stage applied to |f| (QUADPACK "resabs"): sets the scale of round-off
			Function& _func;

			TrapIntegratorCallable(Function& func, Real a, Real b)
				: _func(func), _a(a), _b(b), _currSum(0.0), _currAbsSum(0.0) { }

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
				Real x, sum, absSum, del;
				int subDivNum, j;

				_currStep++;
				if (_currStep == 1) {
					const Real fa = std::invoke(_func, _a);
					const Real fb = std::invoke(_func, _b);
					_currAbsSum = 0.5 * std::abs(_b - _a) * (std::abs(fa) + std::abs(fb));
					return (_currSum = 0.5 * (_b - _a) * (fa + fb));
				}
				else {
					for (subDivNum = 1, j = 1; j < _currStep - 1; j++)
						subDivNum *= 2;

					del = (_b - _a) / subDivNum;
					x = _a + 0.5 * del;

					for (sum = 0.0, absSum = 0.0, j = 0; j < subDivNum; j++, x += del) {
						const Real fx = std::invoke(_func, x);
						sum += fx;
						absSum += std::abs(fx);
					}

					_currSum = 0.5 * (_currSum + (_b - _a) * sum / subDivNum);
					_currAbsSum = 0.5 * (_currAbsSum + std::abs(_b - _a) * absSum / subDivNum);

					return _currSum;
				}
			}
		};

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateTrapImpl(Function& func, Real a, Real b,
			const Real eps = Defaults::TrapezoidIntegrationEPS)
		{
			int j;
			Real currSum, oldSum = 0.0;
			Real last_error = 0.0;

			TrapIntegratorCallable<Function> t(func, a, b);
			for (j = 0; j < Defaults::TrapezoidIntegrationMaxSteps; j++)
			{
				currSum = t.next();

				if (j > 5) {
					Real error = std::abs(currSum - oldSum);
					if (error <= t.tolerance(eps))   // <= so an exactly-zero integral converges
					{
						return IntegrationResult(currSum, error, j, true);
					}
				}

				last_error = std::abs(currSum - oldSum);
				oldSum = currSum;
			}
			return IntegrationResult(currSum, last_error, j, false);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateSimpsonImpl(Function& func, Real a, Real b,
			const Real eps = Defaults::SimpsonIntegrationEPS)
		{
			int j;
			Real currSum, st, ost = 0.0, oldSum = 0.0;
			Real last_error = 0.0;

			TrapIntegratorCallable<Function> t(func, a, b);
			for (j = 0; j < Defaults::SimpsonIntegrationMaxSteps; j++)
			{
				st = t.next();

				currSum = (4.0 * st - ost) / 3.0;

				if (j > 5) {
					Real error = std::abs(currSum - oldSum);
					if (error <= t.tolerance(eps))   // <= so an exactly-zero integral converges
					{
						return IntegrationResult(currSum, error, j, true);
					}
				}

				last_error = std::abs(currSum - oldSum);
				oldSum = currSum;
				ost = st;
			}
			return IntegrationResult(currSum, last_error, j, false);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateRombergImpl(Function& func, Real a, Real b,
			const Real eps = Defaults::RombergIntegrationEPS)
		{
			const int JMAX = Defaults::RombergIntegrationMaxSteps;
			const int K = Defaults::RombergIntegrationUsedPnts;

			std::vector<Real> s(JMAX);
			std::vector<Real> h(JMAX + 1);

			h[0] = 1.0;
			TrapIntegratorCallable<Function> trap(func, a, b);

			Real ss = 0.0;
			Real dss = 0.0;

			for (int j = 0; j < JMAX; j++)
			{
				s[j] = trap.next();

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
								return IntegrationResult(ss, std::abs(dss), j + 1, false);
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
						return IntegrationResult(ss, std::abs(dss), j + 1, true);
					}
				}

				h[j + 1] = 0.25 * h[j];
			}

			return IntegrationResult(ss, std::abs(dss), JMAX, false);
		}

		template<RealFunctionCallable Function>
		IntegrationResult IntegrateGauss10Impl(Function& func, const Real a, const Real b)
		{
			static const Real x[] = { 0.1488743389816312, 0.4333953941292472,
															0.6794095682990244, 0.8650633666889845, 0.9739065285171717 };
			static const Real w[] = { 0.2955242247147529, 0.2692667193099963,
															0.2190863625159821, 0.1494513491505806, 0.0666713443086881 };

			Real xm = 0.5 * (b + a);
			Real xr = 0.5 * (b - a);

			Real s = 0;
			for (int j = 0; j < 5; j++)
			{
				Real dx = xr * x[j];
				s += w[j] * (std::invoke(func, xm + dx) + std::invoke(func, xm - dx));
			}
			Real result = s * xr;

			return IntegrationResult::FixedRule(result);
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
	/// Uses a mixed absolute+relative test:
	///   |S_n - S_{n-1}| < eps * |S_{n-1}| + eps * machine_epsilon
	/// The first term (eps * |S|) provides relative accuracy for non-zero integrals.
	/// The second term (eps * machine_epsilon) provides an absolute floor so that
	/// near-zero integrals can still converge rather than demanding |error| < 0.
	static IntegrationResult IntegrateTrapCore(const IRealFunction& func, Real a, Real b,
																		 const Real eps = Defaults::TrapezoidIntegrationEPS)	{
		return IntegrationDetail::IntegrateTrapImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static auto IntegrateTrapDetailed(const IRealFunction& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateTrapCore(func, a, b, eps);
		r.algorithm_name = "IntegrateTrap";
		return r;
	}

	static Real IntegrateTrap(const IRealFunction& func, Real a, Real b,
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
	///   |S_n - S_{n-1}| < eps * |S_{n-1}| + eps * machine_epsilon
	static IntegrationResult IntegrateSimpsonCore(const IRealFunction& func, Real a, Real b,
																			 const Real eps = Defaults::SimpsonIntegrationEPS)
	{
		return IntegrationDetail::IntegrateSimpsonImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static auto IntegrateSimpsonDetailed(const IRealFunction& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateSimpsonCore(func, a, b, eps);
		r.algorithm_name = "IntegrateSimpson";
		return r;
	}

	static Real IntegrateSimpson(const IRealFunction& func, Real a, Real b,
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
	 *   |dss| <= eps * |ss| + eps * machine_epsilon
	 *
	 * @note Reference: Numerical Recipes §4.3 "Romberg Integration"
	 */
	static IntegrationResult IntegrateRombergCore(const IRealFunction& func, Real a, Real b,
																						const Real eps = Defaults::RombergIntegrationEPS)
	{
		return IntegrationDetail::IntegrateRombergImpl(func, a, b, eps);
	}

	/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
	static auto IntegrateRombergDetailed(const IRealFunction& func, Real a, Real b, const Real eps)
	{
		auto r = IntegrateRombergCore(func, a, b, eps);
		r.algorithm_name = "IntegrateRomberg";
		return r;
	}

	static Real IntegrateRomberg(const IRealFunction& func, Real a, Real b,
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
	 * @return IntegrationResult with value, error_estimate = 0 placeholder, error_estimated = false, 1 iteration
	 *
	 * @note Error estimate is 0 because this is a non-adaptive method.
	 *       Use adaptive methods (Romberg, Simpson) when error bounds are needed.
	 * @note Reference: Abramowitz & Stegun, Table 25.4 for nodes and weights
	 */
	static IntegrationResult IntegrateGauss10Core(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrationDetail::IntegrateGauss10Impl(func, a, b);
	}

	static Real IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)
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

	static Real IntegrateGK21(const IRealFunction& func, Real a, Real b)
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
	static auto IntegrateDetailed(const IRealFunction& func, Real a, Real b, IntegrationMethod method, Real eps)
	{
		auto r = IntegrateCore(func, a, b, method, eps);
		r.algorithm_name = "Integrate";
		return r;
	}

	static Real Integrate(const IRealFunction& func, Real a, Real b,
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
	static Real Integrate(const IRealFunction& func, Real a, Real b,
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
				auto r = IntegrateTrapCore(func, a, b, eps);
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
				auto r = IntegrateSimpsonCore(func, a, b, eps);
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
				auto r = IntegrateRombergCore(func, a, b, eps);
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
				auto r = IntegrateGauss10Core(func, a, b);
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
				auto r = IntegrateGK21Core(func, a, b);
				IntegrationDetail::PopulateFromSimple(result, r, "IntegrateGK21");
				result.function_evaluations = 21;
			});
	}

	/// Unified integration with runtime method selection and full diagnostics
	static IntegrationDetailedResult IntegrateDetailed(
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
	static IntegrationDetailedResult IntegrateDetailed(
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
}

#endif // MML_INTEGRATION_1D_H