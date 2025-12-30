///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IntegrationImproper.h                                               ///
///  Description: Improper integrals with infinite bounds                             ///
///               Variable substitution for semi-infinite and infinite intervals      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_INTEGRATION_IMPROPER_H
#define MML_INTEGRATION_IMPROPER_H

#include <cmath>
#include <limits>

#include "MMLBase.h"
#include "interfaces/IFunction.h"
#include "Integration1D.h"

namespace MML
{
	/**
	 * @brief Open midpoint quadrature for use with improper integrals.
	 *
	 * Unlike the trapezoidal rule, the midpoint rule never evaluates the function
	 * at the endpoints. This is essential for improper integrals where the integrand
	 * may be singular or undefined at the boundaries.
	 *
	 * The midpoint rule uses the extended open formula:
	 *   ∫[a,b] f(x) dx ≈ h * Σ f(xᵢ)  where xᵢ = a + (i + 0.5)h
	 *
	 * Refinement is done by tripling the number of points (not doubling),
	 * which allows reuse of previous function evaluations.
	 */
	class MidpointIntegrator : public IQuadrature
	{
	public:
		Real _a, _b, _currSum;
		const IRealFunction& _func;

		MidpointIntegrator(const IRealFunction& func, Real a, Real b)
			: _func(func), _a(a), _b(b), _currSum(0.0) {}

		Real next() override {
			int it, j;
			Real x, tnm, sum, del, ddel;

			_currStep++;
			if (_currStep == 1) {
				return (_currSum = (_b - _a) * _func(0.5 * (_a + _b)));
			}
			else {
				for (it = 1, j = 1; j < _currStep - 1; j++) 
					it *= 3;  // Triple the subdivisions each step
				tnm = it;
				del = (_b - _a) / (3.0 * tnm);
				ddel = del + del;
				x = _a + 0.5 * del;
				sum = 0.0;
				for (j = 0; j < it; j++) {
					sum += _func(x);
					x += ddel;
					sum += _func(x);
					x += del;
				}
				_currSum = (_currSum + (_b - _a) * sum / tnm) / 3.0;
				return _currSum;
			}
		}
	};

	/**
	 * @brief Transformed integrator for ∫[a, ∞) f(x) dx
	 *
	 * Uses the substitution x = a + 1/t - 1, transforming [a, ∞) to (0, 1].
	 * The transformed integrand is f(a + 1/t - 1) / t²
	 */
	class UpperInfTransform : public IRealFunction
	{
	private:
		const IRealFunction& _func;
		Real _a;  // Lower bound
		
	public:
		UpperInfTransform(const IRealFunction& func, Real a) 
			: _func(func), _a(a) {}

		Real operator()(Real t) const override {
			if (t <= 0.0) return 0.0;  // Limit as t → 0⁺
			Real x = _a + 1.0/t - 1.0;
			return _func(x) / (t * t);
		}
	};

	/**
	 * @brief Transformed integrator for ∫(-∞, b] f(x) dx
	 *
	 * Uses the substitution x = b - 1/t + 1, transforming (-∞, b] to (0, 1].
	 * The transformed integrand is f(b - 1/t + 1) / t²
	 */
	class LowerInfTransform : public IRealFunction
	{
	private:
		const IRealFunction& _func;
		Real _b;  // Upper bound
		
	public:
		LowerInfTransform(const IRealFunction& func, Real b) 
			: _func(func), _b(b) {}

		Real operator()(Real t) const override {
			if (t <= 0.0) return 0.0;  // Limit as t → 0⁺
			Real x = _b - 1.0/t + 1.0;
			return _func(x) / (t * t);
		}
	};

	/**
	 * @brief Open Romberg integration using midpoint rule.
	 *
	 * Similar to standard Romberg but uses the midpoint rule which never
	 * evaluates at endpoints. Essential for improper integrals.
	 *
	 * The step size reduction is h → h/9 (not h/4) because midpoint
	 * refinement triples subdivisions (error goes as h², so h²→h²/9).
	 */
	static IntegrationResult IntegrateOpen(const IRealFunction& func, Real a, Real b,
	                                       const Real eps = Defaults::RombergIntegrationEPS)
	{
		const int JMAX = Defaults::RombergIntegrationMaxSteps;
		const int K = Defaults::RombergIntegrationUsedPnts;
		
		std::vector<Real> s(JMAX);
		std::vector<Real> h(JMAX + 1);
		
		h[0] = 1.0;
		MidpointIntegrator mid(func, a, b);
		
		Real ss = 0.0;
		Real dss = 0.0;
		
		for (int j = 0; j < JMAX; j++)
		{
			s[j] = mid.next();
			
			if (j >= K - 1)
			{
				// Neville's algorithm for polynomial extrapolation to h=0
				std::vector<Real> c(K), d(K);
				int ns = 0;
				Real dif = std::abs(h[j - K + 1]);
				
				for (int i = 0; i < K; i++)
				{
					Real dift = std::abs(h[j - K + 1 + i]);
					if (dift < dif) { ns = i; dif = dift; }
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
						if (den == 0.0)
							return IntegrationResult(ss, std::abs(dss), j + 1, false);
						den = w / den;
						d[i] = hp * den;
						c[i] = ho * den;
					}
					dss = (2 * (ns + 1) < (K - m)) ? c[ns + 1] : d[ns--];
					ss += dss;
				}
				
				if (std::abs(dss) <= eps * std::abs(ss) || (ss == 0.0 && dss == 0.0))
					return IntegrationResult(ss, std::abs(dss), j + 1, true);
			}
			
			// h → h/9 for midpoint rule (not h/4 as in trapezoidal)
			h[j + 1] = h[j] / 9.0;
		}
		
		return IntegrationResult(ss, std::abs(dss), JMAX, false);
	}

	/*************************************************************************/
	/*****                  SEMI-INFINITE INTEGRALS                      *****/
	/*************************************************************************/

	/**
	 * @brief Integrate f(x) from a to infinity: ∫[a, ∞) f(x) dx
	 *
	 * Computes the improper integral over a semi-infinite interval extending
	 * to positive infinity. Uses variable substitution to transform the
	 * infinite interval to a finite one.
	 *
	 * MATHEMATICAL TRANSFORMATION:
	 * The substitution t = 1/(x - a + 1) maps [a, ∞) → (0, 1]:
	 *   - x = a + 1/t - 1
	 *   - dx = -dt/t²
	 *   - ∫[a,∞) f(x) dx = ∫(0,1] f(a + 1/t - 1) / t² dt
	 *
	 * CONVERGENCE:
	 * The integral converges if f(x) = O(1/x^(1+ε)) as x → ∞ for some ε > 0.
	 * Faster decay (e.g., exponential) leads to faster numerical convergence.
	 *
	 * @param func Function to integrate (must be well-behaved for x ≥ a)
	 * @param a Lower bound of integration (finite)
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫[0,∞) e^(-x) dx = 1
	 *   RealFunction f{ [](Real x) { return std::exp(-x); } };
	 *   auto result = IntegrateUpperInf(f, 0.0);
	 *
	 *   // ∫[1,∞) 1/x² dx = 1
	 *   RealFunction g{ [](Real x) { return 1.0/(x*x); } };
	 *   auto result = IntegrateUpperInf(g, 1.0);
	 */
	static IntegrationResult IntegrateUpperInf(const IRealFunction& func, Real a,
	                                           const Real eps = Defaults::RombergIntegrationEPS)
	{
		UpperInfTransform transformed(func, a);
		return IntegrateOpen(transformed, 0.0, 1.0, eps);
	}

	/**
	 * @brief Integrate f(x) from negative infinity to b: ∫(-∞, b] f(x) dx
	 *
	 * Computes the improper integral over a semi-infinite interval extending
	 * to negative infinity.
	 *
	 * MATHEMATICAL TRANSFORMATION:
	 * The substitution t = 1/(b - x + 1) maps (-∞, b] → (0, 1]:
	 *   - x = b - 1/t + 1
	 *   - dx = dt/t²
	 *   - ∫(-∞,b] f(x) dx = ∫(0,1] f(b - 1/t + 1) / t² dt
	 *
	 * @param func Function to integrate (must be well-behaved for x ≤ b)
	 * @param b Upper bound of integration (finite)
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫(-∞,0] e^x dx = 1
	 *   RealFunction f{ [](Real x) { return std::exp(x); } };
	 *   auto result = IntegrateLowerInf(f, 0.0);
	 */
	static IntegrationResult IntegrateLowerInf(const IRealFunction& func, Real b,
	                                           const Real eps = Defaults::RombergIntegrationEPS)
	{
		LowerInfTransform transformed(func, b);
		return IntegrateOpen(transformed, 0.0, 1.0, eps);
	}

	/**
	 * @brief Integrate f(x) from negative infinity to positive infinity: ∫(-∞, ∞) f(x) dx
	 *
	 * Computes the improper integral over the entire real line. The integral
	 * is split at x = 0 into two semi-infinite integrals.
	 *
	 * ALGORITHM:
	 * ∫(-∞, ∞) f(x) dx = ∫(-∞, 0] f(x) dx + ∫[0, ∞) f(x) dx
	 *
	 * The split point is x = 0 by default, which works well for symmetric
	 * functions. For asymmetric functions, consider using the overload
	 * with a custom split point.
	 *
	 * CONVERGENCE:
	 * Requires f(x) → 0 fast enough as x → ±∞:
	 *   - f(x) = O(1/|x|^(1+ε)) is sufficient
	 *   - Gaussian decay e^(-x²) converges very fast
	 *   - Lorentzian decay 1/(1+x²) converges well
	 *
	 * @param func Function to integrate
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫(-∞,∞) e^(-x²) dx = √π
	 *   RealFunction gaussian{ [](Real x) { return std::exp(-x*x); } };
	 *   auto result = IntegrateInf(gaussian);
	 *   // result.value ≈ 1.7724538509...
	 *
	 *   // ∫(-∞,∞) 1/(1+x²) dx = π
	 *   RealFunction lorentz{ [](Real x) { return 1.0/(1.0 + x*x); } };
	 *   auto result = IntegrateInf(lorentz);
	 */
	static IntegrationResult IntegrateInf(const IRealFunction& func,
	                                      const Real eps = Defaults::RombergIntegrationEPS)
	{
		// Split at x = 0
		auto left  = IntegrateLowerInf(func, 0.0, eps);
		auto right = IntegrateUpperInf(func, 0.0, eps);
		
		return IntegrationResult(
			left.value + right.value,
			left.error_estimate + right.error_estimate,
			left.iterations + right.iterations,
			left.converged && right.converged
		);
	}

	/**
	 * @brief Integrate f(x) from negative infinity to positive infinity with custom split point
	 *
	 * For functions that are not symmetric about x = 0, or have features
	 * concentrated in a particular region, specifying a custom split point
	 * can improve accuracy and convergence.
	 *
	 * @param func Function to integrate
	 * @param splitPoint The x value at which to split the integral
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // Function peaked near x = 5
	 *   RealFunction peaked{ [](Real x) { return std::exp(-(x-5)*(x-5)); } };
	 *   auto result = IntegrateInfSplit(peaked, 5.0);  // Split near the peak
	 */
	static IntegrationResult IntegrateInfSplit(const IRealFunction& func, Real splitPoint,
	                                           const Real eps = Defaults::RombergIntegrationEPS)
	{
		auto left  = IntegrateLowerInf(func, splitPoint, eps);
		auto right = IntegrateUpperInf(func, splitPoint, eps);
		
		return IntegrationResult(
			left.value + right.value,
			left.error_estimate + right.error_estimate,
			left.iterations + right.iterations,
			left.converged && right.converged
		);
	}

	/*************************************************************************/
	/*****              SINGULAR INTEGRALS - ENDPOINT SINGULARITIES      *****/
	/*************************************************************************/

	/**
	 * @brief Transformed integrator for lower endpoint singularity
	 *
	 * For integrals with an integrable singularity at x = a, like:
	 *   ∫[a, b] f(x)/√(x-a) dx  or  ∫[a, b] f(x)/(x-a)^α dx  (0 < α < 1)
	 *
	 * Uses the substitution t² = x - a, which maps [a, b] → [0, √(b-a)]:
	 *   - x = a + t²
	 *   - dx = 2t dt
	 *   - √(x-a) = t
	 *   - ∫[a,b] g(x) dx = ∫[0,√(b-a)] g(a + t²) · 2t dt
	 *
	 * The factor of 2t cancels the 1/√(x-a) = 1/t singularity, making the
	 * integrand finite and well-behaved at t = 0.
	 *
	 * This transformation is ideal for inverse square-root singularities
	 * and works well for weaker singularities too.
	 */
	class LowerSingularTransform : public IRealFunction
	{
	private:
		const IRealFunction& _func;
		Real _a;  // Singular point (lower bound)
		
	public:
		LowerSingularTransform(const IRealFunction& func, Real a) 
			: _func(func), _a(a) {}

		Real operator()(Real t) const override {
			// x = a + t², dx = 2t dt
			Real x = _a + t * t;
			return _func(x) * 2.0 * t;
		}
	};

	/**
	 * @brief Transformed integrator for upper endpoint singularity
	 *
	 * For integrals with an integrable singularity at x = b, like:
	 *   ∫[a, b] f(x)/√(b-x) dx  or  ∫[a, b] f(x)/(b-x)^α dx  (0 < α < 1)
	 *
	 * Uses the substitution t² = b - x, which maps [a, b] → [0, √(b-a)]:
	 *   - x = b - t²
	 *   - dx = -2t dt
	 *   - √(b-x) = t
	 *   - ∫[a,b] g(x) dx = ∫[0,√(b-a)] g(b - t²) · 2t dt
	 *
	 * The factor of 2t cancels the 1/√(b-x) = 1/t singularity.
	 */
	class UpperSingularTransform : public IRealFunction
	{
	private:
		const IRealFunction& _func;
		Real _b;  // Singular point (upper bound)
		
	public:
		UpperSingularTransform(const IRealFunction& func, Real b) 
			: _func(func), _b(b) {}

		Real operator()(Real t) const override {
			// x = b - t², dx = -2t dt, but we integrate from 0 to √(b-a) so sign absorbed
			Real x = _b - t * t;
			return _func(x) * 2.0 * t;
		}
	};

	/**
	 * @brief Integrate with singularity at lower bound: ∫[a, b] f(x) dx where f has singularity at x = a
	 *
	 * Handles integrals where the integrand has an integrable singularity at the
	 * lower endpoint. The singularity is typically of the form 1/√(x-a) or
	 * 1/(x-a)^α where 0 < α < 1.
	 *
	 * MATHEMATICAL TRANSFORMATION:
	 * Uses the substitution t² = x - a:
	 *   ∫[a,b] f(x) dx = ∫[0,√(b-a)] f(a + t²) · 2t dt
	 *
	 * IMPORTANT: The function f should already include the singular behavior.
	 * The transformation removes the square-root singularity by the Jacobian 2t.
	 *
	 * EXAMPLES:
	 * - ∫[0,1] 1/√x dx = 2  (singularity at x=0)
	 * - ∫[0,1] x/√x dx = ∫[0,1] √x dx = 2/3
	 * - ∫[0,1] ln(x)/√x dx = -4 (double singularity)
	 *
	 * @param func Function to integrate (may be singular at x = a)
	 * @param a Lower bound with singularity
	 * @param b Upper bound (finite, no singularity)
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫[0,1] 1/√x dx = 2
	 *   RealFunction f{ [](Real x) { return 1.0/std::sqrt(x); } };
	 *   auto result = IntegrateLowerSingular(f, 0.0, 1.0);
	 *   // result.value ≈ 2.0
	 */
	static IntegrationResult IntegrateLowerSingular(const IRealFunction& func, Real a, Real b,
	                                                const Real eps = Defaults::RombergIntegrationEPS)
	{
		LowerSingularTransform transformed(func, a);
		Real upper = std::sqrt(b - a);
		// Use open integration (midpoint) to avoid endpoint issues
		return IntegrateOpen(transformed, 0.0, upper, eps);
	}

	/**
	 * @brief Integrate with singularity at upper bound: ∫[a, b] f(x) dx where f has singularity at x = b
	 *
	 * Handles integrals where the integrand has an integrable singularity at the
	 * upper endpoint.
	 *
	 * MATHEMATICAL TRANSFORMATION:
	 * Uses the substitution t² = b - x:
	 *   ∫[a,b] f(x) dx = ∫[0,√(b-a)] f(b - t²) · 2t dt
	 *
	 * @param func Function to integrate (may be singular at x = b)
	 * @param a Lower bound (finite, no singularity)
	 * @param b Upper bound with singularity
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫[0,1] 1/√(1-x) dx = 2
	 *   RealFunction f{ [](Real x) { return 1.0/std::sqrt(1.0 - x); } };
	 *   auto result = IntegrateUpperSingular(f, 0.0, 1.0);
	 *   // result.value ≈ 2.0
	 */
	static IntegrationResult IntegrateUpperSingular(const IRealFunction& func, Real a, Real b,
	                                                const Real eps = Defaults::RombergIntegrationEPS)
	{
		UpperSingularTransform transformed(func, b);
		Real upper = std::sqrt(b - a);
		return IntegrateOpen(transformed, 0.0, upper, eps);
	}

	/**
	 * @brief Integrate with singularities at both endpoints
	 *
	 * Handles integrals where the integrand has integrable singularities at both
	 * endpoints a and b. The integral is split at the midpoint.
	 *
	 * ALGORITHM:
	 * ∫[a,b] f(x) dx = ∫[a,m] f(x) dx + ∫[m,b] f(x) dx
	 * where m = (a+b)/2, and each half is handled with the appropriate
	 * endpoint singularity transformation.
	 *
	 * @param func Function to integrate (may be singular at x = a and x = b)
	 * @param a Lower bound with singularity
	 * @param b Upper bound with singularity
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫[0,1] 1/√(x(1-x)) dx = π  (Beta function B(1/2, 1/2))
	 *   RealFunction f{ [](Real x) { return 1.0/std::sqrt(x*(1.0-x)); } };
	 *   auto result = IntegrateBothSingular(f, 0.0, 1.0);
	 *   // result.value ≈ π
	 */
	static IntegrationResult IntegrateBothSingular(const IRealFunction& func, Real a, Real b,
	                                               const Real eps = Defaults::RombergIntegrationEPS)
	{
		Real mid = 0.5 * (a + b);
		auto left = IntegrateLowerSingular(func, a, mid, eps);
		auto right = IntegrateUpperSingular(func, mid, b, eps);
		
		return IntegrationResult(
			left.value + right.value,
			left.error_estimate + right.error_estimate,
			left.iterations + right.iterations,
			left.converged && right.converged
		);
	}

	/**
	 * @brief Integrate with interior singularity: ∫[a, b] f(x) dx where f has singularity at x = c ∈ (a, b)
	 *
	 * Handles integrals where the integrand has an integrable singularity at an
	 * interior point c. The integral is split at c and each part is handled
	 * with the appropriate endpoint transformation.
	 *
	 * ALGORITHM:
	 * ∫[a,b] f(x) dx = ∫[a,c] f(x) dx + ∫[c,b] f(x) dx
	 * where the first integral has an upper singularity at c, and the second
	 * has a lower singularity at c.
	 *
	 * @param func Function to integrate (may be singular at x = c)
	 * @param a Lower bound (finite, no singularity)
	 * @param b Upper bound (finite, no singularity)
	 * @param c Interior singular point (a < c < b)
	 * @param eps Desired relative accuracy
	 * @return IntegrationResult with value, error estimate, and convergence status
	 *
	 * @example
	 *   // ∫[-1,1] 1/√|x| dx = 4  (singularity at x=0)
	 *   RealFunction f{ [](Real x) { return 1.0/std::sqrt(std::abs(x)); } };
	 *   auto result = IntegrateInteriorSingular(f, -1.0, 1.0, 0.0);
	 *   // result.value ≈ 4.0
	 */
	static IntegrationResult IntegrateInteriorSingular(const IRealFunction& func, Real a, Real b, Real c,
	                                                   const Real eps = Defaults::RombergIntegrationEPS)
	{
		if (c <= a || c >= b)
			throw std::invalid_argument("IntegrateInteriorSingular: singular point c must be in (a, b)");
		
		auto left = IntegrateUpperSingular(func, a, c, eps);
		auto right = IntegrateLowerSingular(func, c, b, eps);
		
		return IntegrationResult(
			left.value + right.value,
			left.error_estimate + right.error_estimate,
			left.iterations + right.iterations,
			left.converged && right.converged
		);
	}

	/*************************************************************************/
	/*****              CONVENIENCE ALIASES & CONSTANTS                  *****/
	/*************************************************************************/

	/// Mathematical infinity constant for readable function calls
	constexpr Real PosInfinity = std::numeric_limits<Real>::infinity();
	constexpr Real NegInfinity = -std::numeric_limits<Real>::infinity();

} // namespace MML

#endif // MML_INTEGRATION_IMPROPER_H