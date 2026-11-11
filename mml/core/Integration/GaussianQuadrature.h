///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GaussianQuadrature.h                                                ///
///  Description: N-point Gaussian quadrature (Legendre, Laguerre, Hermite, Jacobi)   ///
///               Computed nodes and weights for high-accuracy integration            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_GAUSSIAN_QUADRATURE_H
#define MML_GAUSSIAN_QUADRATURE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <mml/MMLBase.h>
#include <mml/base/StandardFunctions.h>
#include "IntegrationBase.h"

#include <mml/interfaces/IFunction.h>

namespace MML
{
	/**
	 * @brief Result structure for Gaussian quadrature computation
	 */
	struct GaussQuadratureRule {
		static size_t CheckedOrder(int order) {
			if (order <= 0) throw ArgumentError("Gaussian quadrature order must be positive");
			return static_cast<size_t>(order);
		}
		std::vector<Real> nodes;    ///< Quadrature nodes (abscissas)
		std::vector<Real> weights;  ///< Quadrature weights
		int n;                      ///< Number of points
		
		GaussQuadratureRule(int numPoints) : nodes(CheckedOrder(numPoints)), weights(numPoints), n(numPoints) {}
	};

	/*************************************************************************/
	/*****                   GAUSS-LEGENDRE                              *****/
	/*************************************************************************/

	/**
	 * @brief Compute Gauss-Legendre quadrature nodes and weights
	 *
	 * Gauss-Legendre quadrature integrates with weight function w(x) = 1 on [-1, 1].
	 * With n points, it exactly integrates polynomials of degree ≤ 2n-1.
	 *
	 * The nodes are roots of the Legendre polynomial Pₙ(x), and the weights are
	 * computed from the derivative of Pₙ at each node.
	 *
	 * ALGORITHM:
	 * 1. Use asymptotic formula for initial node estimates
	 * 2. Newton-Raphson iteration to find exact roots of Pₙ(x)
	 * 3. Compute weights from the derivative Pₙ'(xᵢ)
	 *
	 * @param[out] x Vector of nodes (size determines n)
	 * @param[out] w Vector of weights (same size as x)
	 * @param x1 Lower bound of integration interval (default -1)
	 * @param x2 Upper bound of integration interval (default +1)
	 *
	 * @note For interval [a,b], nodes and weights are automatically transformed.
	 */
	static void GaussLegendre(std::vector<Real>& x, std::vector<Real>& w,
	                          Real x1 = -1.0, Real x2 = 1.0)
	{
		const Real EPS = std::numeric_limits<Real>::epsilon();
		const int MAXIT = 100;
		const Real PI = Constants::PI;
		
		int n = static_cast<int>(x.size());
		if (w.size() != x.size())
			throw VectorDimensionError("GaussLegendre: x and w must have same size", static_cast<int>(x.size()), static_cast<int>(w.size()));
		
		int m = (n + 1) / 2;  // Exploit symmetry
		Real xl = IntegrationDetail::HalfWidth(x1, x2);
		Real xm = std::midpoint(x1, x2);
		
		for (int i = 0; i < m; i++)
		{
			// Initial approximation for i-th root using asymptotic formula
			Real z = std::cos(PI * (i + MML::Real(0.75)) / (n + MML::Real(0.5)));
			
			// Newton-Raphson iteration
			int its = MAXIT;
			Real z1, pp;
			do {
				Real p1 = 1.0;
				Real p2 = 0.0;
				
				// Evaluate Legendre polynomial Pₙ(z) using recurrence
				for (int j = 0; j < n; j++)
				{
					Real p3 = p2;
					p2 = p1;
					p1 = ((MML::Real(2.0) * j + MML::Real(1.0)) * z * p2 - j * p3) / (j + 1);
				}
				
				// Derivative: Pₙ'(z) = n(zPₙ - Pₙ₋₁)/(z² - 1)
				pp = n * (z * p1 - p2) / (z * z - MML::Real(1.0));
				
				z1 = z;
				z = z1 - p1 / pp;  // Newton update
				if (--its <= 0) break;
			} while (std::abs(z - z1) > EPS);
			
			// Store symmetric nodes and weights
			x[i] = xm - xl * z;
			x[n - 1 - i] = xm + xl * z;
			w[i] = IntegrationDetail::Product({2, xl, Real(1) / ((1 - z * z) * pp * pp)});
			w[n - 1 - i] = w[i];
		}
	}

	/**
	 * @brief Compute n-point Gauss-Legendre rule and return as structure
	 */
	static GaussQuadratureRule GaussLegendreRule(int n, Real a = -1.0, Real b = 1.0)
	{
		GaussQuadratureRule rule(n);
		GaussLegendre(rule.nodes, rule.weights, a, b);
		return rule;
	}

	/*************************************************************************/
	/*****                   GAUSS-LAGUERRE                              *****/
	/*************************************************************************/

	/**
	 * @brief Compute Gauss-Laguerre quadrature nodes and weights
	 *
	 * Gauss-Laguerre quadrature integrates with weight function w(x) = x^α e^(-x) on [0, ∞).
	 * This is ideal for integrals of the form:
	 *   ∫₀^∞ x^α e^(-x) f(x) dx ≈ Σ wᵢ f(xᵢ)
	 *
	 * The nodes are roots of the generalized Laguerre polynomial Lₙ^(α)(x).
	 *
	 * COMMON VALUES OF α:
	 * - α = 0: Standard Laguerre (most common)
	 * - α = -0.5: For ∫ f(x)/√x · e^(-x) dx
	 * - α = 0.5: For ∫ √x · f(x) · e^(-x) dx
	 *
	 * @param[out] x Vector of nodes (size determines n)
	 * @param[out] w Vector of weights (same size as x)
	 * @param alpha Exponent parameter α (default 0)
	 */
	static void GaussLaguerre(std::vector<Real>& x, std::vector<Real>& w, Real alpha = 0.0)
	{
		const int MAXIT = 40;
		// Use relaxed convergence: machine epsilon is too tight for higher-order
		// polynomial evaluation which loses precision through the recurrence relation
		const Real EPS = std::max(std::numeric_limits<Real>::epsilon() * 100, 
		                          static_cast<Real>(1e-15));
		
		int n = static_cast<int>(x.size());
		if (w.size() != x.size())
			throw VectorDimensionError("GaussLaguerre: x and w must have same size", static_cast<int>(x.size()), static_cast<int>(w.size()));
		
		for (int i = 0; i < n; i++)
		{
			// Initial approximation for roots
			Real z;
			if (i == 0) {
				z = (MML::Real(1.0) + alpha) * (MML::Real(3.0) + MML::Real(0.92) * alpha) / (MML::Real(1.0) + MML::Real(2.4) * n + MML::Real(1.8) * alpha);
			} else if (i == 1) {
				z = x[0] + (MML::Real(15.0) + MML::Real(6.25) * alpha) / (MML::Real(1.0) + MML::Real(0.9) * alpha + MML::Real(2.5) * n);
			} else {
				Real ai = Real(i - 1);
				z = x[i-1] + ((MML::Real(1.0) + MML::Real(2.55) * ai) / (MML::Real(1.9) * ai) + MML::Real(1.26) * ai * alpha /
				              (MML::Real(1.0) + MML::Real(3.5) * ai)) * (x[i-1] - x[i-2]) / (MML::Real(1.0) + MML::Real(0.3) * alpha);
			}
			
			// Newton-Raphson iteration
			int its;
			Real pp, p2;
			for (its = 0; its < MAXIT; its++)
			{
				Real p1 = 1.0;
				p2 = 0.0;
				
				// Evaluate Laguerre polynomial using recurrence
				for (int j = 0; j < n; j++)
				{
					Real p3 = p2;
					p2 = p1;
					p1 = ((2 * j + 1 + alpha - z) * p2 - (j + alpha) * p3) / (j + 1);
				}
				
				// Derivative
				pp = (n * p1 - (n + alpha) * p2) / z;
				
				Real z1 = z;
				z = z1 - p1 / pp;
					if (std::abs(z - z1) <= EPS * (1.0 + std::abs(z))) break;
			}
			if (its >= MAXIT)
				throw IntegrationTooManySteps("GaussLaguerre: too many iterations", MAXIT);
			
			x[i] = z;
			// Weight formula using gamma functions
			w[i] = -std::exp(Functions::LGamma(alpha + n) - Functions::LGamma(static_cast<Real>(n))) 
			       / (pp * n * p2);
		}
	}

	/**
	 * @brief Compute n-point Gauss-Laguerre rule and return as structure
	 */
	static GaussQuadratureRule GaussLaguerreRule(int n, Real alpha = 0.0)
	{
		GaussQuadratureRule rule(n);
		GaussLaguerre(rule.nodes, rule.weights, alpha);
		return rule;
	}

	/*************************************************************************/
	/*****                    GAUSS-HERMITE                              *****/
	/*************************************************************************/

	/**
	 * @brief Compute Gauss-Hermite quadrature nodes and weights
	 *
	 * Gauss-Hermite quadrature integrates with weight function w(x) = e^(-x²) on (-∞, ∞).
	 * This is ideal for integrals of the form:
	 *   ∫_{-∞}^∞ e^(-x²) f(x) dx ≈ Σ wᵢ f(xᵢ)
	 *
	 * The nodes are roots of the Hermite polynomial Hₙ(x).
	 *
	 * APPLICATIONS:
	 * - Quantum mechanics (harmonic oscillator)
	 * - Probability (Gaussian integrals)
	 * - Statistical mechanics
	 *
	 * @param[out] x Vector of nodes (size determines n)
	 * @param[out] w Vector of weights (same size as x)
	 */
	static void GaussHermite(std::vector<Real>& x, std::vector<Real>& w)
	{
		const Real EPS = std::numeric_limits<Real>::epsilon();
		const Real PIM4 = MML::Real(0.7511255444649425);  // π^(-1/4)
		const int MAXIT = 15;  // Increased for edge cases
		
		int n = static_cast<int>(x.size());
		if (w.size() != x.size())
			throw VectorDimensionError("GaussHermite: x and w must have same size", static_cast<int>(x.size()), static_cast<int>(w.size()));
		
		int m = (n + 1) / 2;  // Exploit symmetry
		Real z = 0.0;  // Moved outside loop - needed for initial guess formulas
		
		for (int i = 0; i < m; i++)
		{
			// Initial approximation for roots (z preserves value from previous iteration)
			if (i == 0) {
				z = std::sqrt(static_cast<Real>(2 * n + 1)) - 
				    MML::Real(1.85575) * std::pow(static_cast<Real>(2 * n + 1), -MML::Real(0.16667));
			} else if (i == 1) {
				z = z - MML::Real(1.14) * std::pow(static_cast<Real>(n), MML::Real(0.426)) / z;
			} else if (i == 2) {
				z = MML::Real(1.86) * z - MML::Real(0.86) * x[0];
			} else if (i == 3) {
				z = MML::Real(1.91) * z - MML::Real(0.91) * x[1];
			} else {
				z = MML::Real(2.0) * z - x[i - 2];
			}
			
			// Newton-Raphson iteration
			int its;
			Real pp;
			for (its = 0; its < MAXIT; its++)
			{
				Real p1 = PIM4;
				Real p2 = 0.0;
				
				// Evaluate Hermite polynomial using recurrence (normalized form)
				for (int j = 0; j < n; j++)
				{
					Real p3 = p2;
					p2 = p1;
					p1 = z * std::sqrt(MML::Real(2.0) / (j + 1)) * p2 -
					     std::sqrt(static_cast<Real>(j) / (j + 1)) * p3;
				}
				
				// Derivative
				pp = std::sqrt(static_cast<Real>(2 * n)) * p2;
				
				Real z1 = z;
				z = z1 - p1 / pp;
				if (std::abs(z - z1) <= EPS * (1.0 + std::abs(z))) break;
			}
			if (its >= MAXIT)
				throw IntegrationTooManySteps("GaussHermite: too many iterations", MAXIT);
			
			// Store symmetric nodes and weights
			x[i] = z;
			x[n - 1 - i] = -z;
			w[i] = MML::Real(2.0) / (pp * pp);
			w[n - 1 - i] = w[i];
		}
	}

	/**
	 * @brief Compute n-point Gauss-Hermite rule and return as structure
	 */
	static GaussQuadratureRule GaussHermiteRule(int n)
	{
		GaussQuadratureRule rule(n);
		GaussHermite(rule.nodes, rule.weights);
		return rule;
	}

	/*************************************************************************/
	/*****                    GAUSS-JACOBI                               *****/
	/*************************************************************************/

	/**
	 * @brief Compute Gauss-Jacobi quadrature nodes and weights
	 *
	 * Gauss-Jacobi quadrature integrates with weight function 
	 * w(x) = (1-x)^α (1+x)^β on [-1, 1].
	 *
	 * This generalizes several other quadrature rules:
	 * - α = β = 0: Gauss-Legendre
	 * - α = β = -0.5: Gauss-Chebyshev of the first kind
	 * - α = β = 0.5: Gauss-Chebyshev of the second kind
	 *
	 * APPLICATIONS:
	 * - Integrals with endpoint singularities
	 * - Weighted polynomial approximation
	 * - Spectral methods
	 *
	 * @param[out] x Vector of nodes (size determines n)
	 * @param[out] w Vector of weights (same size as x)
	 * @param alpha Parameter α (must be > -1)
	 * @param beta Parameter β (must be > -1)
	 */
	static void GaussJacobi(std::vector<Real>& x, std::vector<Real>& w, 
	                        Real alpha, Real beta)
	{
		const int MAXIT = 15;
		const Real EPS = std::numeric_limits<Real>::epsilon();
		
		int n = static_cast<int>(x.size());
		if (w.size() != x.size())
			throw VectorDimensionError("GaussJacobi: x and w must have same size", static_cast<int>(x.size()), static_cast<int>(w.size()));
		if (alpha <= -1.0 || beta <= -1.0)
			throw NumericalMethodError("GaussJacobi: alpha and beta must be > -1");
		
		Real alfbet = alpha + beta;
		Real z = 0.0;  // Moved outside loop - needed for initial guess formulas
		
		for (int i = 0; i < n; i++)
		{
			// Initial approximation for roots (z preserves value from previous iteration)
			if (i == 0) {
				Real an = alpha / n;
				Real bn = beta / n;
				Real r1 = (MML::Real(1.0) + alpha) * (MML::Real(2.78) / (MML::Real(4.0) + n * n) + MML::Real(0.768) * an / n);
				Real r2 = MML::Real(1.0) + MML::Real(1.48) * an + MML::Real(0.96) * bn + MML::Real(0.452) * an * an + MML::Real(0.83) * an * bn;
				z = MML::Real(1.0) - r1 / r2;
			} else if (i == 1) {
				Real r1 = (MML::Real(4.1) + alpha) / ((MML::Real(1.0) + alpha) * (MML::Real(1.0) + MML::Real(0.156) * alpha));
				Real r2 = MML::Real(1.0) + MML::Real(0.06) * (n - MML::Real(8.0)) * (MML::Real(1.0) + MML::Real(0.12) * alpha) / n;
				Real r3 = MML::Real(1.0) + MML::Real(0.012) * beta * (MML::Real(1.0) + MML::Real(0.25) * std::abs(alpha)) / n;
				z = z - (MML::Real(1.0) - z) * r1 * r2 * r3;
			} else if (i == 2) {
				Real r1 = (MML::Real(1.67) + MML::Real(0.28) * alpha) / (MML::Real(1.0) + MML::Real(0.37) * alpha);
				Real r2 = MML::Real(1.0) + MML::Real(0.22) * (n - MML::Real(8.0)) / n;
				Real r3 = MML::Real(1.0) + MML::Real(8.0) * beta / ((MML::Real(6.28) + beta) * n * n);
				z = z - (x[0] - z) * r1 * r2 * r3;
			} else if (i == n - 2) {
				Real r1 = (MML::Real(1.0) + MML::Real(0.235) * beta) / (MML::Real(0.766) + MML::Real(0.119) * beta);
				Real r2 = MML::Real(1.0) / (MML::Real(1.0) + MML::Real(0.639) * (n - MML::Real(4.0)) / (MML::Real(1.0) + MML::Real(0.71) * (n - MML::Real(4.0))));
				Real r3 = MML::Real(1.0) / (MML::Real(1.0) + MML::Real(20.0) * alpha / ((MML::Real(7.5) + alpha) * n * n));
				z = z + (z - x[n - 4]) * r1 * r2 * r3;
			} else if (i == n - 1) {
				Real r1 = (MML::Real(1.0) + MML::Real(0.37) * beta) / (MML::Real(1.67) + MML::Real(0.28) * beta);
				Real r2 = MML::Real(1.0) / (MML::Real(1.0) + MML::Real(0.22) * (n - MML::Real(8.0)) / n);
				Real r3 = MML::Real(1.0) / (MML::Real(1.0) + MML::Real(8.0) * alpha / ((MML::Real(6.28) + alpha) * n * n));
				z = z + (z - x[n - 3]) * r1 * r2 * r3;
			} else {
				z = MML::Real(3.0) * x[i - 1] - MML::Real(3.0) * x[i - 2] + x[i - 3];
			}
			
			// Newton-Raphson iteration
			int its;
			Real pp, p2;
			for (its = 0; its < MAXIT; its++)
			{
				Real temp = MML::Real(2.0) + alfbet;
				Real p1 = (alpha - beta + temp * z) / MML::Real(2.0);
				p2 = 1.0;
				
				// Evaluate Jacobi polynomial using recurrence
				for (int j = 2; j <= n; j++)
				{
					Real p3 = p2;
					p2 = p1;
					temp = 2 * j + alfbet;
					Real a = 2 * j * (j + alfbet) * (temp - MML::Real(2.0));
					Real b = (temp - MML::Real(1.0)) * (alpha * alpha - beta * beta + temp * (temp - MML::Real(2.0)) * z);
					Real c = MML::Real(2.0) * (j - 1 + alpha) * (j - 1 + beta) * temp;
					p1 = (b * p2 - c * p3) / a;
				}
				
				// Derivative
				Real temp2 = 2 * n + alfbet;
				pp = (n * (alpha - beta - temp2 * z) * p1 + 
				      MML::Real(2.0) * (n + alpha) * (n + beta) * p2) / (temp2 * (MML::Real(1.0) - z * z));
				
				Real z1 = z;
				z = z1 - p1 / pp;
				if (std::abs(z - z1) <= EPS) break;
			}
			if (its >= MAXIT)
				throw IntegrationTooManySteps("GaussJacobi: too many iterations", MAXIT);
			
			x[i] = z;
			// Weight formula using gamma functions - computed in log domain to avoid overflow
			// For extreme parameters (large n, α, β), exp(LGamma(...)) can overflow
			Real temp = 2 * n + alfbet;
			Real log_weight = Functions::LGamma(alpha + n) + Functions::LGamma(beta + n) -
			                  Functions::LGamma(n + MML::Real(1.0)) - Functions::LGamma(n + alfbet + MML::Real(1.0)) +
			                  alfbet * std::log(MML::Real(2.0));
			
			// Check for overflow before calling exp() (exp(709) ≈ 1.7e308, near DBL_MAX)
			if (log_weight > 700.0) {
				throw NumericalMethodError("GaussJacobi: weight computation overflow for given parameters (log_weight=" + std::to_string(log_weight) + ")");
			}
			w[i] = std::exp(log_weight) * temp / (pp * p2);
		}
	}

	/**
	 * @brief Compute n-point Gauss-Jacobi rule and return as structure
	 */
	static GaussQuadratureRule GaussJacobiRule(int n, Real alpha, Real beta)
	{
		GaussQuadratureRule rule(n);
		GaussJacobi(rule.nodes, rule.weights, alpha, beta);
		return rule;
	}

	/*************************************************************************/
	/*****              GAUSS-CHEBYSHEV (Special Cases)                  *****/
	/*************************************************************************/

	/**
	 * @brief Compute Gauss-Chebyshev quadrature of the first kind
	 *
	 * Weight function: w(x) = 1/√(1-x²) on [-1, 1]
	 * This is Gauss-Jacobi with α = β = -1/2
	 *
	 * The nodes and weights have closed-form expressions:
	 *   xᵢ = cos((2i+1)π/(2n))
	 *   wᵢ = π/n
	 */
	static void GaussChebyshev1(std::vector<Real>& x, std::vector<Real>& w)
	{
		int n = static_cast<int>(x.size());
		if (w.size() != x.size())
			throw VectorDimensionError("GaussChebyshev1: x and w must have same size", static_cast<int>(x.size()), static_cast<int>(w.size()));
		
		Real wgt = Constants::PI / n;
		for (int i = 0; i < n; i++)
		{
			x[i] = std::cos(Constants::PI * (2 * i + 1) / (MML::Real(2.0) * n));
			w[i] = wgt;
		}
	}

	/**
	 * @brief Compute Gauss-Chebyshev quadrature of the second kind
	 *
	 * Weight function: w(x) = √(1-x²) on [-1, 1]
	 * This is Gauss-Jacobi with α = β = 1/2
	 *
	 * The nodes and weights have closed-form expressions:
	 *   xᵢ = cos((i+1)π/(n+1))
	 *   wᵢ = π/(n+1) · sin²((i+1)π/(n+1))
	 */
	static void GaussChebyshev2(std::vector<Real>& x, std::vector<Real>& w)
	{
		int n = static_cast<int>(x.size());
		if (w.size() != x.size())
			throw VectorDimensionError("GaussChebyshev2: x and w must have same size", static_cast<int>(x.size()), static_cast<int>(w.size()));
		
		for (int i = 0; i < n; i++)
		{
			Real theta = Constants::PI * (i + 1) / (n + 1);
			x[i] = std::cos(theta);
			Real s = std::sin(theta);
			w[i] = Constants::PI / (n + 1) * s * s;
		}
	}

	/**
	 * @brief Compute n-point Gauss-Chebyshev Type 1 rule and return as structure
	 */
	static GaussQuadratureRule GaussChebyshev1Rule(int n)
	{
		GaussQuadratureRule rule(n);
		GaussChebyshev1(rule.nodes, rule.weights);
		return rule;
	}

	/**
	 * @brief Compute n-point Gauss-Chebyshev Type 2 rule and return as structure
	 */
	static GaussQuadratureRule GaussChebyshev2Rule(int n)
	{
		GaussQuadratureRule rule(n);
		GaussChebyshev2(rule.nodes, rule.weights);
		return rule;
	}

	/*************************************************************************/
	/*****              INTEGRATION USING COMPUTED RULES                 *****/
	/*************************************************************************/

	/**
	 * @brief Integrate using a pre-computed Gaussian quadrature rule
	 *
	 * @param func Function to integrate
	 * @param rule Pre-computed quadrature rule (nodes and weights)
	 * @return Approximate integral value
	 */
	static Real IntegrateWithRuleImpl(const IRealFunction& func, const GaussQuadratureRule& rule,
		IntegrationResult* progress, bool checkFinite)
	{
		if (rule.n <= 0 || rule.nodes.size() != static_cast<size_t>(rule.n) || rule.weights.size() != rule.nodes.size())
			throw ArgumentError("Gaussian quadrature: inconsistent rule storage");
		for (int i = 0; i < rule.n; ++i)
			if (!std::isfinite(rule.nodes[i]) || !std::isfinite(rule.weights[i]))
				throw ArgumentError("Gaussian quadrature: rule nodes and weights must be finite");
		IntegrationDetail::Sum sum;
		int evaluations = 0;
		auto counted = [&](Real point) {
			if (progress) progress->function_evaluations = evaluations;
			return func(point);
		};
		for (int i = 0; i < rule.n; i++)
		{
			sum.Add(IntegrationDetail::Product({rule.weights[i],
				IntegrationDetail::Sample(counted, evaluations, checkFinite, rule.nodes[i])}));
		}
		IntegrationDetail::RequireFinite(sum.Value(), "Gaussian quadrature");
		return sum.Value();
	}

	static Real IntegrateWithRule(const IRealFunction& func, const GaussQuadratureRule& rule) {
		return IntegrateWithRuleImpl(func, rule, nullptr, true);
	}

	/**
	 * @brief Integrate f(x) on [a,b] using n-point Gauss-Legendre
	 *
	 * Convenience function that computes the rule and evaluates the integral.
	 *
	 * @param func Function to integrate
	 * @param a Lower bound
	 * @param b Upper bound
	 * @param n Number of quadrature points (default 20)
	 * @return IntegrationResult with value (error_estimate = 0 for non-adaptive)
	 */
	static IntegrationResult IntegrateGaussLegendreCore(const IRealFunction& func, 
	                                                 Real a, Real b, int n = 20)
	{
		auto rule = GaussLegendreRule(n, a, b);
		Real value = IntegrateWithRule(func, rule);
		return IntegrationResult::FixedRule(value, rule.n);
	}

	static inline Real IntegrateGaussLegendre(const IRealFunction& func,
	                                                 Real a, Real b, int n = 20)
	{
		return IntegrateGaussLegendreCore(func, a, b, n).ValueOrThrow("IntegrateGaussLegendre");
	}

	/// @brief Detailed variant with AlgorithmStatus, timing, and diagnostics
	static inline IntegrationDetailedResult IntegrateGaussLegendreDetailed(const IRealFunction& func,
	                                                                 Real a, Real b, int n = 20,
	                                                                 const IntegrationConfig& config = IntegrationConfig())
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"GaussLegendre", config, [&](IntegrationDetailedResult& r) {
				auto rule = GaussLegendreRule(n, a, b);
				r.value = IntegrateWithRuleImpl(func, rule, &r, config.check_finite);
				r.MarkFixedRule(n);
			});
	}

	/**
	 * @brief Integrate f(x)·x^α·e^(-x) on [0,∞) using n-point Gauss-Laguerre
	 *
	 * Note: The weight function x^α·e^(-x) is built into the quadrature.
	 * You supply f(x) without the weight.
	 *
	 * @param func Function f(x) to integrate (not including weight)
	 * @param n Number of quadrature points (default 20)
	 * @param alpha Exponent parameter (default 0)
	 * @return IntegrationResult
	 */
	static IntegrationResult IntegrateGaussLaguerreCore(const IRealFunction& func, 
	                                                 int n = 20, Real alpha = 0.0)
	{
		auto rule = GaussLaguerreRule(n, alpha);
		Real value = IntegrateWithRule(func, rule);
		return IntegrationResult::FixedRule(value, rule.n);
	}

	static inline Real IntegrateGaussLaguerre(const IRealFunction& func,
	                                                 int n = 20, Real alpha = 0.0)
	{
		return IntegrateGaussLaguerreCore(func, n, alpha).ValueOrThrow("IntegrateGaussLaguerre");
	}

	/// @brief Detailed variant with AlgorithmStatus, timing, and diagnostics
	static inline IntegrationDetailedResult IntegrateGaussLaguerreDetailed(const IRealFunction& func,
	                                                                 int n = 20, Real alpha = 0.0,
	                                                                 const IntegrationConfig& config = IntegrationConfig())
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"GaussLaguerre", config, [&](IntegrationDetailedResult& r) {
				auto rule = GaussLaguerreRule(n, alpha);
				r.value = IntegrateWithRuleImpl(func, rule, &r, config.check_finite);
				r.MarkFixedRule(n);
			});
	}

	/**
	 * @brief Integrate f(x)·e^(-x²) on (-∞,∞) using n-point Gauss-Hermite
	 *
	 * Note: The weight function e^(-x²) is built into the quadrature.
	 * You supply f(x) without the weight.
	 *
	 * @param func Function f(x) to integrate (not including weight)
	 * @param n Number of quadrature points (default 20)
	 * @return IntegrationResult
	 */
	static IntegrationResult IntegrateGaussHermiteCore(const IRealFunction& func, int n = 20)
	{
		auto rule = GaussHermiteRule(n);
		Real value = IntegrateWithRule(func, rule);
		return IntegrationResult::FixedRule(value, rule.n);
	}

	static inline Real IntegrateGaussHermite(const IRealFunction& func, int n = 20)
	{
		return IntegrateGaussHermiteCore(func, n).ValueOrThrow("IntegrateGaussHermite");
	}

	/// @brief Detailed variant with AlgorithmStatus, timing, and diagnostics
	static inline IntegrationDetailedResult IntegrateGaussHermiteDetailed(const IRealFunction& func,
	                                                                int n = 20,
	                                                                const IntegrationConfig& config = IntegrationConfig())
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"GaussHermite", config, [&](IntegrationDetailedResult& r) {
				auto rule = GaussHermiteRule(n);
				r.value = IntegrateWithRuleImpl(func, rule, &r, config.check_finite);
				r.MarkFixedRule(n);
			});
	}

	/**
	 * @brief Integrate f(x)·(1-x)^α·(1+x)^β on [-1,1] using n-point Gauss-Jacobi
	 *
	 * Note: The weight function (1-x)^α·(1+x)^β is built into the quadrature.
	 * You supply f(x) without the weight.
	 *
	 * Special cases:
	 * - α = β = 0: Gauss-Legendre
	 * - α = β = -0.5: Gauss-Chebyshev of the first kind
	 * - α = β = 0.5: Gauss-Chebyshev of the second kind
	 *
	 * @param func Function f(x) to integrate (not including weight)
	 * @param alpha Parameter α (must be > -1)
	 * @param beta Parameter β (must be > -1)
	 * @param n Number of quadrature points (default 20)
	 * @return IntegrationResult
	 */
	static IntegrationResult IntegrateGaussJacobiCore(const IRealFunction& func, 
	                                               Real alpha, Real beta, int n = 20)
	{
		auto rule = GaussJacobiRule(n, alpha, beta);
		Real value = IntegrateWithRule(func, rule);
		return IntegrationResult::FixedRule(value, rule.n);
	}

	static inline Real IntegrateGaussJacobi(const IRealFunction& func,
	                                               Real alpha, Real beta, int n = 20)
	{
		return IntegrateGaussJacobiCore(func, alpha, beta, n).ValueOrThrow("IntegrateGaussJacobi");
	}

	/// @brief Detailed variant with AlgorithmStatus, timing, and diagnostics
	static inline IntegrationDetailedResult IntegrateGaussJacobiDetailed(const IRealFunction& func,
	                                                               Real alpha, Real beta, int n = 20,
	                                                               const IntegrationConfig& config = IntegrationConfig())
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"GaussJacobi", config, [&](IntegrationDetailedResult& r) {
				auto rule = GaussJacobiRule(n, alpha, beta);
				r.value = IntegrateWithRuleImpl(func, rule, &r, config.check_finite);
				r.MarkFixedRule(n);
			});
	}

	/**
	 * @brief Integrate f(x)/√(1-x²) on [-1,1] using n-point Gauss-Chebyshev (first kind)
	 *
	 * Note: The weight function 1/√(1-x²) is built into the quadrature.
	 * You supply f(x) without the weight.
	 *
	 * This is particularly efficient because weights and nodes have closed-form expressions.
	 *
	 * @param func Function f(x) to integrate (not including weight)
	 * @param n Number of quadrature points (default 20)
	 * @return IntegrationResult
	 */
	static IntegrationResult IntegrateGaussChebyshev1Core(const IRealFunction& func, int n = 20)
	{
		auto rule = GaussChebyshev1Rule(n);
		Real value = IntegrateWithRule(func, rule);
		return IntegrationResult::FixedRule(value, rule.n);
	}

	static inline Real IntegrateGaussChebyshev1(const IRealFunction& func, int n = 20)
	{
		return IntegrateGaussChebyshev1Core(func, n).ValueOrThrow("IntegrateGaussChebyshev1");
	}

	/// @brief Detailed variant with AlgorithmStatus, timing, and diagnostics
	static inline IntegrationDetailedResult IntegrateGaussChebyshev1Detailed(const IRealFunction& func,
	                                                                   int n = 20,
	                                                                   const IntegrationConfig& config = IntegrationConfig())
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"GaussChebyshev1", config, [&](IntegrationDetailedResult& r) {
				auto rule = GaussChebyshev1Rule(n);
				r.value = IntegrateWithRuleImpl(func, rule, &r, config.check_finite);
				r.MarkFixedRule(n);
			});
	}

	/**
	 * @brief Integrate f(x)·√(1-x²) on [-1,1] using n-point Gauss-Chebyshev (second kind)
	 *
	 * Note: The weight function √(1-x²) is built into the quadrature.
	 * You supply f(x) without the weight.
	 *
	 * This is particularly efficient because weights and nodes have closed-form expressions.
	 *
	 * @param func Function f(x) to integrate (not including weight)
	 * @param n Number of quadrature points (default 20)
	 * @return IntegrationResult
	 */
	static IntegrationResult IntegrateGaussChebyshev2Core(const IRealFunction& func, int n = 20)
	{
		auto rule = GaussChebyshev2Rule(n);
		Real value = IntegrateWithRule(func, rule);
		return IntegrationResult::FixedRule(value, rule.n);
	}

	static inline Real IntegrateGaussChebyshev2(const IRealFunction& func, int n = 20)
	{
		return IntegrateGaussChebyshev2Core(func, n).ValueOrThrow("IntegrateGaussChebyshev2");
	}

	/// @brief Detailed variant with AlgorithmStatus, timing, and diagnostics
	static inline IntegrationDetailedResult IntegrateGaussChebyshev2Detailed(const IRealFunction& func,
	                                                                   int n = 20,
	                                                                   const IntegrationConfig& config = IntegrationConfig())
	{
		return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
			"GaussChebyshev2", config, [&](IntegrationDetailedResult& r) {
				auto rule = GaussChebyshev2Rule(n);
				r.value = IntegrateWithRuleImpl(func, rule, &r, config.check_finite);
				r.MarkFixedRule(n);
			});
	}

} // namespace MML

#endif // MML_GAUSSIAN_QUADRATURE_H
