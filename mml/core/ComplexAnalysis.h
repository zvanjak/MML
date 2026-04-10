///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ComplexAnalysis.h                                                   ///
///  Description: Complex analysis tools: contours, contour integration,              ///
///               winding number, Cauchy integral, residues, argument principle        ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [AH79]   Ahlfors, L.V. (1979). Complex Analysis, 3rd ed. McGraw-Hill          ///
///    [SS03]   Stein & Shakarchi (2003). Complex Analysis, Princeton Univ. Press     ///
///    [NR3]    Press et al., Numerical Recipes 3rd ed.                               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMPLEX_ANALYSIS_H
#define MML_COMPLEX_ANALYSIS_H

#include "MMLBase.h"

#include "core/AlgorithmTypes.h"
#include "interfaces/IComplexFunction.h"
#include "base/ComplexFunction.h"
#include "core/Derivation/DerivationComplex.h"

#include <cmath>
#include <functional>

namespace MML
{
	namespace ComplexAnalysis
	{
		/********************************************************************************************************************/
		/********                                   Result Types                                                     ********/
		/********************************************************************************************************************/

		/// Result of a complex contour integration
		struct ContourIntegrationResult {
			Complex value = Complex(0.0, 0.0);   ///< Computed integral value
			Real error_estimate = 0.0;            ///< Estimated absolute error (on |value|)
			int function_evaluations = 0;         ///< Number of integrand evaluations
			bool converged = true;                ///< True if error < tolerance

			/// Implicit conversion to Complex for convenience
			operator Complex() const { return value; }
		};

		/********************************************************************************************************************/
		/********                                   Contour Classes                                                  ********/
		/********************************************************************************************************************/

		/// @brief Circle contour γ(t) = center + radius · e^(it), t ∈ [0, 2π]
		/// @details The standard contour for residue computation and Cauchy integrals.
		///          Traversed counter-clockwise (positive orientation).
		class CircleContour : public IRealToComplexFunction
		{
			Complex _center;
			Real _radius;
		public:
			CircleContour(Complex center, Real radius) : _center(center), _radius(radius) {}

			/// @brief γ(t) = center + radius · e^(it)
			Complex operator()(Real t) const override {
				return _center + _radius * Complex(std::cos(t), std::sin(t));
			}

			/// @brief γ'(t) = i · radius · e^(it)
			Complex derivative(Real t) const {
				return _radius * Complex(-std::sin(t), std::cos(t));
			}

			Complex center() const { return _center; }
			Real radius() const { return _radius; }
			Real t_start() const { return 0.0; }
			Real t_end() const { return 2.0 * Constants::PI; }
		};

		/// @brief Line segment contour γ(t) = (1-t)·z1 + t·z2, t ∈ [0, 1]
		class LineSegmentContour : public IRealToComplexFunction
		{
			Complex _z1, _z2;
		public:
			LineSegmentContour(Complex z1, Complex z2) : _z1(z1), _z2(z2) {}

			/// @brief γ(t) = (1-t)·z1 + t·z2
			Complex operator()(Real t) const override {
				return (REAL(1.0) - t) * _z1 + t * _z2;
			}

			/// @brief γ'(t) = z2 - z1 (constant)
			Complex derivative(Real /*t*/) const {
				return _z2 - _z1;
			}

			Real t_start() const { return 0.0; }
			Real t_end() const { return 1.0; }
		};

		/// @brief Circular arc contour γ(t) = center + radius · e^(it), t ∈ [t1, t2]
		class ArcContour : public IRealToComplexFunction
		{
			Complex _center;
			Real _radius;
			Real _t1, _t2;
		public:
			ArcContour(Complex center, Real radius, Real t1, Real t2)
				: _center(center), _radius(radius), _t1(t1), _t2(t2) {}

			/// @brief γ(t) = center + radius · e^(it)  for t ∈ [t1, t2]
			Complex operator()(Real t) const override {
				return _center + _radius * Complex(std::cos(t), std::sin(t));
			}

			Complex derivative(Real t) const {
				return _radius * Complex(-std::sin(t), std::cos(t));
			}

			Real t_start() const { return _t1; }
			Real t_end() const { return _t2; }
		};

		/********************************************************************************************************************/
		/********                          Complex Adaptive Simpson Integration                                       ********/
		/********************************************************************************************************************/

		namespace Detail
		{
			/// @brief Adaptive Simpson's rule for complex-valued integrands
			/// @details Recursively subdivides intervals until error estimate < tolerance.
			///          The integrand g: R → C maps parameter t to the complex value f(γ(t))·γ'(t).
			static Complex AdaptiveSimpsonStep(
				const std::function<Complex(Real)>& g,
				Real a, Real b,
				Complex fa, Complex fm, Complex fb,
				Complex whole,
				Real tol, int depth, int max_depth,
				int& evals, bool& converged)
			{
				Real m = (a + b) / 2.0;
				Real m1 = (a + m) / 2.0;
				Real m2 = (m + b) / 2.0;

				Complex f_m1 = g(m1);
				Complex f_m2 = g(m2);
				evals += 2;

				Real h6 = (b - a) / 6.0;
				Real h12 = h6 / 2.0;

				Complex left = h12 * (fa + REAL(4.0) * f_m1 + fm);
				Complex right = h12 * (fm + REAL(4.0) * f_m2 + fb);
				Complex refined = left + right;

				Real error = std::abs(refined - whole) / 15.0;  // Richardson error estimate

				if (depth >= max_depth)
				{
					converged = false;
					return refined;
				}

				if (error < tol)
				{
					// Richardson extrapolation: add the correction term
					return refined + (refined - whole) / REAL(15.0);
				}

				// Recurse on both halves with halved tolerance
				return AdaptiveSimpsonStep(g, a, m, fa, f_m1, fm, left, tol / 2, depth + 1, max_depth, evals, converged)
				     + AdaptiveSimpsonStep(g, m, b, fm, f_m2, fb, right, tol / 2, depth + 1, max_depth, evals, converged);
			}

			/// @brief Top-level adaptive Simpson for complex-valued functions of a real variable
			static ContourIntegrationResult AdaptiveSimpsonComplex(
				const std::function<Complex(Real)>& g,
				Real a, Real b,
				Real tol = Precision::ComplexAnalysisTolerance,
				int max_depth = 30)
			{
				ContourIntegrationResult result;

				Complex fa = g(a);
				Complex fm = g((a + b) / 2.0);
				Complex fb = g(b);
				int evals = 3;
				bool converged = true;

				Real h6 = (b - a) / 6.0;
				Complex whole = h6 * (fa + REAL(4.0) * fm + fb);

				result.value = AdaptiveSimpsonStep(g, a, b, fa, fm, fb, whole, tol, 0, max_depth, evals, converged);
				result.function_evaluations = evals;
				result.converged = converged;
				// Rough error estimate: if converged, error < tol
				result.error_estimate = converged ? tol : std::abs(result.value) * 1e-6;
				return result;
			}
		}

		/********************************************************************************************************************/
		/********                               Contour Integration                                                  ********/
		/********************************************************************************************************************/

		/// @brief Compute ∮_γ f(z) dz along a circle contour
		/// @details Reduces to ∫₀²π f(γ(t)) · γ'(t) dt using adaptive Simpson
		/// @param f Complex function f:C→C
		/// @param contour CircleContour defining the path
		/// @param tol Absolute tolerance for the integration (default: Precision::ComplexAnalysisTolerance)
		/// @return ContourIntegrationResult with complex integral value
		static ContourIntegrationResult ContourIntegral(
			const IComplexFunction& f,
			const CircleContour& contour,
			Real tol = Precision::ComplexAnalysisTolerance)
		{
			auto integrand = [&](Real t) -> Complex {
				return f(contour(t)) * contour.derivative(t);
			};
			return Detail::AdaptiveSimpsonComplex(integrand, contour.t_start(), contour.t_end(), tol);
		}

		/// @brief Compute ∫_γ f(z) dz along a line segment contour
		static ContourIntegrationResult ContourIntegral(
			const IComplexFunction& f,
			const LineSegmentContour& contour,
			Real tol = Precision::ComplexAnalysisTolerance)
		{
			auto integrand = [&](Real t) -> Complex {
				return f(contour(t)) * contour.derivative(t);
			};
			return Detail::AdaptiveSimpsonComplex(integrand, contour.t_start(), contour.t_end(), tol);
		}

		/// @brief Compute ∫_γ f(z) dz along an arc contour
		static ContourIntegrationResult ContourIntegral(
			const IComplexFunction& f,
			const ArcContour& contour,
			Real tol = Precision::ComplexAnalysisTolerance)
		{
			auto integrand = [&](Real t) -> Complex {
				return f(contour(t)) * contour.derivative(t);
			};
			return Detail::AdaptiveSimpsonComplex(integrand, contour.t_start(), contour.t_end(), tol);
		}

		/// @brief Generic contour integral with user-provided parameterization and derivative
		/// @param f Complex function f:C→C
		/// @param gamma Contour parameterization γ: R → C
		/// @param gamma_deriv Derivative γ': R → C
		/// @param t_start Start of parameter range
		/// @param t_end End of parameter range
		/// @param tol Tolerance (default: Precision::ComplexAnalysisTolerance)
		static ContourIntegrationResult ContourIntegral(
			const IComplexFunction& f,
			const std::function<Complex(Real)>& gamma,
			const std::function<Complex(Real)>& gamma_deriv,
			Real t_start, Real t_end,
			Real tol = Precision::ComplexAnalysisTolerance)
		{
			auto integrand = [&](Real t) -> Complex {
				return f(gamma(t)) * gamma_deriv(t);
			};
			return Detail::AdaptiveSimpsonComplex(integrand, t_start, t_end, tol);
		}

		/********************************************************************************************************************/
		/********                               Winding Number                                                       ********/
		/********************************************************************************************************************/

		/// @brief Compute the winding number of contour γ around point z0
		/// @details n(γ, z0) = (1/2πi) ∮_γ dz / (z - z0)
		///          For a circle contour containing z0, this equals 1.
		///          For a circle not containing z0, this equals 0.
		/// @param contour Circle contour
		/// @param z0 Point to compute winding number around
		/// @param tol Integration tolerance (default: Precision::ComplexAnalysisTolerance)
		/// @return Winding number (integer for closed contours, rounded from numerical result)
		static int WindingNumber(const CircleContour& contour, Complex z0, Real tol = Precision::ComplexAnalysisTolerance)
		{
			// Integrand: 1/(z - z0)
			ComplexFunctionFromStdFunc g([z0](Complex z) -> Complex {
				return REAL(1.0) / (z - z0);
			});

			auto result = ContourIntegral(g, contour, tol);

			// n = (1/2πi) · ∮ dz/(z-z0)
			Complex winding = result.value / (REAL(2.0) * Constants::PI * Complex(0.0, 1.0));
			return static_cast<int>(std::round(winding.real()));
		}

		/********************************************************************************************************************/
		/********                           Cauchy Integral Formula                                                   ********/
		/********************************************************************************************************************/

		/// @brief Evaluate f(z0) using the Cauchy integral formula
		/// @details f(z0) = (1/2πi) ∮_γ f(z)/(z - z0) dz
		///          where γ is a simple closed contour enclosing z0.
		///          This is primarily useful for:
		///          - Validating function implementations against integral representations
		///          - Computing values at points where direct evaluation is numerically difficult
		///          - Educational demonstration of complex analysis theorems
		/// @param f Complex function f:C→C (must be analytic inside and on γ)
		/// @param contour Circle contour enclosing z0
		/// @param z0 Point at which to evaluate f
		/// @param tol Integration tolerance (default: Precision::ComplexAnalysisTolerance)
		/// @return f(z0) computed via contour integration
		static Complex CauchyIntegralFormula(
			const IComplexFunction& f,
			const CircleContour& contour,
			Complex z0,
			Real tol = Precision::ComplexAnalysisTolerance)
		{
			// Integrand: f(z) / (z - z0)
			ComplexFunctionFromStdFunc integrand_func([&f, z0](Complex z) -> Complex {
				return f(z) / (z - z0);
			});

			auto result = ContourIntegral(integrand_func, contour, tol);
			return result.value / (REAL(2.0) * Constants::PI * Complex(0.0, 1.0));
		}

		/// @brief Compute the n-th derivative f^(n)(z0) using the generalized Cauchy integral formula
		/// @details f^(n)(z0) = (n! / 2πi) ∮_γ f(z) / (z - z0)^(n+1) dz
		/// @param f Complex function (must be analytic inside and on the contour)
		/// @param contour Circle contour enclosing z0
		/// @param z0 Point at which to evaluate the derivative
		/// @param n Order of derivative (n >= 0)
		/// @param tol Integration tolerance (default: Precision::ComplexAnalysisTolerance)
		/// @return f^(n)(z0)
		static Complex CauchyDerivative(
			const IComplexFunction& f,
			const CircleContour& contour,
			Complex z0,
			int n,
			Real tol = Precision::ComplexAnalysisTolerance)
		{
			// Integrand: f(z) / (z - z0)^(n+1)
			ComplexFunctionFromStdFunc integrand_func([&f, z0, n](Complex z) -> Complex {
				Complex denom = z - z0;
				Complex denom_power = denom;
				for (int k = 0; k < n; ++k)
					denom_power *= denom;
				return f(z) / denom_power;
			});

			auto result = ContourIntegral(integrand_func, contour, tol);

			// Multiply by n! / (2πi)
			Real factorial = 1.0;
			for (int k = 2; k <= n; ++k)
				factorial *= k;

			return factorial * result.value / (REAL(2.0) * Constants::PI * Complex(0.0, 1.0));
		}

		/********************************************************************************************************************/
		/********                              Residue Computation                                                   ********/
		/********************************************************************************************************************/

		/// @brief Compute the residue of f at z0 via contour integration
		/// @details Res(f, z0) = (1/2πi) ∮_γ f(z) dz
		///          where γ is a small circle around z0 not enclosing other singularities.
		/// @param f Complex function (may have an isolated singularity at z0)
		/// @param z0 Isolated singularity
		/// @param radius Radius of integration circle (must be small enough to exclude
		///               other singularities)
		/// @param tol Integration tolerance (default: Precision::ComplexAnalysisTolerance)
		/// @return Residue of f at z0
		static Complex Residue(
			const IComplexFunction& f,
			Complex z0,
			Real radius = 0.1,
			Real tol = PrecisionValues<Real>::DefaultToleranceStrict)
		{
			CircleContour contour(z0, radius);
			auto result = ContourIntegral(f, contour, tol);
			return result.value / (REAL(2.0) * Constants::PI * Complex(0.0, 1.0));
		}

		/// @brief Compute residue of a simple pole: Res(f, z0) = lim_{z→z0} (z-z0)·f(z)
		/// @details For simple poles (order 1), this direct limit is faster and more
		///          accurate than contour integration.
		/// @param f Complex function with a simple pole at z0
		/// @param z0 Location of the simple pole
		/// @param h Small displacement for limit approximation (default: 1e-8)
		/// @return Residue at z0
		static Complex ResidueSimplePole(
			const IComplexFunction& f,
			Complex z0,
			Real h = PrecisionValues<Real>::DefaultToleranceStrict)
		{
			// Approximate lim_{z→z0} (z - z0) · f(z) from multiple directions
			// Average over 4 directions for robustness
			Complex sum(0.0, 0.0);
			const int n_dirs = 4;
			for (int k = 0; k < n_dirs; ++k)
			{
				Real angle = k * Constants::PI / (2.0 * n_dirs) + Constants::PI / 8.0;  // offset to avoid axes
				Complex dz = h * Complex(std::cos(angle), std::sin(angle));
				Complex z = z0 + dz;
				sum += dz * f(z);
			}
			return sum / Real(n_dirs);
		}

		/********************************************************************************************************************/
		/********                            Argument Principle                                                       ********/
		/********************************************************************************************************************/

		/// @brief Apply the argument principle: count zeros minus poles inside a contour
		/// @details N - P = (1/2πi) ∮_γ f'(z)/f(z) dz
		///          where N = number of zeros, P = number of poles of f inside γ
		///          (counted with multiplicity).
		///
		///          For polynomials and entire functions (no poles), this counts zeros.
		///
		/// @param f Complex function (must be meromorphic inside and analytic on γ)
		/// @param contour Circle contour
		/// @param tol Integration tolerance (default: DefaultToleranceStrict)
		/// @return N - P (rounded to nearest integer)
		static int ArgumentPrinciple(
			const IComplexFunction& f,
			const CircleContour& contour,
			Real tol = PrecisionValues<Real>::DefaultToleranceStrict)
		{
			// Integrand: f'(z) / f(z)
			// We compute f'(z) numerically and divide by f(z)
			ComplexFunctionFromStdFunc logarithmic_derivative([&f](Complex z) -> Complex {
				Complex fz = f(z);
				Complex dfz = Derivation::NDer4Complex(f, z);
				return dfz / fz;
			});

			auto result = ContourIntegral(logarithmic_derivative, contour, tol);
			Complex count = result.value / (REAL(2.0) * Constants::PI * Complex(0.0, 1.0));
			return static_cast<int>(std::round(count.real()));
		}

		/// @brief Count the number of zeros of f inside a contour (for functions with no poles)
		/// @details Equivalent to ArgumentPrinciple when f has no poles inside γ.
		///          This is a convenience alias with clearer semantics.
		/// @param f Complex function (must be analytic inside and on γ)
		/// @param contour Circle contour
		/// @param tol Integration tolerance (default: DefaultToleranceStrict)
		/// @return Number of zeros inside γ (counted with multiplicity)
		static int CountZeros(
			const IComplexFunction& f,
			const CircleContour& contour,
			Real tol = PrecisionValues<Real>::DefaultToleranceStrict)
		{
			return ArgumentPrinciple(f, contour, tol);
		}

	} // namespace ComplexAnalysis
} // namespace MML

#endif // MML_COMPLEX_ANALYSIS_H
