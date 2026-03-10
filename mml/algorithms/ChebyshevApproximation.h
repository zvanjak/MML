///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ChebyshevApproximation.h                                            ///
///  Description: Chebyshev polynomial approximation for function representation      ///
///               Clenshaw recurrence, integration, and derivative evaluation         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CHEBYSHEV_APPROXIMATION_H
#define MML_CHEBYSHEV_APPROXIMATION_H

#include "../MMLBase.h"
#include "../interfaces/IFunction.h"
#include "../base/Vector/Vector.h"
#include "../base/Polynom.h"
#include "../base/ChebyshevPolynom.h"

#include <functional>

namespace MML {

	///////////////////////////////////////////////////////////////////////////////////////////
	///                           CHEBYSHEV APPROXIMATION CLASS                             ///
	///////////////////////////////////////////////////////////////////////////////////////////

	// ChebyshevApproximation - Represents a function as a Chebyshev series expansion.
	//
	// A function f(x) on [a,b] is approximated as:
	//   f(x) ≈ Σ_{j=0}^{n-1} c_j * T_j(y)
	// where y = (2x - a - b) / (b - a) maps [a,b] to [-1,1]
	//
	// Key features:
	// - Near-minimax polynomial approximation
	// - Numerically stable Clenshaw evaluation
	// - Spectral convergence for smooth functions
	// - Implements IRealFunction for seamless integration
	//
	// Reference: Numerical Recipes 3rd Edition, Chapter 5
	class ChebyshevApproximation : public IRealFunction {
	private:
		Vector<Real> _coef; // Chebyshev coefficients c_0, c_1, ..., c_{n-1}
		Real _a, _b;		// Domain [a, b]
		int _m;				// Truncation degree (number of terms used in evaluation, ≤ n)

	public:
		// Default constructor - creates empty approximation
		ChebyshevApproximation()
			: _coef()
			, _a(-1.0)
			, _b(1.0)
			, _m(0) {}

		// Construct from a function on interval [a, b] using n Chebyshev coefficients.
		// Uses DCT-like transform at Chebyshev nodes for coefficient computation.
		//
		// Parameters:
		//   func - Function to approximate (implements IRealFunction)
		//   a, b - Domain interval [a, b]
		//   n    - Number of Chebyshev coefficients (default 50)
		//
		// Algorithm:
		// 1. Sample function at Chebyshev nodes: x_k = cos(π(k+0.5)/n)
		// 2. Compute coefficients via discrete cosine transform
		ChebyshevApproximation(const IRealFunction& func, Real a, Real b, int n = 50)
			: _coef(n)
			, _a(a)
			, _b(b)
			, _m(n) {
			if (n <= 0)
				throw std::invalid_argument("ChebyshevApproximation: n must be positive");
			if (a >= b)
				throw std::invalid_argument("ChebyshevApproximation: require a < b");

			const Real pi = Constants::PI;
			Real bma = 0.5 * (b - a); // half-width
			Real bpa = 0.5 * (b + a); // midpoint

			// Sample function at Chebyshev nodes
			Vector<Real> f(n);
			for (int k = 0; k < n; k++) {
				Real y = std::cos(pi * (k + 0.5) / n); // Chebyshev node on [-1,1]
				f[k] = func(y * bma + bpa);			   // Map to [a,b] and evaluate
			}

			// Compute Chebyshev coefficients via DCT-II formula
			// c_k = (2/n) * Σ_{j=0}^{n-1} f[j] * cos(π * k * (j + 0.5) / n)
			Real fac = 2.0 / n;
			for (int k = 0; k < n; k++) {
				Real sum = 0.0;
				for (int j = 0; j < n; j++)
					sum += f[j] * std::cos(pi * k * (j + 0.5) / n);
				_coef[k] = fac * sum;
			}
		}

		// Construct from std::function for convenience (lambdas, etc.)
		ChebyshevApproximation(std::function<Real(Real)> func, Real a, Real b, int n = 50)
			: _coef(n)
			, _a(a)
			, _b(b)
			, _m(n) {
			if (n <= 0)
				throw std::invalid_argument("ChebyshevApproximation: n must be positive");
			if (a >= b)
				throw std::invalid_argument("ChebyshevApproximation: require a < b");

			const Real pi = Constants::PI;
			Real bma = 0.5 * (b - a);
			Real bpa = 0.5 * (b + a);

			Vector<Real> f(n);
			for (int k = 0; k < n; k++) {
				Real y = std::cos(pi * (k + 0.5) / n);
				f[k] = func(y * bma + bpa);
			}

			// Compute Chebyshev coefficients via DCT-II formula
			Real fac = 2.0 / n;
			for (int k = 0; k < n; k++) {
				Real sum = 0.0;
				for (int j = 0; j < n; j++)
					sum += f[j] * std::cos(pi * k * (j + 0.5) / n);
				_coef[k] = fac * sum;
			}
		}

		// Construct from existing Chebyshev coefficients.
		// Useful for creating approximations from known coefficient sequences,
		// or from results of derivative/integral operations.
		//
		// Parameters:
		//   coefficients - Vector of Chebyshev coefficients
		//   a, b         - Domain interval [a, b]
		ChebyshevApproximation(const Vector<Real>& coefficients, Real a, Real b)
			: _coef(coefficients)
			, _a(a)
			, _b(b)
			, _m(coefficients.size()) {
			if (coefficients.size() == 0)
				throw std::invalid_argument("ChebyshevApproximation: coefficients cannot be empty");
			if (a >= b)
				throw std::invalid_argument("ChebyshevApproximation: require a < b");
		}

		// Copy constructor
		ChebyshevApproximation(const ChebyshevApproximation& other) = default;

		// Move constructor
		ChebyshevApproximation(ChebyshevApproximation&& other) noexcept = default;

		// Copy assignment
		ChebyshevApproximation& operator=(const ChebyshevApproximation& other) = default;

		// Move assignment
		ChebyshevApproximation& operator=(ChebyshevApproximation&& other) noexcept = default;

		// Destructor
		virtual ~ChebyshevApproximation() = default;

		///////////////////////////          EVALUATION          ///////////////////////////

		// Evaluate the Chebyshev approximation at point x using Clenshaw recurrence.
		// This is the IRealFunction interface implementation.
		//
		// The Clenshaw algorithm provides numerically stable O(n) evaluation
		// without explicitly computing Chebyshev polynomials.
		//
		// Algorithm:
		// 1. Map x from [a,b] to y in [-1,1]
		// 2. Apply backward Clenshaw recurrence
		// 3. Return final sum
		//
		// Throws if x is outside [a, b] (can be relaxed for extrapolation)
		Real operator()(Real x) const override {
			// Check domain (comment out for extrapolation)
			if ((x - _a) * (x - _b) > 0.0)
				throw std::domain_error("ChebyshevApproximation: x not in range [a, b]");

			if (_m == 0)
				return 0.0;

			// Map x from [a,b] to y in [-1,1]
			Real y = (2.0 * x - _a - _b) / (_b - _a);
			Real y2 = 2.0 * y;

			// Clenshaw recurrence (backward)
			Real d = 0.0, dd = 0.0;
			for (int j = _m - 1; j > 0; j--) {
				Real sv = d;
				d = y2 * d - dd + _coef[j];
				dd = sv;
			}

			return y * d - dd + 0.5 * _coef[0];
		}

		// Evaluate using only the first m terms (m ≤ _m)
		Real Eval(Real x, int m) const {
			if (m <= 0 || m > _m)
				throw std::invalid_argument("ChebyshevApproximation::Eval: invalid m");

			if ((x - _a) * (x - _b) > 0.0)
				throw std::domain_error("ChebyshevApproximation::Eval: x not in range [a, b]");

			Real y = (2.0 * x - _a - _b) / (_b - _a);
			Real y2 = 2.0 * y;

			Real d = 0.0, dd = 0.0;
			for (int j = m - 1; j > 0; j--) {
				Real sv = d;
				d = y2 * d - dd + _coef[j];
				dd = sv;
			}

			return y * d - dd + 0.5 * _coef[0];
		}

		///////////////////////////         ACCESSORS           ///////////////////////////

		// Get the full degree (number of coefficients - 1)
		int Degree() const { return static_cast<int>(_coef.size()) - 1; }

		// Get the truncation degree (number of terms used - 1)
		int TruncatedDegree() const { return _m - 1; }

		// Get number of coefficients stored
		int NumCoefficients() const { return static_cast<int>(_coef.size()); }

		// Get number of terms used in evaluation
		int NumTerms() const { return _m; }

		// Get the coefficient vector (read-only)
		const Vector<Real>& Coefficients() const { return _coef; }

		// Get individual coefficient
		Real Coefficient(int j) const {
			if (j < 0 || j >= static_cast<int>(_coef.size()))
				throw std::out_of_range("ChebyshevApproximation::Coefficient: index out of range");
			return _coef[j];
		}

		// Get domain endpoints
		Real DomainMin() const { return _a; }
		Real DomainMax() const { return _b; }

		///////////////////////////       TRUNCATION            ///////////////////////////

		// Set the number of terms to use in evaluation.
		// Allows trading accuracy for speed.
		void SetNumTerms(int m) {
			if (m <= 0 || m > static_cast<int>(_coef.size()))
				throw std::invalid_argument("ChebyshevApproximation::SetNumTerms: invalid m");
			_m = m;
		}

		// Automatically truncate based on coefficient magnitude threshold.
		// Returns the new number of terms.
		//
		// Chebyshev coefficients of smooth functions decay rapidly.
		// This method finds the smallest m such that |c_j| < threshold for all j ≥ m.
		int Truncate(Real threshold) {
			while (_m > 1 && std::abs(_coef[_m - 1]) < threshold)
				_m--;
			return _m;
		}

		///////////////////////////     CALCULUS OPERATIONS     ///////////////////////////

		// Compute the derivative of this Chebyshev approximation.
		// Returns a new ChebyshevApproximation representing f'(x).
		//
		// The derivative of a Chebyshev series has a simple recurrence formula:
		//   c'_{n-1} = 0
		//   c'_{n-2} = 2(n-1) * c_{n-1}
		//   c'_{j-1} = c'_{j+1} + 2j * c_j   for j = n-2, ..., 1
		// Then scale by 2/(b-a) for domain mapping.
		ChebyshevApproximation Derivative() const {
			int n = static_cast<int>(_coef.size());
			if (n <= 1) {
				// Derivative of constant is zero
				Vector<Real> zero_coef(1);
				zero_coef[0] = 0.0;
				return ChebyshevApproximation(zero_coef, _a, _b);
			}

			Vector<Real> cder(n);
			cder[n - 1] = 0.0;
			cder[n - 2] = 2.0 * (n - 1) * _coef[n - 1];

			for (int j = n - 2; j > 0; j--)
				cder[j - 1] = cder[j + 1] + 2.0 * j * _coef[j];

			// Scale for domain mapping: d/dx = (2/(b-a)) * d/dy
			Real con = 2.0 / (_b - _a);
			for (int j = 0; j < n; j++)
				cder[j] *= con;

			return ChebyshevApproximation(cder, _a, _b);
		}

		// Compute the indefinite integral of this Chebyshev approximation.
		// Returns a new ChebyshevApproximation representing ∫f(x)dx.
		// The constant of integration is chosen so that the integral is zero at x = a.
		//
		// The integral of a Chebyshev series:
		//   c''_j = (con) * (c_{j-1} - c_{j+1}) / j   for j = 1, ..., n-2
		//   c''_{n-1} = (con) * c_{n-2} / (n-1)
		//   c''_0 chosen to make integral(a) = 0
		// where con = (b-a)/4.
		ChebyshevApproximation Integral() const {
			int n = static_cast<int>(_coef.size());

			Vector<Real> cint(n);
			Real sum = 0.0;
			Real fac = 1.0;
			Real con = 0.25 * (_b - _a);

			for (int j = 1; j < n - 1; j++) {
				cint[j] = con * (_coef[j - 1] - _coef[j + 1]) / j;
				sum += fac * cint[j];
				fac = -fac;
			}

			cint[n - 1] = con * _coef[n - 2] / (n - 1);
			sum += fac * cint[n - 1];

			// c_0 chosen so integral(a) = 0
			cint[0] = 2.0 * sum;

			return ChebyshevApproximation(cint, _a, _b);
		}

		///////////////////////////   POLYNOMIAL CONVERSION     ///////////////////////////

		// Convert Chebyshev coefficients to power series (polynomial) coefficients.
		// Returns polynomial on [-1, 1]: p(y) = Σ d_j * y^j
		//
		// Note: For use on original domain [a,b], apply the substitution
		//       y = (2x - a - b) / (b - a)
		//
		// Uses backward recurrence algorithm from Numerical Recipes.
		Polynom<Real> ToPolynomial() const { return ToPolynomial(_m); }

		// Convert using only the first m terms
		Polynom<Real> ToPolynomial(int m) const {
			if (m <= 0 || m > static_cast<int>(_coef.size()))
				throw std::invalid_argument("ChebyshevApproximation::ToPolynomial: invalid m");

			std::vector<Real> d(m, 0.0);
			std::vector<Real> dd(m, 0.0);

			d[0] = _coef[m - 1];

			for (int j = m - 2; j > 0; j--) {
				for (int k = m - j; k > 0; k--) {
					Real sv = d[k];
					d[k] = 2.0 * d[k - 1] - dd[k];
					dd[k] = sv;
				}
				Real sv = d[0];
				d[0] = -dd[0] + _coef[j];
				dd[0] = sv;
			}

			for (int j = m - 1; j > 0; j--)
				d[j] = d[j - 1] - dd[j];

			d[0] = -dd[0] + 0.5 * _coef[0];

			return Polynom<Real>(d);
		}

		// Create ChebyshevApproximation from a polynomial.
		// Evaluates the polynomial at Chebyshev nodes and computes coefficients.
		static ChebyshevApproximation FromPolynomial(const Polynom<Real>& poly, Real a, Real b, int n = 0) {
			// If n not specified, use polynomial degree + 1
			if (n <= 0)
				n = poly.degree() + 1;

			// Create a wrapper function and use standard constructor
			return ChebyshevApproximation([&poly](Real x) { return poly(x); }, a, b, n);
		}

		///////////////////////////         UTILITIES           ///////////////////////////

		// Compute maximum absolute error compared to a reference function
		// by sampling at many points
		Real MaxError(const IRealFunction& func, int numSamples = 1000) const {
			Real maxErr = 0.0;
			for (int i = 0; i <= numSamples; i++) {
				Real x = _a + i * (_b - _a) / numSamples;
				Real err = std::abs((*this)(x)-func(x));
				if (err > maxErr)
					maxErr = err;
			}
			return maxErr;
		}

		Real MaxError(std::function<Real(Real)> func, int numSamples = 1000) const {
			Real maxErr = 0.0;
			for (int i = 0; i <= numSamples; i++) {
				Real x = _a + i * (_b - _a) / numSamples;
				Real err = std::abs((*this)(x)-func(x));
				if (err > maxErr)
					maxErr = err;
			}
			return maxErr;
		}

		// Print coefficients for debugging
		void PrintCoefficients(std::ostream& os = std::cout, int precision = 10) const {
			os << "ChebyshevApproximation on [" << _a << ", " << _b << "]" << std::endl;
			os << "Coefficients (" << _coef.size() << " stored, " << _m << " used):" << std::endl;
			os << std::scientific << std::setprecision(precision);
			for (int j = 0; j < static_cast<int>(_coef.size()); j++) {
				os << "  c[" << j << "] = " << _coef[j];
				if (j >= _m)
					os << " (truncated)";
				os << std::endl;
			}
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                           UTILITY FUNCTIONS                                         ///
	///////////////////////////////////////////////////////////////////////////////////////////

	// Compute the zeros (roots) of Chebyshev polynomial T_n
	// These are optimal interpolation nodes.
	// x_k = cos((2k + 1)π / (2n))  for k = 0, 1, ..., n-1
	inline Vector<Real> ChebyshevRoots(int n) {
		if (n <= 0)
			throw std::invalid_argument("ChebyshevRoots: n must be positive");

		Vector<Real> roots(n);
		const Real pi = Constants::PI;

		for (int k = 0; k < n; k++)
			roots[k] = std::cos(pi * (2 * k + 1) / (2.0 * n));

		return roots;
	}

	// Compute the extrema of Chebyshev polynomial T_n (Chebyshev-Lobatto points)
	// These include the endpoints and are used in Clenshaw-Curtis quadrature.
	// x_k = cos(kπ / n)  for k = 0, 1, ..., n
	inline Vector<Real> ChebyshevExtrema(int n) {
		if (n <= 0)
			throw std::invalid_argument("ChebyshevExtrema: n must be positive");

		Vector<Real> extrema(n + 1);
		const Real pi = Constants::PI;

		for (int k = 0; k <= n; k++)
			extrema[k] = std::cos(pi * k / n);

		return extrema;
	}

} // namespace MML

#endif // MML_CHEBYSHEV_APPROXIMATION_H
