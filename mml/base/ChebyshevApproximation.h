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
#include "../MMLNumericValidation.h"
#include "../interfaces/IFunction.h"
#include "../base/Vector/Vector.h"
#include "../base/Polynom.h"
#include "../base/ChebyshevPolynom.h"
#include "AlgorithmTypes.h"

#include <algorithm>
#include <functional>
#include <string>

namespace MML {

	///////////////////////////////////////////////////////////////////////////////////////////
	///                           CHEBYSHEV APPROXIMATION CLASS                             ///
	///////////////////////////////////////////////////////////////////////////////////////////

	// ChebyshevApproximation - Represents a function as a Chebyshev series expansion.
	//
	// A function f(x) on [a,b] is approximated as:
	//   f(x) ≈ c_0/2 + Σ_{j=1}^{n-1} c_j * T_j(y)
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

		static int ValidateConstruction(int n, Real a, Real b) {
			if (n <= 0)
				throw ArgumentError("ChebyshevApproximation: number of coefficients must be positive");
			ValidateBounds(a, b, "ChebyshevApproximation");
			if (!(a < b) || !std::isfinite(b - a))
				throw ArgumentError("ChebyshevApproximation: require a < b with a representable span");
			return n;
		}

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
			: ChebyshevApproximation(std::function<Real(Real)>(
				[&func](Real x) { return func(x); }), a, b, n) {}

		// Construct from std::function for convenience (lambdas, etc.)
		ChebyshevApproximation(std::function<Real(Real)> func, Real a, Real b, int n = 50)
			: _coef(ValidateConstruction(n, a, b))
			, _a(a)
			, _b(b)
			, _m(n) {
			if (!func)
				throw ArgumentError("ChebyshevApproximation: function must not be empty");

			const Real pi = Constants::PI;
			Real bma = MML::Real(0.5) * (b - a);
			Real bpa = a + bma;

			Vector<Real> f(n);
			for (int k = 0; k < n; k++) {
				Real y = std::cos(pi * (k + MML::Real(0.5)) / n);
				f[k] = func(y * bma + bpa);
				ValidateFunctionValue(f[k], "Chebyshev sample");
			}

			// Compute Chebyshev coefficients via DCT-II formula
			Real fac = MML::Real(2.0) / n;
			for (int k = 0; k < n; k++) {
				Real sum = 0.0;
				for (int j = 0; j < n; j++)
					sum += (fac * f[j]) * std::cos(pi * k * (j + MML::Real(0.5)) / n);
				ValidateFunctionValue(sum, "Chebyshev coefficient");
				_coef[k] = sum;
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
			, _m(ValidateConstruction(coefficients.size(), a, b)) {
			for (int i = 0; i < _m; ++i)
				ValidateFunctionValue(_coef[i], "Chebyshev coefficient");
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
			ValidateFinite(x, "Chebyshev query");
			if (x < _a || x > _b)
				throw DomainError("ChebyshevApproximation: x not in range [a, b]");

			if (_m == 0)
				return 0.0;

			// Map x from [a,b] to y in [-1,1]
			Real y = MML::Real(2.0) * ((x - _a) / (_b - _a)) - MML::Real(1.0);
			Real y2 = MML::Real(2.0) * y;

			// Clenshaw recurrence (backward)
			Real d = 0.0, dd = 0.0;
			for (int j = _m - 1; j > 0; j--) {
				Real sv = d;
				d = y2 * d - dd + _coef[j];
				dd = sv;
			}

			return y * d - dd + MML::Real(0.5) * _coef[0];
		}

		// Evaluate using only the first m terms (m ≤ _m)
		Real Eval(Real x, int m) const {
			if (m <= 0 || m > _m)
				throw ArgumentError("ChebyshevApproximation::Eval: invalid m");

			ValidateFinite(x, "Chebyshev query");
			if (x < _a || x > _b)
				throw DomainError("ChebyshevApproximation::Eval: x not in range [a, b]");

			Real y = MML::Real(2.0) * ((x - _a) / (_b - _a)) - MML::Real(1.0);
			Real y2 = MML::Real(2.0) * y;

			Real d = 0.0, dd = 0.0;
			for (int j = m - 1; j > 0; j--) {
				Real sv = d;
				d = y2 * d - dd + _coef[j];
				dd = sv;
			}

			return y * d - dd + MML::Real(0.5) * _coef[0];
		}

		Vector<Real> Evaluate(const Vector<Real>& points) const {
			Vector<Real> values(points.size());
			for (int i = 0; i < points.size(); ++i) values[i] = (*this)(points[i]);
			return values;
		}

		Vector<Real> UniformGrid(int numPoints) const {
			if (numPoints < 2)
				throw ArgumentError("ChebyshevApproximation::UniformGrid requires at least two points");
			Vector<Real> points(numPoints);
			for (int i = 0; i < numPoints; ++i)
				points[i] = i == numPoints - 1 ? _b : _a + (_b - _a) * (Real(i) / (numPoints - 1));
			return points;
		}

		void SampleUniform(int numPoints, Vector<Real>& points, Vector<Real>& values) const {
			points = UniformGrid(numPoints);
			values = Evaluate(points);
		}

		///////////////////////////         ACCESSORS           ///////////////////////////

		// Get the full degree (number of coefficients - 1)
		int degree() const { return static_cast<int>(_coef.size()) - 1; }

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
				throw IndexError("ChebyshevApproximation::Coefficient: index out of range");
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
				throw ArgumentError("ChebyshevApproximation::SetNumTerms: invalid m");
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
			int n = _m;
			if (n <= 1) {
				// Derivative of constant is zero
				Vector<Real> zero_coef(1);
				zero_coef[0] = 0.0;
				return ChebyshevApproximation(zero_coef, _a, _b);
			}

			Vector<Real> cder(n);
			cder[n - 1] = 0.0;
			cder[n - 2] = MML::Real(2.0) * (n - 1) * _coef[n - 1];

			for (int j = n - 2; j > 0; j--)
				cder[j - 1] = cder[j + 1] + MML::Real(2.0) * j * _coef[j];

			// Scale for domain mapping: d/dx = (2/(b-a)) * d/dy
			Real con = MML::Real(2.0) / (_b - _a);
			for (int j = 0; j < n; j++)
				cder[j] *= con;

			return ChebyshevApproximation(cder, _a, _b);
		}

		// Compute the indefinite integral of this Chebyshev approximation.
		// Returns a new ChebyshevApproximation representing ∫f(x)dx.
		// The constant of integration is chosen so that the integral is zero at x = a.
		//
		// The integral of a Chebyshev series:
		//   c''_j = con * (c_{j-1} - c_{j+1}) / j, for j = 1, ..., n
		// Missing coefficients are zero; n active terms require n+1 integral terms.
		//   c''_0 chosen to make integral(a) = 0
		// where con = (b-a)/4.
		ChebyshevApproximation Integral() const {
			int n = _m;
			if (n == 0)
				return ChebyshevApproximation(Vector<Real>{0}, _a, _b);
			if (n == std::numeric_limits<int>::max())
				throw ArgumentError("ChebyshevApproximation::Integral: coefficient count exceeds integer range");

			Vector<Real> cint(n + 1);
			Real sum = 0.0;
			Real fac = 1.0;
			Real con = MML::Real(0.25) * (_b - _a);

			for (int j = 1; j <= n; j++) {
				const Real next = j + 1 < n ? _coef[j + 1] : 0;
				cint[j] = con * (_coef[j - 1] - next) / j;
				sum += fac * cint[j];
				fac = -fac;
			}

			// c_0 chosen so integral(a) = 0
			cint[0] = MML::Real(2.0) * sum;

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
				throw ArgumentError("ChebyshevApproximation::ToPolynomial: invalid m");

			std::vector<Real> d(m, 0.0);
			std::vector<Real> dd(m, 0.0);

			d[0] = _coef[m - 1];

			for (int j = m - 2; j > 0; j--) {
				for (int k = m - j; k > 0; k--) {
					Real sv = d[k];
					d[k] = MML::Real(2.0) * d[k - 1] - dd[k];
					dd[k] = sv;
				}
				Real sv = d[0];
				d[0] = -dd[0] + _coef[j];
				dd[0] = sv;
			}

			for (int j = m - 1; j > 0; j--)
				d[j] = d[j - 1] - dd[j];

			d[0] = -dd[0] + MML::Real(0.5) * _coef[0];

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
			return MaxError(std::function<Real(Real)>(
				[&func](Real x) { return func(x); }), numSamples);
		}

		Real MaxError(std::function<Real(Real)> func, int numSamples = 1000) const {
			if (numSamples <= 0 || !func)
				throw ArgumentError("ChebyshevApproximation::MaxError: require positive sampling count and a function");
			Real maxErr = 0.0;
			for (int i = 0; ; ++i) {
				Real x = i == numSamples ? _b : _a + (_b - _a) * (Real(i) / numSamples);
				Real reference = func(x);
				ValidateFunctionValue(reference, "Chebyshev error reference");
				Real err = std::abs((*this)(x) - reference);
				ValidateFunctionValue(err, "Chebyshev error estimate");
				maxErr = std::max(maxErr, err);
				if (i == numSamples) break;
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

	struct ChebyshevApproximationConfig {
		Real tolerance = PrecisionValues<Real>::DefaultTolerance;
		int initial_terms = 8;
		int max_terms = 256;
		int validation_samples = 257;
		int tail_terms = 4;
	};

	struct ChebyshevApproximationResult : public EvaluationResultBase {
		ChebyshevApproximation approximation;
		Vector<Real> coefficients;
		int degree = -1;
		int terms_used = 0;
		int samples_used = 0;
		Real coefficient_tail_error = NoEstimate;
		Real max_error_estimate = NoEstimate;
		Real rms_error_estimate = NoEstimate;
		bool converged = false; // Only adaptive construction tests error estimates against a tolerance.
	};

	namespace ChebyshevApproximationDetail {
		inline Real CoefficientTailError(const ChebyshevApproximation& approximation, int tailTerms) {
			int first = std::max(1, approximation.NumCoefficients() - std::max(1, tailTerms));
			Real tail = 0.0;
			for (int i = first; i < approximation.NumCoefficients(); ++i)
				tail += std::abs(approximation.Coefficient(i));
			ValidateFunctionValue(tail, "Chebyshev coefficient tail");
			return tail;
		}

		template<typename Function>
		void EstimateErrors(const ChebyshevApproximation& approximation, Function&& function,
			int sampleCount, Real& maxError, Real& rmsError) {
			maxError = 0.0;
			rmsError = 0.0;
			int count = sampleCount;
			const Real rootCount = std::sqrt(Real(count));
			for (int i = 0; i < count; ++i) {
				Real x = i == count - 1 ? approximation.DomainMax() : approximation.DomainMin()
					+ (approximation.DomainMax() - approximation.DomainMin()) * (Real(i) / (count - 1));
				Real reference = function(x);
				ValidateFunctionValue(reference, "Chebyshev validation sample");
				Real error = approximation(x) - reference;
				ValidateFunctionValue(error, "Chebyshev validation error");
				maxError = std::max(maxError, std::abs(error));
				rmsError = std::hypot(rmsError, error / rootCount);
			}
			ValidateFunctionValue(rmsError, "Chebyshev RMS error");
		}

		template<typename Function>
		ChebyshevApproximationResult BuildResult(Function&& function, Real a, Real b, int terms,
			int validationSamples, int tailTerms) {
			if (validationSamples < 2 || tailTerms < 1
				|| terms > std::numeric_limits<int>::max() - validationSamples)
				throw ArgumentError("Chebyshev diagnostics: invalid sampling counts");
			AlgorithmTimer timer;
			ChebyshevApproximationResult result;
			result.algorithm_name = "ChebyshevApproximation";
			result.approximation = ChebyshevApproximation(
				std::function<Real(Real)>(function), a, b, terms);
			result.coefficients = result.approximation.Coefficients();
			result.degree = result.approximation.degree();
			result.terms_used = result.approximation.NumTerms();
			result.samples_used = terms + validationSamples;
			result.function_evaluations = result.samples_used;
			result.coefficient_tail_error = CoefficientTailError(result.approximation, tailTerms);
			EstimateErrors(result.approximation, std::forward<Function>(function), validationSamples,
				result.max_error_estimate, result.rms_error_estimate);
			result.status = AlgorithmStatus::Completed;
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	}

	inline ChebyshevApproximationResult ApproximateChebyshevDetailed(
		const std::function<Real(Real)>& function, Real a, Real b, int terms,
		int validationSamples = 257, int tailTerms = 4) {
		auto result = ChebyshevApproximationDetail::BuildResult(
			function, a, b, terms, validationSamples, tailTerms);
		return result;
	}

	inline ChebyshevApproximationResult ApproximateChebyshevAdaptive(
		const std::function<Real(Real)>& function, Real a, Real b,
		const ChebyshevApproximationConfig& config = {}) {
		ValidateTolerance(config.tolerance, "ApproximateChebyshevAdaptive");
		if (config.initial_terms < 2 || config.max_terms < config.initial_terms
			|| config.validation_samples < 2 || config.tail_terms < 1)
			throw ArgumentError("ApproximateChebyshevAdaptive: invalid configuration");

		int terms = config.initial_terms;
		int totalSamples = 0;
		ChebyshevApproximationResult result;
		while (true) {
			result = ChebyshevApproximationDetail::BuildResult(function, a, b, terms,
				config.validation_samples, config.tail_terms);
			if (result.samples_used > std::numeric_limits<int>::max() - totalSamples)
				throw NumericalMethodError("Chebyshev adaptive sample count exceeds integer range");
			totalSamples += result.samples_used;
			result.samples_used = totalSamples;
			result.function_evaluations = totalSamples;
			result.converged = result.coefficient_tail_error <= config.tolerance
				&& result.max_error_estimate <= config.tolerance;
			if (result.converged) {
				result.status = AlgorithmStatus::Success;
				return result;
			}
			if (terms >= config.max_terms) {
				result.status = AlgorithmStatus::ToleranceUnachievable;
				result.error_message = "Chebyshev approximation reached max_terms before tolerance";
				return result;
			}
			terms += std::min(terms, config.max_terms - terms);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	///                           UTILITY FUNCTIONS                                         ///
	///////////////////////////////////////////////////////////////////////////////////////////

	// Compute the zeros (roots) of Chebyshev polynomial T_n
	// These are optimal interpolation nodes.
	// x_k = cos((2k + 1)π / (2n))  for k = 0, 1, ..., n-1
	inline Vector<Real> ChebyshevRoots(int n) {
		if (n <= 0)
			throw ArgumentError("ChebyshevRoots: n must be positive");

		Vector<Real> roots(n);
		const Real pi = Constants::PI;

		for (int k = 0; k < n; k++)
			roots[k] = std::cos(pi * (2 * k + 1) / (MML::Real(2.0) * n));

		return roots;
	}

	// Compute the extrema of Chebyshev polynomial T_n (Chebyshev-Lobatto points)
	// These include the endpoints and are used in Clenshaw-Curtis quadrature.
	// x_k = cos(kπ / n)  for k = 0, 1, ..., n
	inline Vector<Real> ChebyshevExtrema(int n) {
		if (n <= 0)
			throw ArgumentError("ChebyshevExtrema: n must be positive");

		Vector<Real> extrema(n + 1);
		const Real pi = Constants::PI;

		for (int k = 0; k <= n; k++)
			extrema[k] = std::cos(pi * k / n);

		return extrema;
	}

} // namespace MML

#endif // MML_CHEBYSHEV_APPROXIMATION_H
