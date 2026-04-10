///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedRealFunction.h                                          ///
///  Description: 1D interpolation classes (Linear, Polynomial, Spline, Rational)     ///
///               Table-driven function approximation for real-valued functions       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedRealFunction.h
/// @brief 1D interpolation classes for function approximation from tabulated data.
/// @section interp1d_classes Classes
/// - RealFunctionInterpolated - Abstract base class with binary search location
/// - LinearInterpRealFunc - Piecewise linear interpolation
/// - PolynomInterpRealFunc - Polynomial interpolation (Neville's algorithm)
/// - SplineInterpRealFunc - Cubic spline with optional derivatives
/// - RationalInterpRealFunc - Rational function interpolation (Bulirsch-Stoer)
/// - BarycentricRationalInterp - Floater-Hormann barycentric interpolation (no poles)
/// - MonotoneCubicInterpRealFunc - Fritsch-Carlson monotone cubic Hermite (no overshoots)
/// @see Interpolation2DFunction.h for 2D interpolation
/// @see InterpolationParametricCurve.h for parametric curve interpolation
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_REAL_FUNCTION_H
#define MML_INTERPOLATED_REAL_FUNCTION_H

#include <algorithm>
#include <cmath>

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

namespace MML {

	/// @brief Policy for handling queries outside the interpolation data range [MinX, MaxX].
	/// @ingroup Interpolation
	enum class ExtrapolationPolicy {
		Allow,  ///< Silently extrapolate beyond data range (default, backward-compatible)
		Throw,  ///< Throw RealFuncInterpRuntimeError if x is outside [MinX, MaxX]
		Clamp   ///< Clamp x to [MinX, MaxX] before interpolating
	};

	/// @brief Abstract base class for interpolating real functions from tabulated data.
	/// RealFunctionInterpolated provides the foundation for all 1D interpolation methods.
	/// It stores the data points (x, y) and provides efficient binary search to locate
	/// the interval containing a query point.
	/// Derived classes implement `calcInterpValue()` to perform the actual interpolation
	/// using different algorithms (linear, polynomial, spline, rational).
	/// @note Data points are copied internally, so the original vectors can be discarded.
	/// @see LinearInterpRealFunc, PolynomInterpRealFunc, SplineInterpRealFunc, RationalInterpRealFunc
	/// @ingroup Interpolation

	class RealFunctionInterpolated : public IRealFunction {
	private:
		int _numPoints, _usedPoints;
		Vector<Real> _x, _y; // we are storing copies of given values!
		ExtrapolationPolicy _extrapolationPolicy = ExtrapolationPolicy::Allow;

	public:
		/// @brief Construct an interpolated function from data points.
		/// @param x Vector of abscissas (x-values), must be sorted
		/// @param y Vector of ordinates (y-values), same size as x
		/// @param usedPointsInInterpolation Number of points to use in local interpolation (mm)
		/// @throws RealFuncInterpInitError if sizes are invalid

		RealFunctionInterpolated(const Vector<Real>& x, const Vector<Real>& y, int usedPointsInInterpolation)
			: _x(x)
			, _y(y)
			, _numPoints(x.size())
			, _usedPoints(usedPointsInInterpolation) {
			// throw if not enough points
			if (_numPoints < 2 || _usedPoints < 2 || _usedPoints > _numPoints)
				throw RealFuncInterpInitError("RealFunctionInterpolated size error");
			// Check for duplicate x-values (causes division by zero in interpolation)
			for (int i = 0; i < _numPoints; i++)
				for (int j = i + 1; j < _numPoints; j++)
					if (_x[i] == _x[j])
						throw RealFuncInterpInitError("RealFunctionInterpolated: duplicate x-value at indices " + std::to_string(i) + " and " + std::to_string(j));
		}

		virtual ~RealFunctionInterpolated() {}

		/// @brief Calculate the interpolated value at x using points starting at startInd.
		/// @param startInd Starting index in the data arrays
		/// @param x The point at which to interpolate
		/// @return The interpolated value
		/// This pure virtual method must be implemented by derived classes to perform
		/// the specific interpolation algorithm.

		Real virtual calcInterpValue(int startInd, Real x) const = 0;

		/// @name Data Range Accessors
		/// @{

		/// /** @brief Get the minimum x-value in the data. */

		inline Real MinX() const { return X(0); }

		/// /** @brief Get the maximum x-value in the data. */

		inline Real MaxX() const { return X(_numPoints - 1); }

		/// @}

		/// @name Data Point Accessors
		/// @{

		/// /** @brief Get the i-th x-value. */

		inline Real X(int i) const { return _x[i]; }

		/// /** @brief Get the i-th y-value. */

		inline Real Y(int i) const { return _y[i]; }

		/// /** @brief Get the total number of data points. */

		inline int getNumPoints() const { return _numPoints; }

		/// /** @brief Get the number of points used in each local interpolation (mm). */

		inline int getInterpOrder() const { return _usedPoints; }

		/// @}

		/// @name Extrapolation Policy
		/// @{

		/// @brief Check if x is within the interpolation data range [MinX, MaxX].
		bool isInRange(Real x) const { return x >= MinX() && x <= MaxX(); }

		/// @brief Set the extrapolation policy for out-of-range queries.
		void setExtrapolationPolicy(ExtrapolationPolicy policy) { _extrapolationPolicy = policy; }

		/// @brief Get the current extrapolation policy.
		ExtrapolationPolicy getExtrapolationPolicy() const { return _extrapolationPolicy; }

		/// @}

		/// @brief Evaluate the interpolated function at point x.
		/// @param x The point at which to evaluate
		/// @return The interpolated value f(x)
		/// This operator locates the appropriate interval using binary search,
		/// then calls the derived class's `calcInterpValue()` method.

		Real operator()(Real x) const {
			if (_extrapolationPolicy != ExtrapolationPolicy::Allow && !isInRange(x)) {
				if (_extrapolationPolicy == ExtrapolationPolicy::Throw) {
					throw RealFuncInterpRuntimeError(
						"Interpolation: query point x=" + std::to_string(x) +
						" outside data range [" + std::to_string(MinX()) + ", " + std::to_string(MaxX()) + "]");
				}
				// Clamp policy
				x = std::clamp(x, MinX(), MaxX());
			}
			int startInd = locate(x);
			return calcInterpValue(startInd, x);
		}

		/// @brief Locate the interval containing x using binary search.
		/// @param x The query point
		/// @return Starting index j such that x is centered in [j, j+mm-1]
		/// The returned index is clamped to ensure the mm-point interpolation
		/// stencil stays within bounds. Works for both ascending and descending data.

		int locate(const Real x) const {
			int indUpper, indMid, indLower;
			bool isAscending = (_x[_numPoints - 1] >= _x[0]);

			indLower = 0;
			indUpper = _numPoints - 1;

			while (indUpper - indLower > 1) {
				indMid = (indUpper + indLower) / 2;
				if ((x >= _x[indMid]) == isAscending)
					indLower = indMid;
				else
					indUpper = indMid;
			}
			return std::max(0, std::min(_numPoints - _usedPoints, indLower - ((_usedPoints - 2) / 2)));
		}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                         LINEAR INTERPOLATION                                  ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Piecewise linear interpolation of tabulated data.
	/// LinearInterpRealFunc provides simple, robust interpolation by connecting
	/// adjacent data points with straight lines. It is C⁰ continuous (continuous
	/// but with discontinuous first derivative at data points).
	/// **Advantages:**
	/// - Simple and fast O(log n) lookup + O(1) evaluation
	/// - No oscillation or overshoot
	/// - Works well with noisy data
	/// - Optional extrapolation control
	/// **Limitations:**
	/// - Discontinuous derivatives at data points
	/// - May miss curvature in smooth data
	/// @section linear_formula Formula
	/// For @f$ x_i \le x \le x_{i+1} @f$:
	/// @f[
	/// f(x) = y_i + \frac{y_{i+1} - y_i}{x_{i+1} - x_i}(x - x_i)
	/// @f]
	/// @see SplineInterpRealFunc for smoother interpolation
	/// @ingroup Interpolation

	class LinearInterpRealFunc : public RealFunctionInterpolated {
		bool _extrapolateOutsideOfRange;

	public:
		/// @brief Construct a linear interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param extrapolateOutsideOfRange If true, extrapolate beyond data range; if false, return 0

		LinearInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, bool extrapolateOutsideOfRange = false)
			: RealFunctionInterpolated(xv, yv, 2)
			, _extrapolateOutsideOfRange(extrapolateOutsideOfRange) {}

		/// @brief Calculate linearly interpolated value.
		/// @param j Starting index of the interval
		/// @param x Point at which to interpolate
		/// @return Interpolated value (or 0 if out of range and extrapolation disabled)

		Real calcInterpValue(int j, Real x) const override {
			if (_extrapolateOutsideOfRange == false && (x < MinX() || x > MaxX()))
				return 0.0;

			if (X(j) == X(j + 1))
				return Y(j);
			else
				return Y(j) + ((x - X(j)) / (X(j + 1) - X(j)) * (Y(j + 1) - Y(j)));
		}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                      POLYNOMIAL INTERPOLATION                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Polynomial interpolation using Neville's algorithm.
	/// PolynomInterpRealFunc constructs an interpolating polynomial of degree mm-1
	/// through mm consecutive data points. The polynomial passes exactly through
	/// all used points and provides smooth, infinitely differentiable interpolation.
	/// **Algorithm:** Neville's algorithm builds the interpolating polynomial
	/// iteratively, computing successively higher-degree approximations. It also
	/// provides an error estimate as a byproduct.
	/// **Advantages:**
	/// - Exact interpolation through all data points
	/// - Built-in error estimation
	/// - C^∞ smooth within each local region
	/// **Limitations:**
	/// - Runge phenomenon: High-degree polynomials can oscillate wildly
	/// - Not recommended for more than ~10 points
	/// - Error can be large between data points for high degrees
	/// @warning For many data points, use SplineInterpRealFunc or BarycentricRationalInterp instead.
	/// @see BarycentricRationalInterp for oscillation-free alternative
	/// @ingroup Interpolation

	class PolynomInterpRealFunc : public RealFunctionInterpolated {
	private:
		mutable Real _errorEst;

	public:
		/// @brief Construct a polynomial interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param m Number of points to use in interpolation (polynomial degree = m-1)

		PolynomInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, int m)
			: RealFunctionInterpolated(xv, yv, m)
			, _errorEst(0.) {}

		/// /** @brief Get the error estimate from the last interpolation. */

		Real getLastErrorEst() const { return _errorEst; }

		/// @brief Calculate polynomial interpolated value using Neville's algorithm.
		/// @param startInd Starting index in the data array
		/// @param x Point at which to interpolate
		/// @return Interpolated value; error estimate stored in _errorEst
		/// Uses mm-point polynomial interpolation on the subrange [startInd, startInd+mm-1].

		Real calcInterpValue(int startInd, Real x) const override {
			int i, m, ns = 0;
			Real y, den, dif, dift, ho, hp, w, dy;
			Vector<Real> c(getInterpOrder()), d(getInterpOrder());

			dif = std::abs(x - X(startInd));

			// First we find the index ns of the closest table entry
			for (i = 0; i < getInterpOrder(); i++) {
				if ((dift = std::abs(x - X(startInd + i))) < dif) {
					ns = i;
					dif = dift;
				}
				c[i] = Y(startInd + i);
				d[i] = Y(startInd + i);
			}
			y = Y(startInd + ns--); // This is the initial approximation to y

			for (m = 1; m < getInterpOrder(); m++) // For each column of the tableau
			{
				for (i = 0; i < getInterpOrder() - m; i++) {
					// we loop over the current c's and d's and update
					ho = X(startInd + i) - x;
					hp = X(startInd + i + m) - x;
					w = c[i + 1] - d[i];

					if ((den = ho - hp) == 0.0)
						throw RealFuncInterpRuntimeError("Polynomial interpolation: identical x-values encountered");

					den = w / den;
					d[i] = hp * den;
					c[i] = ho * den;
				}
				dy = (2 * (ns + 1) < (getInterpOrder() - m) ? c[ns + 1] : d[ns--]);

				y += dy;
			}
			_errorEst = dy;
			return y;
		}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                     CUBIC SPLINE INTERPOLATION                                ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Cubic spline interpolation with optional derivatives and integration.
	/// SplineInterpRealFunc constructs a cubic spline that passes through all data
	/// points with continuous first and second derivatives (C² continuity). This
	/// is the recommended general-purpose interpolation method for smooth data.
	/// **Types of Splines:**
	/// - **Natural spline:** Second derivatives are zero at endpoints (default)
	/// - **Clamped spline:** First derivatives specified at endpoints
	/// **Algorithm:** Solves a tridiagonal system for the second derivatives at
	/// each data point, then uses these to construct cubic polynomials in each
	/// interval.
	/// @section spline_formula Spline Formula
	/// For @f$ x_i \le x \le x_{i+1} @f$:
	/// @f[
	/// S(x) = \frac{(x_{i+1}-x)^3 y_i'' + (x-x_i)^3 y_{i+1}''}{6h}
	/// + \frac{(x_{i+1}-x) y_i + (x-x_i) y_{i+1}}{h}
	/// - \frac{h[(x_{i+1}-x) y_i'' + (x-x_i) y_{i+1}'']}{6}
	/// @f]
	/// where @f$ h = x_{i+1} - x_i @f$
	/// @section spline_features Features
	/// - `Derivative(x)`: First derivative at any point
	/// - `SecondDerivative(x)`: Second derivative at any point
	/// - `Integrate(a, b)`: Definite integral over [a, b]
	/// @note For data with poles or rapid variations, consider RationalInterpRealFunc.
	/// @see LinearInterpRealFunc for faster but rougher interpolation
	/// @ingroup Interpolation

	class SplineInterpRealFunc : public RealFunctionInterpolated {
	private:
		Vector<Real> _secDerY; ///< Second derivatives at each data point

	public:
		/// @brief Construct a cubic spline interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param yp1 First derivative at first point (1e99 = natural spline)
		/// @param ypn First derivative at last point (1e99 = natural spline)

		SplineInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, Real yp1 = 1.e99, Real ypn = 1.e99)
			: RealFunctionInterpolated(xv, yv, 2)
			, _secDerY(xv.size()) {
			initSecDerivs(&xv[0], &yv[0], yp1, ypn);
		}

		/// @brief Initialize second derivatives by solving the tridiagonal system.
		/// @param xv Pointer to x-values
		/// @param yv Pointer to y-values
		/// @param yp1 First derivative at first point (1e99 = natural)
		/// @param ypn First derivative at last point (1e99 = natural)
		/// Solves the tridiagonal system for second derivatives y'' at each point.
		/// If yp1 and/or ypn are ≥ 1e99, uses natural boundary conditions
		/// (zero second derivative at that boundary).

		void initSecDerivs(const Real* xv, const Real* yv, Real yp1, Real ypn) {
			int i, k;
			Real p, qn, sig, un;

			int numPoints = (int)_secDerY.size();
			Vector<Real> u(numPoints - 1);

			if (yp1 > 0.99e99)
				_secDerY[0] = u[0] = 0.0;
			else {
				_secDerY[0] = -0.5;
				u[0] = (3.0 / (xv[1] - xv[0])) * ((yv[1] - yv[0]) / (xv[1] - xv[0]) - yp1);
			}
			for (i = 1; i < numPoints - 1; i++) {
				sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
				p = sig * _secDerY[i - 1] + 2.0;

				_secDerY[i] = (sig - 1.0) / p;

				u[i] = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
				u[i] = (6.0 * u[i] / (xv[i + 1] - xv[i - 1]) - sig * u[i - 1]) / p;
			}
			if (ypn > 0.99e99)
				qn = un = 0.0;
			else {
				qn = 0.5;
				un = (3.0 / (xv[numPoints - 1] - xv[numPoints - 2])) *
					 (ypn - (yv[numPoints - 1] - yv[numPoints - 2]) / (xv[numPoints - 1] - xv[numPoints - 2]));
			}

			_secDerY[numPoints - 1] = (un - qn * u[numPoints - 2]) / (qn * _secDerY[numPoints - 2] + 1.0);

			for (k = numPoints - 2; k >= 0; k--)
				_secDerY[k] = _secDerY[k] * _secDerY[k + 1] + u[k];
		}

		/// @name Spline Evaluation
		/// @{

		/// @brief Calculate cubic spline interpolated value.
		/// @param startInd Starting index of the interval
		/// @param x Point at which to interpolate
		/// @return Interpolated value using cubic spline formula

		Real calcInterpValue(int startInd, Real x) const override {
			int indLow = startInd, indUpp = startInd + 1;
			Real y, h, b, a;
			h = X(indUpp) - X(indLow);

			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Spline interpolation: zero interval width (identical x-values)");

			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			y = a * Y(indLow) + b * Y(indUpp) + ((POW3(a) - a) * _secDerY[indLow] + (POW3(b) - b) * _secDerY[indUpp]) * (h * h) / 6.0;

			return y;
		}

		/// @}

		/// @name Derivative Evaluation
		/// @{

		/// @brief Evaluate the first derivative of the spline at point x.
		/// @param x Point at which to evaluate the derivative
		/// @return The first derivative dy/dx

		Real Derivative(Real x) const {
			int startInd = locate(x);
			int indLow = startInd, indUpp = startInd + 1;
			Real h, b, a;
			h = X(indUpp) - X(indLow);

			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Spline derivative: zero interval width (identical x-values)");

			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			// Derivative of cubic spline:
			// dy/dx = (y[hi] - y[lo])/h - (3a² - 1)/6 * h * y2[lo] + (3b² - 1)/6 * h * y2[hi]
			Real dydx = (Y(indUpp) - Y(indLow)) / h - (3.0 * a * a - 1.0) / 6.0 * h * _secDerY[indLow] +
						(3.0 * b * b - 1.0) / 6.0 * h * _secDerY[indUpp];

			return dydx;
		}

		/// @brief Evaluate the second derivative of the spline at point x.
		/// @param x Point at which to evaluate
		/// @return The second derivative d²y/dx²
		/// The second derivative is linearly interpolated between the stored
		/// values at the data points (which define the cubic polynomials).

		Real SecondDerivative(Real x) const {
			int startInd = locate(x);
			int indLow = startInd, indUpp = startInd + 1;
			Real h, b, a;
			h = X(indUpp) - X(indLow);

			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Spline second derivative: zero interval width (identical x-values)");

			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			// Second derivative is linear interpolation of stored second derivatives
			return a * _secDerY[indLow] + b * _secDerY[indUpp];
		}

		/// @}

		/// @name Integration
		/// @{

		/// @brief Compute the definite integral of the spline from a to b.
		/// @param a Lower integration bound
		/// @param b Upper integration bound
		/// @return The definite integral ∫[a,b] S(x) dx
		/// @throws If bounds are outside the data range
		/// Uses exact integration of the cubic polynomial in each segment,
		/// summing contributions from all segments intersected by [a, b].

		Real Integrate(Real a, Real b) const {
			if (a > b)
				return -Integrate(b, a);
			if (a < MinX() || b > MaxX())
				throw DomainError("Spline integration: bounds exceed data domain");

			Real integral = 0.0;
			int iLow = locate(a);
			int iHigh = locate(b);

			// Helper lambda to integrate one segment from x0 to x1 within [X(i), X(i+1)]
			auto integrateSegment = [this](int i, Real x0, Real x1) -> Real {
				Real h = X(i + 1) - X(i);
				if (h == 0.0)
					return 0.0;

				// Compute a and b at both endpoints
				Real a0 = (X(i + 1) - x0) / h;
				Real b0 = (x0 - X(i)) / h;
				Real a1 = (X(i + 1) - x1) / h;
				Real b1 = (x1 - X(i)) / h;

				// Integral of spline segment: ∫[a*y_lo + b*y_hi + ((a³-a)*y2_lo + (b³-b)*y2_hi)*h²/6] dx
				// Using substitution and exact integration
				Real y_lo = Y(i), y_hi = Y(i + 1);
				Real y2_lo = _secDerY[i], y2_hi = _secDerY[i + 1];

				// Antiderivative evaluated at x1 minus at x0
				// For term a*y_lo: integral is -h/2 * a² * y_lo
				// For term b*y_hi: integral is h/2 * b² * y_hi
				// For (a³-a)*h²/6*y2_lo: integral is -h/6 * (a⁴/4 - a²/2) * h * y2_lo = -h³/6 * (a⁴/4 - a²/2) * y2_lo
				// For (b³-b)*h²/6*y2_hi: integral is h/6 * (b⁴/4 - b²/2) * h * y2_hi = h³/6 * (b⁴/4 - b²/2) * y2_hi

				Real term1 = h * y_lo * 0.5 * (a0 * a0 - a1 * a1);
				Real term2 = h * y_hi * 0.5 * (b1 * b1 - b0 * b0);
				Real term3 =
					h * h * h / 6.0 * y2_lo * ((a0 * a0 * a0 * a0 / 4.0 - a0 * a0 / 2.0) - (a1 * a1 * a1 * a1 / 4.0 - a1 * a1 / 2.0));
				Real term4 =
					h * h * h / 6.0 * y2_hi * ((b1 * b1 * b1 * b1 / 4.0 - b1 * b1 / 2.0) - (b0 * b0 * b0 * b0 / 4.0 - b0 * b0 / 2.0));

				return term1 + term2 + term3 + term4;
			};

			if (iLow == iHigh) {
				// Both endpoints in same segment
				integral = integrateSegment(iLow, a, b);
			} else {
				// Multiple segments
				integral += integrateSegment(iLow, a, X(iLow + 1));
				for (int i = iLow + 1; i < iHigh; i++) {
					integral += integrateSegment(i, X(i), X(i + 1));
				}
				integral += integrateSegment(iHigh, X(iHigh), b);
			}

			return integral;
		}

		/// @}

		/// @name Data Access
		/// @{

		/// @brief Get the stored second derivative at node i.
		/// @param i Index of the data point
		/// @return The second derivative y''[i]
		/// @throws IndexError if i is out of bounds

		Real GetSecondDerivative(int i) const {
			if (i < 0 || i >= getNumPoints())
				throw IndexError("Index out of range in GetSecondDerivative");
			return _secDerY[i];
		}

		/// @}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                      RATIONAL INTERPOLATION                                   ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Rational function interpolation using Bulirsch-Stoer algorithm.
	/// RationalInterpRealFunc constructs a diagonal rational function (equal degree
	/// numerator and denominator) through mm data points. This is superior to
	/// polynomial interpolation when the underlying function has poles or
	/// asymptotic behavior.
	/// **Algorithm:** Uses the Bulirsch-Stoer algorithm which builds up a continued
	/// fraction representation of the rational interpolant. Provides error estimation.
	/// **Advantages:**
	/// - Handles functions with poles or asymptotes
	/// - Better extrapolation behavior than polynomials
	/// - Built-in error estimation
	/// **Limitations:**
	/// - May detect poles during interpolation (throws exception)
	/// - More complex than polynomial methods
	/// @warning Will throw if a pole is detected during interpolation.
	/// @see PolynomInterpRealFunc for simpler polynomial alternative
	/// @see BarycentricRationalInterp for pole-free rational interpolation
	/// @ingroup Interpolation

	class RationalInterpRealFunc : public RealFunctionInterpolated {
	private:
		mutable Real _errorEst;
		static constexpr Real TINY = 1.0e-99; ///< Small number to prevent 0/0

	public:
		/// @brief Construct a rational interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param m Number of points to use (rational function has degree (m-1)/2 in each)

		RationalInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, int m)
			: RealFunctionInterpolated(xv, yv, m)
			, _errorEst(0.) {}

		/// /** @brief Get the error estimate from the last interpolation. */

		Real getLastErrorEst() const { return _errorEst; }

		/// @brief Calculate rational interpolated value using Bulirsch-Stoer algorithm.
		/// @param startInd Starting index in the data array
		/// @param x Point at which to interpolate
		/// @return Interpolated value; error estimate stored in _errorEst
		/// @throws RealFuncInterpRuntimeError if a pole is detected

		Real calcInterpValue(int startInd, Real x) const override {
			int m, i, ns = 0;
			Real y, w, t, hh, h, dd;
			int mm = getInterpOrder();
			Vector<Real> c(mm), d(mm);

			hh = std::abs(x - X(startInd));
			for (i = 0; i < mm; i++) {
				h = std::abs(x - X(startInd + i));
				if (h == 0.0) {
					_errorEst = 0.0;
					return Y(startInd + i); // Exact hit on a data point
				} else if (h < hh) {
					ns = i;
					hh = h;
				}
				c[i] = Y(startInd + i);
				d[i] = Y(startInd + i) + TINY; // TINY prevents rare zero-over-zero condition
			}
			y = Y(startInd + ns--);

			for (m = 1; m < mm; m++) {
				for (i = 0; i < mm - m; i++) {
					w = c[i + 1] - d[i];
					h = X(startInd + i + m) - x;
					if (std::abs(h) < TINY)
						throw RealFuncInterpRuntimeError("Error in RationalInterpRealFunc: x coincides with a node");
					t = (X(startInd + i) - x) * d[i] / h;
					dd = t - c[i + 1];
					if (dd == 0.0)
						throw RealFuncInterpRuntimeError("Error in RationalInterpRealFunc: pole detected");
					dd = w / dd;
					d[i] = c[i + 1] * dd;
					c[i] = t * dd;
				}
				_errorEst = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]);
				y += _errorEst;
			}
			return y;
		}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                BARYCENTRIC RATIONAL INTERPOLATION                             ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Barycentric rational interpolation with no poles on the real line.
	/// BarycentricRationalInterp implements the Floater-Hormann algorithm for
	/// rational interpolation that is guaranteed to have no poles on the real
	/// axis. This avoids the Runge phenomenon that plagues high-degree polynomial
	/// interpolation.
	/// **Algorithm:** Uses barycentric form with carefully computed weights that
	/// blend local polynomial interpolants. The blending parameter d controls
	/// the trade-off between smoothness and accuracy.
	/// **Parameter d:**
	/// - d = 0: Piecewise constant (step function)
	/// - d = 1: Piecewise linear
	/// - d = n-1: Polynomial interpolation (Lagrange form)
	/// - Recommended: d = 3 to 5 for most applications
	/// **Advantages:**
	/// - No poles on the real line (stable for any input)
	/// - Avoids Runge phenomenon
	/// - O(n) evaluation after O(n²) setup
	/// - Exact at data points
	/// @note This class implements IRealFunction directly, not RealFunctionInterpolated.
	/// @see PolynomInterpRealFunc for pure polynomial alternative
	/// @see RationalInterpRealFunc for Bulirsch-Stoer rational interpolation
	/// @ingroup Interpolation

	class BarycentricRationalInterp : public IRealFunction {
	private:
		int _n;			 ///< Number of data points
		int _d;			 ///< Blending parameter (0 to n-1)
		Vector<Real> _x; ///< Abscissas
		Vector<Real> _y; ///< Ordinates
		Vector<Real> _w; ///< Barycentric weights
		ExtrapolationPolicy _extrapolationPolicy; ///< Policy for out-of-range queries

		/// /** @brief Compute the barycentric weights using Floater-Hormann formula. */

		void computeWeights() {
			for (int k = 0; k < _n; k++) {
				int imin = std::max(k - _d, 0);
				int imax = k >= _n - _d ? _n - _d - 1 : k;
				Real temp = (imin & 1) ? -1.0 : 1.0; // Sign alternates
				Real sum = 0.0;

				for (int i = imin; i <= imax; i++) {
					int jmax = std::min(i + _d, _n - 1);
					Real term = 1.0;
					for (int j = i; j <= jmax; j++) {
						if (j == k)
							continue;
						term *= (_x[k] - _x[j]);
					}
					term = temp / term;
					temp = -temp;
					sum += term;
				}
				_w[k] = sum;
			}
		}

	public:
		/// @brief Construct a barycentric rational interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param d Blending parameter (0 to n-1), default 3
		/// @throws RealFuncInterpInitError if d >= n or d < 0

		BarycentricRationalInterp(const Vector<Real>& xv, const Vector<Real>& yv, int d = 3)
			: _n(xv.size())
			, _d(d)
			, _x(xv)
			, _y(yv)
			, _w(_n)
			, _extrapolationPolicy(ExtrapolationPolicy::Allow) {
			if (_n <= _d)
				throw RealFuncInterpInitError("BarycentricRationalInterp: d too large for number of points");
			if (_d < 0)
				throw RealFuncInterpInitError("BarycentricRationalInterp: d must be non-negative");
			computeWeights();
		}

		/// @name Data Range Accessors
		/// @{

		/// /** @brief Get the minimum x-value. */

		Real MinX() const { return _x[0]; }

		/// /** @brief Get the maximum x-value. */

		Real MaxX() const { return _x[_n - 1]; }

		/// /** @brief Get the number of data points. */

		int getNumPoints() const { return _n; }

		/// /** @brief Get the blending parameter d. */

		int getBlendingParameter() const { return _d; }

		/// @brief Check if x is within the interpolation data range [MinX, MaxX].
		bool isInRange(Real x) const { return x >= MinX() && x <= MaxX(); }

		/// @brief Set the extrapolation policy for out-of-range queries.
		void setExtrapolationPolicy(ExtrapolationPolicy policy) { _extrapolationPolicy = policy; }

		/// @brief Get the current extrapolation policy.
		ExtrapolationPolicy getExtrapolationPolicy() const { return _extrapolationPolicy; }

		/// @}

		/// @brief Evaluate the interpolated function at x.
		/// @param x Point at which to evaluate
		/// @return Interpolated value using barycentric formula

		Real operator()(Real x) const override {
			if (_extrapolationPolicy != ExtrapolationPolicy::Allow && !isInRange(x)) {
				if (_extrapolationPolicy == ExtrapolationPolicy::Throw) {
					throw RealFuncInterpRuntimeError(
						"BarycentricRationalInterp: query point x=" + std::to_string(x) +
						" outside data range [" + std::to_string(MinX()) + ", " + std::to_string(MaxX()) + "]");
				}
				x = std::clamp(x, MinX(), MaxX());
			}
			Real num = 0.0, den = 0.0;
			for (int i = 0; i < _n; i++) {
				Real h = x - _x[i];
				if (h == 0.0) {
					return _y[i]; // Exact hit on a data point
				}
				Real temp = _w[i] / h;
				num += temp * _y[i];
				den += temp;
			}
			return num / den;
		}

		/// @brief Get the barycentric weight at index i.
		/// @param i Index of the data point
		/// @return The barycentric weight w[i]
		/// @throws IndexError if i is out of bounds

		Real getWeight(int i) const {
			if (i < 0 || i >= _n)
				throw IndexError("Index out of range in getWeight");
			return _w[i];
		}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                MONOTONE CUBIC INTERPOLATION                                   ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Monotone cubic interpolation using the Fritsch-Carlson algorithm.
	/// MonotoneCubicInterpRealFunc constructs a piecewise cubic Hermite interpolant
	/// that preserves the monotonicity of the input data. Unlike standard cubic
	/// splines, this method guarantees no spurious oscillations or overshoots
	/// between data points when the data is monotone.
	/// **Algorithm:** Fritsch-Carlson (1980) computes initial tangent estimates
	/// from centered differences, then adjusts them to satisfy monotonicity
	/// constraints using the α-β condition (α² + β² ≤ 9).
	/// **Advantages:**
	/// - Preserves monotonicity of input data
	/// - No oscillations or overshoots between data points
	/// - C¹ continuous (continuous first derivative)
	/// - Exact at data points
	/// **Limitations:**
	/// - C¹ only (second derivative may be discontinuous at knots)
	/// - Slightly less accurate than unconstrained cubic spline for smooth data
	/// @section monotone_formula Hermite Basis
	/// For @f$ x_i \le x \le x_{i+1} @f$ with @f$ t = (x - x_i)/h_i @f$:
	/// @f[
	/// f(x) = (2t^3 - 3t^2 + 1)y_i + (t^3 - 2t^2 + t)h_i d_i
	/// + (-2t^3 + 3t^2)y_{i+1} + (t^3 - t^2)h_i d_{i+1}
	/// @f]
	/// @see SplineInterpRealFunc for C² cubic spline (may overshoot)
	/// @see LinearInterpRealFunc for simpler shape-preserving interpolation
	/// @ingroup Interpolation

	class MonotoneCubicInterpRealFunc : public RealFunctionInterpolated {
	private:
		Vector<Real> _d; ///< Tangent derivatives at each data point

	public:
		/// @brief Construct a monotone cubic interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @throws RealFuncInterpInitError if fewer than 2 points provided

		MonotoneCubicInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv)
			: RealFunctionInterpolated(xv, yv, 2)
			, _d(xv.size())
		{
			initDerivatives();
		}

		/// @brief Initialize tangent derivatives using the Fritsch-Carlson algorithm.
		/// Computes interval slopes δ_k, initial tangent estimates from centered
		/// differences, then enforces the monotonicity constraint α² + β² ≤ 9.

		void initDerivatives()
		{
			int n = getNumPoints();

			if (n == 2)
			{
				// Only one interval — use linear slope for both endpoints
				Real delta = (Y(1) - Y(0)) / (X(1) - X(0));
				_d[0] = delta;
				_d[1] = delta;
				return;
			}

			// Step 1: Compute interval slopes δ_k
			Vector<Real> delta(n - 1);
			for (int k = 0; k < n - 1; k++)
				delta[k] = (Y(k + 1) - Y(k)) / (X(k + 1) - X(k));

			// Step 2: Initial tangent estimates from three-point formula
			_d[0] = delta[0];
			for (int k = 1; k < n - 1; k++)
			{
				if (delta[k - 1] * delta[k] <= 0.0)
					_d[k] = 0.0;  // Different signs or zero — flat tangent
				else
					_d[k] = (delta[k - 1] + delta[k]) / 2.0;
			}
			_d[n - 1] = delta[n - 2];

			// Step 3: Fritsch-Carlson monotonicity correction
			for (int k = 0; k < n - 1; k++)
			{
				if (delta[k] == 0.0)
				{
					// Flat segment — both endpoint tangents must be zero
					_d[k] = 0.0;
					_d[k + 1] = 0.0;
				}
				else
				{
					Real alpha = _d[k] / delta[k];
					Real beta = _d[k + 1] / delta[k];
					Real r2 = alpha * alpha + beta * beta;
					if (r2 > 9.0)
					{
						// Rescale to satisfy α² + β² ≤ 9
						Real tau = 3.0 / std::sqrt(r2);
						_d[k] = tau * alpha * delta[k];
						_d[k + 1] = tau * beta * delta[k];
					}
				}
			}
		}

		/// @brief Calculate monotone cubic Hermite interpolated value.
		/// @param startInd Starting index of the interval
		/// @param x Point at which to interpolate
		/// @return Interpolated value using cubic Hermite basis functions

		Real calcInterpValue(int startInd, Real x) const override
		{
			int lo = startInd, hi = startInd + 1;
			Real h = X(hi) - X(lo);

			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Monotone cubic interpolation: zero interval width (identical x-values)");

			Real t = (x - X(lo)) / h;
			Real t2 = t * t;
			Real t3 = t2 * t;

			// Hermite basis functions
			Real h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
			Real h10 = t3 - 2.0 * t2 + t;
			Real h01 = -2.0 * t3 + 3.0 * t2;
			Real h11 = t3 - t2;

			return h00 * Y(lo) + h10 * h * _d[lo] + h01 * Y(hi) + h11 * h * _d[hi];
		}

		/// @brief Evaluate the first derivative at point x.
		/// @param x Point at which to evaluate the derivative
		/// @return The first derivative dy/dx

		Real Derivative(Real x) const
		{
			int startInd = locate(x);
			int lo = startInd, hi = startInd + 1;
			Real h = X(hi) - X(lo);

			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Monotone cubic derivative: zero interval width (identical x-values)");

			Real t = (x - X(lo)) / h;
			Real t2 = t * t;

			// Derivatives of Hermite basis functions (w.r.t. t), divided by h for dx
			Real dh00 = (6.0 * t2 - 6.0 * t) / h;
			Real dh10 = (3.0 * t2 - 4.0 * t + 1.0) / h;
			Real dh01 = (-6.0 * t2 + 6.0 * t) / h;
			Real dh11 = (3.0 * t2 - 2.0 * t) / h;

			return dh00 * Y(lo) + dh10 * h * _d[lo] + dh01 * Y(hi) + dh11 * h * _d[hi];
		}

		/// @brief Get the stored tangent derivative at node i.
		/// @param i Index of the data point
		/// @return The tangent derivative d[i]
		/// @throws IndexError if i is out of bounds

		Real GetDerivative(int i) const
		{
			if (i < 0 || i >= getNumPoints())
				throw IndexError("Index out of range in GetDerivative");
			return _d[i];
		}
	};

} // namespace MML

#endif // MML_INTERPOLATED_REAL_FUNCTION_H
