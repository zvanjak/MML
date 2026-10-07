///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedFunctionSpline.h                                        ///
///  Description: Spline and monotone cubic interpolation                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedFunctionSpline.h
/// @brief Spline and monotone cubic interpolation.
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_FUNCTION_SPLINE_H
#define MML_INTERPOLATED_FUNCTION_SPLINE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h>

namespace MML {
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

		static bool IsNaturalBoundaryDerivative(Real derivative) {
			if constexpr (std::is_same_v<Real, float>)
				return !std::isfinite(derivative) || derivative >= std::numeric_limits<Real>::max() / REAL(2.0);
			else
				return !std::isfinite(derivative) || derivative > static_cast<Real>(0.99e99L);
		}

	public:
		static constexpr Real NaturalBoundaryDerivative = std::numeric_limits<Real>::max();

		/// @brief Construct a cubic spline interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param yp1 First derivative at first point (NaturalBoundaryDerivative = natural spline)
		/// @param ypn First derivative at last point (NaturalBoundaryDerivative = natural spline)

		SplineInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv,
			Real yp1 = NaturalBoundaryDerivative, Real ypn = NaturalBoundaryDerivative)
			: RealFunctionInterpolated(xv, yv, 2)
			, _secDerY(xv.size()) {
			initSecDerivs(&xv[0], &yv[0], yp1, ypn);
		}

		/// @brief Initialize second derivatives by solving the tridiagonal system.
		/// @param xv Pointer to x-values
		/// @param yv Pointer to y-values
		/// @param yp1 First derivative at first point (NaturalBoundaryDerivative = natural)
		/// @param ypn First derivative at last point (NaturalBoundaryDerivative = natural)
		/// Solves the tridiagonal system for second derivatives y'' at each point.
		/// If yp1 and/or ypn use the natural-boundary sentinel, uses natural boundary conditions
		/// (zero second derivative at that boundary).

		void initSecDerivs(const Real* xv, const Real* yv, Real yp1, Real ypn) {
			int i, k;
			Real p, qn, sig, un;

			int numPoints = (int)_secDerY.size();
			Vector<Real> u(numPoints - 1);

			if (IsNaturalBoundaryDerivative(yp1))
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
			if (IsNaturalBoundaryDerivative(ypn))
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

		const char* InterpolationMethodName() const override { return "CubicSpline"; }

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

		const char* InterpolationMethodName() const override { return "MonotoneCubicHermite"; }

		/// @brief Evaluate the first derivative at point x.
		/// @param x Point at which to evaluate the derivative
		/// @return The first derivative dy/dx

		Real Derivative(Real x) const
		{
			int startInd = locate(x);
			return detail::DifferentiateCubicHermiteSegment(x, X(startInd), X(startInd + 1),
				Y(startInd), Y(startInd + 1), _d[startInd], _d[startInd + 1]);
		}

		Real Integrate(Real from, Real to) const {
			if (from > to) return -Integrate(to, from);
			if (from < MinX() || to > MaxX())
				throw DomainError("Monotone cubic interpolation: integration bounds exceed data domain");
			Real result = 0.0;
			for (int i = 0; i < getNumPoints() - 1; ++i) {
				Real segmentMin = std::min(X(i), X(i + 1));
				Real segmentMax = std::max(X(i), X(i + 1));
				Real left = std::max(from, segmentMin);
				Real right = std::min(to, segmentMax);
				if (left < right)
					result += detail::IntegrateCubicHermiteSegment(left, right, X(i), X(i + 1),
						Y(i), Y(i + 1), _d[i], _d[i + 1]);
			}
			return result;
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

#endif // MML_INTERPOLATED_FUNCTION_SPLINE_H