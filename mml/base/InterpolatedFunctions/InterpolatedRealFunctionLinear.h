///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedRealFunctionLinear.h                                    ///
///  Description: Shared 1D interpolation base and linear interpolation               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedRealFunctionLinear.h
/// @brief Shared 1D interpolation base and piecewise linear interpolation.
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_REAL_FUNCTION_LINEAR_H
#define MML_INTERPOLATED_REAL_FUNCTION_LINEAR_H

#include <algorithm>
#include <cmath>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/InterpolatedFunctions/InterpolationTypes.h>

namespace MML {
	namespace detail {
		inline Real EvaluateCubicHermiteSegment(Real x, Real x0, Real x1,
			Real y0, Real y1, Real derivative0, Real derivative1) {
			Real h = x1 - x0;
			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Cubic Hermite interpolation: zero interval width");
			Real t = (x - x0) / h;
			Real t2 = t * t;
			Real t3 = t2 * t;
			return (REAL(2.0) * t3 - REAL(3.0) * t2 + REAL(1.0)) * y0
				+ (t3 - REAL(2.0) * t2 + t) * h * derivative0
				+ (-REAL(2.0) * t3 + REAL(3.0) * t2) * y1
				+ (t3 - t2) * h * derivative1;
		}

		inline Real DifferentiateCubicHermiteSegment(Real x, Real x0, Real x1,
			Real y0, Real y1, Real derivative0, Real derivative1) {
			Real h = x1 - x0;
			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Cubic Hermite derivative: zero interval width");
			Real t = (x - x0) / h;
			Real t2 = t * t;
			return (REAL(6.0) * t2 - REAL(6.0) * t) * y0 / h
				+ (REAL(3.0) * t2 - REAL(4.0) * t + REAL(1.0)) * derivative0
				+ (-REAL(6.0) * t2 + REAL(6.0) * t) * y1 / h
				+ (REAL(3.0) * t2 - REAL(2.0) * t) * derivative1;
		}

		inline Real IntegrateCubicHermiteSegment(Real from, Real to, Real x0, Real x1,
			Real y0, Real y1, Real derivative0, Real derivative1) {
			Real h = x1 - x0;
			if (h == 0.0)
				throw RealFuncInterpRuntimeError("Cubic Hermite integration: zero interval width");
			auto primitive = [&](Real x) {
				Real t = (x - x0) / h;
				Real t2 = t * t;
				Real t3 = t2 * t;
				Real t4 = t3 * t;
				Real h00 = REAL(0.5) * t4 - t3 + t;
				Real h10 = REAL(0.25) * t4 - REAL(2.0 / 3.0) * t3 + REAL(0.5) * t2;
				Real h01 = -REAL(0.5) * t4 + t3;
				Real h11 = REAL(0.25) * t4 - REAL(1.0 / 3.0) * t3;
				return h * (y0 * h00 + h * derivative0 * h10 + y1 * h01 + h * derivative1 * h11);
			};
			return primitive(to) - primitive(from);
		}
	}

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
		virtual Real calcDetailedValue(int startInd, Real x) const { return calcInterpValue(startInd, x); }

		virtual const char* InterpolationMethodName() const { return "Interpolation"; }

		virtual Real DetailedErrorEstimate() const { return 0.0; }

		/// @name Data Range Accessors
		/// @{

		/// /** @brief Get the minimum x-value in the data. */

		inline Real MinX() const { return std::min(X(0), X(_numPoints - 1)); }

		/// /** @brief Get the maximum x-value in the data. */

		inline Real MaxX() const { return std::max(X(0), X(_numPoints - 1)); }

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

		InterpolationResult EvaluateDetailed(Real x, const InterpolationConfig& config = {}) const {
			AlgorithmTimer timer;
			InterpolationResult result;
			result.algorithm_name = InterpolationMethodName();
			result.query_point = x;
			result.evaluated_point = x;

			try {
				if (!isInRange(x)) {
					if (config.extrapolation_policy == ExtrapolationPolicy::Throw) {
						result.status = AlgorithmStatus::InvalidInput;
						result.interpolation_status = InterpolationStatus::OutOfRange;
						result.error_message = "Interpolation query is outside the data range";
						if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
							throw RealFuncInterpRuntimeError(result.error_message);
						result.elapsed_time_ms = timer.elapsed_ms();
						return result;
					}
					if (config.extrapolation_policy == ExtrapolationPolicy::Clamp) {
						result.evaluated_point = std::clamp(x, MinX(), MaxX());
						result.interpolation_status = InterpolationStatus::Clamped;
					} else {
						result.interpolation_status = InterpolationStatus::Extrapolated;
					}
				}

				result.interval_index = locate(result.evaluated_point);
				result.value = calcDetailedValue(result.interval_index, result.evaluated_point);
				result.error_estimate = config.estimate_error ? DetailedErrorEstimate() : 0.0;
				result.function_evaluations = 1;
				if (config.check_finite && (!std::isfinite(result.value) || !std::isfinite(result.error_estimate))) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = "Interpolation produced a non-finite result";
				}
			}
			catch (const RealFuncInterpRuntimeError& error) {
				if (config.exception_policy == EvaluationExceptionPolicy::Propagate) throw;
				result.status = AlgorithmStatus::AlgorithmSpecificFailure;
				result.interpolation_status = InterpolationStatus::SingularStencil;
				result.error_message = error.what();
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
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

		Real calcDetailedValue(int j, Real x) const override {
			if (X(j) == X(j + 1)) return Y(j);
			return Y(j) + ((x - X(j)) / (X(j + 1) - X(j)) * (Y(j + 1) - Y(j)));
		}

		const char* InterpolationMethodName() const override { return "Linear"; }
	};

} // namespace MML

#endif // MML_INTERPOLATED_REAL_FUNCTION_LINEAR_H