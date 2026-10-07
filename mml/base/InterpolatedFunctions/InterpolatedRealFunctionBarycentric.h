///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedRealFunctionBarycentric.h                               ///
///  Description: Barycentric polynomial and rational interpolation                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedRealFunctionBarycentric.h
/// @brief Barycentric polynomial and rational interpolation with Chebyshev helpers.
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_REAL_FUNCTION_BARYCENTRIC_H
#define MML_INTERPOLATED_REAL_FUNCTION_BARYCENTRIC_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/InterpolatedFunctions/InterpolationTypes.h>

namespace MML {
	/////////////////////////////////////////////////////////////////////////////////////
	///                    BARYCENTRIC POLYNOMIAL INTERPOLATION                       ///
	/////////////////////////////////////////////////////////////////////////////////////

	class BarycentricPolynomialInterp : public IRealFunction {
		int _numPoints;
		Vector<Real> _x;
		Vector<Real> _y;
		Vector<Real> _weights;
		ExtrapolationPolicy _extrapolationPolicy = ExtrapolationPolicy::Allow;

	public:
		BarycentricPolynomialInterp(const Vector<Real>& xv, const Vector<Real>& yv)
			: _numPoints(xv.size()), _x(xv.size()), _y(yv.size()), _weights(xv.size()) {
			if (xv.size() != yv.size() || xv.size() < 2)
				throw RealFuncInterpInitError("BarycentricPolynomialInterp: invalid data sizes");

			std::vector<std::pair<Real, Real>> data(_numPoints);
			for (int i = 0; i < _numPoints; ++i) data[i] = {xv[i], yv[i]};
			std::sort(data.begin(), data.end(), [](const auto& left, const auto& right) {
				return left.first < right.first;
			});
			for (int i = 0; i < _numPoints; ++i) {
				if (i > 0 && data[i].first == data[i - 1].first)
					throw RealFuncInterpInitError("BarycentricPolynomialInterp: duplicate x-value");
				_x[i] = data[i].first;
				_y[i] = data[i].second;
			}
			computeWeights();
		}

		Real operator()(Real x) const override {
			if (_extrapolationPolicy != ExtrapolationPolicy::Allow && !isInRange(x)) {
				if (_extrapolationPolicy == ExtrapolationPolicy::Throw)
					throw RealFuncInterpRuntimeError("BarycentricPolynomialInterp: query outside data range");
				x = std::clamp(x, MinX(), MaxX());
			}
			return evaluate(x);
		}

		InterpolationResult EvaluateDetailed(Real x, const InterpolationConfig& config = {}) const {
			AlgorithmTimer timer;
			InterpolationResult result;
			result.algorithm_name = "BarycentricPolynomial";
			result.query_point = x;
			result.evaluated_point = x;
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
			result.value = evaluate(result.evaluated_point);
			result.interval_index = locateInterval(result.evaluated_point);
			result.function_evaluations = 1;
			if (config.check_finite && !std::isfinite(result.value)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Barycentric polynomial interpolation produced a non-finite result";
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		Real MinX() const { return _x[0]; }
		Real MaxX() const { return _x[_numPoints - 1]; }
		int getNumPoints() const { return _numPoints; }
		bool isInRange(Real x) const { return x >= MinX() && x <= MaxX(); }
		void setExtrapolationPolicy(ExtrapolationPolicy policy) { _extrapolationPolicy = policy; }
		ExtrapolationPolicy getExtrapolationPolicy() const { return _extrapolationPolicy; }

		Real getWeight(int index) const {
			if (index < 0 || index >= _numPoints)
				throw IndexError("BarycentricPolynomialInterp weight index out of range");
			return _weights[index];
		}

	private:
		void computeWeights() {
			Real midpoint = REAL(0.5) * (MinX() + MaxX());
			Real scale = REAL(0.5) * (MaxX() - MinX());
			Real maxWeight = 0.0;
			for (int i = 0; i < _numPoints; ++i) {
				Real normalizedXi = (_x[i] - midpoint) / scale;
				Real denominator = 1.0;
				for (int j = 0; j < _numPoints; ++j) {
					if (i == j) continue;
					Real normalizedXj = (_x[j] - midpoint) / scale;
					denominator *= normalizedXi - normalizedXj;
				}
				_weights[i] = REAL(1.0) / denominator;
				maxWeight = std::max(maxWeight, std::abs(_weights[i]));
			}
			for (int i = 0; i < _numPoints; ++i) _weights[i] /= maxWeight;
		}

		Real evaluate(Real x) const {
			Real numerator = 0.0;
			Real denominator = 0.0;
			for (int i = 0; i < _numPoints; ++i) {
				Real difference = x - _x[i];
				if (difference == 0.0) return _y[i];
				Real term = _weights[i] / difference;
				numerator += term * _y[i];
				denominator += term;
			}
			return numerator / denominator;
		}

		int locateInterval(Real x) const {
			int upper = static_cast<int>(std::upper_bound(&_x[0], &_x[0] + _numPoints, x) - &_x[0]);
			return std::clamp(upper - 1, 0, _numPoints - 2);
		}
	};

	inline Vector<Real> ChebyshevInterpolationNodes(Real min, Real max, int numPoints) {
		if (numPoints < 2 || min >= max)
			throw RealFuncInterpInitError("ChebyshevInterpolationNodes: require numPoints >= 2 and min < max");
		Vector<Real> nodes(numPoints);
		Real midpoint = REAL(0.5) * (min + max);
		Real halfWidth = REAL(0.5) * (max - min);
		for (int i = 0; i < numPoints; ++i) {
			Real angle = Constants::PI * (REAL(2.0) * i + REAL(1.0)) / (REAL(2.0) * numPoints);
			nodes[i] = midpoint - halfWidth * std::cos(angle);
		}
		return nodes;
	}

	inline BarycentricPolynomialInterp MakeChebyshevNodeInterpolator(
		const std::function<Real(Real)>& function, Real min, Real max, int numPoints) {
		Vector<Real> nodes = ChebyshevInterpolationNodes(min, max, numPoints);
		Vector<Real> values(numPoints);
		for (int i = 0; i < numPoints; ++i) values[i] = function(nodes[i]);
		return BarycentricPolynomialInterp(nodes, values);
	}

	inline BarycentricPolynomialInterp MakeChebyshevNodeInterpolator(
		const IRealFunction& function, Real min, Real max, int numPoints) {
		return MakeChebyshevNodeInterpolator(
			[&function](Real x) { return function(x); }, min, max, numPoints);
	}


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

		Real MinX() const { return std::min(_x[0], _x[_n - 1]); }

		/// /** @brief Get the maximum x-value. */

		Real MaxX() const { return std::max(_x[0], _x[_n - 1]); }

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

		InterpolationResult EvaluateDetailed(Real x, const InterpolationConfig& config = {}) const {
			AlgorithmTimer timer;
			InterpolationResult result;
			result.algorithm_name = "FloaterHormannRational";
			result.query_point = x;
			result.evaluated_point = x;
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
			result.value = evaluate(result.evaluated_point);
			result.interval_index = locateInterval(result.evaluated_point);
			result.function_evaluations = 1;
			if (config.check_finite && !std::isfinite(result.value)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Interpolation produced a non-finite result";
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
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

	private:
		Real evaluate(Real x) const {
			Real num = 0.0, den = 0.0;
			for (int i = 0; i < _n; ++i) {
				Real h = x - _x[i];
				if (h == 0.0) return _y[i];
				Real temp = _w[i] / h;
				num += temp * _y[i];
				den += temp;
			}
			return num / den;
		}

		int locateInterval(Real x) const {
			int lower = 0;
			int upper = _n - 1;
			const bool ascending = _x[_n - 1] >= _x[0];
			while (upper - lower > 1) {
				int middle = (lower + upper) / 2;
				if ((x >= _x[middle]) == ascending) lower = middle;
				else upper = middle;
			}
			return std::clamp(lower, 0, _n - 2);
		}
	};

} // namespace MML

#endif // MML_INTERPOLATED_REAL_FUNCTION_BARYCENTRIC_H