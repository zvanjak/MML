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
#include <limits>
#include <utility>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/InterpolatedFunctions/InterpolationTypes.h>

namespace MML {
	namespace detail {
		// Floater-Hormann window terms have a common sign at each node. Log-sum-exp
		// and a common normalization avoid overflowing their products (Floater & Hormann, 2007).
		inline Vector<Real> ComputeNormalizedBarycentricWeights(const Vector<Real>& nodes, int d) {
			const int count = nodes.size();
			Vector<Real> weights(count);
			if (count == 1) {
				weights[0] = 1;
				return weights;
			}
			const bool ascending = nodes[count - 1] > nodes[0];
			const Real span = std::abs(nodes[count - 1] - nodes[0]);
			const Real scale = std::isfinite(span) ? span : std::max(std::abs(nodes[0]), std::abs(nodes[count - 1]));
			Real maximumLog = -std::numeric_limits<Real>::infinity();
			for (int k = 0; k < count; ++k) {
				Real logWeight = -std::numeric_limits<Real>::infinity();
				for (int i = std::max(0, k - d); i <= std::min(k, count - d - 1); ++i) {
					Real logTerm = 0;
					for (int j = i; j <= i + d; ++j) {
						if (j == k) continue;
						const Real difference = nodes[k] - nodes[j];
						const Real normalized = std::isfinite(difference) ? std::abs(difference) / scale
							: std::abs(nodes[k] / scale - nodes[j] / scale);
						logTerm -= normalized > 0 ? std::log(normalized)
							: std::log(std::abs(difference)) - std::log(scale);
					}
					if (logWeight == -std::numeric_limits<Real>::infinity())
						logWeight = logTerm;
					else {
						const Real larger = std::max(logWeight, logTerm);
						logWeight = larger + std::log1p(std::exp(std::min(logWeight, logTerm) - larger));
					}
				}
				weights[k] = logWeight;
				maximumLog = std::max(maximumLog, logWeight);
			}
			for (int k = 0; k < count; ++k) {
				const Real magnitude = std::exp(weights[k] - maximumLog);
				if (!std::isfinite(magnitude) || magnitude == 0)
					throw RealFuncInterpInitError("Barycentric interpolation: relative weight range is not representable");
				const bool negative = ascending ? (d % 2 != k % 2) : (k % 2 != 0);
				weights[k] = negative ? -magnitude : magnitude;
			}
			return weights;
		}

		inline Real EvaluateBarycentric(const Vector<Real>& nodes, const Vector<Real>& values,
			const Vector<Real>& weights, Real x) {
			if (!std::isfinite(x))
				throw RealFuncInterpRuntimeError("Barycentric interpolation: query must be finite");
			int nearest = 0;
			Real distance = std::abs(x - nodes[0]);
			Real valueScale = 0;
			for (int i = 0; i < nodes.size(); ++i) {
				if (x == nodes[i]) return values[i];
				const Real current = std::abs(x - nodes[i]);
				if (current < distance) {
					distance = current;
					nearest = i;
				}
			}
			for (int i = 0; i < values.size(); ++i) {
				if (!std::isfinite(values[i]))
					throw RealFuncInterpRuntimeError("Barycentric interpolation: ordinates must be finite");
				valueScale = std::max(valueScale, std::abs(values[i]));
			}
			if (valueScale == 0) return 0;
			Real halfDistance = std::abs(x / 2 - nodes[nearest] / 2);
			if (!std::isfinite(distance)) {
				for (int i = 1; i < nodes.size(); ++i)
					halfDistance = std::min(halfDistance, std::abs(x / 2 - nodes[i] / 2));
			}
			Real numerator = 0, denominator = 0, numeratorCorrection = 0, denominatorCorrection = 0;
			auto accumulate = [](Real term, Real& sum, Real& correction) {
				const Real next = sum + term;
				correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
				sum = next;
			};
			for (int i = 0; i < nodes.size(); ++i) {
				const Real difference = x - nodes[i];
				const Real ratio = std::isfinite(difference) ? distance / difference
					: (std::isfinite(distance) ? distance / 2 : halfDistance) / (x / 2 - nodes[i] / 2);
				const Real term = weights[i] * ratio;
				accumulate(term * (values[i] / valueScale), numerator, numeratorCorrection);
				accumulate(term, denominator, denominatorCorrection);
			}
			const Real total = denominator + denominatorCorrection;
			if (!std::isfinite(total) || total == 0)
				throw RealFuncInterpRuntimeError("Barycentric interpolation: vanishing or nonfinite denominator");
			const Real result = ((numerator + numeratorCorrection) / total) * valueScale;
			if (!std::isfinite(result))
				throw RealFuncInterpRuntimeError("Barycentric interpolation: result is not representable");
			return result;
		}
	}

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
			for (int i = 0; i < _numPoints; ++i) {
				if (!std::isfinite(xv[i]))
					throw RealFuncInterpInitError("BarycentricPolynomialInterp: nodes must be finite");
				data[i] = {xv[i], yv[i]};
			}
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
			result.function_evaluations = 1;
			try {
				result.value = evaluate(result.evaluated_point);
			} catch (const RealFuncInterpRuntimeError& error) {
				if (config.exception_policy == EvaluationExceptionPolicy::Propagate) throw;
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = error.what();
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}
			result.interval_index = locateInterval(result.evaluated_point);
			result.status = AlgorithmStatus::Completed;
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

		/// Relative weight; the largest absolute weight is normalized to one.
		Real getWeight(int index) const {
			if (index < 0 || index >= _numPoints)
				throw IndexError("BarycentricPolynomialInterp weight index out of range");
			return _weights[index];
		}

	private:
		void computeWeights() {
			_weights = detail::ComputeNormalizedBarycentricWeights(_x, _numPoints - 1);
		}

		Real evaluate(Real x) const {
			return detail::EvaluateBarycentric(_x, _y, _weights, x);
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
		Real midpoint = static_cast<MML::Real>(0.5) * (min + max);
		Real halfWidth = static_cast<MML::Real>(0.5) * (max - min);
		for (int i = 0; i < numPoints; ++i) {
			Real angle = Constants::PI * (static_cast<MML::Real>(2.0) * i + static_cast<MML::Real>(1.0)) / (static_cast<MML::Real>(2.0) * numPoints);
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

		static int ValidateData(const Vector<Real>& x, const Vector<Real>& y, int d) {
			if (x.size() != y.size())
				throw RealFuncInterpInitError("BarycentricRationalInterp: x and y must have the same size");
			if (d < 0 || d >= x.size())
				throw RealFuncInterpInitError("BarycentricRationalInterp: d must be in [0, number of points)");
			detail::ValidateInterpolationNodes(x, 1, "BarycentricRationalInterp");
			return x.size();
		}

		/// /** @brief Compute the barycentric weights using Floater-Hormann formula. */

		void computeWeights() {
			_w = detail::ComputeNormalizedBarycentricWeights(_x, _d);
		}

	public:
		/// @brief Construct a barycentric rational interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param d Blending parameter (0 to n-1), default 3
		/// @throws RealFuncInterpInitError for mismatched sizes, invalid d, or nonfinite/nonmonotonic nodes
		/// @note A single node with d=0 defines a constant interpolant.

		BarycentricRationalInterp(const Vector<Real>& xv, const Vector<Real>& yv, int d = 3)
			: _n(ValidateData(xv, yv, d))
			, _d(d)
			, _x(xv)
			, _y(yv)
			, _w(_n)
			, _extrapolationPolicy(ExtrapolationPolicy::Allow) {
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
			return evaluate(x);
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
			result.function_evaluations = 1;
			try {
				result.value = evaluate(result.evaluated_point);
			} catch (const RealFuncInterpRuntimeError& error) {
				if (config.exception_policy == EvaluationExceptionPolicy::Propagate) throw;
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = error.what();
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			}
			result.interval_index = locateInterval(result.evaluated_point);
			result.status = AlgorithmStatus::Completed;
			if (config.check_finite && !std::isfinite(result.value)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "Interpolation produced a non-finite result";
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		/// @brief Get the barycentric weight at index i.
		/// @param i Index of the data point
		/// @return Relative barycentric weight w[i]; largest absolute weight is one
		/// @throws IndexError if i is out of bounds

		Real getWeight(int i) const {
			if (i < 0 || i >= _n)
				throw IndexError("Index out of range in getWeight");
			return _w[i];
		}

	private:
		Real evaluate(Real x) const {
			if (_n == 1)
				return _y[0];
			return detail::EvaluateBarycentric(_x, _y, _w, x);
		}

		int locateInterval(Real x) const {
			if (_n == 1)
				return -1;
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
