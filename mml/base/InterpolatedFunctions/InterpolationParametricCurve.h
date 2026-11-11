///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolationParametricCurve.h                                      ///
///  Description: Parametric curve interpolation classes                              ///
///               Linear and spline interpolation for N-dimensional curves           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolationParametricCurve.h
/// @brief Parametric curve interpolation classes for N-dimensional curves.
/// @section interp_param_classes Classes
/// - LinInterpParametricCurve<N> - Linear interpolation of N-dimensional curves
/// - SplineInterpParametricCurve<N> - Spline interpolation of N-dimensional curves
/// @see InterpolatedFunctionSpline.h for spline interpolation
/// @see Interpolation2DFunction.h for 2D grid interpolation
/// @ingroup Interpolation

#if !defined MML_INTERPOLATION_PARAMETRIC_CURVE_H
#define MML_INTERPOLATION_PARAMETRIC_CURVE_H

#include <cmath>
#include <vector>
#include <limits>

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <mml/base/InterpolatedFunctions/InterpolatedFunctionSpline.h>

namespace MML {
	namespace detail {
		template<int N>
		int ValidateParametricCurvePoints(const Matrix<Real>& points, Real minT = 0, Real maxT = 1, bool extendClosed = false) {
			if (points.cols() != N)
				throw MatrixDimensionError("Parametric curve: template dimension must match point matrix column count");
			if (points.rows() < 2)
				throw RealFuncInterpInitError("Parametric curve: at least two control points are required");
			if (!std::isfinite(minT) || !std::isfinite(maxT) || !(minT < maxT) || !std::isfinite(maxT - minT))
				throw RealFuncInterpInitError("Parametric curve: parameter range must be finite, increasing and representable");
			if (extendClosed && points.rows() > std::numeric_limits<int>::max() / 2)
				throw RealFuncInterpInitError("Parametric curve: extended closed-curve point count exceeds integer range");
			for (int row = 0; row < points.rows(); ++row)
				for (int col = 0; col < N; ++col)
					if (!std::isfinite(points(row, col)))
						throw RealFuncInterpInitError("Parametric curve: control points must be finite");
			return points.rows();
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                    PARAMETRIC CURVE INTERPOLATION                             ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Linear interpolation of an N-dimensional parametric curve.
	/// LinInterpParametricCurve interpolates a curve in N-dimensional space
	/// given as a sequence of control points. The curve is parameterized by
	/// normalized arc length in [0, 1].
	/// **Parameterization:** Uses chord-length parameterization where the
	/// parameter t is proportional to the cumulative Euclidean distance
	/// along the polyline.
	/// **Features:**
	/// - Works for any dimension N (2D, 3D, or higher)
	/// - Supports both open and closed curves
	/// - Simple and fast (linear segments)
	/// @tparam N Dimension of the space (e.g., 2, 3)
	/// @see SplineInterpParametricCurve for smoother curve interpolation
	/// @see IParametricCurve Base interface for parametric curves
	/// @ingroup Interpolation

	template<int N>
	class LinInterpParametricCurve : public IParametricCurve<N> {
		int _numPoints;			   ///< Number of control points
		Matrix<Real> _curvePoints; ///< Control points [numPoints x N]
		Vector<Real> _arcLen;	   ///< Normalized arc length parameter [0,1]
		bool _isCurveClosed;	   ///< True for closed curves

		static int PointCount(const Matrix<Real>& points, bool close) {
			const int count = detail::ValidateParametricCurvePoints<N>(points);
			if (close && count == std::numeric_limits<int>::max())
				throw RealFuncInterpInitError("Linear parametric curve: closed point count exceeds integer range");
			return count + (close ? 1 : 0);
		}

	public:
		/// @brief Construct a linearly interpolated parametric curve.
		/// @param ptsin Matrix of control points [numPoints x N]
		/// @param close If true, the curve is closed (connects last to first point)
		/// @throws MatrixDimensionError if the matrix does not have exactly N columns
		/// @throws RealFuncInterpInitError for fewer than two points, nonfinite coordinates or unrepresentable chord lengths

		LinInterpParametricCurve(const Matrix<Real>& ptsin, bool close = false)
			: _numPoints(PointCount(ptsin, close))
			, _curvePoints(_numPoints, N)
			, _arcLen(_numPoints)
			, _isCurveClosed(close) {
			for (int i = 0; i < _numPoints; ++i)
				for (int j = 0; j < N; ++j)
					_curvePoints(i, j) = ptsin(i % ptsin.rows(), j);
			// Compute normalized arc length parameterization
			_arcLen[0] = 0.0;
			Real totalLen = 0.0;
			for (int i = 1; i < _numPoints; ++i) {
				Real segLen = 0.0;
				for (int j = 0; j < N; ++j)
					segLen = std::hypot(segLen, _curvePoints[i][j] - _curvePoints[i - 1][j]);
				totalLen += segLen;
				if (!std::isfinite(totalLen))
					throw RealFuncInterpInitError("Linear parametric curve: chord length is not representable");
				_arcLen[i] = totalLen;
			}
			// Normalize to [0,1]
			for (int i = 0; i < _numPoints; ++i)
				_arcLen[i] /= totalLen > 0 ? totalLen : MML::Real(1.0);
		}

		/// /** @brief Get the minimum parameter value (0). */

		Real getMinT() const override { return 0.0; }

		/// /** @brief Get the maximum parameter value (1). */

		Real getMaxT() const override { return 1.0; }

		/// @brief Evaluate the curve at parameter t.
		/// @param t Parameter value in [0, 1]
		/// @return Point on the curve at parameter t
		/// For open curves, t is clamped to [0, 1].
		/// For closed curves, t wraps periodically.

		VectorN<Real, N> operator()(Real t) const override {
			// Clamp or wrap t
			if (_isCurveClosed) {
				t = t - std::floor(t);
			} else {
				if (t <= 0.0)
					t = 0.0;
				if (t >= 1.0)
					t = 1.0;
			}

			// Find segment
			int i = 0;
			while (i < _numPoints - 2 && t > _arcLen[i + 1])
				++i;

			// Linear interpolation
			Real t0 = _arcLen[i], t1 = _arcLen[i + 1];
			Real alpha = (t1 - t0 > 0) ? (t - t0) / (t1 - t0) : MML::Real(0.0);

			VectorN<Real, N> result;
			for (int j = 0; j < N; ++j)
				result[j] = (1 - alpha) * _curvePoints[i][j] + alpha * _curvePoints[i + 1][j];
			return result;
		}
	};


	/// @brief Spline interpolation of an N-dimensional parametric curve.
	/// SplineInterpParametricCurve provides smooth interpolation of a curve
	/// in N-dimensional space using cubic splines for each coordinate.
	/// The curve is parameterized by normalized arc length.
	/// **Algorithm:** Constructs independent cubic splines for each coordinate
	/// dimension, all parameterized by the same arc-length parameter.
	/// **Features:**
	/// - C² continuous (smooth curves)
	/// - Supports both open and closed curves
	/// - Automatic derivative estimation at endpoints
	/// - Custom parameter range [minT, maxT] or default [0, 1]
	/// @note For closed curves, do not duplicate the first point as the last.
	/// The algorithm will connect them automatically.
	/// @tparam N Dimension of the space
	/// @see LinInterpParametricCurve for faster piecewise linear alternative
	/// @see SplineInterpRealFunc for 1D spline details
	/// @ingroup Interpolation

	template<int N>
	class SplineInterpParametricCurve : public IParametricCurve<N> {
		Real _minT, _maxT;	  ///< Parameter range
		int _dim, _numPoints; ///< Dimension and number of input points
		int _bemba;			  ///< Extended point count for closed curves
		bool _isCurveClosed;  ///< True for closed curves

		Matrix<Real> _curvePoints; ///< Rearranged control points [dim x bemba]
		Vector<Real> s;			   ///< Arc-length parameters
		Vector<Real> ans;		   ///< Temporary result storage

		std::vector<SplineInterpRealFunc> srp; ///< Spline for each coordinate

	public:
		/// @brief Construct a spline-interpolated parametric curve.
		/// @param minT Minimum parameter value
		/// @param maxT Maximum parameter value
		/// @param ptsin Matrix of control points [numPoints x dim]
		/// @param close If true, curve is closed (last connects to first)
		/// @throws MatrixDimensionError if the matrix does not have exactly N columns
		/// @throws RealFuncInterpInitError for invalid point counts, nonfinite coordinates or invalid parameter bounds
		/// @note For closed curves, do not duplicate the first point as the last.

		SplineInterpParametricCurve(Real minT, Real maxT, const Matrix<Real>& ptsin, bool close = 0)
			: _minT(minT)
			, _maxT(maxT)
			, _dim(ptsin.cols())
			, _numPoints(detail::ValidateParametricCurvePoints<N>(ptsin, minT, maxT, close))
			, _bemba(close ? 2 * _numPoints : _numPoints)
			, _isCurveClosed(close)
			, _curvePoints(_dim, _bemba)
			, s(_bemba)
			, ans(_dim) {
			int i, ii, im, j, ofs;
			Real ss, soff, db, de;

			ofs = close ? _numPoints / 2 : 0;
			s[0] = 0.;
			for (i = 0; i < _bemba; i++) {
				// Wrap without forming i + _numPoints, which can overflow for large inputs.
				ii = i % _numPoints;
				ii = ii < ofs ? _numPoints - (ofs - ii) : ii - ofs;
				im = ii == 0 ? _numPoints - 1 : ii - 1;

				for (j = 0; j < _dim; j++)
					_curvePoints[j][i] = ptsin[ii][j];

				if (i > 0) {
					s[i] = s[i - 1] + rad(&ptsin[ii][0], &ptsin[im][0]);
					if (!std::isfinite(s[i]))
						throw RealFuncInterpInitError("Spline parametric curve: chord length is not representable");

					if (s[i] == s[i - 1])
						throw RealFuncInterpRuntimeError("Parametric curve interpolation: consecutive identical points not allowed");
					// Consecutive points may not be identical. For a closed curve, the last data
					// point should not duplicate the first.
				}
			}
			ss = close ? s[ofs + _numPoints] - s[ofs] : s[_numPoints - 1] - s[0];
			soff = s[ofs];

			for (i = 0; i < _bemba; i++)
				s[i] = (s[i] - soff) / ss;

			srp.reserve(_dim);
			for (j = 0; j < _dim; j++) {
				db = _bemba < 4 ? MML::Real(1.e99) : fprime(&s[0], &_curvePoints[j][0], 1);
				de = _bemba < 4 ? MML::Real(1.e99) : fprime(&s[_bemba - 1], &_curvePoints[j][_bemba - 1], -1);

				Vector<Real> vec = _curvePoints.VectorFromRow(j);

				srp.emplace_back(s, vec, db, de);
			}
		}

		/// @brief Construct with default parameter range [0, 1].
		/// @param ptsin Matrix of control points
		/// @param close If true, curve is closed

		SplineInterpParametricCurve(const Matrix<Real>& ptsin, bool close = 0)
			: SplineInterpParametricCurve(0.0, 1.0, ptsin, close) {}

		/// /** @brief Get the minimum parameter value. */

		Real getMinT() const { return _minT; }

		/// /** @brief Get the maximum parameter value. */

		Real getMaxT() const { return _maxT; }

		/// @brief Evaluate the curve at parameter t.
		/// @param t Parameter value in [minT, maxT]
		/// @return Point on the curve at parameter t
		/// @throws For closed curves, if t is outside [minT, maxT]
		/// For open curves, values outside [minT, maxT] extrapolate (use with caution).
		/// For closed curves, t is treated as periodic.

		VectorN<Real, N> operator()(Real t) const {
			if (_isCurveClosed == true && (t < _minT || t > _maxT))
				throw DomainError("SplineInterpParametricCurve: parameter t is outside [minT, maxT] for closed curve");

			VectorN<Real, N> ans;

			// we have to map t from [minT, maxT] to [0, 1]
			t = (t - _minT) / (_maxT - _minT);
			if (_isCurveClosed)
				t = t - std::floor(t);
			for (int j = 0; j < _dim; j++)
				ans[j] = srp[j](t);

			return ans;
		}

	private:
		/// @brief Estimate derivative at an endpoint using neighboring points.
		/// @param x Pointer to arc-length parameters
		/// @param y Pointer to coordinate values
		/// @param pm Direction: +1 for left endpoint, -1 for right endpoint
		/// @return Estimated first derivative

		Real fprime(Real* x, Real* y, int pm) {
			Real s1 = x[0] - x[pm * 1], s2 = x[0] - x[pm * 2], s3 = x[0] - x[pm * 3], s12 = s1 - s2, s13 = s1 - s3, s23 = s2 - s3;
			return -(s1 * s2 / (s13 * s23 * s3)) * y[pm * 3] + (s1 * s3 / (s12 * s2 * s23)) * y[pm * 2] -
				   (s2 * s3 / (s1 * s12 * s13)) * y[pm * 1] + (MML::Real(1.) / s1 + MML::Real(1.) / s2 + MML::Real(1.) / s3) * y[0];
		}

		/// @brief Compute Euclidean distance between two points.
		/// @param p1 First point
		/// @param p2 Second point
		/// @return Distance ||p1 - p2||

		Real rad(const Real* p1, const Real* p2) {
			Real sum = 0.;
			for (int i = 0; i < _dim; i++)
				sum = std::hypot(sum, p1[i] - p2[i]);
			return sum;
		}
	};

} // namespace MML

#endif // MML_INTERPOLATION_PARAMETRIC_CURVE_H
