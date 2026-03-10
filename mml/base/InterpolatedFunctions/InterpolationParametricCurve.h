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
/// @see InterpolatedRealFunction.h for 1D interpolation
/// @see Interpolation2DFunction.h for 2D grid interpolation
/// @ingroup Interpolation

#if !defined MML_INTERPOLATION_PARAMETRIC_CURVE_H
#define MML_INTERPOLATION_PARAMETRIC_CURVE_H

#include <cmath>
#include <vector>

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include "base/InterpolatedFunctions/InterpolatedRealFunction.h"

namespace MML {

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

	public:
		/// @brief Construct a linearly interpolated parametric curve.
		/// @param ptsin Matrix of control points [numPoints x N]
		/// @param close If true, the curve is closed (connects last to first point)

		LinInterpParametricCurve(const Matrix<Real>& ptsin, bool close = false)
			: _numPoints(ptsin.rows())
			, _curvePoints(ptsin)
			, _arcLen(_numPoints)
			, _isCurveClosed(close) {
			// Compute normalized arc length parameterization
			_arcLen[0] = 0.0;
			Real totalLen = 0.0;
			for (int i = 1; i < _numPoints; ++i) {
				Real segLen = 0.0;
				for (int j = 0; j < N; ++j)
					segLen += std::pow(_curvePoints[i][j] - _curvePoints[i - 1][j], 2);
				totalLen += std::sqrt(segLen);
				_arcLen[i] = totalLen;
			}
			// Normalize to [0,1]
			for (int i = 0; i < _numPoints; ++i)
				_arcLen[i] /= totalLen > 0 ? totalLen : 1.0;
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
			Real alpha = (t1 - t0 > 0) ? (t - t0) / (t1 - t0) : 0.0;

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

		std::vector<SplineInterpRealFunc*> srp; ///< Spline for each coordinate

	public:
		/// @brief Construct a spline-interpolated parametric curve.
		/// @param minT Minimum parameter value
		/// @param maxT Maximum parameter value
		/// @param ptsin Matrix of control points [numPoints x dim]
		/// @param close If true, curve is closed (last connects to first)
		/// @throws If N doesn't match matrix column count
		/// @note For closed curves, do not duplicate the first point as the last.

		SplineInterpParametricCurve(Real minT, Real maxT, const Matrix<Real>& ptsin, bool close = 0)
			: _numPoints(ptsin.rows())
			, _dim(ptsin.cols())
			, _bemba(close ? 2 * _numPoints : _numPoints)
			, _isCurveClosed(close)
			, _curvePoints(_dim, _bemba)
			, s(_bemba)
			, ans(_dim)
			, srp(_dim)
			, _minT(minT)
			, _maxT(maxT) {
			// check N == dim
			if (N != _dim)
				throw MatrixDimensionError("SplineInterpParametricCurve: template dimension N must match point matrix column count");

			int i, ii, im, j, ofs;
			Real ss, soff, db, de;

			ofs = close ? _numPoints / 2 : 0;
			s[0] = 0.;
			for (i = 0; i < _bemba; i++) {
				ii = (i - ofs + _numPoints) % _numPoints;
				im = (ii - 1 + _numPoints) % _numPoints;

				for (j = 0; j < _dim; j++)
					_curvePoints[j][i] = ptsin[ii][j];

				if (i > 0) {
					s[i] = s[i - 1] + rad(&ptsin[ii][0], &ptsin[im][0]);

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

			for (j = 0; j < _dim; j++) {
				db = _bemba < 4 ? 1.e99 : fprime(&s[0], &_curvePoints[j][0], 1);
				de = _bemba < 4 ? 1.e99 : fprime(&s[_bemba - 1], &_curvePoints[j][_bemba - 1], -1);

				Vector<Real> vec = _curvePoints.VectorFromRow(j);

				srp[j] = new SplineInterpRealFunc(s, vec, db, de);
			}
		}

		/// @brief Construct with default parameter range [0, 1].
		/// @param ptsin Matrix of control points
		/// @param close If true, curve is closed

		SplineInterpParametricCurve(const Matrix<Real>& ptsin, bool close = 0)
			: SplineInterpParametricCurve(0.0, 1.0, ptsin, close) {}

		~SplineInterpParametricCurve() {
			for (int j = 0; j < _dim; j++)
				delete srp[j];
		}

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

			if (_isCurveClosed)
				t = t - floor(t);

			// we have to map t from [minT, maxT] to [0, 1]
			t = (t - _minT) / (_maxT - _minT);
			for (int j = 0; j < _dim; j++)
				ans[j] = (*srp[j])(t);

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
				   (s2 * s3 / (s1 * s12 * s13)) * y[pm * 1] + (1. / s1 + 1. / s2 + 1. / s3) * y[0];
		}

		/// @brief Compute Euclidean distance between two points.
		/// @param p1 First point
		/// @param p2 Second point
		/// @return Distance ||p1 - p2||

		Real rad(const Real* p1, const Real* p2) {
			Real sum = 0.;
			for (int i = 0; i < _dim; i++)
				sum += POW2(p1[i] - p2[i]);
			return sqrt(sum);
		}
	};

} // namespace MML

#endif // MML_INTERPOLATION_PARAMETRIC_CURVE_H
