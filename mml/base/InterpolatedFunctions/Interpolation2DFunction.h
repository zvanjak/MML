///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Interpolation2DFunction.h                                           ///
///  Description: 2D interpolation classes for grid-based function approximation     ///
///               Bilinear and bicubic spline interpolation on regular grids         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Interpolation2DFunction.h
/// @brief 2D interpolation classes for function approximation on regular grids.
/// @section interp2d_classes Classes
/// - BilinearInterp2D - Bilinear interpolation on regular grids
/// - BicubicSplineInterp2D - Bicubic spline interpolation on regular grids
/// @see InterpolatedRealFunction.h for 1D interpolation
/// @see InterpolationParametricCurve.h for parametric curve interpolation
/// @ingroup Interpolation

#if !defined MML_INTERPOLATION_2D_FUNCTION_H
#define MML_INTERPOLATION_2D_FUNCTION_H

#include <vector>

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include "base/InterpolatedFunctions/InterpolatedRealFunction.h"

namespace MML {

	/////////////////////////////////////////////////////////////////////////////////////
	///                       2D INTERPOLATION CLASSES                                ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Bilinear interpolation on a regular 2D grid.
	/// BilinearInterp2D interpolates values on a 2D rectangular grid using
	/// bilinear interpolation. Given function values z(x1[i], x2[j]) at grid
	/// points, it computes z at any interior point.
	/// **Formula:** For a unit square with corners (0,0), (1,0), (0,1), (1,1):
	/// @f[
	/// z(t, u) = (1-t)(1-u)z_{00} + t(1-u)z_{10} + (1-t)u z_{01} + tu z_{11}
	/// @f]
	/// The method also provides partial derivatives `dz/dx1` and `dz/dx2`.
	/// **Continuity:** C⁰ (continuous but with discontinuous derivatives at grid lines)
	/// @see BicubicSplineInterp2D for smoother interpolation with continuous derivatives
	/// @ingroup Interpolation

	class BilinearInterp2D {
	private:
		int _m, _n;				///< Grid dimensions (m x n)
		Vector<Real> _x1, _x2;	///< Grid coordinates
		const Matrix<Real>& _z; ///< Reference to data matrix z[m][n]

		/// @brief Binary search to find interval containing x.
		/// @param arr Array to search
		/// @param x Value to locate
		/// @return Index i such that arr[i] <= x < arr[i+1]

		int locate(const Vector<Real>& arr, Real x) const {
			int n = arr.size();
			int lo = 0, hi = n - 1;
			bool ascending = (arr[n - 1] >= arr[0]);

			while (hi - lo > 1) {
				int mid = (hi + lo) / 2;
				if ((x >= arr[mid]) == ascending)
					lo = mid;
				else
					hi = mid;
			}
			return std::max(0, std::min(n - 2, lo));
		}

	public:
		/// @brief Construct a bilinear interpolation object.
		/// @param x1v Vector of x1 grid coordinates (m values)
		/// @param x2v Vector of x2 grid coordinates (n values)
		/// @param zm Matrix of function values z[m][n]
		/// @throws RealFuncInterpInitError if dimensions don't match

		BilinearInterp2D(const Vector<Real>& x1v, const Vector<Real>& x2v, const Matrix<Real>& zm)
			: _m(x1v.size())
			, _n(x2v.size())
			, _x1(x1v)
			, _x2(x2v)
			, _z(zm) {
			if (zm.rows() != _m || zm.cols() != _n)
				throw RealFuncInterpInitError("BilinearInterp2D: matrix dimensions don't match grid");
		}

		/// @brief Evaluate the bilinear interpolation at (x1, x2).
		/// @param x1 First coordinate
		/// @param x2 Second coordinate
		/// @return Interpolated value z(x1, x2)

		Real operator()(Real x1, Real x2) const {
			int i = locate(_x1, x1);
			int j = locate(_x2, x2);

			// Compute normalized coordinates t, u in [0, 1]
			Real t = (_x1[i + 1] - _x1[i] != 0.0) ? (x1 - _x1[i]) / (_x1[i + 1] - _x1[i]) : 0.0;
			Real u = (_x2[j + 1] - _x2[j] != 0.0) ? (x2 - _x2[j]) / (_x2[j + 1] - _x2[j]) : 0.0;

			// Bilinear interpolation formula
			return (1.0 - t) * (1.0 - u) * _z(i, j) + t * (1.0 - u) * _z(i + 1, j) + (1.0 - t) * u * _z(i, j + 1) +
				   t * u * _z(i + 1, j + 1);
		}

		/// @brief Evaluate interpolation and partial derivatives at (x1, x2).
		/// @param x1 First coordinate
		/// @param x2 Second coordinate
		/// @param[out] z Interpolated value
		/// @param[out] dz_dx1 Partial derivative ∂z/∂x1
		/// @param[out] dz_dx2 Partial derivative ∂z/∂x2

		void interpWithDerivatives(Real x1, Real x2, Real& z, Real& dz_dx1, Real& dz_dx2) const {
			int i = locate(_x1, x1);
			int j = locate(_x2, x2);

			Real dx1 = _x1[i + 1] - _x1[i];
			Real dx2 = _x2[j + 1] - _x2[j];

			Real t = (dx1 != 0.0) ? (x1 - _x1[i]) / dx1 : 0.0;
			Real u = (dx2 != 0.0) ? (x2 - _x2[j]) / dx2 : 0.0;

			// Function value
			z = (1.0 - t) * (1.0 - u) * _z(i, j) + t * (1.0 - u) * _z(i + 1, j) + (1.0 - t) * u * _z(i, j + 1) + t * u * _z(i + 1, j + 1);

			// Partial derivatives
			if (dx1 != 0.0) {
				dz_dx1 = ((1.0 - u) * (_z(i + 1, j) - _z(i, j)) + u * (_z(i + 1, j + 1) - _z(i, j + 1))) / dx1;
			} else {
				dz_dx1 = 0.0;
			}

			if (dx2 != 0.0) {
				dz_dx2 = ((1.0 - t) * (_z(i, j + 1) - _z(i, j)) + t * (_z(i + 1, j + 1) - _z(i + 1, j))) / dx2;
			} else {
				dz_dx2 = 0.0;
			}
		}

		/// @name Grid Information
		/// @{

		/// /** @brief Get the first grid dimension (number of x1 values). */

		int getGridDim1() const { return _m; }

		/// /** @brief Get the second grid dimension (number of x2 values). */

		int getGridDim2() const { return _n; }

		/// @}
	};


	/// @brief Bicubic spline interpolation on a regular 2D grid.
	/// BicubicSplineInterp2D provides smooth 2D interpolation using cubic splines
	/// in both coordinate directions. It offers continuous first derivatives
	/// (C¹ continuity) and is more accurate than bilinear for smooth data.
	/// **Algorithm:** Constructs 1D cubic splines along each row (fixed x1),
	/// then interpolates along x2 first and along x1 second.
	/// **Continuity:** C¹ (continuous first derivatives)
	/// @note Uses dynamically allocated SplineInterpRealFunc objects internally.
	/// Copy and assignment are disabled to prevent double-free.
	/// @see BilinearInterp2D for faster but rougher alternative
	/// @see SplineInterpRealFunc for 1D cubic spline interpolation
	/// @ingroup Interpolation

	class BicubicSplineInterp2D {
	private:
		int _m, _n;										///< Grid dimensions
		Vector<Real> _x1, _x2;							///< Grid coordinates
		Matrix<Real> _z;								///< Local copy of data
		std::vector<SplineInterpRealFunc*> _rowSplines; ///< Splines along x2 for each row

	public:
		/// @brief Construct a bicubic spline interpolation object.
		/// @param x1v Vector of x1 grid coordinates (m values)
		/// @param x2v Vector of x2 grid coordinates (n values)
		/// @param zm Matrix of function values z[m][n]
		/// @throws RealFuncInterpInitError if dimensions don't match

		BicubicSplineInterp2D(const Vector<Real>& x1v, const Vector<Real>& x2v, const Matrix<Real>& zm)
			: _m(x1v.size())
			, _n(x2v.size())
			, _x1(x1v)
			, _x2(x2v)
			, _z(zm)
			, _rowSplines(_m) {
			if (zm.rows() != _m || zm.cols() != _n)
				throw RealFuncInterpInitError("BicubicSplineInterp2D: matrix dimensions don't match grid");

			// Create a spline along x2 for each row (fixed x1[i])
			for (int i = 0; i < _m; i++) {
				Vector<Real> row(_n);
				for (int j = 0; j < _n; j++)
					row[j] = _z(i, j);
				_rowSplines[i] = new SplineInterpRealFunc(_x2, row);
			}
		}

		~BicubicSplineInterp2D() {
			for (int i = 0; i < _m; i++)
				delete _rowSplines[i];
		}

		// Disable copy to avoid double-free
		BicubicSplineInterp2D(const BicubicSplineInterp2D&) = delete;
		BicubicSplineInterp2D& operator=(const BicubicSplineInterp2D&) = delete;

		/// @brief Evaluate the bicubic spline interpolation at (x1, x2).
		/// @param x1 First coordinate
		/// @param x2 Second coordinate
		/// @return Interpolated value z(x1, x2)

		Real operator()(Real x1, Real x2) const {
			// First interpolate along x2 for each x1[i] to get values at fixed x2
			Vector<Real> yv(_m);
			for (int i = 0; i < _m; i++)
				yv[i] = (*_rowSplines[i])(x2);

			// Then interpolate along x1
			SplineInterpRealFunc colSpline(_x1, yv);
			return colSpline(x1);
		}

		/// @brief Evaluate interpolation and partial derivatives at (x1, x2).
		/// @param x1 First coordinate
		/// @param x2 Second coordinate
		/// @param[out] z Interpolated value
		/// @param[out] dz_dx1 Partial derivative ∂z/∂x1
		/// @param[out] dz_dx2 Partial derivative ∂z/∂x2

		void interpWithDerivatives(Real x1, Real x2, Real& z, Real& dz_dx1, Real& dz_dx2) const {
			// Get values and x2-derivatives along rows
			Vector<Real> yv(_m), dyv_dx2(_m);
			for (int i = 0; i < _m; i++) {
				yv[i] = (*_rowSplines[i])(x2);
				dyv_dx2[i] = _rowSplines[i]->Derivative(x2);
			}

			// Interpolate values along x1
			SplineInterpRealFunc colSpline(_x1, yv);
			z = colSpline(x1);
			dz_dx1 = colSpline.Derivative(x1);

			// Interpolate x2-derivatives along x1
			SplineInterpRealFunc dSpline(_x1, dyv_dx2);
			dz_dx2 = dSpline(x1);
		}

		/// @name Grid Information
		/// @{

		/// /** @brief Get the first grid dimension. */

		int getGridDim1() const { return _m; }

		/// /** @brief Get the second grid dimension. */

		int getGridDim2() const { return _n; }

		/// @}
	};

} // namespace MML

#endif // MML_INTERPOLATION_2D_FUNCTION_H
