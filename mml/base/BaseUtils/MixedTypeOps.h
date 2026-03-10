///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MixedTypeOps.h                                                      ///
///  Description: Mixed Complex/Real vector and matrix operations                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_MIXED_TYPE_OPS_H
#define MML_MIXED_TYPE_OPS_H

#include <complex>

#include "mml/MMLBase.h"
#include "mml/MMLExceptions.h"
#include "mml/base/Vector/Vector.h"
#include "mml/base/Matrix/Matrix.h"

namespace MML
{
	namespace Utils
	{
		// ============================================================================
		// Vector<Complex> - Vector<Real> Operations
		// ============================================================================

		/// @brief Adds complex and real vectors element-wise.
		/// @param a Complex vector
		/// @param b Real vector
		/// @return Complex vector with a[i] + b[i]
		/// @throws VectorDimensionError if vectors have different sizes
		static Vector<Complex> AddVec(const Vector<Complex>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("AddVec(Complex, Real) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] + b[i];
			return ret;
		}

		/// @brief Adds real and complex vectors element-wise.
		/// @param a Real vector
		/// @param b Complex vector
		/// @return Complex vector with a[i] + b[i]
		/// @throws VectorDimensionError if vectors have different sizes
		static Vector<Complex> AddVec(const Vector<Real>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("AddVec(Real, Complex) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] + b[i];
			return ret;
		}

		/// @brief Subtracts real vector from complex vector element-wise.
		/// @param a Complex vector
		/// @param b Real vector
		/// @return Complex vector with a[i] - b[i]
		/// @throws VectorDimensionError if vectors have different sizes
		static Vector<Complex> SubVec(const Vector<Complex>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("SubVec(Complex, Real) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] - b[i];
			return ret;
		}

		/// @brief Subtracts complex vector from real vector element-wise.
		/// @param a Real vector
		/// @param b Complex vector
		/// @return Complex vector with a[i] - b[i]
		/// @throws VectorDimensionError if vectors have different sizes
		static Vector<Complex> SubVec(const Vector<Real>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("SubVec(Real, Complex) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] - b[i];
			return ret;
		}

		/// @brief Multiplies complex scalar by real vector.
		/// @param a Complex scalar
		/// @param b Real vector
		/// @return Complex vector with a * b[i]
		static Vector<Complex> MulVec(const Complex& a, const Vector<Real>& b)
		{
			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a * b[i];
			return ret;
		}

		/// @brief Multiplies real vector by complex scalar.
		/// @param a Real vector
		/// @param b Complex scalar
		/// @return Complex vector with b * a[i]
		static Vector<Complex> MulVec(const Vector<Real>& a, const Complex& b)
		{
			Vector<Complex> ret(a.size());;
			for (int i = 0; i < a.size(); i++)
				ret[i] = b * a[i];
			return ret;
		}

		// ============================================================================
		// Matrix<Complex> - Matrix<Real> Operations
		// ============================================================================

		/// @brief Adds complex and real matrices element-wise.
		/// @param a Complex matrix
		/// @param b Real matrix
		/// @return Complex matrix with a[i][j] + b[i][j]
		/// @throws MatrixDimensionError if matrices have different dimensions
		static Matrix<Complex> AddMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.rows() != b.rows() || a.cols() != b.cols())
				throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.rows(), a.cols(), b.rows(), b.cols());

			Matrix<Complex> ret(a);
			for (int i = 0; i < b.rows(); i++)
				for (int j = 0; j < b.cols(); j++)
					ret[i][j] += b[i][j];
			return ret;
		}

		/// @brief Adds real and complex matrices element-wise.
		/// @param a Real matrix
		/// @param b Complex matrix
		/// @return Complex matrix with a[i][j] + b[i][j]
		/// @throws MatrixDimensionError if matrices have different dimensions
		static Matrix<Complex> AddMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.rows() != b.rows() || a.cols() != b.cols())
				throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.rows(), a.cols(), b.rows(), b.cols());

			Matrix<Complex> ret(b);
			for (int i = 0; i < a.rows(); i++)
				for (int j = 0; j < a.cols(); j++)
					ret[i][j] += a[i][j];
			return ret;
		}

		/// @brief Subtracts real matrix from complex matrix element-wise.
		/// @param a Complex matrix
		/// @param b Real matrix
		/// @return Complex matrix with a[i][j] - b[i][j]
		/// @throws MatrixDimensionError if matrices have different dimensions
		static Matrix<Complex> SubMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.rows() != b.rows() || a.cols() != b.cols())
				throw MatrixDimensionError("SubMat(Complex, Real) - must be same dim", a.rows(), a.cols(), b.rows(), b.cols());

			Matrix<Complex> ret(a);
			for (int i = 0; i < b.rows(); i++)
				for (int j = 0; j < b.cols(); j++)
					ret[i][j] -= b[i][j];
			return ret;
		}

		/// @brief Subtracts complex matrix from real matrix element-wise.
		/// @param a Real matrix
		/// @param b Complex matrix
		/// @return Complex matrix with a[i][j] - b[i][j]
		/// @throws MatrixDimensionError if matrices have different dimensions
		static Matrix<Complex> SubMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.rows() != b.rows() || a.cols() != b.cols())
				throw MatrixDimensionError("SubMat(Real, Complex) - must be same dim", a.rows(), a.cols(), b.rows(), b.cols());

			Matrix<Complex> ret(b);
			for (int i = 0; i < a.rows(); i++)
				for (int j = 0; j < a.cols(); j++)
					ret[i][j] = a[i][j] - b[i][j];
			return ret;
		}

		/// @brief Multiplies complex scalar by real matrix.
		/// @param a Complex scalar
		/// @param b Real matrix
		/// @return Complex matrix with a * b[i][j]
		static Matrix<Complex> MulMat(const Complex& a, const Matrix<Real>& b)
		{
			Matrix<Complex> ret(b.rows(), b.cols());

			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret[i][j] = a * b[i][j];
				}

			return ret;
		}

		/// @brief Multiplies real matrix by complex scalar.
		/// @param a Real matrix
		/// @param b Complex scalar
		/// @return Complex matrix with b * a[i][j]
		static Matrix<Complex> MulMat(const Matrix<Real>& a, const Complex& b)
		{
			Matrix<Complex> ret(a.rows(), a.cols());

			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret[i][j] = b * a[i][j];
				}

			return ret;
		}
		
		/// @brief Multiplies complex matrix by real matrix.
		/// @param a Complex matrix (m×n)
		/// @param b Real matrix (n×p)
		/// @return Complex matrix (m×p)
		/// @throws MatrixDimensionError if inner dimensions don't match
		static Matrix<Complex> MulMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.cols() != b.rows())
				throw MatrixDimensionError("MulMat(Complex, Real) - a.colNum must be equal to b.rowNum", a.rows(), a.cols(), b.rows(), b.cols());

			Matrix<Complex> ret(a.rows(), b.cols());
			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret[i][j] = 0.0;
					for (int k = 0; k < a.cols(); k++)
						ret[i][j] += a[i][k] * b[k][j];
				}

			return ret;
		}

		/// @brief Multiplies real matrix by complex matrix.
		/// @param a Real matrix (m×n)
		/// @param b Complex matrix (n×p)
		/// @return Complex matrix (m×p)
		/// @throws MatrixDimensionError if inner dimensions don't match
		static Matrix<Complex> MulMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.cols() != b.rows())
				throw MatrixDimensionError("MulMat(Real, Complex) - a.colNum must be equal to b.rowNum", a.rows(), a.cols(), b.rows(), b.cols());

			Matrix<Complex> ret(a.rows(), b.cols());
			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret[i][j] = 0.0;
					for (int k = 0; k < a.cols(); k++)
						ret[i][j] += a[i][k] * b[k][j];
				}

			return ret;
		}

		// ============================================================================
		// Matrix-Vector Mixed Operations
		// ============================================================================

		/// @brief Multiplies real matrix by complex vector.
		/// @param a Real matrix (m×n)
		/// @param b Complex vector (size n)
		/// @return Complex vector (size m)
		/// @throws MatrixDimensionError if matrix columns don't match vector size
		static Vector<Complex> MulMatVec(const Matrix<Real>& a, const Vector<Complex>& b)
		{
			if (a.cols() != b.size())
				throw MatrixDimensionError("MulMatVec(Real, Complex) - a.colNum must be equal to vector size", a.rows(), a.cols(), (int)b.size(), -1);

			Vector<Complex> ret(a.rows());
			for (int i = 0; i < a.rows(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.cols(); j++)
					ret[i] += a[i][j] * b[j];
			}
			return ret;
		}

		/// @brief Multiplies complex matrix by real vector.
		/// @param a Complex matrix (m×n)
		/// @param b Real vector (size n)
		/// @return Complex vector (size m)
		/// @throws MatrixDimensionError if matrix columns don't match vector size
		static Vector<Complex> MulMatVec(const Matrix<Complex>& a, const Vector<Real>& b)
		{
			if (a.cols() != b.size())
				throw MatrixDimensionError("MulMatVec(Complex, Real) - a.colNum must be equal to vector size", a.rows(), a.cols(), (int)b.size(), -1);

			Vector<Complex> ret(a.rows());
			for (int i = 0; i < a.rows(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.cols(); j++)
					ret[i] += a[i][j] * b[j];
			}
			return ret;
		}

		/// @brief Multiplies complex row vector by real matrix.
		/// @param a Complex vector (size m)
		/// @param b Real matrix (m×n)
		/// @return Complex vector (size n)
		/// @throws MatrixDimensionError if vector size doesn't match matrix rows
		static Vector<Complex> MulVecMat(const Vector<Complex>& a, const Matrix<Real>& b)
		{
			if (a.size() != b.rows())
				throw MatrixDimensionError("MulVecMat(Complex, Real) - vector size must be equal to b.rowNum", (int)a.size(), -1, b.rows(), b.cols());

			Vector<Complex> ret(b.cols());
			for (int i = 0; i < b.cols(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.rows(); j++)
					ret[i] += a[j] * b[j][i];
			}
			return ret;
		}

		/// @brief Multiplies real row vector by complex matrix.
		/// @param a Real vector (size m)
		/// @param b Complex matrix (m×n)
		/// @return Complex vector (size n)
		/// @throws MatrixDimensionError if vector size doesn't match matrix rows
		static Vector<Complex> MulVecMat(const Vector<Real>& a, const Matrix<Complex>& b)
		{
			if (a.size() != b.rows())
				throw MatrixDimensionError("MulVecMat(Real, Complex) - vector size must be equal to b.rowNum", (int)a.size(), -1, b.rows(), b.cols());

			Vector<Complex> ret(b.cols());
			for (int i = 0; i < b.cols(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.rows(); j++)
					ret[i] += a[j] * b[j][i];
			}
			return ret;
		}

	} // namespace Utils
} // namespace MML

#endif // MML_MIXED_TYPE_OPS_H
