///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixBandDiag.h                                                    ///
///  Description: Band diagonal matrix with compact storage                           ///
///               Efficient for tridiagonal and banded systems                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_BAND_DIAG_H
#define MML_MATRIX_BAND_DIAG_H

#include "MMLBase.h"

#include "base/Matrix.h"

namespace MML
{
	class BandDiagonalMatrix
	{
		// The array a[0..n-1][0..m1+m2] stores A as follows: The diagonal elements are in a[0..n-1][m1].
		// Subdiagonal elements are in a[j..n-1][0..m1-1] with j > 0 appropriate to the number of
		// elements on each subdiagonal. Superdiagonal elements are in a[0..j][m1+1..m1+m2] with
		// j < n-1 appropriate to the number of elements on each superdiagonal.        
		int _dim;
		int _m1, _m2;
		Matrix<Real> _data;

	public:
		int RowNum() const { return _dim; }
		int ColNum() const { return _dim; }

		BandDiagonalMatrix(int dim, int m1, int m2, const Matrix<Real>& data) : _dim(dim), _m1(m1), _m2(m2), _data(data)
		{
			if (data.RowNum() != dim || data.ColNum() != m1 + m2 + 1)
				throw("Error in BandDiagonalMatrix constructor: wrong dimensions");
		}

		Real  operator()(int i, int j) const {
			if (i > j + _m1 || j > i + _m2)
				return 0.0;
			else
				return _data[i][j - i + _m1];
		}
		Real& operator()(int i, int j) {
			if (i > j + _m1 || j > i + _m2)
				throw MatrixAccessBoundsError("BandDiagonalMatrix::operator()", i, j, _dim, _dim);
			else
				return _data[i][j - i + _m1];
		}


	// Matrix-vector multiplication for band diagonal matrix
	// Optimized to only multiply non-zero elements
	Vector<Real> operator*(const Vector<Real>& x) const
	{
		if (x.size() != _dim)
			throw std::invalid_argument("BandDiagonalMatrix::operator* - Vector size mismatch");

		Vector<Real> result(_dim);
		for (int i = 0; i < _dim; i++)
		{
			Real sum = 0.0;
			// Only iterate over the band (non-zero elements)
			int j_start = std::max(0, i - _m1);
			int j_end = std::min(_dim - 1, i + _m2);
			
			for (int j = j_start; j <= j_end; j++)
			{
				sum += (*this)(i, j) * x[j];
			}
			result[i] = sum;
		}
		return result;
	}

	// Get full matrix representation (for debugging/testing)
	Matrix<Real> ToFullMatrix() const
	{
		Matrix<Real> full(_dim, _dim);
		for (int i = 0; i < _dim; i++)
		{
			for (int j = 0; j < _dim; j++)
			{
				full[i][j] = (*this)(i, j);
			}
		}
		return full;
	}

	// Get bandwidth information
	int GetLowerBandwidth() const { return _m1; }
	int GetUpperBandwidth() const { return _m2; }
	int GetDimension() const { return _dim; }

	// Check if element is within the band
	bool IsInBand(int i, int j) const
	{
		return !(i > j + _m1 || j > i + _m2);
	}

	///////////////////     ARITHMETIC OPERATIONS     ///////////////////

	// Scalar multiplication
	BandDiagonalMatrix operator*(Real scalar) const
	{
		Matrix<Real> new_data = _data * scalar;
		return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
	}

	// Scalar division
	BandDiagonalMatrix operator/(Real scalar) const
	{
		if (std::abs(scalar) < PrecisionValues<Real>::MatrixElementZeroThreshold)
			throw std::invalid_argument("BandDiagonalMatrix::operator/ - Division by zero");
		Matrix<Real> new_data = _data / scalar;
		return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
	}

	// Unary negation
	BandDiagonalMatrix operator-() const
	{
		Matrix<Real> new_data = -_data;
		return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
	}

	// Matrix addition (requires compatible bandwidth)
	BandDiagonalMatrix operator+(const BandDiagonalMatrix& other) const
	{
		if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
			throw std::invalid_argument("BandDiagonalMatrix::operator+ - Incompatible dimensions or bandwidth");
		
		Matrix<Real> new_data = _data + other._data;
		return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
	}

	// Matrix subtraction (requires compatible bandwidth)
	BandDiagonalMatrix operator-(const BandDiagonalMatrix& other) const
	{
		if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
			throw std::invalid_argument("BandDiagonalMatrix::operator- - Incompatible dimensions or bandwidth");
		
		Matrix<Real> new_data = _data - other._data;
		return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
	}

	// In-place scalar multiplication
	BandDiagonalMatrix& operator*=(Real scalar)
	{
		_data *= scalar;
		return *this;
	}

	// In-place scalar division
	BandDiagonalMatrix& operator/=(Real scalar)
	{
		if (std::abs(scalar) < PrecisionValues<Real>::MatrixElementZeroThreshold)
			throw std::invalid_argument("BandDiagonalMatrix::operator/= - Division by zero");
		_data /= scalar;
		return *this;
	}

	// In-place addition
	BandDiagonalMatrix& operator+=(const BandDiagonalMatrix& other)
	{
		if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
			throw std::invalid_argument("BandDiagonalMatrix::operator+= - Incompatible dimensions or bandwidth");
		
		_data += other._data;
		return *this;
	}

	// In-place subtraction
	BandDiagonalMatrix& operator-=(const BandDiagonalMatrix& other)
	{
		if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
			throw std::invalid_argument("BandDiagonalMatrix::operator-= - Incompatible dimensions or bandwidth");
		
		_data -= other._data;
		return *this;
	}

	///////////////////     MATRIX PROPERTIES     ///////////////////

	// Trace (sum of diagonal elements)
	Real Trace() const
	{
		Real trace = 0.0;
		for (int i = 0; i < _dim; i++)
			trace += (*this)(i, i);
		return trace;
	}

	// Frobenius norm (sqrt of sum of squares of all elements)
	Real NormFrobenius() const
	{
		Real sum = 0.0;
		for (int i = 0; i < _dim; i++)
		{
			int j_start = std::max(0, i - _m1);
			int j_end = std::min(_dim - 1, i + _m2);
			for (int j = j_start; j <= j_end; j++)
			{
				Real val = (*this)(i, j);
				sum += val * val;
			}
		}
		return std::sqrt(sum);
	}

	// Maximum absolute value of elements
	Real NormInf() const
	{
		Real max_val = 0.0;
		for (int i = 0; i < _dim; i++)
		{
			Real row_sum = 0.0;
			int j_start = std::max(0, i - _m1);
			int j_end = std::min(_dim - 1, i + _m2);
			for (int j = j_start; j <= j_end; j++)
				row_sum += std::abs((*this)(i, j));
			max_val = std::max(max_val, row_sum);
		}
		return max_val;
	}

	// Transpose (swaps upper and lower bandwidth)
	BandDiagonalMatrix Transpose() const
	{
		// For transpose, m1 and m2 are swapped
		Matrix<Real> trans_data(_dim, _m1 + _m2 + 1);
		
		for (int i = 0; i < _dim; i++)
		{
			int j_start = std::max(0, i - _m1);
			int j_end = std::min(_dim - 1, i + _m2);
			for (int j = j_start; j <= j_end; j++)
			{
				// A^T[j][i] = A[i][j]
				// In band storage: trans_data[j][i - j + m2] = data[i][j - i + m1]
				trans_data[j][i - j + _m2] = _data[i][j - i + _m1];
			}
		}
		
		return BandDiagonalMatrix(_dim, _m2, _m1, trans_data);
	}

	// Check if matrix is symmetric
	bool IsSymmetric(Real tolerance = 1e-10) const
	{
		if (_m1 != _m2)
			return false;  // Symmetric requires equal bandwidths
		
		for (int i = 0; i < _dim; i++)
		{
			int j_start = std::max(0, i - _m1);
			int j_end = std::min(_dim - 1, i + _m2);
			for (int j = j_start; j <= j_end; j++)
			{
				if (std::abs((*this)(i, j) - (*this)(j, i)) > tolerance)
					return false;
			}
		}
		return true;
	}

	// Equality comparison
	bool IsEqual(const BandDiagonalMatrix& other, Real tolerance = 1e-10) const
	{
		if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
			return false;
		
		return _data.IsEqualTo(other._data, tolerance);
	}

	// Get diagonal as vector
	Vector<Real> GetDiagonal() const
	{
		Vector<Real> diag(_dim);
		for (int i = 0; i < _dim; i++)
			diag[i] = (*this)(i, i);
		return diag;
	}

	// Get subdiagonal k (k=0 is main diagonal, k>0 is below, k<0 is above)
	Vector<Real> GetDiagonal(int k) const
	{
		if (k > _m1 || k < -_m2)
			throw std::invalid_argument("BandDiagonalMatrix::GetDiagonal - Diagonal outside band");
		
		int diag_length;
		int start_row, start_col;
		
		if (k >= 0)
		{
			// Below or on main diagonal
			diag_length = _dim - k;
			start_row = k;
			start_col = 0;
		}
		else
		{
			// Above main diagonal
			diag_length = _dim + k;
			start_row = 0;
			start_col = -k;
		}
		
		Vector<Real> diag(diag_length);
		for (int i = 0; i < diag_length; i++)
			diag[i] = (*this)(start_row + i, start_col + i);
		
		return diag;
	}

		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Dim: " << _dim << std::endl;
			for (int i = 0; i < _dim; i++)
			{
				stream << "[ ";
				for (int j = 0; j < _dim; j++) {
					stream << std::setw(width) << std::setprecision(precision) << (*this)(i, j) << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
	};

	// Friend function for scalar * matrix
	inline BandDiagonalMatrix operator*(Real scalar, const BandDiagonalMatrix& mat)
	{
		return mat * scalar;
	}
}

#endif