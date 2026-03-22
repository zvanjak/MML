///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixBandDiag.h                                                    ///
///  Description: Band diagonal matrix with compact storage                           ///
///               Efficient for tridiagonal and banded systems                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
/// @file MatrixBandDiag.h
/// @brief Band diagonal matrix with compact storage for sparse banded systems.
/// @ingroup Base
/// @section banddiag_overview Overview
/// BandDiagonalMatrix provides efficient storage and operations for matrices where
/// non-zero elements are confined to a diagonal band. Storage complexity is O(n·(m₁+m₂+1))
/// instead of O(n²), making it ideal for:
/// - Tridiagonal systems (m₁ = m₂ = 1)
/// - Pentadiagonal systems (m₁ = m₂ = 2)
/// - General banded systems from discretization schemes
/// @section banddiag_structure Matrix Structure
/// For an n×n matrix with lower bandwidth m₁ and upper bandwidth m₂:
/// @f[
/// A = \begin{pmatrix}
/// a_{0,0} & a_{0,1} & \cdots & a_{0,m_2} & 0 & \cdots \\
/// a_{1,0} & a_{1,1} & a_{1,2} & \cdots & a_{1,m_2+1} & 0 & \cdots \\
/// \vdots & & \ddots & & & \ddots & \\
/// a_{m_1,0} & \cdots & & & & & a_{m_1,m_1+m_2} & 0 \\
/// 0 & a_{m_1+1,1} & & & & & & \ddots \\
/// \vdots & & \ddots & & & & & \\
/// 0 & \cdots & 0 & a_{n-1,n-m_1-1} & \cdots & & & a_{n-1,n-1}
/// \end{pmatrix}
/// @f]
/// @section banddiag_storage Storage Layout
/// Elements are stored in a dense n × (m₁ + m₂ + 1) matrix where:
/// - Column m₁ contains the main diagonal
/// - Columns 0 to m₁-1 contain subdiagonals
/// - Columns m₁+1 to m₁+m₂ contain superdiagonals
/// Mapping: A[i][j] → _data[i][j - i + m₁]
/// @section banddiag_usage Usage Examples
/// @code{.cpp}
/// // Create a tridiagonal matrix: A[i][j] = data for |i-j| <= 1
/// int n = 5;
/// int m1 = 1, m2 = 1;  // Lower and upper bandwidth
/// Matrix<Real> data(n, 3);  // n rows, m1+m2+1=3 columns
/// // Fill: column 0 = subdiag, column 1 = main diag, column 2 = superdiag
/// for (int i = 0; i < n; i++) {
/// data[i][0] = -1.0;   // Subdiagonal
/// data[i][1] = 2.0;    // Main diagonal
/// data[i][2] = -1.0;   // Superdiagonal
/// }
/// BandDiagonalMatrix A(n, m1, m2, data);
/// // Element access
/// Real val = A(2, 3);      // Get element (returns 0 outside band)
/// A(1, 1) = 5.0;           // Set element (throws if outside band)
/// // Matrix-vector product (O(n) operations)
/// Vector<Real> x(n, 1.0);
/// Vector<Real> y = A * x;
/// @endcode
/// @see MatrixTriDiag Specialized tridiagonal matrix
/// @see LinAlgSolver Linear system solvers (may include banded solver)

#if !defined MML_MATRIX_BAND_DIAG_H
#define MML_MATRIX_BAND_DIAG_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Matrix/Matrix.h"

// Standard headers - include what we use
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace MML {
	/// @brief Band diagonal matrix with compact storage.
	/// @ingroup Base
	/// Stores only the non-zero diagonal band, reducing memory from O(n²) to O(n·(m₁+m₂+1)).
	/// Supports efficient O(n) matrix-vector multiplication.
	/// @par Storage Convention
	/// Data is stored in a Matrix<Real> of size n × (m₁ + m₂ + 1):
	/// - Diagonal at column m₁
	/// - Subdiagonals at columns 0 to m₁-1
	/// - Superdiagonals at columns m₁+1 to m₁+m₂

	class BandDiagonalMatrix {
		// The array a[0..n-1][0..m1+m2] stores A as follows: The diagonal elements are in a[0..n-1][m1].
		// Subdiagonal elements are in a[j..n-1][0..m1-1] with j > 0 appropriate to the number of
		// elements on each subdiagonal. Superdiagonal elements are in a[0..j][m1+1..m1+m2] with
		// j < n-1 appropriate to the number of elements on each superdiagonal.
		int _dim;			///< Matrix dimension (n×n)
		int _m1, _m2;		///< Lower and upper bandwidth
		Matrix<Real> _data; ///< Compact storage of band elements

	public:
		/// /** @name Dimension Accessors
		/// @{ */


		/// @brief Returns the number of rows.
		int rows() const { return _dim; }

		/// @brief Returns the number of columns.
		int cols() const { return _dim; }
		/// /** @} */


		/// @brief Constructs a band diagonal matrix.
		/// @param dim Matrix dimension n (creates n×n matrix).
		/// @param m1 Lower bandwidth (number of subdiagonals).
		/// @param m2 Upper bandwidth (number of superdiagonals).
		/// @param data Compact storage matrix of size dim × (m1+m2+1).
		/// @throws Runtime error if data dimensions don't match.

		BandDiagonalMatrix(int dim, int m1, int m2, const Matrix<Real>& data)
			: _dim(dim)
			, _m1(m1)
			, _m2(m2)
			, _data(data) {
			if (data.rows() != dim || data.cols() != m1 + m2 + 1)
				throw MatrixDimensionError("BandDiagonalMatrix constructor: data dimensions must be dim x (m1+m2+1)");
		}

		/// /** @name Element Access
		/// @{ */


		/// @brief Const element access - returns 0 for elements outside the band.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Element value (0 if outside band).

		Real operator()(int i, int j) const {
			if (i > j + _m1 || j > i + _m2)
				return 0.0;
			else
				return _data[i][j - i + _m1];
		}

		/// @brief Mutable element access - throws for elements outside the band.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element.
		/// @throws MatrixAccessBoundsError If (i,j) is outside the band.

		Real& operator()(int i, int j) {
			if (i > j + _m1 || j > i + _m2)
				throw MatrixAccessBoundsError("BandDiagonalMatrix::operator()", i, j, _dim, _dim);
			else
				return _data[i][j - i + _m1];
		}
		/// /** @} */


		/// /** @name Matrix-Vector Operations
		/// @{ */


		/// @brief Matrix-vector multiplication optimized for band structure.
		/// Complexity: O(n · (m₁ + m₂ + 1)) instead of O(n²).
		/// @param x Input vector of size n.
		/// @return Result vector y = A·x.
		/// @throws ArgumentError If vector size doesn't match.

		Vector<Real> operator*(const Vector<Real>& x) const {
			if (x.size() != _dim)
				throw ArgumentError("BandDiagonalMatrix::operator* - Vector size mismatch");

			Vector<Real> result(_dim);
			for (int i = 0; i < _dim; i++) {
				Real sum = 0.0;
				// Only iterate over the band (non-zero elements)
				int j_start = std::max(0, i - _m1);
				int j_end = std::min(_dim - 1, i + _m2);

				for (int j = j_start; j <= j_end; j++) {
					sum += (*this)(i, j) * x[j];
				}
				result[i] = sum;
			}
			return result;
		}
		/// /** @} */


		/// /** @name Conversion and Query
		/// @{ */


		/// @brief Converts to full matrix representation.
		/// @return Dense n×n matrix with all elements (including zeros).
		/// @note Useful for debugging or interfacing with dense algorithms.

		Matrix<Real> ToFullMatrix() const {
			Matrix<Real> full(_dim, _dim);
			for (int i = 0; i < _dim; i++) {
				for (int j = 0; j < _dim; j++) {
					full[i][j] = (*this)(i, j);
				}
			}
			return full;
		}

		/// @brief Returns the lower bandwidth (number of subdiagonals).
		int GetLowerBandwidth() const { return _m1; }

		/// @brief Returns the upper bandwidth (number of superdiagonals).
		int GetUpperBandwidth() const { return _m2; }

		/// @brief Returns the matrix dimension.
		int GetDimension() const { return _dim; }

		/// @brief Checks if an element position is within the band.
		/// @param i Row index.
		/// @param j Column index.
		/// @return True if element can be non-zero.

		bool IsInBand(int i, int j) const { return !(i > j + _m1 || j > i + _m2); }
		/// /** @} */


		/// /** @name Arithmetic Operations
		/// @{ */


		/// @brief Scalar multiplication.
		BandDiagonalMatrix operator*(Real scalar) const {
			Matrix<Real> new_data = _data * scalar;
			return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
		}

		/// @brief Scalar division.
		/// @throws DivisionByZeroError If scalar is near zero.

		BandDiagonalMatrix operator/(Real scalar) const {
			if (std::abs(scalar) < PrecisionValues<Real>::MatrixElementZeroThreshold)
				throw DivisionByZeroError("BandDiagonalMatrix::operator/ - Division by zero");
			Matrix<Real> new_data = _data / scalar;
			return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
		}

		/// @brief Unary negation.
		BandDiagonalMatrix operator-() const {
			Matrix<Real> new_data = -_data;
			return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
		}

		/// @brief Matrix addition (requires identical structure).
		/// @throws MatrixDimensionError If dimensions or bandwidth don't match.

		BandDiagonalMatrix operator+(const BandDiagonalMatrix& other) const {
			if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
				throw MatrixDimensionError("BandDiagonalMatrix::operator+ - Incompatible dimensions or bandwidth");

			Matrix<Real> new_data = _data + other._data;
			return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
		}

		/// @brief Matrix subtraction (requires identical structure).
		/// @throws MatrixDimensionError If dimensions or bandwidth don't match.

		BandDiagonalMatrix operator-(const BandDiagonalMatrix& other) const {
			if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
				throw MatrixDimensionError("BandDiagonalMatrix::operator- - Incompatible dimensions or bandwidth");

			Matrix<Real> new_data = _data - other._data;
			return BandDiagonalMatrix(_dim, _m1, _m2, new_data);
		}

		/// @brief In-place scalar multiplication.
		BandDiagonalMatrix& operator*=(Real scalar) {
			_data *= scalar;
			return *this;
		}

		/// @brief In-place scalar division.
		/// @throws DivisionByZeroError If scalar is near zero.

		BandDiagonalMatrix& operator/=(Real scalar) {
			if (std::abs(scalar) < PrecisionValues<Real>::MatrixElementZeroThreshold)
				throw DivisionByZeroError("BandDiagonalMatrix::operator/= - Division by zero");
			_data /= scalar;
			return *this;
		}

		/// @brief In-place addition.
		/// @throws MatrixDimensionError If dimensions or bandwidth don't match.

		BandDiagonalMatrix& operator+=(const BandDiagonalMatrix& other) {
			if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
				throw MatrixDimensionError("BandDiagonalMatrix::operator+= - Incompatible dimensions or bandwidth");

			_data += other._data;
			return *this;
		}

		/// @brief In-place subtraction.
		/// @throws MatrixDimensionError If dimensions or bandwidth don't match.

		BandDiagonalMatrix& operator-=(const BandDiagonalMatrix& other) {
			if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
				throw MatrixDimensionError("BandDiagonalMatrix::operator-= - Incompatible dimensions or bandwidth");

			_data -= other._data;
			return *this;
		}
		/// /** @} */


		/// /** @name Matrix Properties
		/// @{ */


		/// @brief Returns the trace (sum of diagonal elements).
		Real Trace() const {
			Real trace = 0.0;
			for (int i = 0; i < _dim; i++)
				trace += (*this)(i, i);
			return trace;
		}

		/// @brief Computes the Frobenius norm.
		/// @return √(Σ|aᵢⱼ|²) over all band elements.

		Real NormFrobenius() const {
			Real sum = 0.0;
			for (int i = 0; i < _dim; i++) {
				int j_start = std::max(0, i - _m1);
				int j_end = std::min(_dim - 1, i + _m2);
				for (int j = j_start; j <= j_end; j++) {
					Real val = (*this)(i, j);
					sum += val * val;
				}
			}
			return std::sqrt(sum);
		}

		/// @brief Computes the infinity norm (max row sum).
		/// @return max_i Σⱼ|aᵢⱼ|

		Real NormInf() const {
			Real max_val = 0.0;
			for (int i = 0; i < _dim; i++) {
				Real row_sum = 0.0;
				int j_start = std::max(0, i - _m1);
				int j_end = std::min(_dim - 1, i + _m2);
				for (int j = j_start; j <= j_end; j++)
					row_sum += std::abs((*this)(i, j));
				max_val = std::max(max_val, row_sum);
			}
			return max_val;
		}

		/// @brief Computes the transpose matrix.
		/// @return Transposed matrix (swaps m₁ and m₂).

		BandDiagonalMatrix Transpose() const {
			// For transpose, m1 and m2 are swapped
			Matrix<Real> trans_data(_dim, _m1 + _m2 + 1);

			for (int i = 0; i < _dim; i++) {
				int j_start = std::max(0, i - _m1);
				int j_end = std::min(_dim - 1, i + _m2);
				for (int j = j_start; j <= j_end; j++) {
					// A^T[j][i] = A[i][j]
					// In band storage: trans_data[j][i - j + m2] = data[i][j - i + m1]
					trans_data[j][i - j + _m2] = _data[i][j - i + _m1];
				}
			}

			return BandDiagonalMatrix(_dim, _m2, _m1, trans_data);
		}

		/// @brief Checks if the matrix is symmetric within tolerance.
		/// @param tolerance Maximum allowed difference between A[i][j] and A[j][i].
		/// @return True if m₁ == m₂ and |A[i][j] - A[j][i]| ≤ tolerance for all i,j.

		bool IsSymmetric(Real tolerance = 1e-10) const {
			if (_m1 != _m2)
				return false; // Symmetric requires equal bandwidths

			for (int i = 0; i < _dim; i++) {
				int j_start = std::max(0, i - _m1);
				int j_end = std::min(_dim - 1, i + _m2);
				for (int j = j_start; j <= j_end; j++) {
					if (std::abs((*this)(i, j) - (*this)(j, i)) > tolerance)
						return false;
				}
			}
			return true;
		}

		/// @brief Checks equality with another band matrix.
		/// @param other Matrix to compare with.
		/// @param tolerance Maximum allowed element difference.
		/// @return True if dimensions, bandwidth, and elements match.

		bool IsEqualTo(const BandDiagonalMatrix& other, Real tolerance = Defaults::MatrixIsEqualTolerance) const {
			if (_dim != other._dim || _m1 != other._m1 || _m2 != other._m2)
				return false;

			return _data.IsEqualTo(other._data, tolerance);
		}
		/// /** @} */


		/// /** @name Diagonal Extraction
		/// @{ */


		/// @brief Extracts the main diagonal as a vector.
		Vector<Real> GetDiagonal() const {
			Vector<Real> diag(_dim);
			for (int i = 0; i < _dim; i++)
				diag[i] = (*this)(i, i);
			return diag;
		}

		/// @brief Extracts a specific diagonal.
		/// @param k Diagonal index: 0 = main, k > 0 = below main, k < 0 = above main.
		/// @return Vector of diagonal elements.
		/// @throws IndexError If diagonal is outside the band.

		Vector<Real> GetDiagonal(int k) const {
			if (k > _m1 || k < -_m2)
				throw IndexError("BandDiagonalMatrix::GetDiagonal - Diagonal outside band");

			int diag_length;
			int start_row, start_col;

			if (k >= 0) {
				// Below or on main diagonal
				diag_length = _dim - k;
				start_row = k;
				start_col = 0;
			} else {
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
		/// /** @} */


		/// /** @name Output
		/// @{ */


		/// @brief Prints the full matrix representation to a stream.
		/// @param stream Output stream.
		/// @param width Field width for each element.
		/// @param precision Decimal precision for floating-point output.

		void Print(std::ostream& stream, int width, int precision) const {
			stream << "Dim: " << _dim << std::endl;
			for (int i = 0; i < _dim; i++) {
				stream << "[ ";
				for (int j = 0; j < _dim; j++) {
					stream << std::setw(width) << std::setprecision(precision) << (*this)(i, j) << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
		/// /** @} */
	};

	/// @brief Scalar-matrix multiplication (commutative).
	inline BandDiagonalMatrix operator*(Real scalar, const BandDiagonalMatrix& mat) { return mat * scalar; }
} // namespace MML

#endif
