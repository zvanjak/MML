/// @file MatrixSym.h
/// @brief Symmetric matrix class with optimized storage using only triangular portion.
/// @details This file provides a symmetric matrix class that exploits the symmetry
/// property (A[i][j] = A[j][i]) to store only the lower triangular elements, reducing
/// memory usage by approximately 50% compared to full matrix storage.
/// **Key Features:**
/// - **Compact storage**: Only n(n+1)/2 elements instead of n²
/// - **Memory efficient**: ~50% reduction for large matrices
/// - **Automatic symmetry**: Setting (i,j) automatically sets (j,i)
/// - **Full matrix interface**: Access via (i,j) for any valid indices
/// - **Conversion utilities**: Convert to/from regular Matrix class
/// **Storage Layout:**
/// Elements are stored row-wise in a 1D array (lower triangular):
/// @code
/// [a₀₀, a₁₀, a₁₁, a₂₀, a₂₁, a₂₂, a₃₀, a₃₁, a₃₂, a₃₃, ...]
/// @endcode
/// For element (i, j) where i ≥ j: index = i(i+1)/2 + j
/// **Memory Comparison (n×n matrices):**
/// | n | Full Matrix | MatrixSym | Savings |
/// |---|-------------|-----------|---------|
/// | 100 | 10,000 | 5,050 | 49.5% |
/// | 1000 | 1,000,000 | 500,500 | 50.0% |
/// | 5000 | 25,000,000 | 12,502,500 | 50.0% |
/// **When to Use:**
/// - Covariance matrices
/// - Hessian matrices (optimization)
/// - Graph adjacency matrices (undirected)
/// - Mass/stiffness matrices (FEM)
/// - Quadratic form matrices
/// **Important Note:**
/// Matrix-matrix multiplication of two symmetric matrices does NOT generally
/// produce a symmetric result, so `operator*` returns a regular `Matrix`.
/// **Usage Example:**
/// @code
/// using namespace MML;
/// // Create 3×3 symmetric matrix from lower triangular elements
/// MatrixSym<double> S(3, {
/// 1,        // a₀₀
/// 2, 3,     // a₁₀, a₁₁
/// 4, 5, 6   // a₂₀, a₂₁, a₂₂
/// });
/// // Full matrix view:
/// // [1, 2, 4]
/// // [2, 3, 5]
/// // [4, 5, 6]
/// double trace = S.Trace();  // = 1 + 3 + 6 = 10
/// double norm = S.NormFrobenius();
/// // Create from regular matrix (symmetrizes: (A + A^T)/2)
/// Matrix<double> M(3, 3);
/// // ... fill M ...
/// auto S2 = MatrixSym<double>::FromFullMatrix(M);
/// @endcode
/// @see Matrix For general non-symmetric matrices
/// @see MatrixBandDiag For banded matrices
/// @see MatrixTriDiag For tridiagonal matrices
/// @author Zvonimir Vanjak
/// @date 2024-2025

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixSym.h                                                         ///
///  Description: Symmetric matrix class with optimized storage and operations        ///
///               Stores only upper/lower triangular portion                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_SYM_H
#define MML_MATRIX_SYM_H

#include "MMLBase.h"
#include "base/MatrixPrintFormat.h"

#include "Matrix.h"

#include <vector>
#include <algorithm>
#include <stdexcept>

// Standard headers - include what we use
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>

namespace MML 
{
	/// @brief Allocation limits for symmetric matrices.
	/// @details Prevents accidental allocation of excessively large matrices.

	struct MatrixSymLimits {
		static constexpr int MAX_DIMENSION = 10000;				 ///< Maximum matrix dimension
		static constexpr size_t MAX_ELEMENTS = 50'000'000; ///< Maximum storage elements
	};

	/// @brief Symmetric matrix with optimized compact storage.
	/// @details MatrixSym stores only the lower triangular elements of a symmetric
	/// matrix, automatically providing symmetric access through operator(). This
	/// reduces storage to n(n+1)/2 elements for an n×n matrix.
	/// **Indexing:** For elements (i,j), the class automatically handles the
	/// symmetry: accessing (i,j) or (j,i) returns the same stored value.
	/// @tparam Type Element type (typically Real, float, or Complex)

	template<class Type>
	class MatrixSym 
	{
	private:
		int _dim;								 ///< Matrix dimension (rows = cols = _dim)
		std::vector<Type> _data; ///< Contiguous storage for lower triangular elements

		/// @brief Computes number of stored elements for given dimension.
		/// @param dim Matrix dimension.
		/// @return n(n+1)/2 elements for n×n matrix.

		static constexpr size_t numElements(int dim) noexcept { return static_cast<size_t>(dim) * (dim + 1) / 2; }

		/// @brief Computes linear index into storage array.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Index into _data array (automatically handles i < j by swapping).
		/// @details Storage layout: row-wise lower triangular
		/// [a₀₀, a₁₀, a₁₁, a₂₀, a₂₁, a₂₂, ...]

		static constexpr size_t linearIndex(int i, int j) noexcept {
			// Ensure i >= j for lower triangular access
			if (i < j)
				std::swap(i, j);
			return static_cast<size_t>(i) * (i + 1) / 2 + j;
		}

		/// @brief Validates dimension against limits.
		void validateDimension(int dim) const {
			if (dim < 0)
				throw MatrixDimensionError("MatrixSym: dimension cannot be negative", dim, -1, -1, -1);
			if (dim > MatrixSymLimits::MAX_DIMENSION)
				throw MatrixDimensionError("MatrixSym: dimension exceeds maximum allowed", dim, MatrixSymLimits::MAX_DIMENSION, -1, -1);

			size_t requiredElements = numElements(dim);
			if (requiredElements > MatrixSymLimits::MAX_ELEMENTS)
				throw MatrixAllocationError("MatrixSym: allocation would exceed maximum elements", dim, -1);
		}

	public:
		typedef Type value_type; ///< Element type alias for STL compatibility

		/// @brief Default constructor - creates empty (0×0) matrix.
		MatrixSym() noexcept
				: _dim(0)
				, _data() {}

		/// @brief Dimension constructor - creates dim×dim matrix initialized to zero.
		/// @param dim Matrix dimension.
		/// @throws MatrixDimensionError If dim < 0 or exceeds limits.
		explicit MatrixSym(int dim)
				: _dim(0)
				, _data() {
			if (dim == 0)
				return;

			validateDimension(dim);
			_dim = dim;
			_data.resize(numElements(dim), Type{0});
		}

		/// @brief Dimension + value constructor - creates matrix initialized to val.
		/// @param dim Matrix dimension.
		/// @param val Initial value for all stored elements.
		MatrixSym(int dim, Type val)
				: _dim(0)
				, _data() {
			if (dim == 0)
				return;

			validateDimension(dim);
			_dim = dim;
			_data.resize(numElements(dim), val);
		}

		/// @brief Initializer list constructor - takes lower triangular elements row-wise.
		/// @param dim Matrix dimension.
		/// @param values Lower triangular elements: {a₀₀, a₁₀, a₁₁, a₂₀, a₂₁, a₂₂, ...}
		/// @throws MatrixDimensionError If values.size() != dim(dim+1)/2.
		/// @code
		/// MatrixSym<double> S(3, {1, 2, 3, 4, 5, 6});
		/// // Creates symmetric matrix:
		/// // [1, 2, 4]
		/// // [2, 3, 5]
		/// // [4, 5, 6]
		/// @endcode

		MatrixSym(int dim, std::initializer_list<Type> values)
				: _dim(0)
				, _data() {
			if (dim == 0) {
				if (values.size() != 0)
					throw MatrixDimensionError("MatrixSym: non-empty initializer for zero dimension", 0, -1, -1, -1);
				return;
			}

			validateDimension(dim);

			size_t expected = numElements(dim);
			if (values.size() != expected)
				throw MatrixDimensionError("MatrixSym: initializer list size mismatch - expected " + std::to_string(expected) + " but got " +
																		 std::to_string(values.size()),
																	 dim, -1, static_cast<int>(values.size()), -1);

			_dim = dim;
			_data.assign(values.begin(), values.end());
		}

		/// @brief Copy constructor.
		MatrixSym(const MatrixSym& other)
				: _dim(other._dim)
				, _data(other._data) {}

		/// @brief Move constructor.
		MatrixSym(MatrixSym&& other) noexcept
				: _dim(other._dim)
				, _data(std::move(other._data)) {
			other._dim = 0;
		}

		/// @brief Destructor - default is fine with std::vector.
		~MatrixSym() = default;
		/// /** @} */


		/// /** @name Dimension Accessors
		/// @{ */


		/// @brief Returns matrix dimension (rows = cols = Dim()).
		int Dim() const noexcept { return _dim; }

		/// @brief Returns number of rows (preferred API).
		int rows() const noexcept { return _dim; }

		/// @brief Returns number of columns (preferred API).
		int cols() const noexcept { return _dim; }

		/// @brief Returns number of stored elements: n(n+1)/2.
		size_t size() const noexcept { return _data.size(); }

		/// @brief Checks if matrix is empty (0×0).
		bool empty() const noexcept { return _dim == 0; }

		/// @brief Raw data pointer for interoperability with libraries.
		Type* data() noexcept { return _data.data(); }

		/// @brief Raw data pointer (const) for interoperability.
		const Type* data() const noexcept { return _data.data(); }


		/// @brief Resizes the matrix, destroying existing data.
		/// @param newDim New matrix dimension.
		void Resize(int newDim) {
			if (newDim == _dim)
				return;

			if (newDim == 0) {
				_dim = 0;
				_data.clear();
				return;
			}

			validateDimension(newDim);
			_dim = newDim;
			_data.assign(numElements(newDim), Type{0});
		}

		/// @brief Resizes the matrix with optional data preservation.
		/// @param newDim New matrix dimension.
		/// @param preserveData If true, preserves overlapping elements.
		void Resize(int newDim, bool preserveData) {
			if (!preserveData || _dim == 0) {
				Resize(newDim);
				return;
			}

			if (newDim == _dim)
				return;

			if (newDim == 0) {
				_dim = 0;
				_data.clear();
				return;
			}

			validateDimension(newDim);

			// Create new storage
			std::vector<Type> newData(numElements(newDim), Type{0});

			// Copy existing elements where dimensions overlap
			int minDim = std::min(_dim, newDim);
			for (int i = 0; i < minDim; ++i) {
				for (int j = 0; j <= i; ++j) {
					size_t oldIdx = linearIndex(i, j);
					size_t newIdx = static_cast<size_t>(i) * (i + 1) / 2 + j;
					newData[newIdx] = _data[oldIdx];
				}
			}

			_dim = newDim;
			_data = std::move(newData);
		}

		/// @brief Sets all elements to zero.
		void SetToZero() noexcept { std::fill(_data.begin(), _data.end(), Type{0}); }

		/// @brief Sets all elements to a specified value.
		void SetToValue(Type val) noexcept { std::fill(_data.begin(), _data.end(), val); }

		/// @brief Tolerance-based equality comparison.
		/// @param b Matrix to compare with.
		/// @param eps Maximum allowed element difference.
		/// @return True if dimensions match and all elements are within tolerance.
		bool IsEqualTo(const MatrixSym& b, Type eps = Defaults::MatrixIsEqualTolerance) const {
			if (_dim != b._dim)
				return false;

			for (size_t i = 0; i < _data.size(); ++i) {
				if (Abs(_data[i] - b._data[i]) > eps)
					return false;
			}
			return true;
		}

		/// @brief Static helper for tolerance-based comparison.
		static bool AreEqual(const MatrixSym& a, const MatrixSym& b, Type eps = Defaults::MatrixIsEqualTolerance) {
			return a.IsEqualTo(b, eps);
		}

		/// @brief Converts to full Matrix representation.
		/// @return n×n Matrix with symmetric values filled in.
		Matrix<Type> GetAsMatrix() const {
			if (_dim == 0)
				return Matrix<Type>();

			Matrix<Type> ret(_dim, _dim);
			for (int i = 0; i < _dim; ++i) {
				for (int j = 0; j < _dim; ++j) {
					ret(i, j) = (*this)(i, j);
				}
			}
			return ret;
		}

		/// @brief Creates symmetric matrix from lower triangular part of a matrix.
		/// @param m Source matrix (must be square).
		/// @return MatrixSym containing m's lower triangular elements.
		/// @throws MatrixDimensionError If m is not square.
		static MatrixSym FromLower(const Matrix<Type>& m) {
			if (m.rows() != m.cols())
				throw MatrixDimensionError("MatrixSym::FromLower - must be square matrix", m.rows(), m.cols(), -1, -1);

			int dim = m.rows();
			if (dim == 0)
				return MatrixSym();

			MatrixSym ret(dim);
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ret(i, j) = m(i, j);
				}
			}
			return ret;
		}

		/// @brief Creates symmetric matrix from upper triangular part of a matrix.
		/// @param m Source matrix (must be square).
		/// @return MatrixSym containing m's upper triangular elements (transposed).
		/// @throws MatrixDimensionError If m is not square.
		static MatrixSym FromUpper(const Matrix<Type>& m) {
			if (m.rows() != m.cols())
				throw MatrixDimensionError("MatrixSym::FromUpper - must be square matrix", m.rows(), m.cols(), -1, -1);

			int dim = m.rows();
			if (dim == 0)
				return MatrixSym();

			MatrixSym ret(dim);
			for (int i = 0; i < dim; ++i) {
				for (int j = i; j < dim; ++j) {
					// a(i,j) from upper = a(j,i) in storage = a(i,j) symmetric
					ret(j, i) = m(i, j); // Transpose to lower
				}
			}
			return ret;
		}

		/// @brief Creates symmetric matrix as (M + Mᵀ) / 2.
		/// @param m Source matrix (must be square).
		/// @return MatrixSym containing symmetrized version of m.
		/// @throws MatrixDimensionError If m is not square.
		/// @note This is the standard way to symmetrize an arbitrary matrix.
		static MatrixSym FromFullMatrix(const Matrix<Type>& m) {
			if (m.rows() != m.cols())
				throw MatrixDimensionError("MatrixSym::FromFullMatrix - must be square matrix", m.rows(), m.cols(), -1, -1);

			int dim = m.rows();
			if (dim == 0)
				return MatrixSym();

			MatrixSym ret(dim);
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ret(i, j) = (m(i, j) + m(j, i)) / Type{2};
				}
			}
			return ret;
		}

		/// @brief Creates identity matrix of given dimension.
		/// @param dim Matrix dimension.
		/// @return dim×dim identity matrix (1s on diagonal).
		static MatrixSym Identity(int dim) {
			MatrixSym ret(dim); // Initializes to zero
			for (int i = 0; i < dim; ++i) {
				ret(i, i) = Type{1};
			}
			return ret;
		}

		/// @brief Creates diagonal matrix from vector.
		/// @param diag Vector of diagonal elements.
		/// @return dim×dim diagonal matrix where dim = diag.size().
		static MatrixSym GetDiagonalMatrix(const Vector<Type>& diag) {
			int dim = static_cast<int>(diag.size());
			MatrixSym ret(dim); // Initializes to zero
			for (int i = 0; i < dim; ++i) {
				ret(i, i) = diag[i];
			}
			return ret;
		}

		/// @brief Extracts a row as a vector.
		/// @param rowInd Row index (0-based).
		/// @return Vector of _dim elements.
		/// @throws MatrixAccessBoundsError If rowInd out of range.
		Vector<Type> VectorFromRow(int rowInd) const {
			if (rowInd < 0 || rowInd >= _dim)
				throw MatrixAccessBoundsError("VectorFromRow - row index out of bounds", rowInd, 0, _dim, _dim);

			Vector<Type> ret(_dim);
			for (int j = 0; j < _dim; ++j) {
				ret[j] = (*this)(rowInd, j);
			}
			return ret;
		}

		/// @brief Extracts a column as a vector.
		/// @param colInd Column index (0-based).
		/// @return Vector of _dim elements.
		/// @throws MatrixAccessBoundsError If colInd out of range.
		/// @note For symmetric matrices, VectorFromColumn(i) == VectorFromRow(i).
		Vector<Type> VectorFromColumn(int colInd) const {
			if (colInd < 0 || colInd >= _dim)
				throw MatrixAccessBoundsError("VectorFromColumn - column index out of bounds", 0, colInd, _dim, _dim);

			Vector<Type> ret(_dim);
			for (int i = 0; i < _dim; ++i) {
				ret[i] = (*this)(i, colInd);
			}
			return ret;
		}

		/// @brief Extracts the main diagonal as a vector.
		Vector<Type> VectorFromDiagonal() const {
			Vector<Type> ret(_dim);
			for (int i = 0; i < _dim; ++i) {
				ret[i] = (*this)(i, i);
			}
			return ret;
		}

		/// @brief Computes trace (sum of diagonal elements).
		/// @return @f$ \text{tr}(A) = \sum_{i=0}^{n-1} a_{ii} @f$
		Type Trace() const noexcept {
			Type sum{0};
			for (int i = 0; i < _dim; ++i) {
				sum += (*this)(i, i);
			}
			return sum;
		}

		/// @brief Computes Frobenius norm.
		/// @return @f$ \|A\|_F = \sqrt{\sum_{i,j} |a_{ij}|^2} @f$
		/// @note Off-diagonal elements are counted twice (once from each triangle).
		Real NormFrobenius() const {
			Real sum{0}; // Norm is always real, even for complex matrices
			for (int i = 0; i < _dim; ++i) {
				// Diagonal element: use |a|^2
				if constexpr (std::is_same_v<Type, Complex> || std::is_same_v<Type, std::complex<float>> ||
											std::is_same_v<Type, std::complex<long double>>) {
					sum += std::norm((*this)(i, i));
					// Off-diagonal elements (counted twice in full matrix)
					for (int j = 0; j < i; ++j) {
						sum += Real{2} * std::norm((*this)(i, j));
					}
				} else {
					sum += (*this)(i, i) * (*this)(i, i);
					// Off-diagonal elements (counted twice in full matrix)
					for (int j = 0; j < i; ++j) {
						sum += Real{2} * (*this)(i, j) * (*this)(i, j);
					}
				}
			}
			return std::sqrt(sum);
		}

		/// @brief Computes infinity norm (maximum absolute row sum).
		/// @return @f$ \|A\|_\infty = \max_i \sum_j |a_{ij}| @f$
		Real NormInf() const {
			Real maxSum{0}; // Norm is always real
			for (int i = 0; i < _dim; ++i) {
				Real rowSum{0};
				for (int j = 0; j < _dim; ++j) {
					rowSum += Abs((*this)(i, j));
				}
				if (rowSum > maxSum)
					maxSum = rowSum;
			}
			return maxSum;
		}

		/// @brief Computes 1-norm (maximum absolute column sum).
		/// @return @f$ \|A\|_1 = \max_j \sum_i |a_{ij}| @f$
		/// @note For symmetric matrices, Norm1() == NormInf().

		Real Norm1() const {
			return NormInf(); // Symmetric property: row sum = column sum
		}

		/// @brief Copy assignment.
		MatrixSym& operator=(const MatrixSym& other) {
			if (this != &other) {
				_dim = other._dim;
				_data = other._data;
			}
			return *this;
		}

		/// @brief Move assignment.
		MatrixSym& operator=(MatrixSym&& other) noexcept {
			if (this != &other) {
				_dim = other._dim;
				_data = std::move(other._data);
				other._dim = 0;
			}
			return *this;
		}

		/// @brief Element access (const) without bounds checking.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Element at (i, j), automatically handling symmetry.
		Type operator()(int i, int j) const noexcept { return _data[linearIndex(i, j)]; }

		/// @brief Element access (non-const) without bounds checking.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to stored element (setting (i,j) also sets (j,i)).
		Type& operator()(int i, int j) noexcept { return _data[linearIndex(i, j)]; }

		/// @brief Bounds-checked element access (const).
		/// @param i Row index.
		/// @param j Column index.
		/// @return Element at (i, j).
		/// @throws MatrixAccessBoundsError If indices out of range.

		const Type& at(int i, int j) const {
			if (i < 0 || i >= _dim || j < 0 || j >= _dim)
				throw MatrixAccessBoundsError("MatrixSym::at", i, j, _dim, _dim);
			return _data[linearIndex(i, j)];
		}

		/// @brief Bounds-checked element access (non-const).
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to stored element.
		/// @throws MatrixAccessBoundsError If indices out of range.

		Type& at(int i, int j) {
			if (i < 0 || i >= _dim || j < 0 || j >= _dim)
				throw MatrixAccessBoundsError("MatrixSym::at", i, j, _dim, _dim);
			return _data[linearIndex(i, j)];
		}

		/// @brief Legacy alias for at() (const).
		Type ElemAt(int i, int j) const { return at(i, j); }

		/// @brief Legacy alias for at() (non-const).
		Type& ElemAt(int i, int j) { return at(i, j); }


		/// /** @name Arithmetic Operators

		/// @brief Matrix addition (returns symmetric result).
		MatrixSym operator+(const MatrixSym& b) const {
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator+ - dimensions must match", _dim, -1, b._dim, -1);

			MatrixSym ret(_dim);
			for (size_t i = 0; i < _data.size(); ++i) {
				ret._data[i] = _data[i] + b._data[i];
			}
			return ret;
		}

		/// @brief In-place matrix addition.
		MatrixSym& operator+=(const MatrixSym& b) {
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator+= - dimensions must match", _dim, -1, b._dim, -1);

			for (size_t i = 0; i < _data.size(); ++i) {
				_data[i] += b._data[i];
			}
			return *this;
		}

		/// @brief Unary negation.
		MatrixSym operator-() const {
			MatrixSym ret(_dim);
			for (size_t i = 0; i < _data.size(); ++i) {
				ret._data[i] = -_data[i];
			}
			return ret;
		}

		/// @brief Matrix subtraction (returns symmetric result).
		MatrixSym operator-(const MatrixSym& b) const {
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator- - dimensions must match", _dim, -1, b._dim, -1);

			MatrixSym ret(_dim);
			for (size_t i = 0; i < _data.size(); ++i) {
				ret._data[i] = _data[i] - b._data[i];
			}
			return ret;
		}

		/// @brief In-place matrix subtraction.
		MatrixSym& operator-=(const MatrixSym& b) {
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator-= - dimensions must match", _dim, -1, b._dim, -1);

			for (size_t i = 0; i < _data.size(); ++i) {
				_data[i] -= b._data[i];
			}
			return *this;
		}

		/// @brief Symmetric-symmetric matrix multiplication.
		/// @param b Right-hand symmetric matrix.
		/// @return Regular Matrix (product is NOT symmetric in general).
		/// @throws MatrixDimensionError If dimensions don't match.
		/// @note Returns Matrix, not MatrixSym, because A*B is not symmetric even if A and B are.
		Matrix<Type> operator*(const MatrixSym& b) const {
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator*(MatrixSym) - dimensions must match", _dim, _dim, b._dim, b._dim);

			Matrix<Type> ret(_dim, _dim);
			for (int i = 0; i < _dim; ++i) {
				for (int j = 0; j < _dim; ++j) {
					Type sum{0};
					for (int k = 0; k < _dim; ++k) {
						sum += (*this)(i, k) * b(k, j);
					}
					ret(i, j) = sum;
				}
			}
			return ret;
		}

		/// @brief Symmetric-matrix multiplication.
		/// @param b Right-hand general matrix.
		/// @return Regular Matrix of size (dim × b.cols()).
		/// @throws MatrixDimensionError If dimensions don't match.
		Matrix<Type> operator*(const Matrix<Type>& b) const {
			if (_dim != b.rows())
				throw MatrixDimensionError("MatrixSym::operator*(Matrix) - dimension mismatch", _dim, _dim, b.rows(), b.cols());

			Matrix<Type> ret(_dim, b.cols());
			for (int i = 0; i < _dim; ++i) {
				for (int j = 0; j < b.cols(); ++j) {
					Type sum{0};
					for (int k = 0; k < _dim; ++k) {
						sum += (*this)(i, k) * b(k, j);
					}
					ret(i, j) = sum;
				}
			}
			return ret;
		}

		/// @brief Scalar multiplication (A * scalar).
		friend MatrixSym operator*(const MatrixSym& a, Type scalar) {
			MatrixSym ret(a._dim);
			for (size_t i = 0; i < a._data.size(); ++i) {
				ret._data[i] = a._data[i] * scalar;
			}
			return ret;
		}

		/// @brief Scalar multiplication (scalar * A, commutative).
		friend MatrixSym operator*(Type scalar, const MatrixSym& a) { return a * scalar; }

		/// @brief In-place scalar multiplication.
		MatrixSym& operator*=(Type scalar) {
			for (auto& elem : _data) {
				elem *= scalar;
			}
			return *this;
		}

		/// @brief Scalar division (A / scalar).
		friend MatrixSym operator/(const MatrixSym& a, Type scalar) {
			MatrixSym ret(a._dim);
			for (size_t i = 0; i < a._data.size(); ++i) {
				ret._data[i] = a._data[i] / scalar;
			}
			return ret;
		}

		/// @brief In-place scalar division.
		MatrixSym& operator/=(Type scalar) {
			for (auto& elem : _data) {
				elem /= scalar;
			}
			return *this;
		}

		/// @brief Matrix-vector multiplication (A * v).
		/// @param a Symmetric matrix.
		/// @param v Column vector.
		/// @return Result vector A·v.
		/// @throws MatrixDimensionError If dimensions don't match.
		friend Vector<Type> operator*(const MatrixSym& a, const Vector<Type>& v) {
			if (a._dim != static_cast<int>(v.size()))
				throw MatrixDimensionError("operator*(MatrixSym, Vector) - dimension mismatch", a._dim, a._dim, static_cast<int>(v.size()), -1);

			Vector<Type> ret(a._dim);
			for (int i = 0; i < a._dim; ++i) {
				Type sum{0};
				for (int j = 0; j < a._dim; ++j) {
					sum += a(i, j) * v[j];
				}
				ret[i] = sum;
			}
			return ret;
		}

		/// @brief Row vector times matrix (vᵀ * A).
		/// @param v Row vector.
		/// @param a Symmetric matrix.
		/// @return Result vector vᵀ·A.
		/// @throws MatrixDimensionError If dimensions don't match.
		/// @note For symmetric A: v*A == A*v.
		friend Vector<Type> operator*(const Vector<Type>& v, const MatrixSym& a) {
			if (static_cast<int>(v.size()) != a._dim)
				throw MatrixDimensionError("operator*(Vector, MatrixSym) - dimension mismatch", static_cast<int>(v.size()), -1, a._dim, a._dim);

			Vector<Type> ret(a._dim);
			for (int j = 0; j < a._dim; ++j) {
				Type sum{0};
				for (int i = 0; i < a._dim; ++i) {
					sum += v[i] * a(i, j);
				}
				ret[j] = sum;
			}
			return ret;
		}

		/// @brief Computes the inverse of the matrix.
		/// @return Inverse as a general Matrix (inverse of symmetric is symmetric, but returned as Matrix).
		/// @throws MatrixDimensionError If matrix is empty.
		/// @throws SingularMatrixError If matrix is singular.
		/// @note Converts to Matrix, inverts, and returns result.
		Matrix<Type> GetInverse() const {
			if (_dim == 0)
				throw MatrixDimensionError("MatrixSym::GetInverse - cannot invert empty matrix", 0, 0, -1, -1);

			Matrix<Type> a = GetAsMatrix();
			a.Invert();
			return a;
		}

		/// @brief Converts matrix to formatted string.
		/// @param width Field width for each element (default 10).
		/// @param precision Decimal precision (default 3).
		/// @return Formatted string representation.
		std::string to_string(int width = 10, int precision = 3) const {
			std::stringstream str;
			Print(str, width, precision);
			return str.str();
		}

		/// @brief Prints matrix with configurable formatting.
		/// @param stream Output stream.
		/// @param fmt MatrixPrintFormat specifying width, precision, brackets, etc.
		/// @see MatrixPrintFormat For format options.
		void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const {
			if (fmt.showHeader) {
				stream << "Rows: " << _dim << " Cols: " << _dim << " (Symmetric)" << std::endl;
			}

			std::ios_base::fmtflags oldFlags = stream.flags();
			if (fmt.scientific)
				stream << std::scientific;
			else if (fmt.fixed)
				stream << std::fixed;

			// Compact mode for small matrices
			if (fmt.compactMode && _dim <= 3) {
				if (fmt.showBrackets)
					stream << "[";
				for (int i = 0; i < _dim; ++i) {
					if (i > 0)
						stream << "; ";
					for (int j = 0; j < _dim; ++j) {
						if (j > 0)
							stream << fmt.delimiter;
						stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << (*this)(i, j);
					}
				}
				if (fmt.showBrackets)
					stream << "]";
				stream << std::endl;
			} else {
				// Normal multi-line mode
				for (int i = 0; i < _dim; ++i) {
					if (fmt.showBrackets)
						stream << "[ ";
					for (int j = 0; j < _dim; ++j) {
						stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << (*this)(i, j);
						if (j < _dim - 1)
							stream << fmt.delimiter;
					}
					if (fmt.showBrackets)
						stream << " ]";
					if (i < _dim - 1)
						stream << std::endl;
				}
			}

			stream.flags(oldFlags);
		}

		/// @brief Legacy print method.
		/// @param stream Output stream.
		/// @param width Field width.
		/// @param precision Decimal precision.

		void Print(std::ostream& stream, int width, int precision) const {
			MatrixPrintFormat fmt = MatrixPrintFormat::Default();
			fmt.width = width;
			fmt.precision = precision;
			Print(stream, fmt);
		}

		/// @brief Stream output operator (default formatting).
		friend std::ostream& operator<<(std::ostream& stream, const MatrixSym& a) {
			a.Print(stream, 10, 3);
			return stream;
		}
	};
} // namespace MML
#endif
