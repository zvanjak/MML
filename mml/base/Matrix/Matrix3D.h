///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Matrix3D.h                                                          ///
///  Description: Three-dimensional matrix (tensor) class for volumetric data         ///
///               Contiguous memory layout for cache efficiency                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
/// @file Matrix3D.h
/// @brief Three-dimensional matrix (tensor) class for volumetric data.
/// @ingroup Base
/// @section matrix3d_overview Overview
/// Matrix3D is a template class for storing 3D array data with:
/// - **Contiguous memory layout**: Single allocation for cache efficiency
/// - **Dual access patterns**: Both bracket [i][j][k] and parenthesis (i,j,k)
/// - **Value semantics**: Copy/move constructors and assignment
/// - **Bounds checking**: Optional via MML_DEBUG macro
/// @section matrix3d_memory Memory Layout
/// The data is stored in a single contiguous block with row-major ordering:
/// @code
/// Linear index = i * (dim2 * dim3) + j * dim3 + k
/// @endcode
/// A pointer structure enables natural [i][j][k] indexing without pointer
/// arithmetic at each access:
/// @code
/// _v[i]       -> pointer to slice i
/// _v[i][j]    -> pointer to row j in slice i
/// _v[i][j][k] -> element at (i,j,k)
/// @endcode
/// @section matrix3d_usage Usage Examples
/// @code{.cpp}
/// // Create a 10x10x10 tensor initialized to zero
/// Matrix3D<double> tensor(10, 10, 10, 0.0);
/// // Element access
/// tensor[0][0][0] = 1.0;     // Bracket notation
/// tensor(1, 2, 3) = 2.5;     // Parenthesis notation
/// double val = tensor.At(5, 5, 5);  // Always bounds-checked
/// // Fill and apply
/// tensor.Fill(1.0);
/// tensor.Apply([](double x) { return x * 2.0; });
/// // Scalar operations
/// tensor *= 0.5;
/// auto scaled = tensor * 2.0;
/// // Element-wise operations
/// Matrix3D<double> other(10, 10, 10, 1.0);
/// auto sum = tensor + other;
/// // Direct data access for algorithms
/// double* data = tensor.Data();
/// int total = tensor.Size();  // = dim1 * dim2 * dim3
/// @endcode
/// @section matrix3d_performance Performance Considerations
/// - Use bracket notation [i][j][k] for nested loops with good locality
/// - Use Data() for bulk operations (SIMD, parallel algorithms)
/// - Iterate in k (innermost) → j → i order for cache efficiency
/// @see Tensor General tensor class (N-dimensional)
/// @see MatrixNM 2D matrix with similar memory layout

#if !defined MML_MATRIX_3D_H
#define MML_MATRIX_3D_H

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "MMLBase.h"
#include "MMLExceptions.h"

namespace MML {
	/// @brief Three-dimensional matrix (tensor) with contiguous memory layout.
	/// @ingroup Base
	/// @tparam Type Element type (typically double, float, int, or Complex).
	/// Memory is allocated as a single contiguous std::vector for cache efficiency
	/// and automatic memory management (no manual new/delete).
	/// Supports both [i][j][k] and (i,j,k) access patterns.
	/// @par Thread Safety
	/// Read operations are thread-safe. Write operations to different elements
	/// are safe from different threads. Structural changes (Resize, assignment)
	/// require external synchronization.

	template<class Type>
	class Matrix3D {
	private:
		int _dim1 = 0, _dim2 = 0, _dim3 = 0; ///< Dimensions (depth, rows, columns)
		std::vector<Type>    _data;   ///< Contiguous data block
		std::vector<Type*>   _rows;   ///< Row pointers for [i][j][k] access
		std::vector<Type**>  _slices; ///< Slice pointers (one per slice)

		/// @brief Builds pointer structure for [i][j][k] access from _data.
		void BuildPointers() {
			_rows.resize(_dim1 * _dim2);
			_slices.resize(_dim1);

			for (int i = 0; i < _dim1; ++i) {
				_slices[i] = _rows.data() + i * _dim2;
				for (int j = 0; j < _dim2; ++j)
					_rows[i * _dim2 + j] = _data.data() + (i * _dim2 + j) * _dim3;
			}
		}

		/// @brief Validates dimensions and returns total element count.
		/// @throws ArgumentError If any dimension ≤ 0.
		static int ValidateDims(int d1, int d2, int d3) {
			if (d1 <= 0 || d2 <= 0 || d3 <= 0)
				throw ArgumentError("Matrix3D dimensions must be positive");
			return d1 * d2 * d3;
		}

	public:
		/// /** @name Constructors and Destructor
		/// @{ */

		/// @brief Constructs an empty (0×0×0) tensor.
		Matrix3D() = default;

		/// @brief Constructs a tensor with given dimensions.
		/// @param d1 First dimension (depth/slices).
		/// @param d2 Second dimension (rows per slice).
		/// @param d3 Third dimension (columns per row).
		/// @throws ArgumentError If any dimension ≤ 0.
		/// @note Elements are value-initialized (zero for numeric types).
		Matrix3D(int d1, int d2, int d3)
			: _dim1(d1), _dim2(d2), _dim3(d3)
			, _data(ValidateDims(d1, d2, d3)) {
			BuildPointers();
		}

		/// @brief Constructs a tensor with given dimensions and initial value.
		/// @param d1 First dimension.
		/// @param d2 Second dimension.
		/// @param d3 Third dimension.
		/// @param initVal Value to initialize all elements.
		Matrix3D(int d1, int d2, int d3, const Type& initVal)
			: _dim1(d1), _dim2(d2), _dim3(d3)
			, _data(ValidateDims(d1, d2, d3), initVal) {
			BuildPointers();
		}

		/// @brief Copy constructor (deep copy).
		/// @note Data is deep-copied via std::vector. Pointers are rebuilt.
		Matrix3D(const Matrix3D& other)
			: _dim1(other._dim1), _dim2(other._dim2), _dim3(other._dim3)
			, _data(other._data) {
			if (!_data.empty())
				BuildPointers();
		}

		/// @brief Move constructor (transfers ownership).
		Matrix3D(Matrix3D&& other) noexcept
			: _dim1(other._dim1), _dim2(other._dim2), _dim3(other._dim3)
			, _data(std::move(other._data))
			, _rows(std::move(other._rows))
			, _slices(std::move(other._slices)) {
			other._dim1 = other._dim2 = other._dim3 = 0;
		}

		/// @brief Destructor (default — std::vector handles memory).
		~Matrix3D() = default;
		/// /** @} */


		/// /** @name Assignment Operators
		/// @{ */

		/// @brief Copy assignment (deep copy).
		Matrix3D& operator=(const Matrix3D& other) {
			if (this != &other) {
				_dim1 = other._dim1;
				_dim2 = other._dim2;
				_dim3 = other._dim3;
				_data = other._data;
				if (!_data.empty())
					BuildPointers();
				else {
					_rows.clear();
					_slices.clear();
				}
			}
			return *this;
		}

		/// @brief Move assignment (transfers ownership).
		Matrix3D& operator=(Matrix3D&& other) noexcept {
			if (this != &other) {
				_dim1 = other._dim1;
				_dim2 = other._dim2;
				_dim3 = other._dim3;
				_data = std::move(other._data);
				_rows = std::move(other._rows);
				_slices = std::move(other._slices);
				other._dim1 = other._dim2 = other._dim3 = 0;
			}
			return *this;
		}
		/// /** @} */


		/// /** @name Dimension Accessors
		/// @{ */

		/// @brief Gets the first dimension (depth/number of slices).
		int Dim1() const { return _dim1; }

		/// @brief Gets the second dimension (rows per slice).
		int Dim2() const { return _dim2; }

		/// @brief Gets the third dimension (columns per row).
		int Dim3() const { return _dim3; }

		/// @brief Gets the total number of elements.
		int Size() const { return _dim1 * _dim2 * _dim3; }

		/// @brief Returns true if the tensor has no elements.
		bool IsEmpty() const { return _data.empty(); }

		/// @brief Gets pointer to raw data (for algorithms needing contiguous memory).
		Type* Data() { return _data.data(); }

		/// @brief Gets const pointer to raw data.
		const Type* Data() const { return _data.data(); }
		/// /** @} */


		/// /** @name Element Access
		/// @{ */

		/// @brief Bracket notation access: m[i][j][k].
		/// @param i Slice index (first dimension).
		/// @return Pointer to slice i (allows chained [j][k] access).
		/// @note No bounds checking (use At() for checked access).
		Type** operator[](int i) { return _slices[i]; }

		/// @brief Const bracket access.
		const Type* const* operator[](int i) const { return _slices[i]; }

		/// @brief Parenthesis notation access: m(i, j, k).
		/// @param i Slice index.
		/// @param j Row index.
		/// @param k Column index.
		/// @return Reference to element at (i,j,k).
		/// @note Bounds-checked only if MML_DEBUG is defined.
		Type& operator()(int i, int j, int k) {
#ifdef MML_DEBUG
			if (i < 0 || i >= _dim1 || j < 0 || j >= _dim2 || k < 0 || k >= _dim3)
				throw IndexError("Matrix3D index out of bounds");
#endif
			return _data[(i * _dim2 + j) * _dim3 + k];
		}

		/// @brief Const parenthesis access.
		const Type& operator()(int i, int j, int k) const {
#ifdef MML_DEBUG
			if (i < 0 || i >= _dim1 || j < 0 || j >= _dim2 || k < 0 || k >= _dim3)
				throw IndexError("Matrix3D index out of bounds");
#endif
			return _data[(i * _dim2 + j) * _dim3 + k];
		}

		/// @brief Always bounds-checked element access.
		/// @param i Slice index.
		/// @param j Row index.
		/// @param k Column index.
		/// @return Reference to element at (i,j,k).
		/// @throws IndexError If any index is out of bounds.
		Type& At(int i, int j, int k) {
			if (i < 0 || i >= _dim1 || j < 0 || j >= _dim2 || k < 0 || k >= _dim3)
				throw IndexError("Matrix3D::At index out of bounds");
			return _data[(i * _dim2 + j) * _dim3 + k];
		}

		/// @brief Const bounds-checked access.
		const Type& At(int i, int j, int k) const {
			if (i < 0 || i >= _dim1 || j < 0 || j >= _dim2 || k < 0 || k >= _dim3)
				throw IndexError("Matrix3D::At index out of bounds");
			return _data[(i * _dim2 + j) * _dim3 + k];
		}
		/// /** @} */


		/// /** @name Utility Operations
		/// @{ */

		/// @brief Fills all elements with a value.
		/// @param value The value to fill.
		void Fill(const Type& value) { std::fill(_data.begin(), _data.end(), value); }

		/// @brief Applies a function to all elements in-place.
		/// @param func Unary function Type → Type.
		void Apply(std::function<Type(Type)> func) {
			for (auto& elem : _data)
				elem = func(elem);
		}

		/// @brief Resizes the tensor (discards existing data).
		/// @param d1 New first dimension.
		/// @param d2 New second dimension.
		/// @param d3 New third dimension.
		/// @note If dimensions match current, does nothing.
		void Resize(int d1, int d2, int d3) {
			if (d1 == _dim1 && d2 == _dim2 && d3 == _dim3)
				return;
			if (d1 > 0 && d2 > 0 && d3 > 0) {
				_dim1 = d1; _dim2 = d2; _dim3 = d3;
				_data.assign(d1 * d2 * d3, Type{});
				BuildPointers();
			} else {
				_dim1 = _dim2 = _dim3 = 0;
				_data.clear();
				_rows.clear();
				_slices.clear();
			}
		}

		/// @brief Resizes the tensor and initializes to a value.
		/// @param d1 New first dimension.
		/// @param d2 New second dimension.
		/// @param d3 New third dimension.
		/// @param initVal Value to initialize all elements.
		void Resize(int d1, int d2, int d3, const Type& initVal) {
			if (d1 > 0 && d2 > 0 && d3 > 0) {
				_dim1 = d1; _dim2 = d2; _dim3 = d3;
				_data.assign(d1 * d2 * d3, initVal);
				BuildPointers();
			} else {
				_dim1 = _dim2 = _dim3 = 0;
				_data.clear();
				_rows.clear();
				_slices.clear();
			}
		}
		/// /** @} */


		/// /** @name Comparison Operators
		/// @{ */

		/// @brief Checks equality (dimensions and all elements).
		bool operator==(const Matrix3D& other) const {
			if (_dim1 != other._dim1 || _dim2 != other._dim2 || _dim3 != other._dim3)
				return false;
			return _data == other._data;
		}

		/// @brief Checks inequality.
		bool operator!=(const Matrix3D& other) const { return !(*this == other); }
		/// /** @} */


		/// /** @name Scalar Operations
		/// @{ */

		/// @brief In-place scalar multiplication.
		Matrix3D& operator*=(const Type& scalar) {
			for (auto& elem : _data) elem *= scalar;
			return *this;
		}

		/// @brief In-place scalar division.
		Matrix3D& operator/=(const Type& scalar) {
			for (auto& elem : _data) elem /= scalar;
			return *this;
		}

		/// @brief In-place scalar addition.
		Matrix3D& operator+=(const Type& scalar) {
			for (auto& elem : _data) elem += scalar;
			return *this;
		}

		/// @brief Scalar multiplication returning new tensor.
		Matrix3D operator*(const Type& scalar) const {
			Matrix3D r(*this);
			r *= scalar;
			return r;
		}

		/// @brief Scalar division returning new tensor.
		Matrix3D operator/(const Type& scalar) const {
			Matrix3D r(*this);
			r /= scalar;
			return r;
		}

		/// @brief Scalar addition returning new tensor.
		Matrix3D operator+(const Type& scalar) const {
			Matrix3D r(*this);
			r += scalar;
			return r;
		}
		/// /** @} */


		/// /** @name Element-wise Operations
		/// @{ */

		/// @brief In-place element-wise addition.
		/// @param other Tensor with matching dimensions.
		/// @throws MatrixDimensionError If dimensions don't match.
		Matrix3D& operator+=(const Matrix3D& other) {
			if (_dim1 != other._dim1 || _dim2 != other._dim2 || _dim3 != other._dim3)
				throw MatrixDimensionError("Matrix3D dimensions must match for addition");
			for (int i = 0, n = Size(); i < n; ++i)
				_data[i] += other._data[i];
			return *this;
		}

		/// @brief In-place element-wise subtraction.
		/// @param other Tensor with matching dimensions.
		/// @throws MatrixDimensionError If dimensions don't match.
		Matrix3D& operator-=(const Matrix3D& other) {
			if (_dim1 != other._dim1 || _dim2 != other._dim2 || _dim3 != other._dim3)
				throw MatrixDimensionError("Matrix3D dimensions must match for subtraction");
			for (int i = 0, n = Size(); i < n; ++i)
				_data[i] -= other._data[i];
			return *this;
		}

		/// @brief Element-wise addition returning new tensor.
		Matrix3D operator+(const Matrix3D& other) const {
			Matrix3D r(*this);
			r += other;
			return r;
		}

		/// @brief Element-wise subtraction returning new tensor.
		Matrix3D operator-(const Matrix3D& other) const {
			Matrix3D r(*this);
			r -= other;
			return r;
		}
		/// /** @} */


		/// /** @name String Output
		/// @{ */

		/// @brief Returns brief string representation "Matrix3D[d1×d2×d3]".
		std::string ToString() const {
			std::ostringstream oss;
			oss << "Matrix3D[" << _dim1 << "×" << _dim2 << "×" << _dim3 << "]";
			return oss.str();
		}

		/// @brief Returns detailed string representation with all elements.
		/// @param precision Decimal places for floating-point output.
		std::string ToStringFull(int precision = 4) const {
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(precision);
			oss << "Matrix3D[" << _dim1 << "×" << _dim2 << "×" << _dim3 << "]:\n";
			for (int i = 0; i < _dim1; ++i) {
				oss << "Slice " << i << ":\n";
				for (int j = 0; j < _dim2; ++j) {
					oss << "  [";
					for (int k = 0; k < _dim3; ++k) {
						if (k > 0)
							oss << ", ";
						oss << _data[(i * _dim2 + j) * _dim3 + k];
					}
					oss << "]\n";
				}
			}
			return oss.str();
		}

		/// @brief Stream output operator.
		friend std::ostream& operator<<(std::ostream& os, const Matrix3D& m) { return os << m.ToString(); }
		/// /** @} */
	};

	/// /** @name Type Aliases
	/// @{ */

	using Matrix3Dd = Matrix3D<double>; ///< Double-precision 3D tensor
	using Matrix3Df = Matrix3D<float>;	///< Single-precision 3D tensor
	using Matrix3Di = Matrix3D<int>;	///< Integer 3D tensor
										/// /** @} */

} // namespace MML

#endif
