/// @file MatrixNM.h
/// @brief Fixed-size N×M matrix template with compile-time dimensions.
/// @details This file provides a statically-allocated matrix class optimized for
/// small matrices where dimensions are known at compile time. The template parameters
/// N and M specify the number of rows and columns respectively.
/// **Key Features:**
/// - Stack allocation: No dynamic memory, all storage in fixed array `Type _vals[N][M]`
/// - Compile-time optimization: Dimensions known at compile time enable loop unrolling
/// - Full value semantics: Copy constructors and assignment operators
/// - Rich linear algebra: Inverse, transpose, trace, matrix multiplication
/// - Type flexibility: Works with float, double, Complex, or any numeric type
/// **Storage Layout:**
/// - Row-major order: `_vals[i][j]` stores element at row i, column j
/// - Contiguous memory: `N*M` elements in a 2D C-style array
/// **When to Use:**
/// - Small matrices (2×2, 3×3, 4×4) in performance-critical code
/// - Transformation matrices in graphics/physics
/// - Rotation matrices, Jacobians, small linear systems
/// - When dimensions are fixed and known at compile time
/// **Comparison with Other Matrix Types:**
/// | Type | Allocation | Dimensions | Best For |
/// |------|------------|------------|----------|
/// | MatrixNM | Stack | Fixed (N×M) | Small matrices, transformations |
/// | Matrix | Heap | Dynamic | General computation |
/// | MatrixSym | Heap | Dynamic square | Symmetric matrices |
/// | MatrixBandDiag | Heap | Dynamic | Sparse banded systems |
/// **Usage Example:**
/// @code
/// using namespace MML;
/// // 3×3 rotation matrix
/// MatrixNM<double, 3, 3> R = {
/// {cos(θ), -sin(θ), 0},
/// {sin(θ),  cos(θ), 0},
/// {    0,       0,  1}
/// };
/// // Vector transformation
/// VectorN<double, 3> v = {1, 0, 0};
/// VectorN<double, 3> v_rot = R * v;
/// // Matrix operations
/// auto R_inv = R.GetInverse();
/// auto R_T = R.transpose();
/// double tr = R.Trace();
/// @endcode
/// **Type Aliases (Convenience):**
/// - `Mat22D`, `Mat33D`, `Mat44D` - Double precision
/// - `Mat22F`, `Mat33F`, `Mat44F` - Single precision
/// - `Mat22C`, `Mat33C`, `Mat44C` - Complex
/// @see Matrix For dynamic-size matrices
/// @see VectorN For fixed-size vectors
/// @see MatrixSym For symmetric matrices
/// @author Zvonimir Vanjak
/// @date 2024-2025

///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixNM.h                                                          ///
///  Description: Fixed-size NxM matrix template for compile-time dimensions          ///
///               Static allocation, optimized operations for small matrices          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIXNM_H
#define MML_MATRIXNM_H

#include "MMLBase.h"
#include "base/MatrixPrintFormat.h"
#include "base/Vector/VectorN.h"

#include <initializer_list>
#include <iomanip>
#include <iostream>

// Standard headers - include what we use
#include <sstream>
#include <string>
#include <vector>

namespace MML {
	/// @brief Fixed-size N×M matrix with compile-time dimensions.
	/// @details MatrixNM provides a statically-allocated matrix where dimensions N (rows)
	/// and M (columns) are template parameters. This enables:
	/// - **Zero heap allocation**: All storage on the stack
	/// - **Compile-time size checking**: Dimension mismatches caught at compile time
	/// - **Optimized operations**: Compiler can unroll loops for small matrices
	/// The matrix uses row-major storage in a 2D array `_vals[N][M]`.
	/// **Matrix Multiplication Dimensions:**
	/// For A (N×M) * B (M×K) = C (N×K), the template deduces the result type automatically.
	/// @tparam Type Element type (Real, float, Complex, etc.)
	/// @tparam N Number of rows (compile-time constant)
	/// @tparam M Number of columns (compile-time constant)

	template<class Type, int N, int M>
	class MatrixNM {
		static_assert(N > 0 && M > 0, "MatrixNM dimensions must be positive");
	private:
		Type _vals[N][M] = {{0}}; ///< Row-major storage array (stack allocated)

		template<class U, int P, int Q> friend class MatrixNM;

	public:
		typedef Type value_type; ///< Element type alias for STL compatibility

		/// /** @name Constructors
		/// @{ */


		/// @brief Default constructor, initializes all elements to zero.
		MatrixNM() {}

		/// @brief Constructs from flat initializer list (row-major order).
		/// @param values Elements in row-major order.
		/// @code
		/// MatrixNM<double, 2, 3> A = {1, 2, 3, 4, 5, 6};
		/// // A = [1, 2, 3]
		/// //     [4, 5, 6]
		/// @endcode

		MatrixNM(std::initializer_list<Type> values) {
			auto val = values.begin();
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					if (val != values.end())
						_vals[i][j] = *val++;
					else
						_vals[i][j] = 0.0;
		}

		/// @brief Constructs from nested initializer list (row-wise).
		/// @param rows Initializer list of row initializer lists.
		/// @code
		/// MatrixNM<double, 2, 3> A = {
		/// {1, 2, 3},
		/// {4, 5, 6}
		/// };
		/// @endcode

		MatrixNM(std::initializer_list<std::initializer_list<Type>> rows) {
			size_t i = 0;
			for (auto rowIt = rows.begin(); rowIt != rows.end() && i < N; ++rowIt, ++i) {
				size_t j = 0;
				for (auto colIt = rowIt->begin(); colIt != rowIt->end() && j < M; ++colIt, ++j) {
					_vals[i][j] = *colIt;
				}
				// Fill remaining columns with zero if not enough elements
				for (; j < M; ++j) {
					_vals[i][j] = Type{0};
				}
			}
			// Fill remaining rows with zero if not enough rows
			for (; i < N; ++i) {
				for (size_t j = 0; j < M; ++j) {
					_vals[i][j] = Type{0};
				}
			}
		}

		/// @brief Constructs from flat C-style array.
		/// @param arr Pointer to array of elements (row-major).
		/// @param len Number of elements in array.

		MatrixNM(const Type* arr, size_t len) {
			size_t idx = 0;
			for (size_t i = 0; i < RowNum(); ++i) {
				for (size_t j = 0; j < ColNum(); ++j) {
					if (idx < len)
						_vals[i][j] = arr[idx++];
					else
						_vals[i][j] = Type{0};
				}
			}
		}

		/// @brief Copy constructor.
		MatrixNM(const MatrixNM& m) {
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m._vals[i][j];
		}

		/// @brief Constructs diagonal matrix from scalar.
		/// @param m Scalar value for all diagonal elements.
		/// @code
		/// MatrixNM<double, 3, 3> I(1.0);  // Identity matrix
		/// @endcode

		MatrixNM(const Type& m) {
			for (int i = 0; i < N; i++)
				_vals[i][i] = Type{m};
		}
		/// /** @} */


		/// /** @name Dimension Accessors
		/// @{ */


		/// @brief Returns number of rows (compile-time constant N) - preferred API.
		int rows() const noexcept { return N; }

		/// @brief Returns number of columns (compile-time constant M) - preferred API.
		int cols() const noexcept { return M; }

		/// @brief Returns number of rows (compile-time constant N).
		int RowNum() const { return N; }

		/// @brief Returns number of columns (compile-time constant M).
		int ColNum() const { return M; }
		/// /** @} */


		/// /** @name Factory Methods
		/// @{ */


		/// @brief Creates an identity matrix.
		/// @return N×M matrix with 1s on diagonal, 0s elsewhere.

		static MatrixNM Identity() {
			MatrixNM unitMat;

			for (int i = 0; i < N; i++)
				unitMat._vals[i][i] = 1.0;

			return unitMat;
		}

		/// @brief Converts this matrix to identity (in-place).
		/// @throws MatrixDimensionError If matrix is not square.

		void MakeUnitMatrix(void) {
			if (RowNum() == ColNum()) {
				for (int i = 0; i < RowNum(); i++)
					for (int j = 0; j < ColNum(); j++)
						if (i == j)
							_vals[i][j] = 1;
						else
							_vals[i][j] = 0;
			} else
				throw MatrixDimensionError("MatrixNM::MakeUnitMatrix - must be square matrix", N, M, -1, -1);
		}

		/// @brief Extracts lower triangular part.
		/// @param includeDiagonal If true, includes main diagonal; if false, strictly lower.
		/// @return Lower triangular matrix with zeros above.
		/// @throws MatrixDimensionError If not square.

		MatrixNM GetLower(bool includeDiagonal = true) const {
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetLower - must be square matrix", N, M, -1, -1);

			MatrixNM ret;
			for (int i = 0; i < RowNum(); i++) {
				if (includeDiagonal)
					for (int j = 0; j < i + 1; j++)
						ret[i][j] = _vals[i][j];
				else
					for (int j = 0; j < i; j++)
						ret[i][j] = _vals[i][j];
			}

			return ret;
		}

		/// @brief Extracts upper triangular part.
		/// @param includeDiagonal If true, includes main diagonal; if false, strictly upper.
		/// @return Upper triangular matrix with zeros below.
		/// @throws MatrixDimensionError If not square.

		MatrixNM GetUpper(bool includeDiagonal = true) const {
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetUpper - must be square matrix", N, M, -1, -1);

			MatrixNM ret;
			for (int i = 0; i < RowNum(); i++) {
				if (includeDiagonal)
					for (int j = i; j < ColNum(); j++)
						ret[i][j] = _vals[i][j];
				else
					for (int j = i + 1; j < ColNum(); j++)
						ret[i][j] = _vals[i][j];
			}

			return ret;
		}
		/// /** @} */


		/// /** @name Vector-Matrix Conversion
		/// @{ */


		/// @brief Extracts a row as a vector.
		/// @param rowInd Row index (0-based).
		/// @return Vector of M elements from the specified row.

		VectorN<Type, M> VectorFromRow(int rowInd) {
			VectorN<Type, M> ret;
			for (int j = 0; j < M; j++)
				ret[j] = _vals[rowInd][j];

			return ret;
		}

		/// @brief Extracts a column as a vector.
		/// @param colInd Column index (0-based).
		/// @return Vector of N elements from the specified column.

		VectorN<Type, N> VectorFromColumn(int colInd) {
			VectorN<Type, N> ret;
			for (int i = 0; i < N; i++)
				ret[i] = _vals[i][colInd];

			return ret;
		}

		/// @brief Extracts the main diagonal as a vector.
		/// @return Vector of min(N,M) elements from the diagonal.

		VectorN<Type, N> VectorFromDiagonal() {
			VectorN<Type, N> ret;
			for (int i = 0; i < N; i++)
				ret[i] = _vals[i][i];

			return ret;
		}
		/// /** @} */


		/// /** @name Assignment Operators
		/// @{ */


		/// @brief Copy assignment operator.
		MatrixNM& operator=(const MatrixNM& m) {
			if (this == &m)
				return *this;

			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m._vals[i][j];

			return *this;
		}

		/// @brief Scalar broadcast assignment (sets all elements to scalar).
		MatrixNM& operator=(const Type& m) {
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m;

			return *this;
		}
		/// /** @} */


		/// /** @name Element Access
		/// @{ */


		/// @brief Row access operator (const).
		/// @param i Row index.
		/// @return Pointer to row i (enables `mat[i][j]` syntax).
		/// @warning No bounds checking.

		inline const Type* operator[](int i) const { return _vals[i]; }

		/// @brief Row access operator (non-const).
		/// @param i Row index.
		/// @return Pointer to row i (enables `mat[i][j] = value` syntax).
		/// @warning No bounds checking.

		inline Type* operator[](int i) { return _vals[i]; }

		/// @brief Element access by indices (const).
		/// @param i Row index.
		/// @param j Column index.
		/// @return Element at (i, j).
		/// @warning No bounds checking.

		inline Type operator()(int i, int j) const { return _vals[i][j]; }

		/// @brief Element access by indices (non-const).
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element at (i, j).
		/// @warning No bounds checking.

		inline Type& operator()(int i, int j) { return _vals[i][j]; }

		/// @brief Bounds-checked element access (const).
		/// @param i Row index.
		/// @param j Column index.
		/// @return Element at (i, j).
		/// @throws MatrixAccessBoundsError If indices out of range.

		Type ElemAt(int i, int j) const {
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

			return _vals[i][j];
		}

		/// @brief Bounds-checked element access (non-const).
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element at (i, j).
		/// @throws MatrixAccessBoundsError If indices out of range.

		Type& ElemAt(int i, int j) {
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

			return _vals[i][j];
		}
		/// /** @} */


		/// /** @name Arithmetic Operators
		/// @{ */


		/// @brief Unary negation (returns -A).
		[[nodiscard]] MatrixNM operator-() {
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = -_vals[i][j];
			return temp;
		}

		/// @brief Matrix addition (A + B).
		[[nodiscard]] MatrixNM operator+(const MatrixNM& b) const {
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = b._vals[i][j] + _vals[i][j];
			return temp;
		}

		/// @brief Matrix subtraction (A - B).
		[[nodiscard]] MatrixNM operator-(const MatrixNM& b) const {
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = _vals[i][j] - b._vals[i][j];
			return temp;
		}

		/// @brief Matrix multiplication with compile-time dimension checking.
		/// @details For A (N×M) * B (M×K), returns C (N×K).
		/// @tparam K Number of columns in the right-hand matrix.
		/// @param b Right-hand matrix (M×K).
		/// @return Result matrix (N×K).
		/// @note Dimension compatibility is checked at compile time.

		template<int K>
		[[nodiscard]] MatrixNM<Type, N, K> operator*(const MatrixNM<Type, M, K>& b) const {
			MatrixNM<Type, N, K> ret;

			for (int i = 0; i < ret.rows(); i++)
				for (int j = 0; j < ret.cols(); j++) {
					ret._vals[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret._vals[i][j] += _vals[i][k] * b._vals[k][j];
				}

			return ret;
		}

		/// @brief Scalar multiplication (A * b).
		[[nodiscard]] MatrixNM operator*(const Type& b) const {
			int i, j;
			MatrixNM ret(*this);

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret._vals[i][j] *= b;

			return ret;
		}

		/// @brief In-place scalar multiplication (A *= b).
		MatrixNM& operator*=(const Type& b) {
			int i, j;

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					_vals[i][j] *= b;

			return *this;
		}

		/// @brief Scalar division (A / b).
		[[nodiscard]] MatrixNM operator/(const Type& b) const {
			int i, j;
			MatrixNM ret(*this);

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret._vals[i][j] /= b;

			return ret;
		}

		/// @brief Matrix-vector multiplication (A * v).
		/// @param b Column vector of M elements.
		/// @return Column vector of N elements.

		VectorN<Type, N> operator*(const VectorN<Type, M>& b) const {
			int i, j;
			VectorN<Type, N> ret;

			for (i = 0; i < N; i++) {
				ret[i] = 0;
				for (j = 0; j < M; j++)
					ret[i] += _vals[i][j] * b[j];
			}

			return ret;
		}

		/// @brief Scalar-matrix multiplication (a * B, commutative).
		friend MatrixNM operator*(Type a, const MatrixNM& b) {
			int i, j;
			MatrixNM ret;

			for (i = 0; i < b.rows(); i++)
				for (j = 0; j < b.cols(); j++)
					ret._vals[i][j] = a * b._vals[i][j];

			return ret;
		}

		/// @brief Row vector times matrix (v^T * A).
		/// @param a Row vector of N elements.
		/// @param b Matrix (N×M).
		/// @return Row vector of M elements.

		friend VectorN<Type, M> operator*(const VectorN<Type, N>& a, const MatrixNM<Type, N, M>& b) {
			int i, j;
			VectorN<Type, M> ret;

			for (i = 0; i < M; i++) {
				ret[i] = 0;
				for (j = 0; j < N; j++)
					ret[i] += a[j] * b._vals[j][i];
			}

			return ret;
		}
		/// /** @} */


		/// /** @name Equality Comparisons
		/// @{ */


		/// @brief Exact equality comparison.
		bool operator==(const MatrixNM& b) const {
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
					if (_vals[i][j] != b._vals[i][j])
						return false;

			return true;
		}

		/// @brief Exact inequality comparison.
		bool operator!=(const MatrixNM& b) const { return !(*this == b); }

		/// @brief Tolerance-based equality comparison.
		/// @param b Matrix to compare with.
		/// @param eps Maximum allowed element difference.
		/// @return True if |A[i][j] - B[i][j]| ≤ eps for all i, j.

		bool IsEqualTo(const MatrixNM& b, Type eps = Defaults::MatrixIsEqualTolerance) const {
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (std::abs(_vals[i][j] - b._vals[i][j]) > eps)
						return false;

			return true;
		}

		/// @brief Static helper for tolerance-based comparison.
		static bool AreEqual(const MatrixNM& a, const MatrixNM& b, Type eps = Defaults::MatrixIsEqualTolerance) { return a.IsEqualTo(b, eps); }
		/// /** @} */


		/// /** @name Matrix Operations
		/// @{ */


		/// @brief Computes the trace (sum of diagonal elements).
		/// @return @f$ \text{tr}(A) = \sum_{i=0}^{n-1} a_{ii} @f$
		/// @throws MatrixDimensionError If matrix is not square.

		Type Trace() const {
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Trace - must be square matrix", N, M, -1, -1);

			Type sum = 0;
			for (int i = 0; i < RowNum(); i++)
				sum += _vals[i][i];

			return sum;
		}

		/// @brief Inverts this matrix in-place using Gauss-Jordan elimination.
		/// @throws MatrixDimensionError If matrix is not square.
		/// @throws SingularMatrixError If matrix is singular.
		/// @note Modifies this matrix to contain its inverse.

		void Invert() {
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Invert - must be square matrix", N, M, -1, -1);

			MatrixNM& a = *this;
			MatrixNM<Type, N, 1> b; // dummy rhs

			b(0, 0) = 1.0;

			int i, icol, irow, j, k, l, ll;
			Type big, dum, pivinv;

			int n = RowNum();
			int m = b.cols();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++)
				ipiv[j] = 0;

			// Compute infinity norm for norm-scaled singularity threshold
			Real norm_a = 0.0;
			for (int ii = 0; ii < n; ii++)
				for (int jj = 0; jj < n; jj++) {
					Real abs_val = std::abs(a._vals[ii][jj]);
					if (abs_val > norm_a) norm_a = abs_val;
				}
			Real singularity_threshold = std::numeric_limits<Real>::epsilon() * norm_a * n;

			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (std::abs(a._vals[j][k]) >= big) {
									big = std::abs(a._vals[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++)
						std::swap(a._vals[irow][l], a._vals[icol][l]);
					for (l = 0; l < m; l++)
						std::swap(b._vals[irow][l], b._vals[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (std::abs(a._vals[icol][icol]) < singularity_threshold)
					throw SingularMatrixError("MatrixNM::Invert, gaussj: Singular Matrix", std::abs(a._vals[icol][icol]));

				pivinv = 1.0 / a._vals[icol][icol];
				a._vals[icol][icol] = 1.0;
				for (l = 0; l < n; l++)
					a._vals[icol][l] *= pivinv;
				for (l = 0; l < m; l++)
					b._vals[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a._vals[ll][icol];
						a._vals[ll][icol] = 0.0;
						for (l = 0; l < n; l++)
							a._vals[ll][l] -= a._vals[icol][l] * dum;
						for (l = 0; l < m; l++)
							b._vals[ll][l] -= b._vals[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a._vals[k][indxr[l]], a._vals[k][indxc[l]]);
			}
		}

		/// @brief Returns the inverse without modifying this matrix.
		/// @return A⁻¹ such that A * A⁻¹ = I.
		/// @throws MatrixDimensionError If matrix is not square.
		/// @throws SingularMatrixError If matrix is singular.

		MatrixNM GetInverse() const {
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::GetInverse - must be square matrix", N, M, -1, -1);

			MatrixNM a(*this); // making a copy, where inverse will be stored at the end

			a.Invert();

			return a;
		}

		/// @brief Transposes this matrix in-place.
		/// @throws MatrixDimensionError If matrix is not square (cannot transpose in-place).

		void Transpose() {
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Transpose - inplace Transpose possible only for square  matrix", N, M, -1, -1);

			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = i + 1; j < ColNum(); j++)
					std::swap(_vals[i][j], _vals[j][i]);
		}

		/// @brief Returns the transpose without modifying this matrix (preferred API).
		/// @return M×N matrix Aᵀ where Aᵀ[j][i] = A[i][j].
		/// @note Returns different type for non-square matrices.
		[[nodiscard]] MatrixNM<Type, M, N> transpose() const {
			MatrixNM<Type, M, N> ret;
			for (size_t i = 0; i < cols(); i++)
				for (size_t j = 0; j < rows(); j++)
					ret._vals[i][j] = _vals[j][i];
			return ret;
		}

		/// /** @} */


		/// /** @name Output Methods
		/// @{ */


		/// @brief Converts matrix to formatted string.
		/// @param width Field width for each element.
		/// @param precision Decimal precision.
		/// @return Formatted string representation.

		std::string to_string(int width, int precision) const {
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}

		/// @brief Prints matrix with configurable formatting.
		/// @param stream Output stream.
		/// @param fmt MatrixPrintFormat specifying width, precision, brackets, etc.
		/// @see MatrixPrintFormat For format options.

		void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const {
			// Show header if requested
			if (fmt.showHeader) {
				stream << "Rows: " << RowNum() << " Cols: " << ColNum() << std::endl;
			}

			// Set formatting flags
			std::ios_base::fmtflags oldFlags = stream.flags();
			if (fmt.scientific)
				stream << std::scientific;
			else if (fmt.fixed)
				stream << std::fixed;

			// Compact mode - print on one line for small matrices
			if (fmt.compactMode && RowNum() <= 3 && ColNum() <= 3) {
				if (fmt.showBrackets)
					stream << "[";
				for (int i = 0; i < RowNum(); i++) {
					if (i > 0)
						stream << "; ";
					for (int j = 0; j < ColNum(); j++) {
						if (j > 0)
							stream << fmt.delimiter;
						stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << _vals[i][j];
					}
				}
				if (fmt.showBrackets)
					stream << "]";
				stream << std::endl;
			} else {
				// Normal multi-line mode
				for (int i = 0; i < RowNum(); i++) {
					if (fmt.showBrackets)
						stream << "[ ";
					for (int j = 0; j < ColNum(); j++) {
						stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << _vals[i][j];
						if (j < ColNum() - 1)
							stream << fmt.delimiter;
					}
					if (fmt.showBrackets)
						stream << " ]";
					if (i < RowNum() - 1)
						stream << std::endl;
				}
			}

			// Restore original flags
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

		/// @brief Print with zero threshold (small values displayed as zero).
		/// @param stream Output stream.
		/// @param width Field width.
		/// @param precision Decimal precision.
		/// @param zeroThreshold Values with |value| ≤ threshold shown as 0.

		void Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const {
			stream << "Rows: " << RowNum() << " Cols: " << ColNum();

			for (int i = 0; i < RowNum(); i++) {
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++) {
					Type value{0};
					if (Abs(_vals[i][j]) > zeroThreshold)
						value = _vals[i][j];

					if (j == ColNum() - 1)
						stream << std::setw(width) << std::setprecision(precision) << value;
					else
						stream << std::setw(width) << std::setprecision(precision) << value << ", ";
				}
				if (i == RowNum() - 1)
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}

		/// @brief Stream output operator (default formatting).
		friend std::ostream& operator<<(std::ostream& stream, const MatrixNM& a) {
			a.Print(stream, 10, 3);

			return stream;
		}
		/// /** @} */
	};

	/// /** @name MatrixNM Type Aliases
	/// Convenience typedefs for common matrix sizes and element types.
	/// @{


	/// @name Float Matrices (Legacy Names)
	/// @{
	typedef MatrixNM<float, 2, 2> Matrix22Flt; ///< 2×2 single precision
	typedef MatrixNM<float, 3, 3> Matrix33Flt; ///< 3×3 single precision
	typedef MatrixNM<float, 4, 4> Matrix44Flt; ///< 4×4 single precision
	/// @}

	/// @name Double Matrices (Legacy Names)
	/// @{
	typedef MatrixNM<Real, 2, 2> Matrix22Dbl; ///< 2×2 double precision
	typedef MatrixNM<Real, 3, 3> Matrix33Dbl; ///< 3×3 double precision
	typedef MatrixNM<Real, 4, 4> Matrix44Dbl; ///< 4×4 double precision
	/// @}

	/// @name Complex Matrices (Legacy Names)
	/// @{
	typedef MatrixNM<Complex, 2, 2> Matrix22Complex; ///< 2×2 complex
	typedef MatrixNM<Complex, 3, 3> Matrix33Complex; ///< 3×3 complex
	typedef MatrixNM<Complex, 4, 4> Matrix44Complex; ///< 4×4 complex
	/// @}

	/// @name Short Float Aliases
	/// @{
	typedef MatrixNM<float, 2, 2> Mat22F; ///< 2×2 float (short name)
	typedef MatrixNM<float, 3, 3> Mat33F; ///< 3×3 float (short name)
	typedef MatrixNM<float, 4, 4> Mat44F; ///< 4×4 float (short name)
	/// @}

	/// @name Short Double Aliases
	/// @{
	typedef MatrixNM<Real, 2, 2> Mat22D; ///< 2×2 double (short name)
	typedef MatrixNM<Real, 3, 3> Mat33D; ///< 3×3 double (short name)
	typedef MatrixNM<Real, 4, 4> Mat44D; ///< 4×4 double (short name)
	/// @}

	/// @name Short Complex Aliases
	/// @{
	typedef MatrixNM<Complex, 2, 2> Mat22C; ///< 2×2 complex (short name)
	typedef MatrixNM<Complex, 3, 3> Mat33C; ///< 3×3 complex (short name)
	typedef MatrixNM<Complex, 4, 4> Mat44C; ///< 4×4 complex (short name)
											/// @}
											/// /** @} */

} // namespace MML

#endif