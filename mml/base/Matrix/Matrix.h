///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Matrix.h                                                            ///
///  Description: Generic matrix class with arithmetic, decompositions, utilities     ///
///               Supports real and complex types, row/column operations              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_H
#define MML_MATRIX_H

#include <memory>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <complex>
#include <type_traits>

#include "MMLBase.h"
#include "base/MatrixPrintFormat.h"
#include "base/Vector/Vector.h"

namespace MML
{
	// Forward declarations
	template<class Type> class Matrix;
	template<class Type> class MatrixViewNew;

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               MatrixViewNew - Lightweight view                       ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/// @brief Lightweight non-owning view into a Matrix block
	/// @details Provides strided access to a rectangular subregion of a Matrix without copying.
	///          Used for efficient submatrix operations.
	/// @warning This view stores a raw pointer into the parent Matrix's internal storage.
	///          The view becomes invalid (dangling) if the parent Matrix is destroyed,
	///          resized, or otherwise reallocated. Do not store views beyond the
	///          parent's lifetime or across Resize() calls.
	/// @tparam Type Element type
	template<class Type>
	class MatrixViewNew {
		Type* _data;
		int _rows, _cols, _stride;
	public:
		/// @brief Construct view into matrix data
		/// @param data Pointer to first element of the view
		/// @param rows Number of rows in view
		/// @param cols Number of columns in view
		/// @param stride Row stride (number of columns in parent matrix)
		MatrixViewNew(Type* data, int rows, int cols, int stride)
			: _data(data), _rows(rows), _cols(cols), _stride(stride) {}

		/// @brief Element access (row-major with stride)
		Type& operator()(int i, int j) noexcept { return _data[i * _stride + j]; }
		/// @brief Element access (const)
		const Type& operator()(int i, int j) const noexcept { return _data[i * _stride + j]; }

		/// @brief Get number of rows (preferred API)
		int rows() const noexcept { return _rows; }
		/// @brief Get number of columns (preferred API)
		int cols() const noexcept { return _cols; }
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                                      Matrix                                          ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/// @brief Generic matrix class template for numerical computations
	/// @details Modern implementation using contiguous std::vector storage (row-major).
	///
	///          Key features:
	///          - std::vector for automatic memory management (RAII)
	///          - Exception-safe (no manual new/delete)
	///          - Better cache locality (flat storage vs Type**)
	///          - Move semantics automatic and efficient
	///          - Full arithmetic operations (+, -, *, /)
	///          - Matrix properties (trace, transpose, inverse)
	///          - Row/column extraction and manipulation
	///          - Iterator support for STL algorithms
	///
	/// @threadsafety Thread-safe for const operations. Non-const operations (including Resize)
	///               require external synchronization. See docs/THREADING.md for details.
	///
	/// @tparam Type Element type (typically Real, Complex, or arithmetic types)
	template<class Type>
	class Matrix
	{
	private:
		// Allocation safety limits
		static constexpr int MAX_DIMENSION = 100000;
		static constexpr size_t MAX_ELEMENTS = 100000000;  ///< ~800MB for doubles

		int _rows;
		int _cols;
		std::vector<Type> _data;  // Flat row-major storage

		// Index calculation (row-major)
		inline size_t idx(int i, int j) const { return static_cast<size_t>(i) * _cols + j; }

		// Validate dimensions for allocation
		void ValidateDimensions(int rows, int cols) const
		{
			if (rows < 0 || cols < 0)
				throw MatrixDimensionError("Matrix - dimensions cannot be negative", rows, cols, -1, -1);
			
			if (rows > MAX_DIMENSION || cols > MAX_DIMENSION)
				throw MatrixDimensionError("Matrix - dimensions exceed maximum limit", rows, cols, MAX_DIMENSION, MAX_DIMENSION);
			
			// Integer overflow protection
			size_t numElements = static_cast<size_t>(rows) * static_cast<size_t>(cols);
			if (cols > 0 && numElements / static_cast<size_t>(cols) != static_cast<size_t>(rows))
				throw std::overflow_error("Matrix - size calculation overflow");
			
			if (numElements > MAX_ELEMENTS)
				throw std::bad_alloc();
		}

		bool IsMatrixTypeComplex() const
		{
			return std::is_same_v<Type, std::complex<double>> ||
			       std::is_same_v<Type, std::complex<float>> ||
			       std::is_same_v<Type, std::complex<long double>>;
		}

	public:
		typedef Type value_type;  ///< Element type alias for STL compatibility

		///////////////////////          Constructors and destructor       //////////////////////
		
		/// @brief Default constructor - creates empty matrix
		Matrix() : _rows(0), _cols(0), _data() {}
		
		/// @brief Construct matrix of given dimensions
		/// @param rows Number of rows (must be non-negative)
		/// @param cols Number of columns (must be non-negative)
		/// @throws MatrixDimensionError if dimensions are negative or exceed limits
		/// @details Elements initialized to zero for numeric types
		explicit Matrix(int rows, int cols) : _rows(rows), _cols(cols)
		{
			if (rows < 0 || cols < 0) {
				throw MatrixDimensionError("Matrix: negative dimensions not allowed");
			}
			if (rows == 0 || cols == 0) {
				_rows = 0;
				_cols = 0;
				return;
			}
			ValidateDimensions(rows, cols);
			
			if constexpr (is_MML_simple_numeric<Type>) {
				_data.resize(static_cast<size_t>(rows) * cols, Type{0});
			} else {
				_data.resize(static_cast<size_t>(rows) * cols);
			}
		}
		
		/// @brief Construct matrix with uniform fill value
		/// @param rows Number of rows
		/// @param cols Number of columns
		/// @param val Value to fill all elements with
		/// @throws MatrixDimensionError if dimensions are invalid
		explicit Matrix(int rows, int cols, const Type& val) : _rows(rows), _cols(cols)
		{
			if (rows < 0 || cols < 0) {
				throw MatrixDimensionError("Matrix: negative dimensions not allowed");
			}
			if (rows == 0 || cols == 0) {
				_rows = 0;
				_cols = 0;
				return;
			}
			ValidateDimensions(rows, cols);
			_data.resize(static_cast<size_t>(rows) * cols, val);
		}
		
		/// @brief Construct from nested std::vector (2D array)
		/// @param values 2D vector where values[i][j] is element at row i, column j
		/// @throws MatrixDimensionError if rows have inconsistent sizes
		explicit Matrix(const std::vector<std::vector<Type>>& values)
		{
			if (values.empty() || values[0].empty()) {
				_rows = 0;
				_cols = 0;
				return;
			}
			_rows = static_cast<int>(values.size());
			_cols = static_cast<int>(values[0].size());
			ValidateDimensions(_rows, _cols);
			
			_data.reserve(static_cast<size_t>(_rows) * _cols);
			for (int i = 0; i < _rows; ++i) {
				if (static_cast<int>(values[i].size()) != _cols)
					throw MatrixDimensionError("Matrix - inconsistent row sizes", _rows, _cols, i, static_cast<int>(values[i].size()));
				for (int j = 0; j < _cols; ++j)
					_data.push_back(values[i][j]);
			}
		}
		
		/// @brief Construct from std::array (compile-time sized)
		/// @tparam Rows Number of rows (compile-time)
		/// @tparam Cols Number of columns (compile-time)
		template <std::size_t Rows, std::size_t Cols>
		explicit Matrix(const std::array<std::array<Type, Cols>, Rows>& arr)
			: _rows(static_cast<int>(Rows)), _cols(static_cast<int>(Cols))
		{
			ValidateDimensions(_rows, _cols);
			_data.reserve(static_cast<size_t>(_rows) * _cols);
			for (size_t i = 0; i < Rows; ++i)
				for (size_t j = 0; j < Cols; ++j)
					_data.push_back(arr[i][j]);
		}
		
		/// @brief Construct from raw pointer to contiguous data
		/// @param rows Number of rows
		/// @param cols Number of columns
		/// @param val Pointer to data array
		/// @param isRowWise If true, data is row-major; if false, column-major
		explicit Matrix(int rows, int cols, Type* val, bool isRowWise = true)
			: _rows(rows), _cols(cols)
		{
			if (rows < 0 || cols < 0) {
				throw MatrixDimensionError("Matrix: negative dimensions not allowed");
			}
			if (rows == 0 || cols == 0) {
				_rows = 0;
				_cols = 0;
				return;
			}
			ValidateDimensions(rows, cols);
			_data.resize(static_cast<size_t>(rows) * cols);
			
			if (isRowWise) {
				std::copy(val, val + static_cast<size_t>(rows) * cols, _data.begin());
			} else {
				// Column-wise: need to transpose during copy
				for (int j = 0; j < cols; ++j)
					for (int i = 0; i < rows; ++i)
						_data[idx(i, j)] = *val++;
			}
		}
		
		/// @brief Construct from initializer list (row-major order)
		/// @param rows Number of rows
		/// @param cols Number of columns
		/// @param values Initializer list with rows*cols elements
		/// @param strictMode If true, throws if list size doesn't match dimensions
		/// @throws MatrixDimensionError if strictMode and size mismatch
		explicit Matrix(int rows, int cols, std::initializer_list<Type> values, bool strictMode = true)
			: _rows(rows), _cols(cols)
		{
			if (rows < 0 || cols < 0) {
				throw MatrixDimensionError("Matrix: negative dimensions not allowed");
			}
			if (rows == 0 || cols == 0) {
				_rows = 0;
				_cols = 0;
				return;
			}
			ValidateDimensions(rows, cols);
			
			size_t totalSize = static_cast<size_t>(rows) * cols;
			if (strictMode && values.size() != totalSize)
				throw MatrixDimensionError("Matrix - initializer list size mismatch", rows, cols, -1, -1);
			
			if constexpr (is_MML_simple_numeric<Type>) {
				_data.resize(totalSize, Type{0});
			} else {
				_data.resize(totalSize);
			}
			
			auto it = values.begin();
			for (size_t i = 0; i < totalSize && it != values.end(); ++i, ++it)
				_data[i] = *it;
		}
		
		/// @brief Copy constructor
		Matrix(const Matrix& m) = default;
		
		/// @brief Move constructor
		/// @details Source matrix is left in valid empty state
		Matrix(Matrix&& m) noexcept 
			: _data(std::move(m._data)), _rows(m._rows), _cols(m._cols)
		{
			m._rows = 0;
			m._cols = 0;
		}
		
		/// @brief Submatrix constructor - extract rectangular region
		/// @param m Source matrix
		/// @param ind_row Starting row index
		/// @param ind_col Starting column index
		/// @param row_num Number of rows to extract
		/// @param col_num Number of columns to extract
		/// @throws MatrixDimensionError if region is out of bounds
		Matrix(const Matrix& m, int ind_row, int ind_col, int row_num, int col_num)
		{
			if (ind_row < 0 || ind_row >= m._rows || ind_col < 0 || ind_col >= m._cols)
				throw MatrixDimensionError("Matrix submatrix - invalid start indices", m._rows, m._cols, ind_row, ind_col);
			if (row_num <= 0 || col_num <= 0)
				throw MatrixDimensionError("Matrix submatrix - dimensions must be positive", row_num, col_num, -1, -1);
			if (ind_row + row_num > m._rows || ind_col + col_num > m._cols)
				throw MatrixDimensionError("Matrix submatrix - out of bounds", m._rows, m._cols, ind_row + row_num, ind_col + col_num);
			
			_rows = row_num;
			_cols = col_num;
			_data.reserve(static_cast<size_t>(row_num) * col_num);
			
			for (int i = 0; i < row_num; ++i)
				for (int j = 0; j < col_num; ++j)
					_data.push_back(m(ind_row + i, ind_col + j));
		}
		
		/// @brief Destructor (default - std::vector handles cleanup)
		~Matrix() = default;

		///////////////////////          Assignment operators              //////////////////////
		
		/// @brief Copy assignment
		Matrix& operator=(const Matrix& m) = default;
		
		/// @brief Move assignment
		Matrix& operator=(Matrix&& m) noexcept
		{
			if (this != &m) {
				_data = std::move(m._data);
				_rows = m._rows;
				_cols = m._cols;
				m._rows = 0;
				m._cols = 0;
			}
			return *this;
		}

		///////////////////////              Resize operations             //////////////////////
		
		/// @brief Resize the matrix
		/// @param rows New number of rows
		/// @param cols New number of columns
		/// @param preserveElements If true, existing elements are preserved (up to new dimensions)
		/// @throws MatrixDimensionError if dimensions are invalid
		void Resize(int rows, int cols, bool preserveElements = false)
		{
			if (rows == _rows && cols == _cols)
				return;
			
			if (rows <= 0 || cols <= 0)
				throw MatrixDimensionError("Matrix::Resize - dimensions must be positive", rows, cols, -1, -1);
			
			ValidateDimensions(rows, cols);
			
			if (preserveElements) {
				std::vector<Type> newData(static_cast<size_t>(rows) * cols, Type{0});
				int minRows = std::min(_rows, rows);
				int minCols = std::min(_cols, cols);
				for (int i = 0; i < minRows; ++i)
					for (int j = 0; j < minCols; ++j)
						newData[static_cast<size_t>(i) * cols + j] = _data[idx(i, j)];
				_data = std::move(newData);
			} else {
				if constexpr (is_MML_simple_numeric<Type>) {
					_data.assign(static_cast<size_t>(rows) * cols, Type{0});
				} else {
					_data.resize(static_cast<size_t>(rows) * cols);
				}
			}
			_rows = rows;
			_cols = cols;
		}
		
		/// @brief Convert to identity matrix (in-place)
		/// @throws MatrixDimensionError if matrix is not square
		void MakeUnitMatrix()
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::MakeUnitMatrix - must be square", _rows, _cols, -1, -1);
			
			std::fill(_data.begin(), _data.end(), Type{0});
			for (int i = 0; i < _rows; ++i)
				_data[idx(i, i)] = Type{1};
		}

		///////////////////////          Static factory methods            //////////////////////
		
		/// @brief Create identity matrix
		/// @param dim Matrix dimension (creates dim × dim identity)
		/// @return Identity matrix I where I(i,i) = 1, I(i,j) = 0 for i≠j
		static Matrix Identity(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("Matrix::Identity - dimension must be positive", dim, dim, -1, -1);
			
			Matrix mat(dim, dim);
			mat.MakeUnitMatrix();
			return mat;
		}
		
		/// @brief Create diagonal matrix from vector
		/// @param diagValues Vector of diagonal elements
		/// @return Diagonal matrix D where D(i,i) = diagValues[i]
		static Matrix Diagonal(const Vector<Type>& diagValues)
		{
			int n = diagValues.size();
			if (n <= 0)
				throw MatrixDimensionError("Matrix::Diagonal - vector size must be positive", n, n, -1, -1);
			
			Matrix mat(n, n);
			for (int i = 0; i < n; ++i)
				mat(i, i) = diagValues[i];
			return mat;
		}

		///////////////////////              Standard stuff                //////////////////////
		
		/// @brief Get number of rows (preferred API)
		inline int rows() const noexcept { return _rows; }
		/// @brief Get number of columns (preferred API)
		inline int cols() const noexcept { return _cols; }
		/// @brief Check if matrix is empty (zero dimensions)
		inline bool isEmpty() const noexcept { return _rows == 0 || _cols == 0; }
		
		/// @brief Direct access to underlying data (for algorithms that need it)
		Type* data() noexcept { return _data.data(); }
		/// @brief Direct access to underlying data (const)
		const Type* data() const noexcept { return _data.data(); }

		//////////////////         Accessors         /////////////////
		
		/// @brief Check if matrix is square.
		inline bool isSquare() const noexcept { return _rows == _cols; }
		
		/// @brief Returns string representation of matrix (compact format).
		std::string to_string() const
		{
			std::ostringstream oss;
			oss << "[";
			for (int i = 0; i < _rows; ++i) {
				if (i > 0) oss << "; ";
				oss << "[";
				for (int j = 0; j < _cols; ++j) {
					if (j > 0) oss << ", ";
					oss << (*this)(i, j);
				}
				oss << "]";
			}
			oss << "]";
			return oss.str();
		}

		///////////////////////          Vector extraction                 //////////////////////
		
		/// @brief Extract a row as vector
		/// @param rowInd Row index (0-based)
		/// @return Vector containing row elements
		/// @throws MatrixAccessBoundsError if index out of bounds
		Vector<Type> VectorFromRow(int rowInd) const
		{
			if (rowInd < 0 || rowInd >= _rows)
				throw MatrixAccessBoundsError("VectorFromRow - invalid row index", rowInd, 0, _rows, _cols);
			
			Vector<Type> ret(_cols);
			for (int j = 0; j < _cols; ++j)
				ret[j] = (*this)(rowInd, j);
			return ret;
		}
		
		/// @brief Extract a column as vector
		/// @param colInd Column index (0-based)
		/// @return Vector containing column elements
		/// @throws MatrixAccessBoundsError if index out of bounds
		Vector<Type> VectorFromColumn(int colInd) const
		{
			if (colInd < 0 || colInd >= _cols)
				throw MatrixAccessBoundsError("VectorFromColumn - invalid column index", 0, colInd, _rows, _cols);
			
			Vector<Type> ret(_rows);
			for (int i = 0; i < _rows; ++i)
				ret[i] = (*this)(i, colInd);
			return ret;
		}
		
		/// @brief Extract diagonal as vector (for square matrices)
		/// @return Vector containing diagonal elements
		/// @throws MatrixDimensionError if matrix is not square
		/// @brief Extract diagonal elements as vector (preferred API)
		/// @return Vector containing diagonal elements
		/// @throws MatrixDimensionError if matrix is not square
		Vector<Type> diagonal() const
		{
			if (_rows != _cols)
				throw MatrixDimensionError("diagonal - must be square", _rows, _cols, -1, -1);
			
			Vector<Type> ret(_rows);
			for (int i = 0; i < _rows; ++i)
				ret[i] = (*this)(i, i);
			return ret;
		}

		///////////////////////          Matrix extraction                 //////////////////////
		
		/// @brief Extract lower triangular part
		/// @param includeDiagonal If true, includes diagonal; if false, strictly lower
		/// @return Lower triangular matrix
		/// @throws MatrixDimensionError if matrix is not square
		[[nodiscard]] Matrix GetLower(bool includeDiagonal = true) const
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::GetLower - must be square", _rows, _cols, -1, -1);
			
			Matrix ret(_rows, _cols);
			for (int i = 0; i < _rows; ++i) {
				int jEnd = includeDiagonal ? i + 1 : i;
				for (int j = 0; j < jEnd; ++j)
					ret(i, j) = (*this)(i, j);
			}
			return ret;
		}
		
		/// @brief Extract upper triangular part
		/// @param includeDiagonal If true, includes diagonal; if false, strictly upper
		/// @return Upper triangular matrix
		/// @throws MatrixDimensionError if matrix is not square
		[[nodiscard]] Matrix GetUpper(bool includeDiagonal = true) const
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::GetUpper - must be square", _rows, _cols, -1, -1);
			
			Matrix ret(_rows, _cols);
			for (int i = 0; i < _rows; ++i) {
				int jStart = includeDiagonal ? i : i + 1;
				for (int j = jStart; j < _cols; ++j)
					ret(i, j) = (*this)(i, j);
			}
			return ret;
		}
		
		/// @brief Extract rectangular submatrix
		/// @param start_row Starting row index
		/// @param start_col Starting column index
		/// @param row_num Number of rows
		/// @param col_num Number of columns
		/// @return Submatrix copy
		Matrix GetSubmatrix(int start_row, int start_col, int row_num, int col_num) const
		{
			return Matrix(*this, start_row, start_col, row_num, col_num);
		}

		///////////////////////          Row/column operations             //////////////////////
		
		/// @brief Initialize a row from a vector
		/// @param rowInd Target row index
		/// @param vec Source vector (must have size == ColNum())
		/// @throws MatrixAccessBoundsError if row index invalid
		/// @throws MatrixDimensionError if vector size doesn't match columns
		void InitRowWithVector(int rowInd, const Vector<Type>& vec)
		{
			if (rowInd < 0 || rowInd >= _rows)
				throw MatrixAccessBoundsError("InitRowWithVector - invalid row index", rowInd, 0, _rows, _cols);
			if (vec.size() != _cols)
				throw MatrixDimensionError("InitRowWithVector - vector size must match columns", _rows, _cols, vec.size(), -1);
			
			for (int j = 0; j < _cols; ++j)
				_data[idx(rowInd, j)] = vec[j];
		}
		
		/// @brief Initialize a column from a vector
		/// @param colInd Target column index
		/// @param vec Source vector (must have size == RowNum())
		/// @throws MatrixAccessBoundsError if column index invalid
		/// @throws MatrixDimensionError if vector size doesn't match rows
		void InitColWithVector(int colInd, const Vector<Type>& vec)
		{
			if (colInd < 0 || colInd >= _cols)
				throw MatrixAccessBoundsError("InitColWithVector - invalid column index", 0, colInd, _rows, _cols);
			if (vec.size() != _rows)
				throw MatrixDimensionError("InitColWithVector - vector size must match rows", _rows, _cols, vec.size(), -1);
			
			for (int i = 0; i < _rows; ++i)
				_data[idx(i, colInd)] = vec[i];
		}
		
		/// @brief Swap two rows
		/// @param k First row index
		/// @param l Second row index
		void SwapRows(int k, int l)
		{
			if (k < 0 || k >= _rows || l < 0 || l >= _rows)
				throw MatrixDimensionError("Matrix::SwapRows - invalid row index", _rows, _cols, k, l);
			if (k == l) return;
			
			for (int j = 0; j < _cols; ++j)
				std::swap(_data[idx(k, j)], _data[idx(l, j)]);
		}
		
		/// @brief Swap two columns
		/// @param k First column index
		/// @param l Second column index
		void SwapCols(int k, int l)
		{
			if (k < 0 || k >= _cols || l < 0 || l >= _cols)
				throw MatrixDimensionError("Matrix::SwapCols - invalid column index", _rows, _cols, k, l);
			if (k == l) return;
			
			for (int i = 0; i < _rows; ++i)
				std::swap(_data[idx(i, k)], _data[idx(i, l)]);
		}

		///////////////////////              Creating views                //////////////////////
		
		/// @brief Create mutable view of a rectangular block
		/// @warning The returned view is invalidated by Resize() or destruction of this matrix.
		/// @param startRow Starting row index
		/// @param startCol Starting column index
		/// @param numRows Number of rows in block
		/// @param numCols Number of columns in block
		/// @return MatrixViewNew providing strided access to block
		MatrixViewNew<Type> block(int startRow, int startCol, int numRows, int numCols)
		{
			if (startRow < 0 || startCol < 0 || numRows <= 0 || numCols <= 0 ||
			    startRow + numRows > _rows || startCol + numCols > _cols)
				throw MatrixDimensionError("Matrix::block - invalid parameters", _rows, _cols, startRow, startCol);
			return MatrixViewNew<Type>(&_data[idx(startRow, startCol)], numRows, numCols, _cols);
		}
		
		/// @brief Create const view of a rectangular block
		/// @warning The returned view is invalidated by Resize() or destruction of this matrix.
		const MatrixViewNew<Type> block(int startRow, int startCol, int numRows, int numCols) const
		{
			if (startRow < 0 || startCol < 0 || numRows <= 0 || numCols <= 0 ||
			    startRow + numRows > _rows || startCol + numCols > _cols)
				throw MatrixDimensionError("Matrix::block - invalid parameters", _rows, _cols, startRow, startCol);
			return MatrixViewNew<Type>(const_cast<Type*>(&_data[idx(startRow, startCol)]), numRows, numCols, _cols);
		}

		///////////////////////               Matrix properties            //////////////////////
		
		/// @brief Check if matrix is identity matrix
		/// @param eps Tolerance for floating-point comparison
		/// @return true if |M(i,i)-1| < eps and |M(i,j)| < eps for i≠j
		bool isIdentity(double eps = Defaults::IsMatrixUnitTolerance) const
		{
			if (_rows != _cols) return false;
			
			for (int i = 0; i < _rows; ++i) {
				for (int j = 0; j < _cols; ++j) {
					if (i == j) {
						if (Abs((*this)(i, j) - Type{1}) > eps) return false;
					} else {
						if (Abs((*this)(i, j)) > eps) return false;
					}
				}
			}
			return true;
		}

		/// @brief Check if matrix is diagonal
		/// @param eps Tolerance for off-diagonal elements
		/// @return true if |M(i,j)| < eps for all i≠j
		bool isDiagonal(double eps = Defaults::IsMatrixDiagonalTolerance) const
		{
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					if (i != j && Abs((*this)(i, j)) > eps)
						return false;
			return true;
		}

		/// @brief Check if matrix is diagonally dominant
		/// @return true if |M(i,i)| ≥ Σ_{j≠i} |M(i,j)| for all rows
		bool isDiagonallyDominant() const
		{
			for (int i = 0; i < _rows; ++i) {
				Real sum{0};  // Use Real for magnitude accumulation
				for (int j = 0; j < _cols; ++j)
					if (i != j)
						sum += Abs((*this)(i, j));
				if (Abs((*this)(i, i)) < sum)
					return false;
			}
			return true;
		}

		/// @brief Check if matrix is symmetric (M = Mᵀ)
		/// @return true if M(i,j) = M(j,i) for all i,j
		bool isSymmetric() const
		{
			if (_rows != _cols) return false;
			
			for (int i = 0; i < _rows; ++i)
				for (int j = i + 1; j < _cols; ++j)
					if ((*this)(i, j) != (*this)(j, i))
						return false;
			return true;
		}

		/// @brief Check if matrix is anti-symmetric (skew-symmetric)
		/// @return true if M(i,j) = -M(j,i) (or -conj(M(j,i)) for complex)
		bool isAntiSymmetric() const
		{
			if (_rows != _cols) return false;
			
			for (int i = 0; i < _rows; ++i) {
				for (int j = i + 1; j < _cols; ++j) {
					if (IsMatrixTypeComplex()) {
						if ((*this)(i, j) != -std::conj((*this)(j, i)))
							return false;
					} else {
						if ((*this)(i, j) != -(*this)(j, i))
							return false;
					}
				}
			}
			return true;
		}

		/// @brief Check if matrix is a zero matrix (all elements are zero within tolerance)
		/// @param eps Maximum allowed absolute value for an element to be considered zero
		/// @return true if all elements have absolute value <= eps
		bool isZero(double eps = Defaults::IsMatrixZeroTolerance) const
		{
			for (const auto& elem : _data)
				if (Abs(elem) > eps)
					return false;
			return true;
		}

		///////////////////////             Matrix norm calculations       //////////////////////
		
		/// @brief L1 norm (sum of absolute values of all elements)
		Real NormL1() const
		{
			Real norm{0};
			for (const auto& elem : _data)
				norm += Abs(elem);
			return norm;
		}
		
		/// @brief Frobenius norm: ||M||_F = √(Σ |M(i,j)|²)
		/// @details For complex matrices, uses std::norm(z) = |z|² = real² + imag²
		Real NormL2() const
		{
			Real norm{0};
			for (const auto& elem : _data) {
				if constexpr (std::is_same_v<Type, Complex> || 
				              std::is_same_v<Type, std::complex<float>> ||
				              std::is_same_v<Type, std::complex<long double>>) {
					norm += std::norm(elem);  // |z|² for complex
				} else {
					norm += elem * elem;  // x² for real
				}
			}
			return std::sqrt(norm);
		}
		
		/// @brief L∞ norm (maximum absolute value of any element)
		Real NormLInf() const
		{
			Real norm{0};
			for (const auto& elem : _data)
				norm = std::max(norm, static_cast<Real>(Abs(elem)));
			return norm;
		}

		///////////////////////               Access operators             //////////////////////
		
		/// @brief Primary element access - M(i,j)
		inline Type  operator()(int i, int j) const noexcept { return _data[idx(i, j)]; }
		/// @brief Primary element access (mutable)
		inline Type& operator()(int i, int j)       noexcept { return _data[idx(i, j)]; }
		
		/// @brief Legacy row access - M[i][j] syntax
		/// @details Returns pointer to row start for compatibility with mat[i][j] code
		inline Type* operator[](int i)             noexcept { return &_data[idx(i, 0)]; }
		/// @brief Legacy row access (const)
		inline const Type* operator[](int i) const noexcept { return &_data[idx(i, 0)]; }
		
		/// @brief Bounds-checked element access
		/// @throws MatrixAccessBoundsError if indices out of bounds
		const Type& at(int i, int j) const
		{
			if (i < 0 || i >= _rows || j < 0 || j >= _cols)
				throw MatrixAccessBoundsError("Matrix::at", i, j, _rows, _cols);
			return _data[idx(i, j)];
		}
		
		/// @brief Bounds-checked element access (mutable)
		/// @throws MatrixAccessBoundsError if indices out of bounds
		Type& at(int i, int j)
		{
			if (i < 0 || i >= _rows || j < 0 || j >= _cols)
				throw MatrixAccessBoundsError("Matrix::at", i, j, _rows, _cols);
			return _data[idx(i, j)];
		}

		///////////////////////             Iterator support               //////////////////////
		
		/// @brief Flat iterators - traverse all elements in row-major order
		auto begin()        { return _data.begin(); }
		auto end()          { return _data.end(); }
		auto begin() const  { return _data.begin(); }
		auto end() const    { return _data.end(); }
		auto cbegin() const { return _data.cbegin(); }
		auto cend() const   { return _data.cend(); }

		/// @brief Row iterator proxy for range-based for loops
		class row_range {
			Type* _ptr;
			int _cols;
		public:
			row_range(Type* ptr, int cols) : _ptr(ptr), _cols(cols) {}
			Type* begin() { return _ptr; }
			Type* end()   { return _ptr + _cols; }
		};
		
		class const_row_range {
			const Type* _ptr;
			int _cols;
		public:
			const_row_range(const Type* ptr, int cols) : _ptr(ptr), _cols(cols) {}
			const Type* begin() const { return _ptr; }
			const Type* end() const   { return _ptr + _cols; }
		};
		
		/// @brief Get iterable row range for range-based for
		/// @param i Row index
		/// @return row_range object for iteration
		row_range row(int i)
		{
			if (i < 0 || i >= _rows)
				throw MatrixAccessBoundsError("Matrix::row", i, 0, _rows, _cols);
			return row_range(&_data[idx(i, 0)], _cols);
		}
		
		/// @brief Get iterable row range (const)
		const_row_range row(int i) const
		{
			if (i < 0 || i >= _rows)
				throw MatrixAccessBoundsError("Matrix::row", i, 0, _rows, _cols);
			return const_row_range(&_data[idx(i, 0)], _cols);
		}

		/// @brief Column iterator (strided, forward-only)
		class col_iterator {
			Type* _ptr;
			int _stride;
		public:
			using iterator_category = std::forward_iterator_tag;
			using value_type = Type;
			using difference_type = std::ptrdiff_t;
			using pointer = Type*;
			using reference = Type&;

			col_iterator(Type* ptr, int stride) : _ptr(ptr), _stride(stride) {}
			col_iterator& operator++() { _ptr += _stride; return *this; }
			col_iterator operator++(int) { auto tmp = *this; _ptr += _stride; return tmp; }
			bool operator==(const col_iterator& other) const { return _ptr == other._ptr; }
			bool operator!=(const col_iterator& other) const { return _ptr != other._ptr; }
			Type& operator*() { return *_ptr; }
		};
		
		class col_range {
			Type* _ptr;
			int _rows, _stride;
		public:
			col_range(Type* ptr, int rows, int stride) : _ptr(ptr), _rows(rows), _stride(stride) {}
			col_iterator begin() { return col_iterator(_ptr, _stride); }
			col_iterator end()   { return col_iterator(_ptr + _rows * _stride, _stride); }
		};
		
		/// @brief Get iterable column range
		/// @param j Column index
		/// @return col_range object for strided iteration
		col_range col(int j)
		{
			if (j < 0 || j >= _cols)
				throw MatrixAccessBoundsError("Matrix::col", 0, j, _rows, _cols);
			return col_range(&_data[j], _rows, _cols);
		}

		///////////////////////              Equality operations           //////////////////////
		
		/// @brief Exact equality comparison
		bool operator==(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				return false;
			return _data == b._data;
		}
		
		/// @brief Exact inequality
		bool operator!=(const Matrix& b) const { return !(*this == b); }
		
		/// @brief Approximate equality with tolerance
		/// @param b Matrix to compare with
		/// @param eps Maximum allowed difference per element
		bool IsEqualTo(const Matrix& b, Type eps = Defaults::MatrixIsEqualTolerance) const
		{
			if (_rows != b._rows || _cols != b._cols)
				return false;
			
			for (size_t i = 0; i < _data.size(); ++i)
				if (Abs(_data[i] - b._data[i]) > eps)
					return false;
			return true;
		}
		
		/// @brief Static version of IsEqualTo
		static bool AreEqual(const Matrix& a, const Matrix& b, Type eps = Defaults::MatrixIsEqualTolerance)
		{
			return a.IsEqualTo(b, eps);
		}

		///////////////////////              Arithmetic operators          //////////////////////
		
		/// @brief Unary negation
		[[nodiscard]] Matrix operator-() const
		{
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = -_data[i];
			return ret;
		}
		
		/// @brief Matrix addition
		/// @throws MatrixDimensionError if dimensions don't match
		[[nodiscard]] Matrix operator+(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+ - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] + b._data[i];
			return ret;
		}
		
		/// @brief In-place matrix addition
		Matrix& operator+=(const Matrix& b)
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+= - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] += b._data[i];
			return *this;
		}
		
		/// @brief Matrix subtraction
		/// @throws MatrixDimensionError if dimensions don't match
		[[nodiscard]] Matrix operator-(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator- - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] - b._data[i];
			return ret;
		}
		
		/// @brief In-place matrix subtraction
		Matrix& operator-=(const Matrix& b)
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator-= - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] -= b._data[i];
			return *this;
		}
		
		/// @brief Matrix multiplication: C = A * B
		/// @throws MatrixDimensionError if A.cols != B.rows
		[[nodiscard]] Matrix operator*(const Matrix& b) const
		{
			if (_cols != b._rows)
				throw MatrixDimensionError("Matrix::operator* - a.cols must equal b.rows", _rows, _cols, b._rows, b._cols);
			
			Matrix ret(_rows, b._cols);
			for (int i = 0; i < _rows; ++i) {
				for (int j = 0; j < b._cols; ++j) {
					Type sum{0};
					for (int k = 0; k < _cols; ++k)
						sum += (*this)(i, k) * b(k, j);
					ret(i, j) = sum;
				}
			}
			return ret;
		}
		
		/// @brief Scalar multiplication: M * s
		[[nodiscard]] Matrix operator*(const Type& scalar) const
		{
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] * scalar;
			return ret;
		}
		
		/// @brief In-place scalar multiplication
		Matrix& operator*=(const Type& scalar)
		{
			for (auto& elem : _data)
				elem *= scalar;
			return *this;
		}
		
		/// @brief Scalar division: M / s
		[[nodiscard]] Matrix operator/(const Type& scalar) const
		{
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] / scalar;
			return ret;
		}
		
		/// @brief In-place scalar division
		Matrix& operator/=(const Type& scalar)
		{
			for (auto& elem : _data)
				elem /= scalar;
			return *this;
		}
		
		/// @brief Matrix-vector multiplication: M * v
		/// @throws MatrixDimensionError if cols != vector size
		Vector<Type> operator*(const Vector<Type>& v) const
		{
			if (_cols != v.size())
				throw MatrixDimensionError("Matrix * Vector - cols must equal vector size", _rows, _cols, v.size(), -1);
			
			Vector<Type> ret(_rows);
			for (int i = 0; i < _rows; ++i) {
				Type sum{0};
				for (int j = 0; j < _cols; ++j)
					sum += (*this)(i, j) * v[j];
				ret[i] = sum;
			}
			return ret;
		}
		
		/// @brief Scalar * Matrix (commutative)
		friend Matrix operator*(const Type& scalar, const Matrix& m)
		{
			return m * scalar;
		}
		
		/// @brief Vector * Matrix multiplication: vᵀ * M
		friend Vector<Type> operator*(const Vector<Type>& v, const Matrix& m)
		{
			if (v.size() != m._rows)
				throw MatrixDimensionError("Vector * Matrix - vector size must equal rows", v.size(), -1, m._rows, m._cols);
			
			Vector<Type> ret(m._cols);
			for (int j = 0; j < m._cols; ++j) {
				Type sum{0};
				for (int i = 0; i < m._rows; ++i)
					sum += v[i] * m(i, j);
				ret[j] = sum;
			}
			return ret;
		}

		///////////////////////            Trace, Inverse & Transpose      //////////////////////
		
		/// @brief Calculate matrix trace (sum of diagonal elements) - preferred API
		/// @return tr(M) = Σ M(i,i)
		/// @throws MatrixDimensionError if matrix is not square
		Type trace() const
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::trace - must be square", _rows, _cols, -1, -1);
			
			Type sum{0};
			for (int i = 0; i < _rows; ++i)
				sum += (*this)(i, i);
			return sum;
		}

		/// @brief Invert matrix in-place using Gauss-Jordan elimination
		/// @throws MatrixDimensionError if matrix is not square
		/// @throws SingularMatrixError if matrix is singular
		void Invert()
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::Invert - must be square", _rows, _cols, -1, -1);
			
			// Gauss-Jordan elimination with pivoting
			int n = _rows;
			Matrix& a = *this;
			std::vector<int> indxc(n), indxr(n), ipiv(n, 0);

			// Compute infinity norm for norm-scaled singularity threshold
			Real norm_a = 0.0;
			for (int ii = 0; ii < n; ii++)
				for (int jj = 0; jj < n; jj++) {
					Real abs_val = Abs(a(ii, jj));
					if (abs_val > norm_a) norm_a = abs_val;
				}
			Real singularity_threshold = std::numeric_limits<Real>::epsilon() * norm_a * n;

			for (int i = 0; i < n; ++i) {
				Real big{0};
				int irow = 0, icol = 0;
				
				// Find pivot
				for (int j = 0; j < n; ++j) {
					if (ipiv[j] != 1) {
						for (int k = 0; k < n; ++k) {
							if (ipiv[k] == 0) {
								if (Abs(a(j, k)) >= big) {
									big = Abs(a(j, k));
									irow = j;
									icol = k;
								}
							}
						}
					}
				}
				++ipiv[icol];
				
				// Swap rows if needed
				if (irow != icol) {
					for (int l = 0; l < n; ++l)
						std::swap(a(irow, l), a(icol, l));
				}
				
				indxr[i] = irow;
				indxc[i] = icol;
				
				if (Abs(a(icol, icol)) < singularity_threshold)
					throw SingularMatrixError("Matrix::Invert - Singular Matrix", Abs(a(icol, icol)));
				
				Type pivinv = Type{1} / a(icol, icol);
				a(icol, icol) = Type{1};
				for (int l = 0; l < n; ++l)
					a(icol, l) *= pivinv;
				
				for (int ll = 0; ll < n; ++ll) {
					if (ll != icol) {
						Type dum = a(ll, icol);
						a(ll, icol) = Type{0};
						for (int l = 0; l < n; ++l)
							a(ll, l) -= a(icol, l) * dum;
					}
				}
			}
			
			// Unscramble columns
			for (int l = n - 1; l >= 0; --l) {
				if (indxr[l] != indxc[l]) {
					for (int k = 0; k < n; ++k)
						std::swap(a(k, indxr[l]), a(k, indxc[l]));
				}
			}
		}
		
		/// @brief Get inverse matrix (returns copy) - preferred API
		/// @throws SingularMatrixError if matrix is singular
		[[nodiscard]] Matrix inverse() const
		{
			Matrix ret(*this);
			ret.Invert();
			return ret;
		}

		/// @brief Transpose matrix in-place (only for square matrices)
		/// @throws MatrixDimensionError if matrix is not square
		void Transpose()
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::Transpose - in-place requires square matrix", _rows, _cols, -1, -1);
			
			for (int i = 0; i < _rows; ++i)
				for (int j = i + 1; j < _cols; ++j)
					std::swap(_data[idx(i, j)], _data[idx(j, i)]);
		}
		
		/// @brief Get transpose matrix (returns copy, works for any dimensions) - preferred API
		/// @return Mᵀ where Mᵀ(i,j) = M(j,i)
		[[nodiscard]] Matrix transpose() const
		{
			Matrix ret(_cols, _rows);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					ret(j, i) = (*this)(i, j);
			return ret;
		}

		///////////////////////                    I/O                    //////////////////////
		
		/// @brief Print matrix to stream with formatting options
		/// @param stream Output stream
		/// @param fmt Format options (width, precision, brackets, etc.)
		void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const
		{
			if (isEmpty()) {
				if (fmt.showHeader)
					stream << "Rows: " << _rows << " Cols: " << _cols << " - Empty matrix" << std::endl;
				return;
			}

			if (fmt.showHeader)
				stream << "Rows: " << _rows << " Cols: " << _cols << std::endl;

			std::ios_base::fmtflags oldFlags = stream.flags();
			if (fmt.scientific)
				stream << std::scientific;
			else if (fmt.fixed)
				stream << std::fixed;

			for (int i = 0; i < _rows; ++i) {
				if (fmt.showBrackets) stream << "[ ";
				for (int j = 0; j < _cols; ++j) {
					stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << (*this)(i, j);
					if (j < _cols - 1) stream << fmt.delimiter;
				}
				if (fmt.showBrackets) stream << " ]";
				if (i < _rows - 1) stream << std::endl;
			}
			
			stream.flags(oldFlags);
		}
		
		/// @brief Print with specified width and precision
		void Print(std::ostream& stream, int width, int precision) const
		{
			MatrixPrintFormat fmt = MatrixPrintFormat::Default();
			fmt.width = width;
			fmt.precision = precision;
			fmt.fixed = false;
			Print(stream, fmt);
		}
		
		/// @brief Print with zero threshold (values below threshold printed as 0)
		void Print(std::ostream& stream, int width, int precision, Real threshold) const
		{
			// Create a temporary copy with thresholded values
			Matrix<Type> temp(*this);
			if constexpr (is_MML_simple_numeric<Type>) {
				for (int i = 0; i < _rows; ++i) {
					for (int j = 0; j < _cols; ++j) {
						if (std::abs(static_cast<Real>(temp(i, j))) < threshold) {
							temp(i, j) = Type{0};
						}
					}
				}
			}
			// Use the standard Print function for consistent formatting
			temp.Print(stream, width, precision);
		}
		
		/// @brief Stream output operator
		friend std::ostream& operator<<(std::ostream& stream, const Matrix& a)
		{
			a.Print(stream, 10, 3);
			return stream;
		}
		
		/// @brief Convert to string representation
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;
			Print(str, width, precision);
			return str.str();
		}
	};

	//////////////////////               Default Matrix typedefs       //////////////////////
	
	// Common type aliases
	typedef Matrix<int>     MatrixInt;      ///< Integer matrix
	typedef Matrix<float>   MatrixFlt;      ///< Single-precision matrix
	typedef Matrix<double>  MatrixDbl;      ///< Double-precision matrix
	typedef Matrix<Complex> MatrixComplex;  ///< Complex matrix

	// Short aliases
	typedef Matrix<int>     MatI;   ///< Short alias for MatrixInt
	typedef Matrix<float>   MatF;   ///< Short alias for MatrixFlt
	typedef Matrix<double>  MatD;   ///< Short alias for MatrixDbl
	typedef Matrix<Complex> MatC;   ///< Short alias for MatrixComplex

	// Verify noexcept move operations
	static_assert(std::is_nothrow_move_constructible_v<Matrix<double>>,
	              "Matrix<double> should be nothrow move constructible");
	static_assert(std::is_nothrow_move_assignable_v<Matrix<double>>,
	              "Matrix<double> should be nothrow move assignable");
}


#endif // MML_MATRIX_H
