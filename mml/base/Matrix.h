///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Matrix.h                                                            ///
///  Description: Generic matrix class with arithmetic, decompositions, utilities     ///
///               Supports real and complex types, row/column operations              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
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
#include "Vector.h"

namespace MML
{
	// Forward declarations
	template<class Type> class Matrix;
	template<class Type> class MatrixViewNew;

	///////////////////////////////////////////////////////////////////////////////////////////
	///                               MatrixViewNew - Lightweight view                       ///
	///////////////////////////////////////////////////////////////////////////////////////////
	template<class Type>
	class MatrixViewNew {
		Type* _data;
		int _rows, _cols, _stride;
	public:
		MatrixViewNew(Type* data, int rows, int cols, int stride)
			: _data(data), _rows(rows), _cols(cols), _stride(stride) {}

		// Element access (row-major with stride)
		Type& operator()(int i, int j) { return _data[i * _stride + j]; }
		const Type& operator()(int i, int j) const { return _data[i * _stride + j]; }

		int RowNum() const { return _rows; }
		int ColNum() const { return _cols; }
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	///                                      Matrix                                       ///
	///   Modern Matrix implementation using contiguous std::vector storage (row-major)      ///
	///                                                                                       ///
	///   Key improvements over legacy Matrix:                                               ///
	///   - std::vector for automatic memory management (RAII)                               ///
	///   - Exception-safe (no manual new/delete)                                            ///
	///   - Better cache locality (flat storage vs Type**)                                   ///
	///   - Move semantics automatic and efficient                                           ///
	///   - Compatible interface for drop-in replacement                                     ///
	///////////////////////////////////////////////////////////////////////////////////////////
	template<class Type>
	class Matrix
	{
	private:
		// Allocation safety limits
		static constexpr int MAX_DIMENSION = 100000;
		static constexpr size_t MAX_ELEMENTS = 100000000;  // ~800MB for doubles

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
		typedef Type value_type;

		///////////////////////          Constructors and destructor       //////////////////////
		
		// Default constructor - empty matrix
		Matrix() : _rows(0), _cols(0), _data() {}
		
		// Size constructor - zero-initialized for numeric types
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
		
		// Size + fill value constructor
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
		
		// Constructor from std::vector<std::vector<Type>>
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
		
		// Constructor from std::array<std::array<Type, Cols>, Rows>
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
		
		// Constructor from pointer to continuous data
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
		
		// Constructor from initializer list (row-major)
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
		
		// Copy constructor
		Matrix(const Matrix& m) = default;
		
		// Move constructor - explicitly reset source dimensions
		Matrix(Matrix&& m) noexcept 
			: _data(std::move(m._data)), _rows(m._rows), _cols(m._cols)
		{
			m._rows = 0;
			m._cols = 0;
		}
		
		// Submatrix constructor
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
		
		// Destructor - default is fine, std::vector handles cleanup
		~Matrix() = default;

		///////////////////////          Assignment operators              //////////////////////
		Matrix& operator=(const Matrix& m) = default;
		
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
		
		void MakeUnitMatrix()
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::MakeUnitMatrix - must be square", _rows, _cols, -1, -1);
			
			std::fill(_data.begin(), _data.end(), Type{0});
			for (int i = 0; i < _rows; ++i)
				_data[idx(i, i)] = Type{1};
		}

		///////////////////////          Static factory methods            //////////////////////
		static Matrix GetUnitMatrix(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("Matrix::GetUnitMatrix - dimension must be positive", dim, dim, -1, -1);
			
			Matrix mat(dim, dim);
			mat.MakeUnitMatrix();
			return mat;
		}
		
		static Matrix GetDiagonalMatrix(const Vector<Type>& diagValues)
		{
			int n = diagValues.size();
			if (n <= 0)
				throw MatrixDimensionError("Matrix::GetDiagonalMatrix - vector size must be positive", n, n, -1, -1);
			
			Matrix mat(n, n);
			for (int i = 0; i < n; ++i)
				mat(i, i) = diagValues[i];
			return mat;
		}

		///////////////////////              Standard stuff                //////////////////////
		inline int RowNum() const { return _rows; }
		inline int ColNum() const { return _cols; }
		inline bool IsEmpty() const { return _rows == 0 || _cols == 0; }
		
		// Direct access to underlying data (for algorithms that need it)
		Type* data() { return _data.data(); }
		const Type* data() const { return _data.data(); }

		///////////////////////          Vector extraction                 //////////////////////
		Vector<Type> VectorFromRow(int rowInd) const
		{
			if (rowInd < 0 || rowInd >= _rows)
				throw MatrixAccessBoundsError("VectorFromRow - invalid row index", rowInd, 0, _rows, _cols);
			
			Vector<Type> ret(_cols);
			for (int j = 0; j < _cols; ++j)
				ret[j] = (*this)(rowInd, j);
			return ret;
		}
		
		Vector<Type> VectorFromColumn(int colInd) const
		{
			if (colInd < 0 || colInd >= _cols)
				throw MatrixAccessBoundsError("VectorFromColumn - invalid column index", 0, colInd, _rows, _cols);
			
			Vector<Type> ret(_rows);
			for (int i = 0; i < _rows; ++i)
				ret[i] = (*this)(i, colInd);
			return ret;
		}
		
		Vector<Type> VectorFromDiagonal() const
		{
			if (_rows != _cols)
				throw MatrixDimensionError("VectorFromDiagonal - must be square", _rows, _cols, -1, -1);
			
			Vector<Type> ret(_rows);
			for (int i = 0; i < _rows; ++i)
				ret[i] = (*this)(i, i);
			return ret;
		}
		
		Vector<Type> GetDiagonal() const { return VectorFromDiagonal(); }

		///////////////////////          Matrix extraction                 //////////////////////
		Matrix GetLower(bool includeDiagonal = true) const
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
		
		Matrix GetUpper(bool includeDiagonal = true) const
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
		
		Matrix GetSubmatrix(int start_row, int start_col, int row_num, int col_num) const
		{
			return Matrix(*this, start_row, start_col, row_num, col_num);
		}

		///////////////////////          Row/column operations             //////////////////////
		void InitRowWithVector(int rowInd, const Vector<Type>& vec)
		{
			if (rowInd < 0 || rowInd >= _rows)
				throw MatrixAccessBoundsError("InitRowWithVector - invalid row index", rowInd, 0, _rows, _cols);
			if (vec.size() != _cols)
				throw MatrixDimensionError("InitRowWithVector - vector size must match columns", _rows, _cols, vec.size(), -1);
			
			for (int j = 0; j < _cols; ++j)
				_data[idx(rowInd, j)] = vec[j];
		}
		
		void InitColWithVector(int colInd, const Vector<Type>& vec)
		{
			if (colInd < 0 || colInd >= _cols)
				throw MatrixAccessBoundsError("InitColWithVector - invalid column index", 0, colInd, _rows, _cols);
			if (vec.size() != _rows)
				throw MatrixDimensionError("InitColWithVector - vector size must match rows", _rows, _cols, vec.size(), -1);
			
			for (int i = 0; i < _rows; ++i)
				_data[idx(i, colInd)] = vec[i];
		}
		
		void SwapRows(int k, int l)
		{
			if (k < 0 || k >= _rows || l < 0 || l >= _rows)
				throw MatrixDimensionError("Matrix::SwapRows - invalid row index", _rows, _cols, k, l);
			if (k == l) return;
			
			for (int j = 0; j < _cols; ++j)
				std::swap(_data[idx(k, j)], _data[idx(l, j)]);
		}
		
		void SwapCols(int k, int l)
		{
			if (k < 0 || k >= _cols || l < 0 || l >= _cols)
				throw MatrixDimensionError("Matrix::SwapCols - invalid column index", _rows, _cols, k, l);
			if (k == l) return;
			
			for (int i = 0; i < _rows; ++i)
				std::swap(_data[idx(i, k)], _data[idx(i, l)]);
		}

		///////////////////////              Creating views                //////////////////////
		MatrixViewNew<Type> block(int startRow, int startCol, int numRows, int numCols)
		{
			if (startRow < 0 || startCol < 0 || numRows <= 0 || numCols <= 0 ||
			    startRow + numRows > _rows || startCol + numCols > _cols)
				throw MatrixDimensionError("Matrix::block - invalid parameters", _rows, _cols, startRow, startCol);
			return MatrixViewNew<Type>(&_data[idx(startRow, startCol)], numRows, numCols, _cols);
		}
		
		const MatrixViewNew<Type> block(int startRow, int startCol, int numRows, int numCols) const
		{
			if (startRow < 0 || startCol < 0 || numRows <= 0 || numCols <= 0 ||
			    startRow + numRows > _rows || startCol + numCols > _cols)
				throw MatrixDimensionError("Matrix::block - invalid parameters", _rows, _cols, startRow, startCol);
			return MatrixViewNew<Type>(const_cast<Type*>(&_data[idx(startRow, startCol)]), numRows, numCols, _cols);
		}

		///////////////////////               Matrix properties            //////////////////////
		bool IsUnit(double eps = Defaults::IsMatrixUnitTolerance) const
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
		
		bool IsDiagonal(double eps = Defaults::IsMatrixDiagonalTolerance) const
		{
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					if (i != j && Abs((*this)(i, j)) > eps)
						return false;
			return true;
		}
		
		bool IsDiagDominant() const
		{
			for (int i = 0; i < _rows; ++i) {
				Type sum{0};
				for (int j = 0; j < _cols; ++j)
					if (i != j)
						sum += Abs((*this)(i, j));
				if (Abs((*this)(i, i)) < sum)
					return false;
			}
			return true;
		}
		
		bool IsSymmetric() const
		{
			if (_rows != _cols) return false;
			
			for (int i = 0; i < _rows; ++i)
				for (int j = i + 1; j < _cols; ++j)
					if ((*this)(i, j) != (*this)(j, i))
						return false;
			return true;
		}
		
		bool IsAntiSymmetric() const
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

		///////////////////////             Matrix norm calculations       //////////////////////
		Real NormL1() const
		{
			Real norm{0};
			for (const auto& elem : _data)
				norm += Abs(elem);
			return norm;
		}
		
		Real NormL2() const
		{
			Real norm{0};
			for (const auto& elem : _data)
				norm += elem * elem;
			return std::sqrt(norm);
		}
		
		Real NormLInf() const
		{
			Real norm{0};
			for (const auto& elem : _data)
				norm = std::max(norm, static_cast<Real>(Abs(elem)));
			return norm;
		}

		///////////////////////               Access operators             //////////////////////
		
		// Primary access: operator() - checked bounds in debug, unchecked in release
		inline Type  operator()(int i, int j) const { return _data[idx(i, j)]; }
		inline Type& operator()(int i, int j)       { return _data[idx(i, j)]; }
		
		// Legacy compatibility: operator[] returns pointer to row start
		// This maintains compatibility with code using mat[i][j] syntax
		inline Type* operator[](int i)             { return &_data[idx(i, 0)]; }
		inline const Type* operator[](int i) const { return &_data[idx(i, 0)]; }
		
		// Checked access
		Type at(int i, int j) const
		{
			if (i < 0 || i >= _rows || j < 0 || j >= _cols)
				throw MatrixAccessBoundsError("Matrix::at", i, j, _rows, _cols);
			return _data[idx(i, j)];
		}
		
		Type& at(int i, int j)
		{
			if (i < 0 || i >= _rows || j < 0 || j >= _cols)
				throw MatrixAccessBoundsError("Matrix::at", i, j, _rows, _cols);
			return _data[idx(i, j)];
		}

		///////////////////////             Iterator support               //////////////////////
		
		// Flat iterators - traverse all elements in row-major order
		auto begin()        { return _data.begin(); }
		auto end()          { return _data.end(); }
		auto begin() const  { return _data.begin(); }
		auto end() const    { return _data.end(); }
		auto cbegin() const { return _data.cbegin(); }
		auto cend() const   { return _data.cend(); }

		// Row iterator proxy
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
		
		row_range row(int i)
		{
			if (i < 0 || i >= _rows)
				throw MatrixAccessBoundsError("Matrix::row", i, 0, _rows, _cols);
			return row_range(&_data[idx(i, 0)], _cols);
		}
		
		const_row_range row(int i) const
		{
			if (i < 0 || i >= _rows)
				throw MatrixAccessBoundsError("Matrix::row", i, 0, _rows, _cols);
			return const_row_range(&_data[idx(i, 0)], _cols);
		}

		// Column iterator (strided)
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
		
		col_range col(int j)
		{
			if (j < 0 || j >= _cols)
				throw MatrixAccessBoundsError("Matrix::col", 0, j, _rows, _cols);
			return col_range(&_data[j], _rows, _cols);
		}

		///////////////////////              Equality operations           //////////////////////
		bool operator==(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				return false;
			return _data == b._data;
		}
		
		bool operator!=(const Matrix& b) const { return !(*this == b); }
		
		bool IsEqualTo(const Matrix& b, Type eps = Defaults::MatrixIsEqualTolerance) const
		{
			if (_rows != b._rows || _cols != b._cols)
				return false;
			
			for (size_t i = 0; i < _data.size(); ++i)
				if (Abs(_data[i] - b._data[i]) > eps)
					return false;
			return true;
		}
		
		static bool AreEqual(const Matrix& a, const Matrix& b, Type eps = Defaults::MatrixIsEqualTolerance)
		{
			return a.IsEqualTo(b, eps);
		}

		///////////////////////              Arithmetic operators          //////////////////////
		
		// Unary minus
		Matrix operator-() const
		{
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = -_data[i];
			return ret;
		}
		
		// Matrix + Matrix
		Matrix operator+(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+ - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] + b._data[i];
			return ret;
		}
		
		Matrix& operator+=(const Matrix& b)
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+= - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] += b._data[i];
			return *this;
		}
		
		// Matrix - Matrix
		Matrix operator-(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator- - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] - b._data[i];
			return ret;
		}
		
		Matrix& operator-=(const Matrix& b)
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator-= - dimensions must match", _rows, _cols, b._rows, b._cols);
			
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] -= b._data[i];
			return *this;
		}
		
		// Matrix * Matrix
		Matrix operator*(const Matrix& b) const
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
		
		// Matrix * scalar
		Matrix operator*(const Type& scalar) const
		{
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] * scalar;
			return ret;
		}
		
		Matrix& operator*=(const Type& scalar)
		{
			for (auto& elem : _data)
				elem *= scalar;
			return *this;
		}
		
		// Matrix / scalar
		Matrix operator/(const Type& scalar) const
		{
			Matrix ret(_rows, _cols);
			for (size_t i = 0; i < _data.size(); ++i)
				ret._data[i] = _data[i] / scalar;
			return ret;
		}
		
		Matrix& operator/=(const Type& scalar)
		{
			for (auto& elem : _data)
				elem /= scalar;
			return *this;
		}
		
		// Matrix * Vector
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
		
		// Scalar * Matrix (friend)
		friend Matrix operator*(const Type& scalar, const Matrix& m)
		{
			return m * scalar;
		}
		
		// Vector * Matrix (friend)
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
		Type Trace() const
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::Trace - must be square", _rows, _cols, -1, -1);
			
			Type sum{0};
			for (int i = 0; i < _rows; ++i)
				sum += (*this)(i, i);
			return sum;
		}
		
		void Invert()
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::Invert - must be square", _rows, _cols, -1, -1);
			
			// Gauss-Jordan elimination with pivoting
			int n = _rows;
			Matrix& a = *this;
			std::vector<int> indxc(n), indxr(n), ipiv(n, 0);
			
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
				
				if (a(icol, icol) == Type{0})
					throw SingularMatrixError("Matrix::Invert - Singular Matrix");
				
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
		
		Matrix GetInverse() const
		{
			Matrix ret(*this);
			ret.Invert();
			return ret;
		}
		
		void Transpose()
		{
			if (_rows != _cols)
				throw MatrixDimensionError("Matrix::Transpose - in-place requires square matrix", _rows, _cols, -1, -1);
			
			for (int i = 0; i < _rows; ++i)
				for (int j = i + 1; j < _cols; ++j)
					std::swap(_data[idx(i, j)], _data[idx(j, i)]);
		}
		
		Matrix GetTranspose() const
		{
			Matrix ret(_cols, _rows);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					ret(j, i) = (*this)(i, j);
			return ret;
		}

		///////////////////////                    I/O                    //////////////////////
		void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const
		{
			if (IsEmpty()) {
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
		
		void Print(std::ostream& stream, int width, int precision) const
		{
			MatrixPrintFormat fmt = MatrixPrintFormat::Default();
			fmt.width = width;
			fmt.precision = precision;
			fmt.fixed = false;
			Print(stream, fmt);
		}
		
		// Print with zero threshold - values below threshold printed as 0
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
		
		friend std::ostream& operator<<(std::ostream& stream, const Matrix& a)
		{
			a.Print(stream, 10, 3);
			return stream;
		}
		
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;
			Print(str, width, precision);
			return str.str();
		}

		///////////////////////              File I/O                      //////////////////////
		static bool LoadFromFile(const std::string& filename, Matrix& outMat)
		{
			std::ifstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not open file " << filename << " for reading." << std::endl;
				return false;
			}
			
			int rows, cols;
			file >> rows >> cols;
			outMat.Resize(rows, cols);
			
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					file >> outMat(i, j);
			
			file.close();
			return true;
		}
		
		static bool SaveToFile(const Matrix& mat, const std::string& filename)
		{
			std::ofstream file(filename);
			if (!file.is_open()) {
				std::cerr << "Error: could not create file " << filename << " for writing." << std::endl;
				return false;
			}
			
			file << mat._rows << " " << mat._cols << std::endl;
			for (int i = 0; i < mat._rows; ++i) {
				for (int j = 0; j < mat._cols; ++j)
					file << mat(i, j) << " ";
				file << std::endl;
			}
			
			file.close();
			return true;
		}
		
		static bool LoadFromCSV(const std::string& filename, Matrix& outMat)
		{
			std::ifstream file(filename);
			if (!file.is_open()) return false;
			
			std::vector<std::vector<Type>> data;
			std::string line;
			
			while (std::getline(file, line)) {
				std::vector<Type> row;
				std::stringstream ss(line);
				std::string cell;
				while (std::getline(ss, cell, ',')) {
					std::stringstream cellStream(cell);
					Type value;
					cellStream >> value;
					row.push_back(value);
				}
				if (!row.empty())
					data.push_back(row);
			}
			
			file.close();
			if (data.empty()) return false;
			
			outMat = Matrix(data);
			return true;
		}
		
		static bool SaveToCSV(const Matrix& mat, const std::string& filename)
		{
			std::ofstream file(filename);
			if (!file.is_open()) return false;
			
			for (int i = 0; i < mat._rows; ++i) {
				for (int j = 0; j < mat._cols; ++j) {
					file << mat(i, j);
					if (j < mat._cols - 1) file << ",";
				}
				file << "\n";
			}
			
			file.close();
			return true;
		}
	};

	//////////////////////               Default Matrix typedefs       //////////////////////
	typedef Matrix<int>     MatrixInt;
	typedef Matrix<float>   MatrixFlt;
	typedef Matrix<double>  MatrixDbl;
	typedef Matrix<Complex> MatrixComplex;

	// Legacy typedef aliases (for backward compatibility)
	typedef Matrix<int>     MatI;
	typedef Matrix<float>   MatF;
	typedef Matrix<double>  MatD;
	typedef Matrix<Complex> MatC;

	// Verify noexcept move operations
	static_assert(std::is_nothrow_move_constructible_v<Matrix<double>>,
	              "Matrix<double> should be nothrow move constructible");
	static_assert(std::is_nothrow_move_assignable_v<Matrix<double>>,
	              "Matrix<double> should be nothrow move assignable");
}

#endif // MML_MATRIX_H
