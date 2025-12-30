///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixSym.h                                                         ///
///  Description: Symmetric matrix class with optimized storage and operations        ///
///               Stores only upper/lower triangular portion                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_SYM_H
#define MML_MATRIX_SYM_H

#include "MMLBase.h"

#include "Matrix.h"

#include <vector>
#include <algorithm>
#include <stdexcept>

namespace MML
{
    // Constants for symmetric matrix allocation limits
    struct MatrixSymLimits {
        static constexpr int MAX_DIMENSION = 10000;                 // Maximum dimension
        static constexpr size_t MAX_ELEMENTS = 50'000'000;          // Maximum storage elements
    };

    template<class Type>
    class MatrixSym
    {
    private:
        int _dim;
        std::vector<Type> _data;  // Contiguous storage for lower triangular elements

        //////////////////       Helper functions for linear indexing       //////////////////
        // Storage layout: row-wise lower triangular
        // [a00, a10, a11, a20, a21, a22, a30, a31, a32, a33, ...]
        // For element (i,j) where i >= j: index = i*(i+1)/2 + j
        
        static constexpr size_t numElements(int dim) noexcept {
            return static_cast<size_t>(dim) * (dim + 1) / 2;
        }
        
        static constexpr size_t linearIndex(int i, int j) noexcept {
            // Ensure i >= j for lower triangular access
            if (i < j) std::swap(i, j);
            return static_cast<size_t>(i) * (i + 1) / 2 + j;
        }

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
        typedef Type value_type;      // make Type available externally

        ///////////////////////          Constructors and destructor       //////////////////////
        
        // Default constructor - creates empty (0x0) matrix
        MatrixSym() noexcept : _dim(0), _data() {}
        
        // Dimension constructor - creates dim x dim matrix initialized to zero
        explicit MatrixSym(int dim) : _dim(0), _data() 
        {
            if (dim == 0) return;
            
            validateDimension(dim);
            _dim = dim;
            _data.resize(numElements(dim), Type{0});
        }
        
        // Dimension + value constructor - creates dim x dim matrix initialized to val
        MatrixSym(int dim, Type val) : _dim(0), _data()
        {
            if (dim == 0) return;
            
            validateDimension(dim);
            _dim = dim;
            _data.resize(numElements(dim), val);
        }
        
        // Initializer list constructor - takes lower triangular elements row-wise
        // Example for 3x3: {a00, a10, a11, a20, a21, a22}
        MatrixSym(int dim, std::initializer_list<Type> values) : _dim(0), _data()
        {
            if (dim == 0) {
                if (values.size() != 0)
                    throw MatrixDimensionError("MatrixSym: non-empty initializer for zero dimension", 0, -1, -1, -1);
                return;
            }
            
            validateDimension(dim);
            
            size_t expected = numElements(dim);
            if (values.size() != expected)
                throw MatrixDimensionError("MatrixSym: initializer list size mismatch - expected " + 
                    std::to_string(expected) + " but got " + std::to_string(values.size()), 
                    dim, -1, static_cast<int>(values.size()), -1);
            
            _dim = dim;
            _data.assign(values.begin(), values.end());
        }
        
        // Copy constructor
        MatrixSym(const MatrixSym& other) : _dim(other._dim), _data(other._data) {}
        
        // Move constructor
        MatrixSym(MatrixSym&& other) noexcept 
            : _dim(other._dim), _data(std::move(other._data))
        {
            other._dim = 0;
        }
        
        // Destructor - default is fine with std::vector
        ~MatrixSym() = default;

        ////////////////////////               Standard stuff             ////////////////////////
        
        int Dim() const noexcept { return _dim; }
        int RowNum() const noexcept { return _dim; }
        int ColNum() const noexcept { return _dim; }
        
        // Number of stored elements (lower triangular)
        size_t size() const noexcept { return _data.size(); }
        
        // Check if matrix is empty (0x0)
        bool empty() const noexcept { return _dim == 0; }
        
        // Raw data access (for interop with libraries)
        Type* data() noexcept { return _data.data(); }
        const Type* data() const noexcept { return _data.data(); }
        
        // Resize the matrix (destroys existing data)
        void Resize(int newDim)
        {
            if (newDim == _dim) return;
            
            if (newDim == 0) {
                _dim = 0;
                _data.clear();
                return;
            }
            
            validateDimension(newDim);
            _dim = newDim;
            _data.assign(numElements(newDim), Type{0});
        }
        
        // Resize with preservation of existing elements where possible
        void Resize(int newDim, bool preserveData)
        {
            if (!preserveData || _dim == 0) {
                Resize(newDim);
                return;
            }
            
            if (newDim == _dim) return;
            
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
        
        // Set all elements to zero
        void SetToZero() noexcept
        {
            std::fill(_data.begin(), _data.end(), Type{0});
        }
        
        // Set all elements to a value
        void SetToValue(Type val) noexcept
        {
            std::fill(_data.begin(), _data.end(), val);
        }

        ///////////////////////          Comparison                   //////////////////////
        
        bool IsEqual(const MatrixSym& b, Type eps = Defaults::MatrixIsEqualTolerance) const
        {
            if (_dim != b._dim)
                return false;

            for (size_t i = 0; i < _data.size(); ++i) {
                if (Abs(_data[i] - b._data[i]) > eps)
                    return false;
            }
            return true;
        }
        
        static bool AreEqual(const MatrixSym& a, const MatrixSym& b, Type eps = Defaults::MatrixIsEqualTolerance)
        {
            return a.IsEqual(b, eps);
        }

        ///////////////////////          Conversion to Matrix         //////////////////////
        
        Matrix<Type> GetAsMatrix() const
        {
            if (_dim == 0) return Matrix<Type>();
            
            Matrix<Type> ret(_dim, _dim);
            for (int i = 0; i < _dim; ++i) {
                for (int j = 0; j < _dim; ++j) {
                    ret(i, j) = (*this)(i, j);
                }
            }
            return ret;
        }

        /////////////////////          Init from regular Matrix           /////////////////////
        
        // Create symmetric matrix from lower triangular part of a matrix
        static MatrixSym FromLower(const Matrix<Type>& m)
        {
            if (m.RowNum() != m.ColNum())
                throw MatrixDimensionError("MatrixSym::FromLower - must be square matrix", 
                    m.RowNum(), m.ColNum(), -1, -1);

            int dim = m.RowNum();
            if (dim == 0) return MatrixSym();
            
            MatrixSym ret(dim);
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j <= i; ++j) {
                    ret(i, j) = m(i, j);
                }
            }
            return ret;
        }
        
        // Create symmetric matrix from upper triangular part of a matrix
        static MatrixSym FromUpper(const Matrix<Type>& m)
        {
            if (m.RowNum() != m.ColNum())
                throw MatrixDimensionError("MatrixSym::FromUpper - must be square matrix", 
                    m.RowNum(), m.ColNum(), -1, -1);

            int dim = m.RowNum();
            if (dim == 0) return MatrixSym();
            
            MatrixSym ret(dim);
            for (int i = 0; i < dim; ++i) {
                for (int j = i; j < dim; ++j) {
                    // a(i,j) from upper = a(j,i) in storage = a(i,j) symmetric
                    ret(j, i) = m(i, j);  // Transpose to lower
                }
            }
            return ret;
        }
        
        // Create symmetric matrix as (M + M^T) / 2
        static MatrixSym FromFullMatrix(const Matrix<Type>& m)
        {
            if (m.RowNum() != m.ColNum())
                throw MatrixDimensionError("MatrixSym::FromFullMatrix - must be square matrix", 
                    m.RowNum(), m.ColNum(), -1, -1);

            int dim = m.RowNum();
            if (dim == 0) return MatrixSym();
            
            MatrixSym ret(dim);
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j <= i; ++j) {
                    ret(i, j) = (m(i, j) + m(j, i)) / Type{2};
                }
            }
            return ret;
        }
        
        // Legacy methods (deprecated - use static factory methods instead)
        [[deprecated("Use static MatrixSym::FromLower() instead")]]
        MatrixSym InitFromLower(const Matrix<Type>& b)
        {
            return FromLower(b);
        }
        
        [[deprecated("Use static MatrixSym::FromUpper() instead")]]
        MatrixSym InitFromUpper(const Matrix<Type>& b)
        {
            return FromUpper(b);
        }

        ///////////////////////          Static Factory Methods       //////////////////////
        
        // Create identity matrix of given dimension
        static MatrixSym GetUnitMatrix(int dim)
        {
            MatrixSym ret(dim);  // Initializes to zero
            for (int i = 0; i < dim; ++i) {
                ret(i, i) = Type{1};
            }
            return ret;
        }
        
        // Create diagonal matrix from vector
        static MatrixSym GetDiagonalMatrix(const Vector<Type>& diag)
        {
            int dim = static_cast<int>(diag.size());
            MatrixSym ret(dim);  // Initializes to zero
            for (int i = 0; i < dim; ++i) {
                ret(i, i) = diag[i];
            }
            return ret;
        }

        /////////////////////          Vector-Matrix conversion           /////////////////////
        
        Vector<Type> VectorFromRow(int rowInd) const
        {
            if (rowInd < 0 || rowInd >= _dim)
                throw MatrixAccessBoundsError("VectorFromRow - row index out of bounds", 
                    rowInd, 0, _dim, _dim);

            Vector<Type> ret(_dim);
            for (int j = 0; j < _dim; ++j) {
                ret[j] = (*this)(rowInd, j);
            }
            return ret;
        }
        
        Vector<Type> VectorFromColumn(int colInd) const
        {
            if (colInd < 0 || colInd >= _dim)
                throw MatrixAccessBoundsError("VectorFromColumn - column index out of bounds", 
                    0, colInd, _dim, _dim);

            Vector<Type> ret(_dim);
            for (int i = 0; i < _dim; ++i) {
                ret[i] = (*this)(i, colInd);
            }
            return ret;
        }
        
        Vector<Type> VectorFromDiagonal() const
        {
            Vector<Type> ret(_dim);
            for (int i = 0; i < _dim; ++i) {
                ret[i] = (*this)(i, i);
            }
            return ret;
        }

        ///////////////////////          Matrix Properties            //////////////////////
        
        // Compute trace (sum of diagonal elements)
        Type Trace() const noexcept
        {
            Type sum{0};
            for (int i = 0; i < _dim; ++i) {
                sum += (*this)(i, i);
            }
            return sum;
        }
        
        // Frobenius norm: sqrt(sum of squared elements)
        // Note: off-diagonal elements counted twice in full matrix
        Type NormFrobenius() const
        {
            Type sum{0};
            for (int i = 0; i < _dim; ++i) {
                // Diagonal element
                sum += (*this)(i, i) * (*this)(i, i);
                // Off-diagonal elements (counted twice in full matrix)
                for (int j = 0; j < i; ++j) {
                    sum += Type{2} * (*this)(i, j) * (*this)(i, j);
                }
            }
            return std::sqrt(sum);
        }
        
        // Infinity norm (max row sum of absolute values)
        Type NormInf() const
        {
            Type maxSum{0};
            for (int i = 0; i < _dim; ++i) {
                Type rowSum{0};
                for (int j = 0; j < _dim; ++j) {
                    rowSum += Abs((*this)(i, j));
                }
                if (rowSum > maxSum) maxSum = rowSum;
            }
            return maxSum;
        }
        
        // One norm (max column sum of absolute values)
        // For symmetric matrix, same as infinity norm
        Type Norm1() const
        {
            return NormInf();  // Symmetric property: row sum = column sum
        }

        ///////////////////////////            Operators             ///////////////////////////
        
        // Copy assignment
        MatrixSym& operator=(const MatrixSym& other)
        {
            if (this != &other) {
                _dim = other._dim;
                _data = other._data;
            }
            return *this;
        }
        
        // Move assignment
        MatrixSym& operator=(MatrixSym&& other) noexcept
        {
            if (this != &other) {
                _dim = other._dim;
                _data = std::move(other._data);
                other._dim = 0;
            }
            return *this;
        }

        ////////////////////////           Access operators             ///////////////////////
        
        // Element access without bounds checking - primary access method
        Type operator()(int i, int j) const {
            return _data[linearIndex(i, j)];
        }
        
        Type& operator()(int i, int j) {
            return _data[linearIndex(i, j)];
        }

        // Element access with bounds checking (like std::vector::at)
        Type at(int i, int j) const
        {
            if (i < 0 || i >= _dim || j < 0 || j >= _dim)
                throw MatrixAccessBoundsError("MatrixSym::at", i, j, _dim, _dim);
            return _data[linearIndex(i, j)];
        }
        
        Type& at(int i, int j)
        {
            if (i < 0 || i >= _dim || j < 0 || j >= _dim)
                throw MatrixAccessBoundsError("MatrixSym::at", i, j, _dim, _dim);
            return _data[linearIndex(i, j)];
        }

        // Legacy ElemAt (alias for at)
        Type ElemAt(int i, int j) const { return at(i, j); }
        Type& ElemAt(int i, int j) { return at(i, j); }

        ///////////////////////           Arithmetic operators             //////////////////////
        
        MatrixSym operator+(const MatrixSym& b) const
        {
            if (_dim != b._dim)
                throw MatrixDimensionError("MatrixSym::operator+ - dimensions must match", 
                    _dim, -1, b._dim, -1);

            MatrixSym ret(_dim);
            for (size_t i = 0; i < _data.size(); ++i) {
                ret._data[i] = _data[i] + b._data[i];
            }
            return ret;
        }
        
        MatrixSym& operator+=(const MatrixSym& b)
        {
            if (_dim != b._dim)
                throw MatrixDimensionError("MatrixSym::operator+= - dimensions must match", 
                    _dim, -1, b._dim, -1);

            for (size_t i = 0; i < _data.size(); ++i) {
                _data[i] += b._data[i];
            }
            return *this;
        }
        
        MatrixSym operator-() const
        {
            MatrixSym ret(_dim);
            for (size_t i = 0; i < _data.size(); ++i) {
                ret._data[i] = -_data[i];
            }
            return ret;
        }
        
        MatrixSym operator-(const MatrixSym& b) const
        {
            if (_dim != b._dim)
                throw MatrixDimensionError("MatrixSym::operator- - dimensions must match", 
                    _dim, -1, b._dim, -1);

            MatrixSym ret(_dim);
            for (size_t i = 0; i < _data.size(); ++i) {
                ret._data[i] = _data[i] - b._data[i];
            }
            return ret;
        }
        
        MatrixSym& operator-=(const MatrixSym& b)
        {
            if (_dim != b._dim)
                throw MatrixDimensionError("MatrixSym::operator-= - dimensions must match", 
                    _dim, -1, b._dim, -1);

            for (size_t i = 0; i < _data.size(); ++i) {
                _data[i] -= b._data[i];
            }
            return *this;
        }

        // Matrix-Matrix multiplication (result is NOT symmetric in general)
        Matrix<Type> operator*(const MatrixSym& b) const
        {
            if (_dim != b._dim)
                throw MatrixDimensionError("MatrixSym::operator*(MatrixSym) - dimensions must match", 
                    _dim, _dim, b._dim, b._dim);

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

        Matrix<Type> operator*(const Matrix<Type>& b) const
        {
            if (_dim != b.RowNum())
                throw MatrixDimensionError("MatrixSym::operator*(Matrix) - dimension mismatch", 
                    _dim, _dim, b.RowNum(), b.ColNum());

            Matrix<Type> ret(_dim, b.ColNum());
            for (int i = 0; i < _dim; ++i) {
                for (int j = 0; j < b.ColNum(); ++j) {
                    Type sum{0};
                    for (int k = 0; k < _dim; ++k) {
                        sum += (*this)(i, k) * b(k, j);
                    }
                    ret(i, j) = sum;
                }
            }
            return ret;
        }

        // Scalar multiplication
        friend MatrixSym operator*(const MatrixSym& a, Type scalar)
        {
            MatrixSym ret(a._dim);
            for (size_t i = 0; i < a._data.size(); ++i) {
                ret._data[i] = a._data[i] * scalar;
            }
            return ret;
        }
        
        friend MatrixSym operator*(Type scalar, const MatrixSym& a)
        {
            return a * scalar;
        }
        
        MatrixSym& operator*=(Type scalar)
        {
            for (auto& elem : _data) {
                elem *= scalar;
            }
            return *this;
        }
        
        // Scalar division
        friend MatrixSym operator/(const MatrixSym& a, Type scalar)
        {
            MatrixSym ret(a._dim);
            for (size_t i = 0; i < a._data.size(); ++i) {
                ret._data[i] = a._data[i] / scalar;
            }
            return ret;
        }
        
        MatrixSym& operator/=(Type scalar)
        {
            for (auto& elem : _data) {
                elem /= scalar;
            }
            return *this;
        }

        // Matrix-Vector multiplication
        friend Vector<Type> operator*(const MatrixSym& a, const Vector<Type>& v)
        {
            if (a._dim != static_cast<int>(v.size()))
                throw MatrixDimensionError("operator*(MatrixSym, Vector) - dimension mismatch", 
                    a._dim, a._dim, static_cast<int>(v.size()), -1);

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
        
        friend Vector<Type> operator*(const Vector<Type>& v, const MatrixSym& a)
        {
            if (static_cast<int>(v.size()) != a._dim)
                throw MatrixDimensionError("operator*(Vector, MatrixSym) - dimension mismatch", 
                    static_cast<int>(v.size()), -1, a._dim, a._dim);

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

        ////////////////////////            Inverse              ///////////////////////
        
        Matrix<Type> GetInverse() const
        {
            if (_dim == 0)
                throw MatrixDimensionError("MatrixSym::GetInverse - cannot invert empty matrix", 0, 0, -1, -1);

            Matrix<Type> a = GetAsMatrix();
            a.Invert();
            return a;
        }

        ///////////////////////////               I/O                 ///////////////////////////
        
        std::string to_string(int width = 10, int precision = 3) const
        {
            std::stringstream str;
            Print(str, width, precision);
            return str.str();
        }
        
        // Formatted print with MatrixPrintFormat
        void Print(std::ostream& stream, const MatrixPrintFormat& fmt = MatrixPrintFormat::Default()) const
        {
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
                    if (i > 0) stream << "; ";
                    for (int j = 0; j < _dim; ++j) {
                        if (j > 0) stream << fmt.delimiter;
                        stream << std::setw(fmt.width) << std::setprecision(fmt.precision) << (*this)(i, j);
                    }
                }
                if (fmt.showBrackets)
                    stream << "]";
                stream << std::endl;
            }
            else {
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

        // Legacy print method
        void Print(std::ostream& stream, int width, int precision) const
        {
            MatrixPrintFormat fmt = MatrixPrintFormat::Default();
            fmt.width = width;
            fmt.precision = precision;
            Print(stream, fmt);
        }
        
        friend std::ostream& operator<<(std::ostream& stream, const MatrixSym& a)
        {
            a.Print(stream, 10, 3);
            return stream;
        }
    };
}
#endif
