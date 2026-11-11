///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SparseMatrixCSR.h                                                   ///
///  Description: Compressed Sparse Row (CSR) format - efficient for SpMV and solvers ///
///               Primary format for iterative linear solvers                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_BASE_SPARSE_MATRIX_CSR_H
#define MML_BASE_SPARSE_MATRIX_CSR_H

#include <mml/MMLExceptions.h>

#include "SparseMatrixCOO.h"
#include <mml/base/Vector/Vector.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <complex>
#include <numeric>
#include <utility>
#include <iostream>
#include <iomanip>
#include <limits>
#include <type_traits>

namespace MML::SparseMatrix
{
    ///////////////////////////////////////////////////////////////////////////
    ///                    SparseMatrixCSR Class                            ///
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief Compressed Sparse Row (CSR) format sparse matrix
    /// 
    /// Storage:
    /// - values_: non-zero values in row-major order
    /// - col_indices_: column index for each value
    /// - row_pointers_: index into values_ where each row starts
    /// 
    /// For matrix A with m rows:
    /// - row_pointers_ has size m+1
    /// - Row i has entries from row_pointers_[i] to row_pointers_[i+1]-1
    /// - nnz = row_pointers_[m]
    /// 
    /// Ideal for:
    /// - Matrix-vector products (SpMV)
    /// - Row-wise access
    /// - Iterative solvers (CG, BiCGSTAB, GMRES)
    /// 
    /// @tparam T Element type (typically double or float)
    template<typename T>
    class SparseMatrixCSR {
    private:
        using Magnitude = decltype(std::abs(T{}));
        using SymmetryTolerance = std::conditional_t<std::is_floating_point_v<Magnitude>, Magnitude, Real>;
        int rows_;
        int cols_;
        std::vector<T> values_;
        std::vector<int> col_indices_;
        std::vector<int> row_pointers_;
        
    public:
        // ==================== Constructors ====================
        
        /// @brief Default constructor - creates empty 0x0 matrix
        SparseMatrixCSR() : rows_(0), cols_(0), row_pointers_(1, 0) {}
        
        /// @brief Construct from dimensions (empty matrix)
        SparseMatrixCSR(int rows, int cols) 
            : rows_(rows), cols_(cols) {
            if (rows < 0 || cols < 0) {
                throw MatrixDimensionError("SparseMatrixCSR: matrix dimensions must be non-negative", rows, cols, -1, -1);
            }
            row_pointers_.assign(static_cast<size_t>(rows) + 1, 0);
        }
        
        /// @brief Construct from COO format (main construction path)
        explicit SparseMatrixCSR(SparseMatrixCOO<T>& coo) {
            FromCOO(coo);
        }
        
        /// @brief Construct from raw CSR data (advanced)
        SparseMatrixCSR(int rows, int cols,
                       std::vector<T> values,
                       std::vector<int> col_indices,
                       std::vector<int> row_pointers)
            : rows_(rows), cols_(cols),
              values_(std::move(values)),
              col_indices_(std::move(col_indices)),
              row_pointers_(std::move(row_pointers)) {
            if (!validate()) throw ArgumentError("SparseMatrixCSR: invalid CSR storage");
        }
        
        // ==================== Conversion from COO ====================
        
        /// @brief Build CSR from COO matrix
        void FromCOO(SparseMatrixCOO<T>& coo) {
            // Consolidate COO (sort and sum duplicates)
            coo.consolidate();
            const auto& triplets = coo.triplets();
            if (triplets.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
                throw MatrixAllocationError("SparseMatrixCSR: nonzero count exceeds index range", coo.rows(), coo.cols());
            int nnz = static_cast<int>(triplets.size());
            std::vector<T> values(nnz);
            std::vector<int> indices(nnz);
            std::vector<int> pointers(static_cast<size_t>(coo.rows()) + 1, 0);
            
            // Count entries per row
            for (const auto& t : triplets) {
                pointers[t.row + 1]++;
            }
            
            // Cumulative sum to get row pointers
            for (int i = 0; i < coo.rows(); ++i) {
                pointers[i + 1] += pointers[i];
            }
            
            // Fill values and column indices
            std::vector<int> row_counts(coo.rows(), 0);
            for (const auto& t : triplets) {
                int idx = pointers[t.row] + row_counts[t.row];
                values[idx] = t.value;
                indices[idx] = t.col;
                row_counts[t.row]++;
            }
            values_.swap(values);
            col_indices_.swap(indices);
            row_pointers_.swap(pointers);
            rows_ = coo.rows();
            cols_ = coo.cols();
        }
        
        // ==================== Dimensions ====================
        
        int rows() const { return rows_; }
        int cols() const { return cols_; }
        int nnz() const { return static_cast<int>(values_.size()); }
        
        // ==================== Raw Data Access ====================
        
        const std::vector<T>& values() const { return values_; }
        const std::vector<int>& colIndices() const { return col_indices_; }
        const std::vector<int>& rowPointers() const { return row_pointers_; }
        
        std::vector<T>& values() { return values_; }
        std::vector<int>& colIndices() { return col_indices_; }
        std::vector<int>& rowPointers() { return row_pointers_; }

        /// @brief Check dimensions and canonical CSR storage invariants
        bool validate() const {
            if (rows_ < 0 || cols_ < 0 || values_.size() != col_indices_.size()) return false;
            if (values_.size() > static_cast<size_t>(std::numeric_limits<int>::max())) return false;
            if (row_pointers_.size() != static_cast<size_t>(rows_) + 1 || row_pointers_.front() != 0) return false;
            if (row_pointers_.back() != static_cast<int>(values_.size())) return false;
            for (int row = 0; row < rows_; ++row) {
                if (row_pointers_[row] < 0 || row_pointers_[row] > row_pointers_[row + 1] ||
                    row_pointers_[row + 1] > nnz()) return false;
            }
            for (int row = 0; row < rows_; ++row) {
                int previousColumn = -1;
                for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                    if (col_indices_[k] < 0 || col_indices_[k] >= cols_ || col_indices_[k] <= previousColumn) return false;
                    previousColumn = col_indices_[k];
                }
            }
            return true;
        }
        
        // ==================== Element Access ====================
        
        /// @brief Get element at (row, col) - O(nnz_row) lookup
        T operator()(int row, int col) const {
            if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                throw MatrixAccessBoundsError("SparseMatrixCSR::operator() - index out of bounds", row, col, rows_, cols_);
            }
            
            for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                if (col_indices_[k] == col) {
                    return values_[k];
                }
                if (col_indices_[k] > col) {
                    break;  // Columns are sorted within row
                }
            }
            return T(0);
        }
        
        /// @brief Get pointer to element (nullptr if not present)
        T* findElement(int row, int col) {
            for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                if (col_indices_[k] == col) {
                    return &values_[k];
                }
            }
            return nullptr;
        }
        
        const T* findElement(int row, int col) const {
            for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                if (col_indices_[k] == col) {
                    return &values_[k];
                }
            }
            return nullptr;
        }
        
        // ==================== Row Operations ====================
        
        /// @brief Number of non-zeros in row i
        int rowNnz(int row) const {
            return row_pointers_[row + 1] - row_pointers_[row];
        }
        
        /// @brief Get column indices for row i
        std::pair<const int*, const int*> rowColIndices(int row) const {
            return {col_indices_.data() + row_pointers_[row],
                    col_indices_.data() + row_pointers_[row + 1]};
        }
        
        /// @brief Get values for row i
        std::pair<const T*, const T*> rowValues(int row) const {
            return {values_.data() + row_pointers_[row],
                    values_.data() + row_pointers_[row + 1]};
        }
        
        // ==================== Matrix-Vector Operations ====================
        
        /// @brief Sparse matrix-vector product: y = A * x
        void Multiply(const std::vector<T>& x, std::vector<T>& y) const {
            if (x.size() != static_cast<size_t>(cols_)) {
                throw VectorDimensionError("SparseMatrixCSR::multiply - vector size must match matrix cols", static_cast<int>(x.size()), cols_);
            }
            if (&x == &y) {
                const auto input = x;
                Multiply(input, y);
                return;
            }
            y.resize(rows_);
            
            for (int i = 0; i < rows_; ++i) {
                T sum = T(0);
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    sum += values_[k] * x[col_indices_[k]];
                }
                y[i] = sum;
            }
        }

        void Multiply(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != cols_) {
                throw VectorDimensionError("SparseMatrixCSR::multiply - vector size must match matrix cols", x.size(), cols_);
            }
            if (&x == &y) {
                const auto input = x;
                Multiply(input, y);
                return;
            }
            y.resize(rows_);
            for (int i = 0; i < rows_; ++i) {
                T sum = T(0);
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    sum += values_[k] * x[col_indices_[k]];
                }
                y[i] = sum;
            }
        }
        
        /// @brief Sparse matrix-vector product: returns A * x
        std::vector<T> operator*(const std::vector<T>& x) const {
            std::vector<T> y;
            Multiply(x, y);
            return y;
        }

        Vector<T> operator*(const Vector<T>& x) const {
            Vector<T> y;
            Multiply(x, y);
            return y;
        }
        
        /// @brief Transpose matrix-vector product: y = A^T * x
        void MultiplyTranspose(const std::vector<T>& x, std::vector<T>& y) const {
            if (x.size() != static_cast<size_t>(rows_)) {
                throw VectorDimensionError("SparseMatrixCSR::multiplyTranspose - vector size must match matrix rows", static_cast<int>(x.size()), rows_);
            }
            if (&x == &y) {
                const auto input = x;
                MultiplyTranspose(input, y);
                return;
            }
            y.assign(cols_, T(0));
            
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    y[col_indices_[k]] += values_[k] * x[i];
                }
            }
        }

        void MultiplyTranspose(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::multiplyTranspose - vector size must match matrix rows", x.size(), rows_);
            }
            if (&x == &y) {
                const auto input = x;
                MultiplyTranspose(input, y);
                return;
            }
            y.resize(cols_);
            for (int j = 0; j < cols_; ++j) y[j] = T(0);
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    y[col_indices_[k]] += values_[k] * x[i];
                }
            }
        }
        
        /// @brief y = alpha * A * x + beta * y
        void gemv(T alpha, const std::vector<T>& x, T beta, std::vector<T>& y) const {
            if (x.size() != static_cast<size_t>(cols_)) {
                throw VectorDimensionError("SparseMatrixCSR::gemv - vector size must match matrix cols", static_cast<int>(x.size()), cols_);
            }
            if (y.size() != static_cast<size_t>(rows_)) {
                throw VectorDimensionError("SparseMatrixCSR::gemv - output vector size must match matrix rows", static_cast<int>(y.size()), rows_);
            }
            if (&x == &y) {
                const auto input = x;
                gemv(alpha, input, beta, y);
                return;
            }
            
            for (int i = 0; i < rows_; ++i) {
                T sum = T(0);
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    sum += values_[k] * x[col_indices_[k]];
                }
                y[i] = alpha * sum + beta * y[i];
            }
        }

        void gemv(T alpha, const Vector<T>& x, T beta, Vector<T>& y) const {
            if (x.size() != cols_) {
                throw VectorDimensionError("SparseMatrixCSR::gemv - vector size must match matrix cols", x.size(), cols_);
            }
            if (y.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::gemv - output vector size must match matrix rows", y.size(), rows_);
            }
            if (&x == &y) {
                const auto input = x;
                gemv(alpha, input, beta, y);
                return;
            }
            for (int i = 0; i < rows_; ++i) {
                T sum = T(0);
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    sum += values_[k] * x[col_indices_[k]];
                }
                y[i] = alpha * sum + beta * y[i];
            }
        }
        
        // ==================== Diagonal Operations ====================
        
        /// @brief Extract main diagonal
        std::vector<T> getDiagonal() const {
            int n = std::min(rows_, cols_);
            std::vector<T> diag(n, T(0));
            
            for (int i = 0; i < n; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    if (col_indices_[k] == i) {
                        diag[i] = values_[k];
                        break;
                    }
                }
            }
            return diag;
        }
        
        /// @brief Get diagonal element
        T diagonal(int i) const {
            for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                if (col_indices_[k] == i) {
                    return values_[k];
                }
            }
            return T(0);
        }
        
        // ==================== Norms ====================
        
        // Norms are real-valued even for complex element types
        using NormType = decltype(std::abs(std::declval<T>()));

        /// @brief Frobenius norm: sqrt(sum of |element|^2)
        NormType NormFrobenius() const {
            using AccumulatorType = std::conditional_t<std::is_floating_point_v<NormType>, NormType, Real>;
            StableNormL2AccumulatorT<AccumulatorType> accumulator;
            for (const auto& v : values_) {
                if constexpr (MMLComplex<T>) {
                    accumulator.add(std::abs(v.real()));
                    accumulator.add(std::abs(v.imag()));
                } else {
                    accumulator.add(std::abs(static_cast<AccumulatorType>(v)));
                }
            }
            return static_cast<NormType>(accumulator.value());
        }
        
        /// @brief Infinity norm: max row sum of absolute values
        NormType NormInf() const {
            NormType max_sum = NormType(0);
            for (int i = 0; i < rows_; ++i) {
                NormType row_sum = NormType(0);
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    row_sum += std::abs(values_[k]);
                }
                max_sum = std::max(max_sum, row_sum);
            }
            return max_sum;
        }
        
        /// @brief One norm: max column sum of absolute values
        NormType Norm1() const {
            if (cols_ == 0) return NormType(0);
            std::vector<NormType> col_sums(cols_, NormType(0));
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    col_sums[col_indices_[k]] += std::abs(values_[k]);
                }
            }
            return *std::max_element(col_sums.begin(), col_sums.end());
        }
        
        // ==================== Matrix Properties ====================
        
        /// @brief Check if matrix is square
        bool IsSquare() const { return rows_ == cols_; }
        
        /// @brief Check if matrix has symmetric sparsity pattern
        bool HasSymmetricPattern() const {
            if (!IsSquare()) return false;
            
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    int j = col_indices_[k];
                    if (findElement(j, i) == nullptr) {
                        return false;
                    }
                }
            }
            return true;
        }
        
        /// @brief Check if matrix is symmetric (pattern AND values); tolerance is a real magnitude
        bool IsSymmetric(SymmetryTolerance tol = PrecisionValues<SymmetryTolerance>::IsMatrixSymmetricTolerance) const {
            if (!IsSquare()) return false;
            
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    int j = col_indices_[k];
                    T a_ij = values_[k];
                    T a_ji = (*this)(j, i);
                    const SymmetryTolerance difference = static_cast<SymmetryTolerance>(std::abs(a_ij - a_ji));
                    const SymmetryTolerance magnitude = static_cast<SymmetryTolerance>(std::abs(a_ij));
                    if (difference > tol * (magnitude + SymmetryTolerance(1))) {
                        return false;
                    }
                }
            }
            return true;
        }
        
        // ==================== Triangular Solves ====================
        
        /// @brief Lower triangular solve: L * x = b (L is this matrix)
        /// @param b Right-hand side
        /// @param x Solution (output)
        /// @note Assumes matrix is lower triangular with non-zero diagonal
        void SolveLower(const std::vector<T>& b, std::vector<T>& x) const {
            if (!IsSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveLower - matrix must be square", rows_, cols_, rows_, rows_);
            if (static_cast<int>(b.size()) != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::solveLower - right-hand side size must match matrix rows", static_cast<int>(b.size()), rows_);
            }
            x = b;
            for (int i = 0; i < rows_; ++i) {
                bool foundDiagonal = false;
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    int j = col_indices_[k];
                    if (j < i) {
                        x[i] -= values_[k] * x[j];
                    } else if (j == i) {
                        if (values_[k] == T(0)) throw SingularMatrixError("SparseMatrixCSR::solveLower - zero diagonal");
                        x[i] /= values_[k];
                        foundDiagonal = true;
                        break;
                    }
                }
                if (!foundDiagonal) throw SingularMatrixError("SparseMatrixCSR::solveLower - missing diagonal");
            }
        }

        void SolveLower(const Vector<T>& b, Vector<T>& x) const {
            if (!IsSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveLower - matrix must be square", rows_, cols_, rows_, rows_);
            if (b.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::solveLower - right-hand side size must match matrix rows", b.size(), rows_);
            }
            x = b;
            for (int i = 0; i < rows_; ++i) {
                bool foundDiagonal = false;
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    const int j = col_indices_[k];
                    if (j < i) x[i] -= values_[k] * x[j];
                    else if (j == i) {
                        if (values_[k] == T(0)) throw SingularMatrixError("SparseMatrixCSR::solveLower - zero diagonal");
                        x[i] /= values_[k];
                        foundDiagonal = true;
                        break;
                    }
                }
                if (!foundDiagonal) throw SingularMatrixError("SparseMatrixCSR::solveLower - missing diagonal");
            }
        }
        
        /// @brief Upper triangular solve: U * x = b (U is this matrix)
        void SolveUpper(const std::vector<T>& b, std::vector<T>& x) const {
            if (!IsSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveUpper - matrix must be square", rows_, cols_, rows_, rows_);
            if (static_cast<int>(b.size()) != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::solveUpper - right-hand side size must match matrix rows", static_cast<int>(b.size()), rows_);
            }
            x = b;
            for (int i = rows_ - 1; i >= 0; --i) {
                T diag = T(0);
                bool foundDiagonal = false;
                for (int k = row_pointers_[i + 1] - 1; k >= row_pointers_[i]; --k) {
                    int j = col_indices_[k];
                    if (j > i) {
                        x[i] -= values_[k] * x[j];
                    } else if (j == i) {
                        diag = values_[k];
                        foundDiagonal = true;
                    }
                }
                if (!foundDiagonal) throw SingularMatrixError("SparseMatrixCSR::solveUpper - missing diagonal");
                if (diag == T(0)) throw SingularMatrixError("SparseMatrixCSR::solveUpper - zero diagonal");
                x[i] /= diag;
            }
        }

        void SolveUpper(const Vector<T>& b, Vector<T>& x) const {
            if (!IsSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveUpper - matrix must be square", rows_, cols_, rows_, rows_);
            if (b.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::solveUpper - right-hand side size must match matrix rows", b.size(), rows_);
            }
            x = b;
            for (int i = rows_ - 1; i >= 0; --i) {
                T diagonalValue = T(0);
                bool foundDiagonal = false;
                for (int k = row_pointers_[i + 1] - 1; k >= row_pointers_[i]; --k) {
                    const int j = col_indices_[k];
                    if (j > i) x[i] -= values_[k] * x[j];
                    else if (j == i) {
                        diagonalValue = values_[k];
                        foundDiagonal = true;
                    }
                }
                if (!foundDiagonal) throw SingularMatrixError("SparseMatrixCSR::solveUpper - missing diagonal");
                if (diagonalValue == T(0)) throw SingularMatrixError("SparseMatrixCSR::solveUpper - zero diagonal");
                x[i] /= diagonalValue;
            }
        }
        
        // ==================== Utility ====================
        
        /// @brief Scale all values by a constant
        void scale(T alpha) {
            for (auto& v : values_) {
                v *= alpha;
            }
        }
        
        // ==================== Row Modification (for BCs) ====================
        
        /// @brief Zero out all entries in a row
        /// @note Used for applying boundary conditions
        void ZeroRow(int row) {
            if (row < 0 || row >= rows_) {
                throw IndexError("SparseMatrixCSR::zeroRow - row index out of bounds", row, rows_);
            }
            for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                values_[k] = T(0);
            }
        }
        
        /// @brief Set element at (row, col) to value
        /// @note Only works if entry already exists in sparsity pattern!
        ///       For boundary conditions, use this after assembling the full matrix.
        void set(int row, int col, T value) {
            if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                throw MatrixAccessBoundsError("SparseMatrixCSR::set - index out of bounds", row, col, rows_, cols_);
            }
            
            for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                if (col_indices_[k] == col) {
                    values_[k] = value;
                    return;
                }
            }
            
            // Entry not found - for BC application we need it to exist
            // This can happen if the sparsity pattern doesn't include the diagonal
            // For now, throw an error (user should ensure pattern is complete)
            throw ArgumentError("SparseMatrixCSR::set - entry not in sparsity pattern");
        }
        
        /// @brief Add to element at (row, col)
        /// @note Only works if entry already exists in sparsity pattern
        void add(int row, int col, T value) {
            if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                throw MatrixAccessBoundsError("SparseMatrixCSR::add - index out of bounds", row, col, rows_, cols_);
            }
            
            for (int k = row_pointers_[row]; k < row_pointers_[row + 1]; ++k) {
                if (col_indices_[k] == col) {
                    values_[k] += value;
                    return;
                }
            }
            
            throw ArgumentError("SparseMatrixCSR::add - entry not in sparsity pattern");
        }
        
        /// @brief Memory usage in bytes
        size_t memoryUsage() const {
            return sizeof(*this) + 
                   values_.capacity() * sizeof(T) +
                   col_indices_.capacity() * sizeof(int) +
                   row_pointers_.capacity() * sizeof(int);
        }
        
        /// @brief Print matrix (for debugging)
        void Print(std::ostream& os = std::cout, int precision = 6) const {
            os << "SparseMatrixCSR " << rows_ << "x" << cols_ 
               << " with " << nnz() << " non-zeros:\n";
            os << std::fixed << std::setprecision(precision);
            
            for (int i = 0; i < rows_; ++i) {
                os << "Row " << i << ": ";
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    os << "(" << col_indices_[k] << "," << values_[k] << ") ";
                }
                os << "\n";
            }
        }
        
        /// @brief Convert to dense matrix (for debugging/testing)
        std::vector<std::vector<T>> ToDense() const {
            std::vector<std::vector<T>> dense(rows_, std::vector<T>(cols_, T(0)));
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    dense[i][col_indices_[k]] = values_[k];
                }
            }
            return dense;
        }

        // ---- Deprecated spellings (CONTRACT_POLICY section 11; removed after one minor release) ----
        [[deprecated("use FromCOO()")]]
        void fromCOO(SparseMatrixCOO<T>& coo) { FromCOO(coo); }
        [[deprecated("use NormFrobenius()")]]
        NormType normFrobenius() const { return NormFrobenius(); }
        [[deprecated("use NormInf()")]]
        NormType normInf() const { return NormInf(); }
        [[deprecated("use Norm1()")]]
        NormType norm1() const { return Norm1(); }
        [[deprecated("use IsSquare()")]]
        bool isSquare() const { return IsSquare(); }
        [[deprecated("use HasSymmetricPattern()")]]
        bool hasSymmetricPattern() const { return HasSymmetricPattern(); }
        [[deprecated("use IsSymmetric()")]]
        bool isSymmetric(SymmetryTolerance tol = PrecisionValues<SymmetryTolerance>::IsMatrixSymmetricTolerance) const { return IsSymmetric(tol); }
        [[deprecated("use ZeroRow()")]]
        void zeroRow(int row) { ZeroRow(row); }
        [[deprecated("use Print()")]]
        void print(std::ostream& os = std::cout, int precision = 6) const { Print(os, precision); }
        [[deprecated("use ToDense()")]]
        std::vector<std::vector<T>> toDense() const { return ToDense(); }
        [[deprecated("use Multiply()")]]
        void multiply(const std::vector<T>& x, std::vector<T>& y) const { Multiply(x, y); }
        [[deprecated("use Multiply()")]]
        void multiply(const Vector<T>& x, Vector<T>& y) const { Multiply(x, y); }
        [[deprecated("use MultiplyTranspose()")]]
        void multiplyTranspose(const std::vector<T>& x, std::vector<T>& y) const { MultiplyTranspose(x, y); }
        [[deprecated("use MultiplyTranspose()")]]
        void multiplyTranspose(const Vector<T>& x, Vector<T>& y) const { MultiplyTranspose(x, y); }
        [[deprecated("use SolveLower()")]]
        void solveLower(const std::vector<T>& b, std::vector<T>& x) const { SolveLower(b, x); }
        [[deprecated("use SolveLower()")]]
        void solveLower(const Vector<T>& b, Vector<T>& x) const { SolveLower(b, x); }
        [[deprecated("use SolveUpper()")]]
        void solveUpper(const std::vector<T>& b, std::vector<T>& x) const { SolveUpper(b, x); }
        [[deprecated("use SolveUpper()")]]
        void solveUpper(const Vector<T>& b, Vector<T>& x) const { SolveUpper(b, x); }
    };
    
    // ==================== Free Functions ====================
    
    /// @brief Create identity matrix in CSR format
    template<typename T>
    SparseMatrixCSR<T> IdentityCSR(int n) {
        if (n < 0)
            throw MatrixDimensionError("IdentityCSR: dimension must be non-negative", n, n, -1, -1);
        std::vector<T> values(n, T(1));
        std::vector<int> col_indices(n);
        std::vector<int> row_pointers(static_cast<size_t>(n) + 1);
        
        for (int i = 0; i < n; ++i) {
            col_indices[i] = i;
            row_pointers[i] = i;
        }
        row_pointers[n] = n;
        
        return SparseMatrixCSR<T>(n, n, std::move(values), 
                                  std::move(col_indices), std::move(row_pointers));
    }

    /// @deprecated Use IdentityCSR().
    template<typename T>
    [[deprecated("use IdentityCSR()")]]
    SparseMatrixCSR<T> identityCSR(int n) { return IdentityCSR<T>(n); }
    
    /// @brief Create diagonal matrix in CSR format
    template<typename T>
    SparseMatrixCSR<T> DiagonalCSR(const std::vector<T>& diag) {
        if (diag.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw MatrixAllocationError("DiagonalCSR: dimension exceeds index range", -1, -1);
        int n = static_cast<int>(diag.size());
        std::vector<T> values = diag;
        std::vector<int> col_indices(n);
        std::vector<int> row_pointers(static_cast<size_t>(n) + 1);
        
        for (int i = 0; i < n; ++i) {
            col_indices[i] = i;
            row_pointers[i] = i;
        }
        row_pointers[n] = n;
        
        return SparseMatrixCSR<T>(n, n, std::move(values), 
                                  std::move(col_indices), std::move(row_pointers));
    }

    /// @deprecated Use DiagonalCSR().
    template<typename T>
    [[deprecated("use DiagonalCSR()")]]
    SparseMatrixCSR<T> diagonalCSR(const std::vector<T>& diag) { return DiagonalCSR<T>(diag); }
    
    /// @brief Compute residual r = b - A*x
    template<typename T>
    std::vector<T> residual(const SparseMatrixCSR<T>& A, 
                            const std::vector<T>& x, 
                            const std::vector<T>& b) {
        if (b.size() != static_cast<size_t>(A.rows()))
            throw VectorDimensionError("Sparse residual: RHS size must match matrix rows", static_cast<int>(b.size()), A.rows());
        std::vector<T> r = A * x;
        for (size_t i = 0; i < r.size(); ++i) {
            r[i] = b[i] - r[i];
        }
        return r;
    }
    
    /// @brief Compute L2 norm of vector: sqrt(sum |x_i|^2). Returns a real magnitude also for complex T.
    template<typename T>
    auto norm2(const std::vector<T>& v) {
        using NormValue = std::conditional_t<std::is_integral_v<T>, double, T>;
        using RealT = decltype(std::abs(NormValue{}));
        StableNormL2AccumulatorT<RealT> accumulator;
        for (const auto& x : v)
            accumulator.add(std::abs(NormValue(x)));
        return accumulator.value();
    }
    
    /// @brief Compute (unconjugated, bilinear) dot product sum a_i * b_i
    template<typename T>
    T dot(const std::vector<T>& a, const std::vector<T>& b) {
        if (a.size() != b.size())
            throw VectorDimensionError("Sparse dot: dimensions must match", static_cast<int>(a.size()), static_cast<int>(b.size()));
        T sum = T(0);
        for (size_t i = 0; i < a.size(); ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }
    
    /// @brief axpy: y = alpha * x + y
    template<typename T>
    void axpy(T alpha, const std::vector<T>& x, std::vector<T>& y) {
        if (x.size() != y.size())
            throw VectorDimensionError("Sparse axpy: dimensions must match", static_cast<int>(x.size()), static_cast<int>(y.size()));
        for (size_t i = 0; i < x.size(); ++i) {
            y[i] += alpha * x[i];
        }
    }

		template<typename T>
		Vector<T> residual(const SparseMatrixCSR<T>& A, const Vector<T>& x, const Vector<T>& b) {
			return b - A * x;
		}

		template<typename T>
		auto norm2(const Vector<T>& values) {
			using NormValue = std::conditional_t<std::is_integral_v<T>, double, T>;
			using Magnitude = decltype(std::abs(NormValue{}));
			StableNormL2AccumulatorT<Magnitude> accumulator;
			for (const auto& value : values)
				accumulator.add(std::abs(NormValue(value)));
			return accumulator.value();
		}

		template<typename T>
		T dot(const Vector<T>& a, const Vector<T>& b) {
			if (a.size() != b.size())
				throw VectorDimensionError("Sparse dot: dimensions must match", a.size(), b.size());
			T sum = T(0);
			for (int i = 0; i < a.size(); ++i)
				sum += a[i] * b[i];
			return sum;
		}

		template<typename T>
		void axpy(T alpha, const Vector<T>& x, Vector<T>& y) {
			if (x.size() != y.size())
				throw VectorDimensionError("Sparse axpy: dimensions must match", x.size(), y.size());
			for (int i = 0; i < x.size(); ++i)
				y[i] += alpha * x[i];
		}

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_CSR_H
