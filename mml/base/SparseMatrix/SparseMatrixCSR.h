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
            : rows_(rows), cols_(cols), row_pointers_(rows >= 0 ? rows + 1 : 0, 0) {
            if (rows < 0 || cols < 0) {
                throw MatrixDimensionError("SparseMatrixCSR: matrix dimensions must be non-negative", rows, cols, -1, -1);
            }
        }
        
        /// @brief Construct from COO format (main construction path)
        explicit SparseMatrixCSR(SparseMatrixCOO<T>& coo) {
            fromCOO(coo);
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
        void fromCOO(SparseMatrixCOO<T>& coo) {
            rows_ = coo.rows();
            cols_ = coo.cols();
            
            // Consolidate COO (sort and sum duplicates)
            coo.consolidate();
            const auto& triplets = coo.triplets();
            
            int nnz = static_cast<int>(triplets.size());
            values_.resize(nnz);
            col_indices_.resize(nnz);
            row_pointers_.resize(rows_ + 1);
            
            // Count entries per row
            std::fill(row_pointers_.begin(), row_pointers_.end(), 0);
            for (const auto& t : triplets) {
                row_pointers_[t.row + 1]++;
            }
            
            // Cumulative sum to get row pointers
            for (int i = 0; i < rows_; ++i) {
                row_pointers_[i + 1] += row_pointers_[i];
            }
            
            // Fill values and column indices
            std::vector<int> row_counts(rows_, 0);
            for (const auto& t : triplets) {
                int idx = row_pointers_[t.row] + row_counts[t.row];
                values_[idx] = t.value;
                col_indices_[idx] = t.col;
                row_counts[t.row]++;
            }
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
            if (row_pointers_.size() != static_cast<size_t>(rows_) + 1 || row_pointers_.front() != 0) return false;
            if (row_pointers_.back() != static_cast<int>(values_.size())) return false;
            for (int row = 0; row < rows_; ++row) {
                if (row_pointers_[row] < 0 || row_pointers_[row] > row_pointers_[row + 1]) return false;
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
        void multiply(const std::vector<T>& x, std::vector<T>& y) const {
            if (static_cast<int>(x.size()) != cols_) {
                throw VectorDimensionError("SparseMatrixCSR::multiply - vector size must match matrix cols", static_cast<int>(x.size()), cols_);
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

        void multiply(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != cols_) {
                throw VectorDimensionError("SparseMatrixCSR::multiply - vector size must match matrix cols", x.size(), cols_);
            }
            y.Resize(rows_);
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
            multiply(x, y);
            return y;
        }

        Vector<T> operator*(const Vector<T>& x) const {
            Vector<T> y;
            multiply(x, y);
            return y;
        }
        
        /// @brief Transpose matrix-vector product: y = A^T * x
        void multiplyTranspose(const std::vector<T>& x, std::vector<T>& y) const {
            if (static_cast<int>(x.size()) != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::multiplyTranspose - vector size must match matrix rows", static_cast<int>(x.size()), rows_);
            }
            y.assign(cols_, T(0));
            
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    y[col_indices_[k]] += values_[k] * x[i];
                }
            }
        }

        void multiplyTranspose(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::multiplyTranspose - vector size must match matrix rows", x.size(), rows_);
            }
            y.Resize(cols_);
            for (int j = 0; j < cols_; ++j) y[j] = T(0);
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    y[col_indices_[k]] += values_[k] * x[i];
                }
            }
        }
        
        /// @brief y = alpha * A * x + beta * y
        void gemv(T alpha, const std::vector<T>& x, T beta, std::vector<T>& y) const {
            if (static_cast<int>(x.size()) != cols_) {
                throw VectorDimensionError("SparseMatrixCSR::gemv - vector size must match matrix cols", static_cast<int>(x.size()), cols_);
            }
            if (static_cast<int>(y.size()) != rows_) {
                throw VectorDimensionError("SparseMatrixCSR::gemv - output vector size must match matrix rows", static_cast<int>(y.size()), rows_);
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
        NormType normFrobenius() const {
            NormType sum = NormType(0);
            for (const auto& v : values_) {
                sum += static_cast<NormType>(std::norm(v));  // |v|^2 for complex, v^2 for real
            }
            return std::sqrt(sum);
        }
        
        /// @brief Infinity norm: max row sum of absolute values
        NormType normInf() const {
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
        NormType norm1() const {
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
        bool isSquare() const { return rows_ == cols_; }
        
        /// @brief Check if matrix has symmetric sparsity pattern
        bool hasSymmetricPattern() const {
            if (!isSquare()) return false;
            
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
        
        /// @brief Check if matrix is symmetric (pattern AND values)
        bool isSymmetric(T tol = T(1e-12)) const {
            if (!isSquare()) return false;
            
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    int j = col_indices_[k];
                    T a_ij = values_[k];
                    T a_ji = (*this)(j, i);
                    if (std::abs(a_ij - a_ji) > tol * (std::abs(a_ij) + T(1))) {
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
        void solveLower(const std::vector<T>& b, std::vector<T>& x) const {
            if (!isSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveLower - matrix must be square", rows_, cols_, rows_, rows_);
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

        void solveLower(const Vector<T>& b, Vector<T>& x) const {
            if (!isSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveLower - matrix must be square", rows_, cols_, rows_, rows_);
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
        void solveUpper(const std::vector<T>& b, std::vector<T>& x) const {
            if (!isSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveUpper - matrix must be square", rows_, cols_, rows_, rows_);
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

        void solveUpper(const Vector<T>& b, Vector<T>& x) const {
            if (!isSquare()) throw MatrixDimensionError("SparseMatrixCSR::solveUpper - matrix must be square", rows_, cols_, rows_, rows_);
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
        void zeroRow(int row) {
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
        void print(std::ostream& os = std::cout, int precision = 6) const {
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
        std::vector<std::vector<T>> toDense() const {
            std::vector<std::vector<T>> dense(rows_, std::vector<T>(cols_, T(0)));
            for (int i = 0; i < rows_; ++i) {
                for (int k = row_pointers_[i]; k < row_pointers_[i + 1]; ++k) {
                    dense[i][col_indices_[k]] = values_[k];
                }
            }
            return dense;
        }
    };
    
    // ==================== Free Functions ====================
    
    /// @brief Create identity matrix in CSR format
    template<typename T>
    SparseMatrixCSR<T> identityCSR(int n) {
        std::vector<T> values(n, T(1));
        std::vector<int> col_indices(n);
        std::vector<int> row_pointers(n + 1);
        
        for (int i = 0; i < n; ++i) {
            col_indices[i] = i;
            row_pointers[i] = i;
        }
        row_pointers[n] = n;
        
        return SparseMatrixCSR<T>(n, n, std::move(values), 
                                  std::move(col_indices), std::move(row_pointers));
    }
    
    /// @brief Create diagonal matrix in CSR format
    template<typename T>
    SparseMatrixCSR<T> diagonalCSR(const std::vector<T>& diag) {
        int n = static_cast<int>(diag.size());
        std::vector<T> values = diag;
        std::vector<int> col_indices(n);
        std::vector<int> row_pointers(n + 1);
        
        for (int i = 0; i < n; ++i) {
            col_indices[i] = i;
            row_pointers[i] = i;
        }
        row_pointers[n] = n;
        
        return SparseMatrixCSR<T>(n, n, std::move(values), 
                                  std::move(col_indices), std::move(row_pointers));
    }
    
    /// @brief Compute residual r = b - A*x
    template<typename T>
    std::vector<T> residual(const SparseMatrixCSR<T>& A, 
                            const std::vector<T>& x, 
                            const std::vector<T>& b) {
        std::vector<T> r = A * x;
        for (size_t i = 0; i < r.size(); ++i) {
            r[i] = b[i] - r[i];
        }
        return r;
    }
    
    /// @brief Compute L2 norm of vector
    template<typename T>
    T norm2(const std::vector<T>& v) {
        T sum = T(0);
        for (const auto& x : v) {
            sum += x * x;
        }
        return std::sqrt(sum);
    }
    
    /// @brief Compute dot product
    template<typename T>
    T dot(const std::vector<T>& a, const std::vector<T>& b) {
        T sum = T(0);
        for (size_t i = 0; i < a.size(); ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }
    
    /// @brief axpy: y = alpha * x + y
    template<typename T>
    void axpy(T alpha, const std::vector<T>& x, std::vector<T>& y) {
        for (size_t i = 0; i < x.size(); ++i) {
            y[i] += alpha * x[i];
        }
    }

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_CSR_H
