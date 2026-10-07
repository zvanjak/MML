///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SparseMatrixCSC.h                                                   ///
///  Description: Compressed Sparse Column (CSC) format                               ///
///               Efficient for column-wise access and certain algorithms             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_BASE_SPARSE_MATRIX_CSC_H
#define MML_BASE_SPARSE_MATRIX_CSC_H

#include <mml/MMLExceptions.h>

#include "SparseMatrixCOO.h"
#include "SparseMatrixCSR.h"
#include <mml/base/Vector/Vector.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace MML::SparseMatrix
{
    ///////////////////////////////////////////////////////////////////////////
    ///                    SparseMatrixCSC Class                            ///
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief Compressed Sparse Column (CSC) format sparse matrix
    /// 
    /// Storage:
    /// - values_: non-zero values in column-major order
    /// - row_indices_: row index for each value
    /// - col_pointers_: index into values_ where each column starts
    /// 
    /// Ideal for:
    /// - Column-wise access
    /// - Matrix transpose operations (CSC of A = CSR of A^T)
    /// - Some direct solvers
    /// 
    /// @tparam T Element type (typically double or float)
    template<typename T>
    class SparseMatrixCSC {
    private:
        int rows_;
        int cols_;
        std::vector<T> values_;
        std::vector<int> row_indices_;
        std::vector<int> col_pointers_;
        
    public:
        // ==================== Constructors ====================
        
        /// @brief Default constructor
        SparseMatrixCSC() : rows_(0), cols_(0), col_pointers_(1, 0) {}
        
        /// @brief Construct from dimensions (empty matrix)
        SparseMatrixCSC(int rows, int cols) 
            : rows_(rows), cols_(cols), col_pointers_(cols >= 0 ? cols + 1 : 0, 0) {
            if (rows < 0 || cols < 0) {
                throw MatrixDimensionError("SparseMatrixCSC: matrix dimensions must be non-negative", rows, cols, -1, -1);
            }
        }
        
        /// @brief Construct from COO format
        explicit SparseMatrixCSC(SparseMatrixCOO<T>& coo) {
            fromCOO(coo);
        }
        
        /// @brief Construct from CSR format (transpose structure)
        explicit SparseMatrixCSC(const SparseMatrixCSR<T>& csr) {
            fromCSR(csr);
        }
        
        /// @brief Construct from raw CSC data
        SparseMatrixCSC(int rows, int cols,
                       std::vector<T> values,
                       std::vector<int> row_indices,
                       std::vector<int> col_pointers)
            : rows_(rows), cols_(cols),
              values_(std::move(values)),
              row_indices_(std::move(row_indices)),
                            col_pointers_(std::move(col_pointers)) {
                        if (!validate()) throw ArgumentError("SparseMatrixCSC: invalid CSC storage");
                }
        
        // ==================== Conversion ====================
        
        /// @brief Build CSC from COO matrix
        void fromCOO(SparseMatrixCOO<T>& coo) {
            rows_ = coo.rows();
            cols_ = coo.cols();
            
            coo.consolidate();
            const auto& triplets = coo.triplets();
            
            int nnz = static_cast<int>(triplets.size());
            values_.resize(nnz);
            row_indices_.resize(nnz);
            col_pointers_.resize(cols_ + 1);
            
            // Count entries per column
            std::fill(col_pointers_.begin(), col_pointers_.end(), 0);
            for (const auto& t : triplets) {
                col_pointers_[t.col + 1]++;
            }
            
            // Cumulative sum
            for (int j = 0; j < cols_; ++j) {
                col_pointers_[j + 1] += col_pointers_[j];
            }
            
            // Fill values (need to sort by column first)
            std::vector<Triplet<T>> sorted_triplets = triplets;
            std::sort(sorted_triplets.begin(), sorted_triplets.end(),
                [](const Triplet<T>& a, const Triplet<T>& b) {
                    if (a.col != b.col) return a.col < b.col;
                    return a.row < b.row;
                });
            
            for (int k = 0; k < nnz; ++k) {
                values_[k] = sorted_triplets[k].value;
                row_indices_[k] = sorted_triplets[k].row;
            }
        }
        
        /// @brief Build CSC from CSR matrix
        void fromCSR(const SparseMatrixCSR<T>& csr) {
            rows_ = csr.rows();
            cols_ = csr.cols();
            int nnz = csr.nnz();
            
            values_.resize(nnz);
            row_indices_.resize(nnz);
            col_pointers_.assign(cols_ + 1, 0);
            
            const auto& csr_vals = csr.values();
            const auto& csr_cols = csr.colIndices();
            const auto& csr_rows = csr.rowPointers();
            
            // Count entries per column
            for (int k = 0; k < nnz; ++k) {
                col_pointers_[csr_cols[k] + 1]++;
            }
            
            // Cumulative sum
            for (int j = 0; j < cols_; ++j) {
                col_pointers_[j + 1] += col_pointers_[j];
            }
            
            // Fill values and row indices
            std::vector<int> col_counts(cols_, 0);
            for (int i = 0; i < rows_; ++i) {
                for (int k = csr_rows[i]; k < csr_rows[i + 1]; ++k) {
                    int j = csr_cols[k];
                    int idx = col_pointers_[j] + col_counts[j];
                    values_[idx] = csr_vals[k];
                    row_indices_[idx] = i;
                    col_counts[j]++;
                }
            }
        }
        
        /// @brief Convert to CSR format
        SparseMatrixCSR<T> toCSR() const {
            SparseMatrixCOO<T> coo(rows_, cols_);
            coo.reserve(nnz());
            
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    coo.addEntry(row_indices_[k], j, values_[k]);
                }
            }
            
            return SparseMatrixCSR<T>(coo);
        }
        
        // ==================== Dimensions ====================
        
        int rows() const { return rows_; }
        int cols() const { return cols_; }
        int nnz() const { return static_cast<int>(values_.size()); }
        
        // ==================== Raw Data Access ====================
        
        const std::vector<T>& values() const { return values_; }
        const std::vector<int>& rowIndices() const { return row_indices_; }
        const std::vector<int>& colPointers() const { return col_pointers_; }

        /// @brief Check dimensions and canonical CSC storage invariants
        bool validate() const {
            if (rows_ < 0 || cols_ < 0 || values_.size() != row_indices_.size()) return false;
            if (col_pointers_.size() != static_cast<size_t>(cols_) + 1 || col_pointers_.front() != 0) return false;
            if (col_pointers_.back() != static_cast<int>(values_.size())) return false;
            for (int col = 0; col < cols_; ++col) {
                if (col_pointers_[col] < 0 || col_pointers_[col] > col_pointers_[col + 1]) return false;
                int previousRow = -1;
                for (int k = col_pointers_[col]; k < col_pointers_[col + 1]; ++k) {
                    if (row_indices_[k] < 0 || row_indices_[k] >= rows_ || row_indices_[k] <= previousRow) return false;
                    previousRow = row_indices_[k];
                }
            }
            return true;
        }
        
        // ==================== Element Access ====================
        
        /// @brief Get element at (row, col)
        T operator()(int row, int col) const {
            for (int k = col_pointers_[col]; k < col_pointers_[col + 1]; ++k) {
                if (row_indices_[k] == row) {
                    return values_[k];
                }
            }
            return T(0);
        }
        
        /// @brief Number of non-zeros in column j
        int colNnz(int col) const {
            return col_pointers_[col + 1] - col_pointers_[col];
        }
        
        // ==================== Matrix-Vector Operations ====================
        
        /// @brief Sparse matrix-vector product: y = A * x
        void multiply(const std::vector<T>& x, std::vector<T>& y) const {
            if (static_cast<int>(x.size()) != cols_) {
                throw VectorDimensionError("SparseMatrixCSC::multiply - vector size must match matrix cols", static_cast<int>(x.size()), cols_);
            }
            y.assign(rows_, T(0));
            
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    y[row_indices_[k]] += values_[k] * x[j];
                }
            }
        }

        void multiply(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != cols_) {
                throw VectorDimensionError("SparseMatrixCSC::multiply - vector size must match matrix cols", x.size(), cols_);
            }
            y.Resize(rows_);
            for (int i = 0; i < rows_; ++i) y[i] = T(0);
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    y[row_indices_[k]] += values_[k] * x[j];
                }
            }
        }
        
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
                throw VectorDimensionError("SparseMatrixCSC::multiplyTranspose - vector size must match matrix rows", static_cast<int>(x.size()), rows_);
            }
            y.resize(cols_);
            
            for (int j = 0; j < cols_; ++j) {
                T sum = T(0);
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    sum += values_[k] * x[row_indices_[k]];
                }
                y[j] = sum;
            }
        }

        void multiplyTranspose(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSC::multiplyTranspose - vector size must match matrix rows", x.size(), rows_);
            }
            y.Resize(cols_);
            for (int j = 0; j < cols_; ++j) {
                T sum = T(0);
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    sum += values_[k] * x[row_indices_[k]];
                }
                y[j] = sum;
            }
        }
        
        // ==================== Diagonal ====================
        
        std::vector<T> getDiagonal() const {
            int n = std::min(rows_, cols_);
            std::vector<T> diag(n, T(0));
            
            for (int j = 0; j < n; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    if (row_indices_[k] == j) {
                        diag[j] = values_[k];
                        break;
                    }
                }
            }
            return diag;
        }
        
        // ==================== Utility ====================
        
        size_t memoryUsage() const {
            return sizeof(*this) + 
                   values_.capacity() * sizeof(T) +
                   row_indices_.capacity() * sizeof(int) +
                   col_pointers_.capacity() * sizeof(int);
        }
        
        /// @brief Convert to dense matrix
        std::vector<std::vector<T>> toDense() const {
            std::vector<std::vector<T>> dense(rows_, std::vector<T>(cols_, T(0)));
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    dense[row_indices_[k]][j] = values_[k];
                }
            }
            return dense;
        }
    };

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_CSC_H
