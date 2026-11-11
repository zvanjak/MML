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
            : rows_(rows), cols_(cols) {
            if (rows < 0 || cols < 0) {
                throw MatrixDimensionError("SparseMatrixCSC: matrix dimensions must be non-negative", rows, cols, -1, -1);
            }
            col_pointers_.assign(static_cast<size_t>(cols) + 1, 0);
        }
        
        /// @brief Construct from COO format
        explicit SparseMatrixCSC(SparseMatrixCOO<T>& coo) {
            FromCOO(coo);
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
        void FromCOO(SparseMatrixCOO<T>& coo) {
            coo.consolidate();
            const auto& triplets = coo.triplets();
            if (triplets.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
                throw MatrixAllocationError("SparseMatrixCSC: nonzero count exceeds index range", coo.rows(), coo.cols());
            int nnz = static_cast<int>(triplets.size());
            std::vector<T> values(nnz);
            std::vector<int> indices(nnz);
            std::vector<int> pointers(static_cast<size_t>(coo.cols()) + 1, 0);
            
            // Count entries per column
            for (const auto& t : triplets) {
                pointers[t.col + 1]++;
            }
            
            // Cumulative sum
            for (int j = 0; j < coo.cols(); ++j) {
                pointers[j + 1] += pointers[j];
            }
            
            // Fill values (need to sort by column first)
            std::vector<Triplet<T>> sorted_triplets = triplets;
            std::sort(sorted_triplets.begin(), sorted_triplets.end(),
                [](const Triplet<T>& a, const Triplet<T>& b) {
                    if (a.col != b.col) return a.col < b.col;
                    return a.row < b.row;
                });
            
            for (int k = 0; k < nnz; ++k) {
                values[k] = sorted_triplets[k].value;
                indices[k] = sorted_triplets[k].row;
            }
            values_.swap(values);
            row_indices_.swap(indices);
            col_pointers_.swap(pointers);
            rows_ = coo.rows();
            cols_ = coo.cols();
        }
        
        /// @brief Build CSC from CSR matrix
        void fromCSR(const SparseMatrixCSR<T>& csr) {
            if (!csr.validate())
                throw ArgumentError("SparseMatrixCSC::fromCSR - invalid CSR storage");
            int nnz = csr.nnz();
            
            std::vector<T> values(nnz);
            std::vector<int> indices(nnz);
            std::vector<int> pointers(static_cast<size_t>(csr.cols()) + 1, 0);
            
            const auto& csr_vals = csr.values();
            const auto& csr_cols = csr.colIndices();
            const auto& csr_rows = csr.rowPointers();
            
            // Count entries per column
            for (int k = 0; k < nnz; ++k) {
                pointers[csr_cols[k] + 1]++;
            }
            
            // Cumulative sum
            for (int j = 0; j < csr.cols(); ++j) {
                pointers[j + 1] += pointers[j];
            }
            
            // Fill values and row indices
            std::vector<int> col_counts(csr.cols(), 0);
            for (int i = 0; i < csr.rows(); ++i) {
                for (int k = csr_rows[i]; k < csr_rows[i + 1]; ++k) {
                    int j = csr_cols[k];
                    int idx = pointers[j] + col_counts[j];
                    values[idx] = csr_vals[k];
                    indices[idx] = i;
                    col_counts[j]++;
                }
            }
            values_.swap(values);
            row_indices_.swap(indices);
            col_pointers_.swap(pointers);
            rows_ = csr.rows();
            cols_ = csr.cols();
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
            if (values_.size() > static_cast<size_t>(std::numeric_limits<int>::max())) return false;
            if (col_pointers_.size() != static_cast<size_t>(cols_) + 1 || col_pointers_.front() != 0) return false;
            if (col_pointers_.back() != static_cast<int>(values_.size())) return false;
            for (int col = 0; col < cols_; ++col) {
                if (col_pointers_[col] < 0 || col_pointers_[col] > col_pointers_[col + 1] ||
                    col_pointers_[col + 1] > nnz()) return false;
            }
            for (int col = 0; col < cols_; ++col) {
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
        void Multiply(const std::vector<T>& x, std::vector<T>& y) const {
            if (x.size() != static_cast<size_t>(cols_)) {
                throw VectorDimensionError("SparseMatrixCSC::multiply - vector size must match matrix cols", static_cast<int>(x.size()), cols_);
            }
            if (&x == &y) {
                const auto input = x;
                Multiply(input, y);
                return;
            }
            y.assign(rows_, T(0));
            
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    y[row_indices_[k]] += values_[k] * x[j];
                }
            }
        }

        void Multiply(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != cols_) {
                throw VectorDimensionError("SparseMatrixCSC::multiply - vector size must match matrix cols", x.size(), cols_);
            }
            if (&x == &y) {
                const auto input = x;
                Multiply(input, y);
                return;
            }
            y.resize(rows_);
            for (int i = 0; i < rows_; ++i) y[i] = T(0);
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    y[row_indices_[k]] += values_[k] * x[j];
                }
            }
        }
        
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
                throw VectorDimensionError("SparseMatrixCSC::multiplyTranspose - vector size must match matrix rows", static_cast<int>(x.size()), rows_);
            }
            if (&x == &y) {
                const auto input = x;
                MultiplyTranspose(input, y);
                return;
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

        void MultiplyTranspose(const Vector<T>& x, Vector<T>& y) const {
            if (x.size() != rows_) {
                throw VectorDimensionError("SparseMatrixCSC::multiplyTranspose - vector size must match matrix rows", x.size(), rows_);
            }
            if (&x == &y) {
                const auto input = x;
                MultiplyTranspose(input, y);
                return;
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
        std::vector<std::vector<T>> ToDense() const {
            std::vector<std::vector<T>> dense(rows_, std::vector<T>(cols_, T(0)));
            for (int j = 0; j < cols_; ++j) {
                for (int k = col_pointers_[j]; k < col_pointers_[j + 1]; ++k) {
                    dense[row_indices_[k]][j] = values_[k];
                }
            }
            return dense;
        }

        // ---- Deprecated spellings (CONTRACT_POLICY section 11; removed after one minor release) ----
        [[deprecated("use FromCOO()")]]
        void fromCOO(SparseMatrixCOO<T>& coo) { FromCOO(coo); }
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
    };

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_CSC_H
