///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SparseMatrixCOO.h                                                   ///
///  Description: Coordinate (COO) format sparse matrix - builder pattern             ///
///               Easy construction via triplets, converts to CSR/CSC for computation ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_BASE_SPARSE_MATRIX_COO_H
#define MML_BASE_SPARSE_MATRIX_COO_H

#include <mml/MMLExceptions.h>

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <tuple>

namespace MML::SparseMatrix
{
    ///////////////////////////////////////////////////////////////////////////
    ///                         Triplet Structure                           ///
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief A single (row, col, value) entry for sparse matrix construction
    template<typename T>
    struct Triplet {
        int row;
        int col;
        T value;
        
        Triplet(int r, int c, T v) : row(r), col(c), value(v) {}
        
        // For sorting: row-major order
        bool operator<(const Triplet& other) const {
            if (row != other.row) return row < other.row;
            return col < other.col;
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    ///                    SparseMatrixCOO Class                            ///
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief Coordinate (COO) format sparse matrix
    /// 
    /// This format stores triplets (row, col, value) and is ideal for:
    /// - Incremental matrix construction
    /// - Adding entries in any order
    /// - Assembling from finite difference/element stencils
    /// 
    /// After construction, convert to CSR or CSC for efficient computation.
    /// 
    /// @tparam T Element type (typically double or float)
    template<typename T>
    class SparseMatrixCOO {
    private:
        int rows_;
        int cols_;
        std::vector<Triplet<T>> triplets_;
        bool sorted_ = false;
        
    public:
        /// @brief Construct empty COO matrix with given dimensions
        SparseMatrixCOO(int rows, int cols) : rows_(rows), cols_(cols) {
            if (rows < 0 || cols < 0) {
                throw MatrixDimensionError("SparseMatrixCOO: matrix dimensions must be non-negative", rows, cols, -1, -1);
            }
        }
        
        /// @brief Default constructor - creates 0x0 matrix
        SparseMatrixCOO() : rows_(0), cols_(0) {}
        
        // ==================== Dimensions ====================
        
        int rows() const { return rows_; }
        int cols() const { return cols_; }
        int nnz() const { return static_cast<int>(triplets_.size()); }
        
        /// @brief Resize the matrix (clears all entries)
        void resize(int rows, int cols) {
            if (rows < 0 || cols < 0) {
                throw MatrixDimensionError("SparseMatrixCOO::resize - matrix dimensions must be non-negative", rows, cols, -1, -1);
            }
            rows_ = rows;
            cols_ = cols;
            triplets_.clear();
            sorted_ = false;
        }
        
        // ==================== Entry Addition ====================
        
        /// @brief Add a single entry (row, col, value)
        /// @note Duplicate entries are allowed - they will be summed during conversion
        void addEntry(int row, int col, T value) {
            if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                throw MatrixAccessBoundsError("SparseMatrixCOO::addEntry - index out of bounds", row, col, rows_, cols_);
            }
            if (value != T(0)) {  // Skip explicit zeros
                triplets_.emplace_back(row, col, value);
                sorted_ = false;
            }
        }
        
        /// @brief Add entry using Triplet
        void addEntry(const Triplet<T>& triplet) {
            addEntry(triplet.row, triplet.col, triplet.value);
        }
        
        /// @brief Add multiple entries at once
        void addEntries(const std::vector<Triplet<T>>& entries) {
            for (const auto& t : entries) {
                addEntry(t);
            }
        }
        
        /// @brief Reserve space for expected number of non-zeros
        void reserve(int expected_nnz) {
            triplets_.reserve(expected_nnz);
        }
        
        /// @brief Clear all entries (keep dimensions)
        void clear() {
            triplets_.clear();
            sorted_ = false;
        }
        
        // ==================== Stencil Helpers ====================
        
        /// @brief Add a row of the matrix (useful for equation assembly)
        /// @param row Row index
        /// @param cols Column indices
        /// @param values Corresponding values
        void addRow(int row, const std::vector<int>& cols, const std::vector<T>& values) {
            if (cols.size() != values.size()) {
                throw ArgumentError("SparseMatrixCOO::addRow - cols and values must have same size");
            }
            for (size_t i = 0; i < cols.size(); ++i) {
                addEntry(row, cols[i], values[i]);
            }
        }
        
        /// @brief Add 1D Laplacian stencil [1, -2, 1] / h^2 for interior node
        void addLaplacian1D(int row, int i, T h_squared_inv) {
            if (i > 0) addEntry(row, i - 1, h_squared_inv);
            addEntry(row, i, -2 * h_squared_inv);
            if (i < cols_ - 1) addEntry(row, i + 1, h_squared_inv);
        }
        
        /// @brief Add a bounded 2D Laplacian 5-point stencil
        /// @param row Equation row
        /// @param i x-index in grid
        /// @param j y-index in grid
        /// @param nx Number of x points
        /// @param ny Number of y points
        /// @param h_squared_inv 1/h^2 (assuming uniform grid)
        void addLaplacian2D_5pt(int row, int i, int j, int nx, int ny, T h_squared_inv) {
            if (nx <= 0 || ny <= 0 || i < 0 || i >= nx || j < 0 || j >= ny || nx * ny > cols_) {
                throw ArgumentError("SparseMatrixCOO::addLaplacian2D_5pt - invalid grid dimensions or node index");
            }
            int idx = i + j * nx;
            // Center
            addEntry(row, idx, -4 * h_squared_inv);
            // Left
            if (i > 0) addEntry(row, idx - 1, h_squared_inv);
            // Right
            if (i < nx - 1) addEntry(row, idx + 1, h_squared_inv);
            // Down
            if (j > 0) addEntry(row, idx - nx, h_squared_inv);
            // Up
            if (j < ny - 1) addEntry(row, idx + nx, h_squared_inv);
        }
        
        // ==================== Access ====================
        
        /// @brief Get all triplets (const reference)
        const std::vector<Triplet<T>>& triplets() const { return triplets_; }
        
        /// @brief Sort triplets in row-major order and sum duplicates
        void consolidate() {
            if (triplets_.empty()) {
                sorted_ = true;
                return;
            }
            
            // Sort by (row, col)
            std::sort(triplets_.begin(), triplets_.end());
            
            // Sum duplicates
            std::vector<Triplet<T>> consolidated;
            consolidated.reserve(triplets_.size());
            
            consolidated.push_back(triplets_[0]);
            for (size_t i = 1; i < triplets_.size(); ++i) {
                if (triplets_[i].row == consolidated.back().row &&
                    triplets_[i].col == consolidated.back().col) {
                    consolidated.back().value += triplets_[i].value;
                } else {
                    if (consolidated.back().value == T(0)) {
                        consolidated.pop_back();
                    }
                    consolidated.push_back(triplets_[i]);
                }
            }
            
            // Remove trailing zero if any
            if (!consolidated.empty() && consolidated.back().value == T(0)) {
                consolidated.pop_back();
            }
            
            triplets_ = std::move(consolidated);
            sorted_ = true;
        }
        
        bool isSorted() const { return sorted_; }

        /// @brief Check dimensions, entry bounds, and the consolidated ordering invariant
        bool validate() const {
            if (rows_ < 0 || cols_ < 0) return false;
            for (size_t i = 0; i < triplets_.size(); ++i) {
                const auto& triplet = triplets_[i];
                if (triplet.row < 0 || triplet.row >= rows_ ||
                    triplet.col < 0 || triplet.col >= cols_) return false;
                if (sorted_ && i > 0 && !(triplets_[i - 1] < triplet)) return false;
            }
            return true;
        }
        
        // ==================== Utility ====================
        
        /// @brief Check if matrix is structurally symmetric
        bool isStructurallySymmetric() const {
            if (rows_ != cols_) return false;
            
            // Build a set of (row, col) pairs
            std::vector<std::pair<int, int>> entries;
            entries.reserve(triplets_.size());
            for (const auto& t : triplets_) {
                entries.emplace_back(t.row, t.col);
            }
            std::sort(entries.begin(), entries.end());
            
            // Check each entry has its transpose
            for (const auto& [r, c] : entries) {
                if (!std::binary_search(entries.begin(), entries.end(), std::make_pair(c, r))) {
                    return false;
                }
            }
            return true;
        }
        
        /// @brief Get memory usage in bytes (approximate)
        size_t memoryUsage() const {
            return sizeof(*this) + triplets_.capacity() * sizeof(Triplet<T>);
        }
    };

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_COO_H
