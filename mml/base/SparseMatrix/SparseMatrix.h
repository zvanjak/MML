///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SparseMatrix.h                                                      ///
///  Description: Main header for sparse matrix module                                ///
///               Includes all sparse matrix formats and utilities                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_BASE_SPARSE_MATRIX_H
#define MML_BASE_SPARSE_MATRIX_H

#include "SparseMatrixCOO.h"
#include "SparseMatrixCSR.h"
#include "SparseMatrixCSC.h"
#include <initializer_list>
#include <limits>

namespace MML::SparseMatrix
{
    namespace Detail {
        inline int CheckedGridCount(std::initializer_list<int> dimensions, int stencilWidth) {
            bool empty = false;
            for (int dimension : dimensions) {
                if (dimension < 0)
                    throw MatrixDimensionError("Sparse Laplacian: grid dimensions must be non-negative", dimension, -1, -1, -1);
                empty = empty || dimension == 0;
            }
            if (empty) return 0;
            const int limit = std::numeric_limits<int>::max() / stencilWidth;
            int count = 1;
            for (int dimension : dimensions) {
                if (count > limit / dimension)
                    throw MatrixAllocationError("Sparse Laplacian: grid/stencil exceeds index range", dimension, count);
                count *= dimension;
            }
            return count;
        }

        template<typename T>
        T InverseSquaredSpacing(T spacing) {
            if constexpr (requires { spacing > T(0); }) {
                if (!(spacing > T(0)) || !std::isfinite(std::abs(spacing)))
                    throw ArgumentError("Sparse Laplacian: spacing must be finite and positive");
            } else {
                if (!std::isfinite(std::abs(spacing)) || std::abs(spacing) == 0)
                    throw ArgumentError("Sparse Laplacian: spacing must be finite and nonzero");
            }
            const T inverse = T(1) / spacing;
            const T squared = inverse * inverse;
            if (!std::isfinite(std::abs(squared)) || std::abs(squared) == 0)
                throw MatrixNumericalError("Sparse Laplacian: inverse squared spacing is not representable");
            return squared;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    ///                     Factory Functions                               ///
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief Create 1D Laplacian matrix (-1, 2, -1 stencil)
    /// @param n Number of interior grid points
    /// @param h Grid spacing (default 1.0)
    /// @return Tridiagonal Laplacian matrix in CSR format
    template<typename T = double>
    SparseMatrixCSR<T> Laplacian1D(int n, T h = T(1)) {
        n = Detail::CheckedGridCount({n}, 3);
        const T h2_inv = Detail::InverseSquaredSpacing(h);
        if (!std::isfinite(std::abs(T(2) * h2_inv)))
            throw MatrixNumericalError("Sparse Laplacian: diagonal coefficient is not representable");
        SparseMatrixCOO<T> coo(n, n);
        coo.reserve(3 * n);
        
        for (int i = 0; i < n; ++i) {
            if (i > 0) {
                coo.addEntry(i, i - 1, -h2_inv);
            }
            coo.addEntry(i, i, 2 * h2_inv);
            if (i < n - 1) {
                coo.addEntry(i, i + 1, -h2_inv);
            }
        }
        
        return SparseMatrixCSR<T>(coo);
    }
    
    /// @brief Create 2D Laplacian matrix (5-point stencil)
    /// @param nx Number of x grid points
    /// @param ny Number of y grid points
    /// @param hx Grid spacing in x (default 1.0)
    /// @param hy Grid spacing in y (default 1.0)
    /// @return 5-point Laplacian matrix in CSR format (lexicographic ordering)
    template<typename T = double>
    SparseMatrixCSR<T> Laplacian2D(int nx, int ny, T hx = T(1), T hy = T(1)) {
        const int n = Detail::CheckedGridCount({nx, ny}, 5);
        const T hx2_inv = Detail::InverseSquaredSpacing(hx);
        const T hy2_inv = Detail::InverseSquaredSpacing(hy);
        const T diag = T(2) * hx2_inv + T(2) * hy2_inv;
        if (!std::isfinite(std::abs(diag)))
            throw MatrixNumericalError("Sparse Laplacian: diagonal coefficient is not representable");
        SparseMatrixCOO<T> coo(n, n);
        coo.reserve(5 * n);
        if (n == 0) return SparseMatrixCSR<T>(coo);
        
        
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = i + j * nx;
                
                // Center
                coo.addEntry(idx, idx, diag);
                
                // Left neighbor
                if (i > 0) {
                    coo.addEntry(idx, idx - 1, -hx2_inv);
                }
                
                // Right neighbor
                if (i < nx - 1) {
                    coo.addEntry(idx, idx + 1, -hx2_inv);
                }
                
                // Bottom neighbor
                if (j > 0) {
                    coo.addEntry(idx, idx - nx, -hy2_inv);
                }
                
                // Top neighbor
                if (j < ny - 1) {
                    coo.addEntry(idx, idx + nx, -hy2_inv);
                }
            }
        }
        
        return SparseMatrixCSR<T>(coo);
    }
    
    /// @brief Create 3D Laplacian matrix (7-point stencil)
    template<typename T = double>
    SparseMatrixCSR<T> Laplacian3D(int nx, int ny, int nz, 
                                    T hx = T(1), T hy = T(1), T hz = T(1)) {
        const int n = Detail::CheckedGridCount({nx, ny, nz}, 7);
        const T hx2_inv = Detail::InverseSquaredSpacing(hx);
        const T hy2_inv = Detail::InverseSquaredSpacing(hy);
        const T hz2_inv = Detail::InverseSquaredSpacing(hz);
        const T diag = T(2) * hx2_inv + T(2) * hy2_inv + T(2) * hz2_inv;
        if (!std::isfinite(std::abs(diag)))
            throw MatrixNumericalError("Sparse Laplacian: diagonal coefficient is not representable");
        SparseMatrixCOO<T> coo(n, n);
        coo.reserve(7 * n);
        if (n == 0) return SparseMatrixCSR<T>(coo);
        
        const int nxy = n == 0 ? 0 : nx * ny;
        
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int idx = i + j * nx + k * nxy;
                    
                    coo.addEntry(idx, idx, diag);
                    
                    if (i > 0) coo.addEntry(idx, idx - 1, -hx2_inv);
                    if (i < nx - 1) coo.addEntry(idx, idx + 1, -hx2_inv);
                    if (j > 0) coo.addEntry(idx, idx - nx, -hy2_inv);
                    if (j < ny - 1) coo.addEntry(idx, idx + nx, -hy2_inv);
                    if (k > 0) coo.addEntry(idx, idx - nxy, -hz2_inv);
                    if (k < nz - 1) coo.addEntry(idx, idx + nxy, -hz2_inv);
                }
            }
        }
        
        return SparseMatrixCSR<T>(coo);
    }

    // ---- Deprecated spellings (CONTRACT_POLICY section 11: factories are PascalCase) ----
    template<typename T = double>
    [[deprecated("use Laplacian1D()")]]
    SparseMatrixCSR<T> laplacian1D(int n, T h = T(1)) { return Laplacian1D<T>(n, h); }
    template<typename T = double>
    [[deprecated("use Laplacian2D()")]]
    SparseMatrixCSR<T> laplacian2D(int nx, int ny, T hx = T(1), T hy = T(1)) { return Laplacian2D<T>(nx, ny, hx, hy); }
    template<typename T = double>
    [[deprecated("use Laplacian3D()")]]
    SparseMatrixCSR<T> laplacian3D(int nx, int ny, int nz, T hx = T(1), T hy = T(1), T hz = T(1)) { return Laplacian3D<T>(nx, ny, nz, hx, hy, hz); }

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_H
