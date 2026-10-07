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

namespace MML::SparseMatrix
{
    ///////////////////////////////////////////////////////////////////////////
    ///                     Factory Functions                               ///
    ///////////////////////////////////////////////////////////////////////////
    
    /// @brief Create 1D Laplacian matrix (-1, 2, -1 stencil)
    /// @param n Number of interior grid points
    /// @param h Grid spacing (default 1.0)
    /// @return Tridiagonal Laplacian matrix in CSR format
    template<typename T = double>
    SparseMatrixCSR<T> laplacian1D(int n, T h = T(1)) {
        SparseMatrixCOO<T> coo(n, n);
        coo.reserve(3 * n);
        
        T h2_inv = T(1) / (h * h);
        
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
    SparseMatrixCSR<T> laplacian2D(int nx, int ny, T hx = T(1), T hy = T(1)) {
        int n = nx * ny;
        SparseMatrixCOO<T> coo(n, n);
        coo.reserve(5 * n);
        
        T hx2_inv = T(1) / (hx * hx);
        T hy2_inv = T(1) / (hy * hy);
        T diag = 2 * hx2_inv + 2 * hy2_inv;
        
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
    SparseMatrixCSR<T> laplacian3D(int nx, int ny, int nz, 
                                    T hx = T(1), T hy = T(1), T hz = T(1)) {
        int n = nx * ny * nz;
        SparseMatrixCOO<T> coo(n, n);
        coo.reserve(7 * n);
        
        T hx2_inv = T(1) / (hx * hx);
        T hy2_inv = T(1) / (hy * hy);
        T hz2_inv = T(1) / (hz * hz);
        T diag = 2 * (hx2_inv + hy2_inv + hz2_inv);
        
        int nxy = nx * ny;
        
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

} // namespace MML::SparseMatrix

#endif // MML_BASE_SPARSE_MATRIX_H
