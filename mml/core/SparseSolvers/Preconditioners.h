///////////////////////////////////////////////////////////////////////////////////////////
// Preconditioners.h - Preconditioners for iterative solvers
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Preconditioners accelerate convergence of iterative methods by
// approximately solving M*z = r where M ≈ A. Good preconditioners
// make M⁻¹A close to identity.
//
// Included Preconditioners:
// - JacobiPreconditioner: M = diag(A) - simplest, parallelizable
// - SSORPreconditioner: Symmetric SOR - good for elliptic PDEs
// - ILU0Preconditioner: Incomplete LU - best for general sparse systems
//
// Usage:
//   auto prec = std::make_shared<JacobiPreconditioner<double>>();
//   prec->setup(A);
//   solver.setPreconditioner(prec);
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_CORE_SPARSE_SOLVERS_PRECONDITIONERS_H
#define MML_CORE_SPARSE_SOLVERS_PRECONDITIONERS_H

#include "IterativeSolverBase.h"
#include <algorithm>
#include <stdexcept>

namespace MML::SparseSolvers {

//=============================================================================
// Jacobi (Diagonal) Preconditioner
//=============================================================================

// M = diag(A), so M⁻¹r = r / diag(A)
// Simplest preconditioner, O(n) work, highly parallelizable
template<typename T>
class JacobiPreconditioner : public Preconditioner<T> {
    std::vector<T> invDiag_;  // 1 / diag(A)
    
public:
    void setup(const SparseMatrixCSR<T>& A) override {
        auto diag = A.getDiagonal();
        invDiag_.resize(diag.size());
        
        for (size_t i = 0; i < diag.size(); ++i) {
            // Guard against zero or near-zero diagonal
            if (std::abs(diag[i]) > std::numeric_limits<T>::epsilon() * 100) {
                invDiag_[i] = T(1) / diag[i];
            } else {
                invDiag_[i] = T(1);  // Fall back to identity for this row
            }
        }
    }
    
    void apply(const std::vector<T>& r, std::vector<T>& z) const override {
        if (invDiag_.empty()) {
            throw std::logic_error("JacobiPreconditioner::setup() must be called before apply()");
        }
        if (invDiag_.size() != r.size()) {
            throw std::invalid_argument("JacobiPreconditioner input size does not match setup matrix");
        }
        z.resize(r.size());
        for (size_t i = 0; i < r.size(); ++i) {
            z[i] = r[i] * invDiag_[i];
        }
    }
};

//=============================================================================
// SSOR (Symmetric Successive Over-Relaxation) Preconditioner
//=============================================================================

// M = (D + ωL) D⁻¹ (D + ωU) / (ω(2-ω))
// where A = L + D + U (lower/diagonal/upper)
// 
// Applied as: (D + ωL) z' = r, then (D + ωU) z = (ω(2-ω)) D z'
// Good for elliptic PDEs, ω ∈ (0, 2), typically ω ≈ 1.0-1.5
template<typename T>
class SSORPreconditioner : public Preconditioner<T> {
    const SparseMatrixCSR<T>* A_ = nullptr;
    std::vector<T> diag_;
    T omega_ = T(1.0);  // Relaxation parameter (0 < omega < 2)
    
public:
    explicit SSORPreconditioner(T omega = T(1.0)) : omega_(omega) {}
    
    void setOmega(T omega) { omega_ = omega; }
    T omega() const { return omega_; }
    
    void setup(const SparseMatrixCSR<T>& A) override {
        if (!A.isSquare()) {
            throw MatrixDimensionError("SSORPreconditioner::setup - matrix must be square", A.rows(), A.cols(), A.rows(), A.rows());
        }
        A_ = &A;
        diag_ = A.getDiagonal();
        for (const auto& diagonalValue : diag_) {
            if (diagonalValue == T(0)) {
                A_ = nullptr;
                diag_.clear();
                throw SingularMatrixError("SSORPreconditioner::setup - missing or zero diagonal");
            }
        }
    }
    
    // Apply M⁻¹r where M is the SSOR preconditioner
    // M = (D/ω + L) D⁻¹ (D/ω + U) / (2/ω - 1)
    // Computing z = M⁻¹r:
    //   1. Solve (D/ω + L) y = r  (forward substitution)
    //   2. Solve (D/ω + U) z = D/ω * y  (backward substitution)
    void apply(const std::vector<T>& r, std::vector<T>& z) const override {
        const int n = static_cast<int>(r.size());
        if (!A_) {
            throw std::logic_error("SSORPreconditioner::setup() must be called before apply()");
        }
        if (A_->rows() != n) {
            throw std::invalid_argument("SSORPreconditioner input size does not match setup matrix");
        }
        z.resize(n);
        
        const auto& rowPtr = A_->rowPointers();
        const auto& colIdx = A_->colIndices();
        const auto& values = A_->values();
        
        std::vector<T> y(n);
        
        // Forward sweep: (D/ω + L) y = r
        // y[i] = ω * (r[i] - Σ_{j<i} L[i,j] * y[j]) / D[i]
        for (int i = 0; i < n; ++i) {
            T sum = r[i];
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                int j = colIdx[k];
                if (j < i) {
                    sum -= values[k] * y[j];
                }
            }
            y[i] = omega_ * sum / diag_[i];
        }
        
        // Backward sweep: (D/ω + U) z = D/ω * y
        // z[i] = ω * (D[i]/ω * y[i] - Σ_{j>i} U[i,j] * z[j]) / D[i]
        //      = y[i] - ω * Σ_{j>i} U[i,j] * z[j] / D[i]
        for (int i = n - 1; i >= 0; --i) {
            T sum = y[i];
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                int j = colIdx[k];
                if (j > i) {
                    sum -= omega_ * values[k] * z[j] / diag_[i];
                }
            }
            z[i] = sum;
        }
        
        // Scale by (2-ω)/ω to make M symmetric
        // This is optional but makes the preconditioner more consistent
        // for (int i = 0; i < n; ++i) {
        //     z[i] *= (T(2) - omega_) / omega_;
        // }
    }
};

//=============================================================================
// ILU(0) Preconditioner - Incomplete LU with no fill-in
//=============================================================================

// Computes L and U such that LU ≈ A with the same sparsity pattern as A.
// This is the most effective simple preconditioner for general sparse systems.
//
// Algorithm: IKJ variant
// For k = 1 to n-1:
//   For i = k+1 to n where A[i,k] ≠ 0:
//     A[i,k] = A[i,k] / A[k,k]
//     For j = k+1 to n where A[i,j] ≠ 0:
//       A[i,j] = A[i,j] - A[i,k] * A[k,j]
template<typename T>
class ILU0Preconditioner : public Preconditioner<T> {
    // Store L and U in CSR format (combined, L has unit diagonal)
    SparseMatrixCSR<T> LU_;
    int n_ = 0;
    
public:
    void setup(const SparseMatrixCSR<T>& A) override {
        if (!A.isSquare()) {
            throw MatrixDimensionError("ILU0Preconditioner::setup - matrix must be square", A.rows(), A.cols(), A.rows(), A.rows());
        }
        n_ = A.rows();
        
        // Copy A into working storage
        // We'll modify it in-place to compute L\U
        std::vector<T> values = A.values();
        const auto& rowPtr = A.rowPointers();
        const auto& colIdx = A.colIndices();
        
        // Create maps for fast column lookup within each row
        // diagPtr[i] = index of diagonal in row i
        std::vector<int> diagPtr(n_, -1);
        for (int i = 0; i < n_; ++i) {
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                if (colIdx[k] == i) {
                    diagPtr[i] = k;
                    break;
                }
            }
            if (diagPtr[i] < 0 || values[diagPtr[i]] == T(0)) {
                n_ = 0;
                throw SingularMatrixError("ILU0Preconditioner::setup - missing or zero diagonal");
            }
        }
        
        // IKJ factorization
        for (int i = 1; i < n_; ++i) {
            // For each entry A[i,k] where k < i (in L part)
            for (int pk = rowPtr[i]; pk < rowPtr[i + 1] && colIdx[pk] < i; ++pk) {
                int k = colIdx[pk];
                
                // A[i,k] = A[i,k] / A[k,k]
                if (values[diagPtr[k]] == T(0)) {
                    n_ = 0;
                    throw SingularMatrixError("ILU0Preconditioner::setup - zero factorization pivot");
                }
                values[pk] /= values[diagPtr[k]];
                
                T lik = values[pk];
                
                // For each entry A[k,j] where j > k (in U part of row k)
                // Update A[i,j] -= A[i,k] * A[k,j] if A[i,j] exists
                for (int pj = diagPtr[k] + 1; pj < rowPtr[k + 1]; ++pj) {
                    int j = colIdx[pj];
                    T ukj = values[pj];
                    
                    // Find A[i,j] in row i
                    for (int pi = pk + 1; pi < rowPtr[i + 1]; ++pi) {
                        if (colIdx[pi] == j) {
                            values[pi] -= lik * ukj;
                            break;
                        }
                    }
                }
            }
            if (values[diagPtr[i]] == T(0)) {
                n_ = 0;
                throw SingularMatrixError("ILU0Preconditioner::setup - zero factorization pivot");
            }
        }
        
        // Store the factorization - make copies since we need to pass by value for move
        std::vector<T> values_copy = values;
        std::vector<int> colIdx_copy = colIdx;
        std::vector<int> rowPtr_copy = rowPtr;
        LU_ = SparseMatrixCSR<T>(A.rows(), A.cols(), 
                                  std::move(values_copy),
                                  std::move(colIdx_copy),
                                  std::move(rowPtr_copy));
    }
    
    // Apply ILU: solve L*U*z = r
    // First solve L*y = r (forward substitution)
    // Then solve U*z = y (backward substitution)
    void apply(const std::vector<T>& r, std::vector<T>& z) const override {
        if (n_ == 0) {
            throw std::logic_error("ILU0Preconditioner::setup() must be called before apply()");
        }
        if (n_ != static_cast<int>(r.size())) {
            throw std::invalid_argument("ILU0Preconditioner input size does not match setup matrix");
        }
        z.resize(n_);
        
        const auto& rowPtr = LU_.rowPointers();
        const auto& colIdx = LU_.colIndices();
        const auto& values = LU_.values();
        
        // Find diagonal indices
        std::vector<int> diagPtr(n_, -1);
        for (int i = 0; i < n_; ++i) {
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                if (colIdx[k] == i) {
                    diagPtr[i] = k;
                    break;
                }
            }
        }
        
        // Forward solve: L*y = r (L has unit diagonal)
        std::vector<T> y(n_);
        for (int i = 0; i < n_; ++i) {
            T sum = r[i];
            for (int k = rowPtr[i]; k < rowPtr[i + 1] && colIdx[k] < i; ++k) {
                sum -= values[k] * y[colIdx[k]];
            }
            y[i] = sum;
        }
        
        // Backward solve: U*z = y
        for (int i = n_ - 1; i >= 0; --i) {
            T sum = y[i];
            for (int k = diagPtr[i] + 1; k < rowPtr[i + 1]; ++k) {
                sum -= values[k] * z[colIdx[k]];
            }
            if (diagPtr[i] >= 0 && std::abs(values[diagPtr[i]]) > std::numeric_limits<T>::epsilon() * 100) {
                z[i] = sum / values[diagPtr[i]];
            } else {
                z[i] = sum;  // Fallback for zero diagonal
            }
        }
    }
};

//=============================================================================
// Block Jacobi Preconditioner
//=============================================================================

// Divides the matrix into blocks and inverts each block independently.
// Good for parallel implementations and multi-physics problems.
template<typename T>
class BlockJacobiPreconditioner : public Preconditioner<T> {
    std::vector<std::vector<std::vector<T>>> blockInverses_;  // Dense block inverses
    std::vector<int> blockStarts_;  // Starting index of each block
    int blockSize_ = 4;
    
public:
    explicit BlockJacobiPreconditioner(int blockSize = 4) : blockSize_(blockSize) {}
    
    void setBlockSize(int bs) { blockSize_ = bs; }
    
    void setup(const SparseMatrixCSR<T>& A) override {
        int n = A.rows();
        int numBlocks = (n + blockSize_ - 1) / blockSize_;
        
        blockStarts_.resize(numBlocks + 1);
        blockInverses_.resize(numBlocks);
        
        for (int b = 0; b < numBlocks; ++b) {
            int start = b * blockSize_;
            int end = std::min(start + blockSize_, n);
            int size = end - start;
            
            blockStarts_[b] = start;
            
            // Extract dense block
            std::vector<std::vector<T>> block(size, std::vector<T>(size, T(0)));
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    block[i][j] = A(start + i, start + j);
                }
            }
            
            // Compute inverse (simple Gaussian elimination for small blocks)
            blockInverses_[b] = invertSmallMatrix(block);
        }
        blockStarts_[numBlocks] = n;
    }
    
    void apply(const std::vector<T>& r, std::vector<T>& z) const override {
        if (blockInverses_.empty()) {
            throw std::logic_error("BlockJacobiPreconditioner::setup() must be called before apply()");
        }
        if (blockStarts_.back() != static_cast<int>(r.size())) {
            throw std::invalid_argument("BlockJacobiPreconditioner input size does not match setup matrix");
        }
        z.resize(r.size());
        
        int numBlocks = static_cast<int>(blockInverses_.size());
        for (int b = 0; b < numBlocks; ++b) {
            int start = blockStarts_[b];
            int end = blockStarts_[b + 1];
            int size = end - start;
            
            // z[start:end] = blockInverse[b] * r[start:end]
            for (int i = 0; i < size; ++i) {
                T sum = T(0);
                for (int j = 0; j < size; ++j) {
                    sum += blockInverses_[b][i][j] * r[start + j];
                }
                z[start + i] = sum;
            }
        }
    }
    
private:
    // Simple matrix inversion for small dense matrices
    static std::vector<std::vector<T>> invertSmallMatrix(std::vector<std::vector<T>> A) {
        int n = static_cast<int>(A.size());
        std::vector<std::vector<T>> inv(n, std::vector<T>(n, T(0)));
        
        // Initialize inverse as identity
        for (int i = 0; i < n; ++i) inv[i][i] = T(1);
        
        // Gaussian elimination with partial pivoting
        for (int k = 0; k < n; ++k) {
            // Find pivot
            int maxRow = k;
            for (int i = k + 1; i < n; ++i) {
                if (std::abs(A[i][k]) > std::abs(A[maxRow][k])) {
                    maxRow = i;
                }
            }
            std::swap(A[k], A[maxRow]);
            std::swap(inv[k], inv[maxRow]);
            
            // Check for singular
            if (std::abs(A[k][k]) < std::numeric_limits<T>::epsilon() * 100) {
                // Fallback to identity for this block
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        inv[i][j] = (i == j) ? T(1) : T(0);
                    }
                }
                return inv;
            }
            
            // Scale pivot row
            T pivot = A[k][k];
            for (int j = 0; j < n; ++j) {
                A[k][j] /= pivot;
                inv[k][j] /= pivot;
            }
            
            // Eliminate column
            for (int i = 0; i < n; ++i) {
                if (i != k) {
                    T factor = A[i][k];
                    for (int j = 0; j < n; ++j) {
                        A[i][j] -= factor * A[k][j];
                        inv[i][j] -= factor * inv[k][j];
                    }
                }
            }
        }
        
        return inv;
    }
};

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_PRECONDITIONERS_H
