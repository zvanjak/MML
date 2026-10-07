///////////////////////////////////////////////////////////////////////////////////////////
// GMRES.h - Generalized Minimal Residual solver for general sparse systems
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// GMRES solves Ax = b for general sparse matrices by minimizing the residual
// norm over a Krylov subspace. It's more robust than BiCGSTAB for difficult
// problems but requires more memory.
//
// Algorithm: Restarted GMRES(m) with modified Gram-Schmidt
// 1. r₀ = M⁻¹(b - Ax₀)
// 2. v₁ = r₀ / ||r₀||
// 3. For j = 1, ..., m:
//    a. w = M⁻¹ * A * vⱼ
//    b. For i = 1, ..., j:
//       hᵢⱼ = (w, vᵢ)
//       w = w - hᵢⱼ * vᵢ
//    c. hⱼ₊₁,ⱼ = ||w||
//    d. vⱼ₊₁ = w / hⱼ₊₁,ⱼ
//    e. Apply Givens rotations to H
//    f. Check convergence
// 4. Solve Hy = g (upper triangular)
// 5. x = x₀ + V * y
// 6. If not converged, restart with new x₀ = x
//
// Memory: O(n * m) for basis vectors
// Convergence: Guaranteed to converge in at most n iterations (without restart)
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_CORE_SPARSE_SOLVERS_GMRES_H
#define MML_CORE_SPARSE_SOLVERS_GMRES_H

#include "IterativeSolverBase.h"
#include <cmath>

namespace MML::SparseSolvers {

template<typename T>
class GMRES : public IterativeSolverBase<T> {
    using Base = IterativeSolverBase<T>;
    using Base::config_;
    using Base::callback_;
    using Base::preconditioner_;
    
    int restart_ = 30;  // Restart parameter (Krylov subspace dimension)
    
public:
    GMRES() = default;
    explicit GMRES(int restart) : restart_(restart) {}
    
    void setRestart(int m) { restart_ = m; }
    int restart() const { return restart_; }
    
    // Main solve: Ax = b with initial guess x
    // Uses left preconditioning: solve M⁻¹Ax = M⁻¹b
    SolverResult<T> solve(
        const SparseMatrixCSR<T>& A,
        const std::vector<T>& b,
        std::vector<T>& x) override 
    {
        SolverResult<T> result;
        
        // Validate input
        if (!this->validateInput(A, b, x)) {
            result.status = SolverStatus::InvalidInput;
            result.message = "Matrix and vector dimensions don't match";
            return result;
        }
        
        const int n = A.rows();
        const int m = std::min(restart_, n);  // Krylov subspace dimension
        const T bnorm = norm2(b);
        
        if (config_.verbose) {
            std::cout << "GMRES: n=" << n << ", restart=" << restart_ << ", m=" << m << "\n";
        }
        
        // Handle zero RHS
        if (bnorm < config_.absTolerance) {
            std::fill(x.begin(), x.end(), T(0));
            result.status = SolverStatus::Success;
            result.iterations = 0;
            result.residualNorm = T(0);
            result.relativeResidual = T(0);
            result.message = "Zero RHS";
            return result;
        }
        
        // Allocate Arnoldi vectors: V[0:m+1] vectors of size n
        std::vector<std::vector<T>> V(m + 1, std::vector<T>(n));
        
        // Hessenberg matrix H[0:m+1, 0:m]
        std::vector<std::vector<T>> H(m + 1, std::vector<T>(m, T(0)));
        
        // Givens rotation parameters
        std::vector<T> cs(m);  // Cosines
        std::vector<T> sn(m);  // Sines
        
        // RHS for least squares problem
        std::vector<T> g(m + 1);
        
        // Solution in Krylov basis
        std::vector<T> y(m);
        
        // Work vectors
        std::vector<T> w(n);
        std::vector<T> tmp(n);

        preconditioner_->apply(b, tmp);
        const T preconditionedBnorm = norm2(tmp);
        
        int totalIter = 0;
        T rnorm = T(0);
        
        // Outer restart loop
        for (int restart = 0; restart < config_.maxIterations / m + 1; ++restart) {
            // Compute initial residual: r = b - Ax
            std::vector<T> r = residual(A, x, b);
            rnorm = norm2(r);
            
            this->reportProgress(totalIter, rnorm);
            
            // Check if already converged
            if (this->checkConvergence(rnorm, bnorm)) {
                result.status = SolverStatus::Success;
                result.iterations = totalIter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "Converged";
                return result;
            }

            preconditioner_->apply(r, tmp);
            const T preconditionedRnorm = norm2(tmp);
            if (preconditionedRnorm == T(0)) {
                result.status = SolverStatus::Breakdown;
                result.iterations = totalIter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "GMRES breakdown: preconditioned residual is zero";
                return result;
            }
            
            // v₁ = M⁻¹r / ||M⁻¹r||
            for (int i = 0; i < n; ++i) {
                V[0][i] = tmp[i] / preconditionedRnorm;
            }
            
            // Initialize g = [||M⁻¹r||, 0, 0, ...]
            std::fill(g.begin(), g.end(), T(0));
            g[0] = preconditionedRnorm;
            
            // Reset Hessenberg matrix
            for (auto& row : H) {
                std::fill(row.begin(), row.end(), T(0));
            }
            
            int j;
            for (j = 0; j < m && totalIter < config_.maxIterations; ++j) {
                ++totalIter;
                
                if (config_.verbose) {
                    std::cout << "  j=" << j << ", m=" << m << "\n";
                }
                
                // Arnoldi step for the left-preconditioned operator M⁻¹A
                tmp = A * V[j];
                preconditioner_->apply(tmp, w);
                
                // Modified Gram-Schmidt orthogonalization
                for (int i = 0; i <= j; ++i) {
                    H[i][j] = dot(w, V[i]);
                    axpy(-H[i][j], V[i], w);
                }
                
                H[j+1][j] = norm2(w);
                
                // Debug: Check orthogonality
                if (config_.verbose) {
                    T ortho = T(0);
                    for (int i = 0; i <= j; ++i) {
                        ortho += std::abs(dot(w, V[i]));
                    }
                    std::cout << "  Orthogonality check: " << ortho << "\n";
                    std::cout << "  H[" << j+1 << "][" << j << "] = " << H[j+1][j] << "\n";
                }
                
                // Check for breakdown (H[j+1][j] ≈ 0 means Krylov subspace exhausted)
                // This is "lucky breakdown" - exact solution lies in the current subspace
                const T breakdown_tol = std::numeric_limits<T>::epsilon() * 1000;
                if (H[j+1][j] < breakdown_tol) {
                    if (config_.verbose) {
                        std::cout << "  Breakdown: H[" << j+1 << "][" << j << "] = " 
                                  << H[j+1][j] << " < " << breakdown_tol << "\n";
                    }
                    
                    // Apply previous Givens rotations to column j of H
                    for (int i = 0; i < j; ++i) {
                        T temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
                        H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
                        H[i][j] = temp;
                    }
                    
                    // For breakdown: H[j+1][j] is ~0, so Givens rotation for row j
                    // is essentially identity (c=1, s=0 if H[j][j] > 0)
                    // Just ensure H[j+1][j] is exactly 0
                    H[j+1][j] = T(0);
                    
                    // Compute Givens for potentially small diagonal
                    // This handles case where H[j][j] might also be small
                    if (std::abs(H[j][j]) > breakdown_tol) {
                        cs[j] = T(1);
                        sn[j] = T(0);
                        // g[j+1] stays 0, no update needed
                    } else {
                        // Both H[j][j] and H[j+1][j] are near zero - degenerate case
                        cs[j] = T(1);
                        sn[j] = T(0);
                    }
                    
                    if (config_.verbose) {
                        std::cout << "  Lucky breakdown: including iteration j=" << j << "\n";
                    }
                    ++j;  // Include this iteration
                    break;
                }
                
                // v[j+1] = w / H[j+1][j]
                for (int i = 0; i < n; ++i) {
                    V[j+1][i] = w[i] / H[j+1][j];
                }
                
                // Apply previous Givens rotations to column j of H
                for (int i = 0; i < j; ++i) {
                    T temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
                    H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
                    H[i][j] = temp;
                }
                
                // Compute Givens rotation for row j
                generateGivens(H[j][j], H[j+1][j], cs[j], sn[j]);
                
                // Apply Givens rotation to column j
                H[j][j] = cs[j] * H[j][j] + sn[j] * H[j+1][j];
                H[j+1][j] = T(0);
                
                // Apply Givens rotation to g
                T temp = cs[j] * g[j] + sn[j] * g[j+1];
                g[j+1] = -sn[j] * g[j] + cs[j] * g[j+1];
                g[j] = temp;
                
                // Residual estimate
                rnorm = std::abs(g[j+1]);
                
                if (config_.verbose) {
                    std::cout << "  After Givens: g[" << j+1 << "] = " << g[j+1] << "\n";
                }
                
                this->reportProgress(totalIter, rnorm);
                
                // Check convergence
                if (this->checkConvergence(rnorm, preconditionedBnorm)) {
                    if (config_.verbose) {
                        std::cout << "  CONVERGED at j=" << j << ", rnorm=" << rnorm << "\n";
                    }
                    ++j;  // Include this iteration in solution
                    break;
                }
            }
            
            if (config_.verbose) {
                std::cout << "  Exited inner loop: j=" << j << ", totalIter=" << totalIter << "\n";
            }
            
            // Solve upper triangular system Hy = g
            if (config_.verbose) {
                std::cout << "  H diagonal: [";
                for (int i = 0; i < j; ++i) {
                    std::cout << H[i][i] << (i < j-1 ? ", " : "");
                }
                std::cout << "]\n";
                std::cout << "  g = [";
                for (int i = 0; i < j; ++i) {
                    std::cout << g[i] << (i < j-1 ? ", " : "");
                }
                std::cout << "]\n";
            }
            
            for (int i = j - 1; i >= 0; --i) {
                y[i] = g[i];
                for (int k = i + 1; k < j; ++k) {
                    y[i] -= H[i][k] * y[k];
                }
                if (std::abs(H[i][i]) > std::numeric_limits<T>::epsilon() * 100) {
                    y[i] /= H[i][i];
                }
            }
            
            if (config_.verbose) {
                std::cout << "  y = [";
                for (int i = 0; i < j; ++i) {
                    std::cout << y[i] << (i < j-1 ? ", " : "");
                }
                std::cout << "]\n";
            }
            
            // Update solution: x = x + V * y
            for (int i = 0; i < j; ++i) {
                axpy(y[i], V[i], x);
            }
            
            // Check convergence (true residual)
            r = residual(A, x, b);
            rnorm = norm2(r);
            
            if (this->checkConvergence(rnorm, bnorm)) {
                result.status = SolverStatus::Success;
                result.iterations = totalIter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "Converged";
                return result;
            }
        }
        
        // Did not converge
        result.status = SolverStatus::MaxIterations;
        result.iterations = totalIter;
        result.residualNorm = rnorm;
        result.relativeResidual = rnorm / bnorm;
        result.message = "Maximum iterations reached";
        return result;
    }
    
private:
    // Compute Givens rotation coefficients to zero out b
    // After rotation: [c  s; -s c] * [a; b] = [r; 0] where r = sqrt(a^2 + b^2) >= 0
    static void generateGivens(T a, T b, T& c, T& s) {
        if (b == T(0)) {
            c = (a >= T(0)) ? T(1) : T(-1);
            s = T(0);
        } else if (a == T(0)) {
            c = T(0);
            s = (b >= T(0)) ? T(1) : T(-1);
        } else if (std::abs(b) > std::abs(a)) {
            T tau = a / b;
            T temp = std::sqrt(T(1) + tau * tau);
            s = T(1) / temp;
            if (b < T(0)) s = -s;
            c = s * tau;
        } else {
            T tau = b / a;
            T temp = std::sqrt(T(1) + tau * tau);
            c = T(1) / temp;
            if (a < T(0)) c = -c;
            s = c * tau;
        }
    }
};

// Convenience function
template<typename T>
SolverResult<T> solveGMRES(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    int restart = 30,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    GMRES<T> solver(restart);
    solver.setConfig(config);
    return solver.solve(A, b, x);
}

template<typename T, typename PrecType>
SolverResult<T> solveGMRES(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    std::shared_ptr<PrecType> prec,
    int restart = 30,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    GMRES<T> solver(restart);
    solver.setConfig(config);
    solver.setPreconditioner(prec);
    return solver.solve(A, b, x);
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_GMRES_H
