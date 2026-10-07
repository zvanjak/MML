///////////////////////////////////////////////////////////////////////////////////////////
// BiCGSTAB.h - Biconjugate Gradient Stabilized solver for general sparse systems
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// BiCGSTAB solves Ax = b for general (non-symmetric) sparse matrices.
// It's particularly useful for convection-diffusion problems and other
// non-symmetric PDEs.
//
// Algorithm: Van der Vorst (1992)
// 1. r₀ = b - Ax₀
// 2. Choose r̂₀ (usually r̂₀ = r₀)
// 3. For k = 0, 1, 2, ...
//    a. ρₖ = (r̂₀, rₖ)
//    b. β = (ρₖ/ρₖ₋₁)(αₖ₋₁/ωₖ₋₁)
//    c. pₖ = rₖ + β(pₖ₋₁ - ωₖ₋₁vₖ₋₁)
//    d. p̂ = M⁻¹pₖ
//    e. v = Ap̂
//    f. α = ρₖ/(r̂₀, v)
//    g. s = rₖ - αv
//    h. ŝ = M⁻¹s
//    i. t = Aŝ
//    j. ω = (t, s)/(t, t)
//    k. xₖ₊₁ = xₖ + αp̂ + ωŝ
//    l. rₖ₊₁ = s - ωt
//
// Convergence: Generally converges in O(√κ) iterations for well-conditioned systems
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_CORE_SPARSE_SOLVERS_BICGSTAB_H
#define MML_CORE_SPARSE_SOLVERS_BICGSTAB_H

#include "IterativeSolverBase.h"

namespace MML::SparseSolvers {

template<typename T>
class BiCGSTAB : public IterativeSolverBase<T> {
    using Base = IterativeSolverBase<T>;
    using Base::config_;
    using Base::callback_;
    using Base::preconditioner_;
    
public:
    BiCGSTAB() = default;
    
    // Main solve: Ax = b with initial guess x
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
        const T bnorm = norm2(b);
        
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
        
        // Allocate work vectors
        std::vector<T> r(n);      // Residual
        std::vector<T> r0(n);     // Shadow residual (fixed)
        std::vector<T> p(n);      // Search direction
        std::vector<T> v(n);      // A * p_hat
        std::vector<T> s(n);      // s = r - α*v
        std::vector<T> t(n);      // t = A * s_hat
        std::vector<T> p_hat(n);  // Preconditioned p
        std::vector<T> s_hat(n);  // Preconditioned s
        
        // Initial residual: r = b - Ax
        r = residual(A, x, b);
        T rnorm = norm2(r);
        
        this->reportProgress(0, rnorm);
        
        // Check if already converged
        if (this->checkConvergence(rnorm, bnorm)) {
            result.status = SolverStatus::Success;
            result.iterations = 0;
            result.residualNorm = rnorm;
            result.relativeResidual = rnorm / bnorm;
            result.message = "Converged at initial guess";
            return result;
        }
        
        // Shadow residual: r̂₀ = r₀
        r0 = r;
        
        T rho = T(1);
        T alpha = T(1);
        T omega = T(1);
        
        std::fill(p.begin(), p.end(), T(0));
        std::fill(v.begin(), v.end(), T(0));
        
        // Main BiCGSTAB iteration
        for (int iter = 1; iter <= config_.maxIterations; ++iter) {
            // ρ = (r̂₀, r)
            T rho_new = dot(r0, r);
            
            // Check for breakdown
            if (std::abs(rho_new) < std::numeric_limits<T>::epsilon() * 100) {
                result.status = SolverStatus::Breakdown;
                result.iterations = iter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "BiCGSTAB breakdown: rho ≈ 0";
                return result;
            }
            
            // β = (ρ_new / ρ_old) * (α / ω)
            T beta = (rho_new / rho) * (alpha / omega);
            rho = rho_new;
            
            // p = r + β * (p - ω * v)
            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
            
            // p̂ = M⁻¹p
            preconditioner_->apply(p, p_hat);
            
            // v = A * p̂
            v = A * p_hat;
            
            // α = ρ / (r̂₀, v)
            T r0v = dot(r0, v);
            
            // Check for breakdown
            if (std::abs(r0v) < std::numeric_limits<T>::epsilon() * 100) {
                result.status = SolverStatus::Breakdown;
                result.iterations = iter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "BiCGSTAB breakdown: (r0, v) ≈ 0";
                return result;
            }
            
            alpha = rho / r0v;
            
            // s = r - α * v
            for (int i = 0; i < n; ++i) {
                s[i] = r[i] - alpha * v[i];
            }
            
            // Check if s is small enough (early termination)
            T snorm = norm2(s);
            if (this->checkConvergence(snorm, bnorm)) {
                // x = x + α * p̂
                axpy(alpha, p_hat, x);
                
                result.status = SolverStatus::Success;
                result.iterations = iter;
                result.residualNorm = snorm;
                result.relativeResidual = snorm / bnorm;
                result.message = "Converged (early)";
                return result;
            }
            
            // ŝ = M⁻¹s
            preconditioner_->apply(s, s_hat);
            
            // t = A * ŝ
            t = A * s_hat;
            
            // ω = (t, s) / (t, t)
            T tt = dot(t, t);
            
            // Check for breakdown
            if (std::abs(tt) < std::numeric_limits<T>::epsilon() * 100) {
                result.status = SolverStatus::Breakdown;
                result.iterations = iter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "BiCGSTAB breakdown: (t, t) ≈ 0";
                return result;
            }
            
            omega = dot(t, s) / tt;
            
            // x = x + α * p̂ + ω * ŝ
            axpy(alpha, p_hat, x);
            axpy(omega, s_hat, x);
            
            // r = s - ω * t
            for (int i = 0; i < n; ++i) {
                r[i] = s[i] - omega * t[i];
            }
            
            rnorm = norm2(r);
            this->reportProgress(iter, rnorm);
            
            // Check convergence
            if (this->checkConvergence(rnorm, bnorm)) {
                result.status = SolverStatus::Success;
                result.iterations = iter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "Converged";
                return result;
            }
            
            // Check for stagnation
            if (std::abs(omega) < std::numeric_limits<T>::epsilon() * 100) {
                result.status = SolverStatus::Stagnation;
                result.iterations = iter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "BiCGSTAB stagnation: omega ≈ 0";
                return result;
            }
        }
        
        // Did not converge
        result.status = SolverStatus::MaxIterations;
        result.iterations = config_.maxIterations;
        result.residualNorm = rnorm;
        result.relativeResidual = rnorm / bnorm;
        result.message = "Maximum iterations reached";
        return result;
    }
};

// Convenience function
template<typename T>
SolverResult<T> solveBiCGSTAB(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    BiCGSTAB<T> solver;
    solver.setConfig(config);
    return solver.solve(A, b, x);
}

// Convenience function with preconditioner
template<typename T, typename PrecType>
SolverResult<T> solveBiCGSTAB(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    std::shared_ptr<PrecType> prec,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    BiCGSTAB<T> solver;
    solver.setConfig(config);
    solver.setPreconditioner(prec);
    return solver.solve(A, b, x);
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_BICGSTAB_H
