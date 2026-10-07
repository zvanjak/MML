///////////////////////////////////////////////////////////////////////////////////////////
// ConjugateGradient.h - Conjugate Gradient solver for symmetric positive definite systems
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// The Conjugate Gradient (CG) method solves Ax = b where A is symmetric
// positive definite. It's the method of choice for discrete Laplacians
// and other SPD systems arising from PDEs.
//
// Algorithm: Standard preconditioned CG
// 1. r₀ = b - Ax₀
// 2. z₀ = M⁻¹r₀
// 3. p₀ = z₀
// 4. For k = 0, 1, 2, ...
//    a. αₖ = (rₖ, zₖ) / (pₖ, Apₖ)
//    b. xₖ₊₁ = xₖ + αₖpₖ
//    c. rₖ₊₁ = rₖ - αₖApₖ
//    d. Check convergence
//    e. zₖ₊₁ = M⁻¹rₖ₊₁
//    f. βₖ = (rₖ₊₁, zₖ₊₁) / (rₖ, zₖ)
//    g. pₖ₊₁ = zₖ₊₁ + βₖpₖ
//
// Convergence: O(√κ) iterations where κ = cond(A)
// With preconditioner: O(√κ(M⁻¹A))
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_CORE_SPARSE_SOLVERS_CONJUGATE_GRADIENT_H
#define MML_CORE_SPARSE_SOLVERS_CONJUGATE_GRADIENT_H

#include "IterativeSolverBase.h"

namespace MML::SparseSolvers {

template<typename T>
class ConjugateGradient : public IterativeSolverBase<T> {
    using Base = IterativeSolverBase<T>;
    using Base::config_;
    using Base::callback_;
    using Base::preconditioner_;
    
public:
    ConjugateGradient() = default;
    
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

        if (config_.validateSPD && !isPositiveDefinite(A)) {
            result.status = SolverStatus::InvalidInput;
            result.message = "CG requires a symmetric positive-definite matrix";
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
        std::vector<T> r(n);   // Residual
        std::vector<T> z(n);   // Preconditioned residual
        std::vector<T> p(n);   // Search direction
        std::vector<T> Ap(n);  // Matrix-vector product
        
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
        
        // z = M⁻¹r
        preconditioner_->apply(r, z);
        
        // p = z
        p = z;
        
        // rz = (r, z) = (r, M⁻¹r)
        T rz = dot(r, z);
        
        // Main CG iteration
        for (int iter = 1; iter <= config_.maxIterations; ++iter) {
            // Ap = A * p
            Ap = A * p;
            
            // α = (r, z) / (p, Ap)
            T pAp = dot(p, Ap);
            
            // Check for breakdown
            if (std::abs(pAp) < std::numeric_limits<T>::epsilon() * 100) {
                result.status = SolverStatus::Breakdown;
                result.iterations = iter;
                result.residualNorm = rnorm;
                result.relativeResidual = rnorm / bnorm;
                result.message = "CG breakdown: (p, Ap) ≈ 0";
                return result;
            }
            
            T alpha = rz / pAp;
            
            // x = x + α * p
            axpy(alpha, p, x);
            
            // r = r - α * Ap
            axpy(-alpha, Ap, r);
            
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
            
            // z = M⁻¹r
            preconditioner_->apply(r, z);
            
            // β = (r_new, z_new) / (r_old, z_old)
            T rz_new = dot(r, z);
            T beta = rz_new / rz;
            rz = rz_new;
            
            // p = z + β * p
            for (int i = 0; i < n; ++i) {
                p[i] = z[i] + beta * p[i];
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

private:
    static bool isPositiveDefinite(const SparseMatrixCSR<T>& A) {
        T matrixScale = T(1);
        for (const T value : A.values()) {
            matrixScale = std::max(matrixScale, std::abs(value));
        }

        const T symmetryTolerance = std::sqrt(std::numeric_limits<T>::epsilon()) * matrixScale;
        if (!A.isSymmetric(symmetryTolerance)) {
            return false;
        }

        const int n = A.rows();
        const T pivotTolerance = std::numeric_limits<T>::epsilon() * matrixScale * T(100);
        std::vector<std::vector<T>> lower(n, std::vector<T>(n, T(0)));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                T value = A(i, j);
                for (int k = 0; k < j; ++k) {
                    value -= lower[i][k] * lower[j][k];
                }

                if (i == j) {
                    if (!(value > pivotTolerance)) {
                        return false;
                    }
                    lower[i][j] = std::sqrt(value);
                } else {
                    lower[i][j] = value / lower[j][j];
                }
            }
        }

        return true;
    }
};

// Convenience function
template<typename T>
SolverResult<T> solveCG(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    ConjugateGradient<T> cg;
    cg.setConfig(config);
    return cg.solve(A, b, x);
}

// Convenience function with preconditioner
template<typename T, typename PrecType>
SolverResult<T> solveCG(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    std::shared_ptr<PrecType> prec,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    ConjugateGradient<T> cg;
    cg.setConfig(config);
    cg.setPreconditioner(prec);
    return cg.solve(A, b, x);
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_CONJUGATE_GRADIENT_H
