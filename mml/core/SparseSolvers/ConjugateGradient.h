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
#include <algorithm>

namespace MML::SparseSolvers {

template<typename T>
class ConjugateGradient : public IterativeSolverBase<T> {
    using Base = IterativeSolverBase<T>;
    using Base::config_;
    using Base::callback_;
    using Base::preconditioner_;
    
public:
    ConjugateGradient() = default;

	protected:
		const char* algorithmName() const noexcept override { return "ConjugateGradient"; }
		KrylovResult<T> solveImpl(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x) override {
			KrylovResult<T> result;

			// Validate input
			if (!this->validateInput(A, b, x)) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "Matrix and vector dimensions don't match";
				return result;
			}

			if (config_.validateSPD && !isPositiveDefinite(A)) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "CG requires a symmetric positive-definite matrix";
				return result;
			}

			const int n = A.rows();
			const T bnorm = norm2(b);

			// Handle zero RHS
			if (bnorm == T(0) || bnorm < config_.absTolerance) {
				std::fill(x.begin(), x.end(), T(0));
				result.status = AlgorithmStatus::Success;
				result.iterations_used = 0;
				result.residualNorm = T(0);
				result.relativeResidual = T(0);
				result.error_message = "Zero RHS";
				return result;
			}

			// Allocate work vectors
			Vector<T> r(n);	 // Residual
			Vector<T> z(n);	 // Preconditioned residual
			Vector<T> p(n);	 // Search direction
			Vector<T> Ap(n); // Matrix-vector product

			// Initial residual: r = b - Ax
			r = residual(A, x, b);
			T rnorm = norm2(r);

			this->reportProgress(0, rnorm);

			// Check if already converged
			if (this->checkConvergence(rnorm, bnorm)) {
				result.status = AlgorithmStatus::Success;
				result.iterations_used = 0;
				result.residualNorm = rnorm;
				result.relativeResidual = rnorm / bnorm;
				result.error_message = "Converged at initial guess";
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

				// Breakdown check. For SPD A, (p, Ap) > 0 whenever p != 0, and
				// (p, Ap) / (||p|| ||Ap||) >= 2 sqrt(k)/(k+1) (Kantorovich), so a value that is
				// non-positive or negligible *relative* to ||p|| ||Ap|| signals loss of
				// definiteness. An absolute threshold is wrong here: (p, Ap) is quadratic in
				// the residual and legitimately drops below any fixed epsilon near convergence.
				const T pApScale = norm2(p) * norm2(Ap);
				if (!(pAp > std::numeric_limits<T>::epsilon() * pApScale)) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "CG breakdown: (p, Ap) ≈ 0";
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
					result.status = AlgorithmStatus::Success;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "Converged";
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
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.iterations_used = config_.maxIterations;
			result.residualNorm = rnorm;
			result.relativeResidual = rnorm / bnorm;
			result.error_message = "Maximum iterations reached";
			return result;
		}

private:
    static bool isPositiveDefinite(const SparseMatrixCSR<T>& A) {
        T matrixScale = T(1);
        for (const T value : A.values()) {
            matrixScale = std::max(matrixScale, std::abs(value));
        }

        const T symmetryTolerance = std::sqrt(std::numeric_limits<T>::epsilon()) * matrixScale;
        if (!A.IsSymmetric(symmetryTolerance)) {
            return false;
        }

        const int n = A.rows();
        const T pivotTolerance = std::numeric_limits<T>::epsilon() * matrixScale * T(100);
				std::vector<Vector<T>> lower(n, Vector<T>(n, T(0)));

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
KrylovResult<T> solveCG(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x, const SolverConfig<T>& config = SolverConfig<T>()) {
	ConjugateGradient<T> cg;
	cg.setConfig(config);
	return cg.solve(A, b, x);
}

// Convenience function with preconditioner
template<typename T, typename PrecType>
KrylovResult<T> solveCG(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x, std::shared_ptr<PrecType> prec,
												const SolverConfig<T>& config = SolverConfig<T>()) {
	ConjugateGradient<T> cg;
	cg.setConfig(config);
	cg.setPreconditioner(prec);
	return cg.solve(A, b, x);
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_CONJUGATE_GRADIENT_H
