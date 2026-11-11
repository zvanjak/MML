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
#include <algorithm>

namespace MML::SparseSolvers {

template<MMLReal T>
class BiCGSTAB : public IterativeSolverBase<T> {
    using Base = IterativeSolverBase<T>;
    using Base::config_;
    using Base::callback_;
    using Base::preconditioner_;
    
public:
    BiCGSTAB() = default;

	protected:
		const char* algorithmName() const noexcept override { return "BiCGSTAB"; }
		KrylovResult<T> solveImpl(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x) override {
			KrylovResult<T> result;

			// Validate input
			if (!this->validateInput(A, b, x)) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "Matrix and vector dimensions don't match";
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
			Vector<T> r(n);			// Residual
			Vector<T> r0(n);		// Shadow residual (fixed)
			Vector<T> p(n);			// Search direction
			Vector<T> v(n);			// A * p_hat
			Vector<T> s(n);			// s = r - α*v
			Vector<T> t(n);			// t = A * s_hat
			Vector<T> p_hat(n); // Preconditioned p
			Vector<T> s_hat(n); // Preconditioned s

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

			// Shadow residual: r̂₀ = r₀
			r0 = r;
			T r0norm = rnorm;

			// Breakdown tests are *relative*: ρ = (r̂₀, r), (r̂₀, v) and (t, t) are bilinear in
			// quantities that shrink with the residual, so comparing them with a fixed epsilon
			// falsely reports breakdown just before convergence on any non-toy system.
			const T eps = std::numeric_limits<T>::epsilon();

			T rho = T(1);
			T alpha = T(1);
			T omega = T(1);
			bool restartedLastStep = false;

			std::fill(p.begin(), p.end(), T(0));
			std::fill(v.begin(), v.end(), T(0));
			const auto restartFromTrueResidual = [&]() {
				r = residual(A, x, b);
				rnorm = norm2(r);
				r0 = r;
				r0norm = rnorm;
				rho = alpha = omega = T(1);
				restartedLastStep = false;
				std::fill(p.begin(), p.end(), T(0));
				std::fill(v.begin(), v.end(), T(0));
			};

			// Main BiCGSTAB iteration
			for (int iter = 1; iter <= config_.maxIterations; ++iter) {
				// ρ = (r̂₀, r)
				T rho_new = dot(r0, r);

				// Near-breakdown: r has become (numerically) orthogonal to the arbitrary shadow
				// residual r̂₀. As in Eigen/PETSc, restart the Lanczos process with r̂₀ = r instead
				// of giving up; only a second breakdown immediately after a restart is fatal.
				if (std::abs(rho_new) < eps * r0norm * rnorm) {
					if (restartedLastStep || !(rnorm > T(0))) {
						result.status = AlgorithmStatus::NumericalInstability;
						result.iterations_used = iter;
						result.residualNorm = rnorm;
						result.relativeResidual = rnorm / bnorm;
						result.error_message = "BiCGSTAB breakdown: rho ≈ 0";
						return result;
					}
					r = residual(A, x, b);
					rnorm = norm2(r);
					r0 = r;
					r0norm = rnorm;
					rho_new = rnorm * rnorm;
					p = r;
					std::fill(v.begin(), v.end(), T(0));
					alpha = T(1);
					omega = T(1);
					restartedLastStep = true;
				} else {
					restartedLastStep = false;
					// β = (ρ_new / ρ_old) * (α / ω)
					T beta = (rho_new / rho) * (alpha / omega);

					// p = r + β * (p - ω * v)
					for (int i = 0; i < n; ++i) {
						p[i] = r[i] + beta * (p[i] - omega * v[i]);
					}
				}
				rho = rho_new;

				// p̂ = M⁻¹p
				preconditioner_->apply(p, p_hat);

				// v = A * p̂
				v = A * p_hat;

				// α = ρ / (r̂₀, v)
				T r0v = dot(r0, v);

				// Breakdown: v = A p̂ nearly orthogonal to the shadow residual
				if (std::abs(r0v) < eps * r0norm * norm2(v)) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "BiCGSTAB breakdown: (r0, v) ≈ 0";
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
					restartFromTrueResidual();
					this->reportProgress(iter, rnorm);
					if (!this->checkConvergence(rnorm, bnorm))
						continue;

					result.status = AlgorithmStatus::Success;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "Converged (early)";
					return result;
				}

				// ŝ = M⁻¹s
				preconditioner_->apply(s, s_hat);

				// t = A * ŝ
				t = A * s_hat;

				// ω = (t, s) / (t, t)
				T tt = dot(t, t);

				// Breakdown: t = A ŝ vanished although s did not (singular operator)
				if (!(tt > T(0)) || !std::isfinite(tt)) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "BiCGSTAB breakdown: (t, t) ≈ 0";
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
					restartFromTrueResidual();
					if (!this->checkConvergence(rnorm, bnorm))
						continue;
					result.status = AlgorithmStatus::Success;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "Converged";
					return result;
				}

				// Stagnation: the stabilisation step ω t changed s by less than rounding, and
				// the next β = (ρ_new/ρ)(α/ω) would blow up. Scale-aware: compares ||ω t|| to ||s||.
				if (std::abs(omega) * std::sqrt(tt) <= eps * snorm) {
					result.status = AlgorithmStatus::Stalled;
					result.iterations_used = iter;
					result.residualNorm = rnorm;
					result.relativeResidual = rnorm / bnorm;
					result.error_message = "BiCGSTAB stagnation: omega ≈ 0";
					return result;
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
};

// Convenience function
template<MMLReal T>
KrylovResult<T> solveBiCGSTAB(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x,
															const SolverConfig<T>& config = SolverConfig<T>()) {
	BiCGSTAB<T> solver;
	solver.setConfig(config);
	return solver.solve(A, b, x);
}

// Convenience function with preconditioner
template<MMLReal T, typename PrecType>
KrylovResult<T> solveBiCGSTAB(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x, std::shared_ptr<PrecType> prec,
															const SolverConfig<T>& config = SolverConfig<T>()) {
	BiCGSTAB<T> solver;
	solver.setConfig(config);
	solver.setPreconditioner(prec);
	return solver.solve(A, b, x);
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_BICGSTAB_H
