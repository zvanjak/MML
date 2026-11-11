///////////////////////////////////////////////////////////////////////////////////////////
// IterativeSolverBase.h - Base class for iterative linear solvers
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// This file provides the base infrastructure for iterative linear solvers
// used in finite difference PDE solutions.
//
// Features:
// - Common solver parameters (tolerance, max iterations)
// - Convergence tracking and result structure
// - Callback support for progress monitoring
// - Preconditioner interface
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_CORE_SPARSE_SOLVER_BASE_H
#define MML_CORE_SPARSE_SOLVER_BASE_H

#include "../../base/SparseMatrix/SparseMatrix.h"
#include <mml/base/AlgorithmTypes.h>

#include <functional>
#include <string>
#include <memory>

namespace MML::SparseSolvers {

using namespace ::MML::SparseMatrix;

//=============================================================================
// Solver Configuration
//=============================================================================

template<typename T>
struct SolverConfig {
    T tolerance = PrecisionValues<T>::IterativeSolverTolerance;   // Relative residual tolerance
    T absTolerance = PrecisionValues<T>::RoundoffZeroThreshold;  // Absolute residual floor (roundoff at unit scale)
    int maxIterations = 1000;               // Maximum iterations
    bool verbose = false;                   // Print progress
    bool validateSPD = false;               // Validate CG matrix with a diagnostic Cholesky factorization
    
    // Convergence: ||r|| < tol * ||b|| + absTol
    
    SolverConfig& setTolerance(T tol) { tolerance = tol; return *this; }
    SolverConfig& setAbsTolerance(T atol) { absTolerance = atol; return *this; }
    SolverConfig& setMaxIterations(int n) { maxIterations = n; return *this; }
    SolverConfig& setVerbose(bool v) { verbose = v; return *this; }
    SolverConfig& setValidateSPD(bool validate) { validateSPD = validate; return *this; }
};

//=============================================================================
// Solver Result
//=============================================================================

using SolverStatus [[deprecated("use AlgorithmStatus")]] = AlgorithmStatus;

template<typename T>
struct [[nodiscard]] KrylovResult : IterativeResultBase {
	KrylovResult() {
		status = AlgorithmStatus::InvalidInput;
		error_message = "No Krylov solve performed";
	}
	T residualNorm = std::numeric_limits<T>::quiet_NaN();
	T relativeResidual = std::numeric_limits<T>::quiet_NaN();
};

template<typename T>
using SolverResult [[deprecated("use KrylovResult<T>")]] = KrylovResult<T>;

//=============================================================================
// Preconditioner Interface
//=============================================================================

// Abstract preconditioner: solves M*z = r approximately
template<typename T>
class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    
    // Apply preconditioner: z = M^{-1} * r
		virtual void apply(const Vector<T>& r, Vector<T>& z) const = 0;

		// Optional: setup from matrix (for ILU, etc.)
    virtual void setup(const SparseMatrixCSR<T>& /*A*/) {}
};

// Identity preconditioner (no preconditioning)
template<typename T>
class IdentityPreconditioner : public Preconditioner<T> {
public:
	void apply(const Vector<T>& r, Vector<T>& z) const override { z = r; }
};

//=============================================================================
// Progress Callback
//=============================================================================

template<typename T>
using SolverCallback = std::function<void(int iteration, T residual)>;

//=============================================================================
// Iterative Solver Base Class
//=============================================================================

template<typename T>
class IterativeSolverBase {
protected:
    SolverConfig<T> config_;
    SolverCallback<T> callback_;
    std::shared_ptr<Preconditioner<T>> preconditioner_;
    
public:
    IterativeSolverBase() : preconditioner_(std::make_shared<IdentityPreconditioner<T>>()) {}
    virtual ~IterativeSolverBase() = default;
    
    // Configuration
    void setConfig(const SolverConfig<T>& cfg) { config_ = cfg; }
    SolverConfig<T>& config() { return config_; }
    const SolverConfig<T>& config() const { return config_; }
    
    // Preconditioner
    void setPreconditioner(std::shared_ptr<Preconditioner<T>> prec) {
        preconditioner_ = prec ? prec : std::make_shared<IdentityPreconditioner<T>>();
    }
    
    // Progress callback
    void setCallback(SolverCallback<T> cb) { callback_ = cb; }
    
    // Main solve interface
		[[nodiscard]] KrylovResult<T> solve(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x) {
			AlgorithmTimer timer;
			KrylovResult<T> result;
			result.algorithm_name = algorithmName();
			if (!validateInput(A, b, x) || A.rows() == 0 || !std::isfinite(config_.tolerance) || config_.tolerance <= 0 ||
					!std::isfinite(config_.absTolerance) || config_.absTolerance < 0 || config_.maxIterations < 0) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "Krylov solver requires nonempty matching dimensions, finite positive relative tolerance, "
															 "nonnegative absolute tolerance and iteration budget";
			} else {
				result = solveImpl(A, b, x);
				result.algorithm_name = algorithmName();
				if (result.status != AlgorithmStatus::InvalidInput) {
					result.residualNorm = norm2(residual(A, x, b));
					const T rhsNorm = norm2(b);
					result.relativeResidual = rhsNorm > T(0) ? result.residualNorm / rhsNorm : result.residualNorm;
					result.achieved_tolerance = static_cast<Real>(result.residualNorm);
					if (!std::isfinite(result.residualNorm) || !std::isfinite(rhsNorm)) {
						result.status = AlgorithmStatus::NumericalInstability;
						result.error_message = "Krylov solver returned a non-finite state or true residual";
					} else if (result.status == AlgorithmStatus::Success && !checkConvergence(result.residualNorm, rhsNorm)) {
						result.status = AlgorithmStatus::Stalled;
						result.error_message = "Recursive convergence estimate did not satisfy the true residual criterion";
					}
				}
			}
			result.converged = result.status == AlgorithmStatus::Success;
			if (result.converged)
				result.error_message.clear();
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		// Convenience: solve with zero initial guess
		KrylovResult<T> solve(const SparseMatrixCSR<T>& A, const Vector<T>& b) {
			Vector<T> x(b.size(), T(0));
			return solve(A, b, x);
		}

protected:
	virtual KrylovResult<T> solveImpl(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x) = 0;
	virtual const char* algorithmName() const noexcept = 0;

	// Check convergence: ||r|| < tol * ||b|| + absTol
	bool checkConvergence(T residualNorm, T rhsNorm) const { 	return residualNorm == T(0) || residualNorm < config_.tolerance * rhsNorm + config_.absTolerance; }

	// Report progress
	void reportProgress(int iter, T residual) const {
		if (callback_)
			callback_(iter, residual);
		if (config_.verbose) {
			std::cout << "Iter " << iter << ": residual = " << residual << "\n";
		}
    }
    
    // Validate input dimensions
		bool validateInput(const SparseMatrixCSR<T>& A, const Vector<T>& b, const Vector<T>& x) const {
			return A.rows() == A.cols() && A.rows() == static_cast<int>(b.size()) && A.rows() == static_cast<int>(x.size());
		}
};

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVER_BASE_H
