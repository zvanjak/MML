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
    T tolerance = static_cast<T>(1e-10);   // Relative residual tolerance
    T absTolerance = static_cast<T>(1e-14); // Absolute residual tolerance
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

enum class SolverStatus {
    Success,           // Converged within tolerance
    MaxIterations,     // Did not converge in max iterations
    Stagnation,        // Progress stalled
    Breakdown,         // Method breakdown (e.g., division by zero)
    InvalidInput       // Invalid matrix or vector dimensions
};

template<typename T>
struct SolverResult {
    SolverStatus status = SolverStatus::Success;
    int iterations = 0;
    T residualNorm = T(0);
    T relativeResidual = T(0);
    std::string message;
    
    bool converged() const { return status == SolverStatus::Success; }
    
    operator bool() const { return converged(); }
};

//=============================================================================
// Preconditioner Interface
//=============================================================================

// Abstract preconditioner: solves M*z = r approximately
template<typename T>
class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    
    // Apply preconditioner: z = M^{-1} * r
    virtual void apply(const std::vector<T>& r, std::vector<T>& z) const = 0;
    
    // Optional: setup from matrix (for ILU, etc.)
    virtual void setup(const SparseMatrixCSR<T>& /*A*/) {}
};

// Identity preconditioner (no preconditioning)
template<typename T>
class IdentityPreconditioner : public Preconditioner<T> {
public:
    void apply(const std::vector<T>& r, std::vector<T>& z) const override {
        z = r;
    }
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
    virtual SolverResult<T> solve(
        const SparseMatrixCSR<T>& A,
        const std::vector<T>& b,
        std::vector<T>& x) = 0;
    
    // Convenience: solve with zero initial guess
    SolverResult<T> solve(const SparseMatrixCSR<T>& A, const std::vector<T>& b) {
        std::vector<T> x(b.size(), T(0));
        return solve(A, b, x);
    }
    
protected:
    // Check convergence: ||r|| < tol * ||b|| + absTol
    bool checkConvergence(T residualNorm, T rhsNorm) const {
        return residualNorm < config_.tolerance * rhsNorm + config_.absTolerance;
    }
    
    // Report progress
    void reportProgress(int iter, T residual) const {
        if (callback_) callback_(iter, residual);
        if (config_.verbose) {
            std::cout << "Iter " << iter << ": residual = " << residual << "\n";
        }
    }
    
    // Validate input dimensions
    bool validateInput(const SparseMatrixCSR<T>& A, const std::vector<T>& b, 
                       const std::vector<T>& x) const {
        return A.rows() == A.cols() && 
               A.rows() == static_cast<int>(b.size()) &&
               A.rows() == static_cast<int>(x.size());
    }
};

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVER_BASE_H
