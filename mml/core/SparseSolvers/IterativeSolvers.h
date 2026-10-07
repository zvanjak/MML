///////////////////////////////////////////////////////////////////////////////////////////
// IterativeSolvers.h - Main header for iterative linear solvers
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// This header provides a unified interface to iterative linear solvers
// for sparse systems arising from PDE discretizations.
//
// Included Solvers:
// - ConjugateGradient (CG): For symmetric positive definite systems
// - BiCGSTAB: For general non-symmetric systems
// - GMRES: For general systems (more robust, more memory)
//
// Preconditioners:
// - IdentityPreconditioner: No preconditioning
// - JacobiPreconditioner: Diagonal scaling
// - SSORPreconditioner: Symmetric SOR
// - ILU0Preconditioner: Incomplete LU (no fill-in)
//
// Usage:
//   // Create Laplacian system
//   auto A = laplacian2D<double>(100, 100);
//   std::vector<double> b(10000, 1.0);  // RHS
//   std::vector<double> x(10000, 0.0);  // Initial guess
//
//   // Solve with CG
//   SolverConfig<double> config;
//   config.setTolerance(1e-10).setMaxIterations(1000);
//   auto result = solveCG(A, b, x, config);
//
//   // Or with preconditioner
//   auto prec = std::make_shared<JacobiPreconditioner<double>>();
//   prec->setup(A);
//   result = solveCG(A, b, x, prec, config);
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_CORE_SPARSE_SOLVERS_H
#define MML_CORE_SPARSE_SOLVERS_H

#include "../../base/SparseMatrix/SparseMatrix.h"

#include "IterativeSolverBase.h"

#include "ConjugateGradient.h"
#include "BiCGSTAB.h"
#include "GMRES.h"

#include "Preconditioners.h"

#include <memory>

namespace MML::SparseSolvers {

//=============================================================================
// Solver Selection Helper
//=============================================================================

// Automatic solver selection based on matrix properties
template<typename T>
class AutoSolver {
public:
    // Choose best solver for the given matrix
    static std::unique_ptr<IterativeSolverBase<T>> create(const SparseMatrixCSR<T>& A) {
        // If symmetric and likely SPD (diagonal > 0), use CG
        if (A.isSymmetric() && hasDiagonalDominance(A)) {
            return std::make_unique<ConjugateGradient<T>>();
        }
        // Otherwise, use BiCGSTAB as a good general-purpose choice
        return std::make_unique<BiCGSTAB<T>>();
    }
    
    // Solve with automatic solver selection
    static SolverResult<T> solve(
        const SparseMatrixCSR<T>& A,
        const std::vector<T>& b,
        std::vector<T>& x,
        const SolverConfig<T>& config = SolverConfig<T>()) 
    {
        auto solver = create(A);
        solver->setConfig(config);
        return solver->solve(A, b, x);
    }
    
private:
    static bool hasDiagonalDominance(const SparseMatrixCSR<T>& A) {
        // Quick heuristic: check first few rows
        const int checkRows = std::min(10, A.rows());
        for (int i = 0; i < checkRows; ++i) {
            T diag = A(i, i);
            if (diag <= T(0)) return false;
        }
        return true;
    }
};

//=============================================================================
// Convenience Functions
//=============================================================================

// Auto-select solver and solve
template<typename T>
SolverResult<T> autoSolve(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    const SolverConfig<T>& config = SolverConfig<T>())
{
    return AutoSolver<T>::solve(A, b, x, config);
}

// Solve Poisson equation: -∇²u = f on unit square with zero Dirichlet BC
// Returns solution u on (nx-2) × (ny-2) interior grid
template<typename T>
std::vector<T> solvePoisson2D(
    int nx, int ny,
    const std::vector<T>& f,  // RHS on interior points
    const SolverConfig<T>& config = SolverConfig<T>())
{
    // Create Laplacian matrix
    auto A = laplacian2D<T>(nx - 2, ny - 2);  // Interior points only
    
    // Scale by h² (assuming unit square)
    T h = T(1) / T(nx - 1);
    T h2 = h * h;
    
    // RHS: b = h² * f
    std::vector<T> b(f.size());
    for (size_t i = 0; i < f.size(); ++i) {
        b[i] = h2 * f[i];
    }
    
    // Solve with preconditioned CG
    std::vector<T> u(b.size(), T(0));
    auto prec = std::make_shared<JacobiPreconditioner<T>>();
    prec->setup(A);
    
    auto result = solveCG(A, b, u, prec, config);
    
    if (!result.converged()) {
        // Fallback to BiCGSTAB if CG fails
        std::fill(u.begin(), u.end(), T(0));
        result = solveBiCGSTAB(A, b, u, prec, config);
    }
    
    return u;
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_H
