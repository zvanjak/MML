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
//   auto A = Laplacian2D<double>(100, 100);
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
#include <algorithm>

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
        // CG requires SPD. A symmetric matrix with positive diagonal that is weakly diagonally
        // dominant in every row (strictly in at least one) is positive definite by Gershgorin,
        // so that is the condition we actually check - over ALL rows, not the first ten.
        if (A.IsSymmetric() && isDiagonallyDominantWithPositiveDiagonal(A)) {
            return std::make_unique<ConjugateGradient<T>>();
        }
        // Otherwise, use BiCGSTAB as a good general-purpose choice
        return std::make_unique<BiCGSTAB<T>>();
    }
    
    // Solve with automatic solver selection
		static KrylovResult<T> solve(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x,
																 const SolverConfig<T>& config = SolverConfig<T>()) {
			auto solver = create(A);
			solver->setConfig(config);
			return solver->solve(A, b, x);
		}

		/// Sufficient condition for positive definiteness of a symmetric matrix (Gershgorin):
    /// every diagonal entry positive and >= the sum of |off-diagonal| in its row, with the
    /// inequality strict in at least one row. O(nnz).
    /// @note Until analysis v4 (CORE C4) this was `hasDiagonalDominance`, which only tested
    ///       diag > 0 on the first ten rows and routed indefinite matrices to CG.
    static bool isDiagonallyDominantWithPositiveDiagonal(const SparseMatrixCSR<T>& A) {
        const auto& rowPtr = A.rowPointers();
        const auto& colIdx = A.colIndices();
        const auto& vals   = A.values();
        bool strictSomewhere = false;
        for (int i = 0; i < A.rows(); ++i) {
            T diag = T(0);
            T offSum = T(0);
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                if (colIdx[k] == i) diag += vals[k];
                else                offSum += std::abs(vals[k]);
            }
            if (!(diag > T(0)) || diag < offSum)
                return false;
            if (diag > offSum)
                strictSomewhere = true;
        }
        return strictSomewhere;
    }
};

//=============================================================================
// Convenience Functions
//=============================================================================

// Auto-select solver and solve
template<typename T>
KrylovResult<T> autoSolve(const SparseMatrixCSR<T>& A, const Vector<T>& b, Vector<T>& x,
													const SolverConfig<T>& config = SolverConfig<T>()) {
	return AutoSolver<T>::solve(A, b, x, config);
}

// Solve Poisson equation: -∇²u = f on unit square with zero Dirichlet BC
// Returns solution u on (nx-2) × (ny-2) interior grid
template<typename T>
Vector<T> solvePoisson2D(int nx, int ny,
												 const Vector<T>& f, // RHS on interior points
												 const SolverConfig<T>& config = SolverConfig<T>()) {
	// Create Laplacian matrix
	auto A = Laplacian2D<T>(nx - 2, ny - 2); // Interior points only

	// Scale by h² (assuming unit square)
	T h = T(1) / T(nx - 1);
	T h2 = h * h;

	// RHS: b = h² * f
	Vector<T> b(f.size());
	for (size_t i = 0; i < f.size(); ++i) {
		b[i] = h2 * f[i];
	}

	// Solve with preconditioned CG
	Vector<T> u(b.size(), T(0));
	auto prec = std::make_shared<JacobiPreconditioner<T>>();
	prec->setup(A);

	auto result = solveCG(A, b, u, prec, config);

	if (!result.converged) {
		// Fallback to BiCGSTAB if CG fails
		std::fill(u.begin(), u.end(), T(0));
		result = solveBiCGSTAB(A, b, u, prec, config);
	}
	if (!result.converged)
		throw ConvergenceError(result.error_message, result.iterations_used, static_cast<Real>(result.residualNorm));

	return u;
}

} // namespace MML::SparseSolvers

#endif // MML_CORE_SPARSE_SOLVERS_H
