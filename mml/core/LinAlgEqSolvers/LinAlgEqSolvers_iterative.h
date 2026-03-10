///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgEqSolvers_iterative.h                                         ///
///  Description: Iterative linear system solvers (Jacobi, Gauss-Seidel, SOR, CG)     ///
///               For large sparse systems where direct methods are impractical       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LINEAR_ALG_EQ_SOLVERS_ITERATIVE_H
#define MML_LINEAR_ALG_EQ_SOLVERS_ITERATIVE_H

#include "MMLBase.h"

#include "core/AlgorithmTypes.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                    IterativeSolverConfig                            ///
	///////////////////////////////////////////////////////////////////////////
	/**
	 * @brief Configuration for iterative linear equation solvers
	 * 
	 * Standardized config object following API Standardization Phase 5 pattern.
	 */
	struct IterativeSolverConfig
	{
		Real tolerance = 1e-10;       ///< Convergence tolerance
		int max_iterations = 1000;    ///< Maximum iterations
		Real omega = 1.0;             ///< Relaxation parameter (for SOR, 0 < omega < 2)
		bool verbose = false;         ///< Enable verbose output

		/// Default constructor
		IterativeSolverConfig() = default;

		/// Constructor with tolerance and max iterations
		IterativeSolverConfig(Real tol, int max_iter = 1000)
			: tolerance(tol)
			, max_iterations(max_iter) {}

		/// Factory: High precision configuration
		static IterativeSolverConfig HighPrecision() {
			IterativeSolverConfig cfg;
			cfg.tolerance = 1e-14;
			cfg.max_iterations = 10000;
			return cfg;
		}

		/// Factory: Fast configuration (lower precision, fewer iterations)
		static IterativeSolverConfig Fast() {
			IterativeSolverConfig cfg;
			cfg.tolerance = 1e-6;
			cfg.max_iterations = 200;
			return cfg;
		}

		/// Factory: SOR configuration with optimal omega
		static IterativeSolverConfig SOR(Real omega_val) {
			IterativeSolverConfig cfg;
			cfg.omega = omega_val;
			return cfg;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                    IterativeSolverResult                            ///
	///////////////////////////////////////////////////////////////////////////
	/**
	 * @brief Result structure for iterative linear equation solvers
	 * 
	 * Contains solution vector along with convergence information.
	 * Enhanced with diagnostic fields for API Standardization Phase 5.
	 */
	struct IterativeSolverResult
	{
		Vector<Real> solution;      ///< Solution vector x
		int          iterations;    ///< Number of iterations performed
		Real         residual;      ///< Final residual norm ||Ax - b||
		bool         converged;     ///< True if converged within tolerance
		
		// Enhanced diagnostic fields (API Standardization Phase 5)
		std::string algorithm_name;  ///< Name of the algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;  ///< Algorithm termination status
		std::string error_message;   ///< Error message if failed
		double elapsed_time_ms = 0;  ///< Execution time in milliseconds
		Real initial_residual = 0;   ///< Initial residual norm (before iterations)
		Real omega_used = 0;         ///< Relaxation parameter used (for SOR)
		
		IterativeSolverResult() : iterations(0), residual(0.0), converged(false) {}
		IterativeSolverResult(const Vector<Real>& sol, int iter, Real res, bool conv)
			: solution(sol), iterations(iter), residual(res), converged(conv) {}
	};

	/**
	 * @brief Jacobi iterative method for solving linear systems Ax = b
	 * 
	 * The Jacobi method is a classical iterative technique that updates each component
	 * of the solution using only values from the previous iteration:
	 * 
	 *   x_i^(k+1) = (1/a_ii) * (b_i - Σ_{j≠i} a_ij * x_j^(k))
	 * 
	 * @par Convergence Properties:
	 * - Guaranteed to converge for strictly diagonally dominant matrices
	 * - Converges for symmetric positive definite matrices
	 * - Each iteration can be parallelized (all updates independent)
	 * - Generally slower than Gauss-Seidel but more parallelizable
	 * 
	 * @par Convergence Rate:
	 * The spectral radius ρ(D⁻¹(L+U)) determines convergence rate.
	 * For tridiagonal Poisson: ρ ≈ 1 - π²/(2n²), very slow for large n.
	 * 
	 * @see GaussSeidelSolver, SORSolver
	 */
	class JacobiSolver
	{
	public:
		/**
		 * @brief Solve Ax = b using Jacobi iteration
		 * 
		 * @param A Coefficient matrix (must be square, non-zero diagonal)
		 * @param b Right-hand side vector
		 * @param x0 Initial guess (if empty, uses zero vector)
		 * @param tol Convergence tolerance for relative residual ||Ax-b||/||b||
		 * @param maxIter Maximum number of iterations
		 * @return IterativeSolverResult containing solution and convergence info
		 * 
		 * @throws MatrixDimensionError if dimensions don't match
		 * @throws SingularMatrixError if any diagonal element is zero
		 */
		static IterativeSolverResult Solve(const Matrix<Real>& A, 
		                                    const Vector<Real>& b,
		                                    const Vector<Real>& x0 = Vector<Real>(),
		                                    Real tol = 1e-10,
		                                    int maxIter = 1000)
		{
			int n = A.rows();
			
			// Dimension checks
			if (A.cols() != n)
				throw MatrixDimensionError("JacobiSolver::Solve - matrix must be square", n, A.cols(), -1, -1);
			if (b.size() != n)
				throw MatrixDimensionError("JacobiSolver::Solve - vector b dimension mismatch", n, n, b.size(), -1);
			
			// Check for zero diagonal elements
			for (int i = 0; i < n; i++)
			{
				if (A(i, i) == 0.0)
					throw SingularMatrixError("JacobiSolver::Solve - zero diagonal element at row " + std::to_string(i));
			}
			
			// Initialize solution vector
			Vector<Real> x(n);
			if (x0.size() == n)
				x = x0;
			else
				for (int i = 0; i < n; i++) x[i] = 0.0;
			
			Vector<Real> x_new(n);
			Real b_norm = b.NormL2();
			if (b_norm == 0.0) b_norm = 1.0;  // Avoid division by zero
			
			int iter = 0;
			Real residual = 0.0;
			bool converged = false;
			
			for (iter = 1; iter <= maxIter; iter++)
			{
				// Jacobi iteration: x_new[i] = (b[i] - sum_{j!=i} A[i][j]*x[j]) / A[i][i]
				for (int i = 0; i < n; i++)
				{
					Real sigma = 0.0;
					for (int j = 0; j < n; j++)
					{
						if (j != i)
							sigma += A(i, j) * x[j];
					}
					x_new[i] = (b[i] - sigma) / A(i, i);
				}
				
				// Compute residual ||Ax_new - b||
				Vector<Real> r = A * x_new - b;
				residual = r.NormL2() / b_norm;
				
				// Check convergence
				if (residual < tol)
				{
					converged = true;
					x = x_new;
					break;
				}
				
				// Update for next iteration
				x = x_new;
			}
			
			return IterativeSolverResult(x, iter, residual, converged);
		}
		
		/**
		 * @brief Simple solve returning only solution vector
		 * @throws std::runtime_error if solver doesn't converge
		 */
		static Vector<Real> SolveSimple(const Matrix<Real>& A, 
		                                 const Vector<Real>& b,
		                                 Real tol = 1e-10,
		                                 int maxIter = 1000)
		{
			auto result = Solve(A, b, Vector<Real>(), tol, maxIter);
			if (!result.converged)
				throw ConvergenceError("JacobiSolver::SolveSimple - failed to converge after " + 
				                       std::to_string(result.iterations) + " iterations",
				                       result.iterations, result.residual);
			return result.solution;
		}

		/**
		 * @brief Solve Ax = b using Jacobi iteration with config object
		 * 
		 * @param A Coefficient matrix (must be square, non-zero diagonal)
		 * @param b Right-hand side vector
		 * @param config Solver configuration
		 * @param x0 Initial guess (if empty, uses zero vector)
		 * @return IterativeSolverResult with enhanced diagnostics
		 */
		static IterativeSolverResult Solve(const Matrix<Real>& A, 
		                                    const Vector<Real>& b,
		                                    const IterativeSolverConfig& config,
		                                    const Vector<Real>& x0 = Vector<Real>())
		{
			AlgorithmTimer timer;
			IterativeSolverResult result = Solve(A, b, x0, config.tolerance, config.max_iterations);
			
			result.algorithm_name = "Jacobi";
			result.elapsed_time_ms = timer.elapsed_ms();
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Jacobi method did not converge within " + 
				                       std::to_string(config.max_iterations) + " iterations";
			}
			return result;
		}
	};

	/**
	 * @brief Gauss-Seidel iterative method for solving linear systems Ax = b
	 * 
	 * The Gauss-Seidel method improves on Jacobi by using updated values
	 * immediately as they are computed:
	 * 
	 *   x_i^(k+1) = (1/a_ii) * (b_i - Σ_{j<i} a_ij * x_j^(k+1) - Σ_{j>i} a_ij * x_j^(k))
	 * 
	 * @par Convergence Properties:
	 * - Guaranteed to converge for strictly diagonally dominant matrices
	 * - Guaranteed to converge for symmetric positive definite matrices
	 * - Generally converges faster than Jacobi (about 2x for many problems)
	 * - Cannot be easily parallelized (sequential dependencies)
	 * 
	 * @par Convergence Rate:
	 * The spectral radius ρ((D+L)⁻¹U) determines convergence rate.
	 * For tridiagonal Poisson: ρ_GS ≈ ρ_J², roughly twice as fast as Jacobi.
	 * 
	 * @see JacobiSolver, SORSolver
	 */
	class GaussSeidelSolver
	{
	public:
		/**
		 * @brief Solve Ax = b using Gauss-Seidel iteration
		 * 
		 * @param A Coefficient matrix (must be square, non-zero diagonal)
		 * @param b Right-hand side vector
		 * @param x0 Initial guess (if empty, uses zero vector)
		 * @param tol Convergence tolerance for relative residual ||Ax-b||/||b||
		 * @param maxIter Maximum number of iterations
		 * @return IterativeSolverResult containing solution and convergence info
		 * 
		 * @throws MatrixDimensionError if dimensions don't match
		 * @throws SingularMatrixError if any diagonal element is zero
		 */
		static IterativeSolverResult Solve(const Matrix<Real>& A, 
		                                    const Vector<Real>& b,
		                                    const Vector<Real>& x0 = Vector<Real>(),
		                                    Real tol = 1e-10,
		                                    int maxIter = 1000)
		{
			int n = A.rows();
			
			// Dimension checks
			if (A.cols() != n)
				throw MatrixDimensionError("GaussSeidelSolver::Solve - matrix must be square", n, A.cols(), -1, -1);
			if (b.size() != n)
				throw MatrixDimensionError("GaussSeidelSolver::Solve - vector b dimension mismatch", n, n, b.size(), -1);
			
			// Check for zero diagonal elements
			for (int i = 0; i < n; i++)
			{
				if (A(i, i) == 0.0)
					throw SingularMatrixError("GaussSeidelSolver::Solve - zero diagonal element at row " + std::to_string(i));
			}
			
			// Initialize solution vector
			Vector<Real> x(n);
			if (x0.size() == n)
				x = x0;
			else
				for (int i = 0; i < n; i++) x[i] = 0.0;
			
			Real b_norm = b.NormL2();
			if (b_norm == 0.0) b_norm = 1.0;
			
			int iter = 0;
			Real residual = 0.0;
			bool converged = false;
			
			for (iter = 1; iter <= maxIter; iter++)
			{
				// Gauss-Seidel iteration: use updated values immediately
				for (int i = 0; i < n; i++)
				{
					Real sigma = 0.0;
					// Use new values for j < i (already computed this iteration)
					for (int j = 0; j < i; j++)
						sigma += A(i, j) * x[j];
					// Use old values for j > i
					for (int j = i + 1; j < n; j++)
						sigma += A(i, j) * x[j];
					
					x[i] = (b[i] - sigma) / A(i, i);
				}
				
				// Compute residual ||Ax - b||
				Vector<Real> r = A * x - b;
				residual = r.NormL2() / b_norm;
				
				// Check convergence
				if (residual < tol)
				{
					converged = true;
					break;
				}
			}
			
			return IterativeSolverResult(x, iter, residual, converged);
		}
		
		/**
		 * @brief Simple solve returning only solution vector
		 * @throws std::runtime_error if solver doesn't converge
		 */
		static Vector<Real> SolveSimple(const Matrix<Real>& A, 
		                                 const Vector<Real>& b,
		                                 Real tol = 1e-10,
		                                 int maxIter = 1000)
		{
			auto result = Solve(A, b, Vector<Real>(), tol, maxIter);
			if (!result.converged)
				throw ConvergenceError("GaussSeidelSolver::SolveSimple - failed to converge after " + 
				                       std::to_string(result.iterations) + " iterations",
				                       result.iterations, result.residual);
			return result.solution;
		}

		/**
		 * @brief Solve Ax = b using Gauss-Seidel with config object
		 * 
		 * @param A Coefficient matrix (must be square, non-zero diagonal)
		 * @param b Right-hand side vector
		 * @param config Solver configuration
		 * @param x0 Initial guess (if empty, uses zero vector)
		 * @return IterativeSolverResult with enhanced diagnostics
		 */
		static IterativeSolverResult Solve(const Matrix<Real>& A, 
		                                    const Vector<Real>& b,
		                                    const IterativeSolverConfig& config,
		                                    const Vector<Real>& x0 = Vector<Real>())
		{
			AlgorithmTimer timer;
			IterativeSolverResult result = Solve(A, b, x0, config.tolerance, config.max_iterations);
			
			result.algorithm_name = "GaussSeidel";
			result.elapsed_time_ms = timer.elapsed_ms();
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Gauss-Seidel method did not converge within " + 
				                       std::to_string(config.max_iterations) + " iterations";
			}
			return result;
		}
	};

	/**
	 * @brief Successive Over-Relaxation (SOR) method for solving linear systems Ax = b
	 * 
	 * SOR is a generalization of Gauss-Seidel that introduces a relaxation parameter ω:
	 * 
	 *   x_i^(k+1) = (1-ω)*x_i^(k) + (ω/a_ii) * (b_i - Σ_{j<i} a_ij * x_j^(k+1) - Σ_{j>i} a_ij * x_j^(k))
	 * 
	 * or equivalently:
	 *   x_i^(k+1) = x_i^(k) + ω * (x_i^GS - x_i^(k))
	 * 
	 * where x_i^GS is the Gauss-Seidel update.
	 * 
	 * @par Relaxation Parameter ω:
	 * - ω = 1: Reduces to Gauss-Seidel
	 * - 0 < ω < 1: Under-relaxation (can help convergence for some problems)
	 * - 1 < ω < 2: Over-relaxation (can significantly accelerate convergence)
	 * - ω ≥ 2 or ω ≤ 0: Method diverges
	 * 
	 * @par Optimal ω Selection:
	 * For tridiagonal systems from Poisson equation:
	 *   ω_opt = 2 / (1 + sin(π/n)) ≈ 2 - 2π/n
	 * 
	 * This can reduce iterations from O(n²) to O(n).
	 * 
	 * @par Convergence Properties:
	 * - Converges for 0 < ω < 2 if A is symmetric positive definite
	 * - With optimal ω, convergence can be dramatically faster than Gauss-Seidel
	 * 
	 * @see JacobiSolver, GaussSeidelSolver
	 */
	class SORSolver
	{
	public:
		/**
		 * @brief Solve Ax = b using Successive Over-Relaxation
		 * 
		 * @param A Coefficient matrix (must be square, non-zero diagonal)
		 * @param b Right-hand side vector
		 * @param omega Relaxation parameter (0 < ω < 2 for convergence)
		 * @param x0 Initial guess (if empty, uses zero vector)
		 * @param tol Convergence tolerance for relative residual ||Ax-b||/||b||
		 * @param maxIter Maximum number of iterations
		 * @return IterativeSolverResult containing solution and convergence info
		 * 
		 * @throws MatrixDimensionError if dimensions don't match
		 * @throws SingularMatrixError if any diagonal element is zero
		 * @throws std::invalid_argument if omega not in (0, 2)
		 */
		static IterativeSolverResult Solve(const Matrix<Real>& A, 
		                                    const Vector<Real>& b,
		                                    Real omega,
		                                    const Vector<Real>& x0 = Vector<Real>(),
		                                    Real tol = 1e-10,
		                                    int maxIter = 1000)
		{
			// Validate omega
			if (omega <= 0.0 || omega >= 2.0)
				throw NumericalMethodError("SORSolver::Solve - omega must be in (0, 2), got " + std::to_string(omega));
			
			int n = A.rows();
			
			// Dimension checks
			if (A.cols() != n)
				throw MatrixDimensionError("SORSolver::Solve - matrix must be square", n, A.cols(), -1, -1);
			if (b.size() != n)
				throw MatrixDimensionError("SORSolver::Solve - vector b dimension mismatch", n, n, b.size(), -1);
			
			// Check for zero diagonal elements
			for (int i = 0; i < n; i++)
			{
				if (A(i, i) == 0.0)
					throw SingularMatrixError("SORSolver::Solve - zero diagonal element at row " + std::to_string(i));
			}
			
			// Initialize solution vector
			Vector<Real> x(n);
			if (x0.size() == n)
				x = x0;
			else
				for (int i = 0; i < n; i++) x[i] = 0.0;
			
			Real b_norm = b.NormL2();
			if (b_norm == 0.0) b_norm = 1.0;
			
			int iter = 0;
			Real residual = 0.0;
			bool converged = false;
			
			for (iter = 1; iter <= maxIter; iter++)
			{
				// SOR iteration
				for (int i = 0; i < n; i++)
				{
					Real sigma = 0.0;
					// Use new values for j < i
					for (int j = 0; j < i; j++)
						sigma += A(i, j) * x[j];
					// Use old values for j > i
					for (int j = i + 1; j < n; j++)
						sigma += A(i, j) * x[j];
					
					// Gauss-Seidel update
					Real x_gs = (b[i] - sigma) / A(i, i);
					
					// SOR: weighted average of old value and Gauss-Seidel update
					x[i] = (1.0 - omega) * x[i] + omega * x_gs;
				}
				
				// Compute residual ||Ax - b||
				Vector<Real> r = A * x - b;
				residual = r.NormL2() / b_norm;
				
				// Check convergence
				if (residual < tol)
				{
					converged = true;
					break;
				}
			}
			
			return IterativeSolverResult(x, iter, residual, converged);
		}
		
		/**
		 * @brief Solve with automatic omega selection for tridiagonal systems
		 * 
		 * Uses the optimal ω formula: ω = 2 / (1 + sin(π/n))
		 * This is optimal for tridiagonal systems from discretized Poisson equation.
		 * For other systems, may not be optimal but often gives good results.
		 */
		static IterativeSolverResult SolveOptimal(const Matrix<Real>& A, 
		                                           const Vector<Real>& b,
		                                           const Vector<Real>& x0 = Vector<Real>(),
		                                           Real tol = 1e-10,
		                                           int maxIter = 1000)
		{
			int n = A.rows();
			Real omega = 2.0 / (1.0 + std::sin(Constants::PI / n));
			return Solve(A, b, omega, x0, tol, maxIter);
		}
		
		/**
		 * @brief Simple solve returning only solution vector
		 * @throws std::runtime_error if solver doesn't converge
		 */
		static Vector<Real> SolveSimple(const Matrix<Real>& A, 
		                                 const Vector<Real>& b,
		                                 Real omega,
		                                 Real tol = 1e-10,
		                                 int maxIter = 1000)
		{
			auto result = Solve(A, b, omega, Vector<Real>(), tol, maxIter);
			if (!result.converged)
				throw ConvergenceError("SORSolver::SolveSimple - failed to converge after " + 
				                       std::to_string(result.iterations) + " iterations",
				                       result.iterations, result.residual);
			return result.solution;
		}
		
		/**
		 * @brief Estimate optimal omega empirically by trying several values
		 * 
		 * Performs a small number of iterations with different omega values
		 * and returns the one that gives the smallest residual.
		 * 
		 * @param A Coefficient matrix
		 * @param b Right-hand side vector
		 * @param numSamples Number of omega values to try (default 10)
		 * @param testIter Number of iterations for each test (default 50)
		 * @return Estimated optimal omega value
		 */
		static Real EstimateOptimalOmega(const Matrix<Real>& A,
		                                  const Vector<Real>& b,
		                                  int numSamples = 10,
		                                  int testIter = 50)
		{
			Real bestOmega = 1.0;
			Real bestResidual = std::numeric_limits<Real>::max();
			
			// Try omega values from 1.0 to 1.95
			for (int i = 0; i <= numSamples; i++)
			{
				Real omega = 1.0 + 0.95 * i / numSamples;
				auto result = Solve(A, b, omega, Vector<Real>(), PrecisionValues<Real>::EigenSolverZeroThreshold, testIter);
				if (result.residual < bestResidual)
				{
					bestResidual = result.residual;
					bestOmega = omega;
				}
			}
			
			return bestOmega;
		}

		/**
		 * @brief Solve Ax = b using SOR with config object
		 * 
		 * @param A Coefficient matrix (must be square, non-zero diagonal)
		 * @param b Right-hand side vector
		 * @param config Solver configuration (omega from config.omega)
		 * @param x0 Initial guess (if empty, uses zero vector)
		 * @return IterativeSolverResult with enhanced diagnostics
		 */
		static IterativeSolverResult Solve(const Matrix<Real>& A, 
		                                    const Vector<Real>& b,
		                                    const IterativeSolverConfig& config,
		                                    const Vector<Real>& x0 = Vector<Real>())
		{
			AlgorithmTimer timer;
			IterativeSolverResult result = Solve(A, b, config.omega, x0, config.tolerance, config.max_iterations);
			
			result.algorithm_name = "SOR";
			result.elapsed_time_ms = timer.elapsed_ms();
			result.omega_used = config.omega;
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "SOR method (omega=" + std::to_string(config.omega) + 
				                       ") did not converge within " + 
				                       std::to_string(config.max_iterations) + " iterations";
			}
			return result;
		}
	};

} // namespace MML

#endif // MML_LINEAR_ALG_EQ_SOLVERS_ITERATIVE_H
