///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgEqSolvers_iterative.h                                         ///
///  Description: Dense stationary linear solvers (Jacobi, Gauss-Seidel, SOR)          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LINEAR_ALG_EQ_SOLVERS_ITERATIVE_H
#define MML_LINEAR_ALG_EQ_SOLVERS_ITERATIVE_H

#include <mml/MMLBase.h>

#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                    IterativeSolverConfig                            ///
	///////////////////////////////////////////////////////////////////////////
	/**
	 * @brief Configuration for iterative linear equation solvers
	 * 
	 * Residuals are normalized by ||b||, or absolute when b is zero.
	 * Zero budget checks only the initial guess; negative budgets are invalid.
	 */
	struct IterativeSolverConfig
	{
		Real tolerance = PrecisionValues<Real>::IterativeSolverTolerance;       ///< Finite positive normalized-residual tolerance
		int max_iterations = 1000;    ///< Nonnegative sweep budget; zero checks only the initial guess
		Real omega = 1.0;             ///< Relaxation parameter (for SOR, 0 < omega < 2)
		bool verbose = false;         ///< Reserved; dense solvers do not emit progress output

		/// Default constructor
		IterativeSolverConfig() = default;

		/// Constructor with tolerance and max iterations
		IterativeSolverConfig(Real tol, int max_iter = 1000)
			: tolerance(tol)
			, max_iterations(max_iter) {}

		/// Factory: High precision configuration
		static IterativeSolverConfig HighPrecision() {
			IterativeSolverConfig cfg;
			cfg.tolerance = PrecisionValues<Real>::IterativeSolverTolerance / Real(10000);
			cfg.max_iterations = 10000;
			return cfg;
		}

		/// Factory: Fast configuration (lower precision, fewer iterations)
		static IterativeSolverConfig Fast() {
			IterativeSolverConfig cfg;
			cfg.tolerance = Precision::LooseTolerance;
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
	struct [[nodiscard]] IterativeSolverResult
	{
		Vector<Real> solution;      ///< Solution vector x
		int          iterations;    ///< Sweeps attempted, including a sweep that fails numerically
		Real         residual;      ///< ||Ax-b||/||b|| for returned x; absolute ||Ax|| if b=0
		bool         converged;     ///< True if converged within tolerance
		
		// Enhanced diagnostic fields (API Standardization Phase 5)
		std::string algorithm_name;  ///< Name of the algorithm used
		AlgorithmStatus status = AlgorithmStatus::InvalidInput;  ///< Algorithm termination status
		std::string error_message = "No dense iterative solve performed";   ///< Error message if failed
		double elapsed_time_ms = 0;  ///< Execution time in milliseconds
		Real initial_residual = NoEstimate;   ///< Initial residual, in the same units as residual
		Real omega_used = 0;         ///< Relaxation parameter used (for SOR)
		
		IterativeSolverResult() : iterations(0), residual(NoEstimate), converged(false) {}
		IterativeSolverResult(const Vector<Real>& sol, int iter, Real res, bool conv)
			: solution(sol), iterations(iter), residual(res), converged(conv),
			  status(conv ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded),
			  error_message(conv ? "" : "Iteration budget exhausted") {}
	};

	namespace DenseIterativeDetail
	{
		enum class Method { Jacobi, GaussSeidel, SOR };

		inline bool IsFinite(const Vector<Real>& vector)
		{
			for (int i = 0; i < vector.size(); ++i)
				if (!std::isfinite(vector[i]))
					return false;
			return true;
		}

		inline IterativeSolverResult Solve(const Matrix<Real>& a, const Vector<Real>& b,
			const Vector<Real>& guess, Real tolerance, int budget, Method method, const char* name, Real omega = 0)
		{
			AlgorithmTimer timer;
			const std::string prefix = std::string(name) + "::Solve - ";
			const int n = a.rows();
			if (n == 0 || a.cols() != n)
				throw MatrixDimensionError(prefix + "matrix must be nonempty and square", n, a.cols(), -1, -1);
			if (b.size() != n || (guess.size() != 0 && guess.size() != n))
				throw MatrixDimensionError(prefix + "RHS and nonempty initial guess must match matrix size",
					n, n, b.size(), guess.size());
			if (!std::isfinite(tolerance) || tolerance <= 0 || budget < 0)
				throw ArgumentError(prefix + "tolerance must be finite and positive; iteration budget must be nonnegative");
			if (!IsFinite(b) || !IsFinite(guess))
				throw ArgumentError(prefix + "RHS and initial guess entries must be finite");

			Real maxDiagonal = 0;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j)
					if (!std::isfinite(a(i, j)))
						throw ArgumentError(prefix + "matrix entries must be finite");
				maxDiagonal = std::max(maxDiagonal, std::abs(a(i, i)));
			}
			for (int i = 0; i < n; ++i)
				if (std::abs(a(i, i)) <= maxDiagonal * Precision::Eps)
					throw SingularMatrixError(prefix + "near-zero diagonal element at row " + std::to_string(i));

			IterativeSolverResult result;
			result.solution = guess.size() == 0 ? Vector<Real>(n, Real(0)) : guess;
			result.algorithm_name = name;
			result.omega_used = omega;
			const auto finish = [&]() {
				result.converged = result.status == AlgorithmStatus::Success;
				result.elapsed_time_ms = timer.elapsed_ms();
				return std::move(result);
			};

			const Real rhsNorm = b.NormL2();
			const Real scale = rhsNorm == 0 ? Real(1) : rhsNorm;
			const Real initial = std::isfinite(rhsNorm) ? (a * result.solution - b).NormL2() / scale : NoEstimate;
			if (!std::isfinite(initial)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = prefix + "RHS norm or initial true residual is not representable";
				return finish();
			}
			result.initial_residual = initial;
			result.residual = initial;
			result.status = AlgorithmStatus::MaxIterationsExceeded;

			if (initial >= tolerance && budget > 0) {
				Vector<Real> current = result.solution;
				Vector<Real> next(n);
				for (int iteration = 0; iteration < budget && result.residual >= tolerance; ++iteration) {
					for (int i = 0; i < n; ++i) {
						Real sigma = 0;
						for (int j = 0; j < n; ++j)
							if (j != i)
								sigma += a(i, j) * (method != Method::Jacobi && j < i ? next[j] : current[j]);
						const Real update = (b[i] - sigma) / a(i, i);
						next[i] = method == Method::SOR ? (1 - omega) * current[i] + omega * update : update;
					}
					++result.iterations;
					const Real residual = IsFinite(next) ? (a * next - b).NormL2() / scale : NoEstimate;
					if (!std::isfinite(residual)) {
						result.status = AlgorithmStatus::NumericalInstability;
						result.error_message = prefix + "nonfinite iterate or true residual at sweep " + std::to_string(result.iterations);
						return finish();
					}
					if (residual < result.residual) {
						result.solution = next;
						result.residual = residual;
					}
					std::swap(current, next);
				}
			}
			if (result.residual < tolerance) {
				result.status = AlgorithmStatus::Success;
				result.error_message.clear();
			} else {
				result.error_message = prefix + "failed to converge after " + std::to_string(result.iterations) +
					" sweeps; the normalized residual tolerance was not met";
			}
			return finish();
		}
	}

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
	 * - Symmetric positive definiteness alone does not guarantee convergence
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
		                                    Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
		                                    int maxIter = 1000)
		{
			return DenseIterativeDetail::Solve(A, b, x0, tol, maxIter, DenseIterativeDetail::Method::Jacobi, "Jacobi");
		}
		
		/**
		 * @brief Simple solve returning only solution vector
		 * @throws std::runtime_error if solver doesn't converge
		 */
		static Vector<Real> SolveSimple(const Matrix<Real>& A, 
		                                 const Vector<Real>& b,
		                                 Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
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
			return Solve(A, b, x0, config.tolerance, config.max_iterations);
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
		                                    Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
		                                    int maxIter = 1000)
		{
			return DenseIterativeDetail::Solve(A, b, x0, tol, maxIter, DenseIterativeDetail::Method::GaussSeidel, "GaussSeidel");
		}
		
		/**
		 * @brief Simple solve returning only solution vector
		 * @throws std::runtime_error if solver doesn't converge
		 */
		static Vector<Real> SolveSimple(const Matrix<Real>& A, 
		                                 const Vector<Real>& b,
		                                 Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
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
			return Solve(A, b, x0, config.tolerance, config.max_iterations);
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
		                                    Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
		                                    int maxIter = 1000)
		{
			// Validate omega
			if (!std::isfinite(omega) || omega <= 0.0 || omega >= 2.0)
				throw NumericalMethodError("SORSolver::Solve - omega must be in (0, 2), got " + std::to_string(omega));
			return DenseIterativeDetail::Solve(A, b, x0, tol, maxIter, DenseIterativeDetail::Method::SOR, "SOR", omega);
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
		                                           Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
		                                           int maxIter = 1000)
		{
			int n = A.rows();
			Real omega = n <= 1 ? Real(1) : MML::Real(2.0) / (MML::Real(1.0) + std::sin(Constants::PI / n));
			return Solve(A, b, omega, x0, tol, maxIter);
		}
		
		/**
		 * @brief Simple solve returning only solution vector
		 * @throws std::runtime_error if solver doesn't converge
		 */
		static Vector<Real> SolveSimple(const Matrix<Real>& A, 
		                                 const Vector<Real>& b,
		                                 Real omega,
		                                 Real tol = PrecisionValues<Real>::IterativeSolverTolerance,
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
			if (numSamples <= 0 || testIter <= 0)
				throw ArgumentError("SORSolver::EstimateOptimalOmega - sample count and trial sweep budget must be positive");
			Real bestOmega = 1.0;
			Real bestResidual = std::numeric_limits<Real>::infinity();
			
			// Try omega values from 1.0 to 1.95
			for (int i = 0; ; ++i)
			{
				Real omega = MML::Real(1.0) + MML::Real(0.95) * i / numSamples;
				auto result = Solve(A, b, omega, Vector<Real>(), PrecisionValues<Real>::EigenSolverZeroThreshold, testIter);
				if (result.status != AlgorithmStatus::NumericalInstability && result.residual < bestResidual)
				{
					bestResidual = result.residual;
					bestOmega = omega;
				}
				if (i == numSamples)
					break;
			}
			if (!std::isfinite(bestResidual))
				throw NumericalMethodError("SORSolver::EstimateOptimalOmega - all trial solves failed numerically");
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
			return Solve(A, b, config.omega, x0, config.tolerance, config.max_iterations);
		}
	};

} // namespace MML

#endif // MML_LINEAR_ALG_EQ_SOLVERS_ITERATIVE_H
