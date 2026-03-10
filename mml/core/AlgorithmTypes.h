///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        AlgorithmTypes.h                                                    ///
///  Description: Base types for algorithm configuration and results                  ///
///               Provides consistent patterns across all MML algorithms              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGORITHM_TYPES_H
#define MML_ALGORITHM_TYPES_H

#include "MMLBase.h"

#include <string>
#include <chrono>

namespace MML {

	/******************************************************************************/
	/*****                    Algorithm Status Codes                          *****/
	/******************************************************************************/

	/// Status codes for iterative algorithm outcomes.
	/// 
	/// Provides structured error information beyond simple converged/not-converged.
	/// Use these codes for programmatic error handling.
	enum class AlgorithmStatus {
		/// Algorithm converged successfully
		Success = 0,
		
		/// Algorithm did not converge within max_iterations
		MaxIterationsExceeded = 1,
		
		/// Numerical instability detected (NaN, Inf, or ill-conditioning)
		NumericalInstability = 2,
		
		/// Matrix is singular or nearly singular
		SingularMatrix = 3,
		
		/// Input validation failed (dimensions, constraints)
		InvalidInput = 4,
		
		/// Algorithm stalled (no progress between iterations)
		Stalled = 5,
		
		/// User-requested tolerance is unachievable
		ToleranceUnachievable = 6,
		
		/// Algorithm-specific failure (see error_message for details)
		AlgorithmSpecificFailure = 100
	};

	/// Convert AlgorithmStatus to human-readable string.
	inline std::string ToString(AlgorithmStatus status) {
		switch (status) {
			case AlgorithmStatus::Success:                return "Success";
			case AlgorithmStatus::MaxIterationsExceeded:  return "MaxIterationsExceeded";
			case AlgorithmStatus::NumericalInstability:   return "NumericalInstability";
			case AlgorithmStatus::SingularMatrix:         return "SingularMatrix";
			case AlgorithmStatus::InvalidInput:           return "InvalidInput";
			case AlgorithmStatus::Stalled:                return "Stalled";
			case AlgorithmStatus::ToleranceUnachievable:  return "ToleranceUnachievable";
			case AlgorithmStatus::AlgorithmSpecificFailure: return "AlgorithmSpecificFailure";
			default:                                      return "Unknown";
		}
	}

	/******************************************************************************/
	/*****                    Base Configuration Type                         *****/
	/******************************************************************************/

	/// Base configuration parameters shared by all iterative algorithms.
	/// 
	/// Provides the minimal set of parameters that every iterative algorithm needs.
	/// Derived Config structs should include these fields (via composition or inheritance).
	/// 
	/// @note This is a documentation/pattern struct. Actual algorithm configs
	///       should define their own fields for type safety and discoverability.
	/// 
	/// @example
	/// // In EigenSystemSolvers.h:
	/// struct EigenSolverConfig {
	///     Real tolerance = 1e-10;                         // Standard default
	///     int max_iterations = 100;                       // Standard default
	///     bool verbose = false;                           // From base pattern
	///     bool sort_eigenvalues = true;                   // Algorithm-specific
	///     bool compute_eigenvectors = true;               // Algorithm-specific
	/// };
	struct ConfigBase {
		/// Convergence tolerance (algorithm-dependent interpretation)
		/// - Root finding: |f(x)| < tolerance
		/// - Eigen solvers: max residual < tolerance
		/// - Optimization: |gradient| < tolerance or |x_{n+1} - x_n| < tolerance
		Real tolerance = 1e-10;
		
		/// Maximum number of iterations before giving up
		/// Set to 0 to use algorithm-specific default
		int max_iterations = 100;
		
		/// Enable verbose diagnostic output to std::cout
		/// Useful for debugging convergence issues
		bool verbose = false;
	};

	/******************************************************************************/
	/*****                    Base Result Type                                *****/
	/******************************************************************************/

	/// Base result fields shared by all iterative algorithm results.
	/// 
	/// Provides the minimal diagnostic information that every algorithm should return.
	/// Derived Result structs should include these fields and add algorithm-specific data.
	/// 
	/// @note This is a documentation/pattern struct. Actual algorithm results
	///       should define their own fields for type safety and discoverability.
	/// 
	/// @example
	/// // In EigenSystemSolvers.h:
	/// struct EigenSolverResult {
	///     // === Algorithm output ===
	///     Vector<Real> eigenvalues;
	///     Matrix<Real> eigenvectors;
	///     
	///     // === Convergence info (from IterativeResultBase pattern) ===
	///     bool converged = false;
	///     int iterations_used = 0;
	///     Real achieved_tolerance = 0.0;
	///     AlgorithmStatus status = AlgorithmStatus::Success;
	///     std::string error_message;
	///     
	///     // === Diagnostics ===
	///     std::string algorithm_name;
	///     double elapsed_time_ms = 0.0;
	///     Real max_residual = 0.0;  // Algorithm-specific
	/// };
	struct IterativeResultBase {
		/// True if algorithm converged within tolerance and iteration limits
		bool converged = false;
		
		/// Number of iterations actually performed
		int iterations_used = 0;
		
		/// Actual achieved tolerance (may be better than requested)
		/// Interpretation depends on algorithm (residual, error estimate, etc.)
		Real achieved_tolerance = 0.0;
		
		/// Structured status code for programmatic handling
		AlgorithmStatus status = AlgorithmStatus::Success;
		
		/// Human-readable error message (empty string if successful)
		/// Should be specific and actionable, e.g.:
		/// "Failed to converge after 100 iterations (residual: 1.5e-4)"
		std::string error_message;
		
		/// Name of the algorithm that produced this result
		/// Useful for logging and debugging, e.g. "Jacobi", "Brent", "RKF45"
		std::string algorithm_name;
		
		/// Wall-clock execution time in milliseconds
		/// Measured from algorithm start to completion
		double elapsed_time_ms = 0.0;
		
		/// Number of function evaluations performed
		/// Useful for comparing algorithm efficiency
		int function_evaluations = 0;
	};

	/******************************************************************************/
	/*****                    Algorithm Timer Utility                         *****/
	/******************************************************************************/

	/// RAII-style timer for measuring algorithm execution time.
	/// 
	/// Creates a timer on construction, call elapsed_ms() to get elapsed time.
	/// Designed for use in algorithm implementations to populate elapsed_time_ms.
	/// 
	/// @example
	/// EigenSolverResult Jacobi_Solve(const Matrix<Real>& A, const EigenSolverConfig& config) {
	///     AlgorithmTimer timer;
	///     EigenSolverResult result;
	///     // ... algorithm implementation ...
	///     result.elapsed_time_ms = timer.elapsed_ms();
	///     return result;
	/// }
	class AlgorithmTimer {
	public:
		AlgorithmTimer() : start_(std::chrono::high_resolution_clock::now()) {}
		
		/// Get elapsed time since construction in milliseconds
		[[nodiscard]] double elapsed_ms() const {
			auto now = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - start_);
			return static_cast<double>(duration.count()) / 1000.0;
		}
		
		/// Reset the timer to current time
		void reset() {
			start_ = std::chrono::high_resolution_clock::now();
		}
		
	private:
		std::chrono::time_point<std::chrono::high_resolution_clock> start_;
	};

	/******************************************************************************/
	/*****                    Result Helper Functions                         *****/
	/******************************************************************************/

	/// Create a success result with converged=true and Success status.
	/// 
	/// @tparam ResultType The algorithm-specific result type
	/// @param iterations Number of iterations used
	/// @param tolerance Achieved tolerance
	/// @param algorithm_name Name of the algorithm
	/// @return ResultType with success fields populated
	/// 
	/// @example
	/// auto result = MakeSuccessResult<EigenSolverResult>(sweeps, max_off_diag, "Jacobi");
	/// result.eigenvalues = eigenvalues;  // Add algorithm-specific data
	template<typename ResultType>
	ResultType MakeSuccessResult(int iterations, Real tolerance, const std::string& algorithm_name) {
		ResultType result;
		result.converged = true;
		result.iterations_used = iterations;
		result.achieved_tolerance = tolerance;
		result.status = AlgorithmStatus::Success;
		result.algorithm_name = algorithm_name;
		result.error_message = "";
		return result;
	}

	/// Create a failure result with converged=false and specified status.
	/// 
	/// @tparam ResultType The algorithm-specific result type
	/// @param status The failure status code
	/// @param error_message Descriptive error message
	/// @param iterations Number of iterations performed before failure
	/// @param algorithm_name Name of the algorithm
	/// @return ResultType with failure fields populated
	/// 
	/// @example
	/// if (iter >= config.max_iterations) {
	///     return MakeFailureResult<EigenSolverResult>(
	///         AlgorithmStatus::MaxIterationsExceeded,
	///         "Did not converge after " + std::to_string(iter) + " iterations",
	///         iter, "Jacobi");
	/// }
	template<typename ResultType>
	ResultType MakeFailureResult(AlgorithmStatus status, const std::string& error_message,
	                              int iterations, const std::string& algorithm_name) {
		ResultType result;
		result.converged = false;
		result.iterations_used = iterations;
		result.status = status;
		result.error_message = error_message;
		result.algorithm_name = algorithm_name;
		return result;
	}

} // namespace MML

#endif // MML_ALGORITHM_TYPES_H
