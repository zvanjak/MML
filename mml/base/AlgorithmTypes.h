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

#include <mml/MMLBase.h>

#include <string>
#include <chrono>
#include <cmath>
#include <limits>
#include <optional>

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
		
		/// Algorithm ran to completion but has NO convergence criterion to check
		/// (fixed quadrature rules, direct solvers, closed-form evaluations, geometric
		/// fallbacks). The value is usable; `converged` stays false and estimate fields
		/// that were not computed are NoEstimate. See docs/RESULT_POLICY.md.
		Completed = 7,

		/// A storage or allocation limit prevented further progress.
		ResourceLimitExceeded = 8,

		/// No algorithm has run; a default-constructed result has no usable payload.
		NotComputed = 9,
		
		/// Algorithm-specific failure (see error_message for details)
		AlgorithmSpecificFailure = 100
	};

	/// Convert AlgorithmStatus to human-readable string.
	inline std::string to_string(AlgorithmStatus status) {
		switch (status) {
			case AlgorithmStatus::Success:                return "Success";
			case AlgorithmStatus::MaxIterationsExceeded:  return "MaxIterationsExceeded";
			case AlgorithmStatus::NumericalInstability:   return "NumericalInstability";
			case AlgorithmStatus::SingularMatrix:         return "SingularMatrix";
			case AlgorithmStatus::InvalidInput:           return "InvalidInput";
			case AlgorithmStatus::Stalled:                return "Stalled";
			case AlgorithmStatus::ToleranceUnachievable:  return "ToleranceUnachievable";
			case AlgorithmStatus::Completed:              return "Completed";
			case AlgorithmStatus::ResourceLimitExceeded:  return "ResourceLimitExceeded";
			case AlgorithmStatus::NotComputed:            return "NotComputed";
			case AlgorithmStatus::AlgorithmSpecificFailure: return "AlgorithmSpecificFailure";
			default:                                      return "Unknown";
		}
	}
	/// @deprecated Use to_string() (CONTRACT_POLICY section 11: STL spelling).
	[[deprecated("use to_string()")]]
	inline std::string ToString(AlgorithmStatus status) { return to_string(status); }

	/// True when the result carries a usable value: a checked criterion was met (Success)
	/// or the method has no criterion and ran to completion (Completed).
	[[nodiscard]] constexpr bool IsUsable(AlgorithmStatus status) noexcept {
		return status == AlgorithmStatus::Success || status == AlgorithmStatus::Completed;
	}

	/******************************************************************************/
	/*****                    Estimate Sentinel                               *****/
	/******************************************************************************/

	/// Sentinel for numeric estimate fields (error estimates, achieved tolerances,
	/// residuals) that the algorithm did NOT compute. A literal 0.0 would claim an
	/// exact result; NaN cannot be mistaken for one. Test with HasEstimate().
	inline constexpr Real NoEstimate = std::numeric_limits<Real>::quiet_NaN();

	/// True if an estimate field holds a computed value (i.e. is not NoEstimate).
	[[nodiscard]] inline bool HasEstimate(Real estimate) noexcept {
		return !std::isnan(estimate);
	}

	template<class Type>
	[[nodiscard]] std::enable_if_t<std::is_floating_point_v<Type>, bool>
	HasEstimate(Type estimate) noexcept {
		return !std::isnan(estimate);
	}

	template<class Type>
	[[nodiscard]] bool HasEstimate(const std::optional<Type>& estimate) noexcept {
		return estimate.has_value();
	}

	template<class Type>
	[[nodiscard]] bool HasEstimate(const std::complex<Type>& estimate) noexcept {
		return HasEstimate(estimate.real()) && HasEstimate(estimate.imag());
	}

	template<class Type, int N> class VectorN;
	template<class Type, int N, int M> requires Field<Type> class MatrixNM;

	template<class Type, int N>
	[[nodiscard]] bool HasEstimate(const VectorN<Type, N>& estimate) noexcept {
		for (int i = 0; i < N; ++i)
			if (!HasEstimate(estimate[i])) return false;
		return true;
	}

	template<class Type, int N, int M> requires Field<Type>
	[[nodiscard]] bool HasEstimate(const MatrixNM<Type, N, M>& estimate) noexcept {
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
				if (!HasEstimate(estimate(i, j))) return false;
		return true;
	}

	namespace Detail {
		template<class Scalar, int N>
		struct VectorEstimateShape {
			using ScalarType = Scalar;
			static constexpr int Dimension = N;
		};
		template<class Scalar, int N>
		VectorEstimateShape<Scalar, N> IdentifyVectorEstimate(const VectorN<Scalar, N>*);

		template<class Error, class = void>
		struct EstimateStorage {
			using Type = std::optional<Error>;
			static constexpr bool HasNaNSentinel = false;
			static Type Unknown() { return std::nullopt; }
		};

		template<class Error>
		struct EstimateStorage<Error, std::enable_if_t<std::is_floating_point_v<Error>>> {
			using Type = Error;
			static constexpr bool HasNaNSentinel = true;
			static constexpr Type Unknown() { return std::numeric_limits<Error>::quiet_NaN(); }
		};

		template<class Scalar>
		struct EstimateStorage<std::complex<Scalar>, std::enable_if_t<std::is_floating_point_v<Scalar>>> {
			using Type = std::complex<Scalar>;
			static constexpr bool HasNaNSentinel = true;
			static constexpr Type Unknown() {
				const auto unknown = EstimateStorage<Scalar>::Unknown();
				return Type(unknown, unknown);
			}
		};

		// Deduction through the base also covers Cartesian/spherical/cylindrical vectors.
		template<class Error>
		struct EstimateStorage<Error, std::void_t<
			typename decltype(IdentifyVectorEstimate(static_cast<const Error*>(nullptr)))::ScalarType,
			std::enable_if_t<EstimateStorage<typename decltype(
				IdentifyVectorEstimate(static_cast<const Error*>(nullptr)))::ScalarType>::HasNaNSentinel>>> {
			using Shape = decltype(IdentifyVectorEstimate(static_cast<const Error*>(nullptr)));
			using Type = Error;
			static constexpr bool HasNaNSentinel = true;
			static Type Unknown() {
				Type estimate;
				for (int i = 0; i < Shape::Dimension; ++i)
					estimate[i] = EstimateStorage<typename Shape::ScalarType>::Unknown();
				return estimate;
			}
		};

		template<class Scalar, int N, int M> requires Field<Scalar>
		struct EstimateStorage<MatrixNM<Scalar, N, M>,
			std::enable_if_t<EstimateStorage<Scalar>::HasNaNSentinel>> {
			using Type = MatrixNM<Scalar, N, M>;
			static constexpr bool HasNaNSentinel = true;
			static Type Unknown() {
				Type estimate;
				for (int i = 0; i < N; ++i)
					for (int j = 0; j < M; ++j)
						estimate(i, j) = EstimateStorage<Scalar>::Unknown();
				return estimate;
			}
		};
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
		Real tolerance = Precision::EqualityTolerance;
		
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
	///     Real achieved_tolerance = NoEstimate;
	///     AlgorithmStatus status = AlgorithmStatus::NotComputed;
	///     std::string error_message;
	///     
	///     // === Diagnostics ===
	///     std::string algorithm_name;
	///     double elapsed_time_ms = 0.0;
	///     Real max_residual = NoEstimate;  // Algorithm-specific
	/// };
	struct IterativeResultBase {
		/// True only if a convergence criterion was checked and met (status == Success).
		/// False for Completed results (no criterion exists) and for every failure.
		bool converged = false;
		
		/// Number of iterations actually performed
		int iterations_used = 0;
		
		/// Actual achieved tolerance (may be better than requested)
		/// Interpretation depends on algorithm (residual, error estimate, etc.)
		/// NoEstimate until the algorithm computes it (RESULT_POLICY §3; see HasEstimate()).
		Real achieved_tolerance = NoEstimate;
		
		/// Structured status code for programmatic handling
		AlgorithmStatus status = AlgorithmStatus::NotComputed;
		
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

		/// True when the payload is usable: status is Success or Completed (RESULT_POLICY §2).
		[[nodiscard]] bool IsSuccess() const noexcept { return IsUsable(status); }

		/// Allow `if (result)` style checks; same meaning as IsSuccess().
		explicit operator bool() const noexcept { return IsSuccess(); }
	};

	/******************************************************************************/
	/*****              Non-Iterative Evaluation Contract Types               *****/
	/******************************************************************************/

	/// Policy for handling exceptions raised during evaluation-style algorithms.
	///
	/// Numerical differentiation and similar operations often sit on the boundary
	/// between pure mathematical evaluation and algorithm execution. This policy
	/// lets detailed APIs either preserve the existing exception behavior or map
	/// expected runtime failures into AlgorithmStatus-based results.
	enum class EvaluationExceptionPolicy {
		/// Preserve existing behavior and rethrow exceptions to the caller
		Propagate = 0,

		/// Catch supported exceptions and report them via AlgorithmStatus
		ConvertToStatus = 1
	};

	/// Base configuration for non-iterative evaluation algorithms.
	///
	/// This covers one-shot numerical operations such as differentiation,
	/// interpolation probes, and other evaluation helpers that do not have a
	/// converged/iterations lifecycle but still need configurable diagnostics.
	struct EvaluationConfigBase {
		/// Whether detailed APIs should compute and populate an error estimate
		bool estimate_error = true;

		/// Whether detailed APIs should validate finiteness of intermediate values
		bool check_finite = true;

		/// Whether exceptions should propagate or be converted into status results
		EvaluationExceptionPolicy exception_policy = EvaluationExceptionPolicy::Propagate;
	};

	/// Base result for non-iterative evaluation algorithms.
	///
	/// Unlike IterativeResultBase, this does not model convergence or iteration
	/// counts. It captures the shared diagnostics that still matter for single-shot
	/// numerical operations: structured status, message, timing, and evaluation count.
	struct EvaluationResultBase {
		/// Structured status code for programmatic handling
		AlgorithmStatus status = AlgorithmStatus::NotComputed;

		/// Human-readable error message (empty on success)
		std::string error_message;

		/// Name of the algorithm that produced this result
		std::string algorithm_name;

		/// Wall-clock execution time in milliseconds
		double elapsed_time_ms = 0.0;

		/// Number of function evaluations performed
		int function_evaluations = 0;

		/// True when the payload is usable: status is Success or Completed (RESULT_POLICY §2).
		[[nodiscard]] bool IsSuccess() const noexcept {
			return IsUsable(status);
		}

		/// Allow `if (result)` style checks in user code
		explicit operator bool() const noexcept {
			return IsSuccess();
		}
	};

	/// Generic value + error result for evaluation algorithms.
	///
	/// @tparam TValue Primary output value type
	/// @tparam TError Error estimate type (defaults to TValue)
	template<typename TValue, typename TError = TValue>
	struct EvaluationResult : public EvaluationResultBase {
		/// Primary computed value
		TValue value{};

		/// Unknown until computed: native NaN for floating/fixed-shape estimates,
		/// empty optional for error types without a NaN representation.
		typename Detail::EstimateStorage<TError>::Type error = Detail::EstimateStorage<TError>::Unknown();
	};

	/// Configuration for numerical differentiation algorithms.
	///
	/// Step selection is derivative-specific, so it lives here rather than in the
	/// generic EvaluationConfigBase. A zero step means "choose automatically".
	struct DerivativeConfig : public EvaluationConfigBase {
		/// User-specified step size; 0 means use algorithm-selected default scaling
		Real step = 0.0;
	};

	/// Detailed result for numerical differentiation algorithms.
	///
	/// Extends the generic evaluation result with the effective step size used by
	/// the derivative stencil.
	template<typename TValue, typename TError = TValue>
	struct DerivativeResult : public EvaluationResult<TValue, TError> {
		/// Actual step size used by the derivative evaluation
		Real step_used = 0.0;
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

	/// Create a Completed result: the algorithm has no convergence criterion (fixed rule,
	/// direct method) and ran to the end. converged stays false and achieved_tolerance is
	/// NoEstimate unless the caller supplies a genuinely computed bound afterwards.
	template<typename ResultType>
	ResultType MakeCompletedResult(const std::string& algorithm_name, int iterations = 0) {
		ResultType result;
		result.converged = false;
		result.iterations_used = iterations;
		result.achieved_tolerance = NoEstimate;
		result.status = AlgorithmStatus::Completed;
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

	/// Create a success result for a non-iterative evaluation algorithm.
	///
	/// @tparam ResultType The algorithm-specific evaluation result type
	/// @param algorithm_name Name of the algorithm
	/// @param function_evaluations Number of function evaluations performed
	/// @return ResultType with success fields populated
	template<typename ResultType>
	ResultType MakeEvaluationSuccessResult(const std::string& algorithm_name,
	                                      int function_evaluations = 0) {
		ResultType result;
		result.status = AlgorithmStatus::Success;
		result.error_message.clear();
		result.algorithm_name = algorithm_name;
		result.function_evaluations = function_evaluations;
		return result;
	}

	/// Create a failure result for a non-iterative evaluation algorithm.
	///
	/// @tparam ResultType The algorithm-specific evaluation result type
	/// @param status The failure status code
	/// @param error_message Descriptive error message
	/// @param algorithm_name Name of the algorithm
	/// @param function_evaluations Number of function evaluations performed
	/// @return ResultType with failure fields populated
	template<typename ResultType>
	ResultType MakeEvaluationFailureResult(AlgorithmStatus status,
	                                     const std::string& error_message,
	                                     const std::string& algorithm_name,
	                                     int function_evaluations = 0) {
		ResultType result;
		result.status = status;
		result.error_message = error_message;
		result.algorithm_name = algorithm_name;
		result.function_evaluations = function_evaluations;
		return result;
	}

} // namespace MML

#endif // MML_ALGORITHM_TYPES_H
