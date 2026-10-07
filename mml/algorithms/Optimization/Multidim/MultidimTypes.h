///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Multidim/MultidimTypes.h                                             ///
///  Description: Shared types, configuration, validation helpers, and result objects       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_OPTIMIZATION_MULTIDIM_TYPES_H
#define MML_OPTIMIZATION_MULTIDIM_TYPES_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/interfaces/IFunction.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace MML::Optimization {

	///                     MultidimOptimizationError                       ///
	///////////////////////////////////////////////////////////////////////////
	class MultidimOptimizationError : public std::runtime_error {
	public:
		explicit MultidimOptimizationError(const std::string& message)
			: std::runtime_error("MultidimOptimizationError: " + message) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///                  MultidimOptimizationInputError                     ///
	///////////////////////////////////////////////////////////////////////////
	/// @brief Exception for invalid inputs to multidimensional optimization routines
	class MultidimOptimizationInputError : public std::domain_error {
	public:
		explicit MultidimOptimizationInputError(const std::string& message)
			: std::domain_error("MultidimOptimizationInputError: " + message) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///               Input Validation Helper Functions                     ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Check if a vector contains any NaN or Inf values
	/// @tparam N Dimension of vector
	/// @param v Vector to check
	/// @return true if all elements are finite
	template<int N>
	inline bool IsVectorFinite(const VectorN<Real, N>& v) {
		for (int i = 0; i < N; ++i) {
			if (!std::isfinite(v[i]))
				return false;
		}
		return true;
	}

	/// @brief Validate that a vector contains no NaN or Inf values
	/// @tparam N Dimension of vector
	/// @param v Vector to check
	/// @param context Description for error message
	/// @throws MultidimOptimizationInputError if vector contains NaN or Inf
	template<int N>
	inline void ValidateVectorFinite(const VectorN<Real, N>& v, const char* context = "vector") {
		for (int i = 0; i < N; ++i) {
			if (!std::isfinite(v[i])) {
				throw MultidimOptimizationInputError(std::string(context) + ": component " + 
					std::to_string(i) + " is " + (std::isnan(v[i]) ? "NaN" : "Inf"));
			}
		}
	}

	/// @brief Validate tolerance parameter
	/// @param tol Tolerance to check
	/// @param context Description for error message
	/// @throws MultidimOptimizationInputError if tolerance is invalid
	inline void ValidateMultidimTolerance(Real tol, const char* context = "optimization") {
		if (!std::isfinite(tol) || tol <= 0) {
			throw MultidimOptimizationInputError(std::string(context) + ": tolerance must be positive and finite, got " +
				std::to_string(tol));
		}
	}

	/// @brief Validate a function value returned from evaluation
	/// @param fval Function value to check
	/// @param context Description for error message
	/// @throws MultidimOptimizationInputError if function value is NaN or Inf
	inline void ValidateMultidimFunctionValue(Real fval, const char* context = "function evaluation") {
		if (!std::isfinite(fval)) {
			throw MultidimOptimizationInputError(std::string(context) + ": returned " +
				(std::isnan(fval) ? "NaN" : "Inf"));
		}
	}

	///////////////////////////////////////////////////////////////////////////
	///                  MultidimMinimizationResult                        ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Result structure for multidimensional minimization
     */
	struct MultidimMinimizationResult {
		Vector<Real> xmin; ///< Location of minimum
		Real fmin;		   ///< Function value at minimum
		int iterations;	   ///< Number of iterations (or function evaluations)
		bool converged;	   ///< True if converged within tolerance
		
		// Enhanced diagnostic fields (API Standardization Phase 4)
		std::string algorithm_name;  ///< Name of the algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;  ///< Algorithm termination status
		std::string error_message;   ///< Error message if failed
		double elapsed_time_ms = 0;  ///< Execution time in milliseconds
		int function_evaluations = 0; ///< Number of function evaluations

		MultidimMinimizationResult()
			: fmin(0.0)
			, iterations(0)
			, converged(false) {}

		MultidimMinimizationResult(const Vector<Real>& x, Real f, int iter, bool conv)
			: xmin(x)
			, fmin(f)
			, iterations(iter)
			, converged(conv)
			, function_evaluations(iter) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///                  MultidimOptimizationConfig                        ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Configuration for multidimensional optimization algorithms
     * 
     * Standardized config object following API Standardization Phase 4 pattern.
     */
	struct MultidimOptimizationConfig {
		Real tolerance = 1e-8;           ///< Convergence tolerance
		int max_iterations = 5000;       ///< Maximum iterations
		bool verbose = false;            ///< Enable verbose output
		Real initial_delta = 1.0;        ///< Initial simplex/step size
		int lbfgs_memory_size = 10;      ///< L-BFGS: number of correction pairs (typical: 3-20)

		/// Default constructor
		MultidimOptimizationConfig() = default;

		/// Constructor with key parameters
		MultidimOptimizationConfig(Real tol, int max_iter = 5000)
			: tolerance(tol)
			, max_iterations(max_iter) {}

		/// Factory: High precision configuration
		static MultidimOptimizationConfig HighPrecision() {
			MultidimOptimizationConfig cfg;
			cfg.tolerance = 1e-12;
			cfg.max_iterations = 10000;
			return cfg;
		}

		/// Factory: Fast configuration (lower precision, fewer iterations)
		static MultidimOptimizationConfig Fast() {
			MultidimOptimizationConfig cfg;
			cfg.tolerance = 1e-4;
			cfg.max_iterations = 500;
			return cfg;
		}

		/// Factory: Large-scale configuration (for L-BFGS with many variables)
		static MultidimOptimizationConfig LargeScale(int memory_size = 20) {
			MultidimOptimizationConfig cfg;
			cfg.tolerance = 1e-8;
			cfg.max_iterations = 10000;
			cfg.lbfgs_memory_size = memory_size;
			return cfg;
		}
	};


} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_MULTIDIM_TYPES_H
