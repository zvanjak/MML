///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        NumericValidation.h                                                 ///
///  Description: Shared numeric validation helpers for algorithms                    ///
///               Validates finiteness of bounds, tolerances, and function values     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_NUMERIC_VALIDATION_H
#define MML_NUMERIC_VALIDATION_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include <cmath>
#include <string>

namespace MML {

	///////////////////////////////////////////////////////////////////////////
	///                    NUMERIC VALIDATION HELPERS                       ///
	///////////////////////////////////////////////////////////////////////////
	//
	// These helpers provide consistent input validation across all numerical
	// algorithms. They check for non-finite values (NaN, Inf) and throw
	// NumericInputError with descriptive messages.
	//
	// Usage:
	//   ValidateFinite(x, "initial guess");
	//   ValidateBounds(a, b, "bisection");
	//   ValidateTolerance(tol, "Newton iteration");
	//   ValidateFunctionValue(f(x), "f(x)");
	//
	///////////////////////////////////////////////////////////////////////////

	/// @brief Check if a value is finite (not NaN or Inf)
	/// @param x Value to check
	/// @return true if x is finite, false if NaN or Inf
	inline bool IsFinite(Real x) noexcept {
		return std::isfinite(x);
	}

	/// @brief Check if a value is finite (for generic types)
	/// @tparam T Numeric type
	/// @param x Value to check
	/// @return true if x is finite, false if NaN or Inf
	template<typename T>
	inline bool IsFiniteValue(T x) noexcept {
		return std::isfinite(x);
	}

	/// @brief Validate that a value is finite, throw if not
	/// @param x Value to check
	/// @param name Name of the parameter for error message
	/// @throws NumericInputError if x is NaN or Inf
	inline void ValidateFinite(Real x, const char* name) {
		if (!std::isfinite(x)) {
			throw NumericInputError(std::string(name) + " must be finite (got " + 
				(std::isnan(x) ? "NaN" : "Inf") + ")");
		}
	}

	/// @brief Validate that a value is finite, throw if not (string version)
	/// @param x Value to check
	/// @param name Name of the parameter for error message
	/// @throws NumericInputError if x is NaN or Inf
	inline void ValidateFinite(Real x, const std::string& name) {
		ValidateFinite(x, name.c_str());
	}

	/// @brief Validate bounds for 1D algorithms (root finding, optimization, integration)
	/// @param a Lower bound
	/// @param b Upper bound
	/// @param context Name of the algorithm context for error message
	/// @throws NumericInputError if either bound is NaN or Inf
	inline void ValidateBounds(Real a, Real b, const char* context = "algorithm") {
		if (!std::isfinite(a)) {
			throw NumericInputError(std::string(context) + ": lower bound 'a' must be finite (got " +
				(std::isnan(a) ? "NaN" : "Inf") + ")");
		}
		if (!std::isfinite(b)) {
			throw NumericInputError(std::string(context) + ": upper bound 'b' must be finite (got " +
				(std::isnan(b) ? "NaN" : "Inf") + ")");
		}
	}

	/// @brief Validate bounds (string version)
	inline void ValidateBounds(Real a, Real b, const std::string& context) {
		ValidateBounds(a, b, context.c_str());
	}

	/// @brief Validate tolerance parameter
	/// @param tol Tolerance value
	/// @param context Name of the algorithm context for error message
	/// @throws NumericInputError if tolerance is not positive and finite
	inline void ValidateTolerance(Real tol, const char* context = "algorithm") {
		if (!std::isfinite(tol) || tol <= 0) {
			throw NumericInputError(std::string(context) + ": tolerance must be positive and finite (got " +
				std::to_string(tol) + ")");
		}
	}

	/// @brief Validate tolerance (string version)
	inline void ValidateTolerance(Real tol, const std::string& context) {
		ValidateTolerance(tol, context.c_str());
	}

	/// @brief Validate a function value returned during iteration
	/// @param fval Function value to check
	/// @param context Description of where this was evaluated
	/// @throws NumericInputError if function returned NaN or Inf
	inline void ValidateFunctionValue(Real fval, const char* context = "function evaluation") {
		if (!std::isfinite(fval)) {
			throw NumericInputError(std::string(context) + " returned non-finite value (" +
				(std::isnan(fval) ? "NaN" : "Inf") + ")");
		}
	}

	/// @brief Validate function value (string version)
	inline void ValidateFunctionValue(Real fval, const std::string& context) {
		ValidateFunctionValue(fval, context.c_str());
	}

	/// @brief Check if a function value is finite (non-throwing version)
	/// @param fval Function value to check
	/// @return true if fval is finite, false if NaN or Inf
	inline bool IsFunctionValueValid(Real fval) noexcept {
		return std::isfinite(fval);
	}

	/// @brief Validate maximum iterations parameter
	/// @param maxIter Maximum iterations value
	/// @param context Name of the algorithm context for error message
	/// @throws NumericInputError if maxIter is not positive
	inline void ValidateMaxIterations(int maxIter, const char* context = "algorithm") {
		if (maxIter <= 0) {
			throw NumericInputError(std::string(context) + ": max_iterations must be positive (got " +
				std::to_string(maxIter) + ")");
		}
	}

	/// @brief Validate step size for numerical differentiation/integration
	/// @param h Step size
	/// @param context Name of the algorithm context for error message
	/// @throws NumericInputError if h is not positive and finite
	inline void ValidateStepSize(Real h, const char* context = "algorithm") {
		if (!std::isfinite(h) || h <= 0) {
			throw NumericInputError(std::string(context) + ": step size must be positive and finite (got " +
				std::to_string(h) + ")");
		}
	}

} // namespace MML

#endif // MML_NUMERIC_VALIDATION_H
