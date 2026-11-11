#ifndef MML_INTEGRATION_BASE_H
#define MML_INTEGRATION_BASE_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>

#include <ostream>
#include <string>

namespace MML
{
	/// Configuration for integration Detailed APIs
	struct IntegrationConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
	};

	/******************************************************************************/
	/*****             The one integration result type                         *****/
	/******************************************************************************/

	/// Result of a numerical integration (RESULT_POLICY: one type for every integration layer).
	///
	/// - `converged`        true only if the method's convergence criterion was checked and met
	/// - `error_estimated`  false for fixed rules (Gauss10, Gauss-Legendre, midpoint sums, ...)
	/// - `error_estimate`   a computed bound, or NoEstimate when none exists (HasEstimate())
	/// - `status`           Success (criterion met) / Completed (fixed rule) / MaxIterationsExceeded ...
	/// - `IsSuccess()`      the value is usable (Success or Completed)
	///
	/// Plain integration functions return the value and throw ConvergenceError when the result is
	/// not usable; the `...Detailed` functions return this struct. Conversion to the value is
	/// explicit (`static_cast<TValue>(r)`) or via `.value` so that discarding diagnostics is
	/// visible at the call site.
	template<typename TValue>
	struct [[nodiscard]] BasicIntegrationResult : public EvaluationResultBase {
		TValue value{};                       ///< Computed integral value
		Real error_estimate = NoEstimate;     ///< Estimated absolute error (NoEstimate until computed)
		int iterations = 0;                   ///< Iterations / refinement levels performed (1 for fixed rules)
		bool converged = false;               ///< Criterion checked and met
		bool error_estimated = false;         ///< An error estimate exists (false for fixed rules)

		BasicIntegrationResult() = default;

		/// Field-wise constructor used by the integrator cores. The status is derived:
		/// converged -> Success; no estimate -> Completed; otherwise MaxIterationsExceeded.
		BasicIntegrationResult(TValue val, Real err, int iter, bool conv, bool errEstimated = true)
			: value(val), error_estimate(errEstimated ? err : NoEstimate), iterations(iter)
			, converged(conv), error_estimated(errEstimated)
		{
			if (conv)                 status = AlgorithmStatus::Success;
			else if (!errEstimated)   status = AlgorithmStatus::Completed;
			else                      status = AlgorithmStatus::MaxIterationsExceeded;
		}

		/// Result of a fixed quadrature rule: nothing was checked, there is no estimate (RESULT_POLICY §2).
		static BasicIntegrationResult FixedRule(TValue val, int function_evals = 0)
		{
			BasicIntegrationResult r(val, NoEstimate, 1, /*converged*/ false, /*error_estimated*/ false);
			r.function_evaluations = function_evals;
			return r;
		}

		/// Mark an existing result as produced by a fixed rule.
		void MarkFixedRule(int function_evals)
		{
			error_estimate = NoEstimate;
			error_estimated = false;
			converged = false;
			iterations = 1;
			status = AlgorithmStatus::Completed;
			function_evaluations = function_evals;
		}

		/// True if the value may be used: the criterion was met, or the method has none.
		/// A checked-and-missed criterion (MaxIterationsExceeded, ...) is not usable.
		[[nodiscard]] bool IsUsable() const noexcept { return IsUsable_(status) && (status == AlgorithmStatus::Completed || converged || !error_estimated); }

		/// Returns the value, or throws ConvergenceError when the result is not usable.
		/// This is what the plain (non-Detailed) integration functions return.
		TValue ValueOrThrow(const char* context) const
		{
			if (IsUsable())
				return value;
			if (!error_message.empty())
				throw ConvergenceError(std::string(context) + ": " + error_message, iterations, error_estimate);
			throw ConvergenceError(std::string(context) + " did not converge after " + std::to_string(iterations) +
				" refinements (error estimate " + std::to_string(error_estimate) + ")", iterations, error_estimate);
		}

		/// Explicit conversion to the integral value.
		/// @warning Discards error_estimate, iterations, converged and status - check them first.
		explicit operator TValue() const { return value; }

		/// Prints the value; a result whose convergence criterion was checked and not met is marked.
		friend std::ostream& operator<<(std::ostream& out, const BasicIntegrationResult& r)
		{
			out << r.value;
			if (r.status != AlgorithmStatus::Completed && r.error_estimated && !r.converged)
				out << " (not converged)";
			return out;
		}

	private:
		static bool IsUsable_(AlgorithmStatus s) noexcept { return MML::IsUsable(s); }
	};

	/// Real-valued integration result (1D, multi-dimensional, path, surface, Monte Carlo base)
	using IntegrationResult = BasicIntegrationResult<Real>;

	/// @deprecated Same type as IntegrationResult since the result unification (wewi.3.3).
	using IntegrationDetailedResult = IntegrationResult;
	/******************************************************************************/
	/*****             Integration Detail Execution Helper                    *****/
	/******************************************************************************/

	namespace IntegrationDetail
	{
		/// Copy the integration payload and diagnostics of a core result into a Detailed result
		/// while keeping the wrapper's algorithm_name / elapsed_time_ms, and derive the status
		/// per RESULT_POLICY: converged -> Success, no estimate -> Completed, checked and missed ->
		/// MaxIterationsExceeded with an actionable message.
		template<typename TValue>
		inline void Adopt(BasicIntegrationResult<TValue>& out,
		                  const BasicIntegrationResult<TValue>& r,
		                  const char* algorithm_name)
		{
			out.value = r.value;
			out.iterations = r.iterations;
			out.converged = r.converged;
			out.error_estimated = r.error_estimated;
			out.error_estimate = r.error_estimated ? r.error_estimate : NoEstimate;
			if (r.function_evaluations > 0)
				out.function_evaluations = r.function_evaluations;
			if (!IsUsable(r.status)) {
				out.status = r.status;
				out.error_message = r.error_message;
				return;
			}
			if (r.status == AlgorithmStatus::Completed) {
				out.status = AlgorithmStatus::Completed;
			} else if (r.converged) {
				out.status = AlgorithmStatus::Success;
			} else if (!r.error_estimated) {
				out.status = AlgorithmStatus::Completed;
			} else {
				out.status = AlgorithmStatus::MaxIterationsExceeded;
				out.error_message = std::string(algorithm_name) + " did not converge after " +
					std::to_string(r.iterations) + " refinements (error estimate " +
					std::to_string(r.error_estimate) + ")";
			}
		}

		/// Execute an integration Detailed operation with timing and exception handling.
		///
		/// @tparam ResultType The detailed result type (e.g. IntegrationDetailedResult)
		/// @tparam ComputeFn Lambda that fills in result fields
		/// @param algorithm_name Name of the algorithm for diagnostics
		/// @param config Configuration controlling exception policy
		/// @param compute Lambda receiving ResultType& to populate
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteIntegrationDetailed(const char* algorithm_name,
		                                     const EvaluationConfigBase& config,
		                                     ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;
				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
				compute(result);
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const NumericalMethodError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
			}
			catch (const ResourceLimitError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::ResourceLimitExceeded, ex.what(), algorithm_name);
			}
			catch (const std::invalid_argument& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}
	} // namespace IntegrationDetail
}

#endif  // MML_INTEGRATION_BASE_H
