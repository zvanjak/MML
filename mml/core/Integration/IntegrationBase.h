#ifndef MML_INTEGRATION_BASE_H
#define MML_INTEGRATION_BASE_H

#include "MMLBase.h"
#include "core/AlgorithmTypes.h"

namespace MML
{
	/// Result structure for numerical integration algorithms
	/// Provides convergence status and diagnostics
	/// @note For production code, check error_estimate and converged fields!
	struct IntegrationResult {
		Real value;					 ///< Computed integral value
		Real error_estimate; ///< Estimated absolute error
		int iterations;			 ///< Number of iterations/refinements performed
		bool converged;			 ///< True if convergence criteria met

		/// Implicit conversion to Real for convenience
		/// @note Discards error_estimate, iterations, and converged fields
		operator Real() const { return value; }

		/// Constructor for easy initialization
		IntegrationResult(Real val = 0.0, Real err = 0.0, int iter = 0, bool conv = true)
				: value(val)
				, error_estimate(err)
				, iterations(iter)
				, converged(conv) {}
	};

	/******************************************************************************/
	/*****             Integration Detailed API Types                         *****/
	/******************************************************************************/

	/// Configuration for integration detailed APIs
	struct IntegrationConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
		// No additional integration-specific parameters needed
	};

	/// Detailed result for integration algorithms
	///
	/// Extends EvaluationResultBase with integration-specific fields:
	/// value, error estimate, iteration count, and convergence flag.
	struct IntegrationDetailedResult : public EvaluationResultBase {
		Real value = 0.0;              ///< Computed integral value
		Real error_estimate = 0.0;     ///< Estimated absolute error
		int iterations = 0;            ///< Number of iterations/refinements performed
		bool converged = false;        ///< True if convergence criteria met
	};

	/******************************************************************************/
	/*****             Integration Detail Execution Helper                    *****/
	/******************************************************************************/

	namespace IntegrationDetail
	{
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
