///////////////////////////////////////////////////////////////////////////////////////////
// StatisticsBase.h
//
// Base types for Statistics Detailed API
// Provides StatisticsConfig, detailed result types, and execution helper
// following the EvaluationConfigBase/EvaluationResultBase pattern from AlgorithmTypes.h
//
// Requires: AlgorithmTypes.h
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_STATISTICS_BASE_H
#define MML_STATISTICS_BASE_H

#include "core/AlgorithmTypes.h"

namespace MML
{
	namespace Statistics
	{
		/******************************************************************************/
		/*****             Statistics Detailed API Types                          *****/
		/******************************************************************************/

		/// Configuration for statistics detailed APIs
		struct StatisticsConfig : public EvaluationConfigBase {
			// Inherits: estimate_error, check_finite, exception_policy
		};

		/// Detailed result for hypothesis tests
		///
		/// Extends EvaluationResultBase with all hypothesis-test-specific fields.
		struct HypothesisTestDetailedResult : public EvaluationResultBase {
			Real testStatistic{};        ///< Computed test statistic (t, z, chi-square, F, etc.)
			Real pValue{};               ///< P-value for the test
			Real criticalValue{};        ///< Critical value at the given significance level
			bool rejectNull = false;     ///< Whether to reject the null hypothesis
			Real confidenceLevel{};      ///< Confidence level used (e.g., 0.95 for 95%)
			int degreesOfFreedom = -1;   ///< Degrees of freedom (if applicable, -1 otherwise)
			std::string testName;        ///< Descriptive name of the test
		};

		/// Detailed result for confidence intervals
		///
		/// Extends EvaluationResultBase with confidence-interval-specific fields.
		struct ConfidenceIntervalDetailedResult : public EvaluationResultBase {
			Real estimate{};             ///< Point estimate
			Real lowerBound{};           ///< Lower confidence limit
			Real upperBound{};           ///< Upper confidence limit
			Real marginOfError{};        ///< Half-width of interval
			Real confidenceLevel{};      ///< Confidence level (e.g., 0.95 for 95%)
			std::string parameter;       ///< What we're estimating (e.g., "Mean", "Proportion")
		};

		/// Detailed result for rank correlations
		///
		/// Extends EvaluationResultBase with rank-correlation-specific fields.
		struct RankCorrelationDetailedResult : public EvaluationResultBase {
			Real rho{};                  ///< Correlation coefficient (Spearman or Kendall)
			Real zScore{};               ///< z-score for large sample approximation
			int n{};                     ///< Sample size

			/// Check if correlation is significant at given alpha level
			bool IsSignificant(Real alpha = 0.05) const {
				Real criticalZ = (alpha <= 0.01) ? 2.576 : 1.96;
				return std::abs(zScore) > criticalZ;
			}
		};

		/******************************************************************************/
		/*****             Statistics Detail Execution Helper                     *****/
		/******************************************************************************/

		namespace StatisticsDetail
		{
			/// Execute a statistics Detailed operation with timing and exception handling.
			///
			/// @tparam ResultType The detailed result type
			/// @tparam ComputeFn Lambda that fills in result fields
			/// @param algorithm_name Name of the algorithm for diagnostics
			/// @param config Configuration controlling exception policy
			/// @param compute Lambda receiving ResultType& to populate
			template<typename ResultType, typename ComputeFn>
			ResultType ExecuteStatisticsDetailed(const char* algorithm_name,
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
				catch (const std::runtime_error& ex) {
					return MakeEvaluationFailureResult<ResultType>(
						AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
				}
				catch (const std::exception& ex) {
					return MakeEvaluationFailureResult<ResultType>(
						AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
				}
			}
		} // namespace StatisticsDetail

	} // namespace Statistics
} // namespace MML

#endif // MML_STATISTICS_BASE_H
