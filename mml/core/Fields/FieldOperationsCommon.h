///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FieldOperationsCommon.h                                             ///
///  Description: Shared configuration and helpers for field differential operations  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_FIELD_OPERATIONS_COMMON_H
#define MML_FIELD_OPERATIONS_COMMON_H

#include <mml/MMLBase.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Geometry/Geometry.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Function.h>
#include <mml/base/Tensor.h>

#include <mml/core/Derivation.h>
#include <mml/core/MetricTensor.h>
#include <mml/MMLSingularityHandling.h>
#include <mml/base/AlgorithmTypes.h>

namespace MML
{
	/******************************************************************************/
	/*****               Field Operation Configuration                       *****/
	/******************************************************************************/

	/// Configuration for field differential operations (gradient, divergence, curl, Laplacian).
	///
	/// Extends EvaluationConfigBase with derivative order selection.
	/// Field operations are non-iterative evaluations that internally use numerical
	/// differentiation, so derivative accuracy is the primary tunable parameter.
	struct FieldOperationConfig : public EvaluationConfigBase {
		/// Derivative order for internal numerical differentiation.
		/// 0 = use library default (NDer4), or 1, 2, 4, 6, 8.
		int derivative_order = 0;
	};

	/******************************************************************************/
	/*****               Field Operation Detail Helpers                      *****/
	/******************************************************************************/

	namespace FieldOperationDetail
	{
		/// Execute a field operation with structured result reporting.
		///
		/// Wraps the raw computation in timing, finiteness checking, and
		/// exception handling according to the configured policy.
		///
		/// @tparam ResultType  EvaluationResult specialization for the output
		/// @tparam ComputeFn   Lambda (int& func_evals) -> void that populates result.value/error
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteFieldDetailed(const char* algorithm_name,
		                                const FieldOperationConfig& config,
		                                ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;

				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);

				int func_evals = 0;
				compute(result, func_evals);
				result.function_evaluations = func_evals;

				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const DomainError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const NumericInputError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const NumericalMethodError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
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

		/// Dispatch derivative-order selection for DerivePartialAll.
		/// Returns the gradient and optionally populates per-component error estimates.
		template<int N>
		VectorN<Real, N> DispatchGradient(const IScalarFunction<N>& f, const VectorN<Real, N>& pos,
		                                  int order, bool want_error, VectorN<Real, N>* error,
		                                  int& func_evals)
		{
			// Function pointer type for NDerKPartialByAll(f, pos, error*)
			using DeriveFn = VectorN<Real, N>(*)(const IScalarFunction<N>&,
			                                     const VectorN<Real, N>&,
			                                     VectorN<Real, N>*);

			DeriveFn deriveAll = nullptr;
			int stencil = 0;

			switch (order) {
			case 1:  deriveAll = &Derivation::template NDer1PartialByAll<N>; stencil = 2; break;
			case 2:  deriveAll = &Derivation::template NDer2PartialByAll<N>; stencil = 3; break;
			case 0:  // fall through to default (NDer4)
			case 4:  deriveAll = &Derivation::template NDer4PartialByAll<N>; stencil = 5; break;
			case 6:  deriveAll = &Derivation::template NDer6PartialByAll<N>; stencil = 7; break;
			case 8:  deriveAll = &Derivation::template NDer8PartialByAll<N>; stencil = 9; break;
			default:
				throw ArgumentError("FieldOperation: derivative_order must be 0, 1, 2, 4, 6, or 8");
			}

			func_evals = N * (want_error ? stencil + 1 : stencil);
			return deriveAll(f, pos, want_error ? error : nullptr);
		}
	} // namespace FieldOperationDetail
}

#endif