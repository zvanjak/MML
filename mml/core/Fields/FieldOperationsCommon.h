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
#include <type_traits>

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
		/// 0 = provider/default (NDer4), or 1, 2, 4, 6, 8.
		/// Laplacians support only 0, 2, 4 (including their first partials).
		int derivative_order = 0;
	};

	/******************************************************************************/
	/*****               Field Operation Detail Helpers                      *****/
	/******************************************************************************/

	namespace FieldOperationDetail
	{
		template<typename Value>
		void SetUnavailable(Value& value) {
			if constexpr (std::is_arithmetic_v<Value>)
				value = NoEstimate;
			else
				for (int i = 0; i < value.size(); ++i)
					value[i] = NoEstimate;
		}

		template<typename Value>
		void CheckFinite(const Value& value, bool enabled) {
			if (!enabled) return;
			if constexpr (std::is_arithmetic_v<Value>) {
				if (!std::isfinite(value))
					throw NumericalMethodError("Field operation - nonfinite sample or arithmetic result");
			} else
				for (int i = 0; i < value.size(); ++i)
					CheckFinite(value[i], true);
		}

		inline void ValidateOrder(int order, bool secondDerivatives) {
			if (secondDerivatives) {
				if (order != 0 && order != 2 && order != 4)
					throw ArgumentError("Laplacian derivative_order must be 0, 2, or 4");
			} else if (order != 0 && order != 1 && order != 2 && order != 4 && order != 6 && order != 8)
				throw ArgumentError("Field operation derivative_order must be 0, 1, 2, 4, 6, or 8");
		}

		template<typename Value>
		void CheckEstimate(const Value& value) {
			if constexpr (std::is_arithmetic_v<Value>) {
				if (!std::isnan(value)) {
					CheckFinite(value, true);
					if (value < 0)
						throw NumericalMethodError("Field operation - negative error estimate");
				}
			} else
				for (int i = 0; i < value.size(); ++i)
					CheckEstimate(value[i]);
		}

		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteFieldDetailed(const char* algorithm_name,
		                                const FieldOperationConfig& config,
		                                ComputeFn&& compute, bool secondDerivatives = false)
		{
			AlgorithmTimer timer;
			int func_evals = 0;
			ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
			SetUnavailable(result.error);
			auto execute = [&]() {
				ValidateOrder(config.derivative_order, secondDerivatives);
				compute(result, func_evals);
				CheckFinite(result.value, config.check_finite);
				if (config.check_finite)
					CheckEstimate(result.error);
				result.function_evaluations = func_evals;
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};
			auto failure = [&](AlgorithmStatus status, const std::exception& ex) {
				auto failed = MakeEvaluationFailureResult<ResultType>(status, ex.what(), algorithm_name, func_evals);
				SetUnavailable(failed.value);
				SetUnavailable(failed.error);
				failed.elapsed_time_ms = timer.elapsed_ms();
				return failed;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();
			try {
				return execute();
			}
			catch (const DomainError& ex) {
				return failure(AlgorithmStatus::InvalidInput, ex);
			}
			catch (const NumericInputError& ex) {
				return failure(AlgorithmStatus::InvalidInput, ex);
			}
			catch (const ArgumentError& ex) {
				return failure(AlgorithmStatus::InvalidInput, ex);
			}
			catch (const NumericalMethodError& ex) {
				return failure(AlgorithmStatus::NumericalInstability, ex);
			}
			catch (const std::invalid_argument& ex) {
				return failure(AlgorithmStatus::InvalidInput, ex);
			}
			catch (const std::exception& ex) {
				return failure(AlgorithmStatus::AlgorithmSpecificFailure, ex);
			}
		}

		template<int N, typename Function, typename... Indices>
		Real FirstPartial(const Function& f, const VectorN<Real, N>& pos, int order, Real* error, Indices... indices) {
			switch (order) {
			case 1: return Derivation::NDer1Partial<N>(f, indices..., pos, error);
			case 2: return Derivation::NDer2Partial<N>(f, indices..., pos, error);
			case 0:
			case 4: return Derivation::NDer4Partial<N>(f, indices..., pos, error);
			case 6: return Derivation::NDer6Partial<N>(f, indices..., pos, error);
			case 8: return Derivation::NDer8Partial<N>(f, indices..., pos, error);
			default: throw ArgumentError("Unsupported field derivative order");
			}
		}

		// Scoped adapters let the plain coordinate formulas share configured derivatives
		// without treating derivative-provider calls as field operator() evaluations.
		template<int N>
		class ScalarSampling final : public IScalarFunction<N>, public IScalarDerivativeProvider<N> {
			const IScalarFunction<N>& _source;
			const FieldOperationConfig& _config;
			int& _evaluations;
			const IScalarDerivativeProvider<N>* _provider;
		public:
			mutable VectorN<Real, N> firstErrors{};
			mutable MatrixNM<Real, N, N> secondErrors{};
			ScalarSampling(const IScalarFunction<N>& source, const VectorN<Real, N>& pos,
			               const FieldOperationConfig& config, int& evaluations)
				: _source(source), _config(config), _evaluations(evaluations),
				  _provider(config.derivative_order == 0 ? dynamic_cast<const IScalarDerivativeProvider<N>*>(&source) : nullptr) {
				Singularity::RequireFiniteCoordinates<N>(pos, "Field operation");
			}
			bool HasEstimate() const { return _config.estimate_error && !_provider; }
			Real operator()(const VectorN<Real, N>& p) const override {
				++_evaluations;
				const Real value = _source(p);
				CheckFinite(value, _config.check_finite);
				return value;
			}
			Real EvaluatePartial(const VectorN<Real, N>& p, int axis) const override {
				Real error = 0;
				const Real value = _provider ? _provider->EvaluatePartial(p, axis)
					: FirstPartial<N>(*this, p, _config.derivative_order, HasEstimate() ? &error : nullptr, axis);
				firstErrors[axis] = std::abs(error);
				CheckFinite(value, _config.check_finite);
				CheckFinite(error, _config.check_finite);
				return value;
			}
			Real EvaluateSecondPartial(const VectorN<Real, N>& p, int a, int b) const override {
				Real error = 0;
				Real value;
				if (_provider) value = _provider->EvaluateSecondPartial(p, a, b);
				else if (_config.derivative_order == 2)
					value = Derivation::NSecDer2Partial<N>(*this, a, b, p, HasEstimate() ? &error : nullptr);
				else
					value = Derivation::NSecDer4Partial<N>(*this, a, b, p, HasEstimate() ? &error : nullptr);
				secondErrors(a,b) = std::abs(error);
				CheckFinite(value, _config.check_finite);
				CheckFinite(error, _config.check_finite);
				return value;
			}
		};

		template<int N>
		class VectorSampling final : public IVectorFunction<N>, public IVectorDerivativeProvider<N> {
			const IVectorFunction<N>& _source;
			const FieldOperationConfig& _config;
			int& _evaluations;
			const IVectorDerivativeProvider<N>* _provider;
		public:
			mutable MatrixNM<Real, N, N> errors{};
			VectorSampling(const IVectorFunction<N>& source, const VectorN<Real, N>& pos,
			               const FieldOperationConfig& config, int& evaluations)
				: _source(source), _config(config), _evaluations(evaluations),
				  _provider(config.derivative_order == 0 ? dynamic_cast<const IVectorDerivativeProvider<N>*>(&source) : nullptr) {
				Singularity::RequireFiniteCoordinates<N>(pos, "Field operation");
			}
			bool HasEstimate() const { return _config.estimate_error && !_provider; }
			VectorN<Real, N> operator()(const VectorN<Real, N>& p) const override {
				++_evaluations;
				const auto value = _source(p);
				CheckFinite(value, _config.check_finite);
				return value;
			}
			Real EvaluateComponentPartial(const VectorN<Real, N>& p, int component, int axis) const override {
				Real error = 0;
				const Real value = _provider ? _provider->EvaluateComponentPartial(p, component, axis)
					: FirstPartial<N>(*this, p, _config.derivative_order, HasEstimate() ? &error : nullptr, component, axis);
				errors(component,axis) = std::abs(error);
				CheckFinite(value, _config.check_finite);
				CheckFinite(error, _config.check_finite);
				return value;
			}
		};
	} // namespace FieldOperationDetail
}

#endif