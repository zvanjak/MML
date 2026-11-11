///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationRealFunction.h                                            ///
///  Description: Numerical derivatives of real-valued functions f:R->R               ///
///               Forward, backward, central differences, higher order                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_REAL_FUNCTION_H
#define MML_DERIVATION_REAL_FUNCTION_H

#include <mml/MMLBase.h>
#include <algorithm>

#include <mml/base/AlgorithmTypes.h>
#include <mml/MMLNumericValidation.h>
#include <mml/base/Function.h>

#include "DerivationBase.h"
#include "FirstDerivativeStencil.h"

#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML
{
	namespace Derivation
	{
		namespace Detail
		{
			template<RealFunctionCallable Function, class AdapterCall>
			decltype(auto) WithRealFunctionAdapter(Function&& function, AdapterCall&& call)
				requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
			{
				MML::Detail::RealFunctionCallableAdapter<Function> adapter(function);
				return std::invoke(std::forward<AdapterCall>(call), adapter);
			}

			inline Real ResolveDerivativeStep(const DerivativeConfig& config, Real default_step, Real x)
			{
				return config.step != static_cast<MML::Real>(0.0) ? config.step : ScaleStep(default_step, x);
			}

			inline Real RequireFiniteRichardson(Real value)
			{
				if (!std::isfinite(value))
					throw NumericalMethodError("Richardson derivative: non-finite sample or arithmetic");
				return value;
			}

			template<int DerivativeOrder>
			Real EvaluateRichardsonDerivative(const IRealFunction& f, Real x, Real h,
				Real* error, int maxIterations, Real shrink)
			{
				ValidateFinite(x, "Richardson derivative center");
				ValidateStepSize(h, "Richardson derivative");
				ValidateMaxIterations(maxIterations, "Richardson derivative");
				const Real shrinkSquared = shrink * shrink;
				if (!std::isfinite(shrink) || shrink <= 1 || !std::isfinite(shrinkSquared))
					throw NumericInputError("Richardson derivative: shrink factor and its square must be finite, factor > 1");

				std::vector<std::vector<Real>> table(maxIterations, std::vector<Real>(maxIterations));
				std::optional<Real> centerValue;
				auto estimate = [&](Real step) {
					const Real plus = x + step, minus = x - step;
					if (!std::isfinite(step) || step <= 0 || !std::isfinite(plus) || !std::isfinite(minus) ||
						!(minus < x && x < plus))
						throw NumericalMethodError("Richardson derivative: stencil cannot advance representable coordinates");

					if constexpr (DerivativeOrder == 3) {
						const Real plus2 = x + 2 * step, minus2 = x - 2 * step;
						if (!std::isfinite(plus2) || !std::isfinite(minus2) || !(minus2 < minus && plus < plus2))
							throw NumericalMethodError("Richardson derivative: outer stencil coordinates collapse or overflow");
						const Real fPlus2 = RequireFiniteRichardson(f(plus2));
						const Real fPlus = RequireFiniteRichardson(f(plus));
						const Real fMinus = RequireFiniteRichardson(f(minus));
						const Real fMinus2 = RequireFiniteRichardson(f(minus2));
						return RequireFiniteRichardson(
							((fPlus2 - fMinus2) - 2 * (fPlus - fMinus)) / step / step / step / 2);
					} else {
						if constexpr (DerivativeOrder == 2)
							if (!centerValue) centerValue = RequireFiniteRichardson(f(x));
						const Real fPlus = RequireFiniteRichardson(f(plus));
						const Real fMinus = RequireFiniteRichardson(f(minus));
						if constexpr (DerivativeOrder == 1)
							return RequireFiniteRichardson((fPlus - fMinus) / step / 2);
						else
							return RequireFiniteRichardson(((fPlus - *centerValue) + (fMinus - *centerValue)) / step / step);
					}
				};

				Real step = h;
				table[0][0] = estimate(step);
				Real answer = table[0][0];
				Real bestError = std::numeric_limits<Real>::max();
				for (int i = 1; i < maxIterations; ++i) {
					const Real refinedStep = step / shrink;
					if (refinedStep >= step)
						throw NumericalMethodError("Richardson derivative: step refinement cannot advance");
					step = refinedStep;
					table[0][i] = estimate(step);
					Real factor = shrinkSquared;
					for (int j = 1; j <= i; ++j) {
						const Real difference = RequireFiniteRichardson(table[j - 1][i] - table[j - 1][i - 1]);
						// NR dfridr extrapolation, rearranged to avoid factor * table-entry overflow.
						table[j][i] = RequireFiniteRichardson(table[j - 1][i] + difference / (factor - 1));
						const Real candidateError = std::max(
							RequireFiniteRichardson(std::abs(table[j][i] - table[j - 1][i])),
							RequireFiniteRichardson(std::abs(table[j][i] - table[j - 1][i - 1])));
						if (candidateError <= bestError) {
							bestError = candidateError;
							answer = table[j][i];
						}
						if (j < i) factor = RequireFiniteRichardson(factor * shrinkSquared);
					}
					if (RequireFiniteRichardson(std::abs(table[i][i] - table[i - 1][i - 1])) / 2 >= bestError)
						break;
				}
				if (error) *error = bestError;
				return answer;
			}

			template<typename ComputeFn>
			DerivativeResult<Real> ExecuteDerivativeDetailed(const char* algorithm_name,
			                                              Real h,
			                                              const DerivativeConfig& config,
			                                              int evals_without_error,
			                                              int evals_with_error,
			                                              ComputeFn&& compute)
			{
				auto execute = [&]() {
					ValidateStepSize(h, algorithm_name);
					AlgorithmTimer timer;

					DerivativeResult<Real> result = MakeEvaluationSuccessResult<DerivativeResult<Real>>(
						algorithm_name,
						config.estimate_error ? evals_with_error : evals_without_error);
					result.step_used = h;

					Real error_estimate = NoEstimate;
					result.value = compute(config.estimate_error ? &error_estimate : nullptr);
					result.error = config.estimate_error ? error_estimate : NoEstimate;

					if (config.check_finite)
					{
						if (!IsFinite(result.value))
							throw NumericalMethodError(std::string(algorithm_name) + ": non-finite derivative result");
						if (config.estimate_error && !IsFinite(result.error))
							throw NumericalMethodError(std::string(algorithm_name) + ": non-finite derivative error estimate");
					}

					result.elapsed_time_ms = timer.elapsed_ms();
					return result;
				};

				if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
					return execute();

				try {
					return execute();
				}
				catch (const NumericInputError& ex) {
					auto result = MakeEvaluationFailureResult<DerivativeResult<Real>>(
						AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
					result.step_used = h;
					return result;
				}
				catch (const NumericalMethodError& ex) {
					auto result = MakeEvaluationFailureResult<DerivativeResult<Real>>(
						AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
					result.step_used = h;
					return result;
				}
				catch (const std::exception& ex) {
					auto result = MakeEvaluationFailureResult<DerivativeResult<Real>>(
						AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
					result.step_used = h;
					return result;
				}
			}
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		static DerivativeResult<Real> NDer1Detailed(const IRealFunction& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::One>;
			return Detail::ExecuteDerivativeDetailed("NDer1", h, config, int{Stencil::value_offsets.size()}, int{Stencil::error_offsets.size()},
				[&](Real* error) {
					auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::One>(
						[&](int offset) { return f(x + offset * h); }, [](Real value) { return std::abs(value); }, h, error != nullptr);
					if (error) *error = result.error;
					return result.value;
				});
		}
		static inline DerivativeResult<Real> NDer1Detailed(const IRealFunction& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return NDer1Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer1_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer1Detailed(Function&& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer1Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer1Detailed(Function&& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer1Detailed(adapter, x, config); });
		}
		static Real NDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer1Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static inline Real NDer1(const IRealFunction& f, Real x, Real* error)
		{
			return NDer1(f, x, ScaleStep(NDer1_h, x), error);
		}
		static inline Real NDer1(const IRealFunction& f, Real x)
		{
			return NDer1(f, x, ScaleStep(NDer1_h, x), nullptr);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer1(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer1(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer1(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer1(adapter, x, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer1(Function&& f, Real x)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer1(adapter, x); });
		}
		
		static inline Real NDer1Left(const IRealFunction& f, Real x, Real* error = nullptr)
		{ Real h = ScaleStep(NDer1_h, x); return NDer1(f, x - h, h, error); }
		static inline Real NDer1Right(const IRealFunction& f, Real x, Real* error = nullptr)
		{ Real h = ScaleStep(NDer1_h, x); return NDer1(f, x, h, error); }
		static inline Real NDer1Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{ return NDer1(f, x - h, h, error); }
		static inline Real NDer1Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{ return NDer1(f, x, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                ********/
		/********************************************************************************************************************/
		static DerivativeResult<Real> NDer2Detailed(const IRealFunction& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>;
			return Detail::ExecuteDerivativeDetailed("NDer2", h, config, int{Stencil::value_offsets.size()}, int{Stencil::error_offsets.size()},
				[&](Real* error) {
					auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>(
						[&](int offset) { return f(x + offset * h); }, [](Real value) { return std::abs(value); }, h, error != nullptr);
					if (error) *error = result.error;
					return result.value;
				});
		}
		static inline DerivativeResult<Real> NDer2Detailed(const IRealFunction& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return NDer2Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer2_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer2Detailed(Function&& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer2Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer2Detailed(Function&& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer2Detailed(adapter, x, config); });
		}
		static Real NDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer2Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static inline Real NDer2(const IRealFunction& f, Real x, Real* error)
		{
			return NDer2(f, x, ScaleStep(NDer2_h, x), error);
		}
		static inline Real NDer2(const IRealFunction& f, Real x)
		{
			return NDer2(f, x, ScaleStep(NDer2_h, x), nullptr);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer2(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer2(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer2(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer2(adapter, x, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer2(Function&& f, Real x)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer2(adapter, x); });
		}
		
		// Offset 3h keeps the full stencil (value ±h, error ±2h) strictly on one side of x,
		// matching the NDer4/6/8 family pattern (shift = error-stencil reach + 1).
		static inline Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer2_h, x); return NDer2(f, x - 3 * h, h, error); }
		static inline Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer2_h, x); return NDer2(f, x + 3 * h, h, error); }
		static inline Real NDer2Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x - 3 * h, h, error); }
		static inline Real NDer2Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x + 3 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		static DerivativeResult<Real> NDer4Detailed(const IRealFunction& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Four>;
			return Detail::ExecuteDerivativeDetailed("NDer4", h, config, int{Stencil::value_offsets.size()}, int{Stencil::error_offsets.size()},
				[&](Real* error) {
					auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::Four>(
						[&](int offset) { return f(x + offset * h); }, [](Real value) { return std::abs(value); }, h, error != nullptr);
					if (error) *error = result.error;
					return result.value;
				});
		}
		static inline DerivativeResult<Real> NDer4Detailed(const IRealFunction& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return NDer4Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer4_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer4Detailed(Function&& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer4Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer4Detailed(Function&& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer4Detailed(adapter, x, config); });
		}
		static Real NDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer4Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Real NDer4(const IRealFunction& f, Real x, Real* error)
		{
			return NDer4(f, x, ScaleStep(NDer4_h, x), error);
		}
		static Real NDer4(const IRealFunction& f, Real x)
		{
			return NDer4(f, x, ScaleStep(NDer4_h, x), nullptr);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer4(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer4(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer4(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer4(adapter, x, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer4(Function&& f, Real x)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer4(adapter, x); });
		}

		static inline Real NDer4Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer4_h, x); return NDer4(f, x - 4 * h, h, error); }
		static inline Real NDer4Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer4_h, x); return NDer4(f, x + 4 * h, h, error); }
		static inline Real NDer4Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x - 4 * h, h, error); }
		static inline Real NDer4Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x + 4 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		static DerivativeResult<Real> NDer6Detailed(const IRealFunction& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Six>;
			return Detail::ExecuteDerivativeDetailed("NDer6", h, config, int{Stencil::value_offsets.size()}, int{Stencil::error_offsets.size()},
				[&](Real* error) {
					auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::Six>(
						[&](int offset) { return f(x + offset * h); }, [](Real value) { return std::abs(value); }, h, error != nullptr);
					if (error) *error = result.error;
					return result.value;
				});
		}
		static inline DerivativeResult<Real> NDer6Detailed(const IRealFunction& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return NDer6Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer6_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer6Detailed(Function&& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer6Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer6Detailed(Function&& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer6Detailed(adapter, x, config); });
		}
		static Real NDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer6Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static inline Real NDer6(const IRealFunction& f, Real x, Real* error)
		{
			return NDer6(f, x, ScaleStep(NDer6_h, x), error);
		}
		static inline Real NDer6(const IRealFunction& f, Real x)
		{
			return NDer6(f, x, ScaleStep(NDer6_h, x), nullptr);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer6(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer6(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer6(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer6(adapter, x, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer6(Function&& f, Real x)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer6(adapter, x); });
		}

		static inline Real NDer6Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer6_h, x); return NDer6(f, x - 5 * h, h, error); }
		static inline Real NDer6Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer6_h, x); return NDer6(f, x + 5 * h, h, error); }
		static inline Real NDer6Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x - 5 * h, h, error); }
		static inline Real NDer6Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x + 5 * h, h, error); }

		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		static DerivativeResult<Real> NDer8Detailed(const IRealFunction& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			using Stencil = Detail::FirstDerivativeStencil<Detail::FirstDerivativeOrder::Eight>;
			return Detail::ExecuteDerivativeDetailed("NDer8", h, config, int{Stencil::value_offsets.size()}, int{Stencil::error_offsets.size()},
				[&](Real* error) {
					auto result = Detail::EvaluateFirstDerivativeStencil<Detail::FirstDerivativeOrder::Eight>(
						[&](int offset) { return f(x + offset * h); }, [](Real value) { return std::abs(value); }, h, error != nullptr);
					if (error) *error = result.error;
					return result.value;
				});
		}
		static inline DerivativeResult<Real> NDer8Detailed(const IRealFunction& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return NDer8Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer8_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer8Detailed(Function&& f, Real x, Real h,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer8Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NDer8Detailed(Function&& f, Real x,
		                                           const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer8Detailed(adapter, x, config); });
		}
		static Real NDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer8Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static inline Real NDer8(const IRealFunction& f, Real x, Real* error)
		{
			return NDer8(f, x, ScaleStep(NDer8_h, x), error);
		}
		static inline Real NDer8(const IRealFunction& f, Real x)
		{
			return NDer8(f, x, ScaleStep(NDer8_h, x), nullptr);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer8(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer8(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer8(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer8(adapter, x, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDer8(Function&& f, Real x)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDer8(adapter, x); });
		}

		static inline Real NDer8Left(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer8_h, x); return NDer8(f, x - 6 * h, h, error); }
		static inline Real NDer8Right(const IRealFunction& f, Real x, Real* error = nullptr) { Real h = ScaleStep(NDer8_h, x); return NDer8(f, x + 6 * h, h, error); }
		static inline Real NDer8Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x - 6 * h, h, error); }
		static inline Real NDer8Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x + 6 * h, h, error); }

		/********************************************************************************************************************/
		/********                                      SECOND DERIVATIVES                                            ********/
		/********  NOTE: Direct finite difference formulas (not derivatives of derivatives!)                        ********/
		/********        NSecDer2: 3 function evals (O(h²) accuracy)                                                ********/
		/********        NSecDer4: 5 function evals (O(h⁴) accuracy)                                                ********/
		/********************************************************************************************************************/
		
		// f''(x) ≈ [f(x-h) - 2f(x) + f(x+h)] / h²
		// Second-order accurate (O(h²)), 3 function evaluations
		static DerivativeResult<Real> NSecDer2Detailed(const IRealFunction& f, Real x, Real h,
		                                              const DerivativeConfig& config = {})
		{
			return Detail::ExecuteDerivativeDetailed("NSecDer2", h, config, 3, 5,
				[&](Real* error) {
					Real y0 = f(x);
					Real yh = f(x + h);
					Real ymh = f(x - h);

					Real h2 = h * h;
					Real result = (ymh - MML::Real(2.0) * y0 + yh) / h2;

					if (error)
					{
						Real y2h = f(x + 2 * h);
						Real ym2h = f(x - 2 * h);
						*error = Detail::EstimateSecondTruncation2(ym2h, ymh, y0, yh, y2h, h,
							[](Real value) { return std::abs(value); }) +
						         Constants::Eps * (std::abs(ymh) + 2 * std::abs(y0) + std::abs(yh)) / h2;
					}

					return result;
				});
		}
		static inline DerivativeResult<Real> NSecDer2Detailed(const IRealFunction& f, Real x,
		                                              const DerivativeConfig& config = {})
		{
			return NSecDer2Detailed(f, x, Detail::ResolveDerivativeStep(config, Precision::FourthRootEps, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NSecDer2Detailed(Function&& f, Real x, Real h,
		                                              const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer2Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NSecDer2Detailed(Function&& f, Real x,
		                                              const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer2Detailed(adapter, x, config); });
		}
		static Real NSecDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NSecDer2Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static inline Real NSecDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer2(f, x, ScaleStep(Precision::FourthRootEps, x), error);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDer2(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer2(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDer2(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer2(adapter, x, error); });
		}

		// f''(x) ≈ [-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)] / (12h²)
		// Fourth-order accurate (O(h⁴)), 5 function evaluations
		static DerivativeResult<Real> NSecDer4Detailed(const IRealFunction& f, Real x, Real h,
		                                              const DerivativeConfig& config = {})
		{
			return Detail::ExecuteDerivativeDetailed("NSecDer4", h, config, 5, 7,
				[&](Real* error) {
					Real y0 = f(x);
					Real yh = f(x + h);
					Real ymh = f(x - h);
					Real y2h = f(x + 2 * h);
					Real ym2h = f(x - 2 * h);

					Real h2 = h * h;
					Real result = (-ym2h + MML::Real(16.0) * ymh - MML::Real(30.0) * y0 + MML::Real(16.0) * yh - y2h) / (MML::Real(12.0) * h2);

					if (error)
					{
						Real y3h = f(x + 3 * h);
						Real ym3h = f(x - 3 * h);
						*error = Detail::EstimateSecondTruncation4(ym3h, ym2h, ymh, y0, yh, y2h, y3h, h,
							[](Real value) { return std::abs(value); }) +
						         Constants::Eps * (std::abs(ym2h) + 16 * std::abs(ymh) + 30 * std::abs(y0) +
						                           16 * std::abs(yh) + std::abs(y2h)) / (Real(12) * h2);
					}

					return result;
				});
		}
		static inline DerivativeResult<Real> NSecDer4Detailed(const IRealFunction& f, Real x,
		                                              const DerivativeConfig& config = {})
		{
			return NSecDer4Detailed(f, x, Detail::ResolveDerivativeStep(config, Precision::SixthRootEps, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NSecDer4Detailed(Function&& f, Real x, Real h,
		                                              const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer4Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NSecDer4Detailed(Function&& f, Real x,
		                                              const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer4Detailed(adapter, x, config); });
		}
		static Real NSecDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NSecDer4Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Real NSecDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer4(f, x, ScaleStep(Precision::SixthRootEps, x), error);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDer4(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer4(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDer4(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDer4(adapter, x, error); });
		}

		/********************************************************************************************************************/
		/********                                       THIRD DERIVATIVES                                            ********/
		/********  NOTE: Direct finite difference formulas (not derivatives of derivatives!)                        ********/
		/********        NThirdDer2: 4 function evals (O(h²) accuracy)                                              ********/
		/********        NThirdDer4: 6 function evals (O(h⁴) accuracy)                                              ********/
		/********************************************************************************************************************/
		
		// f'''(x) ≈ [-f(x-2h) + 2f(x-h) - 2f(x+h) + f(x+2h)] / (2h³)
		// Second-order accurate (O(h²)), 4 function evaluations
		static DerivativeResult<Real> NThirdDer2Detailed(const IRealFunction& f, Real x, Real h,
		                                                const DerivativeConfig& config = {})
		{
			return Detail::ExecuteDerivativeDetailed("NThirdDer2", h, config, 4, 6,
				[&](Real* error) {
					Real yh = f(x + h);
					Real ymh = f(x - h);
					Real y2h = f(x + 2 * h);
					Real ym2h = f(x - 2 * h);

					Real h3 = h * h * h;
					Real result = (-ym2h + MML::Real(2.0) * ymh - MML::Real(2.0) * yh + y2h) / (MML::Real(2.0) * h3);

					if (error)
					{
						Real y3h = f(x + 3 * h);
						Real ym3h = f(x - 3 * h);
						*error = Detail::EstimateThirdTruncation2(ym3h, ym2h, ymh, yh, y2h, y3h, h,
							[](Real value) { return std::abs(value); }) +
						         Constants::Eps * (std::abs(ym2h) + 2 * std::abs(ymh) + 2 * std::abs(yh) + std::abs(y2h)) / (Real(2) * h3);
					}

					return result;
				});
		}
		static inline DerivativeResult<Real> NThirdDer2Detailed(const IRealFunction& f, Real x,
		                                                const DerivativeConfig& config = {})
		{
			return NThirdDer2Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer4_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NThirdDer2Detailed(Function&& f, Real x, Real h,
		                                                const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer2Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NThirdDer2Detailed(Function&& f, Real x,
		                                                const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer2Detailed(adapter, x, config); });
		}
		static Real NThirdDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NThirdDer2Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Real NThirdDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer2(f, x, ScaleStep(NDer4_h, x), error);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDer2(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer2(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDer2(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer2(adapter, x, error); });
		}

		// f'''(x) ≈ [f(x-3h) - 8f(x-2h) + 13f(x-h) - 13f(x+h) + 8f(x+2h) - f(x+3h)] / (8h³)
		// Fourth-order accurate (O(h⁴)), 6 function evaluations
		static DerivativeResult<Real> NThirdDer4Detailed(const IRealFunction& f, Real x, Real h,
		                                                const DerivativeConfig& config = {})
		{
			return Detail::ExecuteDerivativeDetailed("NThirdDer4", h, config, 6, 8,
				[&](Real* error) {
					Real yh = f(x + h);
					Real ymh = f(x - h);
					Real y2h = f(x + 2 * h);
					Real ym2h = f(x - 2 * h);
					Real y3h = f(x + 3 * h);
					Real ym3h = f(x - 3 * h);

					Real h3 = h * h * h;
					Real result = (ym3h - MML::Real(8.0) * ym2h + MML::Real(13.0) * ymh - MML::Real(13.0) * yh + MML::Real(8.0) * y2h - y3h) / (MML::Real(8.0) * h3);

					if (error)
					{
						Real y4h = f(x + 4 * h);
						Real ym4h = f(x - 4 * h);
						*error = Detail::EstimateThirdTruncation4(ym4h, ym3h, ym2h, ymh, yh, y2h, y3h, y4h, h,
							[](Real value) { return std::abs(value); }) +
						         Constants::Eps * (std::abs(ym3h) + 8 * std::abs(ym2h) + 13 * std::abs(ymh) +
						                           13 * std::abs(yh) + 8 * std::abs(y2h) + std::abs(y3h)) / (Real(8) * h3);
					}

					return result;
				});
		}
		static inline DerivativeResult<Real> NThirdDer4Detailed(const IRealFunction& f, Real x,
		                                                const DerivativeConfig& config = {})
		{
			return NThirdDer4Detailed(f, x, Detail::ResolveDerivativeStep(config, NDer4_h, x), config);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NThirdDer4Detailed(Function&& f, Real x, Real h,
		                                                const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer4Detailed(adapter, x, h, config); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static DerivativeResult<Real> NThirdDer4Detailed(Function&& f, Real x,
		                                                const DerivativeConfig& config = {})
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer4Detailed(adapter, x, config); });
		}
		static Real NThirdDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NThirdDer4Detailed(f, x, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static inline Real NThirdDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer4(f, x, ScaleStep(NDer4_h, x), error);
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDer4(Function&& f, Real x, Real h, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer4(adapter, x, h, error); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDer4(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDer4(adapter, x, error); });
		}

		/********************************************************************************************************************/
		/********              RICHARDSON EXTRAPOLATION DERIVATIVES (dfridr-style adaptive)                          ********/
		/********  Uses tableau-based Richardson extrapolation to achieve high accuracy automatically               ********/
		/********  Based on Numerical Recipes "dfridr" algorithm with MML Richardson utility                        ********/
		/********************************************************************************************************************/
		
		/// @brief First derivative using Richardson extrapolation (dfridr algorithm).
		/// @details Builds a Richardson tableau using successively smaller step sizes,
		///          automatically finding the optimal step size. Uses central differences.
		///          Typically achieves ~10-12 digits of accuracy for smooth functions.
		///          All three Richardson derivatives reject invalid parameters with NumericInputError
		///          and collapsed/overflowing stencils or nonfinite arithmetic with NumericalMethodError.
		///          With max_iter=1 the error is unavailable (numeric_limits<Real>::max()).
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate derivative
		/// @param h Initial step size (will be refined automatically). Default: 0.1
		/// @param error Output parameter for error estimate (optional)
		/// @param max_iter Maximum tableau depth (default 10)
		/// @param con Step shrinking factor (default 1.4)
		/// @return Derivative estimate with error estimate if requested
		static Real NDerRichardson(const IRealFunction& f, Real x, Real h = 0.1, 
		                           Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			return Detail::EvaluateRichardsonDerivative<1>(f, x, h, error, max_iter, con);
		}

		/// @brief First derivative using Richardson extrapolation with default parameters.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Derivative estimate
		static inline Real NDerRichardson(const IRealFunction& f, Real x, Real* error)
		{
			return NDerRichardson(f, x, MML::Real(0.1), error, 10, MML::Real(1.4));
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDerRichardson(Function&& f, Real x, Real h = 0.1,
		                           Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDerRichardson(adapter, x, h, error, max_iter, con); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDerRichardson(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDerRichardson(adapter, x, error); });
		}

		/// @brief Second derivative using Richardson extrapolation.
		/// @details Uses central difference formula for second derivative with Richardson 
		///          extrapolation to achieve high accuracy automatically.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate second derivative
		/// @param h Initial step size (will be refined automatically). Default: 0.1
		/// @param error Output parameter for error estimate (optional)
		/// @param max_iter Maximum tableau depth (default 10)
		/// @param con Step shrinking factor (default 1.4)
		/// @return Second derivative estimate
		static Real NSecDerRichardson(const IRealFunction& f, Real x, Real h = 0.1,
		                              Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			return Detail::EvaluateRichardsonDerivative<2>(f, x, h, error, max_iter, con);
		}

		/// @brief Second derivative using Richardson extrapolation with default parameters.
		static inline Real NSecDerRichardson(const IRealFunction& f, Real x, Real* error)
		{
			return NSecDerRichardson(f, x, MML::Real(0.1), error, 10, MML::Real(1.4));
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDerRichardson(Function&& f, Real x, Real h = 0.1,
		                              Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDerRichardson(adapter, x, h, error, max_iter, con); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDerRichardson(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDerRichardson(adapter, x, error); });
		}

		/// @brief Third derivative using Richardson extrapolation.
		/// @details Uses central difference formula for third derivative with Richardson
		///          extrapolation. Requires more iterations for convergence.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate third derivative
		/// @param h Initial step size (default 0.2, larger than for lower derivatives)
		/// @param error Output parameter for error estimate (optional)
		/// @param max_iter Maximum tableau depth (default 10)
		/// @param con Step shrinking factor (default 1.4)
		/// @return Third derivative estimate
		static Real NThirdDerRichardson(const IRealFunction& f, Real x, Real h = 0.2,
		                                Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			return Detail::EvaluateRichardsonDerivative<3>(f, x, h, error, max_iter, con);
		}

		/// @brief Third derivative using Richardson extrapolation with default parameters.
		static inline Real NThirdDerRichardson(const IRealFunction& f, Real x, Real* error)
		{
			return NThirdDerRichardson(f, x, MML::Real(0.2), error, 10, MML::Real(1.4));
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDerRichardson(Function&& f, Real x, Real h = 0.2,
		                                Real* error = nullptr, int max_iter = 10, Real con = 1.4)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDerRichardson(adapter, x, h, error, max_iter, con); });
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDerRichardson(Function&& f, Real x, Real* error)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDerRichardson(adapter, x, error); });
		}

		/********************************************************************************************************************/
		/********                          ADAPTIVE STEP SIZE DERIVATIVES (Aliases)                                  ********/
		/********  These are aliases for Richardson extrapolation methods, which automatically                       ********/
		/********  find the optimal step size through tableau-based refinement.                                     ********/
		/********************************************************************************************************************/

		/// @brief Adaptive first derivative (alias for NDerRichardson).
		/// @details Automatically finds optimal step size using Richardson extrapolation.
		///          This solves the fundamental problem of choosing h: too small causes
		///          cancellation error, too large causes truncation error.
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Derivative estimate
		static inline Real NDerAdaptive(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NDerRichardson(f, x, MML::Real(0.1), error, 10, MML::Real(1.4));
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NDerAdaptive(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NDerAdaptive(adapter, x, error); });
		}

		/// @brief Adaptive second derivative (alias for NSecDerRichardson).
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate second derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Second derivative estimate
		static inline Real NSecDerAdaptive(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDerRichardson(f, x, MML::Real(0.1), error, 10, MML::Real(1.4));
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NSecDerAdaptive(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NSecDerAdaptive(adapter, x, error); });
		}

		/// @brief Adaptive third derivative (alias for NThirdDerRichardson).
		/// @param f Function to differentiate
		/// @param x Point at which to evaluate third derivative
		/// @param error Output parameter for error estimate (optional)
		/// @return Third derivative estimate
		static inline Real NThirdDerAdaptive(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDerRichardson(f, x, MML::Real(0.2), error, 10, MML::Real(1.4));
		}
		template<RealFunctionCallable Function>
			requires (!std::derived_from<std::remove_cvref_t<Function>, IRealFunction>)
		static Real NThirdDerAdaptive(Function&& f, Real x, Real* error = nullptr)
		{
			return Detail::WithRealFunctionAdapter(std::forward<Function>(f),
				[&](const IRealFunction& adapter) { return NThirdDerAdaptive(adapter, x, error); });
		}

		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		[[maybe_unused]] static inline Real(*Derive)(const IRealFunction& f, Real x) = Derivation::NDer4;
		[[maybe_unused]] static inline Real(*DeriveErr)(const IRealFunction& f,
																	 Real x, Real* error) = Derivation::NDer4;
		[[maybe_unused]] static inline Real(*DeriveSec)(const IRealFunction& f,
																	 Real x, Real* error) = Derivation::NSecDer4;
		[[maybe_unused]] static inline Real(*DeriveThird)(const IRealFunction& f,
																		 Real x, Real* error) = Derivation::NThirdDer2;
	}
}

#endif // MML_DERIVATION_REAL_FUNCTION_H
