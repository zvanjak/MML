///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationComplex.h                                                 ///
///  Description: Numerical derivatives of complex-valued functions f:C->C            ///
///               Forward, central, and higher-order finite differences               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_COMPLEX_H
#define MML_DERIVATION_COMPLEX_H

#include "MMLBase.h"

#include "core/AlgorithmTypes.h"
#include "core/NumericValidation.h"
#include "interfaces/IComplexFunction.h"

#include "DerivationBase.h"

namespace MML
{
	namespace Derivation
	{
		namespace Detail
		{
			/// @brief Scale step size relative to |z| for complex arguments
			static inline Real ScaleStepComplex(Real h, Complex z) {
				return h * std::max(Real(1), std::abs(z));
			}

			static inline Real ResolveComplexDerivativeStep(const DerivativeConfig& config, Real default_step, Complex z)
			{
				return config.step != REAL(0.0) ? config.step : ScaleStepComplex(default_step, z);
			}

			/// @brief Check if a complex value has finite real and imaginary parts
			static inline bool IsFiniteComplex(Complex z) {
				return std::isfinite(z.real()) && std::isfinite(z.imag());
			}

			template<typename ComputeFn>
			DerivativeResult<Complex, Real> ExecuteComplexDerivativeDetailed(
				const char* algorithm_name,
				Real h,
				const DerivativeConfig& config,
				int evals_without_error,
				int evals_with_error,
				ComputeFn&& compute)
			{
				auto execute = [&]() {
					ValidateStepSize(h, algorithm_name);
					AlgorithmTimer timer;

					DerivativeResult<Complex, Real> result = MakeEvaluationSuccessResult<DerivativeResult<Complex, Real>>(
						algorithm_name,
						config.estimate_error ? evals_with_error : evals_without_error);
					result.step_used = h;

					Real error_estimate = REAL(0.0);
					result.value = compute(config.estimate_error ? &error_estimate : nullptr);
					result.error = config.estimate_error ? error_estimate : REAL(0.0);

					if (config.check_finite)
					{
						if (!IsFiniteComplex(result.value))
							throw NumericalMethodError(std::string(algorithm_name) + ": non-finite complex derivative result");
						if (config.estimate_error && !std::isfinite(result.error))
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
					auto result = MakeEvaluationFailureResult<DerivativeResult<Complex, Real>>(
						AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
					result.step_used = h;
					return result;
				}
				catch (const NumericalMethodError& ex) {
					auto result = MakeEvaluationFailureResult<DerivativeResult<Complex, Real>>(
						AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
					result.step_used = h;
					return result;
				}
				catch (const std::exception& ex) {
					auto result = MakeEvaluationFailureResult<DerivativeResult<Complex, Real>>(
						AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
					result.step_used = h;
					return result;
				}
			}
		}

		/********************************************************************************************************************/
		/********                    Complex derivatives of FIRST order (forward difference)                          ********/
		/********************************************************************************************************************/

		/// @brief First-order forward difference for complex functions f:C→C
		/// @details f'(z) ≈ (f(z+h) - f(z)) / h with real step h
		/// @note Error is O(h); use NDer2Complex for better accuracy
		static DerivativeResult<Complex, Real> NDer1ComplexDetailed(const IComplexFunction& f, Complex z, Real h,
		                                                           const DerivativeConfig& config = {})
		{
			return Detail::ExecuteComplexDerivativeDetailed("NDer1Complex", h, config, 2, 3,
				[&](Real* error) {
					Complex yh = f(z + h);
					Complex y0 = f(z);
					Complex diff = yh - y0;
					if (error)
					{
						Complex ym = f(z - h);
						Real ypph = std::abs(yh - REAL(2.0) * y0 + ym) / h;
						*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Eps / h;
					}
					return diff / h;
				});
		}
		static DerivativeResult<Complex, Real> NDer1ComplexDetailed(const IComplexFunction& f, Complex z,
		                                                           const DerivativeConfig& config = {})
		{
			return NDer1ComplexDetailed(f, z, Detail::ResolveComplexDerivativeStep(config, NDer1_h, z), config);
		}
		static Complex NDer1Complex(const IComplexFunction& f, Complex z, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer1ComplexDetailed(f, z, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Complex NDer1Complex(const IComplexFunction& f, Complex z, Real* error)
		{
			return NDer1Complex(f, z, Detail::ScaleStepComplex(NDer1_h, z), error);
		}
		static Complex NDer1Complex(const IComplexFunction& f, Complex z)
		{
			return NDer1Complex(f, z, Detail::ScaleStepComplex(NDer1_h, z), nullptr);
		}

		/********************************************************************************************************************/
		/********                    Complex derivatives of SECOND order (central difference)                         ********/
		/********************************************************************************************************************/

		/// @brief Second-order central difference for complex functions f:C→C
		/// @details f'(z) ≈ (f(z+h) - f(z-h)) / (2h) with real step h
		/// @note Error is O(h²); recommended default for complex differentiation
		static DerivativeResult<Complex, Real> NDer2ComplexDetailed(const IComplexFunction& f, Complex z, Real h,
		                                                           const DerivativeConfig& config = {})
		{
			return Detail::ExecuteComplexDerivativeDetailed("NDer2Complex", h, config, 2, 4,
				[&](Real* error) {
					Complex yh = f(z + h);
					Complex ymh = f(z - h);
					Complex diff = yh - ymh;

					if (error)
					{
						Complex y2h = f(z + REAL(2.0) * h);
						Complex ym2h = f(z - REAL(2.0) * h);
						*error = Constants::Eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) +
						         std::abs((y2h - ym2h) / REAL(2.0) - diff) / (6 * h);
					}

					return diff / (REAL(2.0) * h);
				});
		}
		static DerivativeResult<Complex, Real> NDer2ComplexDetailed(const IComplexFunction& f, Complex z,
		                                                           const DerivativeConfig& config = {})
		{
			return NDer2ComplexDetailed(f, z, Detail::ResolveComplexDerivativeStep(config, NDer2_h, z), config);
		}
		static Complex NDer2Complex(const IComplexFunction& f, Complex z, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer2ComplexDetailed(f, z, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Complex NDer2Complex(const IComplexFunction& f, Complex z, Real* error)
		{
			return NDer2Complex(f, z, Detail::ScaleStepComplex(NDer2_h, z), error);
		}
		static Complex NDer2Complex(const IComplexFunction& f, Complex z)
		{
			return NDer2Complex(f, z, Detail::ScaleStepComplex(NDer2_h, z), nullptr);
		}

		/********************************************************************************************************************/
		/********                    Complex derivatives of FOURTH order (5-point stencil)                            ********/
		/********************************************************************************************************************/

		/// @brief Fourth-order 5-point stencil for complex functions f:C→C
		/// @details f'(z) ≈ (-f(z+2h) + 8f(z+h) - 8f(z-h) + f(z-2h)) / (12h)
		/// @note Error is O(h⁴); best accuracy/cost tradeoff for smooth functions
		static DerivativeResult<Complex, Real> NDer4ComplexDetailed(const IComplexFunction& f, Complex z, Real h,
		                                                           const DerivativeConfig& config = {})
		{
			return Detail::ExecuteComplexDerivativeDetailed("NDer4Complex", h, config, 4, 6,
				[&](Real* error) {
					Complex yh = f(z + h);
					Complex ymh = f(z - h);
					Complex y2h = f(z + REAL(2.0) * h);
					Complex ym2h = f(z - REAL(2.0) * h);

					Complex y2 = ym2h - y2h;
					Complex y1 = yh - ymh;

					if (error)
					{
						Complex y3h = f(z + REAL(3.0) * h);
						Complex ym3h = f(z - REAL(3.0) * h);

						*error = std::abs((y3h - ym3h) / REAL(2.0) + REAL(2.0) * (ym2h - y2h) +
						                  REAL(5.0) * (yh - ymh) / REAL(2.0)) / (30 * h);
						*error += Constants::Eps * (std::abs(y2h) + std::abs(ym2h) +
						                            8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
					}

					return (y2 + REAL(8.0) * y1) / (REAL(12.0) * h);
				});
		}
		static DerivativeResult<Complex, Real> NDer4ComplexDetailed(const IComplexFunction& f, Complex z,
		                                                           const DerivativeConfig& config = {})
		{
			return NDer4ComplexDetailed(f, z, Detail::ResolveComplexDerivativeStep(config, NDer4_h, z), config);
		}
		static Complex NDer4Complex(const IComplexFunction& f, Complex z, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer4ComplexDetailed(f, z, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Complex NDer4Complex(const IComplexFunction& f, Complex z, Real* error)
		{
			return NDer4Complex(f, z, Detail::ScaleStepComplex(NDer4_h, z), error);
		}
		static Complex NDer4Complex(const IComplexFunction& f, Complex z)
		{
			return NDer4Complex(f, z, Detail::ScaleStepComplex(NDer4_h, z), nullptr);
		}

		/********************************************************************************************************************/
		/********                    Complex derivatives of SIXTH order (7-point stencil)                             ********/
		/********************************************************************************************************************/

		/// @brief Sixth-order 7-point stencil for complex functions f:C→C
		/// @details f'(z) ≈ (f(z+3h) - 9f(z+2h) + 45f(z+h) - 45f(z-h) + 9f(z-2h) - f(z-3h)) / (60h)
		/// @note Error is O(h⁶); highest accuracy, uses 6 (or 8 with error) evaluations
		static DerivativeResult<Complex, Real> NDer6ComplexDetailed(const IComplexFunction& f, Complex z, Real h,
		                                                           const DerivativeConfig& config = {})
		{
			return Detail::ExecuteComplexDerivativeDetailed("NDer6Complex", h, config, 6, 8,
				[&](Real* error) {
					Complex yh = f(z + h);
					Complex ymh = f(z - h);
					Complex y1 = yh - ymh;
					Complex y2 = f(z - REAL(2.0) * h) - f(z + REAL(2.0) * h);
					Complex y3 = f(z + REAL(3.0) * h) - f(z - REAL(3.0) * h);

					if (error)
					{
						Complex y7 = (f(z + REAL(4.0) * h) - f(z - REAL(4.0) * h)
						              - REAL(6.0) * y3 - REAL(14.0) * y1 - REAL(14.0) * y2) / REAL(2.0);
						*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Eps / h;
					}

					return (y3 + REAL(9.0) * y2 + REAL(45.0) * y1) / (REAL(60.0) * h);
				});
		}
		static DerivativeResult<Complex, Real> NDer6ComplexDetailed(const IComplexFunction& f, Complex z,
		                                                           const DerivativeConfig& config = {})
		{
			return NDer6ComplexDetailed(f, z, Detail::ResolveComplexDerivativeStep(config, NDer6_h, z), config);
		}
		static Complex NDer6Complex(const IComplexFunction& f, Complex z, Real h, Real* error = nullptr)
		{
			DerivativeConfig config;
			config.estimate_error = (error != nullptr);
			auto result = NDer6ComplexDetailed(f, z, h, config);
			if (error)
				*error = result.error;
			return result.value;
		}
		static Complex NDer6Complex(const IComplexFunction& f, Complex z, Real* error)
		{
			return NDer6Complex(f, z, Detail::ScaleStepComplex(NDer6_h, z), error);
		}
		static Complex NDer6Complex(const IComplexFunction& f, Complex z)
		{
			return NDer6Complex(f, z, Detail::ScaleStepComplex(NDer6_h, z), nullptr);
		}

	} // namespace Derivation
} // namespace MML

#endif // MML_DERIVATION_COMPLEX_H
