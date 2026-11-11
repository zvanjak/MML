///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingBase.h                                                   ///
///  Description: Configuration and result types for root-finding algorithms         ///
///               Configuration structs, result types, and forward declarations       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_BASE_H
#define MML_ROOTFINDING_BASE_H

#include <mml/MMLBase.h>
#include <algorithm>
#include <limits>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Function.h>

namespace MML {
	namespace RootFinding {
		namespace Detail {
			inline Real EvaluateFunction(const IRealFunction& function, Real x, int& evaluations) {
				if (evaluations == std::numeric_limits<int>::max())
					throw RootFindingError("Root finding: function evaluation counter overflow");
				++evaluations;
				return function(x);
			}

			class CountingRealFunction : public IRealFunction {
			public:
				explicit CountingRealFunction(const IRealFunction& function) : _function(function) {}

				Real operator()(Real x) const override {
					return EvaluateFunction(_function, x, _evaluations);
				}

				[[nodiscard]] int Evaluations() const noexcept { return _evaluations; }

			private:
				const IRealFunction& _function;
				mutable int _evaluations = 0;
			};

			inline bool HaveOppositeSigns(Real first, Real second) {
				return (first < 0 && second > 0) || (first > 0 && second < 0);
			}

			// Callers supply finite values with at least one nonzero endpoint.
			inline Real ScaledSecantFraction(Real numerator, Real other) {
				const Real scale = std::max(std::abs(numerator), std::abs(other));
				const Real normalizedNumerator = numerator / scale;
				const Real difference = normalizedNumerator - other / scale;
				if (difference == 0) return NoEstimate; // tol-ok: exact zero denominator; caller reports the non-finite step.
				return normalizedNumerator / difference;
			}

			// sqrt(abs(first * second)), without forming an overflowing/underflowing product.
			inline Real GeometricMeanMagnitude(Real first, Real second) {
				int firstExponent, secondExponent;
				Real product = std::frexp(std::abs(first), &firstExponent)
					* std::frexp(std::abs(second), &secondExponent);
				int exponent = firstExponent + secondExponent;
				if (exponent % 2 != 0) {
					product *= 2;
					--exponent;
				}
				return std::ldexp(std::sqrt(product), exponent / 2);
			}
		}

		/*********************************************************************/
		/*****           Configuration and Result Types                  *****/
		/*********************************************************************/

		/// Configuration parameters for root-finding algorithms.
		/// 
		/// Provides user control over convergence criteria, iteration limits,
		/// and diagnostic output. All fields have sensible defaults.
		/// 
		/// @example
		/// RootFindingConfig config;
		/// config.tolerance = 1e-12;        // Higher precision
		/// config.max_iterations = 200;     // Allow more iterations
		/// auto result = FindRootBrent(f, 0, 1, config);
		struct RootFindingConfig {
			/// Legacy absolute x tolerance. Used when x_tolerance is zero.
			Real tolerance = Precision::EqualityTolerance;

			/// Absolute tolerance for successive root estimates or bracket width.
			Real x_tolerance = 0.0;

			/// Absolute function-value tolerance. Zero inherits the legacy tolerance.
			Real f_tolerance = 0.0;

			/// Relative x tolerance, scaled by max(1, abs(x)).
			Real relative_tolerance = 0.0;
			
			/// Maximum number of iterations (default: 100)
			/// Set to 0 to use algorithm-specific default from Defaults namespace
			int max_iterations = 100;
			
			/// Reserved compatibility field; scalar refinement solvers do not use it.
			Real initial_step = MML::Real(0.01);
			
			/// Enable verbose output for debugging (default: false)
			/// When true, prints iteration details to *verboseStream (nothing if null)
			bool verbose = false;

			/// Stream for verbose iteration output; null suppresses output (contract policy: no hardcoded std::cout)
			std::ostream* verboseStream = nullptr;

			/// Reserved compatibility field; does not enable configurable stagnation detection.
			/// Must remain nonnegative for configuration validation.
			int max_stagnant_iterations = 3;
		};

		/// Result of a root-finding operation.
		/// 
		/// Contains the root value along with diagnostic information about
		/// the convergence process. Always check `converged` before using `root`.
		/// 
		/// @example
		/// auto result = FindRootBrent(f, 0, 1, config);
		/// if (result.converged) {
		///     std::cout << "Root: " << result.root << " found in " 
		///               << result.iterations_used << " iterations\n";
		/// } else {
		///     std::cerr << "Failed to converge after " << result.iterations_used << " iterations\n";
		/// }
		struct RootFindingResult : public IterativeResultBase {
			/// The computed root value
			Real root = 0.0;
			
			/// Function value at root: f(root), should be near zero if converged
			Real function_value = 0.0;
			
			/// Root estimate change or final bracket width, depending on the method.
			Real x_error = NoEstimate;

			[[nodiscard]] bool IsSuccess() const noexcept {
				return converged && status == AlgorithmStatus::Success;
			}

			explicit operator bool() const noexcept { return IsSuccess(); }
		};

		inline bool IsValidConfig(const RootFindingConfig& config) {
			return std::isfinite(config.tolerance) && config.tolerance > 0.0
				&& std::isfinite(config.x_tolerance) && config.x_tolerance >= 0.0
				&& std::isfinite(config.f_tolerance) && config.f_tolerance >= 0.0
				&& std::isfinite(config.relative_tolerance) && config.relative_tolerance >= 0.0
				&& config.max_iterations >= 0 && config.max_stagnant_iterations >= 0;
		}

		inline Real EffectiveXTolerance(const RootFindingConfig& config, Real x) {
			const Real absoluteTolerance = config.x_tolerance > 0.0
				? config.x_tolerance : config.tolerance;
			return absoluteTolerance
				+ config.relative_tolerance * std::max(Real(1.0), std::abs(x));
		}

		inline Real EffectiveFTolerance(const RootFindingConfig& config) {
			return config.f_tolerance > 0.0 ? config.f_tolerance : config.tolerance;
		}

		inline bool MeetsConvergence(const RootFindingConfig& config, Real x,
			Real xError, Real functionValue) {
			return xError <= EffectiveXTolerance(config, x)
				|| std::abs(functionValue) <= EffectiveFTolerance(config);
		}

		inline bool AreValidEndpoints(Real x1, Real x2) {
			return std::isfinite(x1) && std::isfinite(x2) && x1 != x2;
		}

		inline bool AcceptEndpointRoot(RootFindingResult& result,
			const RootFindingConfig& config, Real x1, Real f1, Real x2, Real f2) {
			if (std::abs(f1) <= EffectiveFTolerance(config)) {
				result.root = x1;
				result.function_value = f1;
			} else if (std::abs(f2) <= EffectiveFTolerance(config)) {
				result.root = x2;
				result.function_value = f2;
			} else {
				return false;
			}
			result.converged = true;
			// The endpoint satisfied the f-tolerance: that residual is the achieved tolerance.
			// No bracket was shrunk, so there is no x-error estimate (RESULT_POLICY §3).
			result.achieved_tolerance = std::abs(result.function_value);
			result.x_error = NoEstimate;
			return true;
		}

		class RootFindingResultFinalizer {
		public:
			RootFindingResultFinalizer(RootFindingResult& result, const char* algorithmName)
				: _result(result) {
				_result.algorithm_name = algorithmName;
			}

			~RootFindingResultFinalizer() {
				_result.elapsed_time_ms = _timer.elapsed_ms();
				if (_result.converged) {
					_result.status = AlgorithmStatus::Success;
					_result.error_message.clear();
				} else if (_result.status == AlgorithmStatus::Success || _result.status == AlgorithmStatus::NotComputed) {
					_result.status = AlgorithmStatus::AlgorithmSpecificFailure;
				}
			}

		private:
			RootFindingResult& _result;
			AlgorithmTimer _timer;
		};

	} // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_BASE_H
