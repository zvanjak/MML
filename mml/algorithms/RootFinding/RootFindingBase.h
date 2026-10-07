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
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Function.h>

namespace MML {
	namespace RootFinding {
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
			Real tolerance = 1e-10;

			/// Absolute tolerance for successive root estimates or bracket width.
			Real x_tolerance = 0.0;

			/// Absolute function-value tolerance. Zero inherits the legacy tolerance.
			Real f_tolerance = 0.0;

			/// Relative x tolerance, scaled by max(1, abs(x)).
			Real relative_tolerance = 0.0;
			
			/// Maximum number of iterations (default: 100)
			/// Set to 0 to use algorithm-specific default from Defaults namespace
			int max_iterations = 100;
			
			/// Initial step size for derivative-free methods (default: 0.01)
			/// Used by some algorithms for initial search or derivative approximation
			Real initial_step = 0.01;
			
			/// Enable verbose output for debugging (default: false)
			/// When true, prints iteration details to *verboseStream (nothing if null)
			bool verbose = false;

			/// Stream for verbose iteration output; null suppresses output (contract policy: no hardcoded std::cout)
			std::ostream* verboseStream = nullptr;

			/// Number of consecutive negligible-progress iterations before reporting Stalled.
			/// Zero disables explicit stagnation detection.
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
			Real x_error = 0.0;

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
			result.achieved_tolerance = 0.0;
			result.x_error = 0.0;
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
				} else if (_result.status == AlgorithmStatus::Success) {
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
