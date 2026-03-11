///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODEAdaptiveIntegrator.h                                             ///
///  Description: Production-ready adaptive ODE integration with dense output        ///
///               Features: FSAL optimization, Hermite interpolation, diagnostics    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                   ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_ADAPTIVE_INTEGRATOR_H
#define MML_ODE_ADAPTIVE_INTEGRATOR_H

#include "mml/MMLBase.h"
#include "mml/core/AlgorithmTypes.h"
#include "mml/interfaces/IODESystem.h"
#include "mml/base/ODESystemSolution.h"
#include "mml/algorithms/ODESolvers/ODESystemSteppers.h"

namespace MML {

	/******************************************************************************/
	/*****            ODE Integrator Configuration                            *****/
	/******************************************************************************/

	/// Configuration for ODE integrators.
	///
	/// Provides user control over step size bounds, tolerances, and integration limits.
	/// All fields have sensible defaults suitable for most problems.
	///
	/// @example
	/// ODEIntegratorConfig config;
	/// config.tolerance = 1e-12;        // Higher precision
	/// config.max_steps = 100000;       // Allow more steps for stiff problems
	/// config.initial_step_size = 0.0;  // Auto-estimate (recommended)
	/// auto result = integrator.integrate(x0, t0, tEnd, config);
	struct ODEIntegratorConfig {
		/// Initial step size (default: 0.0 = auto-estimate)
		/// Set to 0 to use Hairer's automatic step size estimation
		Real initial_step_size = 0.0;

		/// Minimum allowed step size (default: 1e-15)
		/// Integration fails if step size falls below this
		Real min_step_size = 1e-15;

		/// Maximum allowed step size (default: infinity)
		/// Useful to limit step growth in smooth regions
		Real max_step_size = std::numeric_limits<Real>::infinity();

		/// Error tolerance for step acceptance (default: 1e-10)
		/// Used for both absolute and relative error control
		Real tolerance = 1e-10;

		/// Maximum number of accepted steps (default: 100000)
		/// Integration fails if this is exceeded
		int max_steps = 100000;

		/// Output interval for dense output (default: 0.0 = output at integration steps)
		/// If > 0, solutions are interpolated at fixed intervals
		Real output_interval = 0.0;

		/// Enable verbose output for debugging (default: false)
		bool verbose = false;

		/// Factory method for high-precision configuration
		static ODEIntegratorConfig HighPrecision() {
			ODEIntegratorConfig config;
			config.tolerance = 1e-14;
			config.max_steps = 500000;
			return config;
		}

		/// Factory method for fast (lower precision) configuration
		static ODEIntegratorConfig Fast() {
			ODEIntegratorConfig config;
			config.tolerance = 1e-6;
			config.max_steps = 10000;
			return config;
		}

		/// Factory method for stiff problems
		static ODEIntegratorConfig Stiff() {
			ODEIntegratorConfig config;
			config.tolerance = 1e-8;
			config.max_steps = 1000000;
			config.min_step_size = 1e-20;
			return config;
		}
	};

	//===================================================================================
	//                           SolutionStatistics
	//===================================================================================
	/// @brief Diagnostics for adaptive integration with timing information
	struct SolutionStatistics {
		// === Step Counts ===
		int acceptedSteps = 0;	///< Number of accepted steps
		int rejectedSteps = 0;	///< Number of rejected steps
		int totalFuncEvals = 0; ///< Total function evaluations

		// === Step Size Statistics ===
		Real minStepSize = 0;		///< Smallest step size used
		Real maxStepSize = 0;		///< Largest step size used

		// === Diagnostics ===
		std::string algorithm_name = "DormandPrince5";  ///< Stepper algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;  ///< Error description (empty on success)
		double elapsed_time_ms = 0.0;  ///< Total integration time

		void reset() {
			acceptedSteps = rejectedSteps = totalFuncEvals = 0;
			minStepSize = std::numeric_limits<Real>::max();
			maxStepSize = 0;
			status = AlgorithmStatus::Success;
			error_message.clear();
			elapsed_time_ms = 0.0;
		}

		void recordStep(const StepResult& result) {
			if (result.accepted) {
				acceptedSteps++;
				minStepSize = std::min(minStepSize, result.hDone);
				maxStepSize = std::max(maxStepSize, result.hDone);
			} else {
				rejectedSteps++;
			}
			totalFuncEvals += result.funcEvals;
		}

		Real acceptanceRate() const {
			int total = acceptedSteps + rejectedSteps;
			return total > 0 ? static_cast<Real>(acceptedSteps) / total : 1.0;
		}

		int totalSteps() const {
			return acceptedSteps + rejectedSteps;
		}
	};

	//===================================================================================
	//                           ODEAdaptiveIntegrator
	//===================================================================================
	/// @brief Adaptive ODE integrator with dense output support
	///
	/// This integrator uses an adaptive stepper to efficiently solve ODEs while
	/// providing solution output at user-specified fixed intervals (dense output).
	///
	/// Key features:
	/// - Adaptive step size control for efficiency and accuracy
	/// - Dense output: solution at arbitrary times without extra function evals
	/// - FSAL optimization when supported by stepper
	/// - Detailed statistics for performance analysis
	template<typename Stepper = DormandPrince5_Stepper>
	class ODEAdaptiveIntegrator {
	private:
		const IODESystem& _sys;
		Stepper _stepper;
		SolutionStatistics _stats;

	public:
		/// @brief Construct integrator for given ODE system
		explicit ODEAdaptiveIntegrator(const IODESystem& sys)
				: _sys(sys)
				, _stepper(sys) {}

		/// @brief Get integration statistics from last solve
		const SolutionStatistics& getStatistics() const { return _stats; }

		/// @brief Estimate initial step size
		/// Uses Hairer's formula based on local error estimation
		Real estimateInitialStep(Real t0, const Vector<Real>& x0, Real tEnd, Real eps) {
			int n = _sys.getDim();
			Vector<Real> f0(n), f1(n), x1(n);

			_sys.derivs(t0, x0, f0);

			// Estimate scale of problem
			Real d0 = 0, d1 = 0;
			for (int i = 0; i < n; i++) {
				Real scale = std::abs(x0[i]) + 1e-10;
				d0 = std::max(d0, std::abs(x0[i]) / scale);
				d1 = std::max(d1, std::abs(f0[i]) / scale);
			}

			// Initial guess
			Real h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
			h0 = std::min(h0, std::abs(tEnd - t0));

			// Take one explicit Euler step
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] + h0 * f0[i];
			_sys.derivs(t0 + h0, x1, f1);

			// Estimate second derivative
			Real d2 = 0;
			for (int i = 0; i < n; i++) {
				Real scale = std::abs(x0[i]) + 1e-10;
				d2 = std::max(d2, std::abs(f1[i] - f0[i]) / scale / h0);
			}

			// Optimal step for 5th order method
			Real h1;
			if (std::max(d1, d2) <= 1e-15) {
				h1 = std::max(1e-6, h0 * 1e-3);
			} else {
				h1 = std::pow(0.01 / std::max(d1, d2), 1.0 / 5.0);
			}

			return std::min({100 * h0, h1, std::abs(tEnd - t0)});
		}

		/// @brief Integrate ODE system with dense output at fixed intervals
		/// @param x0 Initial conditions
		/// @param t0 Initial time
		/// @param tEnd Final time
		/// @param outputInterval Spacing between output points (dense output)
		/// @param eps Error tolerance (default 1e-10)
		/// @param h0 Initial step size (0 = auto-estimate)
		/// @return Solution with points at regular intervals
		ODESystemSolution integrate(const Vector<Real>& x0, Real t0, Real tEnd, Real outputInterval, Real eps = 1e-10, Real h0 = 0) {
			int n = _sys.getDim();
			int numOutputPoints = static_cast<int>(std::ceil((tEnd - t0) / outputInterval)) + 1;

			ODESystemSolution sol(t0, tEnd, n, numOutputPoints - 1);
			_stats.reset();
			_stepper.resetFSAL();

			// Initialize
			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;

			_sys.derivs(t, x, dxdt);

			// Initial step size
			Real h = (h0 > 0) ? h0 : estimateInitialStep(t0, x0, tEnd, eps);
			h = std::min(h, tEnd - t0);

			// Store initial point
			Real tNextOutput = t0;
			int outputIdx = 0;
			sol.fillValues(outputIdx++, t0, x0);
			tNextOutput += outputInterval;

			// Main integration loop
			while (t < tEnd - Constants::Eps) {
				// Don't step past end
				if (t + h > tEnd) {
					h = tEnd - t;
				}

				StepResult result = _stepper.doStep(t, x, dxdt, h, eps);
				_stats.recordStep(result);

				if (result.accepted) {
					Real tNew = t + result.hDone;

					// Dense output: interpolate at all output points within this step
					while (tNextOutput <= tNew + Constants::Eps && outputIdx < numOutputPoints) {
						if (tNextOutput <= tNew) {
							Vector<Real> xInterp = _stepper.interpolate(tNextOutput);
							sol.fillValues(outputIdx++, tNextOutput, xInterp);
						}
						tNextOutput += outputInterval;
					}

					t = tNew;
					h = result.hNext;

					// FSAL: derivative already updated in x and dxdt by doStep
				} else {
					// Step rejected, try again with smaller step
					h = result.hNext;
				}

				// Safety check
				if (h < Constants::Eps) {
					throw ODESolverError("Step size too small in ODEAdaptiveIntegrator");
				}
			}

			// Ensure final point is captured
			if (outputIdx < numOutputPoints) {
				sol.fillValues(outputIdx++, tEnd, x);
			}

			return sol;
		}

		/// @brief Integrate and return solution at specified times
		/// @param x0 Initial conditions
		/// @param times Vector of output times (must be sorted ascending)
		/// @param eps Error tolerance
		/// @return Solution at specified times
		ODESystemSolution integrateAt(const Vector<Real>& x0, const Vector<Real>& times, Real eps = 1e-10, Real h0 = 0) {
			if (times.size() < 2) {
				throw ODESolverError("Need at least 2 output times");
			}

			int n = _sys.getDim();
			int numPoints = static_cast<int>(times.size());
			Real t0 = times[0];
			Real tEnd = times[numPoints - 1];

			ODESystemSolution sol(t0, tEnd, n, numPoints - 1);
			_stats.reset();
			_stepper.resetFSAL();

			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;

			_sys.derivs(t, x, dxdt);

			Real h = (h0 > 0) ? h0 : estimateInitialStep(t0, x0, tEnd, eps);

			// Store initial point
			int outputIdx = 0;
			sol.fillValues(outputIdx++, t0, x0);

			// Integration loop
			while (outputIdx < numPoints) {
				Real tTarget = times[outputIdx];

				if (t + h > tTarget) {
					h = tTarget - t;
				}

				StepResult result = _stepper.doStep(t, x, dxdt, h, eps);
				_stats.recordStep(result);

				if (result.accepted) {
					Real tNew = t + result.hDone;

					// Output at any times within this step
					while (outputIdx < numPoints && times[outputIdx] <= tNew + Constants::Eps) {
						Real tOut = times[outputIdx];
						if (std::abs(tOut - tNew) < Constants::Eps) {
							sol.fillValues(outputIdx++, tNew, x);
						} else {
							Vector<Real> xInterp = _stepper.interpolate(tOut);
							sol.fillValues(outputIdx++, tOut, xInterp);
						}
					}

					t = tNew;
					h = result.hNext;
				} else {
					h = result.hNext;
				}

				if (h < Constants::Eps && t < tEnd - Constants::Eps) {
					throw ODESolverError("Step size too small in ODEAdaptiveIntegrator");
				}
			}

			return sol;
		}

		/// @brief Integrate ODE system using configuration object
		/// 
		/// @param x0 Initial conditions
		/// @param t0 Initial time
		/// @param tEnd Final time
		/// @param config Integration configuration
		/// @return Solution with timing and diagnostic information in statistics
		ODESystemSolution integrate(const Vector<Real>& x0, Real t0, Real tEnd, const ODEIntegratorConfig& config) {
			AlgorithmTimer timer;  // Starts automatically

			Real outputInterval = (config.output_interval > 0) 
				? config.output_interval 
				: (tEnd - t0) / 100.0;  // Default: 100 output points

			ODESystemSolution sol = integrate(x0, t0, tEnd, outputInterval, 
				config.tolerance, config.initial_step_size);

			// Populate diagnostic fields in statistics
			_stats.elapsed_time_ms = timer.elapsed_ms();
			
			// Check for max steps exceeded (if we tracked it)
			if (_stats.acceptedSteps >= config.max_steps) {
				_stats.status = AlgorithmStatus::MaxIterationsExceeded;
				_stats.error_message = "Maximum steps (" + std::to_string(config.max_steps) + ") exceeded";
			}

			return sol;
		}

		//===================================================================================
		//                           Event Detection Integration
		//===================================================================================

		/**
 * @brief Result of integration with event detection.
 * 
 * Contains both the solution trajectory and information about detected events.
 */
		struct EventResult {
			ODESystemSolution solution;		 ///< Solution trajectory
			std::vector<EventInfo> events; ///< All detected events in chronological order
			bool terminatedByEvent;				 ///< True if integration ended due to EventAction::Stop
			Real finalTime;								 ///< Actual final time (may differ from tEnd if stopped)
			Vector<Real> finalState;			 ///< Final state vector

			EventResult(Real t0, Real tEnd, int dim, int numPoints)
					: solution(t0, tEnd, dim, numPoints)
					, terminatedByEvent(false)
					, finalTime(0) {}
		};

		/**
   * @brief Integrate ODE with event detection using zero-crossing detection.
   * 
   * Integrates the system while monitoring event functions g_i(t, x). When any
   * event function crosses zero in the specified direction, the crossing is
   * located precisely using bisection on dense output interpolation.
   * 
   * @param sys ODE system with event functions (must implement IODESystemWithEvents)
   * @param x0 Initial state vector
   * @param t0 Initial time
   * @param tEnd Final time
   * @param outputInterval Spacing between output points
   * @param eps Error tolerance
   * @param eventTol Tolerance for locating event time (default 1e-12)
   * @param h0 Initial step size (0 = auto-estimate)
   * @return EventResult with solution and event information
   */
		EventResult integrateWithEvents(const IODESystemWithEvents& sys, const Vector<Real>& x0, Real t0, Real tEnd, Real outputInterval,
																		Real eps = 1e-10, Real eventTol = 1e-12, Real h0 = 0) {
			int n = sys.getDim();
			int numEvents = sys.getNumEvents();
			int numOutputPoints = static_cast<int>(std::ceil((tEnd - t0) / outputInterval)) + 1;

			EventResult result(t0, tEnd, n, numOutputPoints - 1);
			_stats.reset();
			_stepper.resetFSAL();

			// Current state
			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;
			sys.derivs(t, x, dxdt);

			// Event function values at current time
			Vector<Real> gPrev(numEvents);
			sys.eventFunctions(t, x, gPrev);

			// Initial step size
			Real h = (h0 > 0) ? h0 : estimateInitialStep(t0, x0, tEnd, eps);
			h = std::min(h, tEnd - t0);

			// Output management
			Real tNextOutput = t0;
			int outputIdx = 0;
			result.solution.fillValues(outputIdx++, t0, x0);
			tNextOutput += outputInterval;

			const int maxEventIterations = 100; // Safety limit for root finding

			// Main integration loop
			while (t < tEnd - Constants::Eps) {
				if (t + h > tEnd) {
					h = tEnd - t;
				}

				StepResult stepResult = _stepper.doStep(t, x, dxdt, h, eps);
				_stats.recordStep(stepResult);

				if (stepResult.accepted) {
					Real tNew = t + stepResult.hDone;

					// Evaluate event functions at new time
					Vector<Real> gNew(numEvents);
					sys.eventFunctions(tNew, x, gNew);

					// Check for zero crossings
					bool eventDetected = false;
					int eventIndex = -1;
					Real eventTime = tNew;

					for (int i = 0; i < numEvents; ++i) {
						bool crossing = false;
						EventDirection actualDir = EventDirection::Both;

						// Check sign change
						if (gPrev[i] * gNew[i] < 0) {
							// Determine crossing direction
							if (gPrev[i] < 0 && gNew[i] > 0) {
								actualDir = EventDirection::Increasing;
							} else {
								actualDir = EventDirection::Decreasing;
							}

							// Check if this direction should trigger
							EventDirection wantDir = sys.getEventDirection(i);
							if (wantDir == EventDirection::Both || wantDir == actualDir) {
								crossing = true;
							}
						}

						if (crossing) {
							// Find exact crossing time using bisection with dense output
							Real tLo = t;
							Real tHi = tNew;
							Real gLo = gPrev[i];
							Real gHi = gNew[i];

							for (int iter = 0; iter < maxEventIterations; ++iter) {
								if (tHi - tLo < eventTol)
									break;

								// Illinois method (modified regula falsi)
								Real tMid = tLo - gLo * (tHi - tLo) / (gHi - gLo);
								// Clamp to interval
								tMid = std::max<Real>(tLo + Real(0.1) * (tHi - tLo), std::min<Real>(tHi - Real(0.1) * (tHi - tLo), tMid));

								Vector<Real> xMid = _stepper.interpolate(tMid);
								Real gMid = sys.eventFunction(i, tMid, xMid);

								if (std::abs(gMid) < eventTol) {
									tLo = tMid;
									break;
								}

								if (gLo * gMid < 0) {
									tHi = tMid;
									gHi = gMid;
								} else {
									tLo = tMid;
									gLo = gMid;
								}
							}

							Real tEvent = tLo;
							if (tEvent < eventTime) {
								eventTime = tEvent;
								eventIndex = i;
								eventDetected = true;
							}
						}
					}

					if (eventDetected && eventIndex >= 0) {
						// Get precise state at event
						Vector<Real> xEvent = _stepper.interpolate(eventTime);
						Real gEvent = sys.eventFunction(eventIndex, eventTime, xEvent);

						EventDirection actualDir;
						if (gPrev[eventIndex] < gNew[eventIndex]) {
							actualDir = EventDirection::Increasing;
						} else {
							actualDir = EventDirection::Decreasing;
						}

						// Record event
						EventInfo evt(eventIndex, eventTime, xEvent, actualDir, gEvent);
						result.events.push_back(evt);

						// Dense output up to event time
						while (tNextOutput <= eventTime + Constants::Eps && outputIdx < numOutputPoints) {
							if (tNextOutput <= eventTime) {
								Vector<Real> xInterp = _stepper.interpolate(tNextOutput);
								result.solution.fillValues(outputIdx++, tNextOutput, xInterp);
							}
							tNextOutput += outputInterval;
						}

						// Determine action
						EventAction action = sys.getEventAction(eventIndex);

						if (action == EventAction::Stop) {
							result.terminatedByEvent = true;
							result.finalTime = eventTime;
							result.finalState = xEvent;
							return result;
						} else if (action == EventAction::Restart) {
							// Apply state update
							sys.handleEvent(eventIndex, eventTime, xEvent);

							// Reset integration from event point
							t = eventTime;
							x = xEvent;
							sys.derivs(t, x, dxdt);
							sys.eventFunctions(t, x, gPrev);
							_stepper.resetFSAL();

							// Re-estimate step size
							h = estimateInitialStep(t, x, tEnd, eps);
							h = std::min(h, tEnd - t);
							continue; // Skip normal update
						}
						// EventAction::Continue: fall through to normal update
					}

					// Dense output for this step
					while (tNextOutput <= tNew + Constants::Eps && outputIdx < numOutputPoints) {
						if (tNextOutput <= tNew) {
							Vector<Real> xInterp = _stepper.interpolate(tNextOutput);
							result.solution.fillValues(outputIdx++, tNextOutput, xInterp);
						}
						tNextOutput += outputInterval;
					}

					// Update for next step
					t = tNew;
					h = stepResult.hNext;
					gPrev = gNew;
				} else {
					h = stepResult.hNext;
				}

				if (h < Constants::Eps) {
					throw ODESolverError("Step size too small in integrateWithEvents");
				}
			}

			// Ensure final point
			if (outputIdx < numOutputPoints) {
				result.solution.fillValues(outputIdx++, tEnd, x);
			}

			result.finalTime = t;
			result.finalState = x;
			return result;
		}
	};

	// Convenient type aliases
	using DormandPrince5Integrator = ODEAdaptiveIntegrator<DormandPrince5_Stepper>;
	using CashKarpIntegrator = ODEAdaptiveIntegrator<CashKarp_Stepper>;
	using DormandPrince8Integrator = ODEAdaptiveIntegrator<DormandPrince8_Stepper>;
	using BulirschStoerIntegrator = ODEAdaptiveIntegrator<BulirschStoer_Stepper>;
	using BulirschStoerRationalIntegrator = ODEAdaptiveIntegrator<BulirschStoerRational_Stepper>;
	using BulirschStoerBulirschSeqIntegrator = ODEAdaptiveIntegrator<BulirschStoerRational_Stepper>; // Alias with clearer name

} // namespace MML

#endif // MML_ODE_ADAPTIVE_INTEGRATOR_H
