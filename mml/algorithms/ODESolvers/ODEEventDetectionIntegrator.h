///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODEEventDetectionIntegrator.h                                       ///
///  Description: Adaptive ODE integration with event detection (zero-crossing)       ///
///               Extends ODEAdaptiveIntegrator with bisection-based event location   ///
///               Supports: state modification, terminal events, direction filtering  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_EVENT_DETECTION_INTEGRATOR_H
#define MML_ODE_EVENT_DETECTION_INTEGRATOR_H

#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "mml/interfaces/IODESystemWithEvents.h"

namespace MML {

	//===================================================================================
	//                           EventResult
	//===================================================================================

	/**
	 * @brief Result of integration with event detection.
	 * 
	 * Contains both the solution trajectory and information about detected events.
	 */
	struct EventResult {
		ODESystemSolution solution;        ///< Solution trajectory
		std::vector<EventInfo> events;     ///< All detected events in chronological order
		bool terminatedByEvent;            ///< True if integration ended due to EventAction::Stop
		Real finalTime;                    ///< Actual final time (may differ from tEnd if stopped)
		Vector<Real> finalState;           ///< Final state vector

		EventResult(Real t0, Real tEnd, int dim, int numPoints)
				: solution(t0, tEnd, dim, numPoints)
				, terminatedByEvent(false)
				, finalTime(0) {}
	};

	//===================================================================================
	//                        ODEEventDetectionIntegrator
	//===================================================================================
	/// @brief Adaptive ODE integrator with event detection support
	///
	/// Extends ODEAdaptiveIntegrator to monitor event functions g_i(t, x) during
	/// integration. When any event function crosses zero in the specified direction,
	/// the crossing is located precisely using bisection on dense output interpolation.
	///
	/// Key features:
	/// - All features of ODEAdaptiveIntegrator (adaptive step, dense output, FSAL)
	/// - Zero-crossing detection with configurable direction filtering
	/// - Bisection-based event time location using dense output
	/// - State modification via event handlers (e.g., velocity reversal)
	/// - Terminal events that stop integration
	///
	/// @tparam Stepper Adaptive stepper type (default: DormandPrince5_Stepper)
	///
	/// @example
	/// BouncingBall ball(9.81, 0.9);
	/// ODEEventDetectionIntegrator<> integrator(ball);
	/// auto result = integrator.integrateWithEvents(ball, x0, 0.0, 10.0, 0.01);
	template<typename Stepper = DormandPrince5_Stepper>
	class ODEEventDetectionIntegrator : public ODEAdaptiveIntegrator<Stepper> {
	public:
		/// @brief Construct event detection integrator for given ODE system
		explicit ODEEventDetectionIntegrator(const IODESystem& sys)
				: ODEAdaptiveIntegrator<Stepper>(sys) {}

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
			this->_stats.reset();
			this->_stepper.resetFSAL();

			// Current state
			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;
			sys.derivs(t, x, dxdt);

			// Event function values at current time
			Vector<Real> gPrev(numEvents);
			sys.eventFunctions(t, x, gPrev);

			// Initial step size
			Real h = (h0 > 0) ? h0 : this->estimateInitialStep(t0, x0, tEnd, eps);
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

				StepResult stepResult = this->_stepper.doStep(t, x, dxdt, h, eps);
				this->_stats.recordStep(stepResult);

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

								Vector<Real> xMid = this->_stepper.interpolate(tMid);
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
						Vector<Real> xEvent = this->_stepper.interpolate(eventTime);
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
								Vector<Real> xInterp = this->_stepper.interpolate(tNextOutput);
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
							this->_stepper.resetFSAL();

							// Re-estimate step size
							h = this->estimateInitialStep(t, x, tEnd, eps);
							h = std::min(h, tEnd - t);
							continue; // Skip normal update
						}
						// EventAction::Continue: fall through to normal update
					}

					// Dense output for this step
					while (tNextOutput <= tNew + Constants::Eps && outputIdx < numOutputPoints) {
						if (tNextOutput <= tNew) {
							Vector<Real> xInterp = this->_stepper.interpolate(tNextOutput);
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

	// Convenient type aliases for event detection integrators
	using DormandPrince5EventIntegrator = ODEEventDetectionIntegrator<DormandPrince5_Stepper>;
	using CashKarpEventIntegrator = ODEEventDetectionIntegrator<CashKarp_Stepper>;
	using DormandPrince8EventIntegrator = ODEEventDetectionIntegrator<DormandPrince8_Stepper>;
	using BulirschStoerEventIntegrator = ODEEventDetectionIntegrator<BulirschStoer_Stepper>;

} // namespace MML

#endif // MML_ODE_EVENT_DETECTION_INTEGRATOR_H
