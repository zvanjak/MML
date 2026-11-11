///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESolverEventDetection.h                                       ///
///  Description: Adaptive ODE integration with event detection (zero-crossing)       ///
///               Extends ODEAdaptiveIntegrator with bisection-based event location   ///
///               Supports: state modification, terminal events, direction filtering  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SOLVER_EVENT_DETECTION_H
#define MML_ODE_SOLVER_EVENT_DETECTION_H

#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#include <algorithm>
#include <vector>
#include <mml/interfaces/IODESystemWithEvents.h>

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
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;

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
																		Real eps = Precision::ODEDefaultTolerance, Real eventTol = Precision::EventTolerance, Real h0 = 0) {
			int n = sys.getDim();
			int numEvents = sys.getNumEvents();
			if (!std::isfinite(tEnd - t0) || !std::isfinite(outputInterval) || outputInterval <= 0)
				throw ODESolverError("Event integration requires finite endpoints and positive output spacing");
			const Real direction = tEnd >= t0 ? Real(1) : Real(-1);
			EventResult result(t0, tEnd, n, 0);
			const int numOutputPoints = static_cast<int>(std::min(
				std::ceil(std::abs(tEnd - t0) / outputInterval) + 1, Real(result.solution.max_size() + 1)));
			this->_stats.reset();
			this->_stepper.resetFSAL();

			// Current state
			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;
			const ODEStepperDetail::CountingODESystem system(sys, this->_stats.totalFuncEvals);
			ODEStepperDetail::EvaluateRHS(system, t, x, dxdt);

			// Event function values at current time
			Vector<Real> gPrev(numEvents);
			sys.eventFunctions(t, x, gPrev);

			// Initial step size
			Real h = (h0 > 0) ? h0 : this->EstimateInitialStepBounded(t0, x0, tEnd, eps, std::numeric_limits<Real>::infinity());
			h = direction * std::min(h, std::abs(tEnd - t0));

			// Output management
			Real tNextOutput = t0;
			int outputIdx = 0;
			result.solution.fillValues(outputIdx++, t0, x0);
			auto storePoint = [&](Real time, const Vector<Real>& state) {
				if (!detail::PrepareTrajectoryPoint(result.solution, outputIdx, result.status, result.error_message)) {
					this->_stats.status = result.status;
					this->_stats.error_message = result.error_message;
					result.terminatedByEvent = false;
					result.finalTime = result.solution.getTValue(result.solution.size() - 1);
					result.finalState = result.solution.getXValuesAtEnd();
					while (!result.events.empty() && direction * (result.events.back().time - result.finalTime) > 0)
						result.events.pop_back();
					return false;
				}
				result.solution.fillValues(outputIdx++, time, state);
				return true;
			};
			tNextOutput += direction * outputInterval;

			const int maxEventIterations = 100; // Safety limit for root finding

			// Main integration loop
			while (direction * (tEnd - t) > 0) {
				if (direction * (t + h - tEnd) > 0) {
					h = tEnd - t;
				}
				if (t + h == t)
					throw ODESolverError("Event step cannot advance floating-point time");

				StepResult stepResult;
				{
					typename ODEAdaptiveIntegrator<Stepper>::StepWorkRecorder work(
						this->_stats, this->_stepper, stepResult, result.solution);
					stepResult = this->_stepper.doStep(t, x, dxdt, h, eps);
					work.recordReturn();
					work.recordStep();
				}

				if (stepResult.accepted) {
					Real tNew = t + stepResult.hDone;

					// Evaluate event functions at new time
					Vector<Real> gNew(numEvents);
					sys.eventFunctions(tNew, x, gNew);

					// Locate every crossing in this step; they are processed in integration order
					// so that several Continue events within one step are all reported.
					struct LocatedCrossing {
						int index;
						Real time;
						EventDirection direction;
					};
					std::vector<LocatedCrossing> crossings;

					for (int i = 0; i < numEvents; ++i) {
						if (!EventSignChange(gPrev[i], gNew[i]))
							continue;

						EventDirection actualDir = gPrev[i] < gNew[i] ? EventDirection::Increasing
																	  : EventDirection::Decreasing;
						EventDirection wantDir = sys.getEventDirection(i);
						if (wantDir != EventDirection::Both && wantDir != actualDir)
							continue;

						// Find exact crossing time using bracketed regula falsi on dense output
						Real tLo = t;
						Real tHi = tNew;
						Real gLo = gPrev[i];
						Real gHi = gNew[i];

						for (int iter = 0; iter < maxEventIterations; ++iter) {
							if (std::abs(tHi - tLo) < eventTol)
								break;

							// Illinois method (modified regula falsi)
							Real tMid = tLo - gLo * (tHi - tLo) / (gHi - gLo);
							// Clamp to interval
							Real inner1 = tLo + Real(0.1) * (tHi - tLo);
							Real inner2 = tHi - Real(0.1) * (tHi - tLo);
							tMid = std::clamp(tMid, std::min(inner1, inner2), std::max(inner1, inner2));

							Vector<Real> xMid = this->_stepper.interpolate(tMid);
							Real gMid = sys.eventFunction(i, tMid, xMid);

							if (std::abs(gMid) < eventTol) {
								tLo = tMid;
								break;
							}

							if (EventSignChange(gLo, gMid)) {
								tHi = tMid;
								gHi = gMid;
							} else {
								tLo = tMid;
								gLo = gMid;
							}
						}

						crossings.push_back({i, tLo, actualDir});
					}

					std::stable_sort(crossings.begin(), crossings.end(),
						[direction](const LocatedCrossing& a, const LocatedCrossing& b) {
							return direction * (a.time - b.time) < 0;
						});

					bool restarted = false;
					for (const LocatedCrossing& crossing : crossings) {
						const int eventIndex = crossing.index;
						const Real eventTime = crossing.time;

						// Get precise state at event
						Vector<Real> xEvent = this->_stepper.interpolate(eventTime);
						Real gEvent = sys.eventFunction(eventIndex, eventTime, xEvent);

						// Record event
						EventInfo evt(eventIndex, eventTime, xEvent, crossing.direction, gEvent);
						result.events.push_back(evt);

						// Dense output up to event time
						while (direction * (tNextOutput - eventTime) <= Constants::Eps && outputIdx < numOutputPoints) {
							if (direction * (tNextOutput - eventTime) <= 0) {
								Vector<Real> xInterp = this->_stepper.interpolate(tNextOutput);
								if (!storePoint(tNextOutput, xInterp)) return result;
							}
							tNextOutput += direction * outputInterval;
						}

						// Determine action
						EventAction action = sys.getEventAction(eventIndex);

						if (action == EventAction::Stop) {
							result.terminatedByEvent = true;
							result.finalTime = eventTime;
							result.finalState = xEvent;
							if (result.solution.getTValue(result.solution.size() - 1) != eventTime)
								if (!storePoint(eventTime, xEvent)) return result;
							return result;
						} else if (action == EventAction::Restart) {
							// Apply state update
							sys.handleEvent(eventIndex, eventTime, xEvent);

							// Reset integration from event point; later crossings of this step
							// belong to the discarded trajectory and are re-detected afterwards.
							t = eventTime;
							x = xEvent;
							ODEStepperDetail::EvaluateRHS(system, t, x, dxdt);
							sys.eventFunctions(t, x, gPrev);
							gPrev[eventIndex] = 0;
							this->_stepper.resetFSAL();

							// Re-estimate step size
							h = this->EstimateInitialStepBounded(t, x, tEnd, eps, std::numeric_limits<Real>::infinity());
							h = direction * std::min(h, std::abs(tEnd - t));
							restarted = true;
							break;
						}
						// EventAction::Continue: proceed to the next crossing of this step
					}
					if (restarted)
						continue; // Skip normal update

					// Dense output for this step
					while (direction * (tNextOutput - tNew) <= Constants::Eps && outputIdx < numOutputPoints) {
						if (direction * (tNextOutput - tNew) <= 0) {
							Vector<Real> xInterp = this->_stepper.interpolate(tNextOutput);
							if (!storePoint(tNextOutput, xInterp)) return result;
						}
						tNextOutput += direction * outputInterval;
					}

					// Update for next step
					t = tNew;
					h = stepResult.hNext;
					gPrev = gNew;
				} else {
					h = stepResult.hNext;
				}

				if (std::abs(h) < Constants::Eps && direction * (tEnd - t) > 0) {
					throw ODESolverError("Step size too small in integrateWithEvents");
				}
			}

			// Ensure final point
			if (result.solution.getTValue(result.solution.size() - 1) != tEnd) {
				if (!storePoint(tEnd, x)) return result;
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

#endif // MML_ODE_SOLVER_EVENT_DETECTION_H
