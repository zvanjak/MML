///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  Backward Euler DAE integration with zero-crossing event handling                ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_EVENT_DETECTION_H
#define MML_DAE_EVENT_DETECTION_H

#include <mml/algorithms/DAESolvers/DAEBackwardEuler.h>
#include <mml/interfaces/IODESystemDAEWithEvents.h>

#include <vector>

namespace MML {

	struct DAEEventConfig {
		Real event_tolerance = PrecisionValues<Real>::IntegrationTolerance;
		int max_root_iterations = 80;
		bool recompute_consistent_ic = true;
	};

	struct DAEEventInfo {
		int event_index = -1;
		Real time = 0.0;
		Vector<Real> differential_state;
		Vector<Real> algebraic_state;
		EventDirection direction = EventDirection::Both;
		Real event_value = 0.0;
	};

	struct DAEEventResult {
		DAESolverResult integration;
		std::vector<DAEEventInfo> events;
		bool terminated_by_event = false;
		Real final_time = 0.0;
		Vector<Real> final_differential_state;
		Vector<Real> final_algebraic_state;

		DAEEventResult(Real t0, Real tEnd, int diffDim, int algDim, int capacity)
			: integration(t0, tEnd, diffDim, algDim, capacity), final_time(t0),
			  final_differential_state(diffDim), final_algebraic_state(algDim) {
			integration.algorithm_name = "DAEBackwardEulerEvents";
		}
	};

	namespace Detail {
		inline bool DAEEventCrossing(Real before, Real after, EventDirection wanted,
			EventDirection& actual) {
			if (!(before * after < REAL(0.0))) return false;
			actual = before < after ? EventDirection::Increasing : EventDirection::Decreasing;
			return wanted == EventDirection::Both || wanted == actual;
		}

		inline void AccumulateDAESegmentDiagnostics(DAESolverResult& target,
			const DAESolverResult& segment) {
			target.newton_iterations += segment.newton_iterations;
			target.jacobian_evaluations += segment.jacobian_evaluations;
			target.last_step_newton_iterations = segment.last_step_newton_iterations;
			target.max_step_newton_iterations = std::max(target.max_step_newton_iterations,
				segment.max_step_newton_iterations);
			target.final_residual_norm = segment.final_residual_norm;
			target.max_residual_norm = std::max(target.max_residual_norm, segment.max_residual_norm);
			target.final_constraint_norm = segment.final_constraint_norm;
			target.max_constraint_violation = std::max(target.max_constraint_violation,
				segment.max_constraint_violation);
		}
	}

	inline DAEEventResult SolveDAEBackwardEulerWithEvents(IODESystemDAEWithEvents& system,
		Real t0, const Vector<Real>& x0, const Vector<Real>& y0, Real tEnd,
		const DAESolverConfig& config = DAESolverConfig(),
		const DAEEventConfig& eventConfig = DAEEventConfig()) {
		const int diffDim = system.getDiffDim();
		const int algDim = system.getAlgDim();
		const int numEvents = system.getNumEvents();
		DAEEventResult result(t0, tEnd, diffDim, algDim, std::min(config.max_steps, 1024));
		Vector<Real> x = x0;
		Vector<Real> y = y0;
		Real t = t0;
		int savedPoint = 0;
		result.integration.solution.fillValues(savedPoint++, t, x, y);
		if (!ValidateDAEInitialState(system, t, x, y, config, result.integration)) {
			result.integration.solution.setFinalSize(0);
			return result;
		}

		Vector<Real> eventBefore(numEvents);
		system.eventFunctions(t, x, y, eventBefore);

		while (t < tEnd && result.integration.accepted_steps + result.integration.rejected_steps < config.max_steps) {
			const Real stepEnd = std::min(t + config.step_size, tEnd);
			DAESolverConfig stepConfig = config;
			stepConfig.step_size = stepEnd - t;
			stepConfig.max_steps = 2;
			DAESolverResult segment = SolveDAEBackwardEuler(system, t, x, y, stepEnd, stepConfig);
			Detail::AccumulateDAESegmentDiagnostics(result.integration, segment);
			if (segment.status != AlgorithmStatus::Success) {
				result.integration.status = segment.status;
				result.integration.failure_reason = segment.failure_reason;
				result.integration.error_message = segment.error_message;
				break;
			}

			Vector<Real> xEnd = segment.solution.getXValuesAtEnd();
			Vector<Real> yEnd = segment.solution.getYValuesAtEnd();
			Vector<Real> eventAfter(numEvents);
			system.eventFunctions(stepEnd, xEnd, yEnd, eventAfter);

			int selectedEvent = -1;
			Real selectedTime = stepEnd;
			Vector<Real> selectedX(diffDim), selectedY(algDim);
			EventDirection selectedDirection = EventDirection::Both;
			Real selectedValue = 0.0;

			for (int eventIndex = 0; eventIndex < numEvents; ++eventIndex) {
				EventDirection actualDirection;
				if (!Detail::DAEEventCrossing(eventBefore[eventIndex], eventAfter[eventIndex],
					system.getEventDirection(eventIndex), actualDirection)) continue;

				Real lo = t;
				Real hi = stepEnd;
				Real valueLo = eventBefore[eventIndex];
				Vector<Real> rootX(diffDim), rootY(algDim);
				Real rootValue = valueLo;
				for (int iteration = 0; iteration < eventConfig.max_root_iterations; ++iteration) {
					const Real mid = (lo + hi) * REAL(0.5);
					const Real fraction = (mid - t) / (stepEnd - t);
					rootX = x + (xEnd - x) * fraction;
					rootY = y + (yEnd - y) * fraction;
					ConsistentICResult ic = ComputeConsistentICDetailed(system, mid, rootX, rootY,
						config.max_newton_iter, config.constraint_tol);
					if (!ic.converged) {
						result.integration.status = AlgorithmStatus::NumericalInstability;
						result.integration.failure_reason = ic.failure_reason;
						result.integration.error_message = "DAE event root projection failed: " + ic.message;
						return result;
					}
					rootValue = system.eventFunction(eventIndex, mid, rootX, rootY);
					if (std::abs(rootValue) <= eventConfig.event_tolerance || hi - lo <= eventConfig.event_tolerance) {
						lo = hi = mid;
						break;
					}
					if (valueLo * rootValue < REAL(0.0)) {
						hi = mid;
					} else {
						lo = mid;
						valueLo = rootValue;
					}
				}
				const Real rootTime = (lo + hi) * REAL(0.5);
				if (rootTime < selectedTime) {
					selectedEvent = eventIndex;
					selectedTime = rootTime;
					selectedX = rootX;
					selectedY = rootY;
					selectedDirection = actualDirection;
					selectedValue = rootValue;
				}
			}

			result.integration.accepted_steps++;
			result.integration.solution.incrementSuccessfulSteps();
			if (selectedEvent < 0) {
				t = stepEnd;
				x = xEnd;
				y = yEnd;
				result.integration.solution.fillValues(savedPoint++, t, x, y);
				eventBefore = eventAfter;
				continue;
			}

			DAEEventInfo info;
			info.event_index = selectedEvent;
			info.time = selectedTime;
			info.differential_state = selectedX;
			info.algebraic_state = selectedY;
			info.direction = selectedDirection;
			info.event_value = selectedValue;
			result.events.push_back(info);

			EventAction action = system.getEventAction(selectedEvent);
			t = selectedTime;
			x = selectedX;
			y = selectedY;
			if (action == EventAction::Restart) {
				system.handleEvent(selectedEvent, t, x, y);
				if (eventConfig.recompute_consistent_ic) {
					ConsistentICResult ic = ComputeConsistentICDetailed(system, t, x, y,
						config.max_newton_iter, config.constraint_tol);
					if (!ic.converged) {
						result.integration.status = AlgorithmStatus::NumericalInstability;
						result.integration.failure_reason = ic.failure_reason;
						result.integration.error_message = "DAE event reinitialization failed: " + ic.message;
						break;
					}
				}
			}
			Vector<Real> constraints(algDim);
			system.algConstraints(t, x, y, constraints);
			if (!ValidateAcceptedDAEState(x, y, constraints, config, result.integration,
				"DAE event state")) break;
			result.integration.solution.fillValues(savedPoint++, t, x, y);

			if (action == EventAction::Stop) {
				result.terminated_by_event = true;
				break;
			}
			system.eventFunctions(t, x, y, eventBefore);
			eventBefore[selectedEvent] = REAL(0.0);
		}

		result.integration.solution.setFinalSize(savedPoint - 1);
		result.integration.total_steps = result.integration.accepted_steps;
		if (t < tEnd && !result.terminated_by_event && result.integration.status == AlgorithmStatus::Success) {
			result.integration.status = AlgorithmStatus::MaxIterationsExceeded;
			result.integration.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.integration.error_message = "Maximum DAE event step count reached before t_end";
		}
		result.final_time = t;
		result.final_differential_state = x;
		result.final_algebraic_state = y;
		return result;
	}

} // namespace MML

#endif // MML_DAE_EVENT_DETECTION_H
