///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        PhaseSpaceAnalysis.h                                               ///
///  Description: Trajectory integration and Poincare-section analysis               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PHASE_SPACE_ANALYSIS_H
#define MML_PHASE_SPACE_ANALYSIS_H

#include <mml/systems/DynamicalSystem/DynamicalAnalysisCommon.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemTypes.h>
#include <mml/algorithms/ODESolvers/ODESolverEventDetection.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace MML::Systems
{
	namespace DynamicalAnalysisDetail
	{
		class PoincareEventSystem final : public IODESystemWithEvents
		{
			IDynamicalSystem& _system;
			const PoincareSection<Real>& _section;

		public:
			PoincareEventSystem(IDynamicalSystem& system, const PoincareSection<Real>& section)
				: _system(system), _section(section) {}

			int getDim() const override { return _system.getDim(); }
			void derivs(Real time, const Vector<Real>& state, Vector<Real>& derivative) const override
			{
				_system.derivs(time, state, derivative);
			}
			int getNumEvents() const override { return 1; }
			Real eventFunction(int, Real, const Vector<Real>& state) const override
			{
				return state[_section.variable] - _section.value;
			}
			EventDirection getEventDirection(int) const override
			{
				if (_section.direction > 0)
					return EventDirection::Increasing;
				if (_section.direction < 0)
					return EventDirection::Decreasing;
				return EventDirection::Both;
			}
		};
	}

	namespace DynamicalAnalysisDetail
	{
		struct PoincareSearchOutcome
		{
			std::vector<Vector<Real>> intersections;
			bool complete = false;	///< numberOfIntersections crossings were found within the budget
			Real searchedTime = static_cast<MML::Real>(0.0); ///< integration time spent after the transient
		};

		/// Bounded search for section crossings: integrates in windows of DefaultRecordTime until
		/// numberOfIntersections crossings are collected or maxSearchTime has been integrated.
		inline PoincareSearchOutcome SearchPoincareCrossings(IDynamicalSystem& system,
			const Vector<Real>& initialState, const PoincareSection<Real>& section,
			int numberOfIntersections, Real initialStep, Real maxSearchTime, Real transientTime)
		{
			if (initialStep <= 0)
				throw ArgumentError("ComputePoincareSection: step size h must be positive");
			if (maxSearchTime <= 0 || !std::isfinite(maxSearchTime))
				throw ArgumentError("ComputePoincareSection: maxSearchTime must be positive and finite");
			if (transientTime < 0 || !std::isfinite(transientTime))
				throw ArgumentError("ComputePoincareSection: transientTime must be non-negative and finite");
			if (section.variable < 0 || section.variable >= system.getDim())
				throw ArgumentError("ComputePoincareSection: section.variable is outside the system dimension");

			PoincareSearchOutcome outcome;
			if (numberOfIntersections <= 0) {
				outcome.complete = true;
				return outcome;
			}

			Vector<Real> state = initialState;
			if (transientTime > 0)
				state = Integrate(system, initialState, static_cast<MML::Real>(0.0), transientTime, transientTime, initialStep).getXValuesAtEnd();

			PoincareEventSystem eventSystem(system, section);
			ODEEventDetectionIntegrator<> integrator(eventSystem);
			Real time = transientTime;
			const Real searchWindow = std::max(DefaultRecordTime, initialStep);
			const Real searchEnd = transientTime + maxSearchTime;

			while (static_cast<int>(outcome.intersections.size()) < numberOfIntersections && time < searchEnd) {
				const Real windowEnd = std::min(time + searchWindow, searchEnd);
				const auto result = integrator.integrateWithEvents(eventSystem, state, time,
					windowEnd, windowEnd - time, Precision::ODEDefaultTolerance,
					Precision::EventTolerance, initialStep);
				for (const auto& event : result.events) {
					outcome.intersections.push_back(event.state);
					if (static_cast<int>(outcome.intersections.size()) == numberOfIntersections)
						break;
				}
				state = result.finalState;
				if (!(result.finalTime > time))
					break; // no progress (should not happen with EventAction::Continue); report incomplete rather than spin
				time = result.finalTime;
			}
			outcome.searchedTime = time - transientTime;
			outcome.complete = static_cast<int>(outcome.intersections.size()) >= numberOfIntersections;
			return outcome;
		}
	}

	class PhaseSpaceAnalyzer
	{
	public:
		/// Collect numberOfIntersections crossings of the Poincaré section after a transient.
		///
		/// The search is bounded: at most maxSearchTime time units are integrated after the
		/// transient. If fewer crossings are found in that budget a ConvergenceError is thrown
		/// (the trajectory does not cross the section, or crosses it in the other direction);
		/// use ComputePoincareSectionDetailed to receive the partial set with a status instead.
		static std::vector<Vector<Real>> ComputePoincareSection(IDynamicalSystem& system,
			const Vector<Real>& initialState, const PoincareSection<Real>& section,
			int numberOfIntersections, Real initialStep = static_cast<MML::Real>(0.01),
			Real maxSearchTime = DefaultPoincareMaxSearchTime, Real transientTime = DefaultTransientTime)
		{
			auto outcome = DynamicalAnalysisDetail::SearchPoincareCrossings(system, initialState, section,
				numberOfIntersections, initialStep, maxSearchTime, transientTime);
			if (!outcome.complete)
				throw ConvergenceError("ComputePoincareSection: found " + std::to_string(outcome.intersections.size()) +
					" of " + std::to_string(numberOfIntersections) + " section crossings within maxSearchTime = " +
					std::to_string(maxSearchTime) + " (trajectory may not cross the section, or not in the requested direction)",
					static_cast<int>(outcome.intersections.size()));
			return std::move(outcome.intersections);
		}

		static std::vector<Vector<Real>> IntegrateTrajectory(IDynamicalSystem& system,
			const Vector<Real>& initialState, Real totalTime, Real outputInterval,
			Real initialStep = static_cast<MML::Real>(0.01))
		{
			if (totalTime < 0)
				throw ArgumentError("IntegrateTrajectory: tTotal must be non-negative (reverse-time integration is not supported)");
			if (initialStep <= 0 || outputInterval <= 0)
				throw ArgumentError("IntegrateTrajectory: step size h and dtOutput must be positive");

			std::vector<Vector<Real>> trajectory;
			if (totalTime == 0) {
				trajectory.push_back(initialState);
				return trajectory;
			}
			const auto solution = DynamicalAnalysisDetail::Integrate(system, initialState,
				static_cast<MML::Real>(0.0), totalTime, outputInterval, initialStep);
			trajectory.reserve(solution.size());
			for (int point = 0; point < solution.size(); ++point)
				trajectory.push_back(DynamicalAnalysisDetail::StateAt(solution, point));
			return trajectory;
		}
	};

	template<typename Type = Real>
	struct PoincareSectionResult : public EvaluationResultBase
	{
		std::vector<Vector<Type>> intersections;
		Type searched_time = Type(0); ///< Integration time spent after the transient looking for crossings
	};

	template<typename Type = Real>
	struct TrajectoryResult : public EvaluationResultBase
	{
		std::vector<Vector<Type>> points;
	};

	/// Detailed variant: never throws for an exhausted search budget. The crossings found so far are
	/// returned with status MaxIterationsExceeded and an explanatory error_message.
	inline PoincareSectionResult<Real> ComputePoincareSectionDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialState, const PoincareSection<Real>& section,
		int numberOfIntersections, Real initialStep = static_cast<MML::Real>(0.01), const DynSysConfig& config = {},
		Real maxSearchTime = DefaultPoincareMaxSearchTime, Real transientTime = DefaultTransientTime)
	{
		return DynSysDetail::ExecuteDynSysDetailed<PoincareSectionResult<Real>>(
			"PoincareSection", config, [&](PoincareSectionResult<Real>& result) {
				auto outcome = DynamicalAnalysisDetail::SearchPoincareCrossings(system, initialState, section,
					numberOfIntersections, initialStep, maxSearchTime, transientTime);
				result.intersections = std::move(outcome.intersections);
				result.searched_time = outcome.searchedTime;
				result.function_evaluations = static_cast<int>(result.intersections.size());
				if (!outcome.complete) {
					result.status = AlgorithmStatus::MaxIterationsExceeded;
					result.error_message = "Found " + std::to_string(result.intersections.size()) + " of " +
						std::to_string(numberOfIntersections) + " section crossings within maxSearchTime = " +
						std::to_string(maxSearchTime);
				}
			});
	}

	inline TrajectoryResult<Real> IntegrateTrajectoryDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialState, Real totalTime, Real outputInterval,
		Real initialStep = static_cast<MML::Real>(0.01), const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<TrajectoryResult<Real>>(
			"TrajectoryIntegration", config, [&](TrajectoryResult<Real>& result) {
				result.points = PhaseSpaceAnalyzer::IntegrateTrajectory(
					system, initialState, totalTime, outputInterval, initialStep);
				result.function_evaluations = static_cast<int>(result.points.size());
			});
	}
}

#endif // MML_PHASE_SPACE_ANALYSIS_H
