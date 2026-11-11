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
			IODESystem& _system;
			const PoincareSection<Real>& _section;
			bool _checkFinite;

		public:
			PoincareEventSystem(IODESystem& system, const PoincareSection<Real>& section, bool checkFinite)
				: _system(system), _section(section), _checkFinite(checkFinite) {}

			int getDim() const override { return _system.getDim(); }
			void derivs(Real time, const Vector<Real>& state, Vector<Real>& derivative) const override
			{
				_system.derivs(time, state, derivative);
			}
			int getNumEvents() const override { return 1; }
			Real eventFunction(int, Real, const Vector<Real>& state) const override
			{
				const Real value = state[_section.variable] - _section.value;
				if (_checkFinite && !std::isfinite(value))
					throw NumericalMethodError("Poincare section function is non-finite");
				return value;
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
			int numberOfIntersections, Real initialStep, Real maxSearchTime, Real transientTime,
			int& evaluations, bool checkFinite)
		{
			DynSysDetail::CheckInitialState(system, initialState);
			DynSysDetail::RequirePositive(initialStep, "Poincare initial step");
			DynSysDetail::RequirePositive(maxSearchTime, "Poincare search budget");
			DynSysDetail::RequireNonnegative(transientTime, "Poincare transient time");
			if (numberOfIntersections < 0)
				throw ArgumentError("ComputePoincareSection: crossing count must be nonnegative");
			if (section.variable < 0 || section.variable >= system.getDim())
				throw ArgumentError("ComputePoincareSection: section.variable is outside the system dimension");
			if (!std::isfinite(section.value) || section.direction < -1 || section.direction > 1)
				throw ArgumentError("ComputePoincareSection: finite section value and direction -1, 0 or 1 required");
			const Real searchEnd = transientTime + maxSearchTime;
			if (!std::isfinite(searchEnd) || !(searchEnd > transientTime))
				throw ArgumentError("ComputePoincareSection: search budget must advance a finite physical clock");

			PoincareSearchOutcome outcome;
			if (numberOfIntersections == 0) {
				outcome.complete = true;
				return outcome;
			}

			Vector<Real> state = initialState;
			if (transientTime > 0)
				state = Integrate(system, initialState, static_cast<MML::Real>(0.0), transientTime,
					transientTime, initialStep, evaluations, checkFinite).getXValuesAtEnd();

			CheckedODESystem checked(system, evaluations, checkFinite);
			PoincareEventSystem eventSystem(checked, section, checkFinite);
			ODEEventDetectionIntegrator<> integrator(eventSystem);
			Real time = transientTime;
			const Real searchWindow = std::max(DefaultRecordTime, initialStep);

			while (static_cast<int>(outcome.intersections.size()) < numberOfIntersections && time < searchEnd) {
				const Real windowEnd = time + std::min(searchWindow, searchEnd - time);
				if (!(windowEnd > time))
					throw NumericalMethodError("Poincare search window cannot advance floating-point time");
				const auto result = integrator.integrateWithEvents(eventSystem, state, time,
					windowEnd, windowEnd - time, Precision::ODEDefaultTolerance,
					Precision::EventTolerance, initialStep);
				if (!IsUsable(result.status)) {
					if (result.status == AlgorithmStatus::ResourceLimitExceeded)
						throw ResourceLimitError(result.error_message);
					throw NumericalMethodError("Poincare event integration failed: " + result.error_message);
				}
				DynSysDetail::CheckState(result.finalState, system.getDim(), "Poincare final state");
				if (!std::isfinite(result.finalTime) || !(result.finalTime > time) || result.finalTime > windowEnd)
					throw NumericalMethodError("Poincare integrator did not make valid physical-time progress");
				for (const auto& event : result.events) {
					DynSysDetail::CheckState(event.state, system.getDim(), "Poincare intersection");
					if (!std::isfinite(event.time) || event.time < time || event.time > result.finalTime)
						throw NumericalMethodError("Poincare intersection has an invalid event time");
					outcome.intersections.push_back(event.state);
					if (static_cast<int>(outcome.intersections.size()) == numberOfIntersections)
						break;
				}
				state = result.finalState;
				time = result.finalTime;
			}
			outcome.searchedTime = time - transientTime;
			outcome.complete = static_cast<int>(outcome.intersections.size()) >= numberOfIntersections;
			return outcome;
		}

		inline std::vector<Vector<Real>> IntegrateTrajectory(IDynamicalSystem& system,
			const Vector<Real>& initialState, Real totalTime, Real outputInterval, Real initialStep,
			int& evaluations, bool checkFinite)
		{
			DynSysDetail::CheckInitialState(system, initialState);
			DynSysDetail::RequireNonnegative(totalTime, "Trajectory total time");
			DynSysDetail::RequirePositive(outputInterval, "Trajectory output interval");
			DynSysDetail::RequirePositive(initialStep, "Trajectory initial step");
			if (totalTime == 0)
				return {initialState};
			const auto solution = Integrate(system, initialState, Real(0), totalTime, outputInterval,
				initialStep, evaluations, checkFinite);
			std::vector<Vector<Real>> trajectory;
			trajectory.reserve(solution.size());
			for (int point = 0; point < solution.size(); ++point)
				trajectory.push_back(StateAt(solution, point));
			return trajectory;
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
			int evaluations = 0;
			auto outcome = DynamicalAnalysisDetail::SearchPoincareCrossings(system, initialState, section,
				numberOfIntersections, initialStep, maxSearchTime, transientTime, evaluations, true);
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
			int evaluations = 0;
			return DynamicalAnalysisDetail::IntegrateTrajectory(system, initialState, totalTime,
				outputInterval, initialStep, evaluations, true);
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
					numberOfIntersections, initialStep, maxSearchTime, transientTime,
					result.function_evaluations, config.check_finite);
				result.intersections = std::move(outcome.intersections);
				result.searched_time = outcome.searchedTime;
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
				result.points = DynamicalAnalysisDetail::IntegrateTrajectory(system, initialState,
					totalTime, outputInterval, initialStep, result.function_evaluations, config.check_finite);
			});
	}
}

#endif // MML_PHASE_SPACE_ANALYSIS_H
