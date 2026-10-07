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

	class PhaseSpaceAnalyzer
	{
	public:
		static std::vector<Vector<Real>> ComputePoincareSection(IDynamicalSystem& system,
			const Vector<Real>& initialState, const PoincareSection<Real>& section,
			int numberOfIntersections, Real initialStep = REAL(0.01))
		{
			if (initialStep <= 0)
				throw ArgumentError("ComputePoincareSection: step size h must be positive");
			std::vector<Vector<Real>> intersections;
			if (numberOfIntersections <= 0)
				return intersections;

			Vector<Real> state = DynamicalAnalysisDetail::Integrate(system, initialState, REAL(0.0),
				DefaultTransientTime, DefaultTransientTime, initialStep).getXValuesAtEnd();
			DynamicalAnalysisDetail::PoincareEventSystem eventSystem(system, section);
			ODEEventDetectionIntegrator<> integrator(eventSystem);
			Real time = DefaultTransientTime;
			const Real searchWindow = std::max(DefaultRecordTime, initialStep);

			while (static_cast<int>(intersections.size()) < numberOfIntersections) {
				const auto result = integrator.integrateWithEvents(eventSystem, state, time,
					time + searchWindow, searchWindow, Precision::ODEDefaultTolerance,
					Precision::EventTolerance, initialStep);
				for (const auto& event : result.events) {
					intersections.push_back(event.state);
					if (static_cast<int>(intersections.size()) == numberOfIntersections)
						break;
				}
				state = result.finalState;
				time = result.finalTime;
			}
			return intersections;
		}

		static std::vector<Vector<Real>> IntegrateTrajectory(IDynamicalSystem& system,
			const Vector<Real>& initialState, Real totalTime, Real outputInterval,
			Real initialStep = REAL(0.01))
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
				REAL(0.0), totalTime, outputInterval, initialStep);
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
	};

	template<typename Type = Real>
	struct TrajectoryResult : public EvaluationResultBase
	{
		std::vector<Vector<Type>> points;
	};

	inline PoincareSectionResult<Real> ComputePoincareSectionDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialState, const PoincareSection<Real>& section,
		int numberOfIntersections, Real initialStep = REAL(0.01), const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<PoincareSectionResult<Real>>(
			"PoincareSection", config, [&](PoincareSectionResult<Real>& result) {
				result.intersections = PhaseSpaceAnalyzer::ComputePoincareSection(
					system, initialState, section, numberOfIntersections, initialStep);
				result.function_evaluations = numberOfIntersections;
			});
	}

	inline TrajectoryResult<Real> IntegrateTrajectoryDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialState, Real totalTime, Real outputInterval,
		Real initialStep = REAL(0.01), const DynSysConfig& config = {})
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