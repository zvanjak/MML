///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        BifurcationAnalysis.h                                              ///
///  Description: Parameter sweeps and bifurcation-diagram generation                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BIFURCATION_ANALYSIS_H
#define MML_BIFURCATION_ANALYSIS_H

#include <mml/systems/DynamicalSystem/DynamicalAnalysisCommon.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemTypes.h>

#include <cmath>
#include <vector>

namespace MML::Systems
{
	namespace DynamicalAnalysisDetail
	{
	class BifurcationComputation
	{
	public:
		static BifurcationDiagram<Real> Sweep(IDynamicalSystem& system, int parameterIndex,
			Real parameterMinimum, Real parameterMaximum, int numberOfSteps, Vector<Real> initialState,
			int component, Real transientTime, Real recordTime, Real initialStep,
			int& evaluations, bool checkFinite)
		{
			DynSysDetail::CheckInitialState(system, initialState);
			DynSysDetail::RequireNonnegative(transientTime, "Bifurcation transient time");
			DynSysDetail::RequireNonnegative(recordTime, "Bifurcation record time");
			DynSysDetail::RequirePositive(initialStep, "Bifurcation initial step");
			DynSysDetail::CheckParameterScan(system, parameterIndex, parameterMinimum, parameterMaximum, numberOfSteps);
			if (numberOfSteps > 1 && parameterMaximum == parameterMinimum)
				throw ArgumentError("BifurcationAnalyzer::Sweep: parameterMaximum must differ from parameterMinimum when numberOfSteps > 1");
			if (component < 0 || component >= system.getDim())
				throw ArgumentError("BifurcationAnalyzer::Sweep: component out of range");
			if (initialState.size() != system.getDim())
				throw VectorDimensionError("BifurcationAnalyzer::Sweep: initialState size must equal system dimension", initialState.size(), system.getDim());
			ValidateClock(numberOfSteps, transientTime, recordTime);

			BifurcationDiagram<Real> diagram;
			diagram.parameterName = system.getParamName(parameterIndex);
			diagram.numTransientSteps = DynSysDetail::NominalIntervalCount(transientTime, initialStep);
			diagram.numRecordedPoints = DynSysDetail::NominalIntervalCount(recordTime, initialStep);
			diagram.parameterValues.reserve(numberOfSteps);
			diagram.attractorValues.reserve(numberOfSteps);
			return DynSysDetail::WithRestoredParameter(system, parameterIndex, [&] {
				Real time = 0;
				for (int step = 0; step < numberOfSteps; ++step) {
					const Real transientEnd = PhaseEnd(time, transientTime);
					const Real recordingEnd = PhaseEnd(transientEnd, recordTime);
					const Real fraction = numberOfSteps > 1 ? Real(step) / (numberOfSteps - 1) : Real(0);
					const Real parameter = std::lerp(parameterMinimum, parameterMaximum, fraction);
					system.setParam(parameterIndex, parameter);
					Vector<Real> state = initialState;
					IntegrateTo(system, state, time, transientEnd, initialStep, evaluations, checkFinite);
					auto maxima = FindLocalMaxima(system, state, component, transientEnd, recordingEnd, initialStep, evaluations, checkFinite);
					diagram.parameterValues.push_back(parameter);
					diagram.attractorValues.push_back(std::move(maxima));
					initialState = state;
					time = recordingEnd;
				}
				return std::move(diagram);
			});
		}

	private:
		static Real PhaseEnd(Real start, Real duration)
		{
			if (duration == 0) return start;
			const Real end = start + duration;
			if (!std::isfinite(end) || !(end > start))
				throw ArgumentError("Bifurcation phase duration must advance a finite physical clock");
			return end;
		}

		static void ValidateClock(int steps, Real transient, Real recording)
		{
			if (transient == 0 && recording == 0) return;
			Real time = 0;
			// A finite total span can still contain a later phase that cannot advance time.
			for (int step = 0; step < steps; ++step)
				time = PhaseEnd(PhaseEnd(time, transient), recording);
		}

		static void IntegrateTo(IDynamicalSystem& system, Vector<Real>& state,
			Real startTime, Real endTime, Real initialStep, int& evaluations, bool checkFinite)
		{
			if (endTime == startTime)
				return;
			state = DynamicalAnalysisDetail::Integrate(
				system, state, startTime, endTime, endTime - startTime, initialStep, evaluations, checkFinite).getXValuesAtEnd();
		}

		static std::vector<Real> FindLocalMaxima(IDynamicalSystem& system, Vector<Real>& state,
			int component, Real startTime, Real endTime, Real initialStep, int& evaluations, bool checkFinite)
		{
			std::vector<Real> maxima;
			if (endTime == startTime)
				return maxima;
			const auto solution = DynamicalAnalysisDetail::Integrate(
				system, state, startTime, endTime, initialStep, initialStep, evaluations, checkFinite);
			for (int point = 1; point + 1 < solution.size(); ++point) {
				const Real previous = solution.getXValue(point - 1, component);
				const Real current = solution.getXValue(point, component);
				const Real next = solution.getXValue(point + 1, component);
				if (current > previous && current > next)
					maxima.push_back(current);
			}
			state = solution.getXValuesAtEnd();
			return maxima;
		}
	};
	}

	class BifurcationAnalyzer
	{
	public:
		/// Continue state and physical time from t=0 through both phases and all parameter samples.
		static BifurcationDiagram<Real> Sweep(IDynamicalSystem& system, int parameterIndex,
			Real parameterMinimum, Real parameterMaximum, int numberOfSteps, Vector<Real> initialState,
			int component, Real transientTime = DefaultTransientTime,
			Real recordTime = DefaultRecordTime, Real initialStep = Real(0.01))
		{
			int evaluations = 0;
			return DynamicalAnalysisDetail::BifurcationComputation::Sweep(system, parameterIndex,
				parameterMinimum, parameterMaximum, numberOfSteps, std::move(initialState),
				component, transientTime, recordTime, initialStep, evaluations, true);
		}
	};

	template<typename Type = Real>
	struct BifurcationAnalysisResult : public EvaluationResultBase
	{
		BifurcationDiagram<Type> diagram;
	};

	inline BifurcationAnalysisResult<Real> SweepBifurcationDetailed(IDynamicalSystem& system,
		int parameterIndex, Real parameterMinimum, Real parameterMaximum, int numberOfSteps,
		Vector<Real> initialState, int component, Real transientTime = DefaultTransientTime,
		Real recordTime = DefaultRecordTime, Real initialStep = static_cast<MML::Real>(0.01),
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<BifurcationAnalysisResult<Real>>(
			"BifurcationAnalyzer", config, [&](BifurcationAnalysisResult<Real>& result) {
				result.diagram = DynamicalAnalysisDetail::BifurcationComputation::Sweep(system, parameterIndex,
					parameterMinimum, parameterMaximum, numberOfSteps, std::move(initialState), component,
					transientTime, recordTime, initialStep, result.function_evaluations, config.check_finite);
			});
	}
}

#endif // MML_BIFURCATION_ANALYSIS_H
