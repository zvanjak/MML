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

#include <vector>

namespace MML::Systems
{
	class BifurcationAnalyzer
	{
	public:
		static BifurcationDiagram<Real> Sweep(IDynamicalSystem& system, int parameterIndex,
			Real parameterMinimum, Real parameterMaximum, int numberOfSteps, Vector<Real> initialState,
			int component, Real transientTime = DefaultTransientTime,
			Real recordTime = DefaultRecordTime, Real initialStep = REAL(0.01))
		{
			BifurcationDiagram<Real> diagram;
			diagram.parameterName = system.getParamName(parameterIndex);
			const Real originalParameter = system.getParam(parameterIndex);
			const Real parameterStep = (parameterMaximum - parameterMinimum) / (numberOfSteps - 1);

			for (int step = 0; step < numberOfSteps; ++step) {
				const Real parameter = parameterMinimum + step * parameterStep;
				system.setParam(parameterIndex, parameter);
				Vector<Real> state = initialState;
				IntegrateTo(system, state, transientTime, initialStep);
				auto maxima = FindLocalMaxima(system, state, component, recordTime, initialStep);
				diagram.parameterValues.push_back(parameter);
				diagram.attractorValues.push_back(maxima);
				initialState = state;
			}

			system.setParam(parameterIndex, originalParameter);
			return diagram;
		}

	private:
		static void IntegrateTo(IDynamicalSystem& system, Vector<Real>& state,
			Real totalTime, Real initialStep)
		{
			if (totalTime < 0)
				throw ArgumentError("IntegrateTo: tTotal must be non-negative (reverse-time integration is not supported)");
			if (initialStep <= 0)
				throw ArgumentError("IntegrateTo: step size h must be positive");
			if (totalTime == 0)
				return;
			state = DynamicalAnalysisDetail::Integrate(
				system, state, REAL(0.0), totalTime, totalTime, initialStep).getXValuesAtEnd();
		}

		static std::vector<Real> FindLocalMaxima(IDynamicalSystem& system, Vector<Real>& state,
			int component, Real totalTime, Real initialStep)
		{
			std::vector<Real> maxima;
			if (totalTime == 0)
				return maxima;
			const auto solution = DynamicalAnalysisDetail::Integrate(
				system, state, REAL(0.0), totalTime, initialStep, initialStep);
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

	template<typename Type = Real>
	struct BifurcationAnalysisResult : public EvaluationResultBase
	{
		BifurcationDiagram<Type> diagram;
	};

	inline BifurcationAnalysisResult<Real> SweepBifurcationDetailed(IDynamicalSystem& system,
		int parameterIndex, Real parameterMinimum, Real parameterMaximum, int numberOfSteps,
		Vector<Real> initialState, int component, Real transientTime = DefaultTransientTime,
		Real recordTime = DefaultRecordTime, Real initialStep = REAL(0.01),
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<BifurcationAnalysisResult<Real>>(
			"BifurcationAnalyzer", config, [&](BifurcationAnalysisResult<Real>& result) {
				result.diagram = BifurcationAnalyzer::Sweep(system, parameterIndex, parameterMinimum,
					parameterMaximum, numberOfSteps, initialState, component, transientTime, recordTime, initialStep);
				result.function_evaluations = numberOfSteps;
			});
	}
}

#endif // MML_BIFURCATION_ANALYSIS_H