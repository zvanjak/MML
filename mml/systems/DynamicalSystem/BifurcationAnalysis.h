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
	class BifurcationAnalyzer
	{
	public:
		static BifurcationDiagram<Real> Sweep(IDynamicalSystem& system, int parameterIndex,
			Real parameterMinimum, Real parameterMaximum, int numberOfSteps, Vector<Real> initialState,
			int component, Real transientTime = DefaultTransientTime,
			Real recordTime = DefaultRecordTime, Real initialStep = static_cast<MML::Real>(0.01))
		{
			if (parameterIndex < 0 || parameterIndex >= system.getNumParam())
				throw ArgumentError("BifurcationAnalyzer::Sweep: parameterIndex out of range");
			if (numberOfSteps < 1)
				throw ArgumentError("BifurcationAnalyzer::Sweep: numberOfSteps must be at least 1");
			if (!std::isfinite(parameterMinimum) || !std::isfinite(parameterMaximum))
				throw ArgumentError("BifurcationAnalyzer::Sweep: parameter range must be finite");
			if (numberOfSteps > 1 && parameterMaximum == parameterMinimum)
				throw ArgumentError("BifurcationAnalyzer::Sweep: parameterMaximum must differ from parameterMinimum when numberOfSteps > 1");
			if (component < 0 || component >= system.getDim())
				throw ArgumentError("BifurcationAnalyzer::Sweep: component out of range");
			if (initialState.size() != system.getDim())
				throw VectorDimensionError("BifurcationAnalyzer::Sweep: initialState size must equal system dimension", initialState.size(), system.getDim());

			BifurcationDiagram<Real> diagram;
			diagram.parameterName = system.getParamName(parameterIndex);
			diagram.numTransientSteps = static_cast<int>(std::lround(transientTime / initialStep));
			diagram.numRecordedPoints = static_cast<int>(std::lround(recordTime / initialStep));

			// The swept parameter is restored even if the integrator throws part-way through
			// (previously the system was left at an arbitrary sweep value on exception)
			struct ParameterRestore {
				IDynamicalSystem& sys; int index; Real original;
				~ParameterRestore() { sys.setParam(index, original); }
			} restore{ system, parameterIndex, system.getParam(parameterIndex) };

			// numberOfSteps == 1: single evaluation at parameterMinimum (step would be 0/0 otherwise)
			const Real parameterStep = numberOfSteps > 1 ? (parameterMaximum - parameterMinimum) / (numberOfSteps - 1) : static_cast<MML::Real>(0.0);

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
				system, state, static_cast<MML::Real>(0.0), totalTime, totalTime, initialStep).getXValuesAtEnd();
		}

		static std::vector<Real> FindLocalMaxima(IDynamicalSystem& system, Vector<Real>& state,
			int component, Real totalTime, Real initialStep)
		{
			std::vector<Real> maxima;
			if (totalTime == 0)
				return maxima;
			const auto solution = DynamicalAnalysisDetail::Integrate(
				system, state, static_cast<MML::Real>(0.0), totalTime, initialStep, initialStep);
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
		Real recordTime = DefaultRecordTime, Real initialStep = static_cast<MML::Real>(0.01),
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
