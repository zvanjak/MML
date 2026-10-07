///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DynamicalAnalysisCommon.h                                          ///
///  Description: Shared integration and instrumentation support for dynamical       ///
///               system analysis                                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_ANALYSIS_COMMON_H
#define MML_DYNAMICAL_ANALYSIS_COMMON_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Vector/Vector.h>
#include <mml/interfaces/IDynamicalSystem.h>
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>

#include <exception>
#include <utility>

namespace MML::Systems
{
	inline constexpr Real DefaultTransientTime = REAL(100.0);
	inline constexpr Real DefaultRecordTime = REAL(50.0);

	namespace DynamicalAnalysisDetail
	{
		inline Vector<Real> StateAt(const ODESystemSolution& solution, int point)
		{
			Vector<Real> state(solution.getSysDim());
			for (int component = 0; component < solution.getSysDim(); ++component)
				state[component] = solution.getXValue(point, component);
			return state;
		}

		inline ODESystemSolution Integrate(IDynamicalSystem& system, const Vector<Real>& initialState,
			Real t0, Real tEnd, Real outputInterval, Real initialStep)
		{
			ODEAdaptiveIntegrator<> integrator(system);
			return integrator.integrate(initialState, t0, tEnd, outputInterval,
				Precision::ODEDefaultTolerance, initialStep);
		}
	}

	struct DynSysConfig : public EvaluationConfigBase
	{
	};

	namespace DynSysDetail
	{
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteDynSysDetailed(const char* algorithmName,
			const DynSysConfig& config, ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;
				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithmName);
				compute(result);
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const std::exception& exception) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, exception.what(), algorithmName);
			}
		}
	}
}

#endif // MML_DYNAMICAL_ANALYSIS_COMMON_H