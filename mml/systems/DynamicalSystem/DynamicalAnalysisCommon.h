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
#include <mml/systems/DynamicalSystem/DynamicalSystemValidation.h>

#include <exception>
#include <utility>

namespace MML::Systems
{
	inline constexpr Real DefaultTransientTime = static_cast<MML::Real>(100.0);
	inline constexpr Real DefaultRecordTime = static_cast<MML::Real>(50.0);
	/// Upper bound on the integration time ComputePoincareSection may spend looking for crossings
	/// (after the transient). Without a bound the search never returns when the trajectory settles
	/// on one side of the section (MinimalMathLibrary-sgho).
	inline constexpr Real DefaultPoincareMaxSearchTime = static_cast<MML::Real>(10000.0);

	namespace DynamicalAnalysisDetail
	{
		class CheckedODESystem final : public IODESystem
		{
			IODESystem& _system;
			int _dimension;
			int& _evaluations;
			bool _checkFinite;
		public:
			CheckedODESystem(IODESystem& system, int& evaluations, bool checkFinite)
				: _system(system), _dimension(system.getDim()), _evaluations(evaluations), _checkFinite(checkFinite) {}
			int getDim() const override { return _dimension; }
			void derivs(Real time, const Vector<Real>& state, Vector<Real>& derivative) const override
			{
				DynSysDetail::CheckState(state, _dimension, "Dynamical integration stage", _checkFinite);
				derivative.resize(_dimension);
				std::fill(derivative.begin(), derivative.end(), NoEstimate);
				if (_evaluations == std::numeric_limits<int>::max())
					throw ResourceLimitError("Dynamical analysis RHS evaluation counter exhausted");
				++_evaluations;
				_system.derivs(time, state, derivative);
				DynSysDetail::CheckState(derivative, _dimension, "Dynamical RHS output");
			}
		};

		inline Vector<Real> StateAt(const ODESystemSolution& solution, int point)
		{
			Vector<Real> state(solution.getSysDim());
			for (int component = 0; component < solution.getSysDim(); ++component)
				state[component] = solution.getXValue(point, component);
			return state;
		}

		inline ODESystemSolution Integrate(IDynamicalSystem& system, const Vector<Real>& initialState,
			Real t0, Real tEnd, Real outputInterval, Real initialStep, int& evaluations, bool checkFinite)
		{
			CheckedODESystem checked(system, evaluations, checkFinite);
			ODEAdaptiveIntegrator<> integrator(checked);
			auto solution = integrator.integrate(initialState, t0, tEnd, outputInterval,
				Precision::ODEDefaultTolerance, initialStep);
			for (int point = 0; point < solution.size(); ++point)
				DynSysDetail::CheckState(StateAt(solution, point), system.getDim(), "Dynamical integration output");
			return solution;
		}
	}

	struct DynSysConfig : public EvaluationConfigBase
	{
		// check_finite controls optional integration-stage checks, never request,
		// callback-output, final-result or convergence/spectrum safety checks.
	};

	namespace DynSysDetail
	{
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteDynSysDetailed(const char* algorithmName,
			const DynSysConfig& config, ComputeFn&& compute)
		{
			AlgorithmTimer timer;
			ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithmName);
			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				compute(result);
			else {
				auto fail = [&](AlgorithmStatus status, const std::exception& error) {
					result.status = status;
					result.error_message = error.what();
				};
				try { compute(result); }
				catch (const SingularMatrixError& error) { fail(AlgorithmStatus::SingularMatrix, error); }
				catch (const std::invalid_argument& error) { fail(AlgorithmStatus::InvalidInput, error); }
				catch (const ResourceLimitError& error) { fail(AlgorithmStatus::ResourceLimitExceeded, error); }
				catch (const ConvergenceError& error) { fail(AlgorithmStatus::MaxIterationsExceeded, error); }
				catch (const MatrixNumericalError& error) { fail(AlgorithmStatus::NumericalInstability, error); }
				catch (const NumericalMethodError& error) { fail(AlgorithmStatus::NumericalInstability, error); }
				catch (const ODESolverError& error) { fail(AlgorithmStatus::NumericalInstability, error); }
				catch (const std::exception& error) { fail(AlgorithmStatus::AlgorithmSpecificFailure, error); }
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	}
}

#endif // MML_DYNAMICAL_ANALYSIS_COMMON_H
