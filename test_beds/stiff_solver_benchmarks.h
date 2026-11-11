#ifndef MML_TEST_BEDS_STIFF_SOLVER_BENCHMARKS_H
#define MML_TEST_BEDS_STIFF_SOLVER_BENCHMARKS_H

#include "stiff_ode_test_bed.h"
#if !defined(MML_SINGLE_HEADER_INCLUDE) && !defined(MML_USE_SINGLE_HEADER)
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace MML::TestBeds {

	struct StiffBenchmarkProblem {
		std::string name;
		std::function<std::unique_ptr<IODESystemWithJacobian>()> createSystem;
		Vector<Real> initial;
		Real endTime;
		Vector<Real> reference;
		Real smokeEndTime;
		Vector<Real> smokeReference;
	};

	inline std::vector<StiffBenchmarkProblem> StiffBenchmarkProblems() {
		// Full endpoints: Bari/CWI IVP test set, rober.f, vdpol.f and hires.f.
		// vdpol uses tau=t/1000 and v=1000*y[1]; convert to MML's unscaled mu=1000 system.
		// Smoke endpoints: independent SciPy Radau/LSODA, verified by generate_stiff_references.py.
		return {
			{"Robertson", [] { return std::make_unique<RobertsonODE>(); },
				RobertsonODE::getInitialCondition(), Real(1e11),
				Vector<Real>{Real(0.2083340149701255e-7), Real(0.8333360770334713e-13), Real(0.9999999791665050)},
				Real(1),
				Vector<Real>{Real(0.9664597373330047), Real(3.0746265785786866e-5), Real(0.033509516401210415)}},
			{"VanDerPol_mu1000", [] { return std::make_unique<VanDerPolStiffODE>(Real(1000)); },
				VanDerPolStiffODE::getInitialCondition(), Real(2000),
				Vector<Real>{Real(1.706167732170483), Real(-0.8928097010247975e-3)}, Real(1),
				Vector<Real>{Real(1.9993333705063092), Real(-0.0006670371231732563)}},
			{"HIRES", [] { return std::make_unique<HIRES_ODE>(); },
				HIRES_ODE::getInitialCondition(), Real(321.8122),
				Vector<Real>{Real(0.7371312573325668e-3), Real(0.1442485726316185e-3),
				 Real(0.5888729740967575e-4), Real(0.1175651343283149e-2),
				 Real(0.2386356198831331e-2), Real(0.6238968252742796e-2),
				 Real(0.2849998395185769e-2), Real(0.2850001604814231e-2)}, Real(1),
				Vector<Real>{Real(0.25549269297154464), Real(0.05690878908653184), Real(0.0194580749770949),
				 Real(0.45851946967112456), Real(0.020147739125070403), Real(0.18228795775952053),
				 Real(0.005499081272420415), Real(0.00020091872757958202)}}
		};
	}

	enum class StiffBenchmarkMethod { BackwardEuler, BDF2, ROS2 };

	inline const char* StiffBenchmarkMethodName(StiffBenchmarkMethod method) {
		switch (method) {
		case StiffBenchmarkMethod::BackwardEuler: return "BackwardEuler";
		case StiffBenchmarkMethod::BDF2: return "BDF2";
		case StiffBenchmarkMethod::ROS2: return "ROS2";
		}
		throw ArgumentError("Unknown stiff benchmark method");
	}

	class CountingStiffSystem : public IODESystemWithJacobian {
		const IODESystemWithJacobian& _system;
	public:
		mutable std::uint64_t rhsCalls = 0;
		mutable std::uint64_t jacobianCalls = 0;

		explicit CountingStiffSystem(const IODESystemWithJacobian& system) : _system(system) {}
		int getDim() const override { return _system.getDim(); }
		void derivs(Real t, const Vector<Real>& y, Vector<Real>& f) const override {
			++rhsCalls;
			_system.derivs(t, y, f);
		}
		void jacobian(Real t, const Vector<Real>& y, Matrix<Real>& J) const override {
			++jacobianCalls;
			_system.jacobian(t, y, J);
		}
	};

	struct StiffBenchmarkMeasurement {
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string message;
		Real reachedTime = NoEstimate;
		Vector<Real> finalState;
		Real maxAbsoluteError = NoEstimate;
		Real maxScaledError = NoEstimate;
		std::uint64_t rhsCalls = 0;
		std::uint64_t jacobianCalls = 0;
		int acceptedSteps = 0;
		int rejectedSteps = 0;
		int luFactorizations = 0;
		int linearSolves = 0;
		double elapsedMilliseconds = 0;
	};

	inline StiffBenchmarkMeasurement RunStiffBenchmark(const StiffBenchmarkProblem& problem,
		StiffBenchmarkMethod method, const StiffSolverConfig& config, bool smoke) {
		StiffBenchmarkMeasurement measurement;
		const Real end = smoke ? problem.smokeEndTime : problem.endTime;
		const auto& reference = smoke ? problem.smokeReference : problem.reference;
		auto original = problem.createSystem();
		CountingStiffSystem system(*original);
		const auto start = std::chrono::steady_clock::now();
		auto result = [&]() {
			switch (method) {
			case StiffBenchmarkMethod::BackwardEuler:
				return SolveBackwardEuler(system, 0, problem.initial, end, config);
			case StiffBenchmarkMethod::BDF2:
				return SolveBDF2(system, 0, problem.initial, end, config);
			case StiffBenchmarkMethod::ROS2: {
				ROS2Solver solver(system);
				return solver.Solve(problem.initial, 0, end, config);
			}
			}
			throw ArgumentError("Unknown stiff benchmark method");
		}();
		measurement.elapsedMilliseconds = std::chrono::duration<double, std::milli>(
			std::chrono::steady_clock::now() - start).count();
		measurement.status = result.status;
		measurement.message = result.error_message;
		measurement.rhsCalls = system.rhsCalls;
		measurement.jacobianCalls = system.jacobianCalls;
		measurement.acceptedSteps = result.solution.getNumStepsOK();
		measurement.rejectedSteps = result.solution.getNumStepsBad();
		measurement.luFactorizations = result.lu_factorizations;
		measurement.linearSolves = result.linear_solves;
		if (!result.solution.isEmpty()) {
			measurement.reachedTime = result.solution.getTValue(result.solution.size() - 1);
			measurement.finalState = result.solution.getXValuesAtEnd();
		}
		else if (measurement.status != AlgorithmStatus::Success) {
			measurement.acceptedSteps = -1;
			measurement.rejectedSteps = -1;
		}
		if (measurement.status != AlgorithmStatus::Success)
			return measurement;
		if (measurement.finalState.size() != reference.size() ||
			!IsNearlyEqual(measurement.reachedTime, end, Precision::EqualityTolerance)) {
			measurement.status = AlgorithmStatus::AlgorithmSpecificFailure;
			measurement.message = "Solver reported success without reaching the reference endpoint";
			return measurement;
		}
		measurement.maxAbsoluteError = 0;
		measurement.maxScaledError = 0;
		for (int i = 0; i < reference.size(); ++i) {
			if (!std::isfinite(measurement.finalState[i])) {
				measurement.status = AlgorithmStatus::NumericalInstability;
				measurement.message = "Non-finite benchmark endpoint state";
				measurement.maxAbsoluteError = NoEstimate;
				measurement.maxScaledError = NoEstimate;
				return measurement;
			}
			Real error = std::abs(measurement.finalState[i] - reference[i]);
			measurement.maxAbsoluteError = std::max(measurement.maxAbsoluteError, error);
			measurement.maxScaledError = std::max(measurement.maxScaledError,
				error / (config.tolerance * (1 + std::abs(reference[i]))));
		}
		return measurement;
	}
}

#endif
