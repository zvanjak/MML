#include "../../TestPrecision.h"
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#include <mml/algorithms/DAESolvers.h>
#include <test_beds/stiff_solver_benchmarks.h>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace MML::Tests::ImplicitRetries {

class CubicODE : public IODESystemWithJacobian {
	Real _onset;
public:
	explicit CubicODE(Real onset = -1) : _onset(onset) {}
	int getDim() const override { return 1; }
	void derivs(Real t, const Vector<Real>& y, Vector<Real>& f) const override {
		f[0] = t > _onset ? -y[0] * y[0] * y[0] : Real(0);
	}
	void jacobian(Real t, const Vector<Real>& y, Matrix<Real>& J) const override {
		J(0, 0) = t > _onset ? -3 * y[0] * y[0] : Real(0);
	}
};

class CubicDAE : public IODESystemDAEWithJacobian {
	Real _onset;
public:
	explicit CubicDAE(Real onset = -1) : _onset(onset) {}
	int getDiffDim() const override { return 1; }
	int getAlgDim() const override { return 1; }
	void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
		f[0] = t > _onset ? -x[0] * x[0] * x[0] : Real(0);
	}
	void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
		g[0] = y[0] - x[0];
	}
	void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& J) const override {
		J(0, 0) = t > _onset ? -3 * x[0] * x[0] : Real(0);
	}
	void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override {
		J(0, 0) = 0;
	}
	void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override {
		J(0, 0) = -1;
	}
	void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override {
		J(0, 0) = 1;
	}
};

class TimePolynomialODE : public IODESystemWithJacobian {
public:
	int getDim() const override { return 1; }
	void derivs(Real t, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 2 * t; }
	void jacobian(Real, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
};

class TimePolynomialDAE : public CubicDAE {
public:
	void diffEqs(Real t, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override {
		f[0] = 2 * t;
	}
	void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override {
		J(0, 0) = 0;
	}
};

class ThrowingODE : public CubicODE {
public:
	void derivs(Real, const Vector<Real>&, Vector<Real>&) const override {
		throw std::runtime_error("problem callback failed");
	}
};

template<class Solution>
void RequireCubicResiduals(const Solution& solution, bool bdf2, Real onset = -1) {
	for (int i = 1; i < solution.size(); ++i) {
		Real t = solution.getTValue(i);
		Real h = t - solution.getTValue(i - 1);
		REQUIRE(h > 0);
		Real current = solution.getXValue(i, 0);
		Real previous = solution.getXValue(i - 1, 0);
		if constexpr (std::is_same_v<Solution, DAESolution>)
			REQUIRE(IsNearlyEqual(current, solution.getYValue(i, 0), Precision::NewtonTolerance));
		Real history = previous;
		Real coefficient = h;
		if (bdf2 && i > 1) {
			Real older = solution.getXValue(i - 2, 0);
			Real ratio = h / (solution.getTValue(i - 1) - solution.getTValue(i - 2));
			history += ratio * ratio / (1 + 2 * ratio) * (previous - older);
			coefficient *= (1 + ratio) / (1 + 2 * ratio);
		}
		Real derivative = t > onset ? -current * current * current : Real(0);
		REQUIRE(std::abs(current - history - coefficient * derivative) <
			Precision::NewtonTolerance * Real(2));
	}
}

TEST_CASE("Implicit retries - ODE BE and BDF2 recover without corrupting history",
	"[StiffSolvers][implicit-retry]")
{
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		CubicODE system;
		StiffSolverConfig config;
		config.step_size = 1;
		config.max_newton_iter = 4;
		config.newton_tol = Precision::NewtonTolerance;
		auto result = bdf2 ? SolveBDF2(system, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEuler(system, 0, Vector<Real>{1}, 1, config);
		INFO(result.error_message);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.getNumStepsBad() > 0);
		REQUIRE(result.solution.getNumStepsOK() > 1);
		REQUIRE(result.total_steps == result.solution.size() - 1);
		REQUIRE(IsNearlyEqual(result.solution.getTValue(result.solution.size() - 1), Real(1)));
		REQUIRE(std::abs(result.solution.getXValuesAtEnd()[0] - 1 / std::sqrt(Real(3))) <
			Precision::CoarseTolerance * Real(1000));
		RequireCubicResiduals(result.solution, bdf2);
	}
}

TEST_CASE("Implicit retries - DAE BE and BDF2 recover with consistent constraints",
	"[dae][implicit-retry]")
{
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		CubicDAE system;
		DAESolverConfig config;
		config.step_size = 1;
		config.max_newton_iter = 4;
		auto result = bdf2 ? SolveDAEBDF2(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config)
			: SolveDAEBackwardEuler(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		INFO(result.error_message);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.rejected_steps > 0);
		REQUIRE(result.accepted_steps == result.solution.size() - 1);
		REQUIRE(result.rejected_steps == result.solution.getNumStepsBad());
		REQUIRE(IsNearlyEqual(result.solution.getTValue(result.solution.size() - 1), Real(1)));
		RequireCubicResiduals(result.solution, bdf2);
	}
}

TEST_CASE("Implicit retries - floor and attempt budget preserve accepted ODE states",
	"[StiffSolvers][implicit-retry]")
{
	for (bool bdf2 : {false, true})
		for (bool floorFailure : {false, true}) {
			CAPTURE(bdf2, floorFailure);
			CubicODE system(Real(0.5));
			StiffSolverConfig config;
			config.step_size = Real(0.25);
			config.min_step_size = Real(0.125);
			config.max_newton_iter = 1;
			config.max_steps = floorFailure ? 20 : 3;
			auto result = bdf2 ? SolveBDF2(system, 0, Vector<Real>{1}, 1, config)
				: SolveBackwardEuler(system, 0, Vector<Real>{1}, 1, config);
			REQUIRE(result.status == (floorFailure ? AlgorithmStatus::NumericalInstability
				: AlgorithmStatus::MaxIterationsExceeded));
			REQUIRE_FALSE(result.error_message.empty());
			REQUIRE(result.solution.size() == 3);
			REQUIRE(result.solution.getTValue(2) == Real(0.5));
			REQUIRE(result.solution.getXValuesAtEnd()[0] == Real(1));
			REQUIRE(result.solution.getNumStepsBad() == (floorFailure ? 2 : 1));
		}
}

TEST_CASE("Implicit retries - floor and attempt budget preserve accepted DAE states",
	"[dae][implicit-retry]")
{
	for (bool bdf2 : {false, true})
		for (bool floorFailure : {false, true}) {
			CAPTURE(bdf2, floorFailure);
			CubicDAE system(Real(0.5));
			DAESolverConfig config;
			config.step_size = Real(0.25);
			config.min_step_size = Real(0.125);
			config.max_newton_iter = 1;
			config.max_steps = floorFailure ? 20 : 3;
			auto result = bdf2 ? SolveDAEBDF2(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config)
				: SolveDAEBackwardEuler(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
			REQUIRE(result.status == (floorFailure ? AlgorithmStatus::NumericalInstability
				: AlgorithmStatus::MaxIterationsExceeded));
			REQUIRE(result.failure_reason == (floorFailure ? DAEFailureReason::NewtonFailure
				: DAEFailureReason::MaxStepsExceeded));
			REQUIRE_FALSE(result.error_message.empty());
			REQUIRE(result.solution.size() == 3);
			REQUIRE(result.solution.getTValue(2) == Real(0.5));
			REQUIRE(result.solution.getXValuesAtEnd()[0] == Real(1));
			REQUIRE(result.solution.getYValuesAtEnd()[0] == Real(1));
			REQUIRE(result.rejected_steps == (floorFailure ? 2 : 1));
		}
}

TEST_CASE("Implicit retries - BDF2 uses actual clipped endpoint step length",
	"[StiffSolvers][dae][implicit-retry]")
{
	// Bootstrap at h=1/4: x1=1/8. Then h=1/8, r=1/2 gives x2=27/128.
	TimePolynomialODE ode;
	TimePolynomialDAE dae;
	StiffSolverConfig odeConfig;
	odeConfig.step_size = Real(0.25);
	DAESolverConfig daeConfig;
	daeConfig.step_size = odeConfig.step_size;
	auto odeResult = SolveBDF2(ode, 0, Vector<Real>{0}, Real(0.375), odeConfig);
	auto daeResult = SolveDAEBDF2(dae, 0, Vector<Real>{0}, Vector<Real>{0}, Real(0.375), daeConfig);
	REQUIRE(odeResult.status == AlgorithmStatus::Success);
	REQUIRE(daeResult.status == AlgorithmStatus::Success);
	REQUIRE(odeResult.solution.size() == 3);
	REQUIRE(daeResult.solution.size() == 3);
	REQUIRE(IsNearlyEqual(odeResult.solution.getXValuesAtEnd()[0], Real(27) / Real(128)));
	REQUIRE(IsNearlyEqual(daeResult.solution.getXValuesAtEnd()[0], Real(27) / Real(128)));
}

TEST_CASE("Implicit retries - rejected later steps preserve ODE and DAE BDF2 history",
	"[StiffSolvers][dae][implicit-retry]")
{
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		CubicODE ode(Real(0.5));
		CubicDAE dae(Real(0.5));
		StiffSolverConfig odeConfig;
		odeConfig.step_size = Real(0.25);
		odeConfig.max_newton_iter = 3;
		DAESolverConfig daeConfig;
		daeConfig.step_size = odeConfig.step_size;
		daeConfig.max_newton_iter = odeConfig.max_newton_iter;
		auto odeResult = bdf2 ? SolveBDF2(ode, 0, Vector<Real>{1}, 1, odeConfig)
			: SolveBackwardEuler(ode, 0, Vector<Real>{1}, 1, odeConfig);
		auto daeResult = bdf2 ? SolveDAEBDF2(dae, 0, Vector<Real>{1}, Vector<Real>{1}, 1, daeConfig)
			: SolveDAEBackwardEuler(dae, 0, Vector<Real>{1}, Vector<Real>{1}, 1, daeConfig);
		INFO(odeResult.error_message);
		INFO(daeResult.error_message);
		REQUIRE(odeResult.status == AlgorithmStatus::Success);
		REQUIRE(daeResult.status == AlgorithmStatus::Success);
		REQUIRE(odeResult.solution.getNumStepsBad() > 0);
		REQUIRE(daeResult.rejected_steps > 0);
		REQUIRE(odeResult.solution.getTValue(2) == Real(0.5));
		REQUIRE(daeResult.solution.getTValue(2) == Real(0.5));
		REQUIRE(odeResult.solution.getTValue(3) - Real(0.5) < odeConfig.step_size);
		REQUIRE(daeResult.solution.getTValue(3) - Real(0.5) < daeConfig.step_size);
		RequireCubicResiduals(odeResult.solution, bdf2, Real(0.5));
		RequireCubicResiduals(daeResult.solution, bdf2, Real(0.5));
	}
}

TEST_CASE("Implicit retries - legacy ODE and detailed APIs share recovery and terminal status",
	"[StiffSolvers][implicit-retry]")
{
	CubicODE system;
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		auto solution = bdf2 ? SolveBDF2(system, 0, Vector<Real>{1}, 1, Real(1), 4)
			: SolveBackwardEuler(system, 0, Vector<Real>{1}, 1, Real(1), 4);
		REQUIRE(solution.getNumStepsBad() > 0);
		REQUIRE(IsNearlyEqual(solution.getTValue(solution.size() - 1), Real(1)));
		ODEStiffConfig config;
		config.solver.step_size = Real(0.25);
		config.solver.min_step_size = Real(0.125);
		config.solver.max_newton_iter = 1;
		CubicODE switched(Real(0.5));
		auto result = bdf2 ? SolveBDF2Detailed(switched, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEulerDetailed(switched, 0, Vector<Real>{1}, 1, config);
		REQUIRE(result.status == AlgorithmStatus::NumericalInstability);
		REQUIRE(result.solution.has_value());
		REQUIRE(result.solution->getTValue(result.solution->size() - 1) == Real(0.5));
		REQUIRE(result.total_steps == 2);
		REQUIRE(result.newton_iterations > 0);
		REQUIRE(result.jacobian_evaluations > 0);
		REQUIRE(result.function_evaluations >= result.jacobian_evaluations);
		REQUIRE_THROWS_AS(bdf2 ? SolveBDF2(switched, 0, Vector<Real>{1}, 1, Real(1), 1,
			Precision::DivisionSafetyThreshold)
			: SolveBackwardEuler(switched, 0, Vector<Real>{1}, 1, Real(1), 1,
			Precision::DivisionSafetyThreshold), ODESolverError);
	}
}

TEST_CASE("Implicit retries - HIRES full reference interval recovers the recorded Newton failure",
	"[StiffSolvers][implicit-retry][HIRES]")
{
	TestBeds::HIRES_ODE system;
	auto initial = TestBeds::HIRES_ODE::getInitialCondition();
	const auto reference = TestBeds::StiffBenchmarkProblems().back().reference;
	StiffSolverConfig config;
	config.step_size = TestBeds::HIRES_ODE::getTypicalEndTime() / Real(200);
	config.newton_tol = Precision::ZeroTolerance;
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		auto result = bdf2 ? SolveBDF2(system, 0, initial, TestBeds::HIRES_ODE::getTypicalEndTime(), config)
			: SolveBackwardEuler(system, 0, initial, TestBeds::HIRES_ODE::getTypicalEndTime(), config);
		INFO(result.error_message);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.getNumStepsBad() > 0);
		REQUIRE(IsNearlyEqual(result.solution.getTValue(result.solution.size() - 1),
			TestBeds::HIRES_ODE::getTypicalEndTime()));
		REQUIRE(IsNearlyEqual(result.solution.getXValuesAtEnd()[6] + result.solution.getXValuesAtEnd()[7],
			initial[6] + initial[7], Precision::IntegrationTolerance));
		for (int i = 0; i < reference.size(); ++i)
			REQUIRE(IsWithinAbs(result.solution.getXValuesAtEnd()[i], reference[i],
				Precision::CoarseTolerance * (bdf2 ? Real(4) : Real(1))));
		auto fullConfig = config;
		fullConfig.reuse_jacobian = false;
		auto full = bdf2 ? SolveBDF2(system, 0, initial, TestBeds::HIRES_ODE::getTypicalEndTime(), fullConfig)
			: SolveBackwardEuler(system, 0, initial, TestBeds::HIRES_ODE::getTypicalEndTime(), fullConfig);
		REQUIRE(full.status == AlgorithmStatus::Success);
		REQUIRE(result.lu_factorizations * 4 < full.lu_factorizations * 3);
		for (int i = 0; i < reference.size(); ++i)
			REQUIRE(IsWithinAbs(result.solution.getXValuesAtEnd()[i], full.solution.getXValuesAtEnd()[i],
				Precision::ZeroTolerance * Real(100)));
	}
}

TEST_CASE("Implicit retries - unrepresentable steps stall explicitly before changing state",
	"[StiffSolvers][dae][implicit-retry]")
{
	CubicODE ode;
	CubicDAE dae;
	Real start = Real(1) / Precision::Eps;
	Real end = start + Real(4);
	StiffSolverConfig odeConfig;
	odeConfig.step_size = Precision::StiffMinStepSize;
	DAESolverConfig daeConfig;
	daeConfig.step_size = odeConfig.step_size;
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		auto odeResult = bdf2 ? SolveBDF2(ode, start, Vector<Real>{1}, end, odeConfig)
			: SolveBackwardEuler(ode, start, Vector<Real>{1}, end, odeConfig);
		auto daeResult = bdf2 ? SolveDAEBDF2(dae, start, Vector<Real>{1}, Vector<Real>{1}, end, daeConfig)
			: SolveDAEBackwardEuler(dae, start, Vector<Real>{1}, Vector<Real>{1}, end, daeConfig);
		REQUIRE(odeResult.status == AlgorithmStatus::Stalled);
		REQUIRE(daeResult.status == AlgorithmStatus::Stalled);
		REQUIRE(daeResult.failure_reason == DAEFailureReason::StepSizeUnderflow);
		REQUIRE_FALSE(odeResult.error_message.empty());
		REQUIRE_FALSE(daeResult.error_message.empty());
		REQUIRE(odeResult.solution.size() == 1);
		REQUIRE(daeResult.solution.size() == 1);
	}
}

TEST_CASE("Implicit retries - failed bootstrap retains only the initial ODE and DAE point",
	"[StiffSolvers][dae][implicit-retry]")
{
	CubicODE ode;
	CubicDAE dae;
	StiffSolverConfig odeConfig;
	odeConfig.step_size = 1;
	odeConfig.min_step_size = Real(0.5);
	odeConfig.max_newton_iter = 1;
	DAESolverConfig daeConfig;
	daeConfig.step_size = odeConfig.step_size;
	daeConfig.min_step_size = odeConfig.min_step_size;
	daeConfig.max_newton_iter = odeConfig.max_newton_iter;
	for (bool bdf2 : {false, true}) {
		auto odeResult = bdf2 ? SolveBDF2(ode, 0, Vector<Real>{1}, 1, odeConfig)
			: SolveBackwardEuler(ode, 0, Vector<Real>{1}, 1, odeConfig);
		auto daeResult = bdf2 ? SolveDAEBDF2(dae, 0, Vector<Real>{1}, Vector<Real>{1}, 1, daeConfig)
			: SolveDAEBackwardEuler(dae, 0, Vector<Real>{1}, Vector<Real>{1}, 1, daeConfig);
		REQUIRE(odeResult.status == AlgorithmStatus::NumericalInstability);
		REQUIRE(daeResult.failure_reason == DAEFailureReason::NewtonFailure);
		REQUIRE(odeResult.solution.size() == 1);
		REQUIRE(daeResult.solution.size() == 1);
		REQUIRE(odeResult.solution.getNumStepsBad() == 2);
		REQUIRE(daeResult.rejected_steps == 2);
		REQUIRE(odeResult.solution.getXValuesAtEnd()[0] == Real(1));
		REQUIRE(daeResult.solution.getYValuesAtEnd()[0] == Real(1));
	}
}

TEST_CASE("Implicit retries - huge intervals obey attempt budgets without span-sized allocation",
	"[StiffSolvers][dae][implicit-retry]")
{
	TimePolynomialODE ode;
	TimePolynomialDAE dae;
	StiffSolverConfig odeConfig;
	odeConfig.step_size = Real(0.25);
	odeConfig.max_steps = 2;
	DAESolverConfig daeConfig;
	daeConfig.step_size = odeConfig.step_size;
	daeConfig.max_steps = odeConfig.max_steps;
	for (bool bdf2 : {false, true}) {
		Real end = Real(1) / Precision::Eps;
		auto odeResult = bdf2 ? SolveBDF2(ode, 0, Vector<Real>{0}, end, odeConfig)
			: SolveBackwardEuler(ode, 0, Vector<Real>{0}, end, odeConfig);
		auto daeResult = bdf2 ? SolveDAEBDF2(dae, 0, Vector<Real>{0}, Vector<Real>{0}, end, daeConfig)
			: SolveDAEBackwardEuler(dae, 0, Vector<Real>{0}, Vector<Real>{0}, end, daeConfig);
		REQUIRE(odeResult.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(daeResult.failure_reason == DAEFailureReason::MaxStepsExceeded);
		REQUIRE(odeResult.solution.size() == 3);
		REQUIRE(daeResult.solution.size() == 3);
		REQUIRE(odeResult.solution.getTValue(2) == Real(0.5));
		REQUIRE(daeResult.solution.getTValue(2) == Real(0.5));
	}
}

TEST_CASE("Implicit retries - callback exceptions obey the detailed ODE exception policy",
	"[StiffSolvers][implicit-retry]")
{
	ThrowingODE system;
	for (bool bdf2 : {false, true}) {
		ODEStiffConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		auto result = bdf2 ? SolveBDF2Detailed(system, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEulerDetailed(system, 0, Vector<Real>{1}, 1, config);
		REQUIRE(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
		REQUIRE(result.error_message == "problem callback failed");
		config.exception_policy = EvaluationExceptionPolicy::Propagate;
		REQUIRE_THROWS_AS(bdf2 ? SolveBDF2Detailed(system, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEulerDetailed(system, 0, Vector<Real>{1}, 1, config), std::runtime_error);
	}
}

TEST_CASE("Implicit retries - ODE sub-step span and non-finite inputs are not successful",
	"[StiffSolvers][implicit-retry]")
{
	CubicODE system;
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		StiffSolverConfig config;
		config.step_size = 1;
		config.max_newton_iter = 10;
		auto result = bdf2 ? SolveBDF2(system, 0, Vector<Real>{1}, Real(0.1), config)
			: SolveBackwardEuler(system, 0, Vector<Real>{1}, Real(0.1), config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.size() == 2);
		REQUIRE(result.solution.getTValue(1) == Real(0.1));
		config.step_size = std::numeric_limits<Real>::quiet_NaN();
		auto invalid = bdf2 ? SolveBDF2(system, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEuler(system, 0, Vector<Real>{1}, 1, config);
		REQUIRE(invalid.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(invalid.error_message.empty());
	}
}

TEST_CASE("Implicit Newton reuse - ODE work counters and roots match full Newton",
	"[StiffSolvers][implicit-newton]")
{
	CubicODE original;
	for (bool bdf2 : {false, true}) {
		CAPTURE(bdf2);
		TestBeds::CountingStiffSystem counted(original);
		StiffSolverConfig config;
		config.step_size = Real(0.125);
		config.max_newton_iter = 20;
		config.newton_tol = Precision::ZeroTolerance;
		auto reused = bdf2 ? SolveBDF2(counted, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEuler(counted, 0, Vector<Real>{1}, 1, config);
		REQUIRE(reused.status == AlgorithmStatus::Success);
		REQUIRE(reused.function_evaluations == counted.rhsCalls);
		REQUIRE(reused.jacobian_evaluations == counted.jacobianCalls);
		config.reuse_jacobian = false;
		auto full = bdf2 ? SolveBDF2(original, 0, Vector<Real>{1}, 1, config)
			: SolveBackwardEuler(original, 0, Vector<Real>{1}, 1, config);
		REQUIRE(full.status == AlgorithmStatus::Success);
		REQUIRE(reused.lu_factorizations < reused.linear_solves);
		REQUIRE(reused.lu_factorizations < full.lu_factorizations);
		REQUIRE(full.lu_factorizations == full.linear_solves);
		REQUIRE(reused.jacobian_evaluations == reused.lu_factorizations);
		REQUIRE(IsWithinAbs(reused.solution.getXValuesAtEnd()[0], full.solution.getXValuesAtEnd()[0],
			Precision::ZeroTolerance * Real(100)));
	}
}

TEST_CASE("Implicit Newton reuse - stiff Van der Pol preserves reference accuracy with less LU work",
	"[StiffSolvers][stiff-benchmark][implicit-newton]")
{
	const auto problem = TestBeds::StiffBenchmarkProblems()[1];
	for (auto method : {TestBeds::StiffBenchmarkMethod::BackwardEuler, TestBeds::StiffBenchmarkMethod::BDF2}) {
		StiffSolverConfig config;
		config.step_size = Real(0.01);
		config.newton_tol = Precision::ZeroTolerance;
		auto reused = TestBeds::RunStiffBenchmark(problem, method, config, true);
		config.reuse_jacobian = false;
		auto full = TestBeds::RunStiffBenchmark(problem, method, config, true);
		REQUIRE(reused.status == AlgorithmStatus::Success);
		REQUIRE(full.status == AlgorithmStatus::Success);
		REQUIRE(reused.luFactorizations < full.luFactorizations);
		REQUIRE(reused.maxAbsoluteError < Precision::NewtonTolerance);
		for (int i = 0; i < problem.initial.size(); ++i)
			REQUIRE(IsWithinAbs(reused.finalState[i], full.finalState[i], Precision::NewtonTolerance));
	}
}

TEST_CASE("Implicit Newton reuse - all Newton DAE drivers preserve constrained solutions",
	"[dae][implicit-newton]")
{
	CubicDAE system;
	auto solve = [&](int method, const DAESolverConfig& config) {
		switch (method) {
		case 0: return SolveDAEBackwardEuler(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		case 1: return SolveDAEBDF2(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		case 2: return SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		case 3: return SolveDAERadauIIA(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		default: return SolveDAEBackwardEulerAdaptive(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		}
	};
	for (int method = 0; method < 5; ++method) {
		CAPTURE(method);
		DAESolverConfig config;
		config.step_size = Real(0.125);
		config.max_newton_iter = 30;
		config.newton_tol = Precision::ZeroTolerance;
		auto reused = solve(method, config);
		config.reuse_jacobian = false;
		auto full = solve(method, config);
		INFO(reused.error_message);
		INFO(full.error_message);
		REQUIRE(reused.status == AlgorithmStatus::Success);
		REQUIRE(full.status == AlgorithmStatus::Success);
		REQUIRE(reused.lu_factorizations < reused.linear_solves);
		REQUIRE(reused.lu_factorizations < full.lu_factorizations);
		REQUIRE(full.lu_factorizations == full.linear_solves);
		REQUIRE(reused.jacobian_evaluations == reused.lu_factorizations);
		REQUIRE(IsWithinAbs(reused.solution.getXValuesAtEnd()[0], full.solution.getXValuesAtEnd()[0],
			Precision::NewtonTolerance));
		for (int i = 0; i < reused.solution.size(); ++i)
			REQUIRE(IsWithinAbs(reused.solution.getXValue(i, 0), reused.solution.getYValue(i, 0),
				config.constraint_tol));
	}
}

TEST_CASE("Implicit Newton reuse - poor contraction refreshes instead of freezing a bad Jacobian",
	"[dae][implicit-newton]")
{
	Vector<Real> x{Real(0.1)}, y;
	DAESolverConfig config;
	config.max_newton_iter = 30;
	config.newton_tol = Precision::ZeroTolerance;
	auto residual = [](const Vector<Real>& xx, const Vector<Real>&, Vector<Real>& value) {
		value[0] = xx[0] * xx[0] - 2;
	};
	auto jacobian = [](const Vector<Real>& xx, const Vector<Real>&, Matrix<Real>& matrix) {
		matrix(0, 0) = 2 * xx[0];
	};
	auto result = SolveDAENewton(x, y, 1, 0, config, residual, jacobian);
	REQUIRE(result.converged);
	REQUIRE(result.lu_factorizations > 1);
	REQUIRE(result.lu_factorizations < result.linear_solves);
	REQUIRE(IsWithinAbs(x[0], std::sqrt(Real(2)), config.newton_tol));
}

TEST_CASE("Implicit Newton reuse - refresh policy respects contraction and remaining iterations",
	"[implicit-newton]")
{
	for (bool shortBudget : {false, true}) {
		Detail::ImplicitNewtonMatrix matrix;
		REQUIRE(matrix.NeedsRefresh(Real(1), Precision::CoarseTolerance, 10, true));
		matrix.Factor(Matrix<Real>::Identity(1));
		REQUIRE(matrix.NeedsRefresh(Real(0.1), Precision::CoarseTolerance,
			shortBudget ? 1 : 10, true) == shortBudget);
		REQUIRE(matrix.NeedsRefresh(Real(0.01), Precision::CoarseTolerance, 10, false));
	}
}

TEST_CASE("Implicit Newton reuse - non-finite and singular matrices remain explicit failures",
	"[dae][implicit-newton]")
{
	for (bool nonfinite : {false, true}) {
		Vector<Real> x{1}, y;
		DAESolverConfig config;
		auto residual = [](const Vector<Real>&, const Vector<Real>&, Vector<Real>& value) {
			value[0] = 1;
		};
		auto jacobian = [&](const Vector<Real>&, const Vector<Real>&, Matrix<Real>& matrix) {
			matrix(0, 0) = nonfinite ? std::numeric_limits<Real>::quiet_NaN() : Real(0);
		};
		auto result = SolveDAENewton(x, y, 1, 0, config, residual, jacobian);
		REQUIRE_FALSE(result.converged);
		REQUIRE(result.failure_reason == (nonfinite ? DAEFailureReason::NonFiniteState
			: DAEFailureReason::SingularJacobian));
		REQUIRE(result.lu_factorizations == 0);
		REQUIRE(result.linear_solves == 0);
		REQUIRE_FALSE(result.message.empty());
		REQUIRE(x[0] == Real(1));
	}
	Vector<Real> x{1, 1}, y;
	DAESolverConfig config;
	auto overflowingResidual = [](const Vector<Real>&, const Vector<Real>&, Vector<Real>& value) {
		value[0] = value[1] = std::numeric_limits<Real>::max();
	};
	auto jacobian = [](const Vector<Real>&, const Vector<Real>&, Matrix<Real>& matrix) {
		matrix = Matrix<Real>::Identity(2);
	};
	auto result = SolveDAENewton(x, y, 2, 0, config, overflowingResidual, jacobian);
	REQUIRE(result.failure_reason == DAEFailureReason::NonFiniteState);
	REQUIRE(result.lu_factorizations == 0);
}

TEST_CASE("Implicit Newton reuse - linearly implicit DAE factors once for both stages",
	"[dae][linearly-implicit][implicit-newton]")
{
	CubicDAE system;
	DAESolverConfig config;
	config.step_size = Real(0.125);
	auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	INFO(result.error_message);
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.accepted_steps == 8);
	REQUIRE(result.lu_factorizations == result.accepted_steps);
	REQUIRE(result.linear_solves == 2 * result.accepted_steps);
	REQUIRE(result.jacobian_evaluations == result.accepted_steps);
}

TEST_CASE("Implicit Newton reuse - new steps and clipped coefficients get fresh matrices",
	"[StiffSolvers][dae][implicit-newton]")
{
	class DecayODE : public CubicODE {
	public:
		void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override { f[0] = -x[0]; }
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
	} ode;
	class DecayDAE : public CubicDAE {
	public:
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
			f[0] = -x[0];
		}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override {
			J(0, 0) = -1;
		}
	} dae;
	StiffSolverConfig odeConfig;
	odeConfig.step_size = Real(0.25);
	DAESolverConfig daeConfig;
	daeConfig.step_size = odeConfig.step_size;
	for (bool bdf2 : {false, true}) {
		auto odeResult = bdf2 ? SolveBDF2(ode, 0, Vector<Real>{1}, Real(0.625), odeConfig)
			: SolveBackwardEuler(ode, 0, Vector<Real>{1}, Real(0.625), odeConfig);
		auto daeResult = bdf2 ? SolveDAEBDF2(dae, 0, Vector<Real>{1}, Vector<Real>{1}, Real(0.625), daeConfig)
			: SolveDAEBackwardEuler(dae, 0, Vector<Real>{1}, Vector<Real>{1}, Real(0.625), daeConfig);
		REQUIRE(odeResult.status == AlgorithmStatus::Success);
		REQUIRE(daeResult.status == AlgorithmStatus::Success);
		REQUIRE(odeResult.lu_factorizations == 3);
		REQUIRE(daeResult.lu_factorizations == 3);
		for (int i = 1; i < odeResult.solution.size(); ++i) {
			Real h = odeResult.solution.getTValue(i) - odeResult.solution.getTValue(i - 1);
			Real previous = odeResult.solution.getXValue(i - 1, 0);
			Real history = previous, coefficient = h;
			if (bdf2 && i > 1) {
				Real ratio = h / (odeResult.solution.getTValue(i - 1) - odeResult.solution.getTValue(i - 2));
				history += (previous - odeResult.solution.getXValue(i - 2, 0)) * ratio * ratio / (1 + 2 * ratio);
				coefficient *= (1 + ratio) / (1 + 2 * ratio);
			}
			Real expected = history / (1 + coefficient);
			REQUIRE(IsWithinAbs(odeResult.solution.getXValue(i, 0), expected, odeConfig.newton_tol));
			REQUIRE(IsWithinAbs(daeResult.solution.getXValue(i, 0), expected, daeConfig.newton_tol));
		}
	}
}

}
