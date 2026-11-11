#include "../../TestPrecision.h"
#include <mml/algorithms/ODESolvers.h>
#include <mml/algorithms/ODESolvers/ODESolverEventDetection.h>
#include <mml/algorithms/DAESolvers.h>
#include <mml/algorithms/DAESolvers/DAEEventDetection.h>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

namespace MML::Tests::ReverseTime {

class PolynomialODE : public IODESystemWithJacobian {
public:
	int getDim() const override { return 1; }
	void derivs(Real t, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 2 * t; }
	void jacobian(Real, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
};

class PolynomialDAE : public IODESystemDAEWithJacobian {
public:
	int getDiffDim() const override { return 1; }
	int getAlgDim() const override { return 1; }
	void diffEqs(Real t, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override {
		f[0] = 2 * t;
	}
	void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
		g[0] = y[0] - 2 * x[0];
	}
	void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
	void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
	void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -2; }
	void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
};

class ExponentialODE : public PolynomialODE {
public:
	void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override { f[0] = x[0]; }
	void jacobian(Real, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
};

class ExponentialDAE : public PolynomialDAE {
public:
	void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
		f[0] = x[0];
	}
	void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
};

class CubicODE : public PolynomialODE {
public:
	void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override { f[0] = x[0] * x[0] * x[0]; }
	void jacobian(Real, const Vector<Real>& x, Matrix<Real>& J) const override { J(0, 0) = 3 * x[0] * x[0]; }
};

class CubicDAE : public PolynomialDAE {
public:
	void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
		f[0] = x[0] * x[0] * x[0];
	}
	void jacobian_fx(Real, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& J) const override {
		J(0, 0) = 3 * x[0] * x[0];
	}
};

class EventODE : public PolynomialODE, public IODESystemWithEvents {
public:
	EventDirection wanted = EventDirection::Decreasing;
	EventAction action = EventAction::Stop;
	bool multiple = false;
	mutable int restarts = 0;
	void derivs(Real, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
	int getNumEvents() const override { return multiple ? 2 : 1; }
	Real eventFunction(int i, Real t, const Vector<Real>&) const override {
		return t - (multiple && i == 0 ? Real(0.25) : Real(0.75));
	}
	EventDirection getEventDirection(int) const override { return wanted; }
	EventAction getEventAction(int) const override { return action; }
	void handleEvent(int, Real, Vector<Real>& x) const override { x[0] += 10; ++restarts; }
};

class EventDAE : public IODESystemDAEWithEvents {
	PolynomialDAE _equations;
public:
	EventDirection wanted = EventDirection::Decreasing;
	EventAction action = EventAction::Stop;
	bool multiple = false;
	mutable int restarts = 0;
	int getDiffDim() const override { return 1; }
	int getAlgDim() const override { return 1; }
	void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
	void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
		_equations.algConstraints(t, x, y, g);
	}
	void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y, Matrix<Real>& J) const override {
		_equations.jacobian_fx(t, x, y, J);
	}
	void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y, Matrix<Real>& J) const override {
		_equations.jacobian_fy(t, x, y, J);
	}
	void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y, Matrix<Real>& J) const override {
		_equations.jacobian_gx(t, x, y, J);
	}
	void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y, Matrix<Real>& J) const override {
		_equations.jacobian_gy(t, x, y, J);
	}
	int getNumEvents() const override { return multiple ? 2 : 1; }
	Real eventFunction(int i, Real t, const Vector<Real>&, const Vector<Real>&) const override {
		return t - (multiple && i == 0 ? Real(0.25) : Real(0.75));
	}
	EventDirection getEventDirection(int) const override { return wanted; }
	EventAction getEventAction(int) const override { return action; }
	void handleEvent(int, Real, Vector<Real>& x, Vector<Real>& y) const override {
		x[0] += 10;
		y[0] = -100;
		++restarts;
	}
};

template<class Solution>
void RequireReverse(const Solution& solution, Real end = 0) {
	REQUIRE(solution.size() >= 2);
	REQUIRE(IsWithinAbs(solution.getTValue(solution.size() - 1), end, Precision::IntegrationTolerance));
	for (int i = 1; i < solution.size(); ++i)
		REQUIRE(solution.getTValue(i) < solution.getTValue(i - 1));
}

template<class Stepper>
void RequireAdaptiveRoundTrip() {
	PolynomialODE system;
	ODEAdaptiveIntegrator<Stepper> integrator(system);
	auto forward = integrator.integrate(Vector<Real>{1}, 0, 1, Real(0.125));
	auto reverse = integrator.integrate(forward.getXValuesAtEnd(), 1, 0, Real(0.125));
	RequireReverse(reverse);
	for (int i = 0; i < reverse.size(); ++i) {
		Real t = reverse.getTValue(i);
		REQUIRE(IsWithinAbs(reverse.getXValue(i, 0), 1 + t * t, Precision::IntegrationTolerance));
	}
	REQUIRE(integrator.getStatistics().minStepSize > 0);
}

TEST_CASE("Reverse time - every adaptive ODE stepper round trips with descending dense output",
	"[reverse-time][ode]")
{
	SECTION("DP5") { RequireAdaptiveRoundTrip<DormandPrince5_Stepper>(); }
	SECTION("DOP853") { RequireAdaptiveRoundTrip<DormandPrince8_Stepper>(); }
	SECTION("Cash Karp") { RequireAdaptiveRoundTrip<CashKarp_Stepper>(); }
	SECTION("Bulirsch Stoer") { RequireAdaptiveRoundTrip<BulirschStoer_Stepper>(); }
	SECTION("Rational Bulirsch Stoer") { RequireAdaptiveRoundTrip<BulirschStoerRational_Stepper>(); }
}

TEST_CASE("Reverse time - implicit ODE solvers reach clipped endpoints and round trip",
	"[reverse-time][StiffSolvers]")
{
	PolynomialODE system;
	for (int method = 0; method < 3; ++method) {
		CAPTURE(method);
		StiffSolverConfig config;
		config.step_size = Real(0.001);
		config.max_steps = 5000;
		config.tolerance = Precision::CoarseTolerance;
		// The float default Newton tolerance exceeds the BDF2 predictor defect.
		if constexpr (std::is_same_v<Real, float>)
			config.newton_tol = Precision::ZeroTolerance;
		auto solve = [&](Real start, Real end, const Vector<Real>& initial) {
			if (method == 0) return SolveBackwardEuler(system, start, initial, end, config);
			if (method == 1) return SolveBDF2(system, start, initial, end, config);
			ROS2Solver solver(system);
			return solver.Solve(initial, start, end, config);
		};
		auto forward = solve(0, Real(0.625), Vector<Real>{1});
		INFO(forward.error_message);
		REQUIRE(forward.status == AlgorithmStatus::Success);
		REQUIRE(forward.solution.getTValue(forward.solution.size() - 1) == Real(0.625));
		auto reverse = solve(Real(0.625), 0, forward.solution.getXValuesAtEnd());
		INFO(reverse.error_message);
		REQUIRE(reverse.status == AlgorithmStatus::Success);
		RequireReverse(reverse.solution);
		REQUIRE(reverse.solution.getTValue(reverse.solution.size() - 1) == Real(0));
		REQUIRE(IsWithinAbs(reverse.solution.getXValuesAtEnd()[0], Real(1), Precision::CoarseTolerance * Real(20)));
	}
}

TEST_CASE("Reverse time - every DAE driver preserves order constraints and round trip",
	"[reverse-time][dae]")
{
	PolynomialDAE system;
	for (int method = 0; method < 7; ++method) {
		CAPTURE(method);
		DAESolverConfig config;
		config.step_size = Real(0.01);
		config.newton_tol = Precision::ZeroTolerance;
		config.abs_tolerance = TOL(PrecisionValues<double>::IntegrationTolerance, Precision::ZeroTolerance);
		config.rel_tolerance = TOL(100 * PrecisionValues<double>::IntegrationTolerance, Precision::ZeroTolerance);
		auto solve = [&](Real start, Real end, const Vector<Real>& x, const Vector<Real>& y) {
			switch (method) {
			case 0: return SolveDAEBackwardEuler(system, start, x, y, end, config);
			case 1: return SolveDAEBDF2(system, start, x, y, end, config);
			case 2: return SolveDAEBDF4(system, start, x, y, end, config);
			case 3: return SolveDAERadauIIA(system, start, x, y, end, config);
			case 4: return SolveDAELinearlyImplicit(system, start, x, y, end, config);
			case 5:
				return SolveDAEBackwardEulerAdaptive(system, start, x, y, end, config);
			default:
				return SolveDAERODAS4(system, start, x, y, end, config);
			}
		};
		auto forward = solve(0, Real(0.625), Vector<Real>{1}, Vector<Real>{2});
		INFO(forward.error_message);
		REQUIRE(forward.status == AlgorithmStatus::Success);
		auto reverse = solve(Real(0.625), 0, forward.solution.getXValuesAtEnd(), forward.solution.getYValuesAtEnd());
		INFO(reverse.error_message);
		REQUIRE(reverse.status == AlgorithmStatus::Success);
		RequireReverse(reverse.solution);
		REQUIRE(IsWithinAbs(reverse.solution.getXValuesAtEnd()[0], Real(1), Precision::CoarseTolerance * Real(200)));
		for (int i = 0; i < reverse.solution.size(); ++i)
			REQUIRE(IsWithinAbs(reverse.solution.getYValue(i, 0), 2 * reverse.solution.getXValue(i, 0),
				config.constraint_tol));
	}
}

TEST_CASE("Reverse time - specified output times must be strictly monotone in either direction",
	"[reverse-time][ode]")
{
	PolynomialODE system;
	DormandPrince8Integrator integrator(system);
	auto solution = integrator.integrateAt(Vector<Real>{2}, Vector<Real>{1, Real(0.75), Real(0.125), 0});
	RequireReverse(solution);
	REQUIRE(IsWithinAbs(solution.getXValuesAtEnd()[0], Real(1), Precision::IntegrationTolerance));
	REQUIRE_THROWS_AS(integrator.integrateAt(Vector<Real>{2}, Vector<Real>{1, Real(0.5), Real(0.75), 0}), ODESolverError);
	REQUIRE_THROWS_AS(integrator.integrateAt(Vector<Real>{2}, Vector<Real>{1, Real(0.5), Real(0.5), 0}), ODESolverError);
}

TEST_CASE("Reverse time - solution interpolation ignores unused capacity and supports descending times",
	"[reverse-time][ode][dae]")
{
	ODESystemSolution ode(1, 0, 1, 20);
	DAESolution dae(1, 0, 1, 1, 20);
	for (int i = 0; i < 3; ++i) {
		Real t = 1 - Real(i) / 2;
		ode.fillValues(i, t, Vector<Real>{1 + t});
		dae.fillValues(i, t, Vector<Real>{1 + t}, Vector<Real>{2 + 2 * t});
	}
	REQUIRE(IsWithinAbs(ode.getSolAsLinInterp(0)(Real(0.25)), Real(1.25), Precision::IntegrationTolerance));
	REQUIRE(IsWithinAbs(ode.getSolAsSplineInterp(0)(Real(0.25)), Real(1.25), Precision::IntegrationTolerance));
	REQUIRE(IsWithinAbs(ode.getSolAsPolyInterp(0, 2)(Real(0.25)), Real(1.25), Precision::IntegrationTolerance));
	REQUIRE(IsWithinAbs(dae.getXAsLinInterp(0)(Real(0.25)), Real(1.25), Precision::IntegrationTolerance));
	REQUIRE(IsWithinAbs(dae.getYAsSplineInterp(0)(Real(0.25)), Real(2.5), Precision::IntegrationTolerance));
}

TEST_CASE("Reverse time - fixed RK and default detailed ODE APIs preserve direction",
	"[reverse-time][ode]")
{
	PolynomialODE system;
	RungeKutta4_StepCalculator rk4;
	ODESystemFixedStepSolver fixed(system, rk4);
	auto forward = fixed.integrate(Vector<Real>{1}, 0, 1, 32);
	REQUIRE(IsWithinAbs(forward.getXValuesAtEnd()[0], Real(2), Precision::IntegrationTolerance));
	auto solution = fixed.integrate(forward.getXValuesAtEnd(), 1, 0, 32);
	RequireReverse(solution);
	REQUIRE(IsWithinAbs(solution.getXValuesAtEnd()[0], Real(1), Precision::IntegrationTolerance));
	auto adaptive = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{2}, 1, 0);
	REQUIRE(adaptive.status == AlgorithmStatus::Success);
	RequireReverse(adaptive.solution.value());
	ODEStiffConfig config;
	config.solver.tolerance = Precision::CoarseTolerance;
	auto ros = SolveROS2Detailed(system, Vector<Real>{2}, 1, 0, config);
	INFO(ros.error_message);
	REQUIRE(ros.status == AlgorithmStatus::Success);
	RequireReverse(ros.solution.value());
}

TEST_CASE("Reverse time - ODE events select first encountered root and filter along the path",
	"[reverse-time][events][ode]")
{
	EventODE system;
	system.multiple = true;
	DormandPrince8EventIntegrator integrator(system);
	auto result = integrator.integrateWithEvents(system, Vector<Real>{1}, 1, 0, Real(0.125),
		Precision::CoarseTolerance, Precision::EventTolerance, 1);
	REQUIRE(result.terminatedByEvent);
	REQUIRE(result.events.size() == 1);
	REQUIRE(result.events[0].eventIndex == 1);
	REQUIRE(result.events[0].direction == EventDirection::Decreasing);
	REQUIRE(IsWithinAbs(result.finalTime, Real(0.75), Precision::IntegrationTolerance));
	RequireReverse(result.solution, result.finalTime);
	system.wanted = EventDirection::Increasing;
	auto filtered = integrator.integrateWithEvents(system, Vector<Real>{1}, 1, 0, Real(0.125));
	REQUIRE(filtered.events.empty());
	REQUIRE(filtered.finalTime == Real(0));
}

TEST_CASE("Reverse time - ODE continue and restart events reach the endpoint once",
	"[reverse-time][events][ode]")
{
	for (auto action : {EventAction::Continue, EventAction::Restart}) {
		EventODE system;
		system.action = action;
		DormandPrince5EventIntegrator integrator(system);
		auto result = integrator.integrateWithEvents(system, Vector<Real>{1}, 1, 0, Real(0.125),
			Precision::CoarseTolerance, Precision::EventTolerance, 1);
		REQUIRE_FALSE(result.terminatedByEvent);
		REQUIRE(result.events.size() == 1);
		REQUIRE(system.restarts == (action == EventAction::Restart ? 1 : 0));
		REQUIRE(result.finalTime == Real(0));
		REQUIRE(result.finalState[0] == (action == EventAction::Restart ? Real(11) : Real(1)));
		RequireReverse(result.solution);
	}
}

TEST_CASE("Reverse time - DAE events select the first root and retain consistent restart states",
	"[reverse-time][events][dae]")
{
	for (auto action : {EventAction::Stop, EventAction::Continue, EventAction::Restart}) {
		EventDAE system;
		system.action = action;
		system.multiple = action == EventAction::Stop;
		DAESolverConfig config;
		config.step_size = 1;
		config.newton_tol = Precision::ZeroTolerance;
		auto result = SolveDAEBackwardEulerWithEvents(system, 1, Vector<Real>{1}, Vector<Real>{2}, 0, config);
		INFO(result.integration.error_message);
		REQUIRE(result.integration.status == AlgorithmStatus::Success);
		REQUIRE(result.events.size() == 1);
		REQUIRE(result.events[0].direction == EventDirection::Decreasing);
		REQUIRE(IsWithinAbs(result.events[0].time, Real(0.75), Precision::IntegrationTolerance));
		RequireReverse(result.integration.solution, action == EventAction::Stop ? result.events[0].time : Real(0));
		REQUIRE(IsWithinAbs(result.final_algebraic_state[0], 2 * result.final_differential_state[0],
			config.constraint_tol));
		REQUIRE(system.restarts == (action == EventAction::Restart ? 1 : 0));
		if (action == EventAction::Stop) REQUIRE(result.events[0].event_index == 1);
	}
}

TEST_CASE("Reverse time - state dependent ODE iteration matrices have the signed step",
	"[reverse-time][StiffSolvers]")
{
	ExponentialODE system;
	for (int method = 0; method < 3; ++method) {
		CAPTURE(method);
		StiffSolverConfig config;
		config.step_size = Real(0.001);
		config.tolerance = Precision::CoarseTolerance;
		if constexpr (std::is_same_v<Real, float>)
			config.newton_tol = Precision::ZeroTolerance;
		auto solve = [&](Real start, Real end, const Vector<Real>& initial) {
			if (method == 0) return SolveBackwardEuler(system, start, initial, end, config);
			if (method == 1) return SolveBDF2(system, start, initial, end, config);
			ROS2Solver solver(system);
			return solver.Solve(initial, start, end, config);
		};
		auto forward = solve(0, Real(0.25), Vector<Real>{1});
		REQUIRE(forward.status == AlgorithmStatus::Success);
		REQUIRE(forward.solution.getTValue(forward.solution.size() - 1) == Real(0.25));
		auto backward = solve(Real(0.25), 0, forward.solution.getXValuesAtEnd());
		REQUIRE(backward.status == AlgorithmStatus::Success);
		REQUIRE(backward.solution.getTValue(backward.solution.size() - 1) == 0);
		RequireReverse(backward.solution);
		REQUIRE(IsWithinAbs(forward.solution.getXValuesAtEnd()[0], std::exp(Real(0.25)),
			Precision::CoarseTolerance * Real(10)));
		REQUIRE(IsWithinAbs(backward.solution.getXValuesAtEnd()[0], Real(1), Precision::CoarseTolerance * Real(10)));
	}
}

TEST_CASE("Reverse time - state dependent DAEs and spans shorter than a requested step",
	"[reverse-time][dae]")
{
	ExponentialDAE system;
	for (int method = 0; method < 6; ++method)
		for (bool shortSpan : {false, true}) {
			CAPTURE(method, shortSpan);
			DAESolverConfig config;
			config.step_size = shortSpan ? Real(1) : Real(0.001);
			config.newton_tol = Precision::ZeroTolerance;
			config.abs_tolerance = TOL(PrecisionValues<double>::IntegrationTolerance, Precision::ZeroTolerance);
			config.rel_tolerance = TOL(100 * PrecisionValues<double>::IntegrationTolerance, Precision::ZeroTolerance);
			auto solve = [&](Real start, Real end, const Vector<Real>& x, const Vector<Real>& y) {
				switch (method) {
				case 0: return SolveDAEBackwardEuler(system, start, x, y, end, config);
				case 1: return SolveDAEBDF2(system, start, x, y, end, config);
				case 2: return SolveDAEBDF4(system, start, x, y, end, config);
				case 3: return SolveDAERadauIIA(system, start, x, y, end, config);
				case 4: return SolveDAELinearlyImplicit(system, start, x, y, end, config);
				default: return SolveDAEBackwardEulerAdaptive(system, start, x, y, end, config);
				}
			};
			Real end = shortSpan ? Real(0.125) : Real(0.25);
			auto forward = solve(0, end, Vector<Real>{1}, Vector<Real>{2});
			INFO(forward.error_message);
			REQUIRE(forward.status == AlgorithmStatus::Success);
			auto backward = solve(end, 0, forward.solution.getXValuesAtEnd(), forward.solution.getYValuesAtEnd());
			INFO(backward.error_message);
			REQUIRE(backward.status == AlgorithmStatus::Success);
			RequireReverse(backward.solution);
			REQUIRE(backward.total_steps == backward.accepted_steps);
			REQUIRE(IsWithinAbs(backward.solution.getXValuesAtEnd()[0], Real(1),
				Precision::CoarseTolerance * (shortSpan ? Real(300) : Real(10))));
			if (shortSpan && method == 2) REQUIRE(backward.solution.size() == 2);
		}
}

TEST_CASE("Reverse time - Newton retries retain accepted history and bounded failure statuses",
	"[reverse-time][StiffSolvers][dae]")
{
	CubicODE ode;
	CubicDAE dae;
	for (bool bdf2 : {false, true}) {
		StiffSolverConfig odeConfig;
		odeConfig.step_size = 1;
		odeConfig.max_newton_iter = 3;
		if constexpr (std::is_same_v<Real, float>)
			odeConfig.newton_tol = Precision::ZeroTolerance;
		DAESolverConfig daeConfig;
		daeConfig.step_size = 1;
		daeConfig.max_newton_iter = 3;
		auto result = bdf2 ? SolveBDF2(ode, 1, Vector<Real>{1}, 0, odeConfig)
			: SolveBackwardEuler(ode, 1, Vector<Real>{1}, 0, odeConfig);
		auto daeResult = bdf2 ? SolveDAEBDF2(dae, 1, Vector<Real>{1}, Vector<Real>{2}, 0, daeConfig)
			: SolveDAEBackwardEuler(dae, 1, Vector<Real>{1}, Vector<Real>{2}, 0, daeConfig);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(daeResult.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.getNumStepsBad() > 0);
		REQUIRE(daeResult.rejected_steps > 0);
		RequireReverse(result.solution);
		RequireReverse(daeResult.solution);
		for (int i = 1; i < result.solution.size(); ++i) {
			Real h = result.solution.getTValue(i) - result.solution.getTValue(i - 1);
			Real value = result.solution.getXValue(i, 0);
			Real previous = result.solution.getXValue(i - 1, 0);
			Real history = previous, coefficient = h;
			if (bdf2 && i > 1) {
				Real ratio = h / (result.solution.getTValue(i - 1) - result.solution.getTValue(i - 2));
				history += (previous - result.solution.getXValue(i - 2, 0)) * ratio * ratio / (1 + 2 * ratio);
				coefficient *= (1 + ratio) / (1 + 2 * ratio);
			}
			REQUIRE(std::abs(value - history - coefficient * value * value * value) < odeConfig.newton_tol * Real(2));
		}
		odeConfig.max_steps = daeConfig.max_steps = 1;
		auto failed = bdf2 ? SolveBDF2(ode, 1, Vector<Real>{1}, 0, odeConfig)
			: SolveBackwardEuler(ode, 1, Vector<Real>{1}, 0, odeConfig);
		auto daeFailed = bdf2 ? SolveDAEBDF2(dae, 1, Vector<Real>{1}, Vector<Real>{2}, 0, daeConfig)
			: SolveDAEBackwardEuler(dae, 1, Vector<Real>{1}, Vector<Real>{2}, 0, daeConfig);
		REQUIRE(failed.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(daeFailed.failure_reason == DAEFailureReason::MaxStepsExceeded);
		REQUIRE(failed.solution.size() == 1);
		REQUIRE(daeFailed.solution.size() == 1);
	}
}

TEST_CASE("Reverse time - parametric solution curves expose ordered parameter bounds",
	"[reverse-time][ode][dae]")
{
	ODESystemSolution ode(1, 0, 3, 10);
	DAESolution dae(1, 0, 3, 1, 10);
	for (int i = 0; i < 3; ++i) {
		Real t = 1 - Real(i) / 2;
		Vector<Real> x{t, 2 * t, 3 * t};
		ode.fillValues(i, t, x);
		dae.fillValues(i, t, x, Vector<Real>{t});
	}
	auto curve = ode.getSolAsParamCurve3D(0, 1, 2);
	auto daeCurve = dae.getXAsParamCurve2D(0, 1);
	REQUIRE(curve.getMinT() == Real(0));
	REQUIRE(curve.getMaxT() == Real(1));
	REQUIRE(IsWithinAbs(curve(Real(0.25))[2], Real(0.75), Precision::IntegrationTolerance));
	REQUIRE(IsWithinAbs(daeCurve(Real(0.25))[1], Real(0.5), Precision::IntegrationTolerance));
}

TEST_CASE("Reverse time - initial step probing stays inside the reverse interval",
	"[reverse-time][ode]")
{
	class DomainODE : public PolynomialODE {
	public:
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& f) const override {
			if (t < 0 || t > 1) throw ArgumentError("Derivative evaluated outside its domain");
			PolynomialODE::derivs(t, x, f);
		}
	} bounded;
	DormandPrince8Integrator integrator(bounded);
	auto solution = integrator.integrate(Vector<Real>{2}, 1, 0, Real(0.125));
	RequireReverse(solution);
	REQUIRE(IsWithinAbs(solution.getXValuesAtEnd()[0], Real(1), Precision::IntegrationTolerance));
	PolynomialODE unbounded;
	DormandPrince5Integrator negative(unbounded);
	auto belowZero = negative.integrate(Vector<Real>{2}, -1, -2, Real(0.125));
	RequireReverse(belowZero, -2);
	REQUIRE(IsWithinAbs(belowZero.getXValuesAtEnd()[0], Real(5), Precision::IntegrationTolerance));
}

TEST_CASE("Reverse time - exact step endpoint events trigger in either direction",
	"[reverse-time][events][ode][dae]")
{
	for (bool backward : {false, true}) {
		EventODE ode;
		ode.wanted = backward ? EventDirection::Decreasing : EventDirection::Increasing;
		DormandPrince5EventIntegrator integrator(ode);
		auto result = integrator.integrateWithEvents(ode, Vector<Real>{1}, backward ? Real(1) : Real(0),
			Real(0.75), Real(0.125), Precision::CoarseTolerance, Precision::EventTolerance, 1);
		REQUIRE(result.terminatedByEvent);
		REQUIRE(result.events.size() == 1);
		REQUIRE(IsWithinAbs(result.finalTime, Real(0.75), Precision::IntegrationTolerance));
		EventDAE dae;
		dae.wanted = ode.wanted;
		DAESolverConfig config;
		config.step_size = 1;
		auto daeResult = SolveDAEBackwardEulerWithEvents(dae, backward ? Real(1) : Real(0),
			Vector<Real>{1}, Vector<Real>{2}, Real(0.75), config);
		REQUIRE(daeResult.terminated_by_event);
		REQUIRE(daeResult.events.size() == 1);
		REQUIRE(IsWithinAbs(daeResult.final_time, Real(0.75), Precision::IntegrationTolerance));
	}
}

TEST_CASE("Reverse time - DAE drivers report stalled clocks without duplicate accepted points",
	"[reverse-time][dae]")
{
	PolynomialDAE system;
	const Real clock = Real(1) / Precision::Eps;
	DAESolverConfig config;
	config.step_size = Real(0.01);
	for (int method = 0; method < 6; ++method) {
		CAPTURE(method);
		auto result = [&]() {
			switch (method) {
			case 0: return SolveDAEBackwardEuler(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
			case 1: return SolveDAEBDF2(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
			case 2: return SolveDAEBDF4(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
			case 3: return SolveDAERadauIIA(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
			case 4: return SolveDAELinearlyImplicit(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
			default: return SolveDAEBackwardEulerAdaptive(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
			}
		}();
		REQUIRE(result.status == AlgorithmStatus::Stalled);
		REQUIRE(result.failure_reason == DAEFailureReason::StepSizeUnderflow);
		REQUIRE(result.solution.size() == 1);
		REQUIRE_FALSE(result.error_message.empty());
	}
}

}
