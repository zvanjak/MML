#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#include <mml/algorithms/ODESolvers/ODESolverEventDetection.h>
#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>
#include <mml/algorithms/DAESolvers.h>
#include <mml/algorithms/DAESolvers/DAEEventDetection.h>
#include <array>

namespace MML::Tests::TrajectoryResources {
	class ConstantODE : public IODESystemWithJacobian {
	public:
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
	};

	class ConstantODEEvents : public IODESystemWithEvents {
	public:
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
		int getNumEvents() const override { return 1; }
		Real eventFunction(int, Real, const Vector<Real>&) const override { return 1; }
	};

	class ConstantDAE : public IODESystemDAEWithEvents {
	public:
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
		int getNumEvents() const override { return 1; }
		Real eventFunction(int, Real, const Vector<Real>&, const Vector<Real>&) const override { return 1; }
	};

	class LastPointStopEvent : public ConstantODEEvents {
	public:
		Real eventFunction(int, Real time, const Vector<Real>&) const override { return time - Real(99999.5); }
		EventAction getEventAction(int) const override { return EventAction::Stop; }
	};

	TEST_CASE("Implicit ODE solvers pass the former trajectory column cap", "[trajectory-storage][ODE]")
	{
		const bool bdf2 = GENERATE(false, true);
		const Real direction = GENERATE(Real(1), Real(-1));
		ConstantODE system;
		StiffSolverConfig config;
		config.step_size = 1;
		config.max_steps = 100001;
		const auto result = bdf2 ? SolveBDF2(system, 0, Vector<Real>{3}, direction * 100001, config)
			: SolveBackwardEuler(system, 0, Vector<Real>{3}, direction * 100001, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.error_message.empty());
		REQUIRE(result.solution.size() == 100002);
		REQUIRE(result.total_steps == 100001);
		REQUIRE(result.solution.getNumStepsOK() == 100001);
		REQUIRE(result.solution.getNumStepsBad() == 0);
		REQUIRE(result.solution.getTValue(100001) == direction * 100001);
		REQUIRE(result.solution.getXValue(100001, 0) == Real(3));
	}

	TEST_CASE("ODE step budgets count intervals rather than stored initial points", "[trajectory-storage][ODE]")
	{
		ConstantODE system;
		StiffSolverConfig config;
		config.step_size = 1;
		config.max_steps = 99999;
		const auto complete = SolveBackwardEuler(system, 0, Vector<Real>{3}, 99999, config);
		REQUIRE(complete.status == AlgorithmStatus::Success);
		REQUIRE(complete.total_steps == 99999);
		REQUIRE(complete.solution.size() == 100000);
		const auto budget = SolveBackwardEuler(system, 0, Vector<Real>{3}, 100000, config);
		REQUIRE(budget.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(budget.solution.size() == 100000);
		config.max_steps = 100000;
		const auto beyond = SolveBackwardEuler(system, 0, Vector<Real>{3}, 100000, config);
		REQUIRE(beyond.status == AlgorithmStatus::Success);
		REQUIRE(beyond.solution.size() == 100001);
	}

	TEST_CASE("Legacy ODE APIs fail explicitly rather than silently dropping points", "[trajectory-storage][ODE]")
	{
		ConstantODE system;
		REQUIRE_THROWS_AS(SolveBackwardEuler(system, 0, Vector<Real>{3}, 100001, 1), ODESolverError);
		REQUIRE_THROWS_AS(SolveBDF2(system, 0, Vector<Real>{3}, 100001, 1), ODESolverError);
		ODESystemFixedStepSolver fixed(system, StepCalculators::RK4_Basic);
		const auto longSolution = fixed.integrate(Vector<Real>{3}, 0, 100000, 100000);
		REQUIRE(longSolution.size() == 100001);
		const int overBudgetSteps = static_cast<int>(Matrix<Real>::MaxElementCount());
		REQUIRE_THROWS_AS(fixed.integrate(Vector<Real>{3}, 0, 1, overBudgetSteps), ResourceLimitError);
	}

	TEST_CASE("Every DAE driver passes the former trajectory column cap", "[trajectory-storage][DAE]")
	{
		using Solver = DAESolverResult (*)(IODESystemDAEWithJacobian&, Real, const Vector<Real>&,
			const Vector<Real>&, Real, const DAESolverConfig&);
		const std::array<Solver, 7> solvers{SolveDAEBackwardEuler, SolveDAEBDF2, SolveDAEBDF4,
			SolveDAERadauIIA, SolveDAELinearlyImplicit, SolveDAERODAS4, SolveDAEBackwardEulerAdaptive};
		const int method = GENERATE(0, 1, 2, 3, 4, 5, 6);
		INFO("DAE method " << method);
		ConstantDAE system;
		DAESolverConfig config;
		config.step_size = config.max_step_size = 1;
		config.max_steps = 100001;
		const auto result = solvers[method](system, 0, Vector<Real>{3}, Vector<Real>{3}, 100001, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.failure_reason == DAEFailureReason::None);
		REQUIRE(result.error_message.empty());
		REQUIRE(result.solution.size() == 100002);
		REQUIRE(result.accepted_steps == 100001);
		REQUIRE(result.total_steps == 100001);
		REQUIRE(result.solution.getNumStepsOK() == 100001);
		REQUIRE(result.rejected_steps == 0);
		REQUIRE(result.solution.getTValue(100001) == Real(100001));
		REQUIRE(result.solution.getXValue(100001, 0) == Real(3));
		REQUIRE(result.solution.getYValue(100001, 0) == Real(3));
	}

	TEST_CASE("Oversized ODE and DAE point reservations retain resource-limit diagnostics", "[trajectory-storage]")
	{
		StiffSolverResult ode(0, 1, 1, 0);
		ode.solution.fillValues(0, 0, Vector<Real>{3});
		REQUIRE_FALSE(detail::PrepareTrajectoryPoint(ode.solution, ode.solution.max_size(), ode.status, ode.error_message));
		REQUIRE(ode.status == AlgorithmStatus::ResourceLimitExceeded);
		REQUIRE_FALSE(ode.error_message.empty());
		REQUIRE(ode.solution.size() == 1);
		REQUIRE(ode.solution.getXValue(0, 0) == Real(3));
		REQUIRE(ode.total_steps == 0);
		DAESolverResult dae(0, 1, 1, 1, 0);
		dae.solution.fillValues(0, 0, Vector<Real>{3}, Vector<Real>{4});
		REQUIRE_FALSE(PrepareDAETrajectoryPoint(dae, dae.solution.max_size()));
		REQUIRE(dae.status == AlgorithmStatus::ResourceLimitExceeded);
		REQUIRE(dae.failure_reason == DAEFailureReason::ResourceLimitExceeded);
		REQUIRE_FALSE(dae.error_message.empty());
		REQUIRE(dae.solution.size() == 1);
		REQUIRE(dae.solution.getYValue(0, 0) == Real(4));
		REQUIRE(dae.accepted_steps == 0);
	}

	TEST_CASE("Fixed DAE budgets permit the requested number of accepted intervals", "[trajectory-storage][DAE]")
	{
		ConstantDAE system;
		DAESolverConfig config;
		config.max_steps = 1;
		config.step_size = 1;
		const auto result = SolveDAERadauIIA(system, 0, Vector<Real>{3}, Vector<Real>{3}, 1, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.accepted_steps == 1);
		REQUIRE(result.solution.size() == 2);
	}

	TEST_CASE("ROS2 step budgets preserve trajectories beyond the former column cap", "[trajectory-storage][ODE]")
	{
		ConstantODE system;
		ODEStiffConfig config;
		config.solver.step_size = 1;
		config.solver.max_steps = 100001;
		const auto result = SolveROS2Detailed(system, Vector<Real>{3}, 0, 100000000, config);
		REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(result.solution.has_value());
		REQUIRE(result.solution->size() == 100002);
		REQUIRE(result.total_steps == 100001);
		REQUIRE(result.solution->getNumStepsOK() == 100001);
		REQUIRE(result.solution->getXValue(100001, 0) == Real(3));
		config.solver.max_steps = 2;
		const auto budget = SolveROS2Detailed(system, Vector<Real>{3}, 0, 100000000, config);
		REQUIRE(budget.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(budget.solution.has_value());
		REQUIRE(budget.total_steps == 2);
		REQUIRE(budget.solution->size() == 3);
	}

	TEST_CASE("Dense output retains points beyond the former trajectory cap", "[trajectory-storage][ODE]")
	{
		ConstantODE system;
		ODEAdaptiveConfig config;
		config.integrator.output_interval = 1;
		config.integrator.initial_step_size = 1;
		const auto result = ODEAdaptiveIntegrateDetailed(system, Vector<Real>{3}, 0, 100001, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.has_value());
		REQUIRE(result.solution->size() == 100002);
		REQUIRE(result.solution->getTValue(100001) == Real(100001));
		REQUIRE(result.solution->getXValue(100001, 0) == Real(3));
		DormandPrince5Integrator legacy(system);
		const auto legacySolution = legacy.integrate(Vector<Real>{3}, 0, 100001, 1);
		REQUIRE(legacySolution.size() == 100002);
	}

	TEST_CASE("Successful adaptive grids preserve compact legacy output storage", "[trajectory-storage][ODE]")
	{
		ConstantODE system;
		DormandPrince5Integrator integrator(system);
		const auto uniform = integrator.integrate(Vector<Real>{3}, 0, 10, 1);
		REQUIRE(uniform.size() == 11);
		REQUIRE(uniform.getTotalSavedSteps() == uniform.size());
		const auto selected = integrator.integrateAt(Vector<Real>{3}, Vector<Real>{0, 1, 2, 3, 4, 5});
		REQUIRE(selected.size() == 6);
		REQUIRE(selected.getTotalSavedSteps() == selected.size());
	}

	TEST_CASE("ODE and DAE event drivers pass the former trajectory column cap", "[trajectory-storage][ODE][DAE]")
	{
		ConstantODEEvents ode;
		DormandPrince5EventIntegrator integrator(ode);
		const auto result = integrator.integrateWithEvents(ode, Vector<Real>{3}, 0, 100001, 1);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.size() == 100002);
		REQUIRE(result.finalTime == result.solution.getTValue(100001));
		REQUIRE(result.finalState[0] == Real(3));
		REQUIRE_FALSE(result.terminatedByEvent);
		ConstantDAE dae;
		DAESolverConfig config;
		config.step_size = 1;
		config.max_steps = 100001;
		const auto daeResult = SolveDAEBackwardEulerWithEvents(dae, 0, Vector<Real>{3}, Vector<Real>{3}, 100001, config);
		REQUIRE(daeResult.integration.status == AlgorithmStatus::Success);
		REQUIRE(daeResult.integration.failure_reason == DAEFailureReason::None);
		REQUIRE(daeResult.integration.solution.size() == 100002);
		REQUIRE(daeResult.integration.accepted_steps == 100001);
	}

	TEST_CASE("A terminal event can append a point beyond the former trajectory cap", "[trajectory-storage][ODE]")
	{
		LastPointStopEvent system;
		DormandPrince5EventIntegrator integrator(system);
		const auto result = integrator.integrateWithEvents(system, Vector<Real>{3}, 0, 100001, 1);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.size() == 100001);
		REQUIRE(result.finalTime > Real(99999));
		REQUIRE(result.finalTime < Real(100000));
		REQUIRE(result.finalTime == result.solution.getTValue(100000));
		REQUIRE(result.terminatedByEvent);
		REQUIRE(result.events.size() == 1);
	}
}
