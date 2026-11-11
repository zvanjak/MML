#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/ODESolvers/BVPShootingMethod.h>
#include <mml/algorithms/ODESolvers/ODESolverStiff.h>
#include <mml/algorithms/DAESolvers.h>
#endif
#include <catch2/catch_test_macros.hpp>

using namespace MML;

namespace {
	class BoundaryODE : public IODESystem {
	public:
		bool stationary = false;
		bool quadratic = false;
		Real scale = 1;
		int getDim() const override { return 2; }
		void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override {
			f[0] = stationary ? Real(0) : scale * (quadratic ? x[1] * x[1] : x[1]);
			f[1] = 0;
		}
	};

	class BoundaryDAE : public IODESystemDAEWithJacobian {
	public:
		bool quadratic = false;
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override { f[0] = -x[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			g[0] = (quadratic ? y[0] * y[0] : y[0]) - x[0];
		}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>& y, Matrix<Real>& J) const override {
			J(0, 0) = quadratic ? 2 * y[0] : Real(1);
		}
	};

	class BoundaryStiffODE : public IODESystemWithJacobian {
	public:
		bool nonfiniteAfterUpdate = false;
		int getDim() const override { return 1; }
		void derivs(Real, const Vector<Real>& x, Vector<Real>& f) const override {
			f[0] = nonfiniteAfterUpdate && x[0] != 1 ? std::numeric_limits<Real>::infinity() : -x[0];
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
	};

	ShootingConfig LinearShootingConfig(Real guess = 0) {
		ShootingConfig config;
		config.unknownIndices = {1};
		config.unknownGuesses = {guess};
		config.targetIndices = {0};
		config.targetValues = {1};
		return config;
	}
}

TEST_CASE("Boundary Newton - scalar shooting failure status and refinement budget are honest", "[boundary-newton][bvp_shooting]") {
	BoundaryODE system;
	SECTION("Unsatisfiable boundary is not Success") {
		system.stationary = true;
		BVPShootingSolver solver(system, 1, Precision::BVPTolerance, 1);
		const auto result = solver.solvePositionBVP(0, 1, 0, 1, -1, 2);
		CHECK_FALSE(result.converged);
		CHECK(result.failureReason == BVPFailureReason::BracketingFailed);
		CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
		CHECK_FALSE(result.failureDetails.empty());
	}
	SECTION("A nonlinear refinement cannot silently exceed its budget") {
		system.quadratic = true;
		BVPShootingSolver limited(system, 1, Precision::BVPTolerance, 1);
		const auto result = limited.solvePositionBVP(0, 1, 0, 2, 0, 2);
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(result.failureReason == BVPFailureReason::MaxIterationsExceeded);
		CHECK(result.iterations == 1);
		CHECK(result.residual > Precision::BVPTolerance);
		BVPShootingSolver sufficient(system, 1, Precision::BVPTolerance, 50);
		const auto success = sufficient.solvePositionBVP(0, 1, 0, 2, 0, 2);
		CHECK(success.success());
		CHECK(success.status == AlgorithmStatus::Success);
		CHECK(success.residual <= Precision::BVPTolerance);
	}
	SECTION("The final Brent update can already satisfy the boundary") {
		BVPShootingSolver solver(system, 1, Precision::BVPTolerance, 1);
		const auto result = solver.solvePositionBVP(0, 1, 0, 1, -1, 2);
		CHECK(result.success());
		CHECK(result.status == AlgorithmStatus::Success);
		CHECK(result.iterations == 1);
		CHECK(result.residual <= Precision::BVPTolerance);
	}
	SECTION("The last allowed bracket expansion is checked") {
		BVPShootingSolver solver(system, 1, Precision::BVPTolerance, 1);
		const auto result = solver.solvePositionBVP(0, 1, 0, 1, 2, 3);
		CHECK(result.success());
		CHECK(result.status == AlgorithmStatus::Success);
		CHECK(result.residual <= Precision::BVPTolerance);
	}
	SECTION("An initial bracket endpoint can be the solution") {
		BVPShootingSolver solver(system, 1, Precision::BVPTolerance, 1);
		const auto result = solver.solvePositionBVP(0, 1, 0, 1, 0, 1);
		CHECK(result.success());
		CHECK(result.status == AlgorithmStatus::Success);
		CHECK(result.iterations == 0);
		CHECK(result.residual == 0);
	}
	SECTION("Invalid budget is not the root finder's default budget") {
		BVPShootingSolver solver(system, 1, Precision::BVPTolerance, 0);
		const auto result = solver.solvePositionBVP(0, 1, 0, 1, -1, 2);
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK(result.failureReason == BVPFailureReason::InvalidConfiguration);
	}
	SECTION("Parameter convergence does not invent boundary convergence") {
		system.scale = 4 / Precision::BVPTolerance;
		BVPShootingSolver solver(system, 1, Precision::BVPTolerance, 50);
		const auto result = solver.solvePositionBVP(0, 1, 0, 1, 0, 1);
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::ToleranceUnachievable);
		CHECK(result.failureReason == BVPFailureReason::RootFindingFailed);
		CHECK(result.residual > Precision::BVPTolerance);
	}
}

TEST_CASE("Boundary Newton - ND shooting checks the final parameter residual", "[boundary-newton][bvp_shooting]") {
	BoundaryODE system;
	BVPShootingSolverND solver(system, 1, Precision::BVPTolerance, 1, Real(1));
	const auto result = solver.solve(0, 1, Vector<Real>{0, 0}, LinearShootingConfig());
	CHECK(result.success());
	CHECK(result.status == AlgorithmStatus::Success);
	CHECK(result.failureReason == BVPFailureReason::None);
	CHECK(result.iterations == 1);
	CHECK(result.residualNorm <= Precision::BVPTolerance);
	CHECK(result.residuals[0] == result.solution.getXValuesAtEnd()[0] - 1);

	BVPShootingSolverND noUpdates(system, 1, Precision::BVPTolerance, 0, Real(1));
	const auto initialRoot = noUpdates.solve(0, 1, Vector<Real>{0, 0}, LinearShootingConfig(1));
	CHECK(initialRoot.success());
	CHECK(initialRoot.iterations == 0);
	const auto exhausted = noUpdates.solve(0, 1, Vector<Real>{0, 0}, LinearShootingConfig());
	CHECK_FALSE(exhausted.converged);
	CHECK(exhausted.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(exhausted.failureReason == BVPFailureReason::MaxIterationsExceeded);
	CHECK(exhausted.residualNorm == 1);

	system.stationary = true;
	const auto singular = solver.solve(0, 1, Vector<Real>{0, 0}, LinearShootingConfig());
	CHECK_FALSE(singular.converged);
	CHECK(singular.status == AlgorithmStatus::SingularMatrix);
	CHECK(singular.failureReason == BVPFailureReason::SingularJacobian);
	const auto invalid = solver.solve(0, 1, Vector<Real>{0, 0}, ShootingConfig{});
	CHECK(invalid.status == AlgorithmStatus::InvalidInput);

	system.stationary = false;
	system.quadratic = true;
	const auto nonlinear = solver.solve(0, 1, Vector<Real>{0, 0}, LinearShootingConfig(2));
	CHECK_FALSE(nonlinear.converged);
	CHECK(nonlinear.iterations == 1);
	CHECK(nonlinear.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(nonlinear.residualNorm > Precision::BVPTolerance);
	CHECK(nonlinear.residualNorm == std::abs(nonlinear.residuals[0]));
}

TEST_CASE("Boundary Newton - consistent IC reports the residual of its returned state", "[boundary-newton][dae]") {
	BoundaryDAE system;
	Vector<Real> y{2};
	const Vector<Real> x{1};
	const auto result = ComputeConsistentICDetailed(system, 0, x, y, 1);
	CHECK(result.converged);
	CHECK(result.failure_reason == DAEFailureReason::None);
	CHECK(result.iterations == 1);
	CHECK(result.initial_residual_norm == 1);
	CHECK(result.final_residual_norm == 0);
	CHECK(y[0] == 1);
	y[0] = 2;
	CHECK(ComputeConsistentIC(system, 0, x, y, 1));
	y[0] = 1;
	const auto initialRoot = ComputeConsistentICDetailed(system, 0, x, y, 0);
	CHECK(initialRoot.converged);
	CHECK(initialRoot.iterations == 0);
	y[0] = 2;
	const auto noUpdates = ComputeConsistentICDetailed(system, 0, x, y, 0);
	CHECK_FALSE(noUpdates.converged);
	CHECK(noUpdates.final_residual_norm == 1);
	CHECK(y[0] == 2);

	system.quadratic = true;
	const auto nonlinear = ComputeConsistentICDetailed(system, 0, Vector<Real>{2}, y, 1);
	CHECK_FALSE(nonlinear.converged);
	CHECK(nonlinear.failure_reason == DAEFailureReason::NewtonFailure);
	CHECK(nonlinear.iterations == 1);
	CHECK(y[0] == Real(1.5));
	CHECK(nonlinear.final_residual_norm == Real(0.25));
}

TEST_CASE("Boundary Newton - shared DAE final check does not grant another correction", "[boundary-newton][dae]") {
	DAESolverConfig config;
	config.max_newton_iter = 1;
	Vector<Real> x{0}, y{0};
	int residualCalls = 0;
	const auto residual = [&](const Vector<Real>& xx, const Vector<Real>& yy, Vector<Real>& r) {
		++residualCalls;
		r[0] = xx[0] - 1;
		r[1] = yy[0] - xx[0];
	};
	const auto jacobian = [](const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) {
		J(0, 0) = 1; J(0, 1) = 0;
		J(1, 0) = -1; J(1, 1) = 1;
	};
	const auto result = SolveDAENewton(x, y, 1, 1, config, residual, jacobian);
	CHECK(result.converged);
	CHECK(result.failure_reason == DAEFailureReason::None);
	CHECK(result.residual_norm == 0);
	CHECK(result.linear_solves == 1);
	CHECK(result.iterations == 1);
	CHECK(residualCalls == 2);
	CHECK(x[0] == 1);
	CHECK(y[0] == 1);
	config.max_newton_iter = 0;
	const auto initialRoot = SolveDAENewton(x, y, 1, 1, config, residual, jacobian);
	CHECK(initialRoot.converged);
	CHECK(initialRoot.iterations == 0);
	CHECK(initialRoot.linear_solves == 0);
	x[0] = 0;
	const auto noUpdates = SolveDAENewton(x, y, 1, 1, config, residual, jacobian);
	CHECK_FALSE(noUpdates.converged);
	CHECK(noUpdates.residual_norm == std::sqrt(Real(2)));
	CHECK(noUpdates.linear_solves == 0);
	CHECK(x[0] == 0);
	config.max_newton_iter = 1;

	x[0] = 2;
	Vector<Real> empty;
	const auto nonlinearResidual = [](const Vector<Real>& xx, const Vector<Real>&, Vector<Real>& r) {
		r[0] = xx[0] * xx[0] - 2;
	};
	const auto nonlinearJacobian = [](const Vector<Real>& xx, const Vector<Real>&, Matrix<Real>& J) { J(0, 0) = 2 * xx[0]; };
	const auto exhausted = SolveDAENewton(x, empty, 1, 0, config, nonlinearResidual, nonlinearJacobian);
	CHECK_FALSE(exhausted.converged);
	CHECK(exhausted.failure_reason == DAEFailureReason::NewtonFailure);
	CHECK(exhausted.residual_norm == Real(0.25));
	CHECK(exhausted.linear_solves == 1);
	CHECK(x[0] == Real(1.5));
}

TEST_CASE("Boundary Newton - final nonfinite residual remains an explicit DAE failure", "[boundary-newton][dae]") {
	DAESolverConfig config;
	config.max_newton_iter = 1;
	Vector<Real> x{0}, y;
	const auto residual = [](const Vector<Real>& xx, const Vector<Real>&, Vector<Real>& r) {
		r[0] = xx[0] == 0 ? Real(-1) : std::numeric_limits<Real>::infinity();
	};
	const auto jacobian = [](const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) { J(0, 0) = 1; };
	const auto result = SolveDAENewton(x, y, 1, 0, config, residual, jacobian);
	CHECK_FALSE(result.converged);
	CHECK(result.failure_reason == DAEFailureReason::NonFiniteState);
	CHECK_FALSE(HasEstimate(result.residual_norm));
	CHECK_FALSE(result.message.empty());
}

TEST_CASE("Boundary Newton - nonfinite final ODE residual does not enter the trajectory", "[boundary-newton][StiffSolvers]") {
	BoundaryStiffODE system;
	system.nonfiniteAfterUpdate = true;
	StiffSolverConfig config;
	config.step_size = Real(0.25);
	config.max_newton_iter = 1;
	const auto result = SolveBackwardEuler(system, 0, Vector<Real>{1}, Real(0.25), config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.linear_solves == 1);
	CHECK(result.function_evaluations == 2);
	CHECK(result.solution.size() == 1);
	CHECK(result.solution.getXValuesAtEnd()[0] == 1);
}

TEST_CASE("Boundary Newton - linear implicit ODE and DAE steps finish with one correction", "[boundary-newton][StiffSolvers][dae]") {
	BoundaryStiffODE ode;
	BoundaryDAE dae;
	for (bool bdf2 : {false, true}) {
		for (Real direction : {Real(-1), Real(1)}) {
			CAPTURE(bdf2, direction);
			StiffSolverConfig odeConfig;
			odeConfig.step_size = Real(0.25);
			odeConfig.min_step_size = odeConfig.step_size;
			odeConfig.max_newton_iter = 1;
			odeConfig.max_steps = 2;
			DAESolverConfig daeConfig;
			daeConfig.step_size = odeConfig.step_size;
			daeConfig.min_step_size = odeConfig.min_step_size;
			daeConfig.max_newton_iter = 1;
			daeConfig.max_steps = 2;
			const auto odeResult = bdf2 ? SolveBDF2(ode, 0, Vector<Real>{1}, direction / 2, odeConfig)
				: SolveBackwardEuler(ode, 0, Vector<Real>{1}, direction / 2, odeConfig);
			const auto daeResult = bdf2 ? SolveDAEBDF2(dae, 0, Vector<Real>{1}, Vector<Real>{1}, direction / 2, daeConfig)
				: SolveDAEBackwardEuler(dae, 0, Vector<Real>{1}, Vector<Real>{1}, direction / 2, daeConfig);
			CHECK(odeResult.status == AlgorithmStatus::Success);
			CHECK(daeResult.status == AlgorithmStatus::Success);
			CHECK(odeResult.total_steps == 2);
			CHECK(daeResult.accepted_steps == 2);
			CHECK(odeResult.solution.getNumStepsBad() == 0);
			CHECK(daeResult.rejected_steps == 0);
			CHECK(odeResult.linear_solves == 2);
			CHECK(daeResult.linear_solves == 2);
			CHECK(daeResult.final_residual_norm < daeConfig.newton_tol);
			CHECK(IsWithinAbs(odeResult.solution.getXValuesAtEnd()[0],
				daeResult.solution.getXValuesAtEnd()[0], Precision::EqualityTolerance));
		}
	}
}
