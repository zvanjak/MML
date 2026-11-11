#include "../../TestPrecision.h"
#include <mml/algorithms/DAESolvers.h>
#include <test_beds/stiff_solver_benchmarks.h>
#include <catch2/catch_test_macros.hpp>
#include <array>
#include <limits>
#include <stdexcept>

namespace MML::Tests::RODAS4 {

	class LinearDAE : public IODESystemDAEWithJacobian {
	public:
		Real fx = -3;
		mutable int jacobians = 0;
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = fx * x[0] + y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - 2 * x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = fx; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -2; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
		void allJacobians(Real t, const Vector<Real>& x, const Vector<Real>& y, Matrix<Real>& a, Matrix<Real>& b, Matrix<Real>& c,
											Matrix<Real>& d) const override {
			++jacobians;
			IODESystemDAEWithJacobian::allJacobians(t, x, y, a, b, c, d);
		}
	};

	class NonlinearDAE : public LinearDAE {
	public:
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0] * x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
		void jacobian_gx(Real, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -2 * x[0]; }
	};

	class TimedDAE : public LinearDAE {
	public:
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			if (t < 0 || t > 1)
				throw std::domain_error("Time probe outside [0,1]");
			f[0] = -3 * x[0] + y[0] + std::cos(t) + std::sin(t) - t * t;
		}
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			g[0] = y[0] - 2 * x[0] - t * t;
		}
	};

	class AlgebraicDAE : public LinearDAE {
	public:
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
		void algConstraints(Real t, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] * y[0] - 1 - t; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>& y, Matrix<Real>& J) const override { J(0, 0) = 2 * y[0]; }
	};

	// Robertson with two differential species and mass conservation as the algebraic equation.
	class RobertsonDAE : public IODESystemDAEWithJacobian {
	public:
		int getDiffDim() const override { return 2; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			f[0] = -Real(0.04) * x[0] + Real(1e4) * x[1] * y[0];
			f[1] = -f[0] - Real(3e7) * x[1] * x[1];
		}
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			g[0] = x[0] + x[1] + y[0] - 1;
		}
		void jacobian_fx(Real, const Vector<Real>& x, const Vector<Real>& y, Matrix<Real>& J) const override {
			J(0, 0) = -Real(0.04);
			J(0, 1) = Real(1e4) * y[0];
			J(1, 0) = Real(0.04);
			J(1, 1) = -Real(1e4) * y[0] - Real(6e7) * x[1];
		}
		void jacobian_fy(Real, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& J) const override {
			J(0, 0) = Real(1e4) * x[1];
			J(1, 0) = -J(0, 0);
		}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override {
			J(0, 0) = 1;
			J(0, 1) = 1;
		}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
	};

	class ODEAsDAE : public IODESystemDAEWithJacobian {
		const IODESystemWithJacobian& _system;

	public:
		explicit ODEAsDAE(const IODESystemWithJacobian& system)
				: _system(system) {}
		int getDiffDim() const override { return _system.getDim(); }
		int getAlgDim() const override { return 0; }
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override { _system.derivs(t, x, f); }
		void algConstraints(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&) const override {}
		void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& J) const override { _system.jacobian(t, x, J); }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
	};

	DAESolverConfig FixedConfig(Real h) {
		DAESolverConfig config;
		config.step_size = config.min_step_size = config.max_step_size = h;
		config.abs_tolerance = 1;
		config.rel_tolerance = 0;
		config.constraint_tol = 1;
		return config;
	}

	void RequireSuccess(const DAESolverResult& result, Real end) {
		INFO(result.error_message);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.failure_reason == DAEFailureReason::None);
		REQUIRE(result.solution.getTValue(result.solution.size() - 1) == end);
		REQUIRE(result.total_steps == result.accepted_steps);
		REQUIRE(result.newton_iterations == 0);
		REQUIRE(result.linear_solves == 6 * result.lu_factorizations);
		REQUIRE(result.lu_factorizations == result.accepted_steps + result.rejected_steps);
		REQUIRE(result.last_error_norm <= 1);
	}

	TEST_CASE("RODAS4 - fourth order on coupled linear and nonlinear index-1 DAEs", "[dae][RODAS4]") {
		LinearDAE linear;
		NonlinearDAE nonlinear;
		for (bool curved : {false, true}) {
			CAPTURE(curved);
			std::array<Real, 3> errors;
			for (int level = 0; level < 3; ++level) {
				Real h = Real(0.125) / (1 << level);
				auto config = FixedConfig(h);
				auto& system = curved ? static_cast<IODESystemDAEWithJacobian&>(nonlinear) : linear;
				auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{curved ? Real(1) : Real(2)}, 1, config);
				RequireSuccess(result, 1);
				REQUIRE(result.rejected_steps == 0);
				errors[level] = std::abs(result.solution.getXValuesAtEnd()[0] - (curved ? Real(0.5) : std::exp(Real(-1))));
			}
			CAPTURE(errors);
			REQUIRE(errors[0] / errors[1] > 12);
			REQUIRE(errors[0] / errors[1] < 20);
			REQUIRE(errors[1] / errors[2] > 12);
			REQUIRE(errors[1] / errors[2] < 20);
		}
	}

	TEST_CASE("RODAS4 - sixth increment is a fourth-power embedded error estimate", "[dae][RODAS4]") {
		LinearDAE system;
		std::array<Real, 2> errors;
		for (int level = 0; level < 2; ++level) {
			Real h = Real(0.125) / (1 << level);
			auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, h, FixedConfig(h));
			RequireSuccess(result, h);
			errors[level] = result.last_error_norm;
			REQUIRE(errors[level] > 0);
		}
		REQUIRE(errors[0] / errors[1] > 12);
		REQUIRE(errors[0] / errors[1] < 20);
	}

	TEST_CASE("RODAS4 - nonautonomous differential and algebraic time terms round trip", "[dae][RODAS4][reverse-time]") {
		TimedDAE system;
		DAESolverConfig config;
		config.step_size = Real(0.125);
		config.abs_tolerance = config.rel_tolerance = Precision::EqualityTolerance;
		config.min_step_size = Precision::ZeroTolerance;
		auto forward = SolveDAERODAS4(system, 0, Vector<Real>{0}, Vector<Real>{0}, 1, config);
		RequireSuccess(forward, 1);
		REQUIRE(IsWithinAbs(forward.solution.getXValuesAtEnd()[0], std::sin(Real(1)), config.abs_tolerance * 10));
		auto reverse = SolveDAERODAS4(system, 1, forward.solution.getXValuesAtEnd(), forward.solution.getYValuesAtEnd(), 0, config);
		RequireSuccess(reverse, 0);
		REQUIRE(IsWithinAbs(reverse.solution.getXValuesAtEnd()[0], 0, config.abs_tolerance * 20));
		for (int i = 1; i < reverse.solution.size(); ++i)
			REQUIRE(reverse.solution.getTValue(i) < reverse.solution.getTValue(i - 1));
		REQUIRE(system.jacobians == forward.jacobian_evaluations + reverse.jacobian_evaluations);
	}

	TEST_CASE("RODAS4 - tolerance refinement rejects trials and preserves nonlinear constraints", "[dae][RODAS4]") {
		NonlinearDAE system;
		std::array<Real, 2> errors;
		int coarseSteps = 0;
		for (int level = 0; level < 2; ++level) {
			DAESolverConfig config;
			config.step_size = 1;
			config.abs_tolerance = config.rel_tolerance = level == 0 ? Precision::LooseTolerance : Precision::EqualityTolerance;
			config.min_step_size = Precision::ZeroTolerance;
			config.constraint_tol = config.abs_tolerance;
			auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
			RequireSuccess(result, 1);
			REQUIRE(result.rejected_steps > 0);
			REQUIRE(result.max_error_norm > 1);
			REQUIRE(result.solution.getNumStepsBad() == result.rejected_steps);
			REQUIRE(result.solution.getNumStepsOK() == result.accepted_steps);
			for (int i = 0; i < result.solution.size(); ++i) {
				Real x = result.solution.getXValue(i, 0), y = result.solution.getYValue(i, 0);
				REQUIRE(std::abs(y - x * x) <= config.constraint_tol);
			}
			errors[level] = std::abs(result.solution.getXValuesAtEnd()[0] - Real(0.5));
			if (level == 0)
				coarseSteps = result.accepted_steps;
			else
				REQUIRE(result.accepted_steps > coarseSteps);
		}
		REQUIRE(errors[1] < errors[0] / 100);
	}

	TEST_CASE("RODAS4 - algebraic error participates when differential error vanishes", "[dae][RODAS4]") {
		AlgebraicDAE system;
		DAESolverConfig config;
		config.step_size = 1;
		config.abs_tolerance = config.rel_tolerance = Precision::EqualityTolerance;
		config.constraint_tol = Precision::EqualityTolerance;
		config.min_step_size = Precision::ZeroTolerance;
		auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		RequireSuccess(result, 1);
		REQUIRE(result.rejected_steps > 0);
		REQUIRE(result.solution.getXValuesAtEnd()[0] == 1);
		REQUIRE(IsWithinAbs(result.solution.getYValuesAtEnd()[0], std::sqrt(Real(2)), Precision::EqualityTolerance * 10));
		REQUIRE(result.max_constraint_violation <= config.constraint_tol);
	}

	TEST_CASE("RODAS4 - L-stability damps unresolved stiff decay", "[dae][RODAS4][StiffSolvers]") {
		LinearDAE system;
		Real previous = 1;
		for (Real rate : {Real(1000), Real(1000000)}) {
			system.fx = -rate - 2;
			auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, FixedConfig(1));
			RequireSuccess(result, 1);
			Real value = std::abs(result.solution.getXValuesAtEnd()[0]);
			REQUIRE(value < Real(0.01));
			REQUIRE(value < previous / 100);
			previous = value;
		}
	}

	TEST_CASE("RODAS4 - Robertson mass-constrained DAE reaches independent short and published full references",
						"[dae][RODAS4][StiffSolvers][stiff-benchmark]") {
		RobertsonDAE system;
		auto problem = TestBeds::StiffBenchmarkProblems().front();
		for (bool full : {false, true}) {
			const Real end = full ? problem.endTime : problem.smokeEndTime;
			const auto& reference = full ? problem.reference : problem.smokeReference;
			DAESolverConfig config;
			config.step_size = Real(1e-6);
			config.min_step_size = Precision::ZeroTolerance;
			config.max_step_size = end;
			config.abs_tolerance = Precision::ZeroTolerance / 100;
			config.rel_tolerance = Precision::EqualityTolerance;
			config.constraint_tol = Precision::EqualityTolerance;
			auto result = SolveDAERODAS4(system, 0, Vector<Real>{1, 0}, Vector<Real>{0}, end, config);
			CAPTURE(full, result.accepted_steps, result.rejected_steps);
			RequireSuccess(result, end);
			auto x = result.solution.getXValuesAtEnd();
			auto y = result.solution.getYValuesAtEnd();
			for (int i = 0; i < 3; ++i) {
				Real value = i < 2 ? x[i] : y[0];
				CAPTURE(i, value, reference[i]);
				REQUIRE(IsWithinRel(value, reference[i], Precision::CoarseTolerance));
			}
			CAPTURE(y[0], reference[2]);
			REQUIRE(IsWithinAbs(y[0], reference[2], Precision::EqualityTolerance));
			for (int i = 0; i < result.solution.size(); ++i)
				REQUIRE(std::abs(result.solution.getXValue(i, 0) + result.solution.getXValue(i, 1) + result.solution.getYValue(i, 0) - 1) <=
								config.constraint_tol);
		}
	}

	TEST_CASE("RODAS4 - zero algebraic dimensions recover classical stiff ODE references", "[dae][RODAS4][StiffSolvers]") {
		for (const auto& problem : TestBeds::StiffBenchmarkProblems()) {
			auto original = problem.createSystem();
			ODEAsDAE system(*original);
			DAESolverConfig config;
			config.step_size = Real(1e-6);
			config.min_step_size = Precision::ZeroTolerance;
			config.abs_tolerance = Precision::ZeroTolerance;
			config.rel_tolerance = Precision::ModerateTolerance;
			auto result = SolveDAERODAS4(system, 0, problem.initial, Vector<Real>{}, problem.smokeEndTime, config);
			CAPTURE(problem.name);
			RequireSuccess(result, problem.smokeEndTime);
			for (int i = 0; i < problem.initial.size(); ++i)
				REQUIRE(IsWithinAbs(result.solution.getXValuesAtEnd()[i], problem.smokeReference[i], Precision::LooseTolerance));
		}
	}

	TEST_CASE("RODAS4 - stage callback errors propagate and nonfinite stages retain accepted state", "[dae][RODAS4]") {
		class BadStage : public LinearDAE {
		public:
			bool throwing = false;
			void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
				if (t > Real(0.25)) {
					if (throwing)
						throw SingularMatrixError("User callback failure, not a stage factorization");
					f[0] = std::numeric_limits<Real>::quiet_NaN();
				} else
					LinearDAE::diffEqs(t, x, y, f);
			}
		} system;
		auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, FixedConfig(1));
		REQUIRE(result.status == AlgorithmStatus::NumericalInstability);
		REQUIRE(result.failure_reason == DAEFailureReason::NonFiniteState);
		REQUIRE(result.accepted_steps == 0);
		REQUIRE(result.rejected_steps == 1);
		REQUIRE(result.solution.size() == 1);
		REQUIRE(result.linear_solves == 1);
		system.throwing = true;
		REQUIRE_THROWS_AS(SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, FixedConfig(1)), SingularMatrixError);
	}

	TEST_CASE("RODAS4 - exhausted budgets keep accepted partial trajectories", "[dae][RODAS4]") {
		LinearDAE system;
		auto config = FixedConfig(Real(0.125));
		config.max_steps = 2;
		auto partial = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
		REQUIRE(partial.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(partial.failure_reason == DAEFailureReason::MaxStepsExceeded);
		REQUIRE(partial.accepted_steps == 2);
		REQUIRE(partial.solution.size() == 3);
		REQUIRE(partial.solution.getTValue(2) == Real(0.25));
		REQUIRE_FALSE(partial.error_message.empty());
		config.step_size = config.max_step_size = 1;
		config.min_step_size = Precision::ZeroTolerance;
		config.max_steps = 1;
		config.abs_tolerance = Precision::ZeroTolerance;
		auto rejected = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
		REQUIRE(rejected.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(rejected.accepted_steps == 0);
		REQUIRE(rejected.rejected_steps == 1);
		REQUIRE(rejected.solution.size() == 1);
		REQUIRE(rejected.solution.getTValue(0) == 0);
	}

	TEST_CASE("RODAS4 - minimum steps never force acceptance of failed error or constraints", "[dae][RODAS4]") {
		LinearDAE linear;
		AlgebraicDAE algebraic;
		for (bool constraint : {false, true}) {
			auto config = FixedConfig(1);
			config.abs_tolerance = constraint ? Real(1) : Precision::ZeroTolerance;
			config.constraint_tol = constraint ? Precision::ZeroTolerance : Real(1);
			auto& system = constraint ? static_cast<IODESystemDAEWithJacobian&>(algebraic) : linear;
			auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{constraint ? Real(1) : Real(2)}, 1, config);
			REQUIRE(result.status == AlgorithmStatus::NumericalInstability);
			REQUIRE(result.failure_reason == (constraint ? DAEFailureReason::ConstraintViolation : DAEFailureReason::StepSizeUnderflow));
			REQUIRE(result.accepted_steps == 0);
			REQUIRE(result.rejected_steps == 1);
			REQUIRE(result.solution.size() == 1);
			REQUIRE_FALSE(result.error_message.empty());
		}
	}

	TEST_CASE("RODAS4 - singular matrices and nonfinite data return explicit failures", "[dae][RODAS4]") {
		LinearDAE system;
		for (Real derivative : {Real(2), std::numeric_limits<Real>::quiet_NaN()}) {
			system.fx = derivative;
			auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, FixedConfig(1));
			REQUIRE(result.status == (std::isfinite(derivative) ? AlgorithmStatus::SingularMatrix : AlgorithmStatus::NumericalInstability));
			REQUIRE(result.lu_factorizations == 0);
			REQUIRE(result.linear_solves == 0);
			REQUIRE(result.rejected_steps == 1);
			REQUIRE(result.solution.size() == 1);
			REQUIRE_FALSE(result.error_message.empty());
		}
		auto initial = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{0}, 1);
		REQUIRE(initial.status == AlgorithmStatus::InvalidInput);
		REQUIRE(initial.failure_reason == DAEFailureReason::InconsistentInitialConditions);
		system.fx = -2;
		auto config = FixedConfig(1);
		config.abs_tolerance = config.rel_tolerance = std::numeric_limits<Real>::max();
		auto overflow = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
		REQUIRE(overflow.status == AlgorithmStatus::NumericalInstability);
		REQUIRE(overflow.accepted_steps == 0);
		REQUIRE(overflow.rejected_steps == 1);
	}

	TEST_CASE("RODAS4 - clipped final steps and stalled reverse clocks", "[dae][RODAS4][reverse-time]") {
		LinearDAE system;
		auto config = FixedConfig(Real(0.125));
		for (bool reverse : {false, true}) {
			auto result = SolveDAERODAS4(system, reverse ? Real(1) : Real(0), Vector<Real>{1}, Vector<Real>{2},
																	 reverse ? Real(0.96875) : Real(0.03125), config);
			RequireSuccess(result, reverse ? Real(0.96875) : Real(0.03125));
			REQUIRE(result.accepted_steps == 1);
			REQUIRE(result.solution.size() == 2);
		}
		const Real clock = 1 / std::numeric_limits<Real>::epsilon();
		auto stalled = SolveDAERODAS4(system, clock, Vector<Real>{1}, Vector<Real>{2}, clock - 4, config);
		REQUIRE(stalled.status == AlgorithmStatus::Stalled);
		REQUIRE(stalled.failure_reason == DAEFailureReason::StepSizeUnderflow);
		REQUIRE(stalled.solution.size() == 1);
	}

	TEST_CASE("RODAS4 - invalid configuration is rejected before solving", "[dae][RODAS4]") {
		LinearDAE system;
		for (int invalid = 0; invalid < 10; ++invalid) {
			DAESolverConfig config;
			switch (invalid) {
			case 0:
				config.step_size = 0;
				break;
			case 1:
				config.min_step_size = -1;
				break;
			case 2:
				config.max_step_size = config.min_step_size / 2;
				break;
			case 3:
				config.abs_tolerance = 0;
				break;
			case 4:
				config.rel_tolerance = -1;
				break;
			case 5:
				config.constraint_tol = std::numeric_limits<Real>::infinity();
				break;
			case 6:
				config.step_safety = 1;
				break;
			case 7:
				config.max_steps = 0;
				break;
			case 8:
				config.step_size = std::numeric_limits<Real>::quiet_NaN();
				break;
			case 9:
				config.max_step_size = std::numeric_limits<Real>::infinity();
				break;
			}
			CAPTURE(invalid);
			REQUIRE_THROWS_AS(SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), ArgumentError);
		}
		REQUIRE_THROWS_AS(SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 0), ArgumentError);
		REQUIRE_THROWS_AS(SolveDAERODAS4(system, 0, Vector<Real>{1, 2}, Vector<Real>{2}, 1), ArgumentError);
	}
} // namespace MML::Tests::RODAS4
