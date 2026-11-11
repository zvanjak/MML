#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/DAESolvers/DAEBDF4.h>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class BDF4LinearDAE : public IODESystemDAEWithJacobian {
	public:
		Real rate = -1;
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override { f[0] = rate * x[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = rate; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
	};

	class BDF4NonlinearDAE : public BDF4LinearDAE {
	public:
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0] * x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gx(Real, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -2 * x[0]; }
	};

	DAESolverConfig BDF4AccuracyConfig(Real h) {
		DAESolverConfig config;
		config.step_size = h;
		config.newton_tol = TOL(1e-13, 1e-5);
		config.constraint_tol = TOL(1e-10, 1e-4);
		config.max_newton_iter = 50;
		return config;
	}

	void CheckBDF4Statistics(const DAESolverResult& result) {
		CHECK(result.algorithm_name == "DAEBDF4");
		CHECK(result.total_steps == result.accepted_steps);
		CHECK(result.solution.size() == result.accepted_steps + 1);
		CHECK(result.solution.getNumStepsOK() == result.accepted_steps);
		CHECK(result.solution.getNumStepsBad() == result.rejected_steps);
	}

	void CheckFourthOrder(const std::vector<Real>& errors) {
		if constexpr (!std::is_same_v<Real, float>) {
			for (size_t i = 1; i < errors.size(); ++i) {
				const Real ratio = errors[i - 1] / errors[i];
				CAPTURE(i, errors[i - 1], errors[i], ratio);
				CHECK(ratio > Real(12));
				CHECK(ratio < Real(20));
			}
		}
	}
}

TEST_CASE("Radau and BDF4 integrate representable endpoint intervals below the time tolerance",
	"[dae][bdf4][Radau][bdf4-contract]") {
	BDF4LinearDAE system;
	for (bool reverse : {false, true}) {
		const Real start = 1;
		const Real end = std::nextafter(start, reverse ? Real(0) : Real(2));
		auto config = BDF4AccuracyConfig(std::abs(end - start));
		config.max_steps = 1;
		REQUIRE(config.step_size < DAETimeTolerance(end, config.step_size));
		for (bool bdf4 : {false, true}) {
			CAPTURE(reverse, bdf4);
			const auto result = bdf4
				? SolveDAEBDF4(system, start, Vector<Real>{1}, Vector<Real>{1}, end, config)
				: SolveDAERadauIIA(system, start, Vector<Real>{1}, Vector<Real>{1}, end, config);
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			REQUIRE(result.accepted_steps == 1);
			REQUIRE(result.solution.size() == 2);
			REQUIRE(result.solution.getTValue(1) == end);
			REQUIRE(result.final_constraint_norm <= config.constraint_tol);
		}
	}
}

TEST_CASE("BDF4 has fourth-order global convergence in both directions", "[dae][bdf4][bdf4-contract]") {
	BDF4LinearDAE system;
	for (bool reverse : {false, true}) {
		std::vector<Real> errors;
		for (Real h : {Real(0.0625), Real(0.03125), Real(0.015625), Real(0.0078125)}) {
			const Real start = reverse ? Real(1) : Real(0), end = reverse ? Real(0) : Real(1);
			const Real initial = std::exp(-start);
			auto result = SolveDAEBDF4(system, start, Vector<Real>{initial}, Vector<Real>{initial}, end, BDF4AccuracyConfig(h));
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			CheckBDF4Statistics(result);
			CHECK(result.solution.getTValue(result.solution.size() - 1) == end);
			CHECK(result.final_constraint_norm <= BDF4AccuracyConfig(h).constraint_tol);
			errors.push_back(std::abs(result.solution.getXValuesAtEnd()[0] - std::exp(-end)));
		}
		CheckFourthOrder(errors);
	}
}

TEST_CASE("BDF4 retains fourth order with coupled nonlinear constraints", "[dae][bdf4][bdf4-contract]") {
	BDF4NonlinearDAE system;
	for (bool reverse : {false, true}) {
		std::vector<Real> errors;
		for (Real h : {Real(0.0625), Real(0.03125), Real(0.015625), Real(0.0078125)}) {
			const Real start = reverse ? Real(1) : Real(0), end = reverse ? Real(0) : Real(1);
			const Real initial = 1 / (1 + start), exact = 1 / (1 + end);
			auto config = BDF4AccuracyConfig(h);
			auto result = SolveDAEBDF4(system, start, Vector<Real>{initial}, Vector<Real>{initial * initial}, end, config);
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			CheckBDF4Statistics(result);
			errors.push_back(std::abs(result.solution.getXValuesAtEnd()[0] - exact));
			CHECK_THAT(result.solution.getYValuesAtEnd()[0], WithinAbs(exact * exact, 3 * errors.back() + config.constraint_tol));
			CHECK(result.final_constraint_norm <= config.constraint_tol);
		}
		CheckFourthOrder(errors);
	}
}

TEST_CASE("BDF4 clipped final steps preserve fourth-order global accuracy", "[dae][bdf4][bdf4-contract]") {
	BDF4LinearDAE system;
	for (bool reverse : {false, true}) {
		std::vector<Real> errors;
		for (Real h : {Real(0.0625), Real(0.03125), Real(0.015625), Real(0.0078125)}) {
			const Real start = reverse ? Real(1) : Real(0), end = reverse ? -h / 2 : 1 + h / 2;
			const Real initial = std::exp(-start);
			auto result = SolveDAEBDF4(system, start, Vector<Real>{initial}, Vector<Real>{initial}, end, BDF4AccuracyConfig(h));
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			CheckBDF4Statistics(result);
			CHECK(result.solution.getTValue(result.solution.size() - 1) == end);
			errors.push_back(std::abs(result.solution.getXValuesAtEnd()[0] - std::exp(-end)));
		}
		CheckFourthOrder(errors);
	}
}

TEST_CASE("BDF4 startup and endpoint steps integrate quartic nonautonomous states", "[dae][bdf4][bdf4-contract]") {
	class QuarticDAE : public BDF4LinearDAE {
	public:
		void diffEqs(Real t, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 4 * t * t * t; }
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0] - t * t; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
	} system;
	for (bool reverse : {false, true}) {
		for (Real span : {Real(0.0625), Real(0.25), Real(0.375), Real(1), Real(1.0625)}) {
			const Real start = reverse ? span : Real(0), end = reverse ? Real(0) : span;
			const Real initial = start * start * start * start;
			auto result = SolveDAEBDF4(system, start, Vector<Real>{initial}, Vector<Real>{initial + start * start}, end,
				BDF4AccuracyConfig(Real(0.125)));
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			CheckBDF4Statistics(result);
			CHECK(result.solution.getTValue(result.solution.size() - 1) == end);
			for (int i = 0; i < result.solution.size(); ++i) {
				const Real t = result.solution.getTValue(i), exact = t * t * t * t;
				CHECK_THAT(result.solution.getXValue(i, 0), WithinAbs(exact, TOL(1e-11, 1e-4)));
				CHECK_THAT(result.solution.getYValue(i, 0), WithinAbs(exact + t * t, TOL(1e-11, 1e-4)));
			}
		}
	}
}

TEST_CASE("BDF4 rejects failed startup main and endpoint candidates transactionally", "[dae][bdf4][bdf4-contract]") {
	class FailingDAE : public BDF4LinearDAE {
	public:
		Real cutoff = 0;
		bool throwCallback = false;
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			if (t > cutoff) {
				if (throwCallback) throw std::logic_error("BDF4 model failure");
				f[0] = std::numeric_limits<Real>::quiet_NaN();
			}
			else BDF4LinearDAE::diffEqs(t, x, y, f);
		}
	};
	for (Real cutoff : {Real(0), Real(0.125), Real(0.375), Real(1)}) {
		FailingDAE system;
		system.cutoff = cutoff;
		auto result = SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, Real(1.0625),
			BDF4AccuracyConfig(Real(0.125)));
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(result.failure_reason == DAEFailureReason::NonFiniteState);
		CHECK(result.accepted_steps == static_cast<int>(cutoff * 8));
		CHECK(result.rejected_steps == 1);
		CheckBDF4Statistics(result);
		CHECK(result.solution.getTValue(result.solution.size() - 1) == cutoff);
		CHECK(result.final_constraint_norm == 0);
		CHECK_FALSE(result.error_message.empty());
	}
	FailingDAE throwing;
	throwing.throwCallback = true;
	CHECK_THROWS_AS(SolveDAEBDF4(throwing, 0, Vector<Real>{1}, Vector<Real>{1}, 1,
		BDF4AccuracyConfig(Real(0.125))), std::logic_error);
}

TEST_CASE("BDF4 keeps stiff damping and budgets accepted main intervals", "[dae][bdf4][bdf4-contract]") {
	BDF4LinearDAE system;
	system.rate = Real(-1e4);
	auto config = BDF4AccuracyConfig(Real(0.125));
	auto stiff = SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	INFO(stiff.error_message);
	REQUIRE(stiff.status == AlgorithmStatus::Success);
	CheckBDF4Statistics(stiff);
	CHECK(std::abs(stiff.solution.getXValuesAtEnd()[0]) < TOL(1e-6, 1e-3));
	for (int i = 1; i < stiff.solution.size(); ++i) {
		CHECK(std::abs(stiff.solution.getXValue(i, 0)) < Real(0.1));
	}
	system.rate = -1;
	for (int budget : {1, 2, 3, 4, 8}) {
		config.max_steps = budget;
		auto limited = SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
		CHECK(limited.status == (budget == 8 ? AlgorithmStatus::Success : AlgorithmStatus::MaxIterationsExceeded));
		CHECK(limited.accepted_steps == budget);
		CHECK(limited.rejected_steps == 0);
		CHECK(limited.solution.getTValue(limited.solution.size() - 1) == budget * config.step_size);
		CheckBDF4Statistics(limited);
	}
}

TEST_CASE("BDF4 preserves short representable intervals and rejects nonadvancing steps", "[dae][bdf4][bdf4-contract]") {
	BDF4LinearDAE system;
	for (bool reverse : {false, true}) {
		const Real tiny = TOL(1e-12, 1e-3);
		for (Real start : {Real(0), Real(1)}) {
			const Real end = reverse ? start - tiny : start + tiny;
			auto result = SolveDAEBDF4(system, start, Vector<Real>{1}, Vector<Real>{1}, end, BDF4AccuracyConfig(1));
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			CHECK(result.accepted_steps == 1);
			CHECK(result.solution.getTValue(1) == end);
			CheckBDF4Statistics(result);
		}
		const Real start = reverse ? Real(1) + tiny : Real(0), end = reverse ? Real(0) : Real(1) + tiny;
		auto result = SolveDAEBDF4(system, start, Vector<Real>{1}, Vector<Real>{1}, end, BDF4AccuracyConfig(Real(0.125)));
		REQUIRE(result.status == AlgorithmStatus::Success);
		CHECK(result.solution.getTValue(result.solution.size() - 1) == end);
		CheckBDF4Statistics(result);
	}
	const Real shortEnd = TOL(1e-18, 1e-9);
	auto shortResult = SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, shortEnd,
		BDF4AccuracyConfig(shortEnd));
	REQUIRE(shortResult.status == AlgorithmStatus::Success);
	CHECK(shortResult.accepted_steps == 1);
	CHECK(shortResult.solution.getTValue(1) == shortEnd);
	CheckBDF4Statistics(shortResult);
	auto stalled = SolveDAEBDF4(system, 1, Vector<Real>{1}, Vector<Real>{1}, 2,
		BDF4AccuracyConfig(Precision::Eps / Real(4)));
	CHECK(stalled.status == AlgorithmStatus::Stalled);
	CHECK(stalled.failure_reason == DAEFailureReason::StepSizeUnderflow);
	CHECK(stalled.accepted_steps == 0);
	CHECK(stalled.rejected_steps == 0);
	CheckBDF4Statistics(stalled);
}

TEST_CASE("BDF4 validates used inputs before callbacks and ignores adaptive controls", "[dae][bdf4][bdf4-contract]") {
	class CountingDAE : public BDF4LinearDAE {
	public:
		mutable int calls = 0;
		int diffDim = 1, algDim = 1;
		int getDiffDim() const override { return diffDim; }
		int getAlgDim() const override { return algDim; }
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			++calls;
			BDF4LinearDAE::diffEqs(t, x, y, f);
		}
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			++calls;
			BDF4LinearDAE::algConstraints(t, x, y, g);
		}
	} system;
	const Real nan = std::numeric_limits<Real>::quiet_NaN(), inf = std::numeric_limits<Real>::infinity();
	DAESolverConfig config;
	for (Real invalid : {Real(0), Real(-1), nan, inf}) {
		for (int field = 0; field < 3; ++field) {
			config = DAESolverConfig{};
			if (field == 0) config.step_size = invalid;
			else if (field == 1) config.newton_tol = invalid;
			else config.constraint_tol = invalid;
			CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
		}
	}
	for (int invalid : {0, -1}) {
		config = DAESolverConfig{};
		config.max_steps = invalid;
		CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
		config = DAESolverConfig{};
		config.max_newton_iter = invalid;
		CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
	}
	config = DAESolverConfig{};
	for (Real invalid : {nan, inf, Real(0)}) {
		CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, invalid, config), ArgumentError);
	}
	CHECK_THROWS_AS(SolveDAEBDF4(system, nan, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
	CHECK_THROWS_AS(SolveDAEBDF4(system, -std::numeric_limits<Real>::max(), Vector<Real>{1},
		Vector<Real>{1}, std::numeric_limits<Real>::max(), config), ArgumentError);
	for (const Vector<Real>& wrong : {Vector<Real>{}, Vector<Real>{1, 2}}) {
		CHECK_THROWS_AS(SolveDAEBDF4(system, 0, wrong, Vector<Real>{1}, 1, config), ArgumentError);
		CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, wrong, 1, config), ArgumentError);
	}
	system.diffDim = 0;
	CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
	system.diffDim = 1;
	system.algDim = -1;
	CHECK_THROWS_AS(SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config), ArgumentError);
	CHECK(system.calls == 0);
	system.algDim = 1;
	config.min_step_size = inf;
	config.max_step_size = -1;
	config.abs_tolerance = nan;
	config.rel_tolerance = -1;
	config.step_safety = nan;
	auto result = SolveDAEBDF4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	REQUIRE(result.status == AlgorithmStatus::Success);
	CHECK(system.calls > 0);
	CheckBDF4Statistics(result);
}

TEST_CASE("BDF4 accounts for Radau starter work and supports no algebraic variables", "[dae][bdf4][bdf4-contract]") {
	class DifferentialOnlyDAE : public BDF4LinearDAE {
	public:
		int getDiffDim() const override { return 2; }
		int getAlgDim() const override { return 0; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
			f[0] = -x[0];
			f[1] = -2 * x[1];
		}
		void algConstraints(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&) const override {}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override {
			j(0, 0) = -1;
			j(1, 1) = -2;
			j(0, 1) = j(1, 0) = 0;
		}
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
	} system;
	auto config = BDF4AccuracyConfig(Real(0.125));
	auto starter = SolveDAEBDF4(system, 0, Vector<Real>{1, 1}, Vector<Real>{}, config.step_size, config);
	auto reference = SolveDAERadauIIA(system, 0, Vector<Real>{1, 1}, Vector<Real>{}, config.step_size, config);
	REQUIRE(starter.status == AlgorithmStatus::Success);
	REQUIRE(reference.status == AlgorithmStatus::Success);
	CheckBDF4Statistics(starter);
	CHECK(starter.accepted_steps == 1);
	CHECK(starter.newton_iterations == reference.newton_iterations);
	CHECK(starter.jacobian_evaluations == reference.jacobian_evaluations);
	CHECK(starter.lu_factorizations == reference.lu_factorizations);
	CHECK(starter.linear_solves == reference.linear_solves);
	CHECK(starter.last_step_newton_iterations == reference.last_step_newton_iterations);
	CHECK(starter.max_step_newton_iterations == reference.max_step_newton_iterations);
	CHECK(starter.final_residual_norm == reference.final_residual_norm);
	CHECK(starter.solution.getXValuesAtEnd() == reference.solution.getXValuesAtEnd());
	config.step_size = Real(0.015625);
	auto result = SolveDAEBDF4(system, 0, Vector<Real>{1, 1}, Vector<Real>{}, 1, config);
	REQUIRE(result.status == AlgorithmStatus::Success);
	CheckBDF4Statistics(result);
	CHECK(result.solution.getYValuesAtEnd().size() == 0);
	CHECK_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(std::exp(Real(-1)), TOL(1e-7, 1e-4)));
	CHECK_THAT(result.solution.getXValuesAtEnd()[1], WithinAbs(std::exp(Real(-2)), TOL(1e-7, 1e-4)));
}

TEST_CASE("BDF4 propagates singular and nonlinear constraint failures without saving candidates", "[dae][bdf4][bdf4-contract]") {
	class SingularDAE : public BDF4LinearDAE {
	public:
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
	} singular;
	auto result = SolveDAEBDF4(singular, 0, Vector<Real>{1}, Vector<Real>{1}, 1, BDF4AccuracyConfig(Real(0.125)));
	CHECK(result.status == AlgorithmStatus::SingularMatrix);
	CHECK(result.failure_reason == DAEFailureReason::SingularJacobian);
	CHECK(result.accepted_steps == 0);
	CHECK(result.rejected_steps == 1);
	CHECK(result.final_constraint_norm == 0);
	CheckBDF4Statistics(result);
	BDF4NonlinearDAE nonlinear;
	auto config = BDF4AccuracyConfig(Real(0.125));
	config.newton_tol = Real(1e-2);
	auto constrained = SolveDAEBDF4(nonlinear, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
	CHECK(constrained.status == AlgorithmStatus::NumericalInstability);
	CHECK(constrained.failure_reason == DAEFailureReason::ConstraintViolation);
	CHECK(constrained.accepted_steps == 0);
	CHECK(constrained.rejected_steps == 1);
	CHECK(constrained.final_constraint_norm == 0);
	CHECK(constrained.max_constraint_violation > config.constraint_tol);
	CheckBDF4Statistics(constrained);
}
