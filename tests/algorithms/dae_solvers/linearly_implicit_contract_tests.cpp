#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/DAESolvers/DAELinearlyImplicit.h>
#include <mml/algorithms/DAESolvers/DAERODAS4.h>
#endif

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	class TrapezoidLinearDAE : public IODESystemDAEWithJacobian {
	public:
		Real rate = -1;
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override { f[0] = rate * x[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - 2 * x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = rate; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -2; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
	};

	class TrapezoidNonlinearDAE : public TrapezoidLinearDAE {
	public:
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0] * x[0]; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
		void jacobian_gx(Real, const Vector<Real>& x, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -2 * x[0]; }
	};

	class TrapezoidTimedDAE : public TrapezoidLinearDAE {
	public:
		void diffEqs(Real t, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = t; }
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0] - t * t; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
	};

	class TrapezoidComplexDAE : public IODESystemDAEWithJacobian {
	public:
		Real a = 0, b = 0;
		int getDiffDim() const override { return 2; }
		int getAlgDim() const override { return 0; }
		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override {
			f[0] = a * x[0] - b * x[1];
			f[1] = b * x[0] + a * x[1];
		}
		void algConstraints(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&) const override {}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override {
			j(0, 0) = j(1, 1) = a;
			j(0, 1) = -b;
			j(1, 0) = b;
		}
		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&) const override {}
	};

	Real TrapezoidAmplification(Real z) { return (1 + z / 2) / (1 - z / 2); }

	void CheckTrapezoidStatistics(const DAESolverResult& result) {
		CHECK(result.total_steps == result.accepted_steps);
		CHECK(result.solution.size() == result.accepted_steps + 1);
		CHECK(result.solution.getNumStepsOK() == result.accepted_steps);
		CHECK(result.solution.getNumStepsBad() == result.rejected_steps);
		CHECK(result.newton_iterations == 0);
	}

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4996)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	DAESolverResult LegacyTrapezoid(IODESystemDAEWithJacobian& system, const DAESolverConfig& config) {
		return SolveDAERODAS(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	}
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
}

TEST_CASE("Linearly implicit DAE has the exact trapezoidal amplification function", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidLinearDAE system;
	DAESolverConfig config;
	config.step_size = 1;
	for (Real z : {Real(-1e6), Real(-100), Real(-10), Real(-2), Real(-1), Real(-0.125), Real(0), Real(0.5), Real(1)}) {
		CAPTURE(z);
		system.rate = z;
		auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		CHECK(result.algorithm_name == "DAELinearlyImplicit");
		CHECK(result.accepted_steps == 1);
		CHECK(result.rejected_steps == 0);
		CHECK(result.lu_factorizations == 1);
		CHECK(result.linear_solves == 2);
		CHECK(result.jacobian_evaluations == 1);
		CheckTrapezoidStatistics(result);
		const Real expected = TrapezoidAmplification(z);
		CHECK_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(expected, TOL(1e-10, 1e-4)));
		CHECK_THAT(result.solution.getYValuesAtEnd()[0], WithinAbs(2 * expected, TOL(1e-10, 1e-4)));
	}
}

TEST_CASE("Linearly implicit DAE converges at order two on linear and nonlinear constraints", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidLinearDAE linear;
	TrapezoidNonlinearDAE nonlinear;
	for (bool curved : {false, true}) {
		IODESystemDAEWithJacobian& system = curved ? static_cast<IODESystemDAEWithJacobian&>(nonlinear) : linear;
		std::vector<Real> errors;
		for (Real h : {Real(0.125), Real(0.0625), Real(0.03125), Real(0.015625)}) {
			DAESolverConfig config;
			config.step_size = h;
			config.constraint_tol = Real(1e-2);
			auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{curved ? Real(1) : Real(2)}, 1, config);
			INFO(result.error_message);
			REQUIRE(result.status == AlgorithmStatus::Success);
			const Real exactX = curved ? Real(0.5) : std::exp(Real(-1));
			const Real exactY = curved ? Real(0.25) : 2 * exactX;
			const Real error = std::abs(result.solution.getXValuesAtEnd()[0] - exactX);
			errors.push_back(error);
			CHECK_THAT(result.solution.getYValuesAtEnd()[0], WithinAbs(exactY, 4 * error + TOL(1e-10, 1e-4)));
			CHECK(result.final_constraint_norm <= config.constraint_tol);
			CheckTrapezoidStatistics(result);
			CHECK(result.lu_factorizations == result.accepted_steps);
			CHECK(result.linear_solves == 2 * result.accepted_steps);
			INFO("curved=" << curved << ", h=" << h << ", error=" << error);
			if (errors.size() > 1) {
				const Real ratio = errors[errors.size() - 2] / errors.back();
				INFO("observed ratio=" << ratio);
				CHECK(ratio > Real(3.6));
				CHECK(ratio < Real(4.4));
			}
		}
	}
}

TEST_CASE("Linearly implicit DAE is A-stable on complex modes including the imaginary axis", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidComplexDAE system;
	DAESolverConfig config;
	config.step_size = 1;
	for (const auto& mode : std::vector<std::pair<Real, Real>>{{Real(0), Real(1)}, {Real(0), Real(10)}, {Real(-0.25), Real(4)}, {Real(-100), Real(25)}}) {
		system.a = mode.first;
		system.b = mode.second;
		const std::complex<Real> z{system.a, system.b};
		const auto expected = (Real(1) + z / Real(2)) / (Real(1) - z / Real(2));
		auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1, 0}, Vector<Real>{}, 1, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		const auto end = result.solution.getXValuesAtEnd();
		CHECK_THAT(end[0], WithinAbs(expected.real(), TOL(1e-10, 1e-4)));
		CHECK_THAT(end[1], WithinAbs(expected.imag(), TOL(1e-10, 1e-4)));
		CHECK(std::hypot(end[0], end[1]) <= 1 + TOL(1e-10, 1e-4));
		if (system.a == 0)
			CHECK_THAT(std::hypot(end[0], end[1]), WithinAbs(1, TOL(1e-10, 1e-4)));
	}
}

TEST_CASE("Linearly implicit trapezoid is explicitly not L-stable while RODAS4 damps stiff decay", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidLinearDAE system;
	system.rate = Real(-1e6);
	DAESolverConfig config;
	config.step_size = config.min_step_size = config.max_step_size = 1;
	config.abs_tolerance = 1;
	config.rel_tolerance = 0;
	auto trapezoid = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	REQUIRE(trapezoid.status == AlgorithmStatus::Success);
	CHECK(trapezoid.solution.getXValuesAtEnd()[0] < Real(-0.99));
	CHECK_THAT(trapezoid.solution.getXValuesAtEnd()[0],
		WithinAbs(TrapezoidAmplification(system.rate), TOL(1e-10, 1e-4)));
	auto rodas = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	REQUIRE(rodas.status == AlgorithmStatus::Success);
	CHECK(std::abs(rodas.solution.getXValuesAtEnd()[0]) < TOL(1e-4, 1e-3));
}

TEST_CASE("Linearly implicit DAE handles nonautonomous equations and nonlinear time constraints in both directions", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidTimedDAE system;
	DAESolverConfig config;
	config.step_size = Real(0.125);
	for (bool reverse : {false, true}) {
		const Real start = reverse ? Real(1) : Real(0);
		const Real end = reverse ? Real(0) : Real(1);
		auto result = SolveDAELinearlyImplicit(system, start, Vector<Real>{start * start / 2}, Vector<Real>{3 * start * start / 2}, end, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		CheckTrapezoidStatistics(result);
		CHECK(result.solution.getTValue(result.solution.size() - 1) == end);
		for (int i = 0; i < result.solution.size(); ++i) {
			const Real t = result.solution.getTValue(i);
			CHECK_THAT(result.solution.getXValue(i, 0), WithinAbs(t * t / 2, TOL(1e-12, 1e-5)));
			CHECK_THAT(result.solution.getYValue(i, 0), WithinAbs(3 * t * t / 2, TOL(1e-12, 1e-5)));
		}
	}
}

TEST_CASE("Linearly implicit DAE clipped steps and the legacy alias preserve method identity", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidLinearDAE system;
	DAESolverConfig config;
	config.step_size = Real(0.125);
	for (bool reverse : {false, true}) {
		const Real start = reverse ? Real(1) : Real(0);
		const Real end = reverse ? Real(0.7) : Real(0.3);
		auto result = SolveDAELinearlyImplicit(system, start, Vector<Real>{1}, Vector<Real>{2}, end, config);
		REQUIRE(result.status == AlgorithmStatus::Success);
		CheckTrapezoidStatistics(result);
		CHECK(result.solution.getTValue(result.solution.size() - 1) == end);
		Real expected = 1;
		for (int i = 1; i < result.solution.size(); ++i) {
			const Real h = result.solution.getTValue(i) - result.solution.getTValue(i - 1);
			expected *= TrapezoidAmplification(system.rate * h);
			CHECK_THAT(result.solution.getXValue(i, 0), WithinAbs(expected, TOL(1e-12, 1e-5)));
			CHECK_THAT(result.solution.getYValue(i, 0), WithinAbs(2 * expected, TOL(1e-12, 1e-5)));
		}
	}
	auto modern = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	auto legacy = LegacyTrapezoid(system, config);
	CHECK(legacy.algorithm_name == modern.algorithm_name);
	CHECK(legacy.status == modern.status);
	CHECK(legacy.accepted_steps == modern.accepted_steps);
	CHECK(legacy.solution.getXValuesAtEnd()[0] == modern.solution.getXValuesAtEnd()[0]);
	CHECK(legacy.solution.getYValuesAtEnd()[0] == modern.solution.getYValuesAtEnd()[0]);
}

TEST_CASE("Linearly implicit DAE fails nonfinite stages without committing the candidate", "[dae][linearly-implicit][trapezoid-contract]") {
	class LateFailureDAE : public TrapezoidLinearDAE {
	public:
		bool throwCallback = false;
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			if (t > Real(0.125)) {
				if (throwCallback) throw std::logic_error("Model callback failure");
				f[0] = std::numeric_limits<Real>::quiet_NaN();
			}
			else TrapezoidLinearDAE::diffEqs(t, x, y, f);
		}
	} system;
	DAESolverConfig config;
	config.step_size = Real(0.125);
	auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.failure_reason == DAEFailureReason::NonFiniteState);
	CHECK(result.accepted_steps == 1);
	CHECK(result.rejected_steps == 1);
	CheckTrapezoidStatistics(result);
	REQUIRE(result.solution.size() == 2);
	CHECK(result.solution.getTValue(1) == config.step_size);
	CHECK_THAT(result.solution.getXValuesAtEnd()[0], WithinAbs(TrapezoidAmplification(-config.step_size), TOL(1e-12, 1e-5)));
	CHECK(result.final_constraint_norm <= config.constraint_tol);
	CHECK_FALSE(result.error_message.empty());
	system.throwCallback = true;
	CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), std::logic_error);
}

TEST_CASE("Linearly implicit DAE constraint failures retain stored-state diagnostics", "[dae][linearly-implicit][trapezoid-contract]") {
	class ConstraintDAE : public TrapezoidLinearDAE {
	public:
		void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& f) const override { f[0] = 0; }
		void algConstraints(Real t, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] * y[0] - 1 - t; }
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>& y, Matrix<Real>& j) const override { j(0, 0) = 2 * y[0]; }
	} system;
	DAESolverConfig config;
	config.step_size = Real(0.125);
	auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{0}, Vector<Real>{1}, config.step_size, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.failure_reason == DAEFailureReason::ConstraintViolation);
	CHECK(result.accepted_steps == 0);
	CHECK(result.rejected_steps == 1);
	CheckTrapezoidStatistics(result);
	REQUIRE(result.solution.size() == 1);
	CHECK(result.solution.getXValuesAtEnd()[0] == 0);
	CHECK(result.solution.getYValuesAtEnd()[0] == 1);
	CHECK(result.final_constraint_norm == 0);
	CHECK(result.max_constraint_violation > config.constraint_tol);
}

TEST_CASE("Linearly implicit DAE validates only its used controls before callbacks", "[dae][linearly-implicit][trapezoid-contract]") {
	class CountingDAE : public TrapezoidLinearDAE {
	public:
		mutable int calls = 0;
		int diffDim = 1, algDim = 1;
		int getDiffDim() const override { return diffDim; }
		int getAlgDim() const override { return algDim; }
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			++calls;
			TrapezoidLinearDAE::diffEqs(t, x, y, f);
		}
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			++calls;
			TrapezoidLinearDAE::algConstraints(t, x, y, g);
		}
	} system;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const Real inf = std::numeric_limits<Real>::infinity();
	DAESolverConfig config;
	for (Real invalid : {Real(0), Real(-1), nan, inf}) {
		config.step_size = invalid;
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), ArgumentError);
		config = DAESolverConfig{};
		config.constraint_tol = invalid;
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), ArgumentError);
		config = DAESolverConfig{};
	}
	for (int invalid : {0, -1}) {
		config.max_steps = invalid;
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), ArgumentError);
	}
	config = DAESolverConfig{};
	for (const auto& endpoints : std::vector<std::pair<Real, Real>>{{Real(0), Real(0)}, {nan, Real(1)}, {Real(0), nan}, {inf, Real(1)}, {Real(0), inf},
		{-std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()}}) {
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, endpoints.first, Vector<Real>{1}, Vector<Real>{2},
			endpoints.second, config), ArgumentError);
	}
	for (const Vector<Real>& wrong : {Vector<Real>{}, Vector<Real>{1, 2}}) {
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, wrong, Vector<Real>{2}, 1, config), ArgumentError);
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, wrong, 1, config), ArgumentError);
	}
	for (int invalid : {0, -1}) {
		system.diffDim = invalid;
		CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), ArgumentError);
	}
	system.diffDim = 1;
	system.algDim = -1;
	CHECK_THROWS_AS(SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config), ArgumentError);
	CHECK(system.calls == 0);
	system.algDim = 1;
	config.max_newton_iter = 0;
	config.newton_tol = nan;
	config.abs_tolerance = nan;
	config.rel_tolerance = -1;
	config.min_step_size = inf;
	config.max_step_size = -1;
	config.step_safety = nan;
	config.reuse_jacobian = false;
	auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	REQUIRE(result.status == AlgorithmStatus::Success);
	CHECK(system.calls > 0);
	CheckTrapezoidStatistics(result);
}

TEST_CASE("Linearly implicit DAE numerical stage failures keep the accepted trajectory", "[dae][linearly-implicit][trapezoid-contract]") {
	enum class Failure { Derivative, PredictorDerivative, Jacobian, Singular, Overflow, FinalConstraint };
	class BrokenStageDAE : public TrapezoidLinearDAE {
	public:
		Failure failure = Failure::Derivative;
		mutable int constraintCalls = 0;
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& f) const override {
			if (failure == Failure::Derivative || (failure == Failure::PredictorDerivative && t > 0))
				f[0] = std::numeric_limits<Real>::quiet_NaN();
			else if (failure == Failure::Overflow)
				f[0] = std::numeric_limits<Real>::max();
			else TrapezoidLinearDAE::diffEqs(t, x, y, f);
		}
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override {
			++constraintCalls;
			if (failure == Failure::FinalConstraint && constraintCalls == 4)
				g[0] = std::numeric_limits<Real>::quiet_NaN();
			else TrapezoidLinearDAE::algConstraints(t, x, y, g);
		}
		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override {
			j(0, 0) = failure == Failure::Jacobian ? std::numeric_limits<Real>::quiet_NaN() :
				(failure == Failure::Singular ? Real(2) : Real(0));
		}
	};
	for (auto failure : {Failure::Derivative, Failure::PredictorDerivative, Failure::Jacobian,
		Failure::Singular, Failure::Overflow, Failure::FinalConstraint}) {
		CAPTURE(static_cast<int>(failure));
		BrokenStageDAE system;
		system.failure = failure;
		DAESolverConfig config;
		config.step_size = failure == Failure::Overflow ? Real(4) : Real(1);
		auto result = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, config.step_size, config);
		CHECK(result.status == (failure == Failure::Singular ?
			AlgorithmStatus::SingularMatrix : AlgorithmStatus::NumericalInstability));
		CHECK(result.failure_reason == (failure == Failure::Singular ?
			DAEFailureReason::SingularJacobian : DAEFailureReason::NonFiniteState));
		CHECK(result.accepted_steps == 0);
		CHECK(result.rejected_steps == 1);
		CheckTrapezoidStatistics(result);
		REQUIRE(result.solution.size() == 1);
		CHECK(result.solution.getTValue(0) == 0);
		CHECK(result.solution.getXValuesAtEnd()[0] == 1);
		CHECK(result.solution.getYValuesAtEnd()[0] == 2);
		CHECK(result.final_constraint_norm == 0);
		CHECK_FALSE(result.error_message.empty());
	}
}

TEST_CASE("Linearly implicit DAE preserves initial-state and fixed-step budget contracts", "[dae][linearly-implicit][trapezoid-contract]") {
	TrapezoidLinearDAE system;
	DAESolverConfig config;
	config.step_size = Real(0.125);
	config.max_steps = 1;
	auto limited = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{2}, 1, config);
	CHECK(limited.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(limited.failure_reason == DAEFailureReason::MaxStepsExceeded);
	CHECK(limited.accepted_steps == 1);
	CHECK(limited.rejected_steps == 0);
	CheckTrapezoidStatistics(limited);
	CHECK(limited.solution.getTValue(1) == config.step_size);
	auto inconsistent = SolveDAELinearlyImplicit(system, 0, Vector<Real>{1}, Vector<Real>{3}, 1, config);
	CHECK(inconsistent.status == AlgorithmStatus::InvalidInput);
	CHECK(inconsistent.failure_reason == DAEFailureReason::InconsistentInitialConditions);
	CheckTrapezoidStatistics(inconsistent);
	CHECK(inconsistent.accepted_steps == 0);
	CHECK(inconsistent.rejected_steps == 0);
	auto nonfinite = SolveDAELinearlyImplicit(system, 0,
		Vector<Real>{std::numeric_limits<Real>::quiet_NaN()}, Vector<Real>{2}, 1, config);
	CHECK(nonfinite.status == AlgorithmStatus::NumericalInstability);
	CHECK(nonfinite.failure_reason == DAEFailureReason::NonFiniteState);
	CheckTrapezoidStatistics(nonfinite);
	CHECK(nonfinite.accepted_steps == 0);
	CHECK(nonfinite.rejected_steps == 0);
}
