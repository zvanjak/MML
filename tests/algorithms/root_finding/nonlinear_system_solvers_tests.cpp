#include <catch2/catch_all.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/RootFinding/NonlinearSystemSolvers.h>
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::NonlinearSystemSolverTests {

	RootFinding::DynamicSystemFunction CircleLineSystem() {
		return [](const Vector<Real>& point) {
			return Vector<Real>({point[0] * point[0] + point[1] * point[1] - REAL(1.0),
				point[0] - point[1]});
		};
	}

	RootFinding::DynamicJacobianFunction CircleLineJacobian() {
		return [](const Vector<Real>& point) {
			return Matrix<Real>(2, 2, {REAL(2.0) * point[0], REAL(2.0) * point[1],
				REAL(1.0), -REAL(1.0)});
		};
	}

	TEST_CASE("Nonlinear Newton solves dynamic system with analytical Jacobian",
		"[RootFinding][NonlinearSystem][AnalyticalJacobian]") {
		RootFinding::NonlinearSystemConfig config;
		config.residual_tolerance = TOL(1e-12, 1e-5);
		config.step_tolerance = TOL(1e-13, 1e-6);
		config.store_trace = true;
		auto result = RootFinding::SolveNonlinearSystemNewton(
			CircleLineSystem(), CircleLineJacobian(), Vector<Real>({REAL(0.8), REAL(0.4)}), config);

		const Real expected = REAL(1.0) / std::sqrt(REAL(2.0));
		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "NonlinearNewtonAnalyticalJacobian");
		REQUIRE_THAT(result.solution[0], WithinAbs(expected, TOL(1e-10, 1e-4)));
		REQUIRE_THAT(result.solution[1], WithinAbs(expected, TOL(1e-10, 1e-4)));
		REQUIRE(result.residual_norm <= config.residual_tolerance);
		REQUIRE(result.iterations_used > 0);
		REQUIRE(result.jacobian_evaluations == result.iterations_used);
		REQUIRE(result.linear_solves == result.iterations_used);
		REQUIRE(result.function_evaluations >= result.iterations_used + 1);
		REQUIRE(result.trace.size() == static_cast<size_t>(result.iterations_used));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		for (size_t index = 1; index < result.trace.size(); ++index)
			REQUIRE(result.trace[index].residual_norm < result.trace[index - 1].residual_norm);
	}

	TEST_CASE("Nonlinear Newton computes dynamic numerical Jacobian",
		"[RootFinding][NonlinearSystem][NumericalJacobian]") {
		RootFinding::NonlinearSystemConfig config;
		config.residual_tolerance = TOL(1e-10, 1e-4);
		auto result = RootFinding::SolveNonlinearSystemNewton(
			CircleLineSystem(), Vector<Real>({REAL(0.8), REAL(0.4)}), config);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "NonlinearNewtonNumericalJacobian");
		REQUIRE(result.jacobian_evaluations == result.iterations_used);
		REQUIRE(result.function_evaluations >= 4 * 2 * result.jacobian_evaluations + 1);
		REQUIRE(result.residual_norm <= config.residual_tolerance);
	}

	TEST_CASE("Nonlinear Newton supports fixed-size callables and Jacobians",
		"[RootFinding][NonlinearSystem][Static][AnalyticalJacobian]") {
		auto function = [](const VectorN<Real, 2>& point) {
			return VectorN<Real, 2>({point[0] * point[0] + point[1] * point[1] - REAL(1.0),
				point[0] - point[1]});
		};
		auto jacobian = [](const VectorN<Real, 2>& point) {
			return MatrixNM<Real, 2, 2>({REAL(2.0) * point[0], REAL(2.0) * point[1],
				REAL(1.0), -REAL(1.0)});
		};
		RootFinding::NonlinearSystemConfig config;
		config.store_trace = true;
		auto result = RootFinding::SolveNonlinearSystemNewton<2>(
			function, jacobian, VectorN<Real, 2>({REAL(0.8), REAL(0.4)}), config);

		const Real expected = REAL(1.0) / std::sqrt(REAL(2.0));
		REQUIRE(result.IsSuccess());
		REQUIRE_THAT(result.solution[0], WithinAbs(expected, TOL(1e-9, 1e-4)));
		REQUIRE_THAT(result.solution[1], WithinAbs(expected, TOL(1e-9, 1e-4)));
		REQUIRE(result.trace.size() == static_cast<size_t>(result.iterations_used));
	}

	class StaticCircleLineSystem : public IVectorFunction<2> {
	public:
		VectorN<Real, 2> operator()(const VectorN<Real, 2>& point) const override {
			return {point[0] * point[0] + point[1] * point[1] - REAL(1.0),
				point[0] - point[1]};
		}
	};

	TEST_CASE("Nonlinear Newton supports fixed-size interface with numerical Jacobian",
		"[RootFinding][NonlinearSystem][Static][NumericalJacobian]") {
		StaticCircleLineSystem function;
		RootFinding::NonlinearSystemConfig config;
		auto result = RootFinding::SolveNonlinearSystemNewton<2>(
			function, VectorN<Real, 2>({REAL(0.8), REAL(0.4)}), config);
		REQUIRE(result.IsSuccess());
		REQUIRE(result.jacobian_evaluations == result.iterations_used);
		REQUIRE(result.residual_norm <= config.residual_tolerance);
	}

	TEST_CASE("Nonlinear Newton accepts an exact initial root",
		"[RootFinding][NonlinearSystem][InitialRoot]") {
		const Real expected = REAL(1.0) / std::sqrt(REAL(2.0));
		auto result = RootFinding::SolveNonlinearSystemNewton(
			CircleLineSystem(), Vector<Real>({expected, expected}));
		REQUIRE(result.IsSuccess());
		REQUIRE(result.iterations_used == 0);
		REQUIRE(result.jacobian_evaluations == 0);
		REQUIRE(result.linear_solves == 0);
		REQUIRE(result.function_evaluations == 1);
	}

	TEST_CASE("Nonlinear Newton backtracking globalizes difficult steps",
		"[RootFinding][NonlinearSystem][Backtracking]") {
		RootFinding::DynamicSystemFunction function = [](const Vector<Real>& point) {
			return Vector<Real>({std::atan(point[0])});
		};
		RootFinding::DynamicJacobianFunction jacobian = [](const Vector<Real>& point) {
			return Matrix<Real>(1, 1, {REAL(1.0) / (REAL(1.0) + point[0] * point[0])});
		};
		RootFinding::NonlinearSystemConfig config;
		config.store_trace = true;
		config.max_iterations = 40;
		config.minimum_damping = REAL(1e-8);
		auto result = RootFinding::SolveNonlinearSystemNewton(
			function, jacobian, Vector<Real>({REAL(10.0)}), config);

		REQUIRE(result.IsSuccess());
		REQUIRE_THAT(result.solution[0], WithinAbs(REAL(0.0), TOL(1e-9, 1e-4)));
		REQUIRE(std::any_of(result.trace.begin(), result.trace.end(),
			[](const auto& entry) { return entry.damping < REAL(1.0); }));
	}

	TEST_CASE("Nonlinear Newton reports structured failures",
		"[RootFinding][NonlinearSystem][Validation]") {
		RootFinding::NonlinearSystemConfig invalid;
		invalid.residual_tolerance = -REAL(1.0);
		auto invalidResult = RootFinding::SolveNonlinearSystemNewton(
			CircleLineSystem(), Vector<Real>({REAL(0.8), REAL(0.4)}), invalid);
		REQUIRE_FALSE(invalidResult.IsSuccess());
		REQUIRE(invalidResult.status == AlgorithmStatus::InvalidInput);

		RootFinding::DynamicSystemFunction empty;
		auto emptyResult = RootFinding::SolveNonlinearSystemNewton(
			empty, Vector<Real>({REAL(0.0)}));
		REQUIRE(emptyResult.status == AlgorithmStatus::InvalidInput);

		RootFinding::DynamicSystemFunction wrongDimension = [](const Vector<Real>&) {
			return Vector<Real>({REAL(1.0)});
		};
		auto dimensionResult = RootFinding::SolveNonlinearSystemNewton(
			wrongDimension, Vector<Real>({REAL(0.0), REAL(0.0)}));
		REQUIRE(dimensionResult.status == AlgorithmStatus::InvalidInput);

		RootFinding::DynamicSystemFunction nonFinite = [](const Vector<Real>& point) {
			return Vector<Real>({std::numeric_limits<Real>::infinity(), point[1]});
		};
		auto nonFiniteResult = RootFinding::SolveNonlinearSystemNewton(
			nonFinite, Vector<Real>({REAL(0.0), REAL(0.0)}));
		REQUIRE(nonFiniteResult.status == AlgorithmStatus::NumericalInstability);

		RootFinding::DynamicJacobianFunction wrongJacobian = [](const Vector<Real>&) {
			return Matrix<Real>(1, 1, {REAL(1.0)});
		};
		auto jacobianResult = RootFinding::SolveNonlinearSystemNewton(
			CircleLineSystem(), wrongJacobian, Vector<Real>({REAL(0.8), REAL(0.4)}));
		REQUIRE(jacobianResult.status == AlgorithmStatus::InvalidInput);

		RootFinding::DynamicSystemFunction singular = [](const Vector<Real>& point) {
			return Vector<Real>({point[0] * point[0] + REAL(1.0), point[1]});
		};
		RootFinding::DynamicJacobianFunction singularJacobian = [](const Vector<Real>& point) {
			return Matrix<Real>(2, 2, {REAL(2.0) * point[0], REAL(0.0), REAL(0.0), REAL(1.0)});
		};
		auto singularResult = RootFinding::SolveNonlinearSystemNewton(
			singular, singularJacobian, Vector<Real>({REAL(0.0), REAL(1.0)}));
		REQUIRE_FALSE(singularResult.IsSuccess());
		REQUIRE(singularResult.status == AlgorithmStatus::SingularMatrix);
		REQUIRE_FALSE(singularResult.error_message.empty());
	}

	TEST_CASE("Nonlinear Newton reports iteration exhaustion",
		"[RootFinding][NonlinearSystem][MaxIterations]") {
		RootFinding::NonlinearSystemConfig config;
		config.max_iterations = 1;
		config.residual_tolerance = TOL(1e-15, 1e-6);
		auto result = RootFinding::SolveNonlinearSystemNewton(
			CircleLineSystem(), CircleLineJacobian(), Vector<Real>({REAL(0.8), REAL(0.4)}), config);
		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(result.iterations_used == 1);
		REQUIRE_FALSE(result.error_message.empty());
	}

} // namespace MML::Tests::Algorithms::NonlinearSystemSolverTests
