///////////////////////////////////////////////////////////////////////////////////////////
///  Tests for Optimization/Constraints/AugmentedLagrangian.h                           ///
///  (repatriated from MML-Packages 2026-10-03, adapted to the MML result contract)     ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h>
#endif

using namespace MML;
using namespace MML::Optimization;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::AugmentedLagrangianTests
{

namespace {
	// Classic Rosenbrock, local to this test file
	struct Rosenbrock2D {
		Real operator()(const VectorN<Real, 2>& x) const {
			Real a = REAL(1.0) - x[0];
			Real b = x[1] - x[0] * x[0];
			return a * a + REAL(100.0) * b * b;
		}
	};
}

///////////////////////////////////////////////////////////////////////////////////////////
// ALMConfig Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: ALMConfig default values", "[ALM][config]")
{
	ALMConfig config;

	REQUIRE(config.max_outer_iterations == 100);
	REQUIRE_THAT(config.feasibility_tolerance, WithinAbs(1e-6, 1e-12));
	REQUIRE_THAT(config.stationarity_tolerance, WithinAbs(1e-6, 1e-12));
	REQUIRE_THAT(config.initial_penalty, WithinAbs(10.0, 1e-12));
	REQUIRE_THAT(config.penalty_multiplier, WithinAbs(10.0, 1e-12));
	REQUIRE_THAT(config.max_penalty, WithinAbs(1e12, 1e6));
	REQUIRE(config.max_inner_iterations == 1000);
	REQUIRE_THAT(config.inner_tolerance, WithinAbs(1e-8, 1e-12));
	REQUIRE(config.verbose_stream == nullptr);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Construction Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: Default construction", "[ALM][construction]")
{
	AugmentedLagrangian<2> alm;
	auto config = alm.GetConfig();
	REQUIRE(config.max_outer_iterations == 100);
}

TEST_CASE("ALM: Construction with config", "[ALM][construction]")
{
	ALMConfig config;
	config.max_outer_iterations = 50;
	config.initial_penalty = REAL(100.0);

	AugmentedLagrangian<2> alm(config);

	auto returnedConfig = alm.GetConfig();
	REQUIRE(returnedConfig.max_outer_iterations == 50);
	REQUIRE_THAT(returnedConfig.initial_penalty, WithinAbs(100.0, 1e-12));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Simple Constrained Optimization Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: Single inequality constraint", "[ALM][optimization]")
{
	// min x^2 + y^2 subject to x + y >= 1 (i.e., -x - y + 1 <= 0)
	// Optimal: x = y = 0.5, f* = 0.5

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	auto constraint = [](const VectorN<Real, 2>& x) -> Real {
		return -x[0] - x[1] + REAL(1.0);  // -x - y + 1 <= 0 means x + y >= 1
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 50;

	AugmentedLagrangian<2> alm(config);
	alm.AddInequalityConstraint(constraint);
	alm.SetBounds(-10.0, 10.0);

	VectorN<Real, 2> x0{ 0.0, 0.0 };
	auto result = alm.Optimize(objective, x0);

	REQUIRE(result.converged);
	REQUIRE_THAT(result.x[0], WithinAbs(0.5, 0.05));
	REQUIRE_THAT(result.x[1], WithinAbs(0.5, 0.05));
	REQUIRE_THAT(result.objective_value, WithinAbs(0.5, 0.05));
}

TEST_CASE("ALM: Single equality constraint", "[ALM][optimization]")
{
	// min x^2 + y^2 subject to x + y = 1
	// Optimal: x = y = 0.5, f* = 0.5

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	auto constraint = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] + x[1] - REAL(1.0);  // x + y = 1
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 50;

	AugmentedLagrangian<2> alm(config);
	alm.AddEqualityConstraint(constraint);
	alm.SetBounds(-10.0, 10.0);

	VectorN<Real, 2> x0{ 2.0, 2.0 };
	auto result = alm.Optimize(objective, x0);

	REQUIRE(result.converged);
	REQUIRE_THAT(result.x[0], WithinAbs(0.5, 0.05));
	REQUIRE_THAT(result.x[1], WithinAbs(0.5, 0.05));
	REQUIRE_THAT(result.objective_value, WithinAbs(0.5, 0.05));
}

TEST_CASE("ALM: Multiple constraints", "[ALM][optimization]")
{
	// min -x - y subject to x^2 + y^2 <= 1, x >= 0, y >= 0
	// Optimal: x = y = 1/sqrt(2) ~ 0.707, f* ~ -1.414

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return -x[0] - x[1];
	};

	auto circleConstraint = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1] - REAL(1.0);
	};
	auto xPositive = [](const VectorN<Real, 2>& x) -> Real { return -x[0]; };
	auto yPositive = [](const VectorN<Real, 2>& x) -> Real { return -x[1]; };

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 100;

	AugmentedLagrangian<2> alm(config);
	alm.AddInequalityConstraint(circleConstraint);
	alm.AddInequalityConstraint(xPositive);
	alm.AddInequalityConstraint(yPositive);
	alm.SetBounds(-2.0, 2.0);

	VectorN<Real, 2> x0{ 0.1, 0.1 };
	auto result = alm.Optimize(objective, x0);

	Real sqrtHalf = REAL(1.0) / std::sqrt(REAL(2.0));

	REQUIRE(result.converged);
	REQUIRE_THAT(result.x[0], WithinAbs(sqrtHalf, REAL(0.1)));
	REQUIRE_THAT(result.x[1], WithinAbs(sqrtHalf, REAL(0.1)));
	REQUIRE_THAT(result.objective_value, WithinAbs(-std::sqrt(REAL(2.0)), REAL(0.1)));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Classic Test Problems
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: Rosenbrock with linear constraint", "[ALM][optimization]")
{
	// min (1-x)^2 + 100(y-x^2)^2 subject to x + y >= 2

	Rosenbrock2D rosenbrock;

	auto constraint = [](const VectorN<Real, 2>& x) -> Real {
		return -(x[0] + x[1] - REAL(2.0));  // x + y >= 2
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 100;
	config.max_inner_iterations = 2000;

	AugmentedLagrangian<2> alm(config);
	alm.AddInequalityConstraint(constraint);
	alm.SetBounds(-5.0, 5.0);

	VectorN<Real, 2> x0{ 0.0, 0.0 };
	auto result = alm.Optimize(rosenbrock, x0);

	REQUIRE(result.x[0] + result.x[1] >= 2.0 - 0.01);
	REQUIRE(result.outer_iterations > 0);
}

TEST_CASE("ALM: Quadratic with ellipsoidal constraint", "[ALM][optimization]")
{
	// min (x-2)^2 + (y-2)^2 subject to x^2/4 + y^2/9 <= 1
	// The unconstrained minimum (2,2) is outside the ellipse

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return (x[0] - REAL(2.0)) * (x[0] - REAL(2.0)) + (x[1] - REAL(2.0)) * (x[1] - REAL(2.0));
	};

	auto ellipse = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] / REAL(4.0) + x[1] * x[1] / REAL(9.0) - REAL(1.0);
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 100;

	AugmentedLagrangian<2> alm(config);
	alm.AddInequalityConstraint(ellipse);
	alm.SetBounds(-5.0, 5.0);

	VectorN<Real, 2> x0{ 0.5, 0.5 };
	auto result = alm.Optimize(objective, x0);

	Real ellipseVal = result.x[0] * result.x[0] / REAL(4.0) + result.x[1] * result.x[1] / REAL(9.0);
	REQUIRE(ellipseVal <= REAL(1.0) + REAL(0.01));
	REQUIRE(result.x[0] > 0);
	REQUIRE(result.x[1] > 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Higher Dimensional Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: 3D optimization with sphere constraint", "[ALM][optimization]")
{
	// min -x - y - z subject to x^2 + y^2 + z^2 <= 1
	// Optimal: x = y = z = 1/sqrt(3), f* = -sqrt(3)

	auto objective = [](const VectorN<Real, 3>& x) -> Real {
		return -x[0] - x[1] - x[2];
	};

	auto sphere = [](const VectorN<Real, 3>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - REAL(1.0);
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-3);
	config.max_outer_iterations = 100;

	AugmentedLagrangian<3> alm(config);
	alm.AddInequalityConstraint(sphere);
	alm.SetBounds(-2.0, 2.0);

	VectorN<Real, 3> x0{ 0.1, 0.1, 0.1 };
	auto result = alm.Optimize(objective, x0);

	Real sphereVal = result.x[0] * result.x[0] + result.x[1] * result.x[1] + result.x[2] * result.x[2];

	REQUIRE(sphereVal <= REAL(1.0) + REAL(0.02));
	REQUIRE_THAT(result.objective_value, WithinAbs(-std::sqrt(REAL(3.0)), REAL(0.2)));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Mixed Equality and Inequality Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: Mixed equality and inequality constraints", "[ALM][optimization]")
{
	// min x^2 + y^2 subject to x + y = 1, x >= 0.3
	// Optimal: x = 0.5, y = 0.5 (equality is tighter)

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	auto equalityConstraint = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] + x[1] - REAL(1.0);
	};
	auto inequalityConstraint = [](const VectorN<Real, 2>& x) -> Real {
		return REAL(0.3) - x[0];  // x >= 0.3 means 0.3 - x <= 0
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 100;

	AugmentedLagrangian<2> alm(config);
	alm.AddEqualityConstraint(equalityConstraint);
	alm.AddInequalityConstraint(inequalityConstraint);
	alm.SetBounds(-10.0, 10.0);

	VectorN<Real, 2> x0{ 0.0, 0.0 };
	auto result = alm.Optimize(objective, x0);

	REQUIRE_THAT(result.x[0] + result.x[1], WithinAbs(REAL(1.0), REAL(0.01)));
	REQUIRE(result.x[0] >= REAL(0.3) - REAL(0.01));
	REQUIRE_THAT(result.x[0], WithinAbs(REAL(0.5), REAL(0.05)));
	REQUIRE_THAT(result.x[1], WithinAbs(REAL(0.5), REAL(0.05)));
}

///////////////////////////////////////////////////////////////////////////////////////////
// Result Structure Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: Result contains multipliers", "[ALM][result]")
{
	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	auto ineqConstraint = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] + x[1] - REAL(1.0);
	};
	auto eqConstraint = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] - x[1];
	};

	ALMConfig config;
	config.max_outer_iterations = 50;

	AugmentedLagrangian<2> alm(config);
	alm.AddInequalityConstraint(ineqConstraint);
	alm.AddEqualityConstraint(eqConstraint);
	alm.SetBounds(-10.0, 10.0);

	VectorN<Real, 2> x0{ 0.0, 0.0 };
	auto result = alm.Optimize(objective, x0);

	REQUIRE(result.inequality_multipliers.size() == 1);
	REQUIRE(result.equality_multipliers.size() == 1);
	REQUIRE(result.outer_iterations > 0);
	REQUIRE(result.inner_iterations > 0);
	REQUIRE(result.function_evaluations > 0);
}

TEST_CASE("ALM: Result follows the IterativeResultBase contract", "[ALM][result]")
{
	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	ALMConfig config;
	config.max_outer_iterations = 5;

	AugmentedLagrangian<2> alm(config);
	alm.SetBounds(-10.0, 10.0);

	// No constraints - converges immediately
	VectorN<Real, 2> x0{ 1.0, 1.0 };
	auto result = alm.Optimize(objective, x0);

	REQUIRE(result.converged);
	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.error_message.empty());               // empty on success (contract policy)
	REQUIRE(result.algorithm_name == "AugmentedLagrangian");
	REQUIRE(result.iterations_used == result.outer_iterations);
	REQUIRE_THAT(result.achieved_tolerance, WithinAbs(result.constraint_violation, REAL(1e-15)));
}

///////////////////////////////////////////////////////////////////////////////////////////
// SolveConstrained Convenience Function Tests
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: SolveConstrained convenience function", "[ALM][convenience]")
{
	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	std::vector<std::function<Real(const VectorN<Real, 2>&)>> ineqConstraints = {
		[](const VectorN<Real, 2>& x) -> Real { return -x[0] - x[1] + REAL(1.0); }  // x + y >= 1
	};

	VectorN<Real, 2> x0{ REAL(0.0), REAL(0.0) };
	auto result = SolveConstrained<2>(objective, x0, ineqConstraints);

	REQUIRE_THAT(result.x[0] + result.x[1], WithinAbs(REAL(1.0), REAL(0.1)));
	REQUIRE(result.outer_iterations > 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Edge Cases
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ALM: Unconstrained optimization (no constraints)", "[ALM][edge-cases]")
{
	// min (x-1)^2 + (y-2)^2 with no constraints; optimal (1, 2)

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return (x[0] - REAL(1.0)) * (x[0] - REAL(1.0)) + (x[1] - REAL(2.0)) * (x[1] - REAL(2.0));
	};

	ALMConfig config;
	config.max_outer_iterations = 10;
	config.max_inner_iterations = 5000;
	config.inner_simplex_delta = REAL(3.0);   // larger simplex to reach (1,2) from (0,0)

	AugmentedLagrangian<2> alm(config);

	VectorN<Real, 2> x0{ REAL(0.0), REAL(0.0) };
	auto result = alm.Optimize(objective, x0);

	REQUIRE(result.converged);
	REQUIRE_THAT(result.x[0], WithinAbs(REAL(1.0), REAL(0.1)));
	REQUIRE_THAT(result.x[1], WithinAbs(REAL(2.0), REAL(0.1)));
}

TEST_CASE("ALM: Infeasible starting point converges to feasible", "[ALM][edge-cases]")
{
	// min x^2 + y^2 subject to x >= 1

	auto objective = [](const VectorN<Real, 2>& x) -> Real {
		return x[0] * x[0] + x[1] * x[1];
	};

	auto constraint = [](const VectorN<Real, 2>& x) -> Real {
		return REAL(1.0) - x[0];  // x >= 1 means 1 - x <= 0
	};

	ALMConfig config;
	config.feasibility_tolerance = REAL(1e-4);
	config.max_outer_iterations = 100;

	AugmentedLagrangian<2> alm(config);
	alm.AddInequalityConstraint(constraint);
	alm.SetBounds(-10.0, 10.0);

	VectorN<Real, 2> x0{ -REAL(5.0), -REAL(5.0) };
	auto result = alm.Optimize(objective, x0);

	REQUIRE(result.converged);
	REQUIRE(result.x[0] >= REAL(1.0) - REAL(0.01));
	REQUIRE_THAT(result.x[1], WithinAbs(REAL(0.0), REAL(0.1)));
}

} // namespace MML::Tests::Algorithms::AugmentedLagrangianTests
