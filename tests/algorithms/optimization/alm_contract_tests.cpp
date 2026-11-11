#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h>
#endif

using namespace MML;
using namespace MML::Optimization;

TEST_CASE("ALM contracts: box bounds constrain evaluations and the returned minimizer", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_outer_iterations = 4;
	config.max_inner_iterations = 256;
	config.inner_simplex_delta = Real(3) / 10;
	AugmentedLagrangian<2> solver(config);
	solver.SetBounds(2, 3);
	int calls = 0;
	bool feasibleCalls = true;
	auto objective = [&](const VectorN<Real, 2>& x) {
		++calls;
		for (int i = 0; i < 2; ++i)
			feasibleCalls = feasibleCalls && x[i] >= 2 && x[i] <= 3;
		return x[0] * x[0] + x[1] * x[1] + 1;
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 2>{Real(5) / 2, Real(12) / 5});
	CHECK(feasibleCalls);
	for (int i = 0; i < 2; ++i) {
		CHECK(result.x[i] >= 2);
		CHECK(result.x[i] <= 3);
	}
	CHECK(std::abs(result.objective_value - 9) <= Precision::LooseTolerance);
	CHECK(result.function_evaluations == calls);
	CHECK(result.converged);
	CHECK(result.last_inner_status == AlgorithmStatus::Success);
	CHECK(result.stationarity_measure <= config.stationarity_tolerance);
}

TEST_CASE("ALM contracts: exhausted inner solve cannot imply outer stationarity", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_outer_iterations = 4;
	config.max_inner_iterations = 2;
	AugmentedLagrangian<1> solver(config);
	int calls = 0;
	auto objective = [&](const VectorN<Real, 1>& x) {
		++calls;
		return x[0] * x[0];
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 1>{1});
	CHECK_FALSE(result.converged);
	CHECK_FALSE(result.IsSuccess());
	CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(result.outer_iterations == 1);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.function_evaluations == calls);
	CHECK(result.function_evaluations > 0);
	CHECK(result.function_evaluations == 3);
	CHECK(result.inner_iterations == 0);
	CHECK(result.last_inner_status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(std::isnan(result.stationarity_measure));
}

TEST_CASE("ALM contracts: invalid configuration is rejected before work", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.initial_penalty = 0;
	CHECK_THROWS_AS(AugmentedLagrangian<1>{config}, MultidimOptimizationInputError);
	config = ALMConfig();
	config.stationarity_tolerance = std::numeric_limits<Real>::quiet_NaN();
	CHECK_THROWS_AS(AugmentedLagrangian<1>{config}, MultidimOptimizationInputError);
	config = ALMConfig();
	config.max_inner_iterations = 0;
	CHECK_THROWS_AS(AugmentedLagrangian<1>{config}, MultidimOptimizationInputError);
}

TEST_CASE("ALM contracts: invalid bound order is rejected by both setters", "[ALM][alm-contract]")
{
	AugmentedLagrangian<1> solver;
	CHECK_THROWS(solver.SetBounds(3, 2));
	CHECK_THROWS(solver.SetBounds(VectorN<Real, 1>{3}, VectorN<Real, 1>{2}));
}

TEST_CASE("ALM contracts: every configuration field is finite and ordered", "[ALM][alm-contract]")
{
	const ALMConfig valid;
	AugmentedLagrangian<1> solver(valid);
	const auto fields = {&ALMConfig::feasibility_tolerance, &ALMConfig::stationarity_tolerance,
		&ALMConfig::initial_penalty, &ALMConfig::penalty_multiplier, &ALMConfig::max_penalty,
		&ALMConfig::multiplier_lower_bound, &ALMConfig::multiplier_upper_bound,
		&ALMConfig::inner_tolerance, &ALMConfig::inner_simplex_delta};
	auto reject = [&](const ALMConfig& invalid) {
		CHECK_THROWS_AS(AugmentedLagrangian<1>{invalid}, MultidimOptimizationInputError);
		CHECK_THROWS_AS(solver.SetConfig(invalid), MultidimOptimizationInputError);
		for (const auto field : fields)
			CHECK(solver.GetConfig().*field == valid.*field);
		CHECK(solver.GetConfig().max_outer_iterations == valid.max_outer_iterations);
		CHECK(solver.GetConfig().max_inner_iterations == valid.max_inner_iterations);
	};
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		for (const auto field : fields) {
			auto config = valid;
			config.*field = invalid;
			reject(config);
		}
	}
	for (const auto field : {&ALMConfig::feasibility_tolerance, &ALMConfig::stationarity_tolerance,
		&ALMConfig::initial_penalty, &ALMConfig::max_penalty, &ALMConfig::inner_tolerance}) {
		for (const Real invalid : {Real(0), Real(-1)}) {
			auto config = valid;
			config.*field = invalid;
			reject(config);
		}
	}
	auto config = valid;
	config.inner_simplex_delta = 0;
	reject(config);
	config = valid;
	config.penalty_multiplier = Real(1) / 2;
	reject(config);
	config = valid;
	config.max_penalty = config.initial_penalty / 2;
	reject(config);
	config = valid;
	config.multiplier_lower_bound = 2;
	config.multiplier_upper_bound = 1;
	reject(config);
	config.multiplier_upper_bound = 3;
	reject(config);
	config = valid;
	config.multiplier_upper_bound = -1;
	reject(config);
	config = valid;
	config.max_outer_iterations = -1;
	reject(config);
	config = valid;
	config.max_inner_iterations = 0;
	reject(config);
	config.max_inner_iterations = -1;
	reject(config);
	config = valid;
	config.penalty_multiplier = 1;
	config.multiplier_lower_bound = 0;
	config.multiplier_upper_bound = 0;
	config.inner_simplex_delta = -1;
	CHECK_NOTHROW(solver.SetConfig(config));
}

TEST_CASE("ALM contracts: projected starts and every callback respect the box", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_inner_iterations = 512;
	config.inner_simplex_delta = Real(3) / 10;
	AugmentedLagrangian<2> solver(config);
	solver.SetBounds(VectorN<Real, 2>{2, 2}, VectorN<Real, 2>{3, 3});
	bool feasible = true;
	auto checkPoint = [&](const VectorN<Real, 2>& x) {
		feasible = feasible && x[0] >= 2 && x[0] <= 3 && x[1] >= 2 && x[1] <= 3;
	};
	solver.AddInequalityConstraint([&](const VectorN<Real, 2>& x) {
		checkPoint(x);
		return 2 - x[0];
	});
	solver.AddEqualityConstraint([&](const VectorN<Real, 2>& x) {
		checkPoint(x);
		return x[0] - x[1];
	});
	int calls = 0;
	VectorN<Real, 2> first;
	auto objective = [&](const VectorN<Real, 2>& x) {
		if (++calls == 1)
			first = x;
		checkPoint(x);
		return 1 + x[0] * x[0] + x[1] * x[1];
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 2>{100, -100});
	CHECK(feasible);
	CHECK(first[0] == 3);
	CHECK(first[1] == 2);
	REQUIRE(result.converged);
	CHECK(result.function_evaluations == calls);
	CHECK(std::abs(result.x[0] - 2) <= Precision::CoarseTolerance);
	CHECK(std::abs(result.x[1] - 2) <= Precision::CoarseTolerance);
	CHECK(std::abs(result.objective_value - 9) <= Precision::CoarseTolerance);
	CHECK(std::abs(result.x[0] - result.x[1]) <= config.feasibility_tolerance);
	CHECK(result.constraint_violation <= config.feasibility_tolerance);
	CHECK(result.stationarity_measure == result.step_norm / std::max({Real(1), std::abs(result.x[0]), std::abs(result.x[1])}));
}

TEST_CASE("ALM contracts: fixed and one-sided bounds have finite feasible evaluations", "[ALM][alm-contract]")
{
	SECTION("all coordinates fixed may converge with zero simplex iterations")
	{
		ALMConfig config;
		config.max_inner_iterations = 1;
		AugmentedLagrangian<2> solver(config);
		solver.SetBounds(2, 2);
		int calls = 0;
		auto objective = [&](const VectorN<Real, 2>& x) {
			++calls;
			CHECK(x[0] == 2);
			CHECK(x[1] == 2);
			return x[0] * x[0] + x[1] * x[1];
		};
		const auto result = solver.Optimize(objective, VectorN<Real, 2>{0, 100});
		CHECK(result.converged);
		CHECK(result.inner_iterations == 0);
		CHECK(result.outer_iterations == 2);
		CHECK(result.last_inner_status == AlgorithmStatus::Success);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 9);
	}
	SECTION("a fixed coordinate and an unbounded side")
	{
		AugmentedLagrangian<2> solver;
		solver.SetBounds(VectorN<Real, 2>{2, -std::numeric_limits<Real>::infinity()}, VectorN<Real, 2>{2, 3});
		bool feasible = true;
		auto objective = [&](const VectorN<Real, 2>& x) {
			feasible = feasible && x[0] == 2 && x[1] <= 3 && std::isfinite(x[1]);
			return x[0] * x[0] + (x[1] - 5) * (x[1] - 5) + 1;
		};
		const auto result = solver.Optimize(objective, VectorN<Real, 2>{-10, 100});
		CHECK(feasible);
		CHECK(result.converged);
		CHECK(result.x[0] == 2);
		CHECK(std::abs(result.x[1] - 3) <= Precision::CoarseTolerance);
	}
}

TEST_CASE("ALM contracts: failed bound setters retain previous valid bounds", "[ALM][alm-contract]")
{
	AugmentedLagrangian<1> solver;
	solver.SetBounds(2, 3);
	const Real inf = std::numeric_limits<Real>::infinity();
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	for (const auto& bounds : {std::pair<Real, Real>{Real(3), Real(2)}, {nan, Real(3)}, {Real(2), nan}, {inf, inf}, {-inf, -inf}})
		CHECK_THROWS_AS(solver.SetBounds(bounds.first, bounds.second), BoundConstraintsError);
	bool feasible = true;
	auto objective = [&](const VectorN<Real, 1>& x) {
		feasible = feasible && x[0] >= 2 && x[0] <= 3;
		return x[0] * x[0] + 1;
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 1>{10});
	CHECK(feasible);
	CHECK(result.converged);
	CHECK(result.x[0] == 2);
}

TEST_CASE("ALM contracts: zero outer budget returns assessed projected input without stationarity", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_outer_iterations = 0;
	AugmentedLagrangian<1> solver(config);
	solver.SetBounds(2, 3);
	solver.AddInequalityConstraint([](const VectorN<Real, 1>& x) { return x[0] - 2; });
	solver.AddEqualityConstraint([](const VectorN<Real, 1>& x) { return x[0] - 4; });
	int calls = 0;
	auto objective = [&](const VectorN<Real, 1>& x) {
		++calls;
		return x[0] * x[0] + 1;
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 1>{10});
	CHECK_FALSE(result.converged);
	CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(result.outer_iterations == 0);
	CHECK(result.iterations_used == 0);
	CHECK(result.inner_iterations == 0);
	CHECK_FALSE(result.last_inner_status.has_value());
	CHECK(std::isnan(result.step_norm));
	CHECK(std::isnan(result.stationarity_measure));
	CHECK(result.x[0] == 3);
	CHECK(result.constraint_violation == 1);
	CHECK(result.achieved_tolerance == 1);
	CHECK(result.objective_value == 10);
	CHECK(result.function_evaluations == 1);
	CHECK(calls == 1);
}

TEST_CASE("ALM contracts: partial inner minimizer is retained without false success", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_inner_iterations = 3;
	AugmentedLagrangian<1> solver(config);
	int calls = 0;
	auto objective = [&](const VectorN<Real, 1>& x) {
		++calls;
		return x[0] * x[0];
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 1>{1});
	CHECK_FALSE(result.converged);
	CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(result.outer_iterations == 1);
	CHECK(result.x[0] == 0);
	CHECK(result.objective_value == 0);
	CHECK(result.function_evaluations == calls);
	CHECK(calls == 5);
	CHECK(std::isnan(result.stationarity_measure));
}

TEST_CASE("ALM contracts: later inner failure invalidates the previous step proxy", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_inner_iterations = 2;
	AugmentedLagrangian<1> solver(config);
	solver.AddEqualityConstraint([](const VectorN<Real, 1>& x) { return x[0] + 1; });
	auto objective = [](const VectorN<Real, 1>& x) { return 1 - 5 * (x[0] + 1) * (x[0] + 1); };
	const auto result = solver.Optimize(objective, VectorN<Real, 1>{0});
	CHECK_FALSE(result.converged);
	CHECK(result.outer_iterations == 2);
	CHECK(result.last_inner_status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(std::isnan(result.stationarity_measure));
	CHECK(result.function_evaluations == 5);
	CHECK(result.error_message.find("inner solve failed") != std::string::npos);
}

TEST_CASE("ALM contracts: bounded inner failure is propagated with exact objective counts", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.max_inner_iterations = 2;
	AugmentedLagrangian<1> solver(config);
	solver.SetBounds(2, 3);
	int calls = 0;
	auto objective = [&](const VectorN<Real, 1>& x) {
		++calls;
		return x[0] * x[0] + 1;
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 1>{Real(5) / 2});
	CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK_FALSE(result.converged);
	CHECK(result.outer_iterations == 1);
	CHECK(result.x[0] == Real(5) / 2);
	CHECK(result.last_inner_status == AlgorithmStatus::MaxIterationsExceeded);
	CHECK(result.function_evaluations == calls);
	CHECK(calls == 4);
}

TEST_CASE("ALM contracts: nonfinite starts callbacks and augmentation are explicit errors", "[ALM][alm-contract]")
{
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		AugmentedLagrangian<1> solver;
		int calls = 0;
		auto objective = [&](const VectorN<Real, 1>& x) { ++calls; return x[0] * x[0]; };
		CHECK_THROWS_AS(solver.Optimize(objective, VectorN<Real, 1>{invalid}), MultidimOptimizationInputError);
		CHECK(calls == 0);
		auto nonfiniteObjective = [&](const VectorN<Real, 1>&) { return invalid; };
		CHECK_THROWS_AS(solver.Optimize(nonfiniteObjective, VectorN<Real, 1>{0}), MultidimOptimizationInputError);
		for (bool equality : {false, true}) {
			solver.ClearConstraints();
			auto constraint = [&](const VectorN<Real, 1>&) { return invalid; };
			if (equality)
				solver.AddEqualityConstraint(constraint);
			else
				solver.AddInequalityConstraint(constraint);
			CHECK_THROWS_AS(solver.Optimize(objective, VectorN<Real, 1>{0}), MultidimOptimizationInputError);
		}
	}
	AugmentedLagrangian<1> solver;
	CHECK_THROWS_AS(solver.AddInequalityConstraint({}), MultidimOptimizationInputError);
	CHECK_THROWS_AS(solver.AddEqualityConstraint({}), MultidimOptimizationInputError);
	solver.AddEqualityConstraint([](const VectorN<Real, 1>&) { return std::numeric_limits<Real>::max(); });
	auto objective = [](const VectorN<Real, 1>&) { return Real(0); };
	CHECK_THROWS_AS(solver.Optimize(objective, VectorN<Real, 1>{0}), MultidimOptimizationInputError);
}

TEST_CASE("ALM contracts: nonfinite constraint updates cannot be hidden by clamping", "[ALM][alm-contract]")
{
	for (bool equality : {false, true}) {
		ALMConfig config;
		config.max_inner_iterations = 1;
		AugmentedLagrangian<1> solver(config);
		solver.SetBounds(2, 2);
		int constraintCalls = 0;
		auto constraint = [&](const VectorN<Real, 1>&) {
			return ++constraintCalls == 6 ? std::numeric_limits<Real>::quiet_NaN() : Real(1);
		};
		if (equality)
			solver.AddEqualityConstraint(constraint);
		else
			solver.AddInequalityConstraint(constraint);
		auto objective = [](const VectorN<Real, 1>&) { return Real(0); };
		CHECK_THROWS_AS(solver.Optimize(objective, VectorN<Real, 1>{0}), MultidimOptimizationInputError);
		CHECK(constraintCalls == 6);
	}
}

TEST_CASE("ALM contracts: nonlinear constraints retain independently known optimum", "[ALM][alm-contract]")
{
	ALMConfig config;
	config.inner_tolerance = Precision::ZeroTolerance;
	AugmentedLagrangian<2> solver(config);
	solver.SetBounds(0, 1);
	solver.AddEqualityConstraint([](const VectorN<Real, 2>& x) { return x[0] + x[1] - 1; });
	solver.AddInequalityConstraint([](const VectorN<Real, 2>& x) { return Real(1) / 4 - x[0]; });
	int calls = 0;
	auto objective = [&](const VectorN<Real, 2>& x) {
		++calls;
		// x0^2+x1^2-1/2, factored to resolve the constrained minimum without an O(1) value floor.
		const Real dx = x[0] - Real(0.5), dy = x[1] - Real(0.5);
		return dx * dx + dy * dy + (x[0] + x[1] - 1);
	};
	const auto result = solver.Optimize(objective, VectorN<Real, 2>{2, 2});
	REQUIRE(result.converged);
	CHECK(std::abs(result.x[0] + result.x[1] - 1) <= config.feasibility_tolerance);
	CHECK(result.x[0] >= Real(1) / 4 - config.feasibility_tolerance);
	// The value-only quadratic minimum has a sqrt(epsilon) coordinate-resolution floor.
	const Real coordinateTolerance = std::max(Precision::CoarseTolerance, std::sqrt(std::numeric_limits<Real>::epsilon()));
	CHECK(std::abs(result.x[0] - Real(1) / 2) <= coordinateTolerance);
	CHECK(std::abs(result.x[1] - Real(1) / 2) <= coordinateTolerance);
	CHECK(std::abs(result.objective_value) <= Precision::LooseTolerance);
	CHECK(result.function_evaluations == calls);
	CHECK(result.last_inner_status == AlgorithmStatus::Success);
	CHECK(result.stationarity_measure <= config.stationarity_tolerance);
}

namespace {
	struct FixedObjective {
		Real operator()(const VectorN<Real, 1>& x) const { return x[0] * x[0]; }
	};
	struct DynamicObjective {
		Real operator()(Vector<Real>& x) const { return x[0] * x[0]; }
	};
	template<class F>
	concept AcceptedObjective = requires(AugmentedLagrangian<1>& solver, F&& function, const VectorN<Real, 1>& x) {
		solver.Optimize(std::forward<F>(function), x);
	};
	template<int N>
	concept AcceptedDimension = requires { typename AugmentedLagrangian<N>; };
}

TEST_CASE("ALM contracts: constrained forwarding accepts native and dynamic callable surfaces", "[ALM][alm-contract]")
{
	STATIC_REQUIRE(AcceptedObjective<FixedObjective>);
	STATIC_REQUIRE(AcceptedObjective<FixedObjective&>);
	STATIC_REQUIRE(AcceptedObjective<DynamicObjective>);
	STATIC_REQUIRE_FALSE(AcceptedObjective<int>);
	STATIC_REQUIRE_FALSE(AcceptedDimension<0>);
	STATIC_REQUIRE_FALSE(AcceptedDimension<-1>);
	ALMConfig config;
	config.max_inner_iterations = 2;
	AugmentedLagrangian<1> solver(config);
	const auto raw = solver.Optimize([](const VectorN<Real, 1>& x) { return x[0] * x[0]; }, VectorN<Real, 1>{1});
	const auto dynamic = solver.Optimize(DynamicObjective(), VectorN<Real, 1>{1});
	const auto convenience = SolveConstrained<1>(FixedObjective(), VectorN<Real, 1>{1}, {}, {}, config);
	int calls = 0;
	ScalarFunctionFromStdFunc<1> native([&](const VectorN<Real, 1>& x) { ++calls; return x[0] * x[0]; });
	IScalarFunction<1>& interface = native;
	const auto throughInterface = solver.Optimize(interface, VectorN<Real, 1>{1});
	CHECK(throughInterface.function_evaluations == calls);
	for (const auto& result : {raw, dynamic, convenience, throughInterface}) {
		CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK_FALSE(result.converged);
		CHECK(result.function_evaluations == 3);
		CHECK(result.x[0] == 1);
	}
}
