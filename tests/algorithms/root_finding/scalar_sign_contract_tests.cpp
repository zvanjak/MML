#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/RootFinding.h>
#endif

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace MML;
using namespace MML::RootFinding;
using Catch::Matchers::WithinAbs;

namespace {
	using DetailedSolver = RootFindingResult (*)(const IRealFunction&, Real, Real, const RootFindingConfig&);
	using PlainSolver = Real (*)(const IRealFunction&, Real, Real, Real);
	struct Method { const char* name; DetailedSolver detailed; PlainSolver plain; };
	const Method methods[] = {
		{"Bisection", FindRootBisectionDetailed, FindRootBisection},
		{"Newton", FindRootNewtonDetailed, FindRootNewton},
		{"Secant", FindRootSecantDetailed, FindRootSecant},
		{"FalsePosition", FindRootFalsePositionDetailed, FindRootFalsePosition},
		{"Ridders", FindRootRiddersDetailed, FindRootRidders},
		{"Brent", FindRootBrentDetailed, FindRootBrent}
	};

	std::vector<Real> Units() {
		const int extreme = std::numeric_limits<Real>::max_exponent * 3 / 4;
		return {std::ldexp(Real(1), -extreme), Real(1), std::ldexp(Real(1), extreme)};
	}

	RootFindingConfig Config(Real unit, int budget = 64) {
		RootFindingConfig config;
		config.tolerance = Precision::EqualityTolerance;
		config.x_tolerance = Precision::EqualityTolerance;
		config.f_tolerance = unit * Precision::EqualityTolerance;
		config.max_iterations = budget;
		return config;
	}
}

TEST_CASE("Brent accepts the reported tiny affine bracket within its requested residual and budget",
	"[RootFinding][scalar-sign-contract]") {
	if constexpr (std::numeric_limits<Real>::max_exponent < 700) {
		SKIP("Reported binary64 witness requires a wider exponent range");
	} else {
		const Real unit = Real(1e-200);
		RootFindingConfig config;
		config.tolerance = Real(1e-210);
		config.f_tolerance = Real(1e-210);
		config.max_iterations = 16;
		int calls = 0;
		RealFunctionFromStdFunc function([&](Real x) { ++calls; return unit * (x - Real(0.25)); });
		const auto result = FindRootBrentDetailed(function, 0, 1, config);
		CHECK(result.IsSuccess());
		CHECK_THAT(result.root, WithinAbs(Real(0.25), config.f_tolerance / unit));
		CHECK(std::abs(result.function_value) <= config.f_tolerance);
		CHECK(result.iterations_used <= 16);
		CHECK(result.function_evaluations == calls);
		CHECK_FALSE(result.status == AlgorithmStatus::InvalidInput);
	}
}

TEST_CASE("All scalar refiners retain analytic roots under finite function-unit rescaling",
	"[RootFinding][scalar-sign-contract]") {
	for (const auto& method : methods)
		for (Real unit : Units())
			for (Real direction : {Real(-1), Real(1)})
				for (bool reversed : {false, true})
					for (bool quadratic : {false, true}) {
						CAPTURE(method.name, unit, direction, reversed, quadratic);
						int calls = 0;
						auto normalized = [=](Real x) { return direction * (quadratic ? x * x - 2 : x - Real(0.25)); };
						RealFunctionFromStdFunc function([&](Real x) { ++calls; return unit * normalized(x); });
						const auto config = Config(unit);
						const Real upper = quadratic ? Real(2) : Real(1);
						const auto result = method.detailed(function, reversed ? upper : Real(0), reversed ? Real(0) : upper, config);
						CHECK(result.IsSuccess());
						CHECK(result.algorithm_name == method.name);
						CHECK_THAT(result.root, WithinAbs(quadratic ? std::sqrt(Real(2)) : Real(0.25),
							Precision::DefaultToleranceStrict));
						CHECK(std::isfinite(result.root));
						CHECK(std::isfinite(result.function_value));
						CHECK(result.function_value == unit * normalized(result.root));
						CHECK(std::abs(result.function_value / unit) <= Precision::DefaultToleranceStrict * 4);
						CHECK(result.iterations_used <= config.max_iterations);
						CHECK(result.function_evaluations == calls);
						CHECK((std::abs(result.function_value) <= config.f_tolerance ||
							(std::isfinite(result.x_error) && result.x_error <= EffectiveXTolerance(config, result.root))));
					}
	for (const auto& method : methods)
		for (Real unit : Units()) {
			RealFunctionFromStdFunc function([=](Real x) { return unit * x; });
			const Real tolerance = std::min(Precision::EqualityTolerance, unit * Precision::EqualityTolerance);
			CHECK_THAT(method.plain(function, -1, 1, tolerance), WithinAbs(Real(0), Precision::DefaultToleranceStrict));
		}
}

TEST_CASE("Same-sign finite brackets are rejected before refinement even when their product underflows",
	"[RootFinding][scalar-sign-contract]") {
	for (const auto& method : methods)
		for (Real unit : Units())
			for (Real direction : {Real(-1), Real(1)}) {
				CAPTURE(method.name, unit, direction);
				int calls = 0;
				RealFunctionFromStdFunc function([&](Real x) { ++calls; return direction * unit * (1 + x * x); });
				const auto result = method.detailed(function, 0, 1, Config(unit));
				CHECK_FALSE(result.IsSuccess());
				CHECK(result.status == AlgorithmStatus::InvalidInput);
				CHECK_FALSE(result.error_message.empty());
				CHECK(result.iterations_used == 0);
				CHECK(result.function_evaluations == 2);
				CHECK(calls == 2);
			}
}

TEST_CASE("Secant and false-position interpolate opposite maximum finite values without overflow",
	"[RootFinding][scalar-sign-contract]") {
	const Real unit = std::numeric_limits<Real>::max();
	for (const auto& method : methods)
		for (Real direction : {Real(-1), Real(1)})
			for (bool reversed : {false, true}) {
				CAPTURE(method.name, direction, reversed);
				RealFunctionFromStdFunc function([=](Real x) { return direction * unit * (2 * x - 1); });
				auto config = Config(unit, 16);
				if (std::string(method.name) == "Secant" || std::string(method.name) == "FalsePosition")
					config.max_iterations = 1;
				const auto result = method.detailed(function, reversed ? Real(1) : Real(0), reversed ? Real(0) : Real(1), config);
				CHECK(result.IsSuccess());
				CHECK(result.root == Real(0.5));
				CHECK(result.function_value == 0);
				CHECK(result.iterations_used <= config.max_iterations);
				CHECK_THAT(method.plain(function, 0, 1, Precision::EqualityTolerance),
					WithinAbs(Real(0.5), Precision::DefaultToleranceStrict));
			}
}

TEST_CASE("Ridders retains one-step affine interpolation at tiny and huge function units",
	"[RootFinding][scalar-sign-contract]") {
	auto units = Units();
	units.push_back(std::numeric_limits<Real>::max());
	for (Real unit : units)
		for (Real direction : {Real(-1), Real(1)})
			for (bool reversed : {false, true}) {
				CAPTURE(unit, direction, reversed);
				RealFunctionFromStdFunc function([=](Real x) { return direction * unit * (x - Real(0.25)); });
				const auto result = FindRootRiddersDetailed(function, reversed ? Real(1) : Real(0), reversed ? Real(0) : Real(1), Config(unit, 1));
				CHECK(result.IsSuccess());
				CHECK_THAT(result.root, WithinAbs(Real(0.25), Precision::DefaultToleranceStrict));
				CHECK(result.iterations_used == 1);
				CHECK(std::abs(result.function_value / unit) <= Precision::DefaultToleranceStrict);
			}
}

TEST_CASE("Newton safeguarding cannot expand a bracket when an imperfect derivative suggests a backward step",
	"[RootFinding][scalar-sign-contract]") {
	for (Real unit : Units())
		for (Real direction : {Real(-1), Real(1)}) {
			CAPTURE(unit, direction);
			std::vector<Real> probes;
			RealFunctionFromStdFunc function([&](Real x) { probes.push_back(x); return direction * unit * (x - Real(0.875)); });
			RealFunctionFromStdFunc derivative([=](Real x) { return direction * unit * (x == Real(0.5) ? 3 : -8); });
			const auto result = FindRootNewtonDetailed(function, derivative, 0, 1, Config(unit, 2));
			CHECK_FALSE(result.IsSuccess());
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK(result.root == Real(0.8125));
			CHECK(probes == std::vector<Real>{0, 1, Real(0.5), Real(0.625), Real(0.8125)});
			CHECK(result.function_evaluations == 5);
			CHECK(result.iterations_used == 2);
		}
}

TEST_CASE("Endpoint tolerance and signed-zero policies precede scalar bracket sign validation",
	"[RootFinding][scalar-sign-contract]") {
	for (const auto& method : methods)
		for (Real unit : Units())
			for (Real direction : {Real(-1), Real(1)})
				for (bool nearEndpoint : {false, true}) {
					const auto config = Config(unit);
					const Real offset = nearEndpoint ? Precision::EqualityTolerance / 2 : Real(0);
					RealFunctionFromStdFunc function([=](Real x) { return direction * unit * (x + offset); });
					const auto result = method.detailed(function, 0, 1, config);
					CHECK(result.IsSuccess());
					CHECK(result.root == 0);
					CHECK(result.function_value == direction * unit * offset);
					CHECK(std::signbit(result.function_value) == (direction < 0));
					CHECK(std::isnan(result.x_error));
					CHECK(result.achieved_tolerance == std::abs(result.function_value));
					CHECK(result.function_evaluations == 2);
				}
}

TEST_CASE("Bracket expansion keeps product-free finite endpoint sign and exact-zero semantics",
	"[RootFinding][scalar-sign-contract]") {
	auto units = Units();
	units.push_back(std::numeric_limits<Real>::denorm_min() * 4);
	units.push_back(std::numeric_limits<Real>::max());
	for (Real unit : units) {
		Real left = 0, right = 1;
		RealFunctionFromStdFunc opposite([=](Real x) { return unit * (2 * x - 1); });
		CHECK(BracketRoot(opposite, left, right, 0));
		CHECK(left == 0);
		CHECK(right == 1);
		RealFunctionFromStdFunc positive([=](Real) { return unit; });
		CHECK_FALSE(BracketRoot(positive, left, right, 0));
		RealFunctionFromStdFunc endpoint([=](Real x) { return -unit * x; });
		CHECK(BracketRoot(endpoint, left, right, 0));
	}
	for (Real unit : Units()) {
		Real left = 0, right = 1;
		int calls = 0;
		RealFunctionFromStdFunc function([&](Real x) { ++calls; return unit * (x - Real(2.5)); });
		CHECK(BracketRoot(function, left, right, 1));
		CHECK(left < Real(2.5));
		CHECK(right > Real(2.5));
		CHECK(calls == 3);
	}
}

TEST_CASE("Ridders handles highly unequal finite endpoint magnitudes without losing their geometric mean",
	"[RootFinding][scalar-sign-contract]") {
	const Real large = std::numeric_limits<Real>::max();
	const Real small = std::numeric_limits<Real>::min();
	for (bool reversed : {false, true}) {
		RealFunctionFromStdFunc function([=](Real x) {
			return x < Real(0.25) ? -large * ((Real(0.25) - x) * 4) : small * ((x - Real(0.25)) / Real(0.75));
		});
		auto config = Config(1, 16);
		config.f_tolerance = std::numeric_limits<Real>::denorm_min();
		const auto result = FindRootRiddersDetailed(function, reversed ? Real(1) : Real(0), reversed ? Real(0) : Real(1), config);
		CHECK(result.IsSuccess());
		CHECK(result.root == Real(0.25));
		CHECK(result.function_value == 0);
		CHECK(result.iterations_used == 2);
	}
}

TEST_CASE("Nonfinite scalar endpoints retain numerical-failure classification before sign validation",
	"[RootFinding][scalar-sign-contract]") {
	for (const auto& method : methods)
		for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
			-std::numeric_limits<Real>::infinity()}) {
			RealFunctionFromStdFunc function([=](Real x) { return x == 0 ? invalid : Real(0); });
			const auto result = method.detailed(function, 0, 1, Config(1));
			CHECK_FALSE(result.IsSuccess());
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			CHECK_FALSE(result.error_message.empty());
			CHECK(result.function_evaluations == 2);
		}

}

TEST_CASE("Subnormal scalar brackets retain sign classification at a representable residual tolerance",
	"[RootFinding][scalar-sign-contract]") {
	const Real unit = std::numeric_limits<Real>::denorm_min() * 64;
	auto config = Config(1);
	config.f_tolerance = std::numeric_limits<Real>::denorm_min();
	for (const auto& method : methods)
		for (Real direction : {Real(-1), Real(1)}) {
			CAPTURE(method.name, direction);
			RealFunctionFromStdFunc function([=](Real x) { return direction * unit * (x - Real(0.25)); });
			const auto result = method.detailed(function, 0, 1, config);
			CHECK(result.IsSuccess());
			CHECK_THAT(result.root, WithinAbs(Real(0.25), (config.f_tolerance * 2) / unit));
			CHECK(std::abs(result.function_value) <= config.f_tolerance);
			CHECK(result.iterations_used <= config.max_iterations);
			CHECK(std::isfinite(result.root));
			RealFunctionFromStdFunc positive([=](Real x) { return direction * unit * (1 + x * x); });
			const auto rejected = method.detailed(positive, 0, 1, config);
			CHECK_FALSE(rejected.IsSuccess());
			CHECK(rejected.status == AlgorithmStatus::InvalidInput);
			CHECK_FALSE(rejected.error_message.empty());
			CHECK(rejected.function_evaluations == 2);
		}
}

TEST_CASE("Scaled secant interpolation still reports an exact zero denominator as numerical failure",
	"[RootFinding][scalar-sign-contract]") {
	RealFunction function([](Real x) { return x == 0 ? Real(-1) : x == 1 ? Real(2) : Real(1); });
	const auto result = FindRootSecantDetailed(function, 0, 1, Config(1));
	CHECK_FALSE(result.IsSuccess());
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.iterations_used == 3);
	CHECK(result.function_evaluations == 4);
}
