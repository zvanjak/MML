#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/RootFinding.h>
#include <mml/base/Function.h>
#include <mml/base/ComplexFunction.h>
#endif

#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Algorithms::RootConvergenceContracts
{
	TEST_CASE("Brent rejects nonpositive fallback iteration budgets before evaluating",
		"[RootFinding][root-convergence-contract]")
	{
		struct RestoreBudget {
			int value = Defaults::BrentMaxSteps;
			~RestoreBudget() { Defaults::BrentMaxSteps = value; }
		} restore;
		for (int budget : {0, -1}) {
			Defaults::BrentMaxSteps = budget;
			int calls = 0;
			auto function = [&](Real x) { ++calls; return x - Real(0.25); };
			RootFinding::RootFindingConfig config;
			config.max_iterations = 0;
			const auto rejected = RootFinding::FindRootBrentDetailed(function, Real(0), Real(1), config);
			CHECK(rejected.status == AlgorithmStatus::InvalidInput);
			CHECK_FALSE(rejected.IsSuccess());
			CHECK_FALSE(rejected.error_message.empty());
			CHECK(rejected.iterations_used == 0);
			CHECK(rejected.function_evaluations == 0);
			CHECK(calls == 0);
			config.max_iterations = 1;
			const auto explicitBudget = RootFinding::FindRootBrentDetailed(function, Real(0), Real(1), config);
			CHECK(explicitBudget.status != AlgorithmStatus::InvalidInput);
			CHECK(calls > 0);
		}
	}

	TEST_CASE("Bisection returns the residual-verified midpoint for either sign and endpoint order",
		"[RootFinding][root-convergence-contract]")
	{
		for (Real offset : {Real(-1e-7), Real(1e-7)}) {
			for (Real sign : {Real(-1), Real(1)}) {
				for (bool reversed : {false, true}) {
					CAPTURE(offset, sign, reversed);
					int calls = 0;
					auto function = [&](Real x) {
						++calls;
						return sign * (x - (Real(0.5) + offset));
					};
					RootFinding::RootFindingConfig config;
					config.tolerance = Real(1e-6);
					config.max_iterations = 16;
					auto result = RootFinding::FindRootBisectionDetailed(
						function, reversed ? Real(1) : Real(0), reversed ? Real(0) : Real(1), config);
					CHECK(result.IsSuccess());
					CHECK(result.root == Real(0.5));
					CHECK(result.function_value == sign * (result.root - (Real(0.5) + offset)));
					CHECK(std::abs(result.function_value) <= config.tolerance);
					CHECK(result.iterations_used == 1);
					CHECK(result.function_evaluations == calls);
					CHECK(calls == 3);
					CHECK(result.x_error == Real(0.5));
				}
			}
		}
	}

	TEST_CASE("Bisection exhaustion and verbose output describe the last evaluated midpoint",
		"[RootFinding][root-convergence-contract]")
	{
		int calls = 0;
		auto function = [&](Real x) { ++calls; return x - Real(0.25); };
		std::ostringstream output;
		RootFinding::RootFindingConfig config;
		config.tolerance = TOL(1e-12, 1e-6);
		config.max_iterations = 1;
		config.verbose = true;
		config.verboseStream = &output;
		auto result = RootFinding::FindRootBisectionDetailed(function, Real(0), Real(1), config);
		CHECK_FALSE(result.IsSuccess());
		CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(result.root == Real(0.5));
		CHECK(result.function_value == result.root - Real(0.25));
		CHECK(result.iterations_used == 1);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 3);
		CHECK(result.x_error == Real(0.5));
		CHECK(output.str().find("x=0.5, f(x)=0.25") != std::string::npos);
	}

	TEST_CASE("Bisection nonfinite midpoint failure identifies the point that was evaluated",
		"[RootFinding][root-convergence-contract]")
	{
		int calls = 0;
		auto function = [&](Real x) {
			++calls;
			return x == Real(0.5) ? std::numeric_limits<Real>::quiet_NaN() : x - Real(0.25);
		};
		auto result = RootFinding::FindRootBisectionDetailed(
			function, Real(0), Real(1), RootFinding::RootFindingConfig{});
		CHECK_FALSE(result.IsSuccess());
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(result.root == Real(0.5));
		CHECK(std::isnan(result.function_value));
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 3);
		CHECK_FALSE(result.error_message.empty());
	}

	TEST_CASE("Safeguarded Newton accepts exact midpoint roots before evaluating a derivative",
		"[RootFinding][root-convergence-contract]")
	{
		int calls = 0;
		int derivativeCalls = 0;
		RealFunctionFromStdFunc function([&](Real x) { ++calls; return x * x * x; });
		RealFunctionFromStdFunc derivative([&](Real x) { ++derivativeCalls; return Real(3) * x * x; });
		RootFinding::RootFindingConfig config;
		config.max_iterations = 1;
		auto result = RootFinding::FindRootNewtonDetailed(function, derivative, Real(-1), Real(1), config);
		CHECK(result.IsSuccess());
		CHECK(result.root == Real(0));
		CHECK(result.function_value == Real(0));
		CHECK(result.iterations_used == 0);
		CHECK(result.achieved_tolerance == Real(0));
		CHECK(std::isnan(result.x_error));
		CHECK(derivativeCalls == 0);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 3);

		calls = 0;
		result = RootFinding::FindRootNewtonDetailed(function, Real(-1), Real(1), config);
		CHECK(result.IsSuccess());
		CHECK(result.iterations_used == 0);
		CHECK(result.function_evaluations == calls);
		CHECK(calls == 3);
		CHECK(RootFinding::FindRootNewton(function, derivative, Real(-1), Real(1), config.tolerance) == Real(0));
	}

	TEST_CASE("Safeguarded Newton bisects for zero and nonfinite derivatives",
		"[RootFinding][root-convergence-contract]")
	{
		for (Real df : {Real(0), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
			CAPTURE(df);
			int calls = 0;
			RealFunctionFromStdFunc function([&](Real x) { ++calls; return x - Real(0.25); });
			RealFunctionFromStdFunc derivative([df](Real) { return df; });
			RootFinding::RootFindingConfig config;
			config.max_iterations = 3;
			auto result = RootFinding::FindRootNewtonDetailed(function, derivative, Real(0), Real(1), config);
			CHECK(result.IsSuccess());
			CHECK(result.root == Real(0.25));
			CHECK(result.function_value == Real(0));
			CHECK(result.function_evaluations == calls);
		}
	}

	TEST_CASE("Safeguarded Newton validates a new residual before accepting a small step",
		"[RootFinding][root-convergence-contract]")
	{
		const Real step = std::max(Real(1e-8), Precision::Eps);
		RootFinding::RootFindingConfig config;
		config.x_tolerance = Real(1e-6);
		config.max_iterations = 1;
		REQUIRE(Real(1) - step < Real(1));
		REQUIRE(step < config.x_tolerance);
		for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
			CAPTURE(invalid, step);
			int calls = 0;
			int derivativeCalls = 0;
			RealFunctionFromStdFunc function([&](Real x) {
				++calls;
				if (x == Real(0)) return Real(-1);
				if (x == Real(1) || x == Real(2)) return Real(1);
				return invalid;
			});
			RealFunctionFromStdFunc derivative([&](Real) { ++derivativeCalls; return Real(1) / step; });
			auto result = RootFinding::FindRootNewtonDetailed(function, derivative, Real(0), Real(2), config);
			CHECK_FALSE(result.IsSuccess());
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			CHECK(result.root == Real(1) - step);
			CHECK_FALSE(std::isfinite(result.function_value));
			if (std::isnan(invalid))
				CHECK(std::isnan(result.function_value));
			else
				CHECK(result.function_value == invalid);
			CHECK(result.x_error < config.x_tolerance);
			CHECK(result.iterations_used == 1);
			CHECK(result.function_evaluations == calls);
			CHECK(calls == 4);
			CHECK(derivativeCalls == 1);
			CHECK_FALSE(result.error_message.empty());
		}
	}

	TEST_CASE("Safeguarded Newton counts numerical samples and does not differentiate after exhaustion",
		"[RootFinding][root-convergence-contract]")
	{
		int calls = 0;
		auto function = [&](Real x) { ++calls; return x * x - Real(2); };
		RootFinding::RootFindingConfig config;
		config.tolerance = TOL(1e-12, 1e-6);
		config.max_iterations = 1;
		auto numeric = RootFinding::FindRootNewtonDetailed(function, Real(1), Real(2), config);
		CHECK_FALSE(numeric.IsSuccess());
		CHECK(numeric.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(numeric.function_evaluations == calls);
		CHECK(calls == 8);

		int derivativeCalls = 0;
		RealFunction realFunction([](Real x) { return x * x - Real(2); });
		RealFunctionFromStdFunc derivative([&](Real x) { ++derivativeCalls; return Real(2) * x; });
		auto analytic = RootFinding::FindRootNewtonDetailed(realFunction, derivative, Real(1), Real(2), config);
		CHECK_FALSE(analytic.IsSuccess());
		CHECK(analytic.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(derivativeCalls == 1);
		CHECK(analytic.function_value == analytic.root * analytic.root - Real(2));
	}

	TEST_CASE("Scalar endpoint acceptance and bracket preconditions remain intact",
		"[RootFinding][root-convergence-contract]")
	{
		RealFunction function([](Real x) { return x - Real(1); });
		RealFunction derivative([](Real) { return Real(1); });
		RootFinding::RootFindingConfig config;
		for (bool reversed : {false, true}) {
			const Real a = reversed ? Real(2) : Real(1);
			const Real b = reversed ? Real(1) : Real(2);
			for (const auto& result : {
				RootFinding::FindRootBisectionDetailed(function, a, b, config),
				RootFinding::FindRootNewtonDetailed(function, derivative, a, b, config)}) {
				CHECK(result.IsSuccess());
				CHECK(result.root == Real(1));
				CHECK(result.iterations_used == 0);
				CHECK(result.function_evaluations == 2);
			}
		}
		CHECK(RootFinding::FindRootNewtonDetailed(function, derivative, Real(2), Real(3), config).status
			== AlgorithmStatus::InvalidInput);
		CHECK(RootFinding::FindRootSecantDetailed(function, Real(2), Real(3), config).status
			== AlgorithmStatus::InvalidInput);
		CHECK(RootFinding::FindRootBisectionDetailed(function, Real(1), Real(1), config).status
			== AlgorithmStatus::InvalidInput);
	}

	TEST_CASE("Safeguarded Newton retains valid tiny derivatives and the scalar step-or-residual rule",
		"[RootFinding][root-convergence-contract]")
	{
		const Real scale = std::sqrt(std::numeric_limits<Real>::min());
		RealFunctionFromStdFunc function([scale](Real x) { return scale * (x - Real(0.25)); });
		RealFunctionFromStdFunc derivative([scale](Real) { return scale; });
		RootFinding::RootFindingConfig config;
		config.tolerance = scale * Precision::EqualityTolerance;
		config.max_iterations = 1;
		auto tiny = RootFinding::FindRootNewtonDetailed(function, derivative, Real(0), Real(1), config);
		CHECK(tiny.IsSuccess());
		CHECK(tiny.root == Real(0.25));
		CHECK(tiny.iterations_used == 1);

		RealFunction linear([](Real x) { return x - Real(0.25); });
		RealFunction steepDerivative([](Real) { return Real(1e8); });
		config = RootFinding::RootFindingConfig{};
		config.x_tolerance = Real(1e-6);
		config.max_iterations = 1;
		auto smallStep = RootFinding::FindRootNewtonDetailed(linear, steepDerivative, Real(0), Real(1), config);
		CHECK(smallStep.IsSuccess());
		CHECK(smallStep.x_error <= config.x_tolerance);
		CHECK(std::abs(smallStep.function_value) > RootFinding::EffectiveFTolerance(config));
		CHECK(smallStep.function_value == linear(smallStep.root));
	}

	TEST_CASE("Reserved scalar root configuration fields do not change solver behavior",
		"[RootFinding][root-convergence-contract]")
	{
		RealFunction function([](Real x) { return x - Real(0.25); });
		RealFunction derivative([](Real) { return Real(0); });
		RootFinding::RootFindingConfig config;
		auto baseline = RootFinding::FindRootNewtonDetailed(function, derivative, Real(0), Real(1), config);
		config.initial_step = Real(12);
		config.max_stagnant_iterations = 1;
		auto reserved = RootFinding::FindRootNewtonDetailed(function, derivative, Real(0), Real(1), config);
		CHECK(reserved.IsSuccess());
		CHECK(reserved.root == baseline.root);
		CHECK(reserved.iterations_used == baseline.iterations_used);
		CHECK(reserved.function_evaluations == baseline.function_evaluations);
	}

	TEST_CASE("Complex Newton accepts a root reached on its final permitted update",
		"[RootFinding][root-convergence-contract]")
	{
		const Real tolerance = TOL(1e-10, Precision::EqualityTolerance);
		int calls = 0;
		Complex lastEvaluated;
		ComplexFunctionFromStdFunc function([&](Complex z) {
			++calls;
			lastEvaluated = z;
			return z - Complex(1);
		});
		const auto noUpdate = RootFinding::FindRootNewtonComplex(function, Complex(0), tolerance, 0);
		REQUIRE_FALSE(noUpdate.converged);
		REQUIRE(noUpdate.status == AlgorithmStatus::MaxIterationsExceeded);
		REQUIRE(std::abs(noUpdate.function_value) > tolerance);
		REQUIRE(calls == 1);
		calls = 0;
		auto result = RootFinding::FindRootNewtonComplex(function, Complex(0), tolerance, 1);
		CHECK(result.converged);
		CHECK(result.status == AlgorithmStatus::Success);
		CHECK(result.iterations_used == 1);
		CHECK(calls == 6);
		CHECK(result.root == lastEvaluated);
		CHECK(std::abs(result.root - Complex(1)) < tolerance);
		CHECK(result.function_value == function(result.root));
		CHECK(result.achieved_tolerance == std::abs(result.function_value));
		CHECK(result.error_message.empty());
	}

	TEST_CASE("Muller accepts an exact quadratic root on its final permitted update",
		"[RootFinding][root-convergence-contract]")
	{
		ComplexFunctionFromStdFunc function([](Complex z) { return z * z + Complex(1); });
		const Real tolerance = TOL(Real(1e-10), Precision::EqualityTolerance);
		for (bool convenience : {false, true}) {
			auto result = convenience
				? RootFinding::FindRootMuller(function, Complex(0.5), tolerance, 1)
				: RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(1), tolerance, 1);
			CHECK(result.converged);
			CHECK(result.status == AlgorithmStatus::Success);
			CHECK(result.iterations_used == 1);
			CHECK(std::abs(result.root.real()) < tolerance);
			CHECK(std::abs(std::abs(result.root.imag()) - Real(1)) < tolerance);
			CHECK(result.function_value == function(result.root));
			CHECK(result.achieved_tolerance == std::abs(result.function_value));
			CHECK(result.error_message.empty());
		}
	}

	TEST_CASE("Complex solvers preserve honest one-update exhaustion",
		"[RootFinding][root-convergence-contract]")
	{
		ComplexFunctionFromStdFunc function([](Complex z) { return z * z * z - Complex(2); });
		for (const auto& result : {
			RootFinding::FindRootNewtonComplex(function, Complex(1), Real(1e-12), 1),
			RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(1), Real(1e-12), 1)}) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK(result.iterations_used == 1);
			CHECK(result.function_value == function(result.root));
			CHECK(result.achieved_tolerance == std::abs(result.function_value));
			CHECK(std::abs(result.function_value) > Real(1e-12));
			CHECK_FALSE(result.error_message.empty());
		}
	}

	TEST_CASE("Complex solvers check initial roots even with no update budget",
		"[RootFinding][root-convergence-contract]")
	{
		ComplexFunctionFromStdFunc function([](Complex z) { return z - Complex(1); });
		for (const auto& result : {
			RootFinding::FindRootNewtonComplex(function, Complex(1), Real(1e-10), 0),
			RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(1), Real(1e-10), 0)}) {
			CHECK(result.converged);
			CHECK(result.status == AlgorithmStatus::Success);
			CHECK(result.root == Complex(1));
			CHECK(result.function_value == Complex(0));
			CHECK(result.iterations_used == 0);
		}
		for (const auto& result : {
			RootFinding::FindRootNewtonComplex(function, Complex(0), Real(1e-10), 0),
			RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(2), Real(1e-10), 0)}) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK(result.iterations_used == 0);
			CHECK(result.function_value == function(result.root));
		}
	}

	TEST_CASE("Complex Newton rejects nonfinite steps while retaining its last evaluated point",
		"[RootFinding][root-convergence-contract]")
	{
		ComplexFunctionFromStdFunc function([](Complex z) {
			return z == Complex(0) ? Complex(std::numeric_limits<Real>::max()) : Real(1e-7) * z;
		});
		auto result = RootFinding::FindRootNewtonComplex(function, Complex(0), Real(1e-10), 1);
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(result.root == Complex(0));
		CHECK(result.function_value == function(result.root));
		CHECK(result.achieved_tolerance == std::abs(result.function_value));
		CHECK(result.iterations_used == 1);
	}

	TEST_CASE("Complex solvers reject a nonfinite residual reached by their final update",
		"[RootFinding][root-convergence-contract]")
	{
		RootFinding::ComplexRootFindingResult result;
		SECTION("Newton") {
			std::vector<Complex> samples;
			ComplexFunctionFromStdFunc function([&](Complex z) {
				samples.push_back(z);
				return z.real() < Real(0.5) ? Complex(std::numeric_limits<Real>::quiet_NaN())
					: z * z + Complex(1);
			});
			result = RootFinding::FindRootNewtonComplex(function, Complex(1, 1), Real(1e-10), 1);
			REQUIRE(samples.size() == 6);
			for (std::size_t i = 0; i + 1 < samples.size(); ++i)
				CHECK(samples[i].real() >= Real(0.5));
			CHECK(samples.back().real() < Real(0.5));
			CHECK(result.root == samples.back());
		}
		SECTION("Muller") {
			ComplexFunctionFromStdFunc function([](Complex z) {
				return std::abs(z.imag()) > Real(0.5) ? Complex(std::numeric_limits<Real>::quiet_NaN())
					: z * z + Complex(1);
			});
			result = RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(1), Real(1e-10), 1);
		}
		CHECK_FALSE(result.converged);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK(std::isnan(result.function_value.real()));
		CHECK(std::isnan(result.achieved_tolerance));
		CHECK(result.iterations_used == 1);
		CHECK_FALSE(result.error_message.empty());
	}

	TEST_CASE("Muller reports degenerate interpolation without corrupting its last point",
		"[RootFinding][root-convergence-contract]")
	{
		ComplexFunctionFromStdFunc function([](Complex) { return Complex(1); });
		for (const auto& result : {
			RootFinding::FindRootMuller(function, Complex(0), Complex(0), Complex(1)),
			RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(1))}) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::Stalled);
			CHECK(result.root == Complex(1));
			CHECK(result.function_value == function(result.root));
			CHECK_FALSE(result.error_message.empty());
		}
	}

	TEST_CASE("Complex solvers reject invalid budgets and nonfinite residuals explicitly",
		"[RootFinding][root-convergence-contract]")
	{
		ComplexFunctionFromStdFunc function([](Complex z) { return z - Complex(1); });
		for (const auto& result : {
			RootFinding::FindRootNewtonComplex(function, Complex(0), Real(1e-10), -1),
			RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(2), Real(1e-10), -1)}) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::InvalidInput);
			CHECK(result.iterations_used == 0);
			CHECK_FALSE(result.error_message.empty());
		}
		for (Real tolerance : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity()}) {
			for (const auto& result : {
				RootFinding::FindRootNewtonComplex(function, Complex(0), tolerance, 1),
				RootFinding::FindRootMuller(function, Complex(0), Complex(0.5), Complex(2), tolerance, 1)}) {
				CHECK_FALSE(result.converged);
				CHECK(result.status == AlgorithmStatus::InvalidInput);
				CHECK(result.iterations_used == 0);
			}
		}
		ComplexFunctionFromStdFunc nonfinite([](Complex) {
			return Complex(std::numeric_limits<Real>::quiet_NaN());
		});
		for (const auto& result : {
			RootFinding::FindRootNewtonComplex(nonfinite, Complex(0), Real(1e-10), 1),
			RootFinding::FindRootMuller(nonfinite, Complex(0), Complex(0.5), Complex(2), Real(1e-10), 1)}) {
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			CHECK(std::isnan(result.function_value.real()));
			CHECK_FALSE(result.error_message.empty());
		}
	}
}
