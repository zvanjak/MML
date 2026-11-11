#include <catch2/catch_all.hpp>
#include <type_traits>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/ComplexFunction.h>
#include <mml/core/Derivation.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/core/ComplexAnalysis.h>
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::ComplexAnalysisTests
{
	// G7/K15 refuses budgets below 50*eps*integral(|integrand|).
	// The unit-circle derivative scales the integral by n!/(2*pi).
	Real CircleBudget(Real outputL1Scale = 1) {
		return std::max(Precision::ComplexAnalysisTolerance,
			Real(2) * Precision::IntegrationRoundoffFloorScale * outputL1Scale);
	}

	Real CountBudget() {
		// NDer4 balances h^4 truncation against eps/h roundoff. For the
		// exp(z^4) adversary, float needs refinement and a 5e-4 count budget;
		// this remains 500 times smaller than the unique-integer window.
		return Testing::Tol(Precision::DefaultToleranceStrict, Real(5e-4));
	}

	TEST_CASE("Contour quadrature resolves harmonics hidden by dyadic samples", "[complex][contour-contract]")
	{
		ComplexFunctionFromStdFunc one([](Complex) { return Complex(1, 0); });
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		const auto fourth = ComplexAnalysis::CauchyDerivative(one, unit, Complex(0, 0), 4, CircleBudget(24));
		CHECK(std::abs(fourth) < std::max(Precision::DefaultTolerance, Real(32 * 24) * Precision::Eps));

		ComplexFunctionFromStdFunc nonvanishing([](Complex z) { return std::exp(z * z * z * z); });
		CHECK(ComplexAnalysis::CountZeros(nonvanishing, unit, CountBudget()) == 0);
	}

	TEST_CASE("Complex contour rejects invalid tolerance before callbacks", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc one([&](Complex) {
			if (++calls > 128) throw ResourceLimitError("Regression probe exhausted");
			return Complex(1, 0);
		});
		ComplexAnalysis::LineSegmentContour line(Complex(0, 0), Complex(1, 0));
		ComplexAnalysis::ContourIntegrationConfig config;
		config.tolerance = 0;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		ComplexAnalysis::ContourIntegrationResult result;
		REQUIRE_NOTHROW(result = ComplexAnalysis::ContourIntegralDetailed(one, line, config));
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK(calls == 0);
	}

	TEST_CASE("Complex contour stops at a nonfinite sample", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc bad([&](Complex) {
			if (++calls > 128) throw ResourceLimitError("Regression probe exhausted");
			return Complex(std::numeric_limits<Real>::quiet_NaN(), 0);
		});
		ComplexAnalysis::LineSegmentContour line(Complex(0, 0), Complex(1, 0));
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		ComplexAnalysis::ContourIntegrationResult result;
		REQUIRE_NOTHROW(result = ComplexAnalysis::ContourIntegralDetailed(bad, line, config));
		CHECK_FALSE(result.IsUsable());
		CHECK(calls == 1);
	}

	TEST_CASE("Simple pole samples cancel analytic backgrounds", "[complex][contour-contract]")
	{
		ComplexFunctionFromStdFunc pole([](Complex z) { return Real(1) / z + Real(100000000); });
		const auto residue = ComplexAnalysis::ResidueSimplePole(pole, Complex(0, 0));
		CHECK(std::abs(residue - Complex(1, 0)) < Precision::DefaultTolerance);
	}

	TEST_CASE("Cauchy derivative rejects negative orders before callbacks", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc one([&](Complex) { ++calls; return Complex(1, 0); });
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		CHECK_THROWS(ComplexAnalysis::CauchyDerivative(one, unit, Complex(0, 0), -1));
		CHECK(calls == 0);
	}

	TEST_CASE("Contour iterations describe actual refinement", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc one([&](Complex) { ++calls; return Complex(1, 0); });
		ComplexAnalysis::LineSegmentContour line(Complex(0, 0), Complex(1, 0));
		const auto result = ComplexAnalysis::ContourIntegral(one, line);
		CHECK(result.converged);
		CHECK(result.function_evaluations == calls);
		CHECK(result.iterations == 0);
	}

	TEST_CASE("Nonuniform contour quadrature matches independent antiderivatives", "[complex][contour-contract]")
	{
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		for (int harmonic : {4, 8, 16, 32, 64}) {
			CAPTURE(harmonic);
			ComplexFunctionFromStdFunc monomial([=](Complex z) { return std::pow(z, Real(harmonic - 1)); });
			const Real budget = CircleBudget(2 * Constants::PI);
			const auto integral = ComplexAnalysis::ContourIntegralDetailed(monomial, unit, budget);
			REQUIRE(integral.IsUsable());
			// pow(z,k) and its oscillatory phase contribute O(k*eps) error.
			CHECK(std::abs(integral.value) < std::max(Precision::DefaultTolerance,
				Real(8 * harmonic) * Precision::Eps * 2 * Constants::PI));
			CHECK(integral.error_estimate <= budget);
		}
		ComplexFunctionFromStdFunc exponential([](Complex z) { return std::exp(z); });
		const Complex start(-1, 2), end(2, -1);
		const auto integral = ComplexAnalysis::ContourIntegralDetailed(exponential, ComplexAnalysis::LineSegmentContour(start, end),
			CircleBudget(std::abs(end - start) * std::exp(end.real())));
		REQUIRE(integral.IsUsable());
		CHECK(std::abs(integral.value - (std::exp(end) - std::exp(start))) < Precision::DefaultToleranceStrict);

		ComplexFunctionFromStdFunc one([](Complex) { return Complex(1, 0); });
		const auto reversed = ComplexAnalysis::ContourIntegralDetailed(one, ComplexAnalysis::ArcContour(Complex(0, 0), 1, Constants::PI, 0),
			CircleBudget(Constants::PI));
		CHECK(std::abs(reversed.ValueOrThrow("arc") - Complex(2, 0)) < Precision::DefaultToleranceStrict);
		const auto gamma = [](Real t) { return Complex(std::cos(t * t), std::sin(t * t)); };
		const auto tangent = [](Real t) { return Real(2) * t * Complex(-std::sin(t * t), std::cos(t * t)); };
		ComplexFunctionFromStdFunc cubic([](Complex z) { return z * z * z; });
		const auto reparameterized = ComplexAnalysis::ContourIntegralDetailed(cubic, gamma, tangent, Real(0), std::sqrt(2 * Constants::PI),
			CircleBudget(2 * Constants::PI));
		CHECK(std::abs(reparameterized.ValueOrThrow("reparameterization")) <
			std::max(Precision::DefaultToleranceStrict, Real(32) * Precision::Eps * 2 * Constants::PI));
	}

	TEST_CASE("Complex contour validates configuration and callbacks before work", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc one([&](Complex) { ++calls; return Complex(1, 0); });
		ComplexAnalysis::LineSegmentContour line(Complex(0, 0), Complex(1, 0));
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		for (Real tol : {Real(-1), Real(0), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
			config.tolerance = tol;
			const auto result = ComplexAnalysis::ContourIntegralDetailed(one, line, config);
			CHECK(result.status == AlgorithmStatus::InvalidInput);
			CHECK(result.function_evaluations == 0);
			CHECK_FALSE(result.error_estimated);
			CHECK_FALSE(result.error_message.empty());
		}
		config.tolerance = Precision::ComplexAnalysisTolerance;
		for (int field = 0; field < 3; ++field) {
			auto invalid = config;
			if (field == 0) invalid.max_depth = -1;
			if (field == 1) invalid.max_evaluations = 0;
			if (field == 2) invalid.initial_panels = 0;
			CHECK(ComplexAnalysis::ContourIntegralDetailed(one, line, invalid).status == AlgorithmStatus::InvalidInput);
		}
		const std::function<Complex(Real)> empty;
		const std::function<Complex(Real)> gamma = [](Real t) { return Complex(t, 0); };
		CHECK(ComplexAnalysis::ContourIntegralDetailed(one, empty, gamma, 0, 1, config).status == AlgorithmStatus::InvalidInput);
		CHECK(ComplexAnalysis::ContourIntegralDetailed(one, gamma, empty, 0, 1, config).status == AlgorithmStatus::InvalidInput);
		CHECK(ComplexAnalysis::ContourIntegralDetailed(one, gamma, gamma, 0, std::numeric_limits<Real>::infinity(), config).status ==
			AlgorithmStatus::InvalidInput);
		const auto zero = ComplexAnalysis::ContourIntegralDetailed(one, gamma, gamma, 1, 1, config);
		CHECK(zero.converged);
		CHECK(zero.value == Complex(0, 0));
		CHECK(zero.error_estimate == 0);
		CHECK(calls == 0);
		config.exception_policy = EvaluationExceptionPolicy::Propagate;
		config.tolerance = 0;
		CHECK_THROWS_AS(ComplexAnalysis::ContourIntegralDetailed(one, line, config), ArgumentError);
	}

	TEST_CASE("Complex contour enforces finite depth and function budgets", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc oscillatory([&](Complex z) { ++calls; return std::pow(z, Real(31)); });
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.max_depth = 0;
		const auto depth = ComplexAnalysis::ContourIntegralDetailed(oscillatory, unit, config);
		CHECK(depth.status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK_FALSE(depth.converged);
		CHECK(depth.error_estimated);
		CHECK(depth.iterations == 0);
		CHECK(depth.function_evaluations == calls);
		config.max_depth = Precision::ComplexAnalysisMaxDepth;
		config.max_evaluations = 32;
		calls = 0;
		const auto budget = ComplexAnalysis::ContourIntegralDetailed(oscillatory, unit, config);
		CHECK(budget.status == AlgorithmStatus::ResourceLimitExceeded);
		CHECK(budget.function_evaluations == calls);
		CHECK(calls == 32);
		CHECK_FALSE(budget.error_message.empty());
		config.max_evaluations = 1;
		calls = 0;
		const auto initial = ComplexAnalysis::ContourIntegralDetailed(oscillatory, unit, config);
		CHECK(initial.status == AlgorithmStatus::ResourceLimitExceeded);
		CHECK(calls == 0);

		ComplexFunctionFromStdFunc constant([](Complex) { return Complex(1, 0); });
		config.max_evaluations = Precision::ComplexAnalysisMaxEvaluations;
		config.max_depth = 2;
		config.tolerance = Precision::Eps * Precision::Eps;
		const auto impossible = ComplexAnalysis::ContourIntegralDetailed(constant,
			ComplexAnalysis::LineSegmentContour(Complex(0, 0), Complex(1, 0)), config);
		CHECK_FALSE(impossible.IsUsable());
		CHECK_FALSE(impossible.converged);
		CHECK(impossible.error_estimate > config.tolerance);
	}

	TEST_CASE("Complex theorem wrappers retain failures instead of returning payloads", "[complex][contour-contract]")
	{
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		ComplexFunctionFromStdFunc exponential([](Complex z) { return std::exp(z); });
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.max_evaluations = 1;
		const auto formula = ComplexAnalysis::CauchyIntegralFormulaDetailed(exponential, unit, Complex(0, 0), config);
		CHECK(formula.status == AlgorithmStatus::ResourceLimitExceeded);
		CHECK(formula.algorithm_name == "CauchyIntegralFormula");
		CHECK_THROWS_AS(ComplexAnalysis::CauchyIntegralFormula(exponential, unit, Complex(0, 0), config), ConvergenceError);
		CHECK_THROWS_AS(ComplexAnalysis::CauchyDerivative(exponential, unit, Complex(0, 0), 4, config), ConvergenceError);
		CHECK_THROWS_AS(ComplexAnalysis::Residue(exponential, Complex(0, 0), 1, config), ConvergenceError);
		CHECK_THROWS_AS(ComplexAnalysis::WindingNumber(unit, Complex(0, 0), config), ConvergenceError);
		CHECK_THROWS_AS(ComplexAnalysis::ArgumentPrinciple(exponential, unit, config), ConvergenceError);
		CHECK_THROWS_AS(ComplexAnalysis::CountZeros(exponential, unit, config), ConvergenceError);
		const auto count = ComplexAnalysis::CountZerosDetailed(exponential, unit, config);
		CHECK_FALSE(count.integer_certified);
		CHECK(count.status == AlgorithmStatus::ResourceLimitExceeded);
		CHECK(count.algorithm_name == "CountZeros");
	}

	TEST_CASE("Complex theorem inputs and transformed error units are validated", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc exponential([&](Complex z) { ++calls; return std::exp(z); });
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		CHECK(ComplexAnalysis::CauchyDerivativeDetailed(exponential, unit, Complex(0, 0), -1, config).status == AlgorithmStatus::InvalidInput);
		CHECK(ComplexAnalysis::CauchyDerivativeDetailed(exponential, unit, Complex(0, 0), std::numeric_limits<int>::max(), config).status ==
			AlgorithmStatus::InvalidInput);
		for (Complex point : {Complex(1, 0), Complex(2, 0), Complex(std::numeric_limits<Real>::quiet_NaN(), 0)})
			CHECK(ComplexAnalysis::CauchyIntegralFormulaDetailed(exponential, unit, point, config).status == AlgorithmStatus::InvalidInput);
		CHECK(ComplexAnalysis::WindingNumberDetailed(unit, Complex(1, 0), config).status == AlgorithmStatus::InvalidInput);
		CHECK(ComplexAnalysis::ResidueDetailed(exponential, Complex(0, 0), -1, config).status == AlgorithmStatus::InvalidInput);
		CHECK(calls == 0);
		config.tolerance = CircleBudget(24 * std::exp(Real(1)));
		const auto derivative = ComplexAnalysis::CauchyDerivativeDetailed(exponential, unit, Complex(0, 0), 4, config);
		REQUIRE(derivative.IsUsable());
		CHECK(std::abs(derivative.value - Complex(1, 0)) < config.tolerance);
		CHECK(derivative.error_estimate <= config.tolerance);
		CHECK(derivative.function_evaluations == calls);
		CHECK(derivative.elapsed_time_ms >= 0);

		config.estimate_error = false;
		config.tolerance = CountBudget();
		const auto hidden = ComplexAnalysis::CountZerosDetailed(exponential, unit, config);
		REQUIRE(hidden.integer_certified);
		CHECK(hidden.integer_value == 0);
		CHECK_FALSE(hidden.error_estimated);
		CHECK_FALSE(HasEstimate(hidden.error_estimate));
	}

	TEST_CASE("Complex count certification rejects unsafe conversions and uncertainty", "[complex][contour-contract]")
	{
		using Result = ComplexAnalysis::IntegerContourResult;
		for (Complex estimate : {Complex(2, 1), Complex(Real(2) + Real(1) / 3, 0),
			// Adding two rounds back to INT_MIN in float, making it a valid
			// boundary rather than an unsafe conversion. Use the next Real.
			Complex(std::nextafter(static_cast<Real>(std::numeric_limits<int>::max()) + Real(2), Constants::PosInf), 0),
			Complex(std::nextafter(static_cast<Real>(std::numeric_limits<int>::min()) - Real(2), Constants::NegInf), 0),
			Complex(0, std::numeric_limits<Real>::quiet_NaN()), Complex(std::numeric_limits<Real>::infinity(), 0)}) {
			Result result;
			static_cast<ComplexAnalysis::ContourIntegrationResult&>(result) = {estimate, 0, 0, true};
			ComplexAnalysis::Detail::CertifyInteger(result, Precision::DefaultTolerance);
			CHECK_FALSE(result.integer_certified);
			CHECK_FALSE(result.IsUsable());
			CHECK_THROWS_AS(result.IntegerOrThrow("regression"), ConvergenceError);
		}
		for (Real error : {Real(1), std::numeric_limits<Real>::quiet_NaN(), Real(-1)}) {
			Result result;
			static_cast<ComplexAnalysis::ContourIntegrationResult&>(result) = {Complex(2, 0), error, 0, true};
			ComplexAnalysis::Detail::CertifyInteger(result, Precision::DefaultTolerance);
			CHECK_FALSE(result.integer_certified);
		}
		// Float rounds INT_MAX upward, outside int's range. Exercise its
		// largest representable in-range integer instead of skipping coverage.
		Real upper = static_cast<Real>(std::numeric_limits<int>::max());
		if (static_cast<long double>(upper) > std::numeric_limits<int>::max())
			upper = std::nextafter(upper, Constants::NegInf);
		for (int integer : {std::numeric_limits<int>::min(), -3, 0, 3, static_cast<int>(upper)}) {
			Result result;
			static_cast<ComplexAnalysis::ContourIntegrationResult&>(result) = {Complex(static_cast<Real>(integer), 0), 0, 0, true};
			ComplexAnalysis::Detail::CertifyInteger(result, Precision::DefaultTolerance);
			REQUIRE(result.integer_certified);
			CHECK(result.IntegerOrThrow("boundary") == integer);
		}
	}

	TEST_CASE("Argument principle budgets include numerical differentiation work", "[complex][contour-contract]")
	{
		int calls = 0;
		ComplexFunctionFromStdFunc nonvanishing([&](Complex z) { ++calls; return std::exp(z * z * z * z); });
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		ComplexAnalysis::ContourIntegrationConfig config;
		config.tolerance = CountBudget();
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto count = ComplexAnalysis::CountZerosDetailed(nonvanishing, unit, config);
		REQUIRE(count.integer_certified);
		CHECK(count.integer_value == 0);
		CHECK(count.function_evaluations == calls);
		CHECK(count.derivative_error_estimate >= 0);
		CHECK(count.error_estimate <= config.tolerance);
		config.max_evaluations = 128;
		calls = 0;
		const auto limited = ComplexAnalysis::ArgumentPrincipleDetailed(nonvanishing, unit, config);
		CHECK(limited.status == AlgorithmStatus::ResourceLimitExceeded);
		CHECK(limited.function_evaluations == calls);
		CHECK(calls <= config.max_evaluations);
		ComplexFunctionFromStdFunc onBoundary([](Complex z) { return z - Real(1); });
		CHECK_FALSE(ComplexAnalysis::CountZerosDetailed(onBoundary, unit, config).IsUsable());
	}

	TEST_CASE("Simple pole result is an honest fixed-step estimate", "[complex][contour-contract]")
	{
		const Complex center(1, 1);
		const Real background = Real(1) / Precision::Eps;
		ComplexFunctionFromStdFunc pole([&](Complex z) { return Real(1) / (z - center) + background; });
		const auto result = ComplexAnalysis::ResidueSimplePoleDetailed(pole, center, Real(3) * Precision::Eps / 4);
		REQUIRE(result.IsUsable());
		CHECK(result.status == AlgorithmStatus::Completed);
		CHECK_FALSE(result.converged);
		CHECK_FALSE(result.error_estimated);
		CHECK_FALSE(HasEstimate(result.error_estimate));
		CHECK(result.function_evaluations == 4);
		CHECK(std::abs(result.value - Complex(1, 0)) < Precision::DefaultToleranceStrict);
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		for (Real h : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), Precision::Eps * Precision::Eps}) {
			const auto bad = ComplexAnalysis::ResidueSimplePoleDetailed(pole, center, h, config);
			CHECK(bad.status == AlgorithmStatus::InvalidInput);
			CHECK(bad.function_evaluations == 0);
		}
		ComplexFunctionFromStdFunc invalid([](Complex) { return Complex(0, std::numeric_limits<Real>::infinity()); });
		CHECK_FALSE(ComplexAnalysis::ResidueSimplePoleDetailed(invalid, Complex(0, 0), Precision::DefaultToleranceStrict, config).IsUsable());
	}

	TEST_CASE("Built-in contours reject invalid or collapsed geometry", "[complex][contour-contract]")
	{
		CHECK_THROWS_AS(ComplexAnalysis::CircleContour(Complex(0, 0), 0), ArgumentError);
		CHECK_THROWS_AS(ComplexAnalysis::CircleContour(Complex(0, 0), -1), ArgumentError);
		CHECK_THROWS_AS(ComplexAnalysis::CircleContour(Complex(0, 0), std::numeric_limits<Real>::infinity()), ArgumentError);
		CHECK_THROWS_AS(ComplexAnalysis::CircleContour(Complex(1, 1), Precision::Eps * Precision::Eps), ArgumentError);
		CHECK_THROWS_AS(ComplexAnalysis::ArcContour(Complex(0, 0), 1, 0, std::numeric_limits<Real>::quiet_NaN()), ArgumentError);
		CHECK_THROWS_AS(ComplexAnalysis::LineSegmentContour(Complex(0, 0), Complex(std::numeric_limits<Real>::infinity(), 0)), ArgumentError);
	}

	TEST_CASE("Interrupted complex theorem results retain scaled partial diagnostics", "[complex][contour-contract]")
	{
		ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
		ComplexAnalysis::ContourIntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.max_evaluations = 32;
		// Normalize the exponential amplitude so an absolute diagnostic budget resolves its rounding error.
		ComplexFunctionFromStdFunc integrand([](Complex z) { return std::exp(Real(20) * (z - Complex(1))) / z; });
		auto rawConfig = config;
		rawConfig.tolerance *= 2 * Constants::PI;
		const auto raw = ComplexAnalysis::ContourIntegralDetailed(integrand, unit, rawConfig);
		REQUIRE(raw.error_estimated);
		REQUIRE_FALSE(raw.IsUsable());
		int calls = 0;
		ComplexFunctionFromStdFunc interrupted([&](Complex z) {
			if (++calls == 33) throw std::runtime_error("Interrupted refinement");
			return std::exp(Real(20) * (z - Complex(1)));
		});
		config.max_evaluations = Precision::ComplexAnalysisMaxEvaluations;
		const auto scaled = ComplexAnalysis::CauchyIntegralFormulaDetailed(interrupted, unit, Complex(0, 0), config);
		CHECK(scaled.status == AlgorithmStatus::AlgorithmSpecificFailure);
		CHECK(scaled.function_evaluations == calls);
		CHECK(calls == 33);
		REQUIRE(scaled.error_estimated);
		CHECK(std::abs(scaled.value - raw.value / Complex(0, 2 * Constants::PI)) < Precision::DefaultToleranceStrict);
		CHECK_THAT(scaled.error_estimate, WithinRel(raw.error_estimate / (2 * Constants::PI), Precision::DefaultToleranceStrict));

		using LegacySignature = ComplexAnalysis::ContourIntegrationResult (*)(const IComplexFunction&,
			const ComplexAnalysis::CircleContour&, Real);
		const LegacySignature legacy = &ComplexAnalysis::ContourIntegral;
		const LegacySignature detailed = &ComplexAnalysis::ContourIntegralDetailed;
		CHECK(legacy != nullptr);
		CHECK(detailed != nullptr);
	}

	/*********************************************************************/
	/*****     IComplexFunction / ComplexFunction wrappers           *****/
	/*********************************************************************/

	TEST_CASE("ComplexFunction_pointer_wrapper", "[complex][interface]")
	{
		// f(z) = z^2
		ComplexFunction f([](Complex z) -> Complex { return z * z; });
		Complex z(3.0, 4.0);
		Complex result = f(z);
		// (3+4i)^2 = 9 + 24i - 16 = -7 + 24i
		REQUIRE_THAT(result.real(), WithinAbs(-7.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.imag(), WithinAbs(24.0, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("ComplexFunctionFromStdFunc_lambda_wrapper", "[complex][interface]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z; });
		Complex z(1.0, 1.0);
		Complex result = f(z);
		// (1+i)^3 = (1+i)(1+2i-1) = (1+i)(2i) = 2i + 2i^2 = -2 + 2i
		REQUIRE_THAT(result.real(), WithinAbs(-2.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.imag(), WithinAbs(2.0, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("RealToComplexFunction_wrapper", "[complex][interface]")
	{
		RealToComplexFunctionFromStdFunc gamma([](Real t) -> Complex {
			return Complex(std::cos(t), std::sin(t));
		});
		Complex result = gamma(Constants::PI / 2.0);
		REQUIRE_THAT(result.real(), WithinAbs(0.0, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(result.imag(), WithinAbs(1.0, TOL(1e-12, 1e-5)));
	}

	/*********************************************************************/
	/*****            Complex Derivatives                            *****/
	/*********************************************************************/

	TEST_CASE("NDer1Complex_z_squared", "[complex][derivation]")
	{
		// f(z) = z^2, f'(z) = 2z
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);  // 2*(2+i) = 4+2i

		Complex der = Derivation::NDer1Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-6, 5e-3)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-6, 5e-3)));
	}

	TEST_CASE("NDer2Complex_z_squared", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);

		// A dyadic step resolves this quadratic exactly, independently of contraction.
		Complex der = Derivation::NDer2Complex(f, z, Real(1) / 32);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-10, 1e-5)));

		DerivativeConfig config;
		config.estimate_error = true;
		const auto sampled = Derivation::NDer2ComplexDetailed(f, z, config);
		const Real h = sampled.step_used;
		// Complex products and rounded sample coordinates contribute O(epsilon/h).
		const Real roundoff = Precision::Eps *
			((std::abs(f(z + h)) + std::abs(f(z - h))) + std::abs(z) * std::abs(expected)) / h;
		REQUIRE(std::isfinite(sampled.error));
		REQUIRE(std::abs(sampled.value - expected) <= roundoff);
		REQUIRE(sampled.error <= roundoff);
		REQUIRE(std::abs(sampled.value - expected) <= sampled.error);
	}

	TEST_CASE("NDer4Complex_z_squared", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);

		Complex der = Derivation::NDer4Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("NDer6Complex_z_squared", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(2.0, 1.0);
		Complex expected(4.0, 2.0);

		Complex der = Derivation::NDer6Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("NDer4Complex_exp", "[complex][derivation]")
	{
		// f(z) = e^z, f'(z) = e^z
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		Complex z(1.0, Constants::PI / 4.0);
		Complex expected = std::exp(z);

		Complex der = Derivation::NDer4Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("NDer2Complex_sin", "[complex][derivation]")
	{
		// f(z) = sin(z), f'(z) = cos(z)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::sin(z); });
		Complex z(1.0, 0.5);
		Complex expected = std::cos(z);

		Complex der = Derivation::NDer2Complex(f, z);
		REQUIRE_THAT(der.real(), WithinAbs(expected.real(), TOL(1e-9, 1e-4)));
		REQUIRE_THAT(der.imag(), WithinAbs(expected.imag(), TOL(1e-9, 1e-4)));
	}

	TEST_CASE("NDer4ComplexDetailed_returns_error", "[complex][derivation]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		Complex z(1.0, 1.0);

		DerivativeConfig config;
		config.estimate_error = true;
		auto result = Derivation::NDer4ComplexDetailed(f, z, config);

		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.error >= 0.0);
		// Error should be very small for a polynomial
		REQUIRE(result.error < TOL(1e-8, 1e-4));
		// Value should be 2z = 2+2i
		REQUIRE_THAT(result.value.real(), WithinAbs(2.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value.imag(), WithinAbs(2.0, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****         Complex Root Finding - Newton                     *****/
	/*********************************************************************/

	TEST_CASE("NewtonComplex_z_squared_plus_1", "[complex][rootfinding]")
	{
		// f(z) = z^2 + 1, roots at z = ±i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + Complex(1.0, 0.0); });

		// Starting near +i
		auto result = RootFinding::FindRootNewtonComplex(f, Complex(MML::Real(0.1), MML::Real(1.2)));
		REQUIRE(result.converged);
		REQUIRE_THAT(result.root.real(), WithinAbs(0.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.root.imag(), WithinAbs(1.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("NewtonComplex_z_squared_plus_1_minus_i", "[complex][rootfinding]")
	{
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + Complex(1.0, 0.0); });

		// Starting near -i
		auto result = RootFinding::FindRootNewtonComplex(f, Complex(-MML::Real(0.1), -MML::Real(0.8)));
		REQUIRE(result.converged);
		REQUIRE_THAT(result.root.real(), WithinAbs(0.0, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.root.imag(), WithinAbs(-1.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("NewtonComplex_cubic", "[complex][rootfinding]")
	{
		// f(z) = z^3 - 1, roots at z = 1, e^(2πi/3), e^(4πi/3)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z - Complex(1.0, 0.0); });

		auto result = RootFinding::FindRootNewtonComplex(f, Complex(MML::Real(1.1), MML::Real(0.1)));
		REQUIRE(result.converged);
		// Should converge to real root z = 1
		Real abs_fz = std::abs(result.function_value);
		REQUIRE(abs_fz < TOL(1e-10, 1e-5));
	}

	/*********************************************************************/
	/*****         Complex Root Finding - Muller                     *****/
	/*********************************************************************/

	TEST_CASE("Muller_z_squared_plus_1", "[complex][rootfinding]")
	{
		// f(z) = z^2 + 1, roots at z = ±i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + Complex(1.0, 0.0); });

		// Muller can find complex roots from real starting points
		auto result = RootFinding::FindRootMuller(f, Complex(0.5, 0.5));
		REQUIRE(result.converged);
		Real abs_fz = std::abs(result.function_value);
		REQUIRE(abs_fz < TOL(1e-8, 1e-4));
		// Root should be +i or -i
		REQUIRE_THAT(std::abs(result.root), WithinAbs(1.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Muller_three_point_cubic", "[complex][rootfinding]")
	{
		// f(z) = z^3 - 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z - Complex(1.0, 0.0); });

		auto result = RootFinding::FindRootMuller(f, Complex(MML::Real(0.5), MML::Real(0.5)), Complex(MML::Real(0.8), MML::Real(0.2)), Complex(MML::Real(1.1), -MML::Real(0.1)));
		REQUIRE(result.converged);
		Real abs_fz = std::abs(result.function_value);
		REQUIRE(abs_fz < TOL(1e-8, 1e-4));
	}

	/*********************************************************************/
	/*****            Contour Classes                                *****/
	/*********************************************************************/

	TEST_CASE("CircleContour_unit_circle", "[complex][contour]")
	{
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);

		// At t=0: γ(0) = 1
		Complex z0 = unit_circle(0.0);
		REQUIRE_THAT(z0.real(), WithinAbs(1.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(z0.imag(), WithinAbs(0.0, TOL(1e-14, 1e-5)));

		// At t=π/2: γ(π/2) = i
		Complex z1 = unit_circle(Constants::PI / 2.0);
		REQUIRE_THAT(z1.real(), WithinAbs(0.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(z1.imag(), WithinAbs(1.0, TOL(1e-14, 1e-5)));

		// Derivative at t=0: γ'(0) = i
		Complex dz0 = unit_circle.derivative(0.0);
		REQUIRE_THAT(dz0.real(), WithinAbs(0.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(dz0.imag(), WithinAbs(1.0, TOL(1e-14, 1e-5)));
	}

	TEST_CASE("LineSegmentContour_basic", "[complex][contour]")
	{
		Complex z1(0.0, 0.0);
		Complex z2(2.0, 2.0);
		ComplexAnalysis::LineSegmentContour seg(z1, z2);

		// At t=0: start
		Complex s = seg(0.0);
		REQUIRE_THAT(s.real(), WithinAbs(0.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(s.imag(), WithinAbs(0.0, TOL(1e-14, 1e-5)));

		// At t=0.5: midpoint
		Complex m = seg(0.5);
		REQUIRE_THAT(m.real(), WithinAbs(1.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(m.imag(), WithinAbs(1.0, TOL(1e-14, 1e-5)));

		// At t=1: end
		Complex e = seg(1.0);
		REQUIRE_THAT(e.real(), WithinAbs(2.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(e.imag(), WithinAbs(2.0, TOL(1e-14, 1e-5)));

		// Derivative is constant: z2 - z1 = 2+2i
		Complex d = seg.derivative(0.5);
		REQUIRE_THAT(d.real(), WithinAbs(2.0, TOL(1e-14, 1e-5)));
		REQUIRE_THAT(d.imag(), WithinAbs(2.0, TOL(1e-14, 1e-5)));
	}

	/*********************************************************************/
	/*****          Contour Integration                              *****/
	/*********************************************************************/

	TEST_CASE("ContourIntegral_z_around_origin", "[complex][integration]")
	{
		// ∮ z dz around unit circle = 0 (z is analytic)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z; });
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);

		auto result = ComplexAnalysis::ContourIntegral(f, unit_circle, CircleBudget(2 * Constants::PI));
		REQUIRE_THAT(std::abs(result.value), WithinAbs(0.0, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("ContourIntegral_one_over_z_around_origin", "[complex][integration]")
	{
		// ∮ (1/z) dz around unit circle = 2πi
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return static_cast<MML::Real>(1.0) / z; });
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);

		auto result = ComplexAnalysis::ContourIntegral(f, unit_circle, CircleBudget(2 * Constants::PI));
		Complex expected(0.0, 2.0 * Constants::PI);

		REQUIRE_THAT(result.value.real(), WithinAbs(expected.real(), TOL(1e-8, 1e-4)));
		REQUIRE_THAT(result.value.imag(), WithinAbs(expected.imag(), TOL(1e-8, 1e-4)));

		// error_estimate is the accumulated Richardson estimate of the accepted leaves, not an
		// invented tol / |value|*1e-6 placeholder (analysis v4, CORE C10); it must bound the true error.
		// (Here f(gamma(t)) gamma'(t) = i is constant, Simpson is exact and the honest estimate is 0.)
		REQUIRE(result.converged);
		REQUIRE(result.error_estimate >= static_cast<MML::Real>(0.0));
		REQUIRE(std::abs(result.value - expected) <= std::max(result.error_estimate, TOL(1e-12, 1e-6)));
		static_assert(!std::is_convertible_v<ComplexAnalysis::ContourIntegrationResult, Complex>, "conversion must be explicit");
		REQUIRE(static_cast<Complex>(result) == result.value);
	}

	TEST_CASE("ContourIntegral_z_squared_analytic", "[complex][integration]")
	{
		// ∮ z^2 dz = 0 (analytic function, Cauchy theorem)
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		ComplexAnalysis::CircleContour circle(Complex(1.0, 0.0), 2.0);

		auto result = ComplexAnalysis::ContourIntegral(f, circle);
		REQUIRE_THAT(std::abs(result.value), WithinAbs(0.0, TOL(1e-7, 1e-3)));
	}

	TEST_CASE("ContourIntegral_line_segment", "[complex][integration]")
	{
		// ∫ z dz from 0 to 1+i along a line segment
		// Antiderivative: z^2/2, so integral = (1+i)^2/2 - 0 = (2i)/2 = i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z; });
		ComplexAnalysis::LineSegmentContour seg(Complex(0.0, 0.0), Complex(1.0, 1.0));

		auto result = ComplexAnalysis::ContourIntegral(f, seg);
		REQUIRE_THAT(result.value.real(), WithinAbs(0.0, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result.value.imag(), WithinAbs(1.0, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****             Winding Number                                *****/
	/*********************************************************************/

	TEST_CASE("WindingNumber_inside", "[complex][winding]")
	{
		// Winding number of unit circle around origin = 1
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);
		int n = ComplexAnalysis::WindingNumber(unit_circle, Complex(0.0, 0.0));
		REQUIRE(n == 1);
	}

	TEST_CASE("WindingNumber_outside", "[complex][winding]")
	{
		// Winding number of unit circle around point outside = 0
		ComplexAnalysis::CircleContour unit_circle(Complex(0.0, 0.0), 1.0);
		int n = ComplexAnalysis::WindingNumber(unit_circle, Complex(5.0, 0.0));
		REQUIRE(n == 0);
	}

	TEST_CASE("WindingNumber_off_center_inside", "[complex][winding]")
	{
		// Circle centered at (2,0) with radius 3, point at origin is inside
		ComplexAnalysis::CircleContour circle(Complex(2.0, 0.0), 3.0);
		int n = ComplexAnalysis::WindingNumber(circle, Complex(0.0, 0.0));
		REQUIRE(n == 1);
	}

	/*********************************************************************/
	/*****          Cauchy Integral Formula                           *****/
	/*********************************************************************/

	TEST_CASE("CauchyIntegralFormula_exp", "[complex][cauchy]")
	{
		// f(z) = e^z, evaluate at z0 = 0 using Cauchy integral
		// Should give f(0) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		Complex result = ComplexAnalysis::CauchyIntegralFormula(f, contour, Complex(0.0, 0.0));
		REQUIRE_THAT(result.real(), WithinAbs(1.0, 1e-6));
		REQUIRE_THAT(result.imag(), WithinAbs(0.0, 1e-6));
	}

	TEST_CASE("CauchyIntegralFormula_z_squared", "[complex][cauchy]")
	{
		// f(z) = z^2, evaluate at z0 = 1+i
		// Should give (1+i)^2 = 2i
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z; });
		ComplexAnalysis::CircleContour contour(Complex(1.0, 1.0), 2.0);

		Complex result = ComplexAnalysis::CauchyIntegralFormula(f, contour, Complex(1.0, 1.0), CircleBudget(8));
		REQUIRE_THAT(result.real(), WithinAbs(0.0, 1e-6));
		REQUIRE_THAT(result.imag(), WithinAbs(2.0, 1e-6));
	}

	TEST_CASE("CauchyDerivative_exp_first", "[complex][cauchy]")
	{
		// f(z) = e^z, f'(0) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		Complex result = ComplexAnalysis::CauchyDerivative(f, contour, Complex(0.0, 0.0), 1);
		REQUIRE_THAT(result.real(), WithinAbs(1.0, 1e-5));
		REQUIRE_THAT(result.imag(), WithinAbs(0.0, 1e-5));
	}

	TEST_CASE("CauchyDerivative_exp_second", "[complex][cauchy]")
	{
		// f(z) = e^z, f''(z) = e^z, f''(0) = 1
		// Uses exp to avoid Simpson aliasing with polynomial integrands
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return std::exp(z); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		Complex result = ComplexAnalysis::CauchyDerivative(f, contour, Complex(0.0, 0.0), 2, CircleBudget(2 * std::exp(Real(1))));
		REQUIRE_THAT(result.real(), WithinAbs(1.0, 1e-4));
		REQUIRE_THAT(result.imag(), WithinAbs(0.0, 1e-4));
	}

	/*********************************************************************/
	/*****            Residue Computation                            *****/
	/*********************************************************************/

	TEST_CASE("Residue_one_over_z", "[complex][residue]")
	{
		// f(z) = 1/z has Res(f, 0) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return static_cast<MML::Real>(1.0) / z; });

		Complex res = ComplexAnalysis::Residue(f, Complex(0.0, 0.0), 0.5, CircleBudget());
		REQUIRE_THAT(res.real(), WithinAbs(1.0, 1e-6));
		REQUIRE_THAT(res.imag(), WithinAbs(0.0, 1e-6));
	}

	TEST_CASE("Residue_one_over_z_squared_plus_1_at_i", "[complex][residue]")
	{
		// f(z) = 1/(z^2 + 1) = 1/((z-i)(z+i))
		// Res(f, i) = 1/(i + i) = 1/(2i) = -i/2
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return static_cast<MML::Real>(1.0) / (z * z + static_cast<MML::Real>(1.0)); });

		Complex res = ComplexAnalysis::Residue(f, Complex(0.0, 1.0), 0.5, CircleBudget());
		REQUIRE_THAT(res.real(), WithinAbs(0.0, 1e-5));
		REQUIRE_THAT(res.imag(), WithinAbs(-0.5, 1e-5));
	}

	TEST_CASE("ResidueSimplePole_one_over_z", "[complex][residue]")
	{
		// f(z) = 1/z, Res(f, 0) = lim_{z→0} z·(1/z) = 1
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return static_cast<MML::Real>(1.0) / z; });

		Complex res = ComplexAnalysis::ResidueSimplePole(f, Complex(0.0, 0.0));
		REQUIRE_THAT(res.real(), WithinAbs(1.0, 1e-5));
		REQUIRE_THAT(res.imag(), WithinAbs(0.0, 1e-5));
	}

	/*********************************************************************/
	/*****           Argument Principle / Count Zeros                *****/
	/*********************************************************************/

	TEST_CASE("ArgumentPrinciple_z_squared_plus_1", "[complex][argument]")
	{
		// f(z) = z^2 + 1, roots at ±i
		// A circle of radius 2 around origin contains both roots → N = 2
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + static_cast<MML::Real>(1.0); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 2.0);

		int n = ComplexAnalysis::ArgumentPrinciple(f, contour, CountBudget());
		REQUIRE(n == 2);
	}

	TEST_CASE("CountZeros_z_cubed", "[complex][argument]")
	{
		// f(z) = z^3 has a triple root at 0
		// Circle radius 1 around origin → N = 3
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z * z; });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 0.0), 1.0);

		int n = ComplexAnalysis::CountZeros(f, contour, CountBudget());
		REQUIRE(n == 3);
	}

	TEST_CASE("CountZeros_polynomial_partial", "[complex][argument]")
	{
		// f(z) = z^2 + 1, roots at ±i
		// Circle of radius 0.5 centered at i should contain exactly 1 root
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + static_cast<MML::Real>(1.0); });
		ComplexAnalysis::CircleContour contour(Complex(0.0, 1.0), 0.5);

		int n = ComplexAnalysis::CountZeros(f, contour, CountBudget());
		REQUIRE(n == 1);
	}

	TEST_CASE("CountZeros_no_zeros_inside", "[complex][argument]")
	{
		// f(z) = z^2 + 1, roots at ±i
		// Circle of radius 0.5 centered at (5, 0) contains no roots
		ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return z * z + static_cast<MML::Real>(1.0); });
		ComplexAnalysis::CircleContour contour(Complex(5.0, 0.0), 0.5);

		int n = ComplexAnalysis::CountZeros(f, contour, CountBudget());
		REQUIRE(n == 0);
	}

} // namespace MML::Tests::Core::ComplexAnalysisTests
