#include <MML.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Single header integration retains whole-domain estimates at the work limit", "[single_header][integration-contract]") {
	using namespace MML;
	int calls = 0;
	const auto result = Integration::IntegrateAdaptive2DCore([&](Real x, Real) {
		++calls;
		return std::pow(x, 14);
	}, 0, 1, 0, 1, Real(1e-10), 5, 226);
	CHECK(result.status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK(result.function_evaluations == calls);
	CHECK(calls == 225);
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(1) / 15, Precision::EqualityTolerance));
	CHECK_FALSE(result.converged);
}

TEST_CASE("Single header integration separates fixed completion and global convergence", "[single_header][integration-contract]") {
	using namespace MML;
	const auto fixed = Integration::IntegrateGK21Core([](Real x) { return x; }, 0, 1);
	CHECK(fixed.status == AlgorithmStatus::Completed);
	CHECK_FALSE(fixed.converged);
	CHECK(fixed.error_estimated);
	int calls = 0;
	const auto adaptive = Integration::IntegrateGKAdaptiveCore([&](Real x) {
		++calls;
		return std::exp(10 * x) - std::sinh(Real(10)) / 10;
	}, -1, 1, Real(1e-8), Real(1e-8), 10);
	CHECK(adaptive.converged);
	CHECK(adaptive.error_estimate <= Real(1e-8));
	CHECK(adaptive.function_evaluations == calls);
}

TEST_CASE("Single header closed-loop work accounts for pointwise cancellation", "[single_header][integration-contract]") {
	using namespace MML;
	const ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{std::cos(t), std::sin(t), 0}; });
	const VectorFunction<3> gradient([](const VectorN<Real, 3>& p) { return VectorN<Real, 3>{2 * p[0], 2 * p[1], 0}; });
	const auto result = PathIntegration::LineIntegralDetailed(gradient, circle, 0, 2 * Constants::PI, Real(1e-5));
	CHECK(result.converged);
	CHECK_THAT(result.value, Catch::Matchers::WithinAbs(Real(0), Real(1e-8)));
	CHECK(result.function_evaluations > 0);
}

TEST_CASE("Single header integration rejects nonfinite fixed-rule payloads", "[single_header][integration-contract]") {
	using namespace MML;
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	RealFunction invalid([](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
	const auto result = IntegrateGauss10Detailed(invalid, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.function_evaluations == 1);
	CHECK_FALSE(result.IsUsable());
}
