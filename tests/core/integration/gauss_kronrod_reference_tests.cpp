#include <catch2/catch_test_macros.hpp>
#include <mml/core/Integration.h>
#include <mml/base/Function.h>

#include <functional>

TEST_CASE("Adaptive Gauss-Kronrod accepts core function references without ADL ambiguity",
          "[gauss-kronrod][integration][reference-wrapper]")
{
	using namespace MML;
	RealFunction polynomial(+[](Real x) { return x * x; });
	const IRealFunction& function = polynomial;
	for (const auto rule : {Integration::GKRule::GK15, Integration::GKRule::GK21, Integration::GKRule::GK31}) {
		const auto result = Integration::IntegrateGKAdaptiveDetailed(
			std::cref(function), -Real(1), Real(1), Precision::IntegrationTolerance,
			Precision::IntegrationTolerance, 8, rule);
		REQUIRE(result.converged);
		REQUIRE(IsNearlyEqual(result.value, Real(2) / 3));
	}
}
