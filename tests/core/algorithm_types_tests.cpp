#include <catch2/catch_test_macros.hpp>

#include "core/AlgorithmTypes.h"

using namespace MML;

namespace MML::Tests::Core::AlgorithmTypesTests {

TEST_CASE("AlgorithmTypes - EvaluationResultBase reports success by default", "[core][api]") {
	EvaluationResultBase result;

	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.error_message.empty());
	REQUIRE(result.algorithm_name.empty());
	REQUIRE(result.function_evaluations == 0);
	REQUIRE(result.IsSuccess());
	REQUIRE(static_cast<bool>(result));
}

TEST_CASE("AlgorithmTypes - DerivativeConfig has derivation-friendly defaults", "[core][api]") {
	DerivativeConfig config;

	REQUIRE(config.step == REAL(0.0));
	REQUIRE(config.estimate_error);
	REQUIRE(config.check_finite);
	REQUIRE(config.exception_policy == EvaluationExceptionPolicy::Propagate);
}

TEST_CASE("AlgorithmTypes - MakeEvaluationSuccessResult populates shared fields", "[core][api]") {
	auto result = MakeEvaluationSuccessResult<DerivativeResult<Real>>("NDer4Detailed", 5);

	REQUIRE(result.status == AlgorithmStatus::Success);
	REQUIRE(result.error_message.empty());
	REQUIRE(result.algorithm_name == "NDer4Detailed");
	REQUIRE(result.function_evaluations == 5);
	REQUIRE(result.IsSuccess());
	REQUIRE(static_cast<bool>(result));
	REQUIRE(result.step_used == REAL(0.0));
	REQUIRE(result.value == REAL(0.0));
	REQUIRE(result.error == REAL(0.0));
}

TEST_CASE("AlgorithmTypes - MakeEvaluationFailureResult populates failure diagnostics", "[core][api]") {
	auto result = MakeEvaluationFailureResult<DerivativeResult<Real>>(
		AlgorithmStatus::NumericalInstability,
		"non-finite derivative sample encountered",
		"NDer4Detailed",
		7);

	REQUIRE(result.status == AlgorithmStatus::NumericalInstability);
	REQUIRE(result.error_message == "non-finite derivative sample encountered");
	REQUIRE(result.algorithm_name == "NDer4Detailed");
	REQUIRE(result.function_evaluations == 7);
	REQUIRE_FALSE(result.IsSuccess());
	REQUIRE_FALSE(static_cast<bool>(result));
}

} // namespace MML::Tests::Core::AlgorithmTypesTests