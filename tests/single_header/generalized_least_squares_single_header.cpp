#include <MML.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Single header fits correlated observation errors with full GLS", "[single_header][GLS]")
{
	using namespace MML;
	const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
	const Matrix<Real> covariance(4, 4, {2, 1, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});
	const Vector<std::function<Real(Real)>> basis{
		[](Real) { return Real(1); }, [](Real value) { return value; }};
	const auto result = GeneralizedLinearFitDetailed(x, y, covariance, basis);
	REQUIRE(result.status == AlgorithmStatus::Completed);
	REQUIRE(result.IsSuccess());
	REQUIRE(result.effective_rank == 2);
	REQUIRE(IsNearlyEqual(result.coefficients[0], Real(7) / 8, Precision::EqualityTolerance));
	REQUIRE(IsNearlyEqual(result.coefficients[1], Real(1), Precision::EqualityTolerance));
	REQUIRE(IsNearlyEqual(result.residual_norm * result.residual_norm, Real(21) / 8, Precision::EqualityTolerance));
	CurveFittingConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	REQUIRE(GeneralizedLinearFitDetailed(x, y, Matrix<Real>(4, 4), basis, config).status == AlgorithmStatus::InvalidInput);
}
