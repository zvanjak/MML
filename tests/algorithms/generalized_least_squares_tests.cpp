#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../TestPrecision.h"
#include <mml/algorithms/CurveFitting.h>
#include <mml/base/Function.h>
#include <limits>
#include <stdexcept>

namespace MML::Tests::GLS {
	using Catch::Matchers::WithinAbs;

	Vector<std::function<Real(Real)>> LinearBasis() {
		return Vector<std::function<Real(Real)>>{[](Real) { return Real(1); }, [](Real x) { return x; }};
	}

	Matrix<Real> CorrelatedCovariance() {
		return Matrix<Real>(4, 4, {2, 1, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});
	}

	TEST_CASE("GLS uses off-diagonal covariance and matches an exact correlated reference", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		const auto covariance = CorrelatedCovariance();
		const auto saved = covariance;
		const auto fit = GeneralizedLinearLeastSquares(x, y, covariance, LinearBasis());
		REQUIRE_THAT(fit.coefficients[0], WithinAbs(Real(7) / 8, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.coefficients[1], WithinAbs(Real(1), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.residual_norm * fit.residual_norm, WithinAbs(Real(21) / 8, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.mean_squared_error, WithinAbs(Real(21) / 32, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.r_squared, WithinAbs(Real(8) / 15, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.adjusted_r_squared, WithinAbs(Real(3) / 10, TOL(1e-12, 1e-5)));
		REQUIRE(fit.effective_rank == 2);
		REQUIRE(fit.num_data_points == 4);
		REQUIRE(fit.num_basis_functions == 2);
		REQUIRE(fit.weight_sum == Real(0));
		REQUIRE(covariance == saved);
		const auto diagonalOnly = WeightedGeneralLinearLeastSquares(x, y, Vector<Real>{Real(0.5), Real(0.5), 1, 1}, LinearBasis());
		REQUIRE_FALSE(IsNearlyEqual(fit.coefficients[0], diagonalOnly.coefficients[0], Precision::EqualityTolerance));
	}

	TEST_CASE("GLS identity and diagonal covariance reduce to existing least squares", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		const auto basis = LinearBasis();
		const auto identity = GeneralizedLinearLeastSquares(x, y, Matrix<Real>::Identity(4), basis);
		const auto ordinary = GeneralLinearLeastSquares(x, y, basis);
		Matrix<Real> covariance(4, 4);
		const Vector<Real> variances{1, 4, 2, 8};
		Vector<Real> weights(4);
		for (int i = 0; i < 4; ++i) {
			covariance(i, i) = variances[i];
			weights[i] = 1 / variances[i];
		}
		const auto diagonal = GeneralizedLinearLeastSquares(x, y, covariance, basis);
		const auto weighted = WeightedGeneralLinearLeastSquares(x, y, weights, basis);
		for (int i = 0; i < 2; ++i) {
			REQUIRE_THAT(identity.coefficients[i], WithinAbs(ordinary.coefficients[i], TOL(1e-12, 1e-5)));
			REQUIRE_THAT(diagonal.coefficients[i], WithinAbs(weighted.coefficients[i], TOL(1e-12, 1e-5)));
		}
		REQUIRE_THAT(identity.residual_norm, WithinAbs(ordinary.residual_norm, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(diagonal.residual_norm, WithinAbs(weighted.residual_norm, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(diagonal.r_squared, WithinAbs(weighted.r_squared, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("GLS is invariant to observation ordering and global covariance units", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		const auto covariance = CorrelatedCovariance();
		const auto reference = GeneralizedLinearLeastSquares(x, y, covariance, LinearBasis());
		for (Real scale : {std::sqrt(std::numeric_limits<Real>::min()), std::sqrt(std::numeric_limits<Real>::max())}) {
			const auto fit = GeneralizedLinearLeastSquares(x, y, covariance * scale, LinearBasis());
			for (int j = 0; j < 2; ++j)
				REQUIRE_THAT(fit.coefficients[j], WithinAbs(reference.coefficients[j], TOL(1e-12, 1e-5)));
			REQUIRE_THAT(fit.residual_norm * std::sqrt(scale), WithinAbs(reference.residual_norm, TOL(1e-12, 1e-5)));
			REQUIRE_THAT(fit.r_squared, WithinAbs(reference.r_squared, TOL(1e-12, 1e-5)));
			REQUIRE_THAT(fit.condition_number, WithinAbs(reference.condition_number, TOL(1e-12, 1e-5)));
		}
		const int order[]{3, 1, 0, 2};
		Vector<Real> shuffledX(4), shuffledY(4);
		Matrix<Real> shuffledCovariance(4, 4);
		for (int i = 0; i < 4; ++i) {
			shuffledX[i] = x[order[i]];
			shuffledY[i] = y[order[i]];
			for (int j = 0; j < 4; ++j) shuffledCovariance(i, j) = covariance(order[i], order[j]);
		}
		const auto shuffled = GeneralizedLinearLeastSquares(shuffledX, shuffledY, shuffledCovariance, LinearBasis());
		for (int j = 0; j < 2; ++j)
			REQUIRE_THAT(shuffled.coefficients[j], WithinAbs(reference.coefficients[j], TOL(1e-12, 1e-5)));
		REQUIRE_THAT(shuffled.residual_norm, WithinAbs(reference.residual_norm, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(shuffled.r_squared, WithinAbs(reference.r_squared, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("GLS satisfies covariance-weighted normal orthogonality", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		const auto covariance = CorrelatedCovariance();
		const auto fit = GeneralizedLinearLeastSquares(x, y, covariance, LinearBasis());
		Vector<Real> residual(4), precisionResidual;
		for (int i = 0; i < 4; ++i)
			residual[i] = y[i] - fit.coefficients[0] - fit.coefficients[1] * x[i];
		CholeskySolver<Real> cholesky(covariance);
		cholesky.Solve(residual, precisionResidual);
		Real constantMoment = 0, linearMoment = 0, quadratic = 0;
		for (int i = 0; i < 4; ++i) {
			constantMoment += precisionResidual[i];
			linearMoment += x[i] * precisionResidual[i];
			quadratic += residual[i] * precisionResidual[i];
		}
		REQUIRE_THAT(constantMoment, WithinAbs(Real(0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(linearMoment, WithinAbs(Real(0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(quadratic, WithinAbs(fit.residual_norm * fit.residual_norm, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("GLS reports minimum-norm rank-deficient fits honestly", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		Vector<std::function<Real(Real)>> basis{
			[](Real) { return Real(1); }, [](Real) { return Real(1); }};
		const auto fit = GeneralizedLinearLeastSquares(x, y, CorrelatedCovariance(), basis);
		REQUIRE(fit.effective_rank == 1);
		REQUIRE(std::isinf(fit.condition_number));
		REQUIRE_THAT(fit.coefficients[0], WithinAbs(Real(23) / 16, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.coefficients[1], WithinAbs(Real(23) / 16, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.r_squared, WithinAbs(Real(0), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.adjusted_r_squared, WithinAbs(Real(0), TOL(1e-12, 1e-5)));
		basis[0] = basis[1] = [](Real) { return Real(0); };
		const auto zero = GeneralizedLinearLeastSquares(x, y, CorrelatedCovariance(), basis);
		REQUIRE(zero.effective_rank == 0);
		REQUIRE(std::isinf(zero.condition_number));
		REQUIRE(zero.coefficients[0] == Real(0));
		REQUIRE(zero.coefficients[1] == Real(0));
		REQUIRE(std::isfinite(zero.residual_norm));
	}

	TEST_CASE("GLS handles exact square fits and a single observation", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1}, y{2, 5};
		const auto fit = GeneralizedLinearLeastSquares(x, y, Matrix<Real>(2, 2, {2, 1, 1, 2}), LinearBasis());
		REQUIRE_THAT(fit.coefficients[0], WithinAbs(Real(2), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.coefficients[1], WithinAbs(Real(3), TOL(1e-12, 1e-5)));
		REQUIRE_THAT(fit.residual_norm, WithinAbs(Real(0), TOL(1e-12, 1e-5)));
		REQUIRE(std::isnan(fit.adjusted_r_squared));
		Vector<std::function<Real(Real)>> constant{[](Real) { return Real(1); }};
		const auto one = GeneralizedLinearLeastSquares(Vector<Real>{0}, Vector<Real>{0}, Matrix<Real>(1, 1, {4}), constant);
		REQUIRE(one.coefficients[0] == Real(0));
		REQUIRE(one.residual_norm == Real(0));
		REQUIRE(std::isnan(one.r_squared));
		REQUIRE(std::isnan(one.adjusted_r_squared));
		const auto constantResponse = GeneralizedLinearLeastSquares(
			Vector<Real>{0, 1, 2, 3}, Vector<Real>{3, 3, 3, 3}, CorrelatedCovariance(), constant);
		REQUIRE_THAT(constantResponse.coefficients[0], WithinAbs(Real(3), TOL(1e-12, 1e-5)));
		REQUIRE(std::isnan(constantResponse.r_squared));
		REQUIRE(std::isnan(constantResponse.adjusted_r_squared));
	}

	TEST_CASE("GLS validates shapes and covariance before evaluating basis functions", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		int calls = 0;
		Vector<std::function<Real(Real)>> basis{[&calls](Real) { ++calls; return Real(1); }};
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(Vector<Real>{}, Vector<Real>{}, Matrix<Real>{}, basis), CurveFittingError);
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, Vector<Real>{1}, CorrelatedCovariance(), basis), CurveFittingError);
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, Matrix<Real>(3, 4), basis), CurveFittingError);
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, Matrix<Real>(4, 3), basis), CurveFittingError);
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(Vector<Real>{0}, Vector<Real>{0}, Matrix<Real>::Identity(1), LinearBasis()), CurveFittingError);
		Vector<std::function<Real(Real)>> empty;
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, CorrelatedCovariance(), empty), CurveFittingError);
		for (const auto& covariance : {
			Matrix<Real>(4, 4), Matrix<Real>(4, 4, Real(1)), Matrix<Real>::Identity(4) * Real(-1),
			Matrix<Real>(4, 4, {1, 2, 0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}),
			Matrix<Real>(4, 4, {2, 1, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1})})
			REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, covariance, basis), CurveFittingError);
		REQUIRE(calls == 0);
	}

	TEST_CASE("GLS rejects nonfinite observations covariance and basis values", "[CurveFitting][GLS]")
	{
		Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		auto covariance = CorrelatedCovariance();
		const auto originalX = x, originalY = y;
		for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			x[0] = invalid;
			REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, covariance, LinearBasis()), CurveFittingError);
			x = originalX;
			y[0] = invalid;
			REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, covariance, LinearBasis()), CurveFittingError);
			y = originalY;
			covariance(0, 1) = invalid;
			REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, covariance, LinearBasis()), CurveFittingError);
			covariance = CorrelatedCovariance();
			Vector<std::function<Real(Real)>> badBasis{[invalid](Real) { return invalid; }};
			REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, covariance, badBasis), CurveFittingError);
		}
		Vector<std::function<Real(Real)>> uncallable(1);
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, covariance, uncallable), CurveFittingError);
	}

	TEST_CASE("GLS detailed API honors statistics and exception policies", "[CurveFitting][GLS][Detailed]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		const auto fit = GeneralizedLinearFitDetailed(x, y, CorrelatedCovariance(), LinearBasis());
		REQUIRE(fit.IsSuccess());
		REQUIRE(fit.status == AlgorithmStatus::Completed);
		REQUIRE(fit.algorithm_name == "GeneralizedLinearLeastSquares");
		REQUIRE(fit.error_message.empty());
		REQUIRE(fit.function_evaluations == 8);
		REQUIRE(fit.elapsed_time_ms >= 0);
		CurveFittingConfig config;
		config.estimate_error = false;
		const auto noStatistics = GeneralizedLinearFitDetailed(x, y, CorrelatedCovariance(), LinearBasis(), config);
		REQUIRE(noStatistics.IsSuccess());
		REQUIRE(std::isfinite(noStatistics.residual_norm));
		REQUIRE(noStatistics.effective_rank == 2);
		REQUIRE(std::isnan(noStatistics.r_squared));
		REQUIRE(std::isnan(noStatistics.mean_squared_error));
		REQUIRE(std::isnan(noStatistics.adjusted_r_squared));
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto bad = GeneralizedLinearFitDetailed(x, y, Matrix<Real>(4, 4), LinearBasis(), config);
		REQUIRE(bad.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(bad.IsSuccess());
		REQUIRE_FALSE(bad.error_message.empty());
		REQUIRE(bad.coefficients.empty());
		REQUIRE_THROWS_AS(GeneralizedLinearFitDetailed(x, y, Matrix<Real>(4, 4), LinearBasis()), CurveFittingError);
	}

	TEST_CASE("GLS pointer basis overloads match callable overloads and reject nulls", "[CurveFitting][GLS]")
	{
		RealFunctionFromStdFunc constant([](Real) { return Real(1); });
		RealFunctionFromStdFunc linear([](Real x) { return x; });
		Vector<const IRealFunction*> basis{&constant, &linear};
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		const auto fit = GeneralizedLinearLeastSquares(x, y, CorrelatedCovariance(), basis);
		REQUIRE_THAT(fit.coefficients[0], WithinAbs(Real(7) / 8, TOL(1e-12, 1e-5)));
		const auto detailed = GeneralizedLinearFitDetailed(x, y, CorrelatedCovariance(), basis);
		REQUIRE(detailed.status == AlgorithmStatus::Completed);
		REQUIRE(detailed.function_evaluations == 8);
		basis[1] = nullptr;
		REQUIRE_THROWS_AS(GeneralizedLinearLeastSquares(x, y, CorrelatedCovariance(), basis), CurveFittingError);
		CurveFittingConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		REQUIRE(GeneralizedLinearFitDetailed(x, y, CorrelatedCovariance(), basis, config).status == AlgorithmStatus::InvalidInput);
	}

	TEST_CASE("GLS retains callback exceptions and reports arithmetic overflow explicitly", "[CurveFitting][GLS]")
	{
		const Vector<Real> x{0, 1, 2, 3}, y{0, 2, 4, 3};
		Vector<std::function<Real(Real)>> throwing(1);
		throwing[0] = [](Real) -> Real { throw std::runtime_error("GLS callback failed"); };
		REQUIRE_THROWS_WITH(GeneralizedLinearFitDetailed(x, y, CorrelatedCovariance(), throwing), "GLS callback failed");
		CurveFittingConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto failure = GeneralizedLinearFitDetailed(x, y, CorrelatedCovariance(), throwing, config);
		REQUIRE(failure.status == AlgorithmStatus::AlgorithmSpecificFailure);
		REQUIRE(failure.error_message == "GLS callback failed");
		const auto tinyCovariance = Matrix<Real>::Identity(4) * std::numeric_limits<Real>::denorm_min();
		const auto overflow = GeneralizedLinearFitDetailed(x, y, tinyCovariance, LinearBasis(), config);
		REQUIRE(overflow.status == AlgorithmStatus::NumericalInstability);
		REQUIRE_FALSE(overflow.error_message.empty());
	}
}
