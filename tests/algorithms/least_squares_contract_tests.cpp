#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/CurveFitting.h>
#include <mml/base/Function.h>
#endif
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>
#include <stdexcept>

using namespace MML;
using Catch::Matchers::WithinAbs;

namespace {
	Vector<std::function<Real(Real)>> ConstantBasis() {
		Vector<std::function<Real(Real)>> basis(1);
		basis[0] = [](Real) { return Real(1); };
		return basis;
	}
}

TEST_CASE("Least squares line fitting preserves finite small-sample conventions and independent unit scales",
	"[CurveFitting][least-squares-contract]") {
	const auto one = LinearLeastSquaresDetailed(Vector<Real>{0}, Vector<Real>{7});
	CHECK(one.a == 0);
	CHECK(one.b == 7);
	CHECK(one.residual_norm == 0);
	CHECK(one.r_squared == 1);
	const Real large = std::numeric_limits<Real>::max() / 4;
	const auto horizontal = LinearLeastSquaresDetailed(Vector<Real>{2, 2}, Vector<Real>{large, large});
	CHECK(horizontal.a == 0);
	CHECK(horizontal.b == large);
	CHECK(horizontal.residual_norm == 0);
	const auto two = LinearLeastSquaresDetailed(Vector<Real>{0, 0}, Vector<Real>{1, 3});
	CHECK(two.a == 0);
	CHECK(two.b == 2);
	CHECK_THAT(two.residual_norm, WithinAbs(std::sqrt(Real(2)), Precision::LooseTolerance));
	CHECK_THAT(two.r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
	const Real tiny = std::scalbn(Real(1), std::numeric_limits<Real>::min_exponent / 2 - 8);
	const Real huge = std::scalbn(Real(1), std::numeric_limits<Real>::max_exponent / 2 + 8);
	for (Real xScale : {tiny, Real(1), huge}) {
		for (Real yScale : {tiny, Real(1), huge}) {
			CAPTURE(xScale, yScale);
			const Vector<Real> x{-2 * xScale, -xScale, 0, xScale, 2 * xScale};
			const Vector<Real> y{-yScale, yScale, 3 * yScale, 5 * yScale, 7 * yScale};
			const Real expectedSlope = 2 * (yScale / xScale);
			CurveFittingConfig config;
			config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
			const auto fit = LinearFitDetailed(x, y, config);
			if (!std::isfinite(expectedSlope)) {
				CHECK(fit.status == AlgorithmStatus::NumericalInstability);
			} else {
				REQUIRE(fit.IsSuccess());
				CHECK_THAT(fit.coefficients[0] / expectedSlope, WithinAbs(Real(1), Precision::LooseTolerance));
				CHECK_THAT(fit.coefficients[1] / yScale, WithinAbs(Real(3), Precision::LooseTolerance));
				CHECK_THAT(fit.r_squared, WithinAbs(Real(1), Precision::LooseTolerance));
				CHECK(std::isfinite(fit.residual_norm));
			}
		}
	}
}

TEST_CASE("Least squares validates every observation and callable before evaluating the basis",
	"[CurveFitting][least-squares-contract]") {
	int calls = 0;
	Vector<std::function<Real(Real)>> basis{[&](Real) { ++calls; return Real(1); }};
	const Vector<Real> x{0, 1, 2}, y{1, 2, 3}, weights{1, 1, 1}, penalty{0};
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		Vector<Real> invalid = x;
		invalid[2] = bad;
		for (bool badX : {false, true}) {
			const auto& xs = badX ? invalid : x;
			const auto& ys = badX ? y : invalid;
			CHECK_THROWS_AS(GeneralLinearLeastSquares(xs, ys, basis), CurveFittingError);
			CHECK_THROWS_AS(WeightedGeneralLinearLeastSquares(xs, ys, weights, basis), CurveFittingError);
			CHECK_THROWS_AS(TikhonovGeneralLinearLeastSquares(xs, ys, basis, penalty), CurveFittingError);
			CHECK_THROWS_AS(LinearLeastSquaresDetailed(xs, ys), CurveFittingError);
			CurveFittingConfig config;
			config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
			for (bool checkFinite : {false, true}) {
				config.check_finite = checkFinite;
				CHECK(GeneralLinearFitDetailed(xs, ys, basis, config).status == AlgorithmStatus::InvalidInput);
				CHECK(WeightedGeneralLinearFitDetailed(xs, ys, weights, basis, config).status == AlgorithmStatus::InvalidInput);
				CHECK(PolynomialFitDetailed(xs, ys, 0, config).status == AlgorithmStatus::InvalidInput);
				CHECK(LinearFitDetailed(xs, ys, config).status == AlgorithmStatus::InvalidInput);
			}
		}
	}
	CHECK(calls == 0);
	Vector<std::function<Real(Real)>> empty{basis[0], std::function<Real(Real)>{}};
	CHECK_THROWS_AS(GeneralLinearLeastSquares(x, y, empty), CurveFittingError);
	CHECK_THROWS_AS(WeightedGeneralLinearLeastSquares(x, y, weights, empty), CurveFittingError);
	CHECK_THROWS_AS(TikhonovGeneralLinearLeastSquares(x, y, empty, Vector<Real>{0, 0}), CurveFittingError);
	RealFunctionFromStdFunc constant([](Real) { return Real(1); });
	const Vector<const IRealFunction*> pointers{&constant, nullptr};
	CHECK_THROWS_AS(GeneralLinearLeastSquares(x, y, pointers), CurveFittingError);
	CHECK(calls == 0);
	for (Real bad : {Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(WeightedGeneralLinearLeastSquares(x, y, Vector<Real>{1, 1, bad}, basis), CurveFittingError);
		CHECK_THROWS_AS(TikhonovGeneralLinearLeastSquares(x, y, basis, Vector<Real>{bad}), CurveFittingError);
	}
	CHECK(calls == 0);
}

TEST_CASE("Least squares reports full-design and retained-subspace condition numbers separately",
	"[CurveFitting][least-squares-contract]") {
	const Vector<Real> x{0, 1, 2, 3}, y{1, 2, 3, 4}, weights{1, 1, 1, 1};
	Vector<std::function<Real(Real)>> duplicate{[](Real) { return Real(1); }, [](Real) { return Real(1); }};
	RealFunctionFromStdFunc constant([](Real) { return Real(1); });
	const Vector<const IRealFunction*> pointers{&constant, &constant};
	for (const auto& fit : {
		GeneralLinearLeastSquares(x, y, duplicate), GeneralLinearLeastSquares(x, y, pointers),
		WeightedGeneralLinearLeastSquares(x, y, weights, duplicate),
		TikhonovGeneralLinearLeastSquares(x, y, duplicate, Vector<Real>{0, 0}),
		GeneralizedLinearLeastSquares(x, y, Matrix<Real>::Identity(4), duplicate)}) {
		CHECK(fit.effective_rank == 1);
		CHECK(std::isinf(fit.condition_number));
		CHECK_THAT(fit.effective_condition_number, WithinAbs(Real(1), Precision::LooseTolerance));
		CHECK_THAT(fit.coefficients[0], WithinAbs(Real(5) / 4, Precision::LooseTolerance));
		CHECK_THAT(fit.coefficients[1], WithinAbs(Real(5) / 4, Precision::LooseTolerance));
		CHECK_THAT(fit.r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
		CHECK_THAT(fit.adjusted_r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
	}
	for (const auto& fit : {GeneralLinearFitDetailed(x, y, duplicate),
		WeightedGeneralLinearFitDetailed(x, y, weights, duplicate),
		GeneralizedLinearFitDetailed(x, y, Matrix<Real>::Identity(4), duplicate)}) {
		CHECK(fit.IsSuccess());
		CHECK(std::isinf(fit.condition_number));
		CHECK_THAT(fit.effective_condition_number, WithinAbs(Real(1), Precision::LooseTolerance));
	}
	duplicate[0] = duplicate[1] = [](Real) { return Real(0); };
	const auto zero = GeneralLinearLeastSquares(x, y, duplicate);
	CHECK(zero.effective_rank == 0);
	CHECK(std::isinf(zero.condition_number));
	CHECK(std::isinf(zero.effective_condition_number));
	CHECK(zero.coefficients[0] == 0);
	CHECK(zero.coefficients[1] == 0);
	duplicate[0] = [](Real) { return Real(1); };
	duplicate[1] = [](Real t) { return t; };
	const auto full = GeneralLinearFitDetailed(x, y, duplicate);
	CHECK(full.effective_rank == 2);
	CHECK(std::isfinite(full.condition_number));
	CHECK(full.condition_number == full.effective_condition_number);
}

TEST_CASE("Least squares optional statistics stay undefined without forcing overflowing squares",
	"[CurveFitting][least-squares-contract]") {
	const Real scale = std::sqrt(std::numeric_limits<Real>::max()) * 16;
	const Vector<Real> x{-1, 0, 1}, y{scale, -scale, scale}, weights{1, 1, 1};
	CurveFittingConfig config;
	config.estimate_error = false;
	for (const auto& fit : {LinearFitDetailed(x, y, config), GeneralLinearFitDetailed(x, y, ConstantBasis(), config),
		PolynomialFitDetailed(x, y, 0, config), WeightedGeneralLinearFitDetailed(x, y, weights, ConstantBasis(), config)}) {
		CHECK(fit.IsSuccess());
		CHECK(std::isfinite(fit.residual_norm));
		CHECK(std::isnan(fit.r_squared));
		CHECK(std::isnan(fit.mean_squared_error));
		CHECK(std::isnan(fit.adjusted_r_squared));
	}
	config.estimate_error = true;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	for (const auto& fit : {LinearFitDetailed(x, y, config), GeneralLinearFitDetailed(x, y, ConstantBasis(), config),
		PolynomialFitDetailed(x, y, 0, config), WeightedGeneralLinearFitDetailed(x, y, weights, ConstantBasis(), config)}) {
		CHECK(fit.status == AlgorithmStatus::NumericalInstability);
		CHECK_FALSE(fit.IsSuccess());
		CHECK_FALSE(fit.error_message.empty());
	}
	const auto square = PolynomialFit(Vector<Real>{0, 1}, Vector<Real>{1, 2}, 1);
	CHECK(std::isnan(square.adjusted_r_squared));
	const auto constantResponse = GeneralLinearLeastSquares(x, Vector<Real>{3, 3, 3}, ConstantBasis());
	CHECK(constantResponse.r_squared == 1);
}

TEST_CASE("Least squares validates nonfinite basis results and retains callback exception policies",
	"[CurveFitting][least-squares-contract]") {
	const Vector<Real> x{0, 1, 2}, y{1, 2, 3}, weights{1, 1, 1};
	CurveFittingConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		const Vector<std::function<Real(Real)>> basis{[bad](Real) { return bad; }};
		CHECK_THROWS_AS(GeneralLinearLeastSquares(x, y, basis), CurveFittingError);
		CHECK_THROWS_AS(WeightedGeneralLinearLeastSquares(x, y, weights, basis), CurveFittingError);
		CHECK_THROWS_AS(TikhonovGeneralLinearLeastSquares(x, y, basis, Vector<Real>{0}), CurveFittingError);
		CHECK(GeneralLinearFitDetailed(x, y, basis, config).status == AlgorithmStatus::InvalidInput);
	}
	const Vector<std::function<Real(Real)>> throwing{[](Real) -> Real { throw std::runtime_error("basis callback failed"); }};
	CHECK_THROWS_WITH(GeneralLinearFitDetailed(x, y, throwing), "basis callback failed");
	const auto failure = GeneralLinearFitDetailed(x, y, throwing, config);
	CHECK(failure.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(failure.error_message == "basis callback failed");
	config.estimate_error = false;
	const Vector<std::function<Real(Real)>> tiny{[](Real) { return std::numeric_limits<Real>::min(); }};
	const Real large = std::numeric_limits<Real>::max() / 4;
	const auto overflow = GeneralLinearFitDetailed(x, Vector<Real>{large, large, large}, tiny, config);
	CHECK(overflow.status == AlgorithmStatus::NumericalInstability);
	CHECK(overflow.coefficients.empty());
}

TEST_CASE("Least squares result evaluation validates borrowed callables and finite arithmetic",
	"[CurveFitting][least-squares-contract]") {
	const auto fit = GeneralLinearLeastSquares(Vector<Real>{0, 1, 2}, Vector<Real>{1, 1, 1}, ConstantBasis());
	const Vector<const IRealFunction*> null{nullptr};
	const Vector<std::function<Real(Real)>> empty{std::function<Real(Real)>{}};
	CHECK_THROWS_AS(fit.evaluate(Real(0), null), CurveFittingError);
	CHECK_THROWS_AS(fit.evaluate(Real(0), empty), CurveFittingError);
	CHECK_THROWS_AS(fit.evaluate(std::numeric_limits<Real>::quiet_NaN(), ConstantBasis()), CurveFittingError);
	const Vector<std::function<Real(Real)>> nonfinite{[](Real) { return std::numeric_limits<Real>::infinity(); }};
	CHECK_THROWS_AS(fit.evaluate(Real(0), nonfinite), CurveFittingError);
	GeneralLinearFitResult<Real> overflow;
	overflow.coefficients = Vector<Real>{std::numeric_limits<Real>::max()};
	const Vector<std::function<Real(Real)>> two{[](Real) { return Real(2); }};
	CHECK_THROWS_AS(overflow.evaluate(Real(0), two), MatrixNumericalError);
}

TEST_CASE("Least squares centering retains variation under a large common offset",
	"[CurveFitting][least-squares-contract]") {
	const Real offset = 1 / std::numeric_limits<Real>::epsilon();
	const Vector<Real> x{0, 1, 2, 3};
	for (Real units : {Real(1), Real(1) / 8}) {
		const Vector<Real> y{offset * units, (offset + 2) * units, (offset + 4) * units, (offset + 6) * units};
		const auto fit = GeneralLinearLeastSquares(x, y, ConstantBasis());
		const Real ratio = fit.residual_norm / units;
		CHECK_THAT(fit.r_squared, WithinAbs(1 - ratio * ratio / 20, Precision::LooseTolerance));
		const auto line = LinearLeastSquaresDetailed(y, Vector<Real>{0, 2 * units, 4 * units, 6 * units});
		CHECK_THAT(line.a, WithinAbs(Real(1), Precision::LooseTolerance));
		CHECK_THAT(line.b / (offset * units), WithinAbs(Real(-1), Precision::LooseTolerance));
		const Real roundingTolerance = std::max(Precision::LooseTolerance,
			Real(16) * std::numeric_limits<Real>::epsilon());
		CHECK(line.residual_norm / units < roundingTolerance);
	}
}

TEST_CASE("Least squares weighted statistics and coefficients retain independent response and weight units",
	"[CurveFitting][least-squares-contract]") {
	const Vector<Real> x{0, 1, 2};
	for (Real units : {std::numeric_limits<Real>::min() * 16, Real(1), std::sqrt(std::numeric_limits<Real>::max()) / 16}) {
		for (Real weightScale : {std::numeric_limits<Real>::denorm_min(), Real(1), std::numeric_limits<Real>::max() / 16}) {
			CAPTURE(units, weightScale);
			const auto fit = WeightedGeneralLinearLeastSquares(x, Vector<Real>{0, units, 2 * units},
				Vector<Real>{weightScale, 2 * weightScale, weightScale}, ConstantBasis());
			CHECK_THAT(fit.coefficients[0] / units, WithinAbs(Real(1), Precision::LooseTolerance));
			CHECK_THAT(fit.r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
			CHECK_THAT(fit.adjusted_r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
			CHECK(std::isfinite(fit.residual_norm));
			CHECK(std::isfinite(fit.mean_squared_error));
			CHECK(fit.weight_sum == 4 * weightScale);
		}
	}
}

TEST_CASE("Least squares rejects the zero-scale singular line instead of successful NaN coefficients",
	"[CurveFitting][least-squares-contract][least-squares-witness]") {
	const Vector<Real> x{0, 0, 0}, y{1, 2, 3};
	Real slope = 0, intercept = 0;
	CHECK_THROWS_AS(LinearLeastSquares(x, y, slope, intercept), SingularMatrixError);
	CurveFittingConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto fit = LinearFitDetailed(x, y, config);
	CHECK(fit.status == AlgorithmStatus::SingularMatrix);
	CHECK_FALSE(fit.IsSuccess());
	CHECK(fit.coefficients.empty());
	CHECK_FALSE(fit.error_message.empty());
}

TEST_CASE("Least squares rejects nonfinite responses instead of successful NaN coefficients",
	"[CurveFitting][least-squares-contract][least-squares-witness]") {
	const Vector<Real> x{0, 1, 2}, y{0, std::numeric_limits<Real>::quiet_NaN(), 2};
	CurveFittingConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto fit = GeneralLinearFitDetailed(x, y, ConstantBasis(), config);
	CHECK(fit.status == AlgorithmStatus::InvalidInput);
	CHECK_FALSE(fit.IsSuccess());
	CHECK(fit.coefficients.empty());
	CHECK_FALSE(fit.error_message.empty());
}

TEST_CASE("Least squares R squared is invariant to response units for constant-only fits",
	"[CurveFitting][least-squares-contract][least-squares-witness]") {
	const Vector<Real> x{0, 1, 2};
	for (Real scale : {Real(1), Real(1e-10), std::sqrt(std::numeric_limits<Real>::min()) * 16,
		std::sqrt(std::numeric_limits<Real>::max()) / 16}) {
		CAPTURE(scale);
		const Vector<Real> y{0, scale, 2 * scale};
		const auto ordinary = GeneralLinearLeastSquares(x, y, ConstantBasis());
		const auto weighted = WeightedGeneralLinearLeastSquares(x, y, Vector<Real>{1, 2, 1}, ConstantBasis());
		const auto regularized = RidgeGeneralLinearLeastSquares(x, y, ConstantBasis(), Real(0));
		for (const auto& fit : {ordinary, weighted, regularized}) {
			CHECK_THAT(fit.r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
			CHECK_THAT(fit.adjusted_r_squared, WithinAbs(Real(0), Precision::LooseTolerance));
			CHECK_THAT(fit.coefficients[0] / scale, WithinAbs(Real(1), Precision::LooseTolerance));
			CHECK(std::isfinite(fit.residual_norm));
		}
	}
}
