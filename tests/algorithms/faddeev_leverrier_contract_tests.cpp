#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/MatrixAlg/Algebra.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#endif

#include <cmath>
#include <limits>

namespace {
	using namespace MML;
	using namespace MML::MatrixAlg;

	template<MMLScalar Scalar>
	void CheckInverse(const Matrix<Scalar>& matrix, const FaddeevLeVerrierResult<Scalar>& result)
	{
		REQUIRE(result.inverse.has_value());
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				CHECK(std::isfinite(std::abs((*result.inverse)(row, col))));
		const auto identity = Matrix<Scalar>::Identity(matrix.rows());
		CHECK((matrix * *result.inverse).IsEqualTo(identity, Precision::MatrixIsEqualTolerance));
		CHECK((*result.inverse * matrix).IsEqualTo(identity, Precision::MatrixIsEqualTolerance));
	}
}

TEST_CASE("Faddeev LeVerrier preserves representable coefficients across units",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * 3 / 4;
	for (Real scale : {Real(2), large, -large, Real(1) / large}) {
		CAPTURE(scale);
		const Matrix<Real> matrix{2, 2, {scale, 0, 0, scale}};
		const auto result = FaddeevLeVerrier(matrix);
		REQUIRE(result.characteristicCoefficients.size() == 3);
		CHECK(result.characteristicCoefficients[0] == Real(1));
		CHECK(std::isfinite(result.characteristicCoefficients[1]));
		CHECK(std::isfinite(result.characteristicCoefficients[2]));
		CHECK(std::abs(result.characteristicCoefficients[1] / scale + Real(2)) <= Precision::EqualityTolerance);
		CHECK(std::abs(result.characteristicCoefficients[2] / (scale * scale) - Real(1)) <= Precision::EqualityTolerance);
		CHECK(std::abs(result.determinant / (scale * scale) - Real(1)) <= Precision::EqualityTolerance);
		CheckInverse(matrix, result);
		REQUIRE(result.inverse.has_value());
		CHECK(std::abs((*result.inverse)(0, 0) * scale - Real(1)) <= Precision::EqualityTolerance);
		CHECK(std::abs((*result.inverse)(1, 1) * scale - Real(1)) <= Precision::EqualityTolerance);
	}
}

TEST_CASE("Faddeev LeVerrier handles scaled complex triangular matrices",
		  "[MatrixAlg][FaddeevLeVerrier][Complex]")
{
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * 3 / 4;
	for (Real scale : {Real(1), large, Real(1) / large}) {
		CAPTURE(scale);
		const Complex diagonal{scale, scale / 2};
		const Matrix<Complex> matrix{2, 2, {diagonal, Complex{scale, -scale}, Complex{}, diagonal}};
		const auto result = FaddeevLeVerrier(matrix);
		CHECK(std::abs(result.characteristicCoefficients[1] / diagonal + Complex{2}) <= Precision::EqualityTolerance);
		CHECK(std::abs(result.determinant / (diagonal * diagonal) - Complex{1}) <= Precision::EqualityTolerance);
		CheckInverse(matrix, result);
	}
}

TEST_CASE("Faddeev LeVerrier inverse does not divide by cancelled coefficients",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const Matrix<Real> matrix{2, 2, {Real(1e20), 0, 0, 1}};
	CheckInverse(matrix, FaddeevLeVerrier(matrix));
	CheckInverse(matrix, MatrixAnalyzer<Real>(matrix).FaddeevLeVerrier());
}

TEST_CASE("Faddeev LeVerrier absent inverse means LU detected singularity",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * 3 / 4;
	for (Real scale : {Real(1), large, Real(1) / large}) {
		const Matrix<Real> matrix{2, 2, {scale, scale, 0, 0}};
		const auto result = FaddeevLeVerrier(matrix);
		CHECK_FALSE(result.inverse.has_value());
		CHECK(result.determinant == Real(0));
		CHECK(result.characteristicCoefficients[2] == Real(0));
	}
	CHECK_FALSE(FaddeevLeVerrier(Matrix<Complex>{2, 2, {
		Complex{1, 1}, Complex{1, 1}, Complex{}, Complex{}}}).inverse.has_value());
}

TEST_CASE("Faddeev LeVerrier rejects nonfinite and unrepresentable results",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{}), MatrixDimensionError);
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{2, 3}), MatrixDimensionError);
	for (Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{1, 1, {invalid}}), MatrixNumericalError);
		REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Complex>{1, 1, {Complex{1, invalid}}}), MatrixNumericalError);
	}
	const Real large = std::sqrt(std::numeric_limits<Real>::max()) * 2;
	const Real small = std::sqrt(std::numeric_limits<Real>::denorm_min()) / 2;
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{2, 2, {large, 0, 0, large}}),
					  MatrixNumericalError);
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{2, 2, {small, 0, 0, small}}),
					  MatrixNumericalError);
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{1, 1, {std::numeric_limits<Real>::denorm_min()}}),
					  MatrixNumericalError);
}

TEST_CASE("Faddeev LeVerrier higher order coefficients retain their units",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const Real large = std::pow(std::numeric_limits<Real>::max(), Real(1) / 8);
	for (Real scale : {Real(1), large, -large, Real(1) / large}) {
		CAPTURE(scale);
		const Matrix<Real> matrix{4, 4, {
			scale, scale, 0, 0,
			0, 2 * scale, scale, 0,
			0, 0, 3 * scale, scale,
			0, 0, 0, 4 * scale}};
		const auto result = FaddeevLeVerrier(matrix);
		const Real coefficients[]{1, -10, 35, -50, 24};
		Real scalePower = 1;
		for (int index = 0; index <= 4; ++index) {
			CHECK(std::isfinite(result.characteristicCoefficients[index]));
			CHECK(std::abs(result.characteristicCoefficients[index] / scalePower - coefficients[index])
				  <= std::abs(coefficients[index]) * Precision::EqualityTolerance);
			if (index != 4)
				scalePower *= scale;
		}
		CHECK(std::abs(result.determinant / scalePower - Real(24)) <= Real(24) * Precision::EqualityTolerance);
		CheckInverse(matrix, result);
	}
}

TEST_CASE("Faddeev LeVerrier rejects normalization that destroys input",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const Real large = std::numeric_limits<Real>::max() / 2;
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Real>{2, 2, {large, tiny, 0, 1}}), MatrixNumericalError);
	REQUIRE_THROWS_AS(FaddeevLeVerrier(Matrix<Complex>{1, 1, {Complex{large, tiny}}}), MatrixNumericalError);
}

TEST_CASE("Faddeev LeVerrier checks finite inverse residuals before publishing",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const int size = std::numeric_limits<Real>::digits <= 24 ? 5 : 10;
	Matrix<Real> hilbert(size, size);
	for (int row = 0; row < size; ++row)
		for (int col = 0; col < size; ++col)
			hilbert(row, col) = Real(1) / (row + col + 1);
	const auto inverse = Inverse(hilbert);
	const auto product = hilbert * inverse;
	Real maxResidual = 0;
	for (int row = 0; row < size; ++row)
		for (int col = 0; col < size; ++col) {
			REQUIRE(std::isfinite(inverse(row, col)));
			maxResidual = std::max(maxResidual, std::abs(product(row, col) - (row == col ? Real(1) : Real(0))));
		}
	REQUIRE(maxResidual > static_cast<Real>(size) * Precision::MatrixIsEqualTolerance);
	REQUIRE_THROWS_WITH(FaddeevLeVerrier(hilbert),
		Catch::Matchers::ContainsSubstring("LU inverse failed two-sided identity validation"));
}

TEST_CASE("Faddeev LeVerrier scalar and zero controls stay usable",
		  "[MatrixAlg][FaddeevLeVerrier]")
{
	const auto zero = FaddeevLeVerrier(Matrix<Real>{3, 3});
	CHECK(zero.characteristicCoefficients == Vector<Real>{1, 0, 0, 0});
	CHECK(zero.determinant == Real(0));
	CHECK_FALSE(zero.inverse.has_value());
	for (Real value : {Real(2), std::numeric_limits<Real>::max(), std::numeric_limits<Real>::min()}) {
		CAPTURE(value);
		const Matrix<Real> matrix{1, 1, {value}};
		const auto result = FaddeevLeVerrier(matrix);
		CHECK(result.characteristicCoefficients[1] == -value);
		CHECK(result.determinant == value);
		CheckInverse(matrix, result);
	}
	const Matrix<Complex> scalar{1, 1, {Complex{2, -1}}};
	CheckInverse(scalar, FaddeevLeVerrier(scalar));
}
