#include <catch2/catch_all.hpp>
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#include <mml/algorithms/MatrixSingularity.h>
#include <mml/core/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/FrameTags.h>

using namespace MML;

TEMPLATE_TEST_CASE("Algorithm-only inversion uses LU for dynamic fixed symmetric and tridiagonal matrices",
                   "[MatrixInversion][MatrixAlg]", Real, Complex)
{
	const Matrix<TestType> matrix(2, 2, {TestType(0), TestType(2), TestType(1), TestType(3)});
	const auto inverse = MatrixAlg::Inverse(matrix);
	REQUIRE((matrix * inverse).IsIdentity(Precision::LooseTolerance));
	MatrixNM<TestType, 2, 2> fixed{TestType(0), TestType(2), TestType(1), TestType(3)};
	const auto fixedInverse = MatrixAlg::Inverse(fixed);
	REQUIRE((fixed * fixedInverse).IsEqualTo(MatrixNM<TestType, 2, 2>::Identity(), Precision::LooseTolerance));
	MatrixSym<TestType> symmetric(2, {TestType(0), TestType(1), TestType(2)});
	REQUIRE((symmetric.ToFullMatrix() * MatrixAlg::Inverse(symmetric)).IsIdentity(Precision::LooseTolerance));
	// A nonsingular tridiagonal system with a zero leading pivot cannot use Thomas elimination.
	TridiagonalMatrix<TestType> tridiagonal(2, {TestType(0), TestType(2), TestType(1), TestType(3)});
	REQUIRE((tridiagonal.ToFullMatrix() * MatrixAlg::Inverse(tridiagonal)).IsIdentity(Precision::LooseTolerance));
	auto mutableMatrix = matrix;
	MatrixAlg::Invert(mutableMatrix);
	REQUIRE(mutableMatrix.IsEqualTo(inverse, Precision::LooseTolerance));
	MatrixAlg::Invert(fixed);
	REQUIRE(fixed.IsEqualTo(fixedInverse, Precision::LooseTolerance));
}

TEST_CASE("Algorithm inversion validates shape finite input and leaves failed mutations unchanged", "[MatrixInversion]")
{
	REQUIRE_THROWS_AS(MatrixAlg::Inverse(Matrix<Real>{}), MatrixDimensionError);
	REQUIRE_THROWS_AS(MatrixAlg::Inverse(Matrix<Real>(2, 3)), MatrixDimensionError);
	Matrix<Real> singular(2, 2, {1, 2, 2, 4});
	const auto before = singular;
	REQUIRE_THROWS_AS(MatrixAlg::Invert(singular), SingularMatrixError);
	REQUIRE(singular == before);
	for (Real value : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		Matrix<Real> invalid(2, 2, {1, value, 0, 1});
		REQUIRE_THROWS_AS(MatrixAlg::Inverse(invalid), MatrixNumericalError);
		Matrix<Complex> invalidComplex(1, 1, {Complex(1, value)});
		REQUIRE_THROWS_AS(MatrixAlg::Inverse(invalidComplex), MatrixNumericalError);
	}
	for (Real scale : {std::sqrt(std::numeric_limits<Real>::min()), std::sqrt(std::numeric_limits<Real>::max())}) {
		Matrix<Real> scaled(2, 2, {2 * scale, 0, 0, scale});
		REQUIRE((scaled * MatrixAlg::Inverse(scaled)).IsIdentity(Precision::LooseTolerance));
	}
}

TEST_CASE("Relocated typed metric and inverse policy preserve their geometry contracts", "[MatrixInversion][metric]")
{
	MatrixNM<Real, 2, 2> matrix{2, 0, 0, 4};
	const Metric<2, Cartesian2> metric(matrix);
	REQUIRE(IsNearlyEqual(metric.contravariantComponents()(0, 0), Real(1) / 2));
	REQUIRE(IsNearlyEqual(metric.contravariantComponents()(1, 1), Real(1) / 4));
	const auto inverse = Singularity::InverseMatrix(matrix);
	REQUIRE(IsNearlyEqual(inverse(1, 1), Real(1) / 4));
	MatrixNM<Real, 2, 2> singular;
	REQUIRE_THROWS_AS(Singularity::InverseMatrix(singular), SingularMatrixError);
	const auto zeros = Singularity::InverseMatrix(singular, SingularityPolicy::ReturnZero);
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			REQUIRE(zeros(i, j) == 0);
}
