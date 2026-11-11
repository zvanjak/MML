#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#if !defined(MML_SINGLE_HEADER_INCLUDE)
#include <mml/base/BaseUtils.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/MMLSingularityHandling.h>
#endif

#include <limits>

using namespace MML;

namespace {
	const Real largeScale = std::pow(Real(10), std::numeric_limits<Real>::max_exponent10 >= 200 ? Real(200) : Real(25));
	const Real smallScale = Real(1) / largeScale;
}

TEST_CASE("Base comparisons never certify nonfinite data", "[BaseUtilityContract][BaseUtilityFinite]")
{
	for (const Real bad : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CHECK_FALSE(Utils::AreEqual(Complex{bad, 0}, Complex{bad, 0}));
		CHECK_FALSE(Utils::AreEqualAbs(Complex{bad, 0}, Complex{}));
		CHECK_FALSE(Utils::AreEqual(Vector<Real>{bad}, Vector<Real>{0}));
		CHECK_FALSE(Utils::AreEqual(Vector<Complex>{Complex{0, bad}}, Vector<Complex>{Complex{}}));
		CHECK_FALSE(Utils::AreEqualAbs(Vector<Complex>{Complex{bad, 0}}, Vector<Complex>{Complex{}}));
		CHECK_FALSE(Vector<Real>{bad}.IsNearZero());
		CHECK_FALSE(VectorN<Real, 1>{bad}.IsEqualTo(VectorN<Real, 1>{0}));
		CHECK_FALSE(VectorN<Real, 1>{bad}.IsNearZero());
		CHECK_FALSE(VectorN<Complex, 1>{Complex{0, bad}}.IsEqualTo(VectorN<Complex, 1>{Complex{}}));
		CHECK_FALSE(Matrix<Real>{1, 1, {bad}}.IsEqualTo(Matrix<Real>{1, 1, {0}}));
		CHECK_FALSE(Matrix<Complex>{1, 1, {Complex{bad, 0}}}.IsEqualTo(Matrix<Complex>{1, 1, {Complex{}}}));
		MatrixNM<Real, 1, 1> fixed;
		fixed(0, 0) = bad;
		CHECK_FALSE(fixed.IsEqualTo(MatrixNM<Real, 1, 1>{}));
		MatrixSym<Real> symmetric(1);
		symmetric(0, 0) = bad;
		CHECK_FALSE(symmetric.IsEqualTo(MatrixSym<Real>(1)));
		CHECK_FALSE(IsOrthonormal(std::vector<VectorN<Real, 1>>{VectorN<Real, 1>{bad}}));
	}
	CHECK(Utils::AreEqual(Complex{1, 2}, Complex{1, 2}, 0));
	CHECK(Vector<Real>{0}.IsNearZero(0));
}

TEST_CASE("Base predicates reject invalid comparison tolerances", "[BaseUtilityContract]")
{
	for (Real tolerance : {-Real(1), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity()}) {
		CHECK_FALSE(Utils::AreEqual(Complex{}, Complex{}, tolerance));
		CHECK_FALSE(Utils::AreEqualAbs(Complex{}, Complex{}, tolerance));
		CHECK_FALSE(Utils::AreEqual(Vector<Real>{}, Vector<Real>{}, tolerance));
		CHECK_FALSE(Utils::AreEqual(Vector<Complex>{}, Vector<Complex>{}, tolerance));
		CHECK_FALSE(Utils::AreEqualAbs(Vector<Complex>{}, Vector<Complex>{}, tolerance));
		CHECK_FALSE(Vector<Real>{}.IsNearZero(tolerance));
		CHECK_FALSE(VectorN<Real, 1>{0}.IsEqualTo(VectorN<Real, 1>{0}, tolerance));
		CHECK_FALSE(VectorN<Real, 1>{0}.IsNearZero(tolerance));
		CHECK_FALSE(Matrix<Real>{}.IsEqualTo(Matrix<Real>{}, tolerance));
		CHECK_FALSE(MatrixNM<Real, 1, 1>{}.IsEqualTo(MatrixNM<Real, 1, 1>{}, tolerance));
		CHECK_FALSE(MatrixSym<Real>(1).IsEqualTo(MatrixSym<Real>(1), tolerance));
		CHECK_FALSE(IsOrthonormal(std::vector<VectorN<Real, 1>>{}, tolerance));
	}
}

TEST_CASE("Matrix utility factories validate shapes before copying", "[BaseUtilityContract]")
{
	const std::vector<Vector<Real>> empty;
	CHECK(Utils::MatrixFromVectorsInRows(empty).rows() == 0);
	CHECK(Utils::MatrixFromVectorsInColumns(empty).cols() == 0);
	for (const auto& ragged : {std::vector<Vector<Real>>{Vector<Real>{1}, Vector<Real>{2, 3}},
		std::vector<Vector<Real>>{Vector<Real>{1, 2}, Vector<Real>{3}}}) {
		CHECK_THROWS_AS(Utils::MatrixFromVectorsInRows(ragged), MatrixDimensionError);
		CHECK_THROWS_AS(Utils::MatrixFromVectorsInColumns(ragged), MatrixDimensionError);
	}
	const std::vector<Vector<Real>> rectangular{Vector<Real>{1, 2}, Vector<Real>{3, 4}};
	CHECK(Utils::MatrixFromVectorsInRows(rectangular).IsEqualTo(Matrix<Real>{2, 2, {1, 2, 3, 4}}));
	CHECK(Utils::MatrixFromVectorsInColumns(rectangular).IsEqualTo(Matrix<Real>{2, 2, {1, 3, 2, 4}}));
	CHECK_THROWS_AS(Utils::SimilarityTransform(Matrix<Real>(2, 3), Matrix<Real>(2, 2)), MatrixDimensionError);
	CHECK_THROWS_AS(Utils::SimilarityTransform(Matrix<Real>(3, 3), Matrix<Real>(2, 2)), MatrixDimensionError);
	CHECK_THROWS_AS(Utils::SimilarityTransform(Matrix<Real>(2, 2), Matrix<Real>(2, 3)), MatrixDimensionError);
	CHECK(Utils::SimilarityTransform(Matrix<Real>{}, Matrix<Real>{}).rows() == 0);
}

TEST_CASE("Maximum matrix difference rejects nonfinite entries", "[BaseUtilityContract]")
{
	const Matrix<Real> zero{1, 1, {0}};
	CHECK_THROWS_AS(Utils::MaxAbsDiff(Matrix<Real>{1, 1, {std::numeric_limits<Real>::quiet_NaN()}}, zero),
		NumericInputError);
	CHECK_THROWS_AS(Utils::MaxAbsDiff(Matrix<Real>{1, 1, {std::numeric_limits<Real>::infinity()}}, zero),
		NumericInputError);
	CHECK_THROWS_AS(Utils::MaxAbsDiff(
		Matrix<Complex>{1, 1, {Complex{0, std::numeric_limits<Real>::quiet_NaN()}}}, Matrix<Complex>{1, 1, {Complex{}}}),
		NumericInputError);
	CHECK_THROWS_AS(Utils::MaxAbsDiff(Matrix<Real>{1, 1, {std::numeric_limits<Real>::max()}},
		Matrix<Real>{1, 1, {-std::numeric_limits<Real>::max()}}), NumericalMethodError);
	CHECK(Utils::MaxAbsDiff(Matrix<Real>{2, 1, {1, 99}}, zero) == 1);
}

TEST_CASE("Angle wrapping has bounded reduction and half-open endpoints", "[BaseUtilityContract]")
{
	for (Real angle : {Real(1e20), -Real(1e20), std::numeric_limits<Real>::max(),
		-std::numeric_limits<Real>::max(), Real(0), -Real(0), Constants::PI, -Constants::PI,
		2 * Constants::PI, -2 * Constants::PI, -std::numeric_limits<Real>::denorm_min()}) {
		const Real positive = Utils::AngleTo2PiRange(angle);
		const Real symmetric = Utils::AngleToPiPiRange(angle);
		CHECK(std::isfinite(positive));
		CHECK(positive >= 0);
		CHECK(positive < 2 * Constants::PI);
		CHECK(symmetric >= -Constants::PI);
		CHECK(symmetric < Constants::PI);
		CHECK(symmetric == normalizeAngle(angle));
	}
	CHECK(Utils::AngleToPiPiRange(Constants::PI) == -Constants::PI);
	CHECK(Utils::AngleToPiPiRange(-Constants::PI) == -Constants::PI);
	CHECK(Utils::AngleTo2PiRange(2 * Constants::PI) == 0);
	CHECK(Utils::AngleToPiPiRange(Precision::Eps) == Precision::Eps);
	const Real maximum = std::numeric_limits<Real>::max();
	const Real reducedDifference = normalizeAngle(normalizeAngle(maximum) - normalizeAngle(-maximum));
	CHECK(AnglesAreEqual(maximum, -maximum, Precision::EqualityTolerance) ==
		(std::abs(reducedDifference) < Precision::EqualityTolerance));
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(Utils::AngleTo2PiRange(bad), DomainError);
		CHECK_THROWS_AS(Utils::AngleToPiPiRange(bad), DomainError);
		CHECK_FALSE(AnglesAreEqual(bad, 0, Precision::EqualityTolerance));
	}
}

TEST_CASE("Cartesian coordinate radii survive huge and tiny finite scales", "[BaseUtilityContract]")
{
	for (Real scale : {largeScale, smallScale}) {
		Real radius, phi, theta, z;
		Utils::CartesianToPolar(3 * scale, 4 * scale, radius, phi);
		CHECK_THAT(radius / scale, Catch::Matchers::WithinRel(Real(5), Precision::EqualityTolerance));
		Utils::CartesianToCylindrical(3 * scale, 4 * scale, 12 * scale, radius, phi, z);
		CHECK_THAT(radius / scale, Catch::Matchers::WithinRel(Real(5), Precision::EqualityTolerance));
		CHECK(z == 12 * scale);
		Utils::CartesianToSpherical(3 * scale, 4 * scale, 12 * scale, radius, theta, phi);
		CHECK_THAT(radius / scale, Catch::Matchers::WithinRel(Real(13), Precision::EqualityTolerance));
		CHECK_THAT(theta, Catch::Matchers::WithinRel(std::atan2(Real(5), Real(12)), Precision::EqualityTolerance));
	}
	Real radius, phi, theta;
	Utils::CartesianToSpherical(0, 0, 0, radius, theta, phi);
	CHECK(radius == 0);
	CHECK(theta == 0);
	CHECK(phi == 0);
	CHECK_THROWS_AS(Utils::CartesianToPolar(std::numeric_limits<Real>::infinity(), 1, radius, phi), NumericInputError);
	CHECK_THROWS_AS(Utils::CartesianToSpherical(1, 1, std::numeric_limits<Real>::quiet_NaN(), radius, theta, phi),
		NumericInputError);
	CHECK_THROWS_AS(Utils::CartesianToPolar(std::numeric_limits<Real>::max(),
		std::numeric_limits<Real>::max(), radius, phi), NumericalMethodError);
}

TEST_CASE("Vector projection preserves representable huge tiny and mixed-scale results", "[BaseUtilityContract]")
{
	for (Real scale : {largeScale, smallScale, std::numeric_limits<Real>::max() * Real(.75)}) {
		const Vector<Real> vector{scale, scale};
		const auto projection = Utils::VectorProjectionParallelTo(vector, vector);
		for (int i = 0; i < vector.size(); ++i)
			CHECK_THAT(projection[i] / scale, Catch::Matchers::WithinRel(Real(1), Precision::EqualityTolerance));
	}
	const Real large = largeScale, small = smallScale;
	const auto mixed = Utils::VectorProjectionParallelTo(Vector<Real>{large, small}, Vector<Real>{small, large});
	CHECK_THAT(mixed[1] / small, Catch::Matchers::WithinRel(Real(2), Precision::EqualityTolerance));
	const auto orthogonal = Utils::VectorProjectionParallelTo(Vector<Real>{large, 0}, Vector<Real>{0, small});
	CHECK(orthogonal.IsZero());
	const Real maximum = std::numeric_limits<Real>::max();
	const auto cancellation = Utils::VectorProjectionParallelTo(Vector<Real>{maximum, maximum, -maximum},
		Vector<Real>{1, 1, 1});
	CHECK_THAT(cancellation[0] / maximum, Catch::Matchers::WithinRel(Real(1) / 3, Precision::EqualityTolerance));
	CHECK_THROWS_AS(Utils::VectorProjectionParallelTo(Vector<Real>{1}, Vector<Real>{0}), ArgumentError);
	CHECK_THROWS_AS(Utils::VectorProjectionParallelTo(Vector<Real>{1}, Vector<Real>{1, 2}), VectorDimensionError);
	CHECK_THROWS_AS(Utils::VectorProjectionParallelTo(Vector<Real>{1},
		Vector<Real>{std::numeric_limits<Real>::infinity()}), NumericInputError);
}

TEST_CASE("SafeDivide validates inputs and preserves signed singular policies", "[BaseUtilityContract]")
{
	using namespace Singularity;
	const Real tolerance = DEFAULT_SINGULARITY_TOL;
	for (Real numerator : {Real(1), -Real(1), Real(0), -Real(0)}) {
		for (Real denominator : {Real(0), -Real(0), tolerance / 2, -tolerance / 2}) {
			const auto value = SafeDivide(numerator, denominator, SingularityPolicy::ReturnInf);
			CHECK(std::isinf(value));
			CHECK(std::signbit(value) == (std::signbit(numerator) != std::signbit(denominator)));
			CHECK(SafeDivide(numerator, denominator, SingularityPolicy::Clamp) ==
				numerator / std::copysign(tolerance, denominator));
		}
	}
	for (Real invalid : {Real(0), -Real(1), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity()})
		CHECK_THROWS_AS(SafeDivide(1, 2, SingularityPolicy::Throw, nullptr, invalid), NumericInputError);
	CHECK_THROWS_AS(SafeDivide(1, 2, static_cast<SingularityPolicy>(-1)), ArgumentError);
	CHECK_THROWS_AS(SafeDivide(1, std::numeric_limits<Real>::infinity()), NumericInputError);
	CHECK_THROWS_AS(SafeDivide(std::numeric_limits<Real>::quiet_NaN(), 0, SingularityPolicy::ReturnZero), NumericInputError);
	CHECK(std::isnan(SafeDivide(1, 0, SingularityPolicy::ReturnNaN)));
	CHECK(SafeDivide(1, 0, SingularityPolicy::ReturnZero) == 0);
	CHECK_THROWS_AS(SafeDivide(1, 0), DomainError);
}

TEST_CASE("Gram Schmidt is scale-aware and keeps documented dependent-column behavior", "[BaseUtilityContract]")
{
	for (Real scale : {largeScale, smallScale}) {
		const Matrix<Real> matrix{3, 3, {scale, scale, 0, scale, 0, 0, 0, scale, 0}};
		const auto basis = GramSchmidt(matrix);
		CHECK_THAT(basis.VectorFromColumn(0).NormL2(), Catch::Matchers::WithinRel(Real(1), Precision::EqualityTolerance));
		CHECK_THAT(basis.VectorFromColumn(1).NormL2(), Catch::Matchers::WithinRel(Real(1), Precision::EqualityTolerance));
		CHECK(std::abs(Utils::ScalarProduct(basis.VectorFromColumn(0), basis.VectorFromColumn(1))) <= Precision::OrthogonalityTolerance);
		CHECK(basis.VectorFromColumn(2).IsZero());
		const std::vector<VectorN<Real, 2>> vectors{{scale, scale}, {scale, 0}, {0, 0}, {2 * scale, 0}};
		const auto extracted = GramSchmidtVectors(vectors);
		REQUIRE(extracted.size() == 2);
		CHECK(IsOrthonormal(extracted));
	}
	const Real delta = std::sqrt(Precision::Eps);
	const auto nearlyParallel = GramSchmidt(Matrix<Real>{2, 2, {1, 1, 1, 1 + delta}});
	CHECK((nearlyParallel.Transposed() * nearlyParallel).IsIdentity(Precision::OrthogonalityTolerance));
	CHECK(GramSchmidt(Matrix<Real>{2, 2, {largeScale, 0, 0, smallScale}}).IsIdentity(Precision::EqualityTolerance));
}

TEST_CASE("Gram Schmidt rejects invalid tolerances and nonfinite data", "[BaseUtilityContract][BaseUtilityFinite]")
{
	for (Real tolerance : {Real(0), -Real(1), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(GramSchmidt(Matrix<Real>{1, 1, {0}}, tolerance), NumericInputError);
		CHECK_THROWS_AS(GramSchmidtVectors(std::vector<VectorN<Real, 1>>{}, tolerance), NumericInputError);
	}
	CHECK_THROWS_AS(GramSchmidt(Matrix<Real>{1, 1, {std::numeric_limits<Real>::quiet_NaN()}}), NumericInputError);
	CHECK_THROWS_AS(GramSchmidtVectors(std::vector<VectorN<Real, 1>>{{std::numeric_limits<Real>::infinity()}}), NumericInputError);
	CHECK(GramSchmidt(Matrix<Real>{1, 1, {0}}).VectorFromColumn(0).IsZero());
}
