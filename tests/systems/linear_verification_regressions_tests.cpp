#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/LinearSystem.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

using namespace MML;
using namespace MML::Systems;
using Catch::Matchers::WithinRel;

TEST_CASE("Linear verification keeps representable ratios across RHS units", "[verification-regression]") {
	const Real large = std::numeric_limits<Real>::max() / 2;
	const Real tiny = std::numeric_limits<Real>::min();
	const Matrix<Real> identity = Matrix<Real>::Identity(4);
	for (const Real magnitude : {tiny, Real(1), large}) {
		const LinearSystem<Real> system(identity, Vector<Real>(4, magnitude));
		const Vector<Real> wrong{0, magnitude, magnitude, magnitude};
		const auto result = system.Verify(wrong);
		CHECK_THAT(system.RelativeResidual(wrong), WithinRel(Real(.5), Precision::EqualityTolerance));
		CHECK_THAT(result.relativeResidual, WithinRel(Real(.5), Precision::EqualityTolerance));
		CHECK_THAT(result.backwardError, WithinRel(Real(1) / std::sqrt(Real(3)), Precision::EqualityTolerance));
		CHECK_THAT(result.estimatedForwardError, WithinRel(Real(.5), Precision::EqualityTolerance));
		CHECK_FALSE(result.isAccurate);
		const auto exact = system.Verify(Vector<Real>(4, magnitude));
		CHECK(exact.relativeResidual == 0);
		CHECK(exact.isAccurate);
	}
	const Real nearMaximum = std::numeric_limits<Real>::max() * Real(.75);
	const LinearSystem<Real> overflowingNorm(identity, Vector<Real>(4, nearMaximum));
	const auto result = overflowingNorm.Verify(Vector<Real>{0, nearMaximum, nearMaximum, nearMaximum});
	CHECK(std::isfinite(result.absoluteResidual));
	CHECK_THAT(result.relativeResidual, WithinRel(Real(.5), Precision::EqualityTolerance));
	CHECK_FALSE(result.isAccurate);
}

TEST_CASE("Linear verification scales products before forming error denominators", "[verification-regression]") {
	const Real large = std::sqrt(std::numeric_limits<Real>::max());
	const Matrix<Real> coefficient(1, 1, large);
	const LinearSystem<Real> system(coefficient, Vector<Real>{1});
	const auto result = system.Verify(Vector<Real>{large * Real(2)});
	CHECK(std::isfinite(result.backwardError));
	CHECK_THAT(result.backwardError, WithinRel(Real(1), Precision::EqualityTolerance));
	CHECK_FALSE(result.isAccurate);
}

TEST_CASE("Linear verification defines homogeneous zero RHS behavior", "[verification-regression]") {
	const LinearSystem<Real> system(Matrix<Real>::Identity(2), Vector<Real>{0, 0});
	const auto exact = system.Verify(Vector<Real>{0, 0});
	CHECK(exact.relativeResidual == 0);
	CHECK(exact.backwardError == 0);
	CHECK(exact.isAccurate);
	const auto wrong = system.Verify(Vector<Real>{1, 0});
	CHECK(std::isinf(wrong.relativeResidual));
	CHECK_FALSE(wrong.isAccurate);
}

TEST_CASE("Linear verification accepts every finite candidate for a zero operator and RHS",
	"[verification-regression]") {
	const LinearSystem<Real> system(Matrix<Real>(2, 2, Real(0)), Vector<Real>{0, 0});
	for (const Vector<Real>& candidate : {Vector<Real>{0, 0}, Vector<Real>{Real(-.25), 1}}) {
		CAPTURE(candidate[0], candidate[1]);
		CHECK(system.ResidualNorm(candidate) == 0);
		CHECK(system.RelativeResidual(candidate) == 0);
		const auto result = system.Verify(candidate);
		CHECK(result.absoluteResidual == 0);
		CHECK(result.relativeResidual == 0);
		CHECK(result.backwardError == 0);
		CHECK(result.isAccurate);
	}
}

TEST_CASE("Linear verification rejects invalid tolerances modes and nonfinite data", "[verification-regression]") {
	const Matrix<Real> identity = Matrix<Real>::Identity(2);
	const LinearSystem<Real> system(identity, Vector<Real>{1, 1});
	for (const Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(),
							 std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(system.Verify(Vector<Real>{1, 1}, invalid), NumericInputError);
	}
	for (const Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(system.Verify(Vector<Real>{invalid, 1}), NumericInputError);
		const LinearSystem<Real> invalidRHS(identity, Vector<Real>{invalid, 1});
		CHECK_THROWS_AS(invalidRHS.RelativeResidual(Vector<Real>{1, 1}), NumericInputError);
		CHECK_THROWS_AS(invalidRHS.ResidualNorm(Vector<Real>{1, 1}), NumericInputError);
		CHECK_THROWS_AS(invalidRHS.Verify(Vector<Real>{1, 1}), NumericInputError);
		Matrix<Real> coefficient = identity;
		coefficient(0, 0) = invalid;
		const LinearSystem<Real> invalidMatrix(coefficient, Vector<Real>{1, 1});
		CHECK_THROWS_AS(invalidMatrix.Verify(Vector<Real>{1, 1}), NumericInputError);
		CHECK_THROWS_AS(invalidMatrix.RelativeResidual(Vector<Real>{1, 1}), NumericInputError);
		CHECK_THROWS_AS(invalidMatrix.ResidualNorm(Vector<Real>{1, 1}), NumericInputError);
	}
	CHECK_THROWS_AS(system.Verify(Vector<Real>{1}), VectorDimensionError);
	const LinearSystem<Real> multiple(identity, Matrix<Real>::Identity(2));
	CHECK_THROWS_AS(multiple.Verify(Vector<Real>{1, 1}), InvalidStateError);
	CHECK_THROWS_AS(multiple.RelativeResidual(Vector<Real>{1, 1}), InvalidStateError);
	CHECK_THROWS_AS(multiple.ResidualNorm(Vector<Real>{1, 1}), InvalidStateError);
}

TEST_CASE("Linear verification retains ratios when absolute residual entries overflow", "[verification-regression]") {
	const Real magnitude = std::numeric_limits<Real>::max() * Real(.75);
	const LinearSystem<Real> system(Matrix<Real>::Identity(1), Vector<Real>{-magnitude});
	const auto result = system.Verify(Vector<Real>{magnitude});
	CHECK(std::isinf(result.absoluteResidual));
	CHECK_THAT(result.relativeResidual, WithinRel(Real(2), Precision::EqualityTolerance));
	CHECK_THAT(result.backwardError, WithinRel(Real(2), Precision::EqualityTolerance));
	CHECK_FALSE(result.isAccurate);
}

TEST_CASE("Linear verification does not silently erase unresolvable cancellation", "[verification-regression]") {
	const Real large = std::numeric_limits<Real>::max() / 2;
	const Real tiny = std::numeric_limits<Real>::min();
	const Matrix<Real> coefficient(1, 3, {large, -large, tiny});
	const LinearSystem<Real> system(coefficient, Vector<Real>{tiny});
	CHECK_THROWS_AS(system.RelativeResidual(Vector<Real>{1, 1, 1}), MatrixNumericalError);
}

TEST_CASE("Linear verification matches an independent non-diagonal residual oracle", "[verification-regression]") {
	// Leave room for coefficients/RHS up to 4*scale; retain the original double scale range.
	const int scaleExponent = std::min(500, std::numeric_limits<Real>::max_exponent - 4);
	for (const int exponent : {-scaleExponent, 0, scaleExponent}) {
		CAPTURE(exponent);
		const Real scale = std::ldexp(Real(1), exponent);
		REQUIRE(std::isnormal(scale));
		REQUIRE(std::isfinite(4 * scale));
		const Matrix<Real> coefficient(2, 2, {2 * scale, scale, scale, 3 * scale});
		const LinearSystem<Real> system(coefficient, Vector<Real>{3 * scale, 4 * scale});
		if constexpr (std::is_same_v<Real, float>) {
			const auto& rhs = system.getRHS();
			const Real naiveSquaredNorm = rhs[0] * rhs[0] + rhs[1] * rhs[1];
			if (exponent < 0)
				CHECK(naiveSquaredNorm == 0);
			else if (exponent > 0)
				CHECK(std::isinf(naiveSquaredNorm));
		}
		for (const Real delta : {Real(-.25), Real(.125), Real(1)}) {
			CAPTURE(delta);
			const Vector<Real> candidate{1 + delta, 1 - delta};
			const Real expected = std::abs(delta) / std::sqrt(Real(5));
			const Real expectedAbsolute = std::sqrt(Real(5)) * std::abs(delta) * scale;
			REQUIRE(std::isnormal(expectedAbsolute));
			CHECK_THAT(system.ResidualNorm(candidate), WithinRel(expectedAbsolute, Precision::EqualityTolerance));
			CHECK_THAT(system.RelativeResidual(candidate), WithinRel(expected, Precision::EqualityTolerance));
			const auto result = system.Verify(candidate);
			CHECK_THAT(result.absoluteResidual, WithinRel(expectedAbsolute, Precision::EqualityTolerance));
			CHECK_THAT(result.relativeResidual, WithinRel(expected, Precision::EqualityTolerance));
			const Real expectedBackward = std::sqrt(Real(5)) * std::abs(delta) /
				(4 * std::hypot(1 + delta, 1 - delta));
			CHECK_THAT(result.backwardError, WithinRel(expectedBackward, Precision::EqualityTolerance));
			CHECK_FALSE(result.isAccurate);
		}
	}
}
