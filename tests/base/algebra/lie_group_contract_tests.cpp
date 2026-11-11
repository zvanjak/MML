#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Algebra/LieGroups.h>
#endif

#include <array>
#include <cmath>
#include <limits>

using namespace MML;
using namespace MML::Algebra;

namespace {
	constexpr Real RotationRoundoff = 64 * Precision::Eps;

	void CheckRotationValue(Real actual, Real expected) {
		CHECK(std::isfinite(actual));
		CHECK_THAT(actual, Catch::Matchers::WithinAbs(expected, RotationRoundoff));
	}

	template<int N>
	void CheckRotationMatrix(const MatrixNM<Real, N, N>& actual, const MatrixNM<Real, N, N>& expected) {
		for (int row = 0; row < N; ++row)
			for (int column = 0; column < N; ++column)
				CheckRotationValue(actual(row, column), expected(row, column));
	}

	MatrixNM<Real, 3, 3> RodriguesMatrix(const VectorN<Real, 3>& axis, Real angle) {
		const Real cosine = std::cos(angle), sine = std::sin(angle);
		const std::array<std::array<Real, 3>, 3> cross{{
			{0, -axis[2], axis[1]}, {axis[2], 0, -axis[0]}, {-axis[1], axis[0], 0}
		}};
		MatrixNM<Real, 3, 3> matrix;
		for (int row = 0; row < 3; ++row)
			for (int column = 0; column < 3; ++column)
				matrix(row, column) = (row == column ? cosine : 0)
					+ (1 - cosine) * axis[row] * axis[column] + sine * cross[row][column];
		return matrix;
	}

	void CheckCanonicalRotation(const SO3& rotation) {
		const auto& q = rotation.quaternion();
		CHECK(q.w() >= 0);
		CheckRotationValue(std::hypot(std::hypot(q.w(), q.x()), std::hypot(q.y(), q.z())), 1);
		const auto logarithm = rotation.log();
		const Real angle = std::hypot(logarithm[0], std::hypot(logarithm[1], logarithm[2]));
		CHECK(std::isfinite(angle));
		CHECK(angle <= Constants::PI + RotationRoundoff);
		CheckRotationMatrix(SO3::Exp(logarithm).matrix(), rotation.matrix());
	}
}

TEST_CASE("Quaternion Euler extraction preserves rotations at both gimbal poles",
	"[quaternion][euler][float-contract]") {
	for (bool xyz : {false, true})
		for (Real sign : {Real(-1), Real(1)})
			for (Real offset : {Real(0), Precision::EqualityTolerance / 4, Real(0.1)})
				for (Real yaw : {Real(-0.7), Real(0.3)})
					for (Real roll : {Real(-0.4), Real(0.2)}) {
						const Real pitch = sign * (Constants::PI / 2 - offset);
						const auto rotation = xyz ? Quaternion::FromEulerXYZ(roll, pitch, yaw) :
							Quaternion::FromEulerZYX(yaw, pitch, roll);
						const Real budget = RotationRoundoff + (offset < Precision::EqualityTolerance ? offset : Real(0));
						for (Real polarity : {Real(-1), Real(1)}) {
							CAPTURE(xyz, sign, offset, yaw, roll, polarity);
							const Quaternion q(polarity * rotation.w(), polarity * rotation.x(),
								polarity * rotation.y(), polarity * rotation.z());
							const auto angles = xyz ? q.ToEulerXYZ() : q.ToEulerZYX();
							const auto restored = xyz ? Quaternion::FromEulerXYZ(angles[0], angles[1], angles[2]) :
								Quaternion::FromEulerZYX(angles[0], angles[1], angles[2]);
							CHECK(IsWithinAbs(angles[1], pitch, budget));
							const auto expected = q.ToRotationMatrix(), actual = restored.ToRotationMatrix();
							for (int i = 0; i < 3; ++i)
								for (int j = 0; j < 3; ++j) {
									CHECK(std::isfinite(angles[i]));
									CHECK(IsWithinAbs(actual(i, j), expected(i, j), budget));
								}
						}
					}
}

TEST_CASE("SO3 principal logs and canonicalization are stable on both sides of pi",
	"[lie-group-contract][SO3]") {
	const std::array<VectorN<Real, 3>, 8> axes{{
		{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
		VectorN<Real, 3>{-1, 2, -3}.Normalized(), VectorN<Real, 3>{1, -2, 3}.Normalized()
	}};
	const Real offset = Precision::DefaultTolerance / 4;
	for (const auto& axis : axes) {
		for (Real angle : {Constants::PI / 3, Constants::PI - offset, Constants::PI + offset,
			Constants::PI - 4 * offset, Constants::PI + 4 * offset}) {
			CAPTURE(axis[0], axis[1], axis[2], angle);
			const Real principal = angle <= Constants::PI ? angle : angle - 2 * Constants::PI;
			for (const auto& rotation : {SO3::Exp(axis * angle), SO3::FromAxisAngle(axis, angle)}) {
				CheckCanonicalRotation(rotation);
				CheckRotationMatrix(rotation.matrix(), RodriguesMatrix(axis, angle));
				for (int i = 0; i < 3; ++i) CheckRotationValue(rotation.log()[i], axis[i] * principal);
				auto repeated = rotation;
				const auto original = rotation.quaternion();
				for (int iteration = 0; iteration < 4; ++iteration) {
					repeated = repeated.compose(SO3::Identity());
					CHECK(repeated.quaternion().w() >= 0);
					CheckRotationValue(repeated.quaternion().w(), original.w());
					CheckRotationValue(repeated.quaternion().x(), original.x());
					CheckRotationValue(repeated.quaternion().y(), original.y());
					CheckRotationValue(repeated.quaternion().z(), original.z());
				}
				const auto restored = SO3::FromMatrix(rotation.matrix());
				CheckCanonicalRotation(restored);
				CheckRotationMatrix(restored.matrix(), rotation.matrix());
			}
		}
	}
}

TEST_CASE("SO3 exact half-turn ties and adjacent angles preserve the hemisphere",
	"[lie-group-contract][SO3]") {
	for (const auto& axis : {VectorN<Real, 3>{1, 0, 0}, VectorN<Real, 3>{0, 1, 0},
		VectorN<Real, 3>{0, 0, 1}, VectorN<Real, 3>{Real(3) / 5, -Real(4) / 5, 0}}) {
		MatrixNM<Real, 3, 3> halfTurn;
		for (int row = 0; row < 3; ++row)
			for (int column = 0; column < 3; ++column)
				halfTurn(row, column) = 2 * axis[row] * axis[column] - (row == column ? 1 : 0);
		const auto rotation = SO3::FromMatrix(halfTurn);
		CHECK(rotation.quaternion().w() == 0);
		for (Real component : {rotation.quaternion().x(), rotation.quaternion().y(), rotation.quaternion().z()}) {
			if (component != 0) {
				CHECK(component > 0);
				break;
			}
		}
		CheckCanonicalRotation(rotation);
		CheckRotationMatrix(rotation.inverse().matrix(), halfTurn);
		CheckRotationMatrix(rotation.compose(rotation).matrix(), MatrixNM<Real, 3, 3>::Identity());
	}
	for (Real sign : {Real(-1), Real(1)}) {
		for (Real angle : {std::nextafter(Constants::PI, Real(0)), Constants::PI,
			std::nextafter(Constants::PI, std::numeric_limits<Real>::infinity())}) {
			const VectorN<Real, 3> axis{sign, 0, 0};
			CheckCanonicalRotation(SO3::Exp(axis * angle));
			CheckCanonicalRotation(SO3::FromAxisAngle(axis, angle));
		}
	}
}

TEST_CASE("SO3 distances and interpolation follow the shortest near-pi rotation",
	"[lie-group-contract][SO3]") {
	const Real offset = Precision::DefaultTolerance / 4;
	const auto left = SO3::Exp({Real(1) / 5, -Real(3) / 10, Real(1) / 10});
	for (Real sign : {Real(-1), Real(1)}) {
		const VectorN<Real, 3> axis{0, sign, 0};
		for (Real angle : {Constants::PI - offset, Constants::PI + offset}) {
			const auto relative = SO3::Exp(axis * angle);
			const auto right = left.compose(relative);
			const Real principal = angle <= Constants::PI ? angle : angle - 2 * Constants::PI;
			CheckRotationValue(left.geodesic_distance(right), std::abs(principal));
			CheckRotationValue(right.geodesic_distance(left), std::abs(principal));
			CHECK(left.geodesic_distance(right) <= Constants::PI + RotationRoundoff);
			for (Real t : {Real(0), Real(1) / 4, Real(1) / 2, Real(3) / 4, Real(1)}) {
				const auto interpolated = SO3::Slerp(SO3::Identity(), relative, t);
				CheckCanonicalRotation(interpolated);
				CheckRotationMatrix(interpolated.matrix(), RodriguesMatrix(axis, principal * t));
			}
		}
	}
}

TEST_CASE("SE2 strict matrix construction rejects every nonfinite entry", "[lie-group-contract][SE2]") {
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		for (int row = 0; row < 3; ++row)
			for (int column = 0; column < 3; ++column) {
				auto matrix = MatrixNM<Real, 3, 3>::Identity();
				matrix(row, column) = bad;
				CAPTURE(row, column, bad);
				CHECK_THROWS_AS(SE2::FromHomogeneousMatrix(matrix), DomainError);
			}
	}
}

TEST_CASE("SE3 strict matrix construction rejects every nonfinite entry", "[lie-group-contract][SE3]") {
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
		-std::numeric_limits<Real>::infinity()}) {
		for (int row = 0; row < 4; ++row)
			for (int column = 0; column < 4; ++column) {
				auto matrix = MatrixNM<Real, 4, 4>::Identity();
				matrix(row, column) = bad;
				CAPTURE(row, column, bad);
				CHECK_THROWS_AS(SE3::FromHomogeneousMatrix(matrix), DomainError);
			}
	}
}

TEST_CASE("Rotation and rigid matrix validation share a finite nonnegative tolerance contract",
	"[lie-group-contract][SE2][SE3][SO3]") {
	const auto planar = MatrixNM<Real, 3, 3>::Identity();
	const auto spatial = MatrixNM<Real, 4, 4>::Identity();
	for (Real bad : {Real(-1), std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(SO3::IsRotationMatrix(planar, bad), NumericInputError);
		CHECK_THROWS_AS(SO3::FromMatrix(planar, bad), NumericInputError);
		CHECK_THROWS_AS(SE2::FromHomogeneousMatrix(planar, bad), NumericInputError);
		CHECK_THROWS_AS(SE3::FromHomogeneousMatrix(spatial, bad), NumericInputError);
	}
	CHECK(SO3::IsRotationMatrix(planar, 0));
	CheckRotationMatrix(SO3::FromMatrix(planar, 0).matrix(), planar);
	CheckRotationMatrix(SE2::FromHomogeneousMatrix(planar, 0).homogeneous_matrix(), planar);
	CheckRotationMatrix(SE3::FromHomogeneousMatrix(spatial, 0).homogeneous_matrix(), spatial);
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		for (int row = 0; row < 3; ++row)
			for (int column = 0; column < 3; ++column) {
				auto matrix = planar;
				matrix(row, column) = bad;
				CHECK_FALSE(SO3::IsRotationMatrix(matrix));
				CHECK_THROWS_AS(SO3::FromMatrix(matrix), DomainError);
			}
	}
}

TEST_CASE("Strict rigid matrices preserve finite transforms and bottom-row tolerance boundaries",
	"[lie-group-contract][SE2][SE3]") {
	const SE2 planar(SO2(Real(1) / 3), {2, -3});
	const SE3 spatial(SO3::Exp({Real(1) / 3, -Real(1) / 5, Real(1) / 7}), {2, -3, 4});
	CheckRotationMatrix(SE2::FromHomogeneousMatrix(planar.homogeneous_matrix()).homogeneous_matrix(),
		planar.homogeneous_matrix());
	CheckRotationMatrix(SE3::FromHomogeneousMatrix(spatial.homogeneous_matrix()).homogeneous_matrix(),
		spatial.homogeneous_matrix());
	const Real tolerance = Precision::OrthogonalityTolerance;
	for (Real factor : {Real(1) / 2, Real(2)}) {
		auto matrix2 = MatrixNM<Real, 3, 3>::Identity();
		auto matrix3 = MatrixNM<Real, 4, 4>::Identity();
		matrix2(2, 0) = matrix3(3, 0) = factor * tolerance;
		if (factor < 1) {
			CHECK_NOTHROW(SE2::FromHomogeneousMatrix(matrix2, tolerance));
			CHECK_NOTHROW(SE3::FromHomogeneousMatrix(matrix3, tolerance));
		} else {
			CHECK_THROWS_AS(SE2::FromHomogeneousMatrix(matrix2, tolerance), DomainError);
			CHECK_THROWS_AS(SE3::FromHomogeneousMatrix(matrix3, tolerance), DomainError);
		}
	}
	auto matrix2 = MatrixNM<Real, 3, 3>::Identity();
	auto matrix3 = MatrixNM<Real, 4, 4>::Identity();
	matrix2(0, 2) = matrix3(0, 3) = std::numeric_limits<Real>::max();
	CHECK(SE2::FromHomogeneousMatrix(matrix2, 0).translation()[0] == std::numeric_limits<Real>::max());
	CHECK(SE3::FromHomogeneousMatrix(matrix3, 0).translation()[0] == std::numeric_limits<Real>::max());
	matrix2(2, 0) = matrix3(3, 0) = std::numeric_limits<Real>::denorm_min();
	CHECK_THROWS_AS(SE2::FromHomogeneousMatrix(matrix2, 0), DomainError);
	CHECK_THROWS_AS(SE3::FromHomogeneousMatrix(matrix3, 0), DomainError);
	matrix2 = MatrixNM<Real, 3, 3>::Identity();
	matrix3 = MatrixNM<Real, 4, 4>::Identity();
	matrix2(0, 0) = matrix3(0, 0) = -1;
	CHECK_THROWS_AS(SE2::FromHomogeneousMatrix(matrix2), DomainError);
	CHECK_THROWS_AS(SE3::FromHomogeneousMatrix(matrix3), DomainError);
}
