#include <catch2/catch_test_macros.hpp>

#include "../../TestPrecision.h"
#include <mml/base/Algebra_base.h>

namespace MML::Tests::Base::AlgebraTests
{
	constexpr Real RigidMotionTolerance = Testing::Tol(REAL(1e-12), REAL(1e-6));

	struct Body2;
	struct World2;
	struct Body3;
	struct World3;

	TEST_CASE("SE2 distinguishes points from vectors and composes", "[Algebra][SE2][Base]")
	{
		const Algebra::SE2 transform(Algebra::SO2(Constants::PI / REAL(2.0)), {2, 3});
		REQUIRE(transform.apply_point({1, 0}).IsEqualTo({2, 4}, RigidMotionTolerance));
		REQUIRE(transform.apply_vector({1, 0}).IsEqualTo({0, 1}, RigidMotionTolerance));
		REQUIRE(transform.inverse().apply_point(transform.apply_point({-1, 4})).IsEqualTo({-1, 4}, RigidMotionTolerance));

		const Algebra::SE2 right(Algebra::SO2(REAL(0.2)), {-1, 2});
		const VectorN<Real, 2> point{3, -2};
		REQUIRE(transform.compose(right).apply_point(point).IsEqualTo(
			transform.apply_point(right.apply_point(point)), RigidMotionTolerance));
	}

	TEST_CASE("SE2 homogeneous and typed actions round trip", "[Algebra][SE2][Base]")
	{
		const Algebra::SE2 transform(Algebra::SO2(REAL(0.4)), {2, -3});
		const auto restored = Algebra::SE2::FromHomogeneousMatrix(transform.homogeneous_matrix());
		Point<Real, 2, Body2> point{1, 2};
		TangentVector<2, Body2> vector{1, 2};

		REQUIRE(restored.homogeneous_matrix().IsEqualTo(transform.homogeneous_matrix(), RigidMotionTolerance));
		const auto typedPoint = transform.apply_point<Body2, World2>(point);
		const auto typedVector = transform.apply_vector<Body2, World2>(vector);
		REQUIRE(typedPoint.coordinates().IsEqualTo(transform.apply_point(point.coordinates()), RigidMotionTolerance));
		REQUIRE(typedVector.components().IsEqualTo(transform.apply_vector(vector.components()), RigidMotionTolerance));

		auto malformed = transform.homogeneous_matrix();
		malformed(0, 0) = REAL(2.0);
		REQUIRE_THROWS_AS(Algebra::SE2::FromHomogeneousMatrix(malformed), DomainError);
	}

	TEST_CASE("SE3 distinguishes points from vectors and composes", "[Algebra][SE3][Base]")
	{
		const Algebra::SE3 transform(
			Algebra::SO3::FromAxisAngle({0, 0, 1}, Constants::PI / REAL(2.0)), {2, 3, 4});
		REQUIRE(transform.apply_point({1, 0, 0}).IsEqualTo({2, 4, 4}, RigidMotionTolerance));
		REQUIRE(transform.apply_vector({1, 0, 0}).IsEqualTo({0, 1, 0}, RigidMotionTolerance));
		REQUIRE(transform.inverse().apply_point(transform.apply_point({-1, 4, 2})).IsEqualTo({-1, 4, 2}, RigidMotionTolerance));

		const Algebra::SE3 right(Algebra::SO3::Exp({REAL(0.1), REAL(-0.2), REAL(0.3)}), {-1, 2, 1});
		const VectorN<Real, 3> point{3, -2, 1};
		REQUIRE(transform.compose(right).apply_point(point).IsEqualTo(
			transform.apply_point(right.apply_point(point)), RigidMotionTolerance));
	}

	TEST_CASE("SE3 homogeneous and typed actions round trip", "[Algebra][SE3][Base]")
	{
		const Algebra::SE3 transform(Algebra::SO3::Exp({REAL(0.2), REAL(0.1), REAL(-0.3)}), {2, -3, 1});
		const auto restored = Algebra::SE3::FromHomogeneousMatrix(transform.homogeneous_matrix());
		Point<Real, 3, Body3> point{1, 2, 3};
		TangentVector<3, Body3> vector{1, 2, 3};

		REQUIRE(restored.homogeneous_matrix().IsEqualTo(transform.homogeneous_matrix(), RigidMotionTolerance));
		const auto typedPoint = transform.apply_point<Body3, World3>(point);
		const auto typedVector = transform.apply_vector<Body3, World3>(vector);
		REQUIRE(typedPoint.coordinates().IsEqualTo(transform.apply_point(point.coordinates()), RigidMotionTolerance));
		REQUIRE(typedVector.components().IsEqualTo(transform.apply_vector(vector.components()), RigidMotionTolerance));

		auto malformed = transform.homogeneous_matrix();
		malformed(3, 0) = REAL(1.0);
		REQUIRE_THROWS_AS(Algebra::SE3::FromHomogeneousMatrix(malformed), DomainError);
	}
}