#include <catch2/catch_test_macros.hpp>

#include "../../TestPrecision.h"
#include <mml/base/Algebra_base.h>

#include <cmath>

namespace MML::Tests::Base::AlgebraTests
{
	constexpr Real GroupTolerance = Testing::Tol(REAL(1e-10), REAL(1e-6));

	bool MatrixNear(const MatrixNM<Real, 3, 3>& left, const MatrixNM<Real, 3, 3>& right, Real tolerance = GroupTolerance)
	{
		return left.IsEqualTo(right, tolerance);
	}

	bool VectorNear(const VectorN<Real, 3>& left, const VectorN<Real, 3>& right, Real tolerance = GroupTolerance)
	{
		return left.IsEqualTo(right, tolerance);
	}

	TEST_CASE("SO2 normalizes composes and follows shortest interpolation", "[Algebra][SO2][Base]")
	{
		const Algebra::SO2 quarterTurn(Constants::PI / REAL(2.0));
		const auto rotated = quarterTurn.apply({1, 0});

		REQUIRE(std::abs(rotated[0]) < GroupTolerance);
		REQUIRE(std::abs(rotated[1] - REAL(1.0)) < GroupTolerance);
		REQUIRE(std::abs(Algebra::SO2(REAL(3.0) * Constants::PI).angle() + Constants::PI) < GroupTolerance);
		REQUIRE(std::abs(quarterTurn.compose(quarterTurn.inverse()).angle()) < GroupTolerance);
		REQUIRE(std::abs(Algebra::SO2::Exp(quarterTurn.log()).angle() - quarterTurn.angle()) < GroupTolerance);
		REQUIRE(std::abs(Algebra::SO2::Interpolate(Algebra::SO2(REAL(0.9) * Constants::PI),
			Algebra::SO2(-REAL(0.9) * Constants::PI), REAL(0.5)).angle() + Constants::PI) < GroupTolerance);
	}

	TEST_CASE("SO3 hat vee exp and log are mutually consistent", "[Algebra][SO3][Base]")
	{
		const VectorN<Real, 3> omega{REAL(0.2), REAL(-0.3), REAL(0.4)};
		const auto rotation = Algebra::SO3::Exp(omega);

		REQUIRE(VectorNear(Algebra::SO3::Vee(Algebra::SO3::Hat(omega)), omega));
		REQUIRE(VectorNear(rotation.log(), omega));
		REQUIRE(VectorNear(Algebra::SO3::Exp({REAL(1e-12), REAL(-2e-12), REAL(3e-12)}).log(),
			VectorN<Real, 3>{REAL(1e-12), REAL(-2e-12), REAL(3e-12)}, REAL(1e-14)));
		REQUIRE_THROWS_AS(Algebra::SO3::Vee(MatrixNM<Real, 3, 3>::Identity()), DomainError);
	}

	TEST_CASE("SO3 logarithm remains stable near pi", "[Algebra][SO3][Base]")
	{
		const VectorN<Real, 3> axis = VectorN<Real, 3>{1, 2, -1}.Normalized();
		const VectorN<Real, 3> omega = axis * (Constants::PI - REAL(1e-8));
		const auto recovered = Algebra::SO3::Exp(omega).log();

		REQUIRE(VectorNear(recovered, omega, REAL(1e-9)));
		REQUIRE(MatrixNear(Algebra::SO3::Exp(recovered).matrix(), Algebra::SO3::Exp(omega).matrix()));
	}

	TEST_CASE("SO3 axis-angle preserves norms and composes", "[Algebra][SO3][Base]")
	{
		const auto rotation = Algebra::SO3::FromAxisAngle({0, 0, 1}, Constants::PI / REAL(2.0));
		const VectorN<Real, 3> vector{1, 0, 0};
		const auto rotated = rotation.apply(vector);

		REQUIRE(VectorNear(rotated, {0, 1, 0}));
		REQUIRE(std::abs(rotated.NormL2() - vector.NormL2()) < GroupTolerance);
		REQUIRE(MatrixNear(rotation.compose(rotation.inverse()).matrix(), MatrixNM<Real, 3, 3>::Identity()));
		REQUIRE(Algebra::SO3::IsRotationMatrix(rotation.matrix()));
		REQUIRE_THROWS_AS(Algebra::SO3::FromAxisAngle({0, 0, 0}, REAL(1.0)), ArgumentError);
	}

	TEST_CASE("SO3 matrix construction is strict and projection repairs drift", "[Algebra][SO3][Base]")
	{
		const auto original = Algebra::SO3::Exp({REAL(0.3), REAL(-0.1), REAL(0.2)});
		const auto roundTrip = Algebra::SO3::FromMatrix(original.matrix());
		MatrixNM<Real, 3, 3> drifted = original.matrix();
		drifted(0, 0) += REAL(1e-4);
		drifted(1, 2) -= REAL(2e-4);

		REQUIRE(MatrixNear(roundTrip.matrix(), original.matrix()));
		REQUIRE_THROWS_AS(Algebra::SO3::FromMatrix(drifted), DomainError);
		const auto projected = Algebra::SO3::Project(drifted);
		REQUIRE(Algebra::SO3::IsRotationMatrix(projected.matrix()));
		REQUIRE(projected.geodesic_distance(original) < REAL(1e-3));

		const MatrixNM<Real, 3, 3> reflection{{1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
		REQUIRE_THROWS_AS(Algebra::SO3::FromMatrix(reflection), DomainError);
		REQUIRE_THROWS_AS((Algebra::SO3::Project({{1, 1, 0}, {0, 0, 1}, {0, 0, 0}})), DomainError);
	}

	TEST_CASE("SO3 geodesic distance and slerp follow shortest rotation", "[Algebra][SO3][Base]")
	{
		const auto identity = Algebra::SO3::Identity();
		const auto halfTurn = Algebra::SO3::FromAxisAngle({0, 0, 1}, Constants::PI);
		const auto midpoint = Algebra::SO3::Slerp(identity, halfTurn, REAL(0.5));

		REQUIRE(std::abs(identity.geodesic_distance(halfTurn) - Constants::PI) < GroupTolerance);
		REQUIRE(std::abs(identity.geodesic_distance(midpoint) - Constants::PI / REAL(2.0)) < GroupTolerance);
		const auto midpointVector = midpoint.apply({1, 0, 0});
		REQUIRE(std::abs(midpointVector[0]) < GroupTolerance);
		REQUIRE(std::abs(std::abs(midpointVector[1]) - REAL(1.0)) < GroupTolerance);
		REQUIRE(std::abs(midpointVector[2]) < GroupTolerance);
	}
}