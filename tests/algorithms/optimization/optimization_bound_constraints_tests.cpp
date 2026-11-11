#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/Constraints/BoundConstraints.h>
#endif

using namespace MML;
using namespace MML::Optimization;

namespace MML::Tests::Algorithms::Optimization
{
	TEST_CASE("BoundConstraints - Construction validates dimensions and bounds", "[Optimization][Bounds]")
	{
		Vector<Real> lower{ static_cast<MML::Real>(-1.0), static_cast<MML::Real>(0.0) };
		Vector<Real> upper{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0) };
		BoundConstraints bounds(lower, upper);

		REQUIRE(bounds.size() == 2);
		REQUIRE(bounds.lower(0) == static_cast<MML::Real>(-1.0));
		REQUIRE(bounds.upper(1) == static_cast<MML::Real>(2.0));

		REQUIRE_THROWS_AS(BoundConstraints(Vector<Real>{ static_cast<MML::Real>(0.0) }, upper), BoundConstraintsError);
		REQUIRE_THROWS_AS(BoundConstraints(Vector<Real>{ static_cast<MML::Real>(2.0) }, Vector<Real>{ static_cast<MML::Real>(1.0) }), BoundConstraintsError);
		REQUIRE_THROWS_AS(BoundConstraints(Vector<Real>{ std::numeric_limits<Real>::quiet_NaN() }, Vector<Real>{ static_cast<MML::Real>(1.0) }), BoundConstraintsError);
	}

	TEST_CASE("BoundConstraints - Project clamps dynamic vectors", "[Optimization][Bounds]")
	{
		BoundConstraints bounds(Vector<Real>{ static_cast<MML::Real>(-1.0), static_cast<MML::Real>(0.0), -std::numeric_limits<Real>::infinity() },
			Vector<Real>{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), std::numeric_limits<Real>::infinity() });

		Vector<Real> point{ static_cast<MML::Real>(-2.5), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0) };
		Vector<Real> projected = bounds.Project(point);

		REQUIRE(projected[0] == static_cast<MML::Real>(-1.0));
		REQUIRE(projected[1] == static_cast<MML::Real>(2.0));
		REQUIRE(projected[2] == static_cast<MML::Real>(4.0));
		REQUIRE(bounds.IsFeasible(projected));
		REQUIRE(bounds.MaxViolation(point) == static_cast<MML::Real>(1.5));
	}

	TEST_CASE("BoundConstraints - Project supports fixed-size vectors", "[Optimization][Bounds]")
	{
		VectorN<Real, 3> lower{ static_cast<MML::Real>(-1.0), static_cast<MML::Real>(-2.0), static_cast<MML::Real>(0.0) };
		VectorN<Real, 3> upper{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0) };
		BoundConstraints bounds(lower, upper);

		VectorN<Real, 3> point{ static_cast<MML::Real>(2.0), static_cast<MML::Real>(-3.0), static_cast<MML::Real>(4.0) };
		VectorN<Real, 3> projected = bounds.Project(point);

		REQUIRE(projected[0] == static_cast<MML::Real>(1.0));
		REQUIRE(projected[1] == -static_cast<MML::Real>(2.0));
		REQUIRE(projected[2] == static_cast<MML::Real>(4.0));
		REQUIRE(bounds.IsFeasible(projected));
		REQUIRE(bounds.MaxViolation(point) == static_cast<MML::Real>(1.0));
	}

	TEST_CASE("BoundConstraints - Active and free variables are detected", "[Optimization][Bounds]")
	{
		BoundConstraints bounds(Vector<Real>{ static_cast<MML::Real>(-1.0), static_cast<MML::Real>(0.0), -std::numeric_limits<Real>::infinity() },
			Vector<Real>{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), std::numeric_limits<Real>::infinity() });

		Vector<Real> point{ static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(5.0) };

		REQUIRE(bounds.IsLowerActive(point, 0));
		REQUIRE_FALSE(bounds.IsUpperActive(point, 0));
		REQUIRE(bounds.IsActive(point, 0));
		REQUIRE(bounds.IsFree(point, 1));
		REQUIRE(bounds.IsFree(point, 2));
	}

	TEST_CASE("BoundConstraints - ProjectInPlace and dimension checks", "[Optimization][Bounds]")
	{
		BoundConstraints bounds(2, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		Vector<Real> point{ -static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) };

		bounds.ProjectInPlace(point);
		REQUIRE(point[0] == static_cast<MML::Real>(0.0));
		REQUIRE(point[1] == static_cast<MML::Real>(1.0));

		Vector<Real> wrongDimension{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0) };
		REQUIRE_THROWS_AS(bounds.Project(wrongDimension), VectorDimensionError);
	}
}
