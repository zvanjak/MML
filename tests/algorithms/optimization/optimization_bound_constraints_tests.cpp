#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#elif defined MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/Constraints/BoundConstraints.h>
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedNelderMead.h>
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedPowell.h>
#include <mml/algorithms/Optimization/Constraints/ProjectedGradient.h>
#endif

#include <cmath>
#include <limits>
#include <string>
#include <utility>

using namespace MML;
using namespace MML::Optimization;

namespace MML::Tests::Algorithms::Optimization
{
	class CountingBoundsQuadratic final : public IDifferentiableScalarFunction<2>
	{
	public:
		mutable int calls = 0;
		mutable int gradients = 0;
		static Real Value(const VectorN<Real, 2>& point)
		{
			return (point[0] - 2) * (point[0] - 2) + (point[1] + 3) * (point[1] + 3);
		}
		Real operator()(const VectorN<Real, 2>& point) const override
		{
			++calls;
			return Value(point);
		}
		void Gradient(const VectorN<Real, 2>& point, VectorN<Real, 2>& gradient) const override
		{
			++gradients;
			gradient[0] = 2 * (point[0] - 2);
			gradient[1] = 2 * (point[1] + 3);
		}
	};

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

	TEST_CASE("Bound feasibility contracts: nonfinite coordinates are errors even without finite bounds",
		"[Optimization][Bounds][bounds-contract]")
	{
		const Real inf = std::numeric_limits<Real>::infinity();
		for (int kind = 0; kind < 3; ++kind) {
			const auto bounds = kind == 0 ? BoundConstraints(3, 0, 1)
				: kind == 1 ? BoundConstraints::Unbounded(3)
				: BoundConstraints(Vector<Real>{0, -inf, 0}, Vector<Real>{1, 1, inf});
			for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), inf, -inf}) {
				for (int index = 0; index < 3; ++index) {
					Vector<Real> dynamic{Real(1) / 2, Real(1) / 2, Real(1) / 2};
					VectorN<Real, 3> fixed{Real(1) / 2, Real(1) / 2, Real(1) / 2};
					dynamic[index] = fixed[index] = invalid;
					CHECK_THROWS_AS(bounds.IsFeasible(dynamic), BoundConstraintsError);
					CHECK_THROWS_AS(bounds.MaxViolation(dynamic), BoundConstraintsError);
					CHECK_THROWS_AS(bounds.IsFeasible(fixed), BoundConstraintsError);
					CHECK_THROWS_AS(bounds.MaxViolation(fixed), BoundConstraintsError);
				}
			}
		}
	}

	TEST_CASE("Bound feasibility contracts: invalid bound data and empty query domains are explicit errors",
		"[Optimization][Bounds][bounds-contract]")
	{
		const Real inf = std::numeric_limits<Real>::infinity(), nan = std::numeric_limits<Real>::quiet_NaN();
		for (const auto& limits : {std::pair<Real, Real>{nan, Real(1)}, {Real(0), nan}, {Real(2), Real(1)}, {inf, inf}, {-inf, -inf}}) {
			CHECK_THROWS_AS(BoundConstraints(1, limits.first, limits.second), BoundConstraintsError);
			CHECK_THROWS_AS(BoundConstraints(Vector<Real>{limits.first}, Vector<Real>{limits.second}), BoundConstraintsError);
			CHECK_THROWS_AS(BoundConstraints(VectorN<Real, 1>{limits.first}, VectorN<Real, 1>{limits.second}), BoundConstraintsError);
		}
		const BoundConstraints empty;
		CHECK_THROWS_AS(empty.IsFeasible(Vector<Real>{}), BoundConstraintsError);
		CHECK_THROWS_AS(empty.MaxViolation(Vector<Real>{}), BoundConstraintsError);
		const BoundConstraints valid(1, 0, 1);
		CHECK_THROWS_AS(valid.IsFeasible(Vector<Real>{}), VectorDimensionError);
		CHECK_THROWS_AS(valid.MaxViolation(Vector<Real>{}), VectorDimensionError);
		CHECK_THROWS_AS(valid.IsFeasible(VectorN<Real, 2>{0, 0}), VectorDimensionError);
		CHECK_THROWS_AS(valid.MaxViolation(VectorN<Real, 2>{0, 0}), VectorDimensionError);
	}

	TEST_CASE("Bound feasibility contracts: tolerances are finite nonnegative and retain exact comparisons",
		"[Optimization][Bounds][bounds-contract]")
	{
		const BoundConstraints bounds(1, 0, 1);
		for (Real invalid : {Real(-1), std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
			CHECK_THROWS_AS(bounds.IsFeasible(Vector<Real>{Real(1) / 2}, invalid), BoundConstraintsError);
			CHECK_THROWS_AS(bounds.IsFeasible(VectorN<Real, 1>{Real(1) / 2}, invalid), BoundConstraintsError);
		}
		CHECK(bounds.IsFeasible(Vector<Real>{1}, 0));
		CHECK(bounds.IsFeasible(VectorN<Real, 1>{0}, 0));
		const Real tolerance = Precision::VectorIsEqualTolerance;
		CHECK(bounds.IsFeasible(Vector<Real>{1 + tolerance / 2}, tolerance));
		CHECK_FALSE(bounds.IsFeasible(Vector<Real>{1 + 2 * tolerance}, tolerance));
		CHECK(bounds.IsFeasible(VectorN<Real, 1>{1 + tolerance / 2}, tolerance));
		CHECK_FALSE(bounds.IsFeasible(VectorN<Real, 1>{1 + 2 * tolerance}, tolerance));
	}

	TEST_CASE("Bound feasibility contracts: late invalid coordinates cannot hide behind overflowed finite violations",
		"[Optimization][Bounds][bounds-contract]")
	{
		const Real maximum = std::numeric_limits<Real>::max();
		const Real nan = std::numeric_limits<Real>::quiet_NaN();
		const BoundConstraints bounds(Vector<Real>{-maximum, 0}, Vector<Real>{-maximum, 1});
		CHECK_THROWS_AS(bounds.MaxViolation(Vector<Real>{maximum, nan}), BoundConstraintsError);
		CHECK_THROWS_AS(bounds.IsFeasible(Vector<Real>{maximum, nan}), BoundConstraintsError);
		CHECK_THROWS_AS(bounds.MaxViolation(VectorN<Real, 2>{maximum, nan}), BoundConstraintsError);
		CHECK_THROWS_AS(bounds.IsFeasible(VectorN<Real, 2>{maximum, nan}), BoundConstraintsError);
		CHECK(std::isinf(bounds.MaxViolation(Vector<Real>{maximum, 0})));
		CHECK(std::isinf(bounds.MaxViolation(VectorN<Real, 2>{maximum, 0})));
		CHECK_FALSE(bounds.IsFeasible(Vector<Real>{maximum, 0}));
		CHECK_FALSE(bounds.IsFeasible(VectorN<Real, 2>{maximum, 0}));
	}

	TEST_CASE("Bound feasibility contracts: projection fixed faces and correctly oriented infinite endpoints are preserved",
		"[Optimization][Bounds][bounds-contract]")
	{
		const Real inf = std::numeric_limits<Real>::infinity(), maximum = std::numeric_limits<Real>::max();
		const BoundConstraints bounds(Vector<Real>{0, -inf, 2}, Vector<Real>{1, inf, 2});
		Vector<Real> dynamic{-2, maximum, 5};
		VectorN<Real, 3> fixed{-2, maximum, 5};
		const auto projectedDynamic = bounds.Project(dynamic);
		const auto projectedFixed = bounds.Project(fixed);
		CHECK_FALSE(bounds.IsFeasible(dynamic));
		CHECK_FALSE(bounds.IsFeasible(fixed));
		CHECK(bounds.MaxViolation(dynamic) == 3);
		CHECK(bounds.MaxViolation(fixed) == 3);
		bounds.ProjectInPlace(dynamic);
		bounds.ProjectInPlace(fixed);
		for (int index = 0; index < 3; ++index) {
			CHECK(dynamic[index] == projectedDynamic[index]);
			CHECK(fixed[index] == projectedFixed[index]);
		}
		CHECK(dynamic[0] == 0);
		CHECK(dynamic[1] == maximum);
		CHECK(dynamic[2] == 2);
		CHECK(bounds.IsFeasible(dynamic));
		CHECK(bounds.IsFeasible(fixed));
		CHECK(bounds.MaxViolation(dynamic) == 0);
		CHECK(bounds.MaxViolation(fixed) == 0);
		const auto unbounded = BoundConstraints::Unbounded(2);
		CHECK(unbounded.IsFeasible(Vector<Real>{-maximum, maximum}));
		CHECK(unbounded.MaxViolation(VectorN<Real, 2>{-maximum, maximum}) == 0);
		const BoundConstraints lowerOnly(1, 0, inf), upperOnly(1, -inf, 1);
		CHECK(lowerOnly.IsFeasible(VectorN<Real, 1>{maximum}));
		CHECK(upperOnly.IsFeasible(Vector<Real>{-maximum}));
	}

	TEST_CASE("Bound feasibility contracts: invalid-coordinate diagnostics name the inspected coordinate",
		"[Optimization][Bounds][bounds-contract]")
	{
		const BoundConstraints bounds(2, 0, 1);
		const Real nan = std::numeric_limits<Real>::quiet_NaN();
		bool rejected = false;
		try {
			const auto unexpected = bounds.MaxViolation(Vector<Real>{0, nan});
			(void)unexpected;
		} catch (const BoundConstraintsError& error) {
			rejected = true;
			CHECK(std::string(error.what()).find("point[1]") != std::string::npos);
			CHECK(std::string(error.what()).find("finite") != std::string::npos);
		}
		CHECK(rejected);
	}

	TEST_CASE("Bound feasibility contracts: optimizer entry points independently reject invalid starts before callbacks",
		"[Optimization][Bounds][bounds-contract]")
	{
		const BoundConstraints bounds(2, 0, 1);
		for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity(),
			-std::numeric_limits<Real>::infinity()}) {
			for (int index = 0; index < 2; ++index) {
				CountingBoundsQuadratic function;
				VectorN<Real, 2> start{Real(1) / 2, Real(1) / 2};
				start[index] = invalid;
				CHECK_THROWS_AS(ProjectedGradient().Minimize(function, start, bounds), MultidimOptimizationInputError);
				CHECK_THROWS_AS(ProjectedGradientMinimize(function, start, bounds), MultidimOptimizationInputError);
				CHECK_THROWS_AS(BoxConstrainedPowell().Minimize(function, start, bounds), MultidimOptimizationInputError);
				CHECK_THROWS_AS(BoxConstrainedPowellMinimize(function, start, bounds), MultidimOptimizationInputError);
				CHECK_THROWS_AS(BoxConstrainedNelderMead().Minimize(function, start, bounds), MultidimOptimizationInputError);
				CHECK_THROWS_AS(BoxConstrainedNelderMeadMinimize(function, start, bounds), MultidimOptimizationInputError);
				const Vector<Real> deltas{Real(1) / 2, Real(1) / 2};
				CHECK_THROWS_AS(BoxConstrainedNelderMead().Minimize(function, start, bounds, deltas), MultidimOptimizationInputError);
				CHECK_THROWS_AS(BoxConstrainedNelderMeadMinimize(function, start, bounds, deltas), MultidimOptimizationInputError);
				CHECK(function.calls == 0);
				CHECK(function.gradients == 0);
			}
		}
	}

	TEST_CASE("Bound feasibility contracts: box wrappers retain projected returned-objective recomputation",
		"[Optimization][Bounds][bounds-contract]")
	{
		const BoundConstraints bounds(2, 0, 1);
		const VectorN<Real, 2> start{2, -2};
		MultidimOptimizationConfig config;
		config.max_iterations = 1;
		for (int method = 0; method < 2; ++method) {
			CountingBoundsQuadratic function;
			const auto result = method == 0 ? BoxConstrainedPowellMinimize(function, start, bounds, config)
				: BoxConstrainedNelderMeadMinimize(function, start, bounds, config);
			CHECK(bounds.IsFeasible(result.xmin));
			CHECK(result.function_evaluations == function.calls);
			VectorN<Real, 2> returned;
			for (int index = 0; index < 2; ++index) returned[index] = result.xmin[index];
			CHECK(result.fmin == CountingBoundsQuadratic::Value(returned));
			CHECK(std::isfinite(result.fmin));
		}
	}
}
