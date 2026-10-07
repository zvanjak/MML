#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedNelderMead.h>
#endif

using namespace MML;
using namespace MML::Optimization;

namespace MML::Tests::Algorithms::Optimization
{
	class BoxNelderQuadratic2D : public IScalarFunction<2>
	{
		VectorN<Real, 2> _center;

	public:
		BoxNelderQuadratic2D(Real x, Real y)
			: _center{ x, y } { }

		Real operator()(const VectorN<Real, 2>& point) const override
		{
			Real dx = point[0] - _center[0];
			Real dy = point[1] - _center[1];
			return dx * dx + dy * dy;
		}
	};

	class BoxNelderRosenbrock2D : public IScalarFunction<2>
	{
	public:
		Real operator()(const VectorN<Real, 2>& point) const override
		{
			Real a = REAL(1.0) - point[0];
			Real b = point[1] - point[0] * point[0];
			return a * a + REAL(100.0) * b * b;
		}
	};

	TEST_CASE("BoxConstrainedNelderMead - Finds interior optimum from infeasible start", "[Optimization][BoxNelderMead][Bounds]")
	{
		BoxNelderQuadratic2D objective(REAL(0.25), -REAL(0.5));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = REAL(1e-10);
		config.max_iterations = 2000;
		config.initial_delta = REAL(0.5);

		auto result = BoxConstrainedNelderMeadMinimize(objective, VectorN<Real, 2>{ REAL(5.0), REAL(5.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(bounds.IsFeasible(result.xmin, REAL(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(REAL(0.25)).margin(REAL(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(-REAL(0.5)).margin(REAL(1e-5)));
		REQUIRE(result.fmin < REAL(1e-10));
	}

	TEST_CASE("BoxConstrainedNelderMead - Projects constrained optimum onto active bounds", "[Optimization][BoxNelderMead][Bounds]")
	{
		BoxNelderQuadratic2D objective(REAL(2.0), -REAL(3.0));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = REAL(1e-10);
		config.max_iterations = 2000;
		config.initial_delta = REAL(0.5);

		auto result = BoxConstrainedNelderMeadMinimize(objective, VectorN<Real, 2>{ REAL(0.0), REAL(0.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(bounds.IsFeasible(result.xmin, REAL(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(REAL(1.0)).margin(REAL(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(-REAL(1.0)).margin(REAL(1e-5)));
		REQUIRE(bounds.IsUpperActive(result.xmin, 0, REAL(1e-5)));
		REQUIRE(bounds.IsLowerActive(result.xmin, 1, REAL(1e-5)));
	}

	TEST_CASE("BoxConstrainedNelderMead - Supports custom repaired simplex deltas", "[Optimization][BoxNelderMead][Bounds]")
	{
		BoxNelderQuadratic2D objective(REAL(0.5), REAL(0.5));
		BoundConstraints bounds(Vector<Real>{ REAL(0.0), REAL(0.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = REAL(1e-10);
		config.max_iterations = 2000;
		Vector<Real> deltas{ REAL(0.25), REAL(0.25) };

		auto result = BoxConstrainedNelderMeadMinimize(objective, VectorN<Real, 2>{ REAL(1.0), REAL(1.0) }, bounds, deltas, config);

		REQUIRE(result.converged);
		REQUIRE(bounds.IsFeasible(result.xmin, REAL(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(REAL(0.5)).margin(REAL(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(REAL(0.5)).margin(REAL(1e-5)));
	}

	TEST_CASE("BoxConstrainedNelderMead - Handles bounded Rosenbrock valley", "[Optimization][BoxNelderMead][Bounds]")
	{
		BoxNelderRosenbrock2D objective;
		BoundConstraints bounds(Vector<Real>{ -REAL(2.0), -REAL(1.0) }, Vector<Real>{ REAL(2.0), REAL(3.0) });
		MultidimOptimizationConfig config;
		config.tolerance = REAL(1e-10);
		config.max_iterations = 5000;
		config.initial_delta = REAL(0.5);

		auto result = BoxConstrainedNelderMeadMinimize(objective, VectorN<Real, 2>{ -REAL(1.2), REAL(1.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(bounds.IsFeasible(result.xmin, REAL(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(REAL(1.0)).margin(REAL(1e-4)));
		REQUIRE(result.xmin[1] == Catch::Approx(REAL(1.0)).margin(REAL(1e-4)));
		REQUIRE(result.fmin < REAL(1e-8));
	}

	TEST_CASE("BoxConstrainedNelderMead - Rejects invalid configuration and deltas", "[Optimization][BoxNelderMead][Bounds]")
	{
		BoxNelderQuadratic2D objective(REAL(0.0), REAL(0.0));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		MultidimOptimizationConfig config;
		config.initial_delta = REAL(0.0);
		REQUIRE_THROWS_AS(BoxConstrainedNelderMead(config), MultidimOptimizationInputError);

		config.initial_delta = REAL(1.0);
		Vector<Real> deltas{ REAL(0.25), REAL(0.0) };
		BoxConstrainedNelderMead optimizer(config);
		REQUIRE_THROWS_AS(optimizer.Minimize(objective, VectorN<Real, 2>{ REAL(0.0), REAL(0.0) }, bounds, deltas), MultidimOptimizationInputError);
	}
}
