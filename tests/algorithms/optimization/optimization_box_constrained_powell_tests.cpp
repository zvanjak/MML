#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/Constraints/BoxConstrainedPowell.h>
#endif

using namespace MML;
using namespace MML::Optimization;

namespace MML::Tests::Algorithms::Optimization
{
	class BoxPowellQuadratic2D : public IScalarFunction<2>
	{
		VectorN<Real, 2> _center;

	public:
		BoxPowellQuadratic2D(Real x, Real y)
			: _center{ x, y } { }

		Real operator()(const VectorN<Real, 2>& point) const override
		{
			Real dx = point[0] - _center[0];
			Real dy = point[1] - _center[1];
			return dx * dx + dy * dy;
		}
	};

	TEST_CASE("BoxConstrainedPowell - Finds interior optimum from infeasible start", "[Optimization][BoxPowell][Bounds]")
	{
		BoxPowellQuadratic2D objective(REAL(0.25), -REAL(0.5));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = REAL(1e-8);
		config.max_iterations = 300;

		auto result = BoxConstrainedPowellMinimize(objective, VectorN<Real, 2>{ REAL(5.0), REAL(5.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(bounds.IsFeasible(result.xmin, REAL(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(REAL(0.25)).margin(REAL(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(-REAL(0.5)).margin(REAL(1e-5)));
		REQUIRE(result.fmin < REAL(1e-10));
	}

	TEST_CASE("BoxConstrainedPowell - Projects constrained optimum onto active bounds", "[Optimization][BoxPowell][Bounds]")
	{
		BoxPowellQuadratic2D objective(REAL(2.0), -REAL(3.0));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = REAL(1e-8);
		config.max_iterations = 300;

		auto result = BoxConstrainedPowellMinimize(objective, VectorN<Real, 2>{ REAL(0.0), REAL(0.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(bounds.IsFeasible(result.xmin, REAL(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(REAL(1.0)).margin(REAL(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(-REAL(1.0)).margin(REAL(1e-5)));
		REQUIRE(bounds.IsUpperActive(result.xmin, 0, REAL(1e-5)));
		REQUIRE(bounds.IsLowerActive(result.xmin, 1, REAL(1e-5)));
	}

	TEST_CASE("BoxConstrainedPowell - Rejects invalid configuration", "[Optimization][BoxPowell][Bounds]")
	{
		MultidimOptimizationConfig config;
		config.tolerance = REAL(0.0);
		REQUIRE_THROWS_AS(BoxConstrainedPowell(config), MultidimOptimizationInputError);

		config.tolerance = REAL(1e-8);
		config.max_iterations = 0;
		REQUIRE_THROWS_AS(BoxConstrainedPowell(config), MultidimOptimizationInputError);
	}
}
