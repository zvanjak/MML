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
		BoxPowellQuadratic2D objective(static_cast<MML::Real>(0.25), -static_cast<MML::Real>(0.5));
		BoundConstraints bounds(Vector<Real>{ -static_cast<MML::Real>(1.0), -static_cast<MML::Real>(1.0) }, Vector<Real>{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = static_cast<MML::Real>(1e-8);
		config.max_iterations = 300;

		auto result = BoxConstrainedPowellMinimize(objective, VectorN<Real, 2>{ static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(bounds.IsFeasible(result.xmin, static_cast<MML::Real>(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(static_cast<MML::Real>(0.25)).margin(static_cast<MML::Real>(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(-static_cast<MML::Real>(0.5)).margin(static_cast<MML::Real>(1e-5)));
		REQUIRE(result.fmin < static_cast<MML::Real>(1e-10));
	}

	TEST_CASE("BoxConstrainedPowell - Projects constrained optimum onto active bounds", "[Optimization][BoxPowell][Bounds]")
	{
		BoxPowellQuadratic2D objective(static_cast<MML::Real>(2.0), -static_cast<MML::Real>(3.0));
		BoundConstraints bounds(Vector<Real>{ -static_cast<MML::Real>(1.0), -static_cast<MML::Real>(1.0) }, Vector<Real>{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0) });
		MultidimOptimizationConfig config;
		config.tolerance = static_cast<MML::Real>(1e-8);
		config.max_iterations = 300;

		auto result = BoxConstrainedPowellMinimize(objective, VectorN<Real, 2>{ static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(bounds.IsFeasible(result.xmin, static_cast<MML::Real>(1e-10)));
		REQUIRE(result.xmin[0] == Catch::Approx(static_cast<MML::Real>(1.0)).margin(static_cast<MML::Real>(1e-5)));
		REQUIRE(result.xmin[1] == Catch::Approx(-static_cast<MML::Real>(1.0)).margin(static_cast<MML::Real>(1e-5)));
		REQUIRE(bounds.IsUpperActive(result.xmin, 0, static_cast<MML::Real>(1e-5)));
		REQUIRE(bounds.IsLowerActive(result.xmin, 1, static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("BoxConstrainedPowell - Rejects invalid configuration", "[Optimization][BoxPowell][Bounds]")
	{
		MultidimOptimizationConfig config;
		config.tolerance = static_cast<MML::Real>(0.0);
		REQUIRE_THROWS_AS(BoxConstrainedPowell(config), MultidimOptimizationInputError);

		config.tolerance = static_cast<MML::Real>(1e-8);
		config.max_iterations = 0;
		REQUIRE_THROWS_AS(BoxConstrainedPowell(config), MultidimOptimizationInputError);
	}
}
