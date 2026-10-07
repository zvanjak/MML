#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Optimization/Constraints/ProjectedGradient.h>
#endif

using namespace MML;
using namespace MML::Optimization;

namespace MML::Tests::Algorithms::Optimization
{
	class OffsetQuadratic2D : public IDifferentiableScalarFunction<2>
	{
		VectorN<Real, 2> _center;

	public:
		OffsetQuadratic2D(Real x, Real y)
			: _center{ x, y } { }

		Real operator()(const VectorN<Real, 2>& point) const override
		{
			Real dx = point[0] - _center[0];
			Real dy = point[1] - _center[1];
			return dx * dx + dy * dy;
		}

		void Gradient(const VectorN<Real, 2>& point, VectorN<Real, 2>& gradient) const override
		{
			gradient[0] = REAL(2.0) * (point[0] - _center[0]);
			gradient[1] = REAL(2.0) * (point[1] - _center[1]);
		}
	};

	TEST_CASE("ProjectedGradient - Finds interior constrained quadratic optimum", "[Optimization][ProjectedGradient][Bounds]")
	{
		OffsetQuadratic2D objective(REAL(0.25), -REAL(0.5));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		ProjectedGradientConfig config;
		config.gradient_tolerance = REAL(1e-9);
		config.max_iterations = 200;

		auto result = ProjectedGradientMinimize(objective, VectorN<Real, 2>{ REAL(0.9), REAL(0.9) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.projected_gradient_norm < REAL(1e-8));
		REQUIRE(result.xmin.IsEqualTo(VectorN<Real, 2>{ REAL(0.25), -REAL(0.5) }, REAL(1e-6)));
		REQUIRE(result.fmin < REAL(1e-12));
	}

	TEST_CASE("ProjectedGradient - Projects optimum outside box to active face", "[Optimization][ProjectedGradient][Bounds]")
	{
		OffsetQuadratic2D objective(REAL(2.0), -REAL(3.0));
		BoundConstraints bounds(Vector<Real>{ -REAL(1.0), -REAL(1.0) }, Vector<Real>{ REAL(1.0), REAL(1.0) });
		ProjectedGradientConfig config;
		config.gradient_tolerance = REAL(1e-9);
		config.max_iterations = 200;

		auto result = ProjectedGradientMinimize(objective, VectorN<Real, 2>{ REAL(0.0), REAL(0.0) }, bounds, config);

		REQUIRE(result.converged);
		REQUIRE(result.xmin.IsEqualTo(VectorN<Real, 2>{ REAL(1.0), -REAL(1.0) }, REAL(1e-6)));
		REQUIRE(bounds.IsUpperActive(result.xmin, 0));
		REQUIRE(bounds.IsLowerActive(result.xmin, 1));
	}

	TEST_CASE("ProjectedGradient - Inactive bounds match unconstrained solution", "[Optimization][ProjectedGradient][Bounds]")
	{
		OffsetQuadratic2D objective(-REAL(2.0), REAL(3.0));
		BoundConstraints bounds(Vector<Real>{ -REAL(10.0), -REAL(10.0) }, Vector<Real>{ REAL(10.0), REAL(10.0) });

		auto result = ProjectedGradientMinimize(objective, VectorN<Real, 2>{ REAL(5.0), -REAL(5.0) }, bounds);

		REQUIRE(result.converged);
		REQUIRE(result.xmin.IsEqualTo(VectorN<Real, 2>{ -REAL(2.0), REAL(3.0) }, REAL(1e-6)));
		REQUIRE(bounds.IsFree(result.xmin, 0));
		REQUIRE(bounds.IsFree(result.xmin, 1));
	}

	TEST_CASE("ProjectedGradient - Rejects invalid configuration", "[Optimization][ProjectedGradient][Bounds]")
	{
		ProjectedGradientConfig config;
		config.initial_step = REAL(0.0);
		REQUIRE_THROWS_AS(ProjectedGradient(config), MultidimOptimizationInputError);

		config.initial_step = REAL(1.0);
		config.contraction = REAL(1.0);
		REQUIRE_THROWS_AS(ProjectedGradient(config), MultidimOptimizationInputError);
	}
}
