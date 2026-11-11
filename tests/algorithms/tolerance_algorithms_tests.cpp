// TOLERANCE_POLICY in mml/algorithms and mml/systems (task MinimalMathLibrary-wewi.4.4): the three
// chaos thresholds are one named constant, config defaults are named, sentinels are NaN, and the
// Defaults:: tolerance mirrors are deprecated in favour of Precision::.
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Function.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/CurveFitting.h>
#include <mml/algorithms/Optimization/LP/LPTypes.h>
#include <mml/algorithms/Optimization/Multidim/MultidimTypes.h>
#include <mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h>
#include <mml/algorithms/MatrixAlg/Subspaces.h>
#endif

#include <cmath>

using namespace MML;

namespace MML::Tests::Algorithms::ToleranceAlgorithmsTests
{
	TEST_CASE("ToleranceAlgorithms_tiers_and_named_config_defaults", "[tolerance][algorithms]")
	{
		static_assert(Precision::LooseTolerance < Precision::CoarseTolerance);
		static_assert(Precision::WellConditionedLimit < Precision::IllConditionedLimit);
		static_assert(Precision::BrentZeroEpsilon < Precision::SqrtEps);

		RootFinding::RootFindingConfig rc;
		REQUIRE(rc.tolerance == Precision::EqualityTolerance);

		Optimization::LPConfig lp;
		REQUIRE(lp.tolerance == Precision::EqualityTolerance);
		REQUIRE(lp.pivotTolerance == Precision::ZeroTolerance);
		REQUIRE(lp.bigM == Precision::BigMPenalty);

		Optimization::MultidimOptimizationConfig mo;
		REQUIRE(mo.tolerance == Precision::OptimizationGradientTolerance);

		Optimization::ALMConfig al;
		REQUIRE(al.feasibility_tolerance == Precision::LooseTolerance);
		REQUIRE(al.inner_tolerance == Precision::ModerateTolerance);
		REQUIRE(al.max_penalty == Precision::PenaltyParameterMax);
	}

	TEST_CASE("ToleranceAlgorithms_Ridders_without_the_minus_1e99_sentinel", "[tolerance][algorithms][RootFinding]")
	{
		RealFunction f(+[](Real x) { return x * x - 2.0; });
		auto r = RootFinding::FindRootRiddersDetailed(f, 1.0, 2.0, 1e-12);
		REQUIRE(r.converged);
		REQUIRE(r.root == Catch::Approx(std::sqrt(2.0)).epsilon(1e-11));
		// a root exactly at a bracket endpoint must still be found (the old sentinel compared against -1e98)
		RealFunction g(+[](Real x) { return x - 1.0; });
		auto r2 = RootFinding::FindRootRiddersDetailed(g, 1.0, 3.0, 1e-12);
		REQUIRE(IsUsable(r2.status));
		REQUIRE(r2.root == Catch::Approx(1.0).margin(1e-12));
	}

	TEST_CASE("ToleranceAlgorithms_curve_fitting_not_applicable_is_NaN", "[tolerance][algorithms]")
	{
		CurveFittingResult<Real> stats;
		REQUIRE(std::isnan(stats.adjusted_r_squared));
		REQUIRE(std::isnan(stats.condition_number));
		REQUIRE(stats.effective_rank == -1);   // integer "no index" stays an int
	}

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	TEST_CASE("ToleranceAlgorithms_Defaults_mirrors_equal_Precision", "[tolerance][deprecated]")
	{
		REQUIRE(Defaults::VectorIsEqualTolerance == Precision::VectorIsEqualTolerance);
		REQUIRE(Defaults::MatrixIsEqualTolerance == Precision::MatrixIsEqualTolerance);
		REQUIRE(Defaults::DefaultTolerance == Precision::DefaultTolerance);
		REQUIRE(Defaults::OrthogonalityTolerance == Precision::OrthogonalityTolerance);
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}
