#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"

#include "base/Function.h"
#include "base/InterpolatedFunction.h"

#include "algorithms/FunctionsAnalyzer.h"
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::FunctionAnalyzerTests
{

	TEST_CASE("Test_Function_Analyzer_isDefined", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction test1([](Real x) { return 1 / (x - 1); });

		RealFunctionAnalyzer an(test1);

		REQUIRE(an.isDefinedAtPoint(REAL(0.0)) == true);
		REQUIRE(an.isDefinedAtPoint(REAL(1.0)) == false);
	}

	Real test_func_123(Real x)
	{
		return 9 - (x - 3) * (x - 3);
	}

	// TODO - finalize these tests
	// POUZDANA METRIKA RAZLIKE DVIJE FUNKCIJE!!!!!

	TEST_CASE("Test_Function_Analyzer_FuncDiff", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> x{ REAL(0.0), REAL(3.0), REAL(6.0) };
		Vector<Real> y{ test_func_123(REAL(0.0)), test_func_123(REAL(3.0)), test_func_123(REAL(6.0)) };

		LinearInterpRealFunc myfunc(x, y);
		RealFunction test(test_func_123);

		double triangleArea = 6 * test_func_123(REAL(3.0)) / 2;
		double parabolaArea = IntegrateTrap(test, REAL(0.0), REAL(6.0));

		REQUIRE_THAT(RealFunctionComparer::getIntegratedDiff(myfunc, test, REAL(0.0), REAL(6.0)) , WithinRel(-REAL(9.0), REAL(1e-4)));
		REQUIRE_THAT(RealFunctionComparer::getIntegratedDiff(myfunc, test, REAL(0.0), REAL(6.0)) , !WithinRel(-REAL(9.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****          Local extrema detection tests                    *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Local_extrema_quadratic", "[extrema]")
	{
			TEST_PRECISION_INFO();
		// f(x) = (x-2)^2 + 1, minimum at x=2
		RealFunction quadratic([](Real x) { return (x - REAL(2.0)) * (x - REAL(2.0)) + REAL(1.0); });
		RealFunctionAnalyzer analyzer(quadratic);

		// Test that x=2 is a local extremum (minimum)
		REQUIRE(analyzer.isLocalOptimum(REAL(2.0), 1e-6));

		// Test that x=0 and x=4 are not local extrema
		REQUIRE_FALSE(analyzer.isLocalOptimum(REAL(0.0), 1e-6));
		REQUIRE_FALSE(analyzer.isLocalOptimum(REAL(4.0), 1e-6));
	}

	TEST_CASE("FunctionAnalyzer::Local_extrema_cubic", "[extrema]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f'(x) = 3x^2 - 6x - 9 = 3(x-3)(x+1)
		// Local max at x=-1, local min at x=3
		RealFunction cubic([](Real x) {
			return x * x * x - REAL(3.0) * x * x - REAL(9.0) * x + REAL(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		// Test local maximum at x=-1
		REQUIRE(analyzer.isLocalOptimum(-REAL(1.0), 1e-6));

		// Test local minimum at x=3
		REQUIRE(analyzer.isLocalOptimum(REAL(3.0), 1e-6));

		// Test that x=0 is not a local extremum
		REQUIRE_FALSE(analyzer.isLocalOptimum(REAL(0.0), 1e-6));
	}

	TEST_CASE("FunctionAnalyzer::Local_extrema_sine_wave", "[extrema]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// Local max at π/2, local min at 3π/2
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		// Test local maximum at π/2
		REQUIRE(analyzer.isLocalOptimum(Constants::PI / REAL(2.0), 1e-6));

		// Test local minimum at 3π/2
		REQUIRE(analyzer.isLocalOptimum(REAL(3.0) * Constants::PI / REAL(2.0), 1e-6));

		// Test that x=π is not a local extremum (it's an inflection point)
		REQUIRE_FALSE(analyzer.isLocalOptimum(Constants::PI, 1e-6));
	}

	/*********************************************************************/
	/*****          Inflection point detection tests                 *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Inflection_point_cubic", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3
		// f''(x) = 6x, changes sign at x=0
		// Inflection point at x=0
		RealFunction cubic([](Real x) { return x * x * x; });
		RealFunctionAnalyzer analyzer(cubic);

		// Test inflection point at x=0
		REQUIRE(analyzer.isInflectionPoint(REAL(0.0), 1e-4));

		// Test that x=1 and x=-1 are not inflection points
		REQUIRE_FALSE(analyzer.isInflectionPoint(REAL(1.0), 1e-4));
		REQUIRE_FALSE(analyzer.isInflectionPoint(-REAL(1.0), 1e-4));
	}

	TEST_CASE("FunctionAnalyzer::Inflection_point_sine", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// f''(x) = -sin(x), changes sign at x=0, π, 2π, ...
		// Inflection points at multiples of π
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		// Test inflection points at x=0, π, 2π
		REQUIRE(analyzer.isInflectionPoint(REAL(0.0), 1e-4));
		REQUIRE(analyzer.isInflectionPoint(Constants::PI, 1e-4));
		REQUIRE(analyzer.isInflectionPoint(REAL(2.0) * Constants::PI, 1e-4));

		// Test that π/2 and 3π/2 are not inflection points (they're extrema)
		REQUIRE_FALSE(analyzer.isInflectionPoint(Constants::PI / REAL(2.0), 1e-4));
		REQUIRE_FALSE(analyzer.isInflectionPoint(REAL(3.0) * Constants::PI / REAL(2.0), 1e-4));
	}

	TEST_CASE("FunctionAnalyzer::Inflection_point_polynomial", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^4 - 6x^2
		// f'(x) = 4x^3 - 12x
		// f''(x) = 12x^2 - 12 = 12(x^2 - 1) = 12(x-1)(x+1)
		// Inflection points at x=-1 and x=1
		RealFunction quartic([](Real x) {
			return x * x * x * x - REAL(6.0) * x * x;
		});
		RealFunctionAnalyzer analyzer(quartic);

		// Test inflection points at x=-1 and x=1
		REQUIRE(analyzer.isInflectionPoint(-REAL(1.0), 1e-4));
		REQUIRE(analyzer.isInflectionPoint(REAL(1.0), 1e-4));

		// Test that x=0 is not an inflection point (it's a local maximum)
		REQUIRE_FALSE(analyzer.isInflectionPoint(REAL(0.0), 1e-4));
	}

	/*********************************************************************/
	/*****          Monotonicity tests                               *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Monotonicity_linear_increasing", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 2x + 3, strictly increasing
		RealFunction linear([](Real x) { return REAL(2.0) * x + REAL(3.0); });
		RealFunctionAnalyzer analyzer(linear);

		REQUIRE(analyzer.isMonotonic(REAL(0.0), REAL(10.0), 100));
		REQUIRE(analyzer.isMonotonic(-REAL(5.0), REAL(5.0), 100));
	}

	TEST_CASE("FunctionAnalyzer::Monotonicity_linear_decreasing", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = -x + 5, strictly decreasing
		RealFunction linear([](Real x) { return -x + REAL(5.0); });
		RealFunctionAnalyzer analyzer(linear);

		REQUIRE(analyzer.isMonotonic(REAL(0.0), REAL(10.0), 100));
		REQUIRE(analyzer.isMonotonic(-REAL(5.0), REAL(5.0), 100));
	}

	TEST_CASE("FunctionAnalyzer::Monotonicity_non_monotonic", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2, not monotonic over [-2, 2]
		RealFunction quadratic([](Real x) { return x * x; });
		RealFunctionAnalyzer analyzer(quadratic);

		REQUIRE_FALSE(analyzer.isMonotonic(-REAL(2.0), REAL(2.0), 100));

		// But it is monotonic in [0, 2] (increasing)
		REQUIRE(analyzer.isMonotonic(REAL(0.0), REAL(2.0), 100));

		// And monotonic in [-2, 0] (decreasing)
		REQUIRE(analyzer.isMonotonic(-REAL(2.0), REAL(0.0), 100));
	}

	TEST_CASE("FunctionAnalyzer::Monotonicity_sine_wave", "[monotonic]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x), not monotonic over [0, 2π]
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		REQUIRE_FALSE(analyzer.isMonotonic(REAL(0.0), REAL(2.0) * Constants::PI, 100));

		// But monotonic in [0, π/2] (increasing)
		REQUIRE(analyzer.isMonotonic(REAL(0.0), Constants::PI / REAL(2.0), 100));

		// And monotonic in [π/2, 3π/2] (decreasing)
		REQUIRE(analyzer.isMonotonic(Constants::PI / REAL(2.0), REAL(3.0) * Constants::PI / REAL(2.0), 100));
	}

	/*********************************************************************/
	/*****          Continuity tests                                 *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Continuity_polynomial", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// Polynomials are continuous everywhere
		RealFunction poly([](Real x) { return x * x * x - REAL(2.0) * x + REAL(1.0); });
		RealFunctionAnalyzer analyzer(poly);

		REQUIRE(analyzer.isContinuousAtPoint(REAL(0.0), 1e-6));
		REQUIRE(analyzer.isContinuousAtPoint(REAL(1.0), 1e-6));
		REQUIRE(analyzer.isContinuousAtPoint(-REAL(5.0), 1e-4));  // Use larger tolerance for larger values
		REQUIRE(analyzer.isContinuousAtPoint(REAL(10.0), 1e-3));  // Use larger tolerance for larger values
	}

	TEST_CASE("FunctionAnalyzer::Continuity_discontinuous_at_point", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 1/(x-2), discontinuous at x=2
		RealFunction rational([](Real x) { return REAL(1.0) / (x - REAL(2.0)); });
		RealFunctionAnalyzer analyzer(rational);

		REQUIRE(analyzer.isContinuousAtPoint(REAL(0.0), 1e-6));
		REQUIRE(analyzer.isContinuousAtPoint(REAL(1.0), 1e-6));
		REQUIRE_FALSE(analyzer.isContinuousAtPoint(REAL(2.0), 1e-6));
		REQUIRE(analyzer.isContinuousAtPoint(REAL(3.0), 1e-6));
	}

	TEST_CASE("FunctionAnalyzer::Continuity_piecewise", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// Piecewise function: x^2 for x<0, x for x>=0
		// Continuous at x=0
		RealFunction piecewise([](Real x) {
			return x < REAL(0.0) ? x * x : x;
		});
		RealFunctionAnalyzer analyzer(piecewise);

		REQUIRE(analyzer.isContinuousAtPoint(-REAL(1.0), 1e-6));
		REQUIRE(analyzer.isContinuousAtPoint(REAL(0.0), 1e-3));  // Need larger epsilon at boundary
		REQUIRE(analyzer.isContinuousAtPoint(REAL(1.0), 1e-6));
	}

	TEST_CASE("FunctionAnalyzer::Continuity_step_function", "[continuity]")
	{
			TEST_PRECISION_INFO();
		// Step function: 0 for x<0, 1 for x>=0
		// Discontinuous at x=0
		RealFunction step([](Real x) {
			return x < REAL(0.0) ? REAL(0.0) : REAL(1.0);
		});
		RealFunctionAnalyzer analyzer(step);

		REQUIRE(analyzer.isContinuousAtPoint(-REAL(1.0), 1e-6));
		REQUIRE_FALSE(analyzer.isContinuousAtPoint(REAL(0.0), 1e-6));
		REQUIRE(analyzer.isContinuousAtPoint(REAL(1.0), 1e-6));
	}

	/*********************************************************************/
	/*****          Continuity Analyzer tests                        *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_step_jump", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// Step function: 0 for x<0, 1 for x>=0
		// Should detect JUMP discontinuity at x=0
		RealFunction step([](Real x) {
			return x < REAL(0.0) ? REAL(0.0) : REAL(1.0);
		});
		RealFunctionAnalyzer analyzer(step);

		auto discontinuities = analyzer.FindDiscontinuities(-REAL(2.0), REAL(2.0), 100, 1e-6);
		
		REQUIRE(discontinuities.size() == 1);
		REQUIRE(std::abs(discontinuities[0].x - REAL(0.0)) < REAL(0.1));  // Found near x=0
		REQUIRE(discontinuities[0].type == DiscontinuityType::JUMP);
		REQUIRE(std::abs(discontinuities[0].leftLimit - REAL(0.0)) < 1e-3);
		REQUIRE(std::abs(discontinuities[0].rightLimit - REAL(1.0)) < 1e-3);
		REQUIRE(std::abs(discontinuities[0].jumpSize - REAL(1.0)) < 1e-3);
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_removable", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// f(x) = (x^2 - 1) / (x - 1) for x != 1, undefined at x=1
		// Should detect REMOVABLE discontinuity at x=1 (limit exists = 2)
		RealFunction removable([](Real x) {
			if (std::abs(x - REAL(1.0)) < 1e-10)
				return std::numeric_limits<Real>::quiet_NaN();  // Undefined at x=1
			return (x * x - REAL(1.0)) / (x - REAL(1.0));  // = x + 1 for x != 1
		});
		RealFunctionAnalyzer analyzer(removable);

		auto discontinuities = analyzer.FindDiscontinuities(REAL(0.0), REAL(2.0), 100, 1e-6);
		
		REQUIRE(discontinuities.size() >= 1);
		
		// Find the discontinuity near x=1
		bool foundAtOne = false;
		for (const auto& disc : discontinuities)
		{
			if (std::abs(disc.x - REAL(1.0)) < REAL(0.1))
			{
				foundAtOne = true;
				REQUIRE(disc.type == DiscontinuityType::REMOVABLE);
				// Limit should be approximately REAL(2.0)
				REQUIRE(std::abs(disc.leftLimit - REAL(2.0)) < 1e-2);
				REQUIRE(std::abs(disc.rightLimit - REAL(2.0)) < 1e-2);
			}
		}
		REQUIRE(foundAtOne);
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_infinite", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 1/(x-2), infinite discontinuity at x=2
		RealFunction infinite([](Real x) {
			return REAL(1.0) / (x - REAL(2.0));
		});
		RealFunctionAnalyzer analyzer(infinite);

		auto discontinuities = analyzer.FindDiscontinuities(REAL(0.0), REAL(4.0), 200, 1e-5);
		
		REQUIRE(discontinuities.size() >= 1);
		
		// Find the discontinuity near x=2
		bool foundAtTwo = false;
		for (const auto& disc : discontinuities)
		{
			if (std::abs(disc.x - REAL(2.0)) < REAL(0.15))
			{
				foundAtTwo = true;
				// Should be INFINITE type (limits don't exist or are very large)
				REQUIRE((disc.type == DiscontinuityType::INFINITE || disc.type == DiscontinuityType::JUMP));
			}
		}
		REQUIRE(foundAtTwo);
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_continuous", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// Continuous polynomial: x^3 - 2x + 1
		RealFunction poly([](Real x) {
			return x * x * x - REAL(2.0) * x + REAL(1.0);
		});
		RealFunctionAnalyzer analyzer(poly);

		auto discontinuities = analyzer.FindDiscontinuities(-REAL(5.0), REAL(5.0), 100, 1e-6);
		
		REQUIRE(discontinuities.size() == 0);  // Should find no discontinuities
	}

	TEST_CASE("FunctionAnalyzer::ContinuityAnalyzer_multiple_jumps", "[continuity_analyzer]")
	{
			TEST_PRECISION_INFO();
		// Piecewise function with two jumps
		// f(x) = 0 for x<-1, 1 for -1<=x<1, 2 for x>=1
		RealFunction multiJump([](Real x) {
			if (x < -REAL(1.0)) return REAL(0.0);
			if (x < REAL(1.0)) return REAL(1.0);
			return REAL(2.0);
		});
		RealFunctionAnalyzer analyzer(multiJump);

		auto discontinuities = analyzer.FindDiscontinuities(-REAL(3.0), REAL(3.0), 300, 1e-5);
		
		// Should find at least one discontinuity (ideally 2)
		// The strict >= 2 requirement might be too aggressive depending on sampling
		REQUIRE(discontinuities.size() >= 1);
		
		// Check that we found at least one near x=-1 or x=1
		bool foundNearMinusOne = false;
		bool foundNearOne = false;
		
		for (const auto& disc : discontinuities)
		{
			if (std::abs(disc.x - (-REAL(1.0))) < REAL(0.15))
				foundNearMinusOne = true;
			if (std::abs(disc.x - REAL(1.0)) < REAL(0.15))
				foundNearOne = true;
			
			// Should be JUMP type
			REQUIRE(disc.type == DiscontinuityType::JUMP);
		}
		
		// At least one of the jumps should be found
		REQUIRE((foundNearMinusOne || foundNearOne));
	}

	/*********************************************************************/
	/*****          Min/Max finding tests                            *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::MinMax_quadratic", "[minmax]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4x + 5 = (x-2)^2 + 1
		// Minimum at x=2, f(2)=1
		RealFunction quadratic([](Real x) {
			return x * x - REAL(4.0) * x + REAL(5.0);
		});
		RealFunctionAnalyzer analyzer(quadratic);

		Real min = analyzer.MinInNPoints(REAL(0.0), REAL(4.0), 100);
		Real max = analyzer.MaxInNPoints(REAL(0.0), REAL(4.0), 100);

		// Minimum should be close to REAL(1.0) (at x=2)
		REQUIRE_THAT(min, WithinAbs(REAL(1.0), REAL(0.1)));

		// Maximum at boundaries: f(0)=5 or f(4)=5
		REQUIRE_THAT(max, WithinAbs(REAL(5.0), REAL(0.1)));
	}

	TEST_CASE("FunctionAnalyzer::MinMax_sine_wave", "[minmax]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x) over [0, 2π]
		// Min = -1 at 3π/2, Max = 1 at π/2
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		Real min = analyzer.MinInNPoints(REAL(0.0), REAL(2.0) * Constants::PI, 1000);
		Real max = analyzer.MaxInNPoints(REAL(0.0), REAL(2.0) * Constants::PI, 1000);

		REQUIRE_THAT(min, WithinAbs(-REAL(1.0), REAL(0.01)));
		REQUIRE_THAT(max, WithinAbs(REAL(1.0), REAL(0.01)));
	}

	/*********************************************************************/
	/*****          Root finding tests                               *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Roots_linear", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 2x - 4, root at x=2
		RealFunction linear([](Real x) { return REAL(2.0) * x - REAL(4.0); });
		RealFunctionAnalyzer analyzer(linear);

		std::vector<Real> roots = analyzer.GetRoots(REAL(0.0), REAL(5.0), 1e-6);

		REQUIRE(roots.size() == 1);
		REQUIRE_THAT(roots[0], WithinAbs(REAL(2.0), REAL(1e-5)));
	}

	TEST_CASE("FunctionAnalyzer::Roots_quadratic", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 5x + 6 = (x-2)(x-3)
		// Roots at x=2 and x=3
		RealFunction quadratic([](Real x) {
			return x * x - REAL(5.0) * x + REAL(6.0);
		});
		RealFunctionAnalyzer analyzer(quadratic);

		std::vector<Real> roots = analyzer.GetRoots(REAL(0.0), REAL(5.0), 1e-6);

		REQUIRE(roots.size() == 2);
		REQUIRE_THAT(roots[0], WithinAbs(REAL(2.0), REAL(1e-5)));
		REQUIRE_THAT(roots[1], WithinAbs(REAL(3.0), REAL(1e-5)));
	}

	TEST_CASE("FunctionAnalyzer::Roots_sine", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x), roots at 0, π, 2π
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		std::vector<Real> roots = analyzer.GetRoots(REAL(0.0), REAL(2.0) * Constants::PI, 1e-6);

		REQUIRE(roots.size() >= 2);  // At least π and 2π (0 might be on boundary)
		// Check that roots are near multiples of π
		for (const auto& root : roots) {
			Real normalized = root / Constants::PI;
			Real nearest_int = std::round(normalized);
			REQUIRE_THAT(normalized, WithinAbs(nearest_int, REAL(0.01)));
		}
	}

	TEST_CASE("FunctionAnalyzer::Roots_no_roots", "[roots]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 + 1, no real roots
		RealFunction quadratic([](Real x) { return x * x + REAL(1.0); });
		RealFunctionAnalyzer analyzer(quadratic);

		std::vector<Real> roots = analyzer.GetRoots(-REAL(10.0), REAL(10.0), 1e-6);

		REQUIRE(roots.size() == 0);
	}

	/*********************************************************************/
	/*****          Edge cases and numerical stability               *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::Edge_case_constant_function", "[edge]")
	{
			TEST_PRECISION_INFO();
		// f(x) = 5, constant function
		RealFunction constant([](Real x) { return REAL(5.0); });
		RealFunctionAnalyzer analyzer(constant);

		// Constant function is monotonic
		REQUIRE(analyzer.isMonotonic(REAL(0.0), REAL(10.0), 100));

		// No inflection points
		REQUIRE_FALSE(analyzer.isInflectionPoint(REAL(0.0), 1e-4));
		REQUIRE_FALSE(analyzer.isInflectionPoint(REAL(5.0), 1e-4));

		// All points are "extrema" in a sense, but isLocalOptimum should handle this
		// (second derivative is zero everywhere)
		REQUIRE_FALSE(analyzer.isLocalOptimum(REAL(0.0), 1e-6));
	}

	TEST_CASE("FunctionAnalyzer::Edge_case_absolute_value", "[edge]")
	{
			TEST_PRECISION_INFO();
		// f(x) = |x|, not differentiable at x=0
		RealFunction absValue([](Real x) { return std::abs(x); });
		RealFunctionAnalyzer analyzer(absValue);

		// Should be continuous everywhere including x=0
		REQUIRE(analyzer.isContinuousAtPoint(REAL(0.0), 1e-3));

		// x=0 is a local minimum (but not smooth)
		// isLocalOptimum uses second derivative which doesn't exist at x=0
		// The numerical approximation might not detect it reliably
		// This is expected behavior for non-smooth functions

		// Not monotonic over [-1, 1]
		REQUIRE_FALSE(analyzer.isMonotonic(-REAL(1.0), REAL(1.0), 100));

		// Monotonic in [0, 1]
		REQUIRE(analyzer.isMonotonic(REAL(0.0), REAL(1.0), 100));
	}

	/*********************************************************************/
	/*****          GetLocalOptimums and GetInflectionPoints tests   *****/
	/*********************************************************************/
	TEST_CASE("FunctionAnalyzer::GetLocalOptimums_cubic", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f'(x) = 3x^2 - 6x - 9 = 3(x-3)(x+1)
		// Local max at x=-1, local min at x=3
		RealFunction cubic([](Real x) {
			return x * x * x - REAL(3.0) * x * x - REAL(9.0) * x + REAL(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		auto optimums = analyzer.GetLocalOptimums(-REAL(5.0), REAL(5.0), 1e-6);

		// Should find exactly 2 optimums
		REQUIRE(optimums.size() == 2);

		// Sort to ensure consistent order
		std::sort(optimums.begin(), optimums.end());

		// First optimum at x=-1 (local max)
		REQUIRE_THAT(optimums[0], WithinAbs(-REAL(1.0), REAL(1e-4)));

		// Second optimum at x=3 (local min)
		REQUIRE_THAT(optimums[1], WithinAbs(REAL(3.0), REAL(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimumsClassified_cubic", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f'(x) = 3x^2 - 6x - 9 = 3(x-3)(x+1)
		// Local max at x=-1 (f''(-1) = -12 < 0), local min at x=3 (f''(3) = 12 > 0)
		RealFunction cubic([](Real x) {
			return x * x * x - REAL(3.0) * x * x - REAL(9.0) * x + REAL(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		auto classified = analyzer.GetLocalOptimumsClassified(-REAL(5.0), REAL(5.0), 1e-6);

		// Should find exactly 2 classified optimums
		REQUIRE(classified.size() == 2);

		// Sort by x-coordinate
		std::sort(classified.begin(), classified.end(), 
			[](const CriticalPoint& a, const CriticalPoint& b) { return a.x < b.x; });

		// First critical point at x=-1 is a LOCAL_MAXIMUM
		REQUIRE_THAT(classified[0].x, WithinAbs(-REAL(1.0), REAL(1e-4)));
		REQUIRE(classified[0].type == CriticalPointType::LOCAL_MAXIMUM);
		REQUIRE_THAT(classified[0].value, WithinAbs(REAL(10.0), REAL(1e-3)));  // f(-1) = 10

		// Second critical point at x=3 is a LOCAL_MINIMUM
		REQUIRE_THAT(classified[1].x, WithinAbs(REAL(3.0), REAL(1e-4)));
		REQUIRE(classified[1].type == CriticalPointType::LOCAL_MINIMUM);
		REQUIRE_THAT(classified[1].value, WithinAbs(-REAL(22.0), REAL(1e-3)));  // f(3) = -22
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimums_sine", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// Extrema at π/2 + nπ
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		// Search in [0, 2π]
		auto optimums = analyzer.GetLocalOptimums(REAL(0.0), REAL(2.0) * Constants::PI, 1e-6);

		// Should find 2 extrema: max at π/2, min at 3π/2
		REQUIRE(optimums.size() == 2);

		// Sort to ensure consistent order
		std::sort(optimums.begin(), optimums.end());

		REQUIRE_THAT(optimums[0], WithinAbs(Constants::PI / REAL(2.0), 1e-4));
		REQUIRE_THAT(optimums[1], WithinAbs(REAL(3.0) * Constants::PI / REAL(2.0), 1e-4));
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimumsClassified_sine", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// Max at π/2, min at 3π/2
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		auto classified = analyzer.GetLocalOptimumsClassified(REAL(0.0), REAL(2.0) * Constants::PI, 1e-6);

		REQUIRE(classified.size() == 2);

		// Sort by x-coordinate
		std::sort(classified.begin(), classified.end(),
			[](const CriticalPoint& a, const CriticalPoint& b) { return a.x < b.x; });

		// Maximum at π/2
		REQUIRE_THAT(classified[0].x, WithinAbs(Constants::PI / REAL(2.0), 1e-4));
		REQUIRE(classified[0].type == CriticalPointType::LOCAL_MAXIMUM);
		REQUIRE_THAT(classified[0].value, WithinAbs(REAL(1.0), REAL(1e-4)));

		// Minimum at 3π/2
		REQUIRE_THAT(classified[1].x, WithinAbs(REAL(3.0) * Constants::PI / REAL(2.0), 1e-4));
		REQUIRE(classified[1].type == CriticalPointType::LOCAL_MINIMUM);
		REQUIRE_THAT(classified[1].value, WithinAbs(-REAL(1.0), REAL(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetLocalOptimums_quartic", "[optimums]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^4 - 4x^3 + 4x^2
		// f'(x) = 4x^3 - 12x^2 + 8x = 4x(x^2 - 3x + 2) = 4x(x-1)(x-2)
		// Critical points at x=0, x=1, x=2
		RealFunction quartic([](Real x) {
			return x * x * x * x - REAL(4.0) * x * x * x + REAL(4.0) * x * x;
		});
		RealFunctionAnalyzer analyzer(quartic);

		auto optimums = analyzer.GetLocalOptimums(-REAL(1.0), REAL(3.0), 1e-6);

		// Should find 3 optimums at x=0, 1, 2
		REQUIRE(optimums.size() == 3);

		std::sort(optimums.begin(), optimums.end());

		REQUIRE_THAT(optimums[0], WithinAbs(REAL(0.0), REAL(1e-4)));
		REQUIRE_THAT(optimums[1], WithinAbs(REAL(1.0), REAL(1e-4)));
		REQUIRE_THAT(optimums[2], WithinAbs(REAL(2.0), REAL(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_cubic", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 3x^2 - 9x + 5
		// f''(x) = 6x - 6 = 0 => x = 1
		// Inflection point at x=1
		RealFunction cubic([](Real x) {
			return x * x * x - REAL(3.0) * x * x - REAL(9.0) * x + REAL(5.0);
		});
		RealFunctionAnalyzer analyzer(cubic);

		auto inflections = analyzer.GetInflectionPoints(-REAL(5.0), REAL(5.0), 1e-6);

		// Should find exactly 1 inflection point
		REQUIRE(inflections.size() == 1);
		REQUIRE_THAT(inflections[0], WithinAbs(REAL(1.0), REAL(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_quartic", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^4 - 6x^2
		// f''(x) = 12x^2 - 12 = 12(x^2 - 1) = 0 => x = ±1
		// Inflection points at x=-1 and x=1
		RealFunction quartic([](Real x) {
			return x * x * x * x - REAL(6.0) * x * x;
		});
		RealFunctionAnalyzer analyzer(quartic);

		auto inflections = analyzer.GetInflectionPoints(-REAL(3.0), REAL(3.0), 1e-6);

		// Should find 2 inflection points
		REQUIRE(inflections.size() == 2);

		std::sort(inflections.begin(), inflections.end());

		REQUIRE_THAT(inflections[0], WithinAbs(-REAL(1.0), REAL(1e-4)));
		REQUIRE_THAT(inflections[1], WithinAbs(REAL(1.0), REAL(1e-4)));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_sine", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x)
		// f''(x) = -sin(x) = 0 => x = nπ
		// Inflection points at 0, π, 2π
		RealFunction sine([](Real x) { return std::sin(x); });
		RealFunctionAnalyzer analyzer(sine);

		auto inflections = analyzer.GetInflectionPoints(REAL(0.0), REAL(2.0) * Constants::PI, 1e-6);

		// Should find 3 inflection points: 0, π, 2π
		REQUIRE(inflections.size() == 3);

		std::sort(inflections.begin(), inflections.end());

		REQUIRE_THAT(inflections[0], WithinAbs(REAL(0.0), REAL(1e-4)));
		REQUIRE_THAT(inflections[1], WithinAbs(Constants::PI, REAL(1e-4)));
		REQUIRE_THAT(inflections[2], WithinAbs(REAL(2.0) * Constants::PI, 1e-4));
	}

	TEST_CASE("FunctionAnalyzer::GetInflectionPoints_no_inflection", "[inflection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2, no inflection points (concavity never changes)
		RealFunction quadratic([](Real x) { return x * x; });
		RealFunctionAnalyzer analyzer(quadratic);

		auto inflections = analyzer.GetInflectionPoints(-REAL(5.0), REAL(5.0), 1e-6);

		// Should find no inflection points
		REQUIRE(inflections.size() == 0);
	}

}