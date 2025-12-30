#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Optimization.h"
#endif

using namespace MML;
using namespace MML::Testing;
// Removed: using Catch::Approx; - now using precision-aware matchers
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::Optimization
{
	using namespace MML;
using namespace MML::Testing;

	/************************************************************************/
	/*****                    Test Functions                            *****/
	/************************************************************************/
	
	// Simple quadratic: f(x) = (x-2)^2 + 1, minimum at x=2, f(2)=1
	class QuadraticFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override { return (x - REAL(2.0)) * (x - REAL(2.0)) + REAL(1.0); }
	};
	
	// Derivative: f'(x) = 2(x-2)
	class QuadraticDerivFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override { return REAL(2.0) * (x - REAL(2.0)); }
	};
	
	// Cubic with local minimum: f(x) = x^3 - 6x^2 + 9x + 1
	// Local min at x=3, f(3)=1; local max at x=1, f(1)=5
	class CubicFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override { 
			return x*x*x - REAL(6.0)*x*x + REAL(9.0)*x + REAL(1.0); 
		}
	};
	
	// Derivative: f'(x) = 3x^2 - 12x + 9 = 3(x-1)(x-3)
	class CubicDerivFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override { 
			return REAL(3.0)*x*x - REAL(12.0)*x + REAL(9.0); 
		}
	};
	
	// Rosenbrock in 1D slice: f(x) = (1-x)^2 + 100*(0-x^2)^2 = (1-x)^2 + 100*x^4
	// Minimum at x ≈ REAL(0.78) (approx), complex shape
	class Rosenbrock1D : public IRealFunction {
	public:
		Real operator()(Real x) const override {
			Real a = REAL(1.0) - x;
			Real b = -x * x;  // y = 0, so this is (y - x^2)
			return a * a + REAL(100.0) * b * b;
		}
	};
	
	// sin(x) for testing local minima
	class SinFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override { return std::sin(x); }
	};
	
	class CosFunc : public IRealFunction {  // derivative of sin
	public:
		Real operator()(Real x) const override { return std::cos(x); }
	};
	
	// Quartic: f(x) = x^4 - 2x^2 + 1 = (x^2 - 1)^2, minima at x = ±1
	class QuarticFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override {
			Real x2 = x * x;
			return x2 * x2 - REAL(2.0) * x2 + REAL(1.0);
		}
	};
	
	class QuarticDerivFunc : public IRealFunction {
	public:
		Real operator()(Real x) const override {
			return REAL(4.0) * x * x * x - REAL(4.0) * x;
		}
	};

	/************************************************************************/
	/*****              REAL(10.1) Bracket Minimum Tests                      *****/
	/************************************************************************/
	
	TEST_CASE("BracketMinimum_QuadraticFunction_FindsValidBracket", "[Optimization][Bracket]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto bracket = Minimization::BracketMinimum(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(bracket.valid);
		// fb should be less than both fa and fc (minimum is bracketed)
		REQUIRE(bracket.fb <= bracket.fa);
		REQUIRE(bracket.fb <= bracket.fc);
		// The minimum at x=2 should be contained in the bracket
		Real minBound = std::min(bracket.ax, bracket.cx);
		Real maxBound = std::max(bracket.ax, bracket.cx);
		REQUIRE(minBound <= REAL(2.0));
		REQUIRE(maxBound >= REAL(2.0));
	}
	
	TEST_CASE("BracketMinimum_InitialPointsReversed_StillWorks", "[Optimization][Bracket]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		// Start with b < a (reversed order)
		auto bracket = Minimization::BracketMinimum(f, REAL(5.0), REAL(1.0));
		
		REQUIRE(bracket.valid);
		REQUIRE(bracket.fb <= bracket.fa);
		REQUIRE(bracket.fb <= bracket.fc);
	}
	
	TEST_CASE("BracketMinimum_MinimumAtOrigin_FindsBracket", "[Optimization][Bracket]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2, minimum at x=0
		RealFunctionFromStdFunc f([](Real x) { return x * x; });
		
		auto bracket = Minimization::BracketMinimum(f, -REAL(1.0), REAL(1.0));
		
		REQUIRE(bracket.valid);
		Real minBound = std::min(bracket.ax, bracket.cx);
		Real maxBound = std::max(bracket.ax, bracket.cx);
		REQUIRE(minBound <= REAL(0.0));
		REQUIRE(maxBound >= REAL(0.0));
	}
	
	TEST_CASE("BracketMinimum_CubicWithLocalMin_FindsBracket", "[Optimization][Bracket]")
	{
			TEST_PRECISION_INFO();
		CubicFunc f;
		
		// Search near the local minimum at x=3
		auto bracket = Minimization::BracketMinimum(f, REAL(2.0), REAL(4.0));
		
		REQUIRE(bracket.valid);
		REQUIRE(bracket.fb <= bracket.fa);
		REQUIRE(bracket.fb <= bracket.fc);
	}

	/************************************************************************/
	/*****          REAL(10.2) Golden Section Search Tests                    *****/
	/************************************************************************/
	
	TEST_CASE("GoldenSection_QuadraticFunction_FindsMinimum", "[Optimization][GoldenSection]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto result = Minimization::GoldenSectionSearch(f, REAL(0.0), REAL(5.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-6));
	}
	
	TEST_CASE("GoldenSection_WithBracket_FindsMinimum", "[Optimization][GoldenSection]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		auto bracket = Minimization::BracketMinimum(f, REAL(0.0), REAL(1.0));
		
		auto result = Minimization::GoldenSectionSearch(f, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-6));
	}
	
	TEST_CASE("GoldenSection_SinFunction_FindsLocalMin", "[Optimization][GoldenSection]")
	{
			TEST_PRECISION_INFO();
		SinFunc f;
		
		// sin(x) has minimum at x = 3π/2 ≈ REAL(4.712)
		// Create a known bracket around 3π/2: a=π, b=3π/2, c=2π
		MinimumBracket bracket;
		bracket.ax = Constants::PI;
		bracket.bx = REAL(3.0) * Constants::PI / REAL(2.0);
		bracket.cx = REAL(2.0) * Constants::PI;
		bracket.fa = f(bracket.ax);  // sin(π) = 0
		bracket.fb = f(bracket.bx);  // sin(3π/2) = -1
		bracket.fc = f(bracket.cx);  // sin(2π) = 0
		bracket.valid = true;
		
		auto result = Minimization::GoldenSectionSearch(f, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(3.0) * Constants::PI / REAL(2.0), 1e-5));
		REQUIRE_THAT(result.fmin, RealWithinRel(-REAL(1.0), 1e-6));
	}
	
	TEST_CASE("GoldenSection_QuarticFunction_FindsMinimumAtOne", "[Optimization][GoldenSection]")
	{
			TEST_PRECISION_INFO();
		QuarticFunc f;
		
		// Search for minimum at x=1 (there's also one at x=-1)
		auto result = Minimization::GoldenSectionSearch(f, REAL(0.5), REAL(2.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(1.0), 1e-5));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(0.0), 1e-6));
	}
	
	TEST_CASE("GoldenSection_QuarticFunction_FindsMinimumAtNegOne", "[Optimization][GoldenSection]")
	{
			TEST_PRECISION_INFO();
		QuarticFunc f;
		
		// Search for minimum at x=-1
		auto result = Minimization::GoldenSectionSearch(f, -REAL(2.0), -REAL(0.5));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(-REAL(1.0), 1e-5));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(0.0), 1e-6));
	}
	
	TEST_CASE("GoldenSection_CustomTolerance_RespectsIt", "[Optimization][GoldenSection]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		// Looser tolerance should require fewer iterations
		auto result1 = Minimization::GoldenSectionSearch(f, REAL(0.0), REAL(5.0), 1e-4);
		auto result2 = Minimization::GoldenSectionSearch(f, REAL(0.0), REAL(5.0), 1e-10);
		
		REQUIRE(result1.iterations < result2.iterations);
	}

	/************************************************************************/
	/*****     REAL(10.3) Brent's Method (Parabolic Interpolation) Tests      *****/
	/************************************************************************/
	
	TEST_CASE("BrentMinimize_QuadraticFunction_FindsMinimum", "[Optimization][Brent]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto result = Minimization::BrentMinimize(f, REAL(0.0), REAL(5.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-7));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentMinimize_WithBracket_FindsMinimum", "[Optimization][Brent]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		auto bracket = Minimization::BracketMinimum(f, -REAL(1.0), REAL(0.0));
		
		auto result = Minimization::BrentMinimize(f, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-7));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentMinimize_CubicLocalMin_FindsIt", "[Optimization][Brent]")
	{
			TEST_PRECISION_INFO();
		CubicFunc f;
		
		// Create a bracket around the local minimum at x=3
		// f(x) = x^3 - 6x^2 + 9x + 1, has min at x=3 with f(3)=1
		MinimumBracket bracket;
		bracket.ax = REAL(2.0);
		bracket.bx = REAL(3.0);
		bracket.cx = REAL(4.0);
		bracket.fa = f(bracket.ax);  // f(2) = 8 - 24 + 18 + 1 = 3
		bracket.fb = f(bracket.bx);  // f(3) = 27 - 54 + 27 + 1 = 1
		bracket.fc = f(bracket.cx);  // f(4) = 64 - 96 + 36 + 1 = 5
		bracket.valid = true;
		
		auto result = Minimization::BrentMinimize(f, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(3.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-6));
	}
	
	TEST_CASE("BrentMinimize_SinFunction_FindsMinimum", "[Optimization][Brent]")
	{
			TEST_PRECISION_INFO();
		SinFunc f;
		
		// Create a known bracket around 3π/2
		MinimumBracket bracket;
		bracket.ax = Constants::PI;
		bracket.bx = REAL(3.0) * Constants::PI / REAL(2.0);
		bracket.cx = REAL(2.0) * Constants::PI;
		bracket.fa = f(bracket.ax);
		bracket.fb = f(bracket.bx);
		bracket.fc = f(bracket.cx);
		bracket.valid = true;
		
		auto result = Minimization::BrentMinimize(f, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(3.0) * Constants::PI / REAL(2.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(-REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentMinimize_FasterThanGoldenSection", "[Optimization][Brent]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto goldenResult = Minimization::GoldenSectionSearch(f, REAL(0.0), REAL(5.0), 1e-8);
		auto brentResult = Minimization::BrentMinimize(f, REAL(0.0), REAL(5.0), 1e-8);
		
		// Brent should converge in fewer iterations for smooth quadratic
		REQUIRE(brentResult.iterations <= goldenResult.iterations);
		// Both should find the same minimum
		REQUIRE_THAT(brentResult.xmin, RealWithinRel(goldenResult.xmin, 1e-6));
	}
	
	TEST_CASE("BrentMinimize_Rosenbrock1D_Converges", "[Optimization][Brent]")
	{
			TEST_PRECISION_INFO();
		Rosenbrock1D f;
		
		auto result = Minimization::BrentMinimize(f, REAL(0.0), REAL(2.0));
		
		REQUIRE(result.converged);
		// The 1D Rosenbrock slice has minimum around x ≈ REAL(0.79) with f ≈ REAL(0.77)
		// The minimum is not 0 for this 1D slice
		REQUIRE(result.fmin < REAL(1.0));
	}

	/************************************************************************/
	/*****    REAL(10.4) Brent with Derivatives Tests                         *****/
	/************************************************************************/
	
	TEST_CASE("BrentWithDeriv_QuadraticFunction_FindsMinimum", "[Optimization][BrentDeriv]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		QuadraticDerivFunc df;
		
		auto result = Minimization::BrentMinimizeWithDeriv(f, df, REAL(0.0), REAL(5.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-7));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentWithDeriv_WithBracket_FindsMinimum", "[Optimization][BrentDeriv]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		QuadraticDerivFunc df;
		auto bracket = Minimization::BracketMinimum(f, -REAL(1.0), REAL(0.0));
		
		auto result = Minimization::BrentMinimizeWithDeriv(f, df, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-7));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentWithDeriv_CubicLocalMin_FindsIt", "[Optimization][BrentDeriv]")
	{
			TEST_PRECISION_INFO();
		CubicFunc f;
		CubicDerivFunc df;
		
		// Create a bracket around the local minimum at x=3
		MinimumBracket bracket;
		bracket.ax = REAL(2.0);
		bracket.bx = REAL(3.0);
		bracket.cx = REAL(4.0);
		bracket.fa = f(bracket.ax);
		bracket.fb = f(bracket.bx);
		bracket.fc = f(bracket.cx);
		bracket.valid = true;
		
		auto result = Minimization::BrentMinimizeWithDeriv(f, df, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(3.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-6));
	}
	
	TEST_CASE("BrentWithDeriv_SinFunction_FindsMinimum", "[Optimization][BrentDeriv]")
	{
			TEST_PRECISION_INFO();
		SinFunc f;
		CosFunc df;  // d/dx sin(x) = cos(x)
		
		// Create a known bracket around 3π/2
		MinimumBracket bracket;
		bracket.ax = Constants::PI;
		bracket.bx = REAL(3.0) * Constants::PI / REAL(2.0);
		bracket.cx = REAL(2.0) * Constants::PI;
		bracket.fa = f(bracket.ax);
		bracket.fb = f(bracket.bx);
		bracket.fc = f(bracket.cx);
		bracket.valid = true;
		
		auto result = Minimization::BrentMinimizeWithDeriv(f, df, bracket);
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(3.0) * Constants::PI / REAL(2.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(-REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentWithDeriv_QuarticFunction_FindsMinimum", "[Optimization][BrentDeriv]")
	{
			TEST_PRECISION_INFO();
		QuarticFunc f;
		QuarticDerivFunc df;
		
		auto result = Minimization::BrentMinimizeWithDeriv(f, df, REAL(0.5), REAL(2.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(1.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinAbs(REAL(0.0), 1e-12));
	}
	
	TEST_CASE("BrentWithDeriv_FasterThanBrentWithoutDeriv", "[Optimization][BrentDeriv]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		QuadraticDerivFunc df;
		
		auto brentResult = Minimization::BrentMinimize(f, REAL(0.0), REAL(5.0), 1e-10);
		auto dbrentResult = Minimization::BrentMinimizeWithDeriv(f, df, REAL(0.0), REAL(5.0), 1e-10);
		
		// With derivatives should converge in fewer or equal iterations
		REQUIRE(dbrentResult.iterations <= brentResult.iterations);
		// Both should find the same minimum
		REQUIRE_THAT(dbrentResult.xmin, RealWithinRel(brentResult.xmin, 1e-8));
	}

	/************************************************************************/
	/*****                 Maximization Tests                           *****/
	/************************************************************************/
	
	TEST_CASE("GoldenSectionMaximize_NegativeQuadratic_FindsMaximum", "[Optimization][Maximize]")
	{
			TEST_PRECISION_INFO();
		// f(x) = -(x-3)^2 + 5, maximum at x=3, f(3)=5
		RealFunctionFromStdFunc f([](Real x) { return -(x - REAL(3.0)) * (x - REAL(3.0)) + REAL(5.0); });
		
		auto result = Minimization::GoldenSectionMaximize(f, REAL(0.0), REAL(6.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(3.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(5.0), 1e-6));
	}
	
	TEST_CASE("BrentMaximize_SinFunction_FindsMaximum", "[Optimization][Maximize]")
	{
			TEST_PRECISION_INFO();
		SinFunc f;
		
		// sin(x) has maximum at x = π/2 ≈ REAL(1.571)
		// For maximize, we negate the function, so need bracket where sin is high
		// Use simple bounds where the maximum is clearly inside
		// The convenience function auto-brackets, which can expand.
		// Use the bracket-based version with negated function directly
		RealFunctionFromStdFunc negF([&f](Real x) { return -f(x); });
		
		MinimumBracket bracket;
		bracket.ax = REAL(0.0);
		bracket.bx = Constants::PI / REAL(2.0);
		bracket.cx = Constants::PI;
		bracket.fa = negF(bracket.ax);  // -sin(0) = 0
		bracket.fb = negF(bracket.bx);  // -sin(π/2) = -1
		bracket.fc = negF(bracket.cx);  // -sin(π) = 0
		bracket.valid = true;
		
		auto result = Minimization::BrentMinimize(negF, bracket);
		result.fmin = -result.fmin;  // Restore original sign
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(Constants::PI / REAL(2.0), 1e-6));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(1.0), 1e-7));
	}
	
	TEST_CASE("BrentMaximize_CubicLocalMax_FindsIt", "[Optimization][Maximize]")
	{
			TEST_PRECISION_INFO();
		CubicFunc f;
		
		// Cubic has local max at x=1, f(1)=5
		// f(x) = x^3 - 6x^2 + 9x + 1
		// f(0) = 1, f(1) = 5, f(2) = 3
		RealFunctionFromStdFunc negF([&f](Real x) { return -f(x); });
		
		MinimumBracket bracket;
		bracket.ax = REAL(0.0);
		bracket.bx = REAL(1.0);
		bracket.cx = REAL(2.0);
		bracket.fa = negF(bracket.ax);  // -1
		bracket.fb = negF(bracket.bx);  // -5 (minimum of negated = maximum of original)
		bracket.fc = negF(bracket.cx);  // -3
		bracket.valid = true;
		
		auto result = Minimization::BrentMinimize(negF, bracket);
		result.fmin = -result.fmin;
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(1.0), 1e-5));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(5.0), 1e-5));
	}

	/************************************************************************/
	/*****                Edge Cases and Error Handling                 *****/
	/************************************************************************/
	
	TEST_CASE("BracketMinimum_VeryNarrowInitialInterval_Works", "[Optimization][EdgeCase]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto bracket = Minimization::BracketMinimum(f, REAL(1.9), REAL(2.1));
		
		REQUIRE(bracket.valid);
	}
	
	TEST_CASE("BracketMinimum_WideInitialInterval_Works", "[Optimization][EdgeCase]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto bracket = Minimization::BracketMinimum(f, -REAL(100.0), REAL(100.0));
		
		REQUIRE(bracket.valid);
	}
	
	TEST_CASE("GoldenSection_VeryTightTolerance_Converges", "[Optimization][EdgeCase]")
	{
			TEST_PRECISION_INFO();
		QuadraticFunc f;
		
		auto result = Minimization::GoldenSectionSearch(f, REAL(0.0), REAL(5.0), 1e-12);
		
		REQUIRE(result.converged);
		// Golden section achieves good but not perfect precision
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(2.0), 1e-7));
	}
	
	TEST_CASE("BrentMinimize_FlatRegion_Converges", "[Optimization][EdgeCase]")
	{
			TEST_PRECISION_INFO();
		// f(x) = max(0, (x-2)^2 - 1) - flat at minimum
		RealFunctionFromStdFunc f([](Real x) { 
			Real val = (x - REAL(2.0)) * (x - REAL(2.0)) - REAL(1.0);
			return val > 0 ? val : REAL(0.0);
		});
		
		auto result = Minimization::BrentMinimize(f, REAL(0.0), REAL(5.0));
		
		REQUIRE(result.converged);
		// Should find something in the flat region [1, 3]
		REQUIRE(result.xmin >= REAL(1.0) - REAL(0.1));
		REQUIRE(result.xmin <= REAL(3.0) + REAL(0.1));
	}
	
	TEST_CASE("Optimization_LambdaFunction_Works", "[Optimization][Lambda]")
	{
			TEST_PRECISION_INFO();
		// Test with lambda wrapped in RealFunctionFromStdFunc
		RealFunctionFromStdFunc f([](Real x) { return (x - REAL(1.5)) * (x - REAL(1.5)) + REAL(0.5); });
		
		auto result = Minimization::BrentMinimize(f, -REAL(5.0), REAL(5.0));
		
		REQUIRE(result.converged);
		REQUIRE_THAT(result.xmin, RealWithinRel(REAL(1.5), 1e-7));
		REQUIRE_THAT(result.fmin, RealWithinRel(REAL(0.5), 1e-7));
	}

}  // namespace MML::Tests::Algorithms::Optimization

/***************************************************************************************************
 * COMPREHENSIVE TEST BED INTEGRATION TESTS
 * 
 * These tests iterate over ALL 1D test cases from optimization_test_bed.h to ensure
 * complete coverage of the 1D optimization algorithms.
 ***************************************************************************************************/

#include "../../test_data/optimization_test_bed.h"

using namespace MML::TestBeds;

TEST_CASE("GoldenSection_All1DUnimodalTestBed", "[Optimization][GoldenSection][TestBed][Comprehensive]")
{
	TEST_PRECISION_INFO();
	
	auto unimodalTests = getUnimodal1DTests();
	
	INFO("Testing " << unimodalTests.size() << " unimodal 1D functions from optimization test bed");
	
	for (const auto& test : unimodalTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Description: " << test.description);
			INFO("Search range: [" << test.searchLow << ", " << test.searchHigh << "]");
			INFO("True minimum at x=" << test.trueMinimumX << ", f=" << test.trueMinimumF);
			
			OptRealFunctionWrapper f(test.func);
			
			auto result = Minimization::GoldenSectionSearch(f, test.searchLow, test.searchHigh);
			
			CHECK(result.converged);
			
			// Adjust tolerance based on difficulty
			Real xTol = 1e-5;
			Real fTol = 1e-5;
			if (test.difficulty >= 2) { xTol = 1e-4; fTol = 1e-4; }
			if (test.difficulty >= 3) { xTol = 1e-3; fTol = 1e-3; }
			
			if (result.converged)
			{
				INFO("Found minimum at x=" << result.xmin << ", f=" << result.fmin);
				CHECK(std::abs(result.xmin - test.trueMinimumX) < xTol);
				CHECK(std::abs(result.fmin - test.trueMinimumF) < fTol);
			}
		}
	}
}

TEST_CASE("Brent_All1DUnimodalTestBed", "[Optimization][Brent][TestBed][Comprehensive]")
{
	TEST_PRECISION_INFO();
	
	auto unimodalTests = getUnimodal1DTests();
	
	INFO("Testing " << unimodalTests.size() << " unimodal 1D functions with Brent's method");
	
	for (const auto& test : unimodalTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Description: " << test.description);
			
			OptRealFunctionWrapper f(test.func);
			
			auto result = Minimization::BrentMinimize(f, test.searchLow, test.searchHigh);
			
			CHECK(result.converged);
			
			// Brent is more accurate than Golden Section
			Real xTol = 1e-6;
			Real fTol = 1e-6;
			if (test.difficulty >= 2) { xTol = 1e-5; fTol = 1e-5; }
			if (test.difficulty >= 3) { xTol = 1e-4; fTol = 1e-4; }
			
			if (result.converged)
			{
				INFO("Found minimum at x=" << result.xmin << ", f=" << result.fmin);
				CHECK(std::abs(result.xmin - test.trueMinimumX) < xTol);
				CHECK(std::abs(result.fmin - test.trueMinimumF) < fTol);
			}
		}
	}
}

TEST_CASE("BrentWithDeriv_All1DWithDerivativeTestBed", "[Optimization][BrentDeriv][TestBed][Comprehensive]")
{
	TEST_PRECISION_INFO();
	
	auto unimodalTests = getUnimodal1DTests();
	
	INFO("Testing functions with derivatives using Brent with derivative");
	
	for (const auto& test : unimodalTests)
	{
		if (!test.hasDerivative || test.derivative == nullptr)
			continue;
		
		DYNAMIC_SECTION(test.name + " (with derivative)")
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			
			OptRealFunctionWrapper f(test.func);
			OptRealFunctionWrapper df(test.derivative);
			
			auto result = Minimization::BrentMinimizeWithDeriv(f, df, test.searchLow, test.searchHigh);
			
			CHECK(result.converged);
			
			// With derivative info, should be even more accurate
			Real xTol = 1e-7;
			Real fTol = 1e-7;
			if (test.difficulty >= 2) { xTol = 1e-6; fTol = 1e-6; }
			if (test.difficulty >= 3) { xTol = 1e-5; fTol = 1e-5; }
			
			if (result.converged)
			{
				INFO("Found minimum at x=" << result.xmin << ", f=" << result.fmin);
				CHECK(std::abs(result.xmin - test.trueMinimumX) < xTol);
				CHECK(std::abs(result.fmin - test.trueMinimumF) < fTol);
			}
		}
	}
}

TEST_CASE("GoldenSection_Multimodal1DTestBed", "[Optimization][GoldenSection][TestBed][Multimodal]")
{
	TEST_PRECISION_INFO();
	
	auto multimodalTests = getMultimodal1DTests();
	
	INFO("Testing " << multimodalTests.size() << " multimodal 1D functions");
	INFO("Note: Golden Section finds LOCAL minima - we verify the found minimum is valid");
	
	for (const auto& test : multimodalTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Description: " << test.description);
			INFO("Global minimum at x=" << test.trueMinimumX << ", f=" << test.trueMinimumF);
			
			OptRealFunctionWrapper f(test.func);
			
			auto result = Minimization::GoldenSectionSearch(f, test.searchLow, test.searchHigh);
			
			CHECK(result.converged);
			
			if (result.converged)
			{
				INFO("Found minimum at x=" << result.xmin << ", f=" << result.fmin);
				
				// For multimodal, verify it's a valid local minimum:
				// f(xmin) should be <= f(xmin ± small delta)
				Real delta = 1e-4;
				Real fLeft = test.func(result.xmin - delta);
				Real fRight = test.func(result.xmin + delta);
				
				CHECK(result.fmin <= fLeft + 1e-8);
				CHECK(result.fmin <= fRight + 1e-8);
				
				// Check if we found the global minimum (not required but interesting)
				bool foundGlobal = std::abs(result.fmin - test.trueMinimumF) < 0.1;
				if (foundGlobal)
					INFO("Found GLOBAL minimum");
				else
					INFO("Found local minimum (expected for multimodal)");
			}
		}
	}
}

TEST_CASE("Brent_Multimodal1DTestBed", "[Optimization][Brent][TestBed][Multimodal]")
{
	TEST_PRECISION_INFO();
	
	auto multimodalTests = getMultimodal1DTests();
	
	INFO("Testing " << multimodalTests.size() << " multimodal 1D functions with Brent");
	
	for (const auto& test : multimodalTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Category: " << test.category);
			INFO("Difficulty: " << test.difficulty);
			INFO("Description: " << test.description);
			
			OptRealFunctionWrapper f(test.func);
			
			auto result = Minimization::BrentMinimize(f, test.searchLow, test.searchHigh);
			
			CHECK(result.converged);
			
			if (result.converged)
			{
				INFO("Found minimum at x=" << result.xmin << ", f=" << result.fmin);
				
				// Verify it's a valid local minimum
				Real delta = 1e-4;
				Real fLeft = test.func(result.xmin - delta);
				Real fRight = test.func(result.xmin + delta);
				
				CHECK(result.fmin <= fLeft + 1e-8);
				CHECK(result.fmin <= fRight + 1e-8);
			}
		}
	}
}

TEST_CASE("BracketMinimum_All1DTestBed", "[Optimization][Bracket][TestBed]")
{
	TEST_PRECISION_INFO();
	
	auto allTests = getAll1DOptimizationTests();
	
	INFO("Testing bracketing on " << allTests.size() << " 1D functions");
	
	for (const auto& test : allTests)
	{
		DYNAMIC_SECTION(test.name)
		{
			INFO("Search range: [" << test.searchLow << ", " << test.searchHigh << "]");
			
			OptRealFunctionWrapper f(test.func);
			
			// Start bracketing from the middle of the search range
			Real startA = (test.searchLow + test.searchHigh) / 2 - 0.5;
			Real startB = (test.searchLow + test.searchHigh) / 2 + 0.5;
			
			MinimumBracket bracket = Minimization::BracketMinimum(f, startA, startB);
			
			if (bracket.valid)
			{
				INFO("Found bracket: ax=" << bracket.ax << ", bx=" << bracket.bx << ", cx=" << bracket.cx);
				
				// Verify bracket property: f(ax) >= f(bx) && f(cx) >= f(bx)
				Real fa = f(bracket.ax);
				Real fb = f(bracket.bx);
				Real fc = f(bracket.cx);
				
				CHECK(fa >= fb - 1e-10);  // Small tolerance for numerical issues
				CHECK(fc >= fb - 1e-10);
				
				// Verify ordering: ax < bx < cx or ax > bx > cx
				bool ordered = (bracket.ax < bracket.bx && bracket.bx < bracket.cx) || 
				               (bracket.ax > bracket.bx && bracket.bx > bracket.cx);
				CHECK(ordered);
			}
			else
			{
				// Some functions may fail to bracket (e.g., monotonic in region)
				INFO("Bracketing failed for this test case");
			}
		}
	}
}

TEST_CASE("Optimization1D_ConvergenceComparison", "[Optimization][TestBed][Performance]")
{
	TEST_PRECISION_INFO();
	
	// Compare Golden Section vs Brent on standard test cases
	auto easyTests = getUnimodal1DTests();
	
	INFO("Comparing Golden Section vs Brent convergence");
	
	int goldenWins = 0, brentWins = 0;
	
	for (const auto& test : easyTests)
	{
		if (test.difficulty > 2) continue;  // Skip hard ones for this comparison
		
		DYNAMIC_SECTION(test.name + " comparison")
		{
			OptRealFunctionWrapper f(test.func);
			
			auto goldenResult = Minimization::GoldenSectionSearch(f, test.searchLow, test.searchHigh);
			auto brentResult = Minimization::BrentMinimize(f, test.searchLow, test.searchHigh);
			
			Real goldenError = std::abs(goldenResult.xmin - test.trueMinimumX);
			Real brentError = std::abs(brentResult.xmin - test.trueMinimumX);
			
			INFO("Golden Section error: " << goldenError);
			INFO("Brent error: " << brentError);
			
			if (brentError < goldenError)
				brentWins++;
			else
				goldenWins++;
			
			// Both should converge
			CHECK(goldenResult.converged);
			CHECK(brentResult.converged);
		}
	}
	
	// Brent should generally be better or equal
	INFO("Brent wins: " << brentWins << ", Golden Section wins: " << goldenWins);
}


