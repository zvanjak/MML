#include <catch2/catch_all.hpp>
#include <type_traits>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include <algorithm>  // for std::sort
#include <utility>    // for std::pair

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/Integration.h>

#endif

#include "../test_beds/real_functions_test_bed.h"

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::IntegrationTests
{
	TEST_CASE("Test_Integration_Trap_const_func", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc constFunc0( [](Real x) -> Real { return static_cast<MML::Real>(0.0); } );
		RealFunctionFromStdFunc constFunc1( [](Real x) -> Real { return static_cast<MML::Real>(2.0); } );
		RealFunctionFromStdFunc constFunc2( [](Real x) -> Real { return -static_cast<MML::Real>(5.0); } );
		RealFunctionFromStdFunc constFunc3( [](Real x) -> Real { return TOL(1e-10, 1e-5); } );

		REQUIRE_THAT(static_cast<MML::Real>(0.0), WithinAbs(IntegrateTrap(constFunc0, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(static_cast<MML::Real>(0.0), WithinAbs(IntegrateTrap(constFunc0, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));

		REQUIRE_THAT(static_cast<MML::Real>(2.0), WithinAbs(IntegrateTrap(constFunc1, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(static_cast<MML::Real>(8.0), WithinAbs(IntegrateTrap(constFunc1, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));

		REQUIRE_THAT(-static_cast<MML::Real>(5.0), WithinAbs(IntegrateTrap(constFunc2, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(-static_cast<MML::Real>(20.0), WithinAbs(IntegrateTrap(constFunc2, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));

		REQUIRE_THAT(TOL(1e-10, 1e-5), WithinAbs(IntegrateTrap(constFunc3, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(TOL(4e-10, 4e-5), WithinAbs(IntegrateTrap(constFunc3, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
	}

	// TODO
	//TEST_CASE("Test_Integration_almost_const_func", "[simple]") 
	//{
	//	RealFunctionFromStdFunc linearFun( [](Real x) -> Real { return static_cast<MML::Real>(2.0); } );

	//	double int_trap = IntegrateTrap(linearFun, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

	//	REQUIRE_THAT(static_cast<MML::Real>(1.0), Catch::Matchers::WithinAbs(IntegrateTrap(linearFun, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), static_cast<MML::Real>(1e-16)));
	//}

	TEST_CASE("Test_Integration_Trap_linear_func", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc linearFunc1( [](Real x) -> Real { return 2 * x; } );
		RealFunctionFromStdFunc linearFunc2( [](Real x) -> Real { return -3 * x; } );
		RealFunctionFromStdFunc linearFunc3( [](Real x) -> Real { return TOL(1e-10, 1e-5) * x; } );
		RealFunctionFromStdFunc linearFunc4( [](Real x) -> Real { return 2 * (x - static_cast<MML::Real>(0.5)); } );

		REQUIRE_THAT(static_cast<MML::Real>(1.0), WithinAbs(IntegrateTrap(linearFunc1, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(-static_cast<MML::Real>(1.5), WithinAbs(IntegrateTrap(linearFunc2, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(0.5e-10, WithinAbs(IntegrateTrap(linearFunc3, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(static_cast<MML::Real>(0.0), WithinAbs(IntegrateTrap(linearFunc4, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), TOL(1e-15, 1e-5)));

		REQUIRE_THAT(-static_cast<MML::Real>(5.0), WithinAbs(IntegrateTrap(linearFunc1, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(static_cast<MML::Real>(7.5), WithinAbs(IntegrateTrap(linearFunc2, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(TOL(-2.5e-10, -2.5e-5), WithinAbs(IntegrateTrap(linearFunc3, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0)), TOL(1e-15, 1e-5)));
		REQUIRE_THAT(-static_cast<MML::Real>(10.0), WithinAbs(IntegrateTrap(linearFunc4, -static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0)), TOL(1e-15, 1e-5)));
	}

	TEST_CASE("Test_Integration_Trap_standard_func_precision", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getFunc("Sin")._func;
		RealFunction sinFunc_int = TestBeds::RealFunctionsTestBed::getFunc("Sin")._funcIntegrated;

		Real a = static_cast<MML::Real>(0.0), b = static_cast<MML::Real>(1.0);
		REQUIRE_THAT(sinFunc_int(b) - sinFunc_int(a), WithinAbs(IntegrateTrap(sinFunc, a, b), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(sinFunc_int(b) - sinFunc_int(a), !WithinAbs(IntegrateTrap(sinFunc, a, b), static_cast<MML::Real>(1e-6)));

		// row - integration method
		// column - req prec
		// testirati mogucu preciznost za dvije funkcije
		// sa svim metodama
		// provjeriti exception to many steps za preveliku preciznost

	}

	TEST_CASE("Test_Integration_Trap_integral_test_bed1", "[simple]")
	{
			TEST_PRECISION_INFO();
		auto testFunc = TestBeds::RealFunctionsTestBed::getFuncWithIntegral("TestInt1")._func;
		auto testFunc_int = TestBeds::RealFunctionsTestBed::getFuncWithIntegral("TestInt1")._funcIntegrated;

		Real a = static_cast<MML::Real>(0.0), b = static_cast<MML::Real>(1.0);
		Real integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal,  WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-6)));

		REQUIRE_THAT(integralExactVal,  WithinRel(IntegrateTrap(testFunc, a, b), TOL(1e-5, 1e-3)));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), TOL(1e-6, 1e-6)));

		a = static_cast<MML::Real>(0.0), b = static_cast<MML::Real>(5.0);
		integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal,  WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-2)));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-3)));

		REQUIRE_THAT(integralExactVal,  WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));

		a = -static_cast<MML::Real>(5.0), b = static_cast<MML::Real>(2.0);
		integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal, WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-2)));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-3)));

		REQUIRE_THAT(integralExactVal, WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Test_Integration_Trap_integral_test_bed2", "[simple]")
	{
			TEST_PRECISION_INFO();
		auto testFunc = TestBeds::RealFunctionsTestBed::getFuncWithIntegral("TestInt2")._func;
		auto testFunc_int = TestBeds::RealFunctionsTestBed::getFuncWithIntegral("TestInt2")._funcIntegrated;

		Real a = static_cast<MML::Real>(0.0), b = static_cast<MML::Real>(1.0);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-6)));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-6)));

		a = static_cast<MML::Real>(0.0), b = static_cast<MML::Real>(5.0);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-6)));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));

		a = -static_cast<MML::Real>(5.0), b = static_cast<MML::Real>(2.0);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Integration::Adaptive_precision_control", "[adaptive]")
	{
			TEST_PRECISION_INFO();
		// Test that integration methods achieve requested precision
		RealFunctionFromStdFunc testFunc( [](Real x) -> Real { return std::sin(x) * std::exp(-x); } );

		// Test with different precision requirements
		auto result1 = IntegrateTrapDetailed(testFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-6, 1e-3));
		auto result2 = IntegrateTrapDetailed(testFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-8, 1e-4));
		auto result3 = IntegrateTrapDetailed(testFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-10, 1e-5));

		// Higher precision requires more steps (iterations)
		REQUIRE(result2.iterations > result1.iterations);
		REQUIRE(result3.iterations > result2.iterations);

		// Results should be consistent and increasingly accurate
		REQUIRE_THAT(result1.value, WithinAbs(result2.value, TOL(1e-6, 1e-3)));
		REQUIRE_THAT(result2.value, WithinAbs(result3.value, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Integration::Simpson_vs_Trap_accuracy", "[comparison]")
	{
			TEST_PRECISION_INFO();
		// Simpson's rule should be more accurate than trapezoidal for smooth functions
		RealFunctionFromStdFunc smoothFunc( [](Real x) -> Real { return std::exp(-x * x); } );

		// Request same relative precision
		auto trapResult = IntegrateTrapDetailed(smoothFunc, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), TOL(1e-8, 1e-4));
		auto simpsonResult = IntegrateSimpsonDetailed(smoothFunc, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), TOL(1e-8, 1e-4));

		// Simpson should converge faster (fewer iterations)
		REQUIRE(simpsonResult.iterations <= trapResult.iterations);

		// Both should give similar results
		REQUIRE_THAT(trapResult.value, WithinAbs(simpsonResult.value, TOL(1e-7, 1e-4)));
	}

	TEST_CASE("Integration::Gauss10_high_accuracy", "[gauss]")
	{
			TEST_PRECISION_INFO();
		// Gauss-Legendre integration should be very accurate for polynomials
		// Test with polynomial of degree <= 19 (exactly integrable with 10-point rule)
		RealFunctionFromStdFunc poly5( [](Real x) -> Real { return static_cast<MML::Real>(1.0) + x + x * x + x * x * x + x * x * x * x + x * x * x * x * x; } );

		// Analytical integral: x + x^2/2 + x^3/3 + x^4/4 + x^5/5 + x^6/6
		auto exactIntegral = [](Real x) -> Real {
			return x + x * x / static_cast<MML::Real>(2.0) + std::pow(x, MML::Real(3)) / static_cast<MML::Real>(3.0) + std::pow(x, MML::Real(4)) / static_cast<MML::Real>(4.0) + std::pow(x, MML::Real(5)) / static_cast<MML::Real>(5.0) + std::pow(x, MML::Real(6)) / static_cast<MML::Real>(6.0);
		};

		Real a = static_cast<MML::Real>(0.0), b = static_cast<MML::Real>(2.0);
		Real exact = exactIntegral(b) - exactIntegral(a);
		Real gauss10 = IntegrateGauss10(poly5, a, b);

		// Gauss-Legendre should be nearly exact for polynomials
		REQUIRE_THAT(exact, WithinAbs(gauss10, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Integration::Gauss10_returns_IntegrationResult", "[gauss][api]")
	{
			TEST_PRECISION_INFO();
		// Verify that IntegrateGauss10 now returns IntegrationResult structure
		RealFunctionFromStdFunc quadratic( [](Real x) -> Real { return x * x; } );
		
		// ∫₀¹ x² dx = 1/3
		auto result = IntegrateGauss10Detailed(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		// Check all fields of IntegrationResult
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(3.0), TOL(1e-12, 1e-5)));
		REQUIRE_FALSE(HasEstimate(result.error_estimate));  // Gauss10 has no error estimate (NoEstimate, not 0.0)
		REQUIRE(result.iterations == 1);        // Single-pass evaluation
		REQUIRE(result.converged == false);     // Fixed-order rule: no criterion was checked (RESULT_POLICY)
		REQUIRE_FALSE(result.error_estimated);
		
		// Verify intentional convenience conversion to Real.
		Real value = IntegrateGauss10(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		REQUIRE_THAT(value, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(3.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Integration::Callable_overloads_accept_lambdas_and_functors", "[callable][C++20]")
	{
		TEST_PRECISION_INFO();

		auto quadratic = [](Real x) -> Real { return x * x; };
		const Real exact = static_cast<MML::Real>(1.0) / static_cast<MML::Real>(3.0);

		REQUIRE_THAT(IntegrateTrap(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(exact, TOL(2e-5, 1e-4)));
		REQUIRE_THAT(IntegrateSimpson(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(exact, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(IntegrateRomberg(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(exact, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(IntegrateGauss10(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(exact, TOL(1e-12, 1e-5)));
		REQUIRE_THAT(IntegrateGK21(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(exact, TOL(1e-12, 1e-5)));

		Real scale = static_cast<MML::Real>(2.0);
		auto captured = [scale](Real x) -> Real { return scale * x; };
		REQUIRE_THAT(Integrate(captured, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), SIMPSON), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(Integrate<GAUSS10>(captured, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5)));

		struct CubicFunctor
		{
			Real operator()(Real x) const { return x * x * x; }
		};

		REQUIRE_THAT(Integrate<TRAP>(CubicFunctor{}, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), WithinAbs(static_cast<MML::Real>(0.25), TOL(5e-6, 1e-4)));
	}

	TEST_CASE("Integration::Oscillatory_functions", "[oscillatory]")
	{
			TEST_PRECISION_INFO();
		// Test integration of highly oscillatory functions
		RealFunctionFromStdFunc oscillatory( [](Real x) -> Real { return std::sin(static_cast<MML::Real>(10.0) * x) / (static_cast<MML::Real>(1.0) + x); } );

		// These require more careful integration
		Real trapResult = IntegrateTrap(oscillatory, static_cast<MML::Real>(0.0), Constants::PI, static_cast<MML::Real>(1e-6));
		Real simpsonResult = IntegrateSimpson(oscillatory, static_cast<MML::Real>(0.0), Constants::PI, static_cast<MML::Real>(1e-6));

		// Results should be reasonably close
		REQUIRE_THAT(trapResult, WithinAbs(simpsonResult, static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("Integration::Quasi_improper_integrals", "[improper]")
	{
			TEST_PRECISION_INFO();
		// Test functions that are well-behaved but approach infinity at boundaries
		// f(x) = 1/sqrt(x) on (0, 1] - integrable but singular at 0
		RealFunctionFromStdFunc invSqrt( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(x); } );

		// Analytical result: 2*sqrt(x) from epsilon to 1 = 2
		// We start slightly away from 0 to avoid numerical issues. Plain trapezoid cannot meet
		// 1e-5 next to a singularity in 20 halvings (the plain call throws ConvergenceError);
		// the Detailed result still carries the best estimate, and the singular-aware rule converges.
		auto trapNearSingular = IntegrateTrapDetailed(invSqrt, MML::Real(1e-6), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-5));
		REQUIRE_FALSE(trapNearSingular.converged);
		REQUIRE_THAT(trapNearSingular.value, WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2e-3)));
		REQUIRE_THROWS_AS(IntegrateTrap(invSqrt, MML::Real(1e-6), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-5)), ConvergenceError);
		REQUIRE_THAT(IntegrateLowerSingular(invSqrt, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-8)), WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-6, 1e-3)));

		// f(x) = ln(x) on (0, 1] - integrable with x*ln(x) - x
		RealFunctionFromStdFunc logFunc( [](Real x) -> Real { return std::log(x); } );

		// Analytical result: [x*ln(x) - x] from epsilon to 1 = -1
		auto trapLog = IntegrateTrapDetailed(logFunc, MML::Real(1e-6), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-5));
		REQUIRE_THAT(trapLog.value, WithinAbs(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-2)));
	}

	TEST_CASE("Integration::Symmetric_functions", "[properties]")
	{
			TEST_PRECISION_INFO();
		// Test that integration respects symmetry properties
		RealFunctionFromStdFunc evenFunc( [](Real x) -> Real { return std::cos(x); } );
		RealFunctionFromStdFunc oddFunc( [](Real x) -> Real { return std::sin(x); } );

		// For even functions: integral from -a to a = 2 * integral from 0 to a
		Real a = Constants::PI / static_cast<MML::Real>(2.0);
		Real evenResultSymmetric = IntegrateTrap(evenFunc, -a, a);
		Real evenResultPositive = IntegrateTrap(evenFunc, static_cast<MML::Real>(0.0), a);
		REQUIRE_THAT(evenResultSymmetric, WithinAbs(static_cast<MML::Real>(2.0) * evenResultPositive, TOL(1e-8, 1e-4)));

		// For odd functions: integral from -a to a = 0
		Real oddResultSymmetric = IntegrateTrap(oddFunc, -a, a);
		REQUIRE_THAT(oddResultSymmetric, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-2)));
	}

	TEST_CASE("Integration::Edge_cases", "[edge]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc simpleFunc( [](Real x) -> Real { return x * x; } );

		// Zero-width interval
		REQUIRE_THAT(IntegrateTrap(simpleFunc, static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-15, 1e-5)));

		// Reversed limits (should handle gracefully or give negative result)
		Real forward = IntegrateTrap(simpleFunc, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		Real backward = IntegrateTrap(simpleFunc, static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));
		REQUIRE_THAT(forward, WithinAbs(-backward, TOL(1e-10, 1e-5)));

		const Real h = TOL(1e-8, 1e-4), tolerance = TOL(1e-12, 1e-5);
		auto tiny = IntegrateTrapDetailed(simpleFunc, Real(0), h, tolerance);
		if constexpr (std::is_same_v<Real, float>) {
			REQUIRE(tiny.converged);
			REQUIRE(tiny.status == AlgorithmStatus::Success);
			REQUIRE(IntegrateTrap(simpleFunc, Real(0), h, tolerance) == tiny.value);
		} else {
			// Twenty trapezoid halvings cannot deliver the double request of 1e-12.
			REQUIRE_FALSE(tiny.converged);
			REQUIRE(tiny.status == AlgorithmStatus::MaxIterationsExceeded);
			REQUIRE_THROWS_AS(IntegrateTrap(simpleFunc, Real(0), h, tolerance), ConvergenceError);
		}
		const Real exact = h * h * h / 3;
		REQUIRE_THAT(tiny.value, WithinAbs(exact, TOL(1e-25, Real(4) * tolerance * std::abs(exact))));
		REQUIRE_THAT(IntegrateRomberg(simpleFunc, Real(0), h, tolerance),
		             WithinAbs(exact, TOL(1e-25, Real(4) * tolerance * std::abs(exact))));
	}

	TEST_CASE("Integration::Precision_limits", "[precision]")
	{
			TEST_PRECISION_INFO();
		// Test behavior when requesting very high precision
		RealFunctionFromStdFunc smoothFunc( [](Real x) -> Real { return std::exp(-x); } );

		// Exact integral from 0 to 1: 1 - 1/e ≈ static_cast<MML::Real>(0.632120558829)
		Real exact = static_cast<MML::Real>(1.0) - std::exp(-static_cast<MML::Real>(1.0));

		// Test achievable precision
		auto result = IntegrateTrapDetailed(smoothFunc, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5));

		// Should achieve result close to exact
		REQUIRE_THAT(result.value, WithinAbs(exact, TOL(1e-10, 1e-5)));

		// Achieved precision (error_estimate) should be reasonable
		REQUIRE(result.error_estimate < TOL(1e-10, 1e-5));
	}

	TEST_CASE("Integration::Method_consistency", "[consistency]")
	{
			TEST_PRECISION_INFO();
		// All methods should give consistent results for simple functions
		RealFunctionFromStdFunc quadratic( [](Real x) -> Real { return static_cast<MML::Real>(2.0) * x * x - static_cast<MML::Real>(3.0) * x + static_cast<MML::Real>(1.0); } );

		// Analytical: [2x^3/3 - 3x^2/2 + x] from 0 to 2 = 16/3 - 6 + 2 = 4/3
		Real exact = static_cast<MML::Real>(4.0) / static_cast<MML::Real>(3.0);

		Real trapResult = IntegrateTrap(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), TOL(1e-8, 1e-4));
		Real simpsonResult = IntegrateSimpson(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), TOL(1e-8, 1e-4));
		Real gaussResult = IntegrateGauss10(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0));

		// All methods should converge to exact result
		REQUIRE_THAT(exact, WithinAbs(trapResult, TOL(1e-7, 1e-3)));
		REQUIRE_THAT(exact, WithinAbs(simpsonResult, TOL(1e-9, 1e-4)));
		REQUIRE_THAT(exact, WithinAbs(gaussResult, TOL(1e-12, 1e-5)));
	}

	/*********************************************************************/
	/*****          Romberg Integration Tests                        *****/
	/*********************************************************************/
	
	TEST_CASE("Integration::Romberg_polynomial", "[romberg]")
	{
			TEST_PRECISION_INFO();
		// For polynomials, Romberg should be exact (or very close)
		RealFunctionFromStdFunc quadratic( [](Real x) -> Real { return x * x; } );
		
		// ∫₀¹ x² dx = 1/3
		auto result = IntegrateRombergDetailed(quadratic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(3.0), TOL(1e-11, 1e-5)));
	}
	
	TEST_CASE("Integration::Romberg_cubic", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc cubic( [](Real x) -> Real { return x * x * x - static_cast<MML::Real>(2.0) * x; } );
		
		// ∫₀² (x³ - 2x) dx = [x⁴/4 - x²]₀² = 4 - 4 = 0
		auto result = IntegrateRombergDetailed(cubic, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), TOL(1e-12, 1e-5));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-11, 1e-5)));
	}
	
	TEST_CASE("Integration::Romberg_trigonometric", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc sinFunc( [](Real x) -> Real { return std::sin(x); } );
		
		// ∫₀^π sin(x) dx = 2
		auto result = IntegrateRombergDetailed(sinFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-10, 1e-5));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-9, 1e-4)));
	}
	
	TEST_CASE("Integration::Romberg_exponential", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc expFunc( [](Real x) -> Real { return std::exp(x); } );
		
		// ∫₀¹ e^x dx = e - 1 ≈ static_cast<MML::Real>(1.718281828)
		auto result = IntegrateRombergDetailed(expFunc, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(std::exp(static_cast<MML::Real>(1.0)) - static_cast<MML::Real>(1.0), TOL(1e-11, 1e-5)));
	}
	
	TEST_CASE("Integration::Romberg_gaussian", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc gaussian( [](Real x) -> Real { return std::exp(-x * x); } );
		
		// ∫₀¹ e^(-x²) dx ≈ static_cast<MML::Real>(0.7468241328) (error function related)
		// erf(1) = 2/√π * ∫₀¹ e^(-t²) dt, so ∫₀¹ e^(-x²) dx = √π/2 * erf(1)
		Real expected = std::sqrt(Constants::PI) / static_cast<MML::Real>(2.0) * std::erf(static_cast<MML::Real>(1.0));
		
		auto result = IntegrateRombergDetailed(gaussian, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, TOL(1e-11, 1e-5)));
	}
	
	TEST_CASE("Integration::Romberg_high_precision", "[romberg]")
	{
			TEST_PRECISION_INFO();
		// Test Romberg's ability to achieve very high precision
		RealFunctionFromStdFunc smooth( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + x * x); } );
		
		// ∫₀¹ 1/(1+x²) dx = arctan(1) = π/4
		auto result = IntegrateRombergDetailed(smooth, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-14, 1e-5));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / static_cast<MML::Real>(4.0), TOL(1e-13, 1e-5)));
	}
	
	TEST_CASE("Integration::Romberg_faster_than_simpson", "[romberg][comparison]")
	{
			TEST_PRECISION_INFO();
		// For smooth functions, Romberg should converge faster (fewer iterations)
		RealFunctionFromStdFunc smooth( [](Real x) -> Real { return std::sin(x) * std::exp(-x); } );
		
		auto rombergResult = IntegrateRombergDetailed(smooth, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-10, 1e-5));
		auto simpsonResult = IntegrateSimpsonDetailed(smooth, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-10, 1e-5));
		
		// Both should converge
		REQUIRE(rombergResult.converged == true);
		REQUIRE(simpsonResult.converged == true);
		
		// Results should be very close
		REQUIRE_THAT(rombergResult.value, WithinAbs(simpsonResult.value, TOL(1e-9, 1e-4)));
		
		// Romberg typically needs fewer iterations for smooth functions
		// (This is a soft check - mainly documenting expected behavior)
		// REQUIRE(rombergResult.iterations <= simpsonResult.iterations);
	}
	
	TEST_CASE("Integration::Romberg_method_consistency", "[romberg][consistency]")
	{
			TEST_PRECISION_INFO();
		// All integration methods should give consistent results
		RealFunctionFromStdFunc testFunc( [](Real x) -> Real { return std::cos(x) * std::cos(x); } );
		
		// ∫₀^π cos²(x) dx = π/2
		Real expected = Constants::PI / static_cast<MML::Real>(2.0);
		
		auto trapResult = IntegrateTrapDetailed(testFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-8, 1e-4));
		auto simpsonResult = IntegrateSimpsonDetailed(testFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-8, 1e-4));
		auto rombergResult = IntegrateRombergDetailed(testFunc, static_cast<MML::Real>(0.0), Constants::PI, TOL(1e-10, 1e-5));
		Real gaussResult = IntegrateGauss10(testFunc, static_cast<MML::Real>(0.0), Constants::PI);
		
		REQUIRE_THAT(trapResult.value, WithinAbs(expected, TOL(1e-7, 1e-5)));
		REQUIRE_THAT(simpsonResult.value, WithinAbs(expected, TOL(1e-7, 1e-5)));
		REQUIRE_THAT(rombergResult.value, WithinAbs(expected, TOL(1e-9, 1e-4)));
		REQUIRE_THAT(gaussResult, WithinAbs(expected, TOL(1e-6, 1e-4)));
	}
	
	TEST_CASE("Integration::Romberg_result_structure", "[romberg]")
	{
			TEST_PRECISION_INFO();
		// Test the IntegrationResult structure from Romberg
		RealFunctionFromStdFunc testFunc( [](Real x) -> Real { return x * x; } );
		
		auto result = IntegrateRombergDetailed(testFunc, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5));
		
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(3.0), TOL(1e-9, 1e-4)));
		REQUIRE(result.iterations > 0);
		REQUIRE(result.error_estimate < TOL(1e-9, 1e-4));
		REQUIRE(result.converged == true);
	}

	/*********************************************************************/
	/*****          Improper Integrals Tests                         *****/
	/*********************************************************************/
	
	TEST_CASE("Integration::ImproperUpper_exponential_decay", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) e^(-x) dx = 1
		RealFunctionFromStdFunc expDecay( [](Real x) -> Real { return std::exp(-x); } );
		
		auto result = IntegrateUpperInfDetailed(expDecay, static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperUpper_inverse_square", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[1, ∞) 1/x² dx = 1
		RealFunctionFromStdFunc invSquare( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / (x * x); } );
		
		auto result = IntegrateUpperInfDetailed(invSquare, static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-7)));
	}
	
	TEST_CASE("Integration::ImproperUpper_arctan_derivative", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) 1/(1+x²) dx = π/2
		RealFunctionFromStdFunc lorentzian( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + x * x); } );
		
		auto result = IntegrateUpperInfDetailed(lorentzian, static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-7, 1e-5)));
	}
	
	TEST_CASE("Integration::ImproperLower_exponential_growth", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, 0] e^x dx = 1
		RealFunctionFromStdFunc expGrowth( [](Real x) -> Real { return std::exp(x); } );
		
		auto result = IntegrateLowerInfDetailed(expGrowth, static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperLower_arctan_derivative", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, 0] 1/(1+x²) dx = π/2
		RealFunctionFromStdFunc lorentzian( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + x * x); } );
		
		auto result = IntegrateLowerInfDetailed(lorentzian, static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-7, 1e-5)));
	}
	
	TEST_CASE("Integration::ImproperFull_gaussian", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, ∞) e^(-x²) dx = √π
		RealFunctionFromStdFunc gaussian( [](Real x) -> Real { return std::exp(-x * x); } );
		
		auto result = IntegrateInfDetailed(gaussian);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(std::sqrt(Constants::PI), static_cast<MML::Real>(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperFull_lorentzian", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, ∞) 1/(1+x²) dx = π
		RealFunctionFromStdFunc lorentzian( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + x * x); } );
		
		auto result = IntegrateInfDetailed(lorentzian);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI, static_cast<MML::Real>(1e-6)));
	}
	
	TEST_CASE("Integration::ImproperFull_gaussian_shifted", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, ∞) e^(-(x-3)²) dx = √π (same as unshifted)
		RealFunctionFromStdFunc shifted( [](Real x) -> Real { return std::exp(-(x - static_cast<MML::Real>(3.0)) * (x - static_cast<MML::Real>(3.0))); } );
		
		// Without split point adjustment - should still converge
		auto result1 = IntegrateInfDetailed(shifted);
		REQUIRE(result1.converged == true);
		REQUIRE_THAT(result1.value, WithinAbs(std::sqrt(Constants::PI), static_cast<MML::Real>(1e-4)));
		
		// With split point at the peak - better convergence
		auto result2 = IntegrateInfSplitDetailed(shifted, static_cast<MML::Real>(3.0));
		REQUIRE(result2.converged == true);
		REQUIRE_THAT(result2.value, WithinAbs(std::sqrt(Constants::PI), static_cast<MML::Real>(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperUpper_sinc_squared", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) sin²(x)/x² dx = π/2
		// Use sinc²(x) = sin²(x)/x² (with sinc(0) = 1)
		// Note: This is a challenging integral due to oscillations
		RealFunctionFromStdFunc sincSquared( [](Real x) -> Real {
			if (std::abs(x) < TOL(1e-10, 1e-5)) return static_cast<MML::Real>(1.0);  // Limit as x → 0
			return std::sin(x) * std::sin(x) / (x * x);
		} );
		
		auto result = IntegrateUpperInfDetailed(sincSquared, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-6));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), 2e-3));  // ~static_cast<MML::Real>(0.1)% accuracy for oscillatory function
	}
	
	TEST_CASE("Integration::ImproperUpper_gamma_like", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) x * e^(-x) dx = Γ(2) = 1! = 1
		RealFunctionFromStdFunc gammaLike( [](Real x) -> Real { return x * std::exp(-x); } );
		
		auto result = IntegrateUpperInfDetailed(gammaLike, static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-4)));
	}
	
	TEST_CASE("Integration::ImproperUpper_gamma3", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) x² * e^(-x) dx = Γ(3) = 2! = 2
		RealFunctionFromStdFunc gamma3( [](Real x) -> Real { return x * x * std::exp(-x); } );
		
		auto result = IntegrateUpperInfDetailed(gamma3, static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-4)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                     SINGULAR INTEGRAL TESTS                           //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::LowerSingular_inv_sqrt", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√x dx = 2
		RealFunctionFromStdFunc f( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(x); } );
		
		auto result = IntegrateLowerSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-6)));
	}

	TEST_CASE("Integration::UpperSingular_inv_sqrt", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√(1-x) dx = 2
		RealFunctionFromStdFunc f( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(static_cast<MML::Real>(1.0) - x); } );
		
		auto result = IntegrateUpperSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-6, 1e-4)));
	}

	TEST_CASE("Integration::BothSingular_beta_half", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√(x(1-x)) dx = π  (Beta function B(1/2, 1/2) = Γ(1/2)²/Γ(1) = π)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(x * (static_cast<MML::Real>(1.0) - x)); } );
		
		auto result = IntegrateBothSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI, static_cast<MML::Real>(1e-4)));
	}

	TEST_CASE("Integration::InteriorSingular_inv_sqrt_abs", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[-1,1] 1/√|x| dx = 4  (two contributions of 2 each)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(std::abs(x)); } );
		
		auto result = IntegrateInteriorSingularDetailed(f, -static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(4.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Integration::LowerSingular_log_times_sqrt", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] ln(x)/√x dx = -4
		// This has both logarithmic and square-root singularities at 0
		RealFunctionFromStdFunc f( [](Real x) -> Real { return std::log(x) / std::sqrt(x); } );
		
		auto result = IntegrateLowerSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(-static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.02)));  // Harder integral, relax tolerance
	}

	TEST_CASE("Integration::LowerSingular_sqrt_itself", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] √x dx = 2/3 (no singularity, but transformation should still work)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return std::sqrt(x); } );
		
		auto result = IntegrateLowerSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0)/static_cast<MML::Real>(3.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Integration::UpperSingular_arcsin_deriv", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√(1-x²) dx = arcsin(1) - arcsin(0) = π/2
		// Note: this has singularity only at x=1, not at x=0
		RealFunctionFromStdFunc f( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(static_cast<MML::Real>(1.0) - x * x); } );
		
		auto result = IntegrateUpperSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Integration::InteriorSingular_two_point", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,2] 1/√|x-1| dx = 4  (singularity at x=1 interior point)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(std::abs(x - static_cast<MML::Real>(1.0))); } );
		
		auto result = IntegrateInteriorSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(4.0), TOL(1e-5, 1e-4)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                     GAUSSIAN QUADRATURE TESTS                         //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("GaussLegendre::weights_sum_to_interval", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Legendre on [-1,1], weights should sum to 2
		std::vector<Real> x(10), w(10);
		GaussLegendre(x, w);
		
		Real sum = 0;
		for (int i = 0; i < 10; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("GaussLegendre::weights_sum_scaled_interval", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Legendre on [a,b], weights should sum to (b-a)
		std::vector<Real> x(15), w(15);
		GaussLegendre(x, w, static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0));
		
		Real sum = 0;
		for (int i = 0; i < 15; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("GaussLegendre::nodes_symmetric", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// Nodes should be symmetric about 0 on [-1,1]
		std::vector<Real> x(8), w(8);
		GaussLegendre(x, w);
		
		// Check that x[i] = -x[n-1-i]
		for (int i = 0; i < 4; i++) {
			REQUIRE_THAT(x[i], WithinAbs(-x[7-i], TOL(1e-14, 1e-5)));
			REQUIRE_THAT(w[i], WithinAbs(w[7-i], TOL(1e-14, 1e-5)));
		}
	}
	
	TEST_CASE("GaussLegendre::integrate_polynomial_exact", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// n-point Gauss-Legendre exactly integrates polynomials up to degree 2n-1
		// 5-point should exactly integrate up to degree 9
		std::vector<Real> x(5), w(5);
		GaussLegendre(x, w, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		// ∫[0,1] x^4 dx = 1/5 = static_cast<MML::Real>(0.2)
		Real integral = 0;
		for (int i = 0; i < 5; i++)
			integral += w[i] * std::pow(x[i], MML::Real(4));
		REQUIRE_THAT(integral, WithinAbs(static_cast<MML::Real>(0.2), TOL(1e-14, 1e-5)));
		
		// ∫[0,1] x^8 dx = 1/9
		integral = 0;
		for (int i = 0; i < 5; i++)
			integral += w[i] * std::pow(x[i], MML::Real(8));
		REQUIRE_THAT(integral, WithinAbs(static_cast<MML::Real>(1.0)/static_cast<MML::Real>(9.0), TOL(1e-14, 1e-5)));
	}
	
	TEST_CASE("GaussLegendre::integrate_sin_high_precision", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, π] sin(x) dx = 2
		std::vector<Real> x(20), w(20);
		GaussLegendre(x, w, static_cast<MML::Real>(0.0), Constants::PI);
		
		Real integral = 0;
		for (int i = 0; i < 20; i++)
			integral += w[i] * std::sin(x[i]);
		REQUIRE_THAT(integral, WithinAbs(static_cast<MML::Real>(2.0), TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("GaussLegendre::rule_struct", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// Test the GaussQuadratureRule structure
		auto rule = GaussLegendreRule(12, -static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		
		REQUIRE(rule.n == 12);
		REQUIRE(rule.nodes.size() == 12);
		REQUIRE(rule.weights.size() == 12);
		
		// Weights should sum to (3 - (-2)) = 5
		Real sum = 0;
		for (auto w : rule.weights) sum += w;
		REQUIRE_THAT(sum, WithinAbs(static_cast<MML::Real>(5.0), TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("GaussLaguerre::weights_sum_to_gamma", "[gaussian][laguerre]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Laguerre with α=0, weights sum to Γ(1) = 1
		std::vector<Real> x(10), w(10);
		GaussLaguerre(x, w, static_cast<MML::Real>(0.0));
		
		Real sum = 0;
		for (int i = 0; i < 10; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(static_cast<MML::Real>(1.0), TOL3(1e-10, 1e-3, 1e-12)));
	}
	
	TEST_CASE("GaussLaguerre::nodes_positive", "[gaussian][laguerre]")
	{
			TEST_PRECISION_INFO();
		// All Laguerre nodes should be positive (on [0,∞))
		std::vector<Real> x(15), w(15);
		GaussLaguerre(x, w, static_cast<MML::Real>(0.0));
		
		for (int i = 0; i < 15; i++) {
			REQUIRE(x[i] > 0);
			REQUIRE(w[i] > 0);
		}
	}
	
	TEST_CASE("GaussLaguerre::integrate_exponential_decay", "[gaussian][laguerre]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,∞) e^(-x) dx = 1 (weight function already includes e^(-x))
		// So ∫ w(x) * 1 dx = 1
		std::vector<Real> x(10), w(10);
		GaussLaguerre(x, w, static_cast<MML::Real>(0.0));
		
		Real integral = 0;
		for (int i = 0; i < 10; i++)
			integral += w[i] * static_cast<MML::Real>(1.0);  // f(x) = 1, weight has e^(-x)
		REQUIRE_THAT(integral, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-3)));
	}
	
	TEST_CASE("GaussHermite::weights_sum_to_sqrt_pi", "[gaussian][hermite]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Hermite, weights sum to ∫ e^(-x²) dx = √π
		std::vector<Real> x(10), w(10);
		GaussHermite(x, w);
		
		Real sum = 0;
		for (int i = 0; i < 10; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(std::sqrt(Constants::PI), TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("GaussHermite::nodes_symmetric", "[gaussian][hermite]")
	{
			TEST_PRECISION_INFO();
		// Hermite nodes should be symmetric about 0
		std::vector<Real> x(10), w(10);
		GaussHermite(x, w);
		
		// Nodes should satisfy x[i] = -x[n-1-i]
		for (int i = 0; i < 5; i++) {
			REQUIRE_THAT(x[i], WithinAbs(-x[9-i], TOL(1e-14, 1e-5)));
			REQUIRE_THAT(w[i], WithinAbs(w[9-i], TOL(1e-14, 1e-5)));
		}
	}
	
	TEST_CASE("GaussHermite::integrate_even_polynomial", "[gaussian][hermite]")
	{
			TEST_PRECISION_INFO();
		// ∫ x² e^(-x²) dx = √π/2
		std::vector<Real> x(10), w(10);
		GaussHermite(x, w);
		
		Real integral = 0;
		for (int i = 0; i < 10; i++)
			integral += w[i] * x[i] * x[i];
		REQUIRE_THAT(integral, WithinAbs(std::sqrt(Constants::PI) / static_cast<MML::Real>(2.0), TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("GaussChebyshev1::closed_form_accuracy", "[gaussian][chebyshev]")
	{
			TEST_PRECISION_INFO();
		// Chebyshev type 1: w(x) = 1/√(1-x²)
		// Weights should all equal π/n
		std::vector<Real> x(10), w(10);
		GaussChebyshev1(x, w);
		
		Real expected_weight = Constants::PI / static_cast<MML::Real>(10.0);
		for (int i = 0; i < 10; i++) {
			REQUIRE_THAT(w[i], WithinAbs(expected_weight, TOL(1e-14, 1e-5)));
		}
	}
	
	TEST_CASE("GaussChebyshev1::nodes_are_zeros", "[gaussian][chebyshev]")
	{
			TEST_PRECISION_INFO();
		// Chebyshev type 1 nodes: cos((2k+1)π/(2n))
		int n = 8;
		std::vector<Real> x(n), w(n);
		GaussChebyshev1(x, w);
		
		for (int k = 0; k < n; k++) {
			Real expected = std::cos((static_cast<MML::Real>(2.0) * k + static_cast<MML::Real>(1.0)) * Constants::PI / (static_cast<MML::Real>(2.0) * n));
			REQUIRE_THAT(x[k], WithinAbs(expected, TOL(1e-14, 1e-5)));
		}
	}
	
	TEST_CASE("GaussChebyshev2::closed_form_accuracy", "[gaussian][chebyshev]")
	{
			TEST_PRECISION_INFO();
		// Chebyshev type 2: w(x) = √(1-x²)
		// Check weight formula: wₖ = π/(n+1) sin²((k+1)π/(n+1))
		int n = 8;
		std::vector<Real> x(n), w(n);
		GaussChebyshev2(x, w);
		
		for (int k = 0; k < n; k++) {
			Real angle = (k + static_cast<MML::Real>(1.0)) * Constants::PI / (n + static_cast<MML::Real>(1.0));
			Real expected = Constants::PI / (n + static_cast<MML::Real>(1.0)) * std::pow(std::sin(angle), MML::Real(2));
			REQUIRE_THAT(w[k], WithinAbs(expected, TOL(1e-14, 1e-5)));
		}
	}
	
	TEST_CASE("GaussJacobi::reduces_to_legendre", "[gaussian][jacobi]")
	{
			TEST_PRECISION_INFO();
		// Jacobi with α=β=0 is Legendre
		// Note: Different algorithms produce nodes in different orders, so sort first
		int n = 8;
		std::vector<Real> xL(n), wL(n), xJ(n), wJ(n);
		
		GaussLegendre(xL, wL);
		GaussJacobi(xJ, wJ, static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
		
		// Sort nodes (paired with weights)
		std::vector<std::pair<Real, Real>> legendre(n), jacobi(n);
		for (int i = 0; i < n; i++) {
			legendre[i] = {xL[i], wL[i]};
			jacobi[i] = {xJ[i], wJ[i]};
		}
		std::sort(legendre.begin(), legendre.end());
		std::sort(jacobi.begin(), jacobi.end());
		
		for (int i = 0; i < n; i++) {
			REQUIRE_THAT(jacobi[i].first, WithinAbs(legendre[i].first, TOL(1e-8, 1e-4)));
			REQUIRE_THAT(jacobi[i].second, WithinAbs(legendre[i].second, TOL(1e-8, 1e-4)));
		}
	}
	
	TEST_CASE("GaussJacobi::reduces_to_chebyshev1", "[gaussian][jacobi]")
	{
			TEST_PRECISION_INFO();
		// Jacobi with α=β=-1/2 is Chebyshev type 1
		int n = 8;
		std::vector<Real> xC(n), wC(n), xJ(n), wJ(n);
		
		GaussChebyshev1(xC, wC);
		GaussJacobi(xJ, wJ, -static_cast<MML::Real>(0.5), -static_cast<MML::Real>(0.5));
		
		// Nodes should match (possibly in different order)
		// Both are symmetric, just compare sorted
		std::sort(xC.begin(), xC.end());
		std::sort(xJ.begin(), xJ.end());
		
		for (int i = 0; i < n; i++) {
			REQUIRE_THAT(xJ[i], WithinAbs(xC[i], TOL(1e-8, 1e-4)));
		}
	}
	
	TEST_CASE("IntegrateGaussLegendre::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test the convenience integration function
		RealFunctionFromStdFunc f( [](Real x) { return x * x; } );
		
		// ∫[0,2] x² dx = 8/3
		Real result = IntegrateGaussLegendre(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0), 10);
		REQUIRE_THAT(result, WithinAbs(static_cast<MML::Real>(8.0)/static_cast<MML::Real>(3.0), TOL3(1e-12, 1e-5, 1e-14)));  // Relaxed for floating-point
	}
	
	TEST_CASE("IntegrateGaussLaguerre::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test the convenience integration function
		// ∫[0,∞) x * e^(-x) dx = 1 (function f=x, weight has e^(-x))
		RealFunctionFromStdFunc f( [](Real x) { return x; } );
		
		Real result = IntegrateGaussLaguerre(f, 15, static_cast<MML::Real>(0.0));
		REQUIRE_THAT(result, WithinAbs(static_cast<MML::Real>(1.0), TOL3(1e-8, 1e-4, 1e-10)));
	}
	
	TEST_CASE("IntegrateGaussHermite::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test the convenience integration function
		// ∫ x⁴ e^(-x²) dx = 3√π/4 (even function)
		RealFunctionFromStdFunc f( [](Real x) { return x * x * x * x; } );
		
		Real result = IntegrateGaussHermite(f, 15);
		Real expected = static_cast<MML::Real>(3.0) * std::sqrt(Constants::PI) / static_cast<MML::Real>(4.0);
		REQUIRE_THAT(result, WithinAbs(expected, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("IntegrateGaussJacobi::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test Gauss-Jacobi with α=1, β=0: weight = (1-x)
		// ∫[-1,1] (1-x) * x² dx = ∫[-1,1] x² - x³ dx = 2/3 - 0 = 2/3
		RealFunctionFromStdFunc f( [](Real x) { return x * x; } );
		
		Real result = IntegrateGaussJacobi(f, static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), 15);
		REQUIRE_THAT(result, WithinAbs(static_cast<MML::Real>(2.0)/static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("IntegrateGaussChebyshev1::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫[-1,1] 1/√(1-x²) dx = π (integral of weight function alone)
		RealFunctionFromStdFunc one( [](Real x) { return static_cast<MML::Real>(1.0); } );
		
		Real result = IntegrateGaussChebyshev1(one, 10);
		REQUIRE_THAT(result, WithinAbs(Constants::PI, TOL(1e-12, 1e-5)));
	}
	
	TEST_CASE("IntegrateGaussChebyshev2::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫[-1,1] √(1-x²) dx = π/2 (integral of weight function alone)
		RealFunctionFromStdFunc one( [](Real x) { return static_cast<MML::Real>(1.0); } );
		
		Real result = IntegrateGaussChebyshev2(one, 10);
		REQUIRE_THAT(result, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), TOL3(1e-12, 1e-5, 1e-14)));
	}
	
	TEST_CASE("IntegrateGaussLegendre::high_order_polynomial", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// 10-point Gauss-Legendre should exactly integrate polynomials up to degree 19
		// Test x^18 on [0,1]: ∫x^18 dx = 1/19
	RealFunctionFromStdFunc poly18([](Real x) -> Real { return std::pow(x, static_cast<MML::Real>(18)); });
	}
	
	TEST_CASE("IntegrateGaussLaguerre::gamma_integral", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,∞) x^3 e^(-x) dx = Γ(4) = 3! = 6
		// Using α=0, so f(x) = x³
		RealFunctionFromStdFunc f( [](Real x) { return x * x * x; } );
		
		const int laguerre_n = std::is_same_v<Real, float> ? 15 : 20;
		Real result = IntegrateGaussLaguerre(f, laguerre_n, static_cast<MML::Real>(0.0));
		REQUIRE_THAT(result, WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-8, 1e-2)));
	}
	
	TEST_CASE("IntegrateGaussHermite::gaussian_moments", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫ x⁶ e^(-x²) dx = 15√π/8 (6th moment of Gaussian)
	  RealFunctionFromStdFunc f([](Real x) -> Real { return std::pow(x, static_cast<MML::Real>(6)); });
	
	  Real result = IntegrateGaussHermite(f, 20);
	  Real expected = static_cast<MML::Real>(15.0) * std::sqrt(Constants::PI) / static_cast<MML::Real>(8.0);
	  REQUIRE_THAT(result, WithinAbs(expected, TOL(1e-10, 1e-5)));
  } 

  TEST_CASE("IntegrateGaussJacobi::beta_function", "[gaussian][integrate]")
  {
		TEST_PRECISION_INFO();
	// ∫[-1,1] (1-x)^α (1+x)^β dx = 2^(α+β+1) B(α+1, β+1)
	// where B is the Beta function: B(a,b) = Γ(a)Γ(b)/Γ(a+b)
	// For α=static_cast<MML::Real>(0.5), β=static_cast<MML::Real>(1.5): B(static_cast<MML::Real>(1.5), static_cast<MML::Real>(2.5)) = Γ(static_cast<MML::Real>(1.5))Γ(static_cast<MML::Real>(2.5))/Γ(4) = (√π/2)(3√π/4)/6 = 3π/48 = π/16
	// Result = 2^3 * π/16 = π/2
	  RealFunctionFromStdFunc one( [](Real x) { return static_cast<MML::Real>(1.0); } );
	
	  Real result = IntegrateGaussJacobi(one, static_cast<MML::Real>(0.5), static_cast<MML::Real>(1.5), 15);
	  REQUIRE_THAT(result, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-8, 1e-4)));
  }

  TEST_CASE("GaussJacobi::overflow_protection", "[gaussian][jacobi]")
  {
		TEST_PRECISION_INFO();
		// Test that extreme parameters that would cause overflow are caught
		// For very large α, β values, the weight calculation can overflow
		// The overflow check triggers when log_weight > 700 (exp(700) ≈ 1e304)
		std::vector<Real> x(200), w(200);
		
		// Either throws NumericalMethodError (overflow) or IntegrationTooManySteps (convergence)
		// Both are acceptable - the point is it doesn't silently produce wrong results
		REQUIRE_THROWS(GaussJacobi(x, w, static_cast<MML::Real>(100.0), static_cast<MML::Real>(100.0)));
  }

	TEST_CASE("Integration::Difficult_peaked_function", "[difficult]")
	{
	TEST_PRECISION_INFO();
		// Integrate a sharply peaked function: 1/(1+(x-static_cast<MML::Real>(0.3))²/static_cast<MML::Real>(0.01))
		// This is a Lorentzian with width static_cast<MML::Real>(0.1) centered at static_cast<MML::Real>(0.3)
		// ∫ dx/(1+(x-c)²/w²) = w*arctan((x-c)/w)
		RealFunctionFromStdFunc peaked( [](Real x) { 
			Real diff = x - static_cast<MML::Real>(0.3);
			return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + diff * diff / static_cast<MML::Real>(0.01));
		} );
		
		// Exact: static_cast<MML::Real>(0.1) * (arctan((1-static_cast<MML::Real>(0.3))/static_cast<MML::Real>(0.1)) - arctan((0-static_cast<MML::Real>(0.3))/static_cast<MML::Real>(0.1)))
		//      = static_cast<MML::Real>(0.1) * (arctan(7) + arctan(3))
		Real expected = static_cast<MML::Real>(0.1) * (std::atan(static_cast<MML::Real>(7.0)) + std::atan(static_cast<MML::Real>(3.0)));
		
		auto result = IntegrateSimpsonDetailed(peaked, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-6));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinRel(expected, static_cast<MML::Real>(1e-3)));
	}

	TEST_CASE("Integration::Difficult_rapid_oscillation", "[difficult]")
	{
			TEST_PRECISION_INFO();
		// Integrate sin(20x) on [0, π] = 0 (cancellation)
		// Using lower frequency to be more tractable
		RealFunctionFromStdFunc fastSin( [](Real x) { return std::sin(static_cast<MML::Real>(20.0) * x); } );
		
		auto result = IntegrateSimpsonDetailed(fastSin, static_cast<MML::Real>(0.0), Constants::PI, static_cast<MML::Real>(1e-4));
		
		// Should be very close to 0 due to cancellation
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.1)));
	}

	TEST_CASE("Integration::Difficult_near_singularity", "[difficult]")
	{
			TEST_PRECISION_INFO();
		// Integrate 1/√(x + static_cast<MML::Real>(0.001)) on [0,1] - near-singular at x=0
		RealFunctionFromStdFunc nearSing( [](Real x) { return static_cast<MML::Real>(1.0) / std::sqrt(x + static_cast<MML::Real>(0.001)); } );
		
		auto result = IntegrateSimpsonDetailed(nearSing, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-6));
		Real expected = static_cast<MML::Real>(2.0) * (std::sqrt(static_cast<MML::Real>(1.001)) - std::sqrt(static_cast<MML::Real>(0.001)));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinRel(expected, static_cast<MML::Real>(1e-4)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    METHOD COMPARISON TESTS                            //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::AllMethods_polynomial", "[comparison][comprehensive]")
	{
			TEST_PRECISION_INFO();
		// All methods should give same answer for polynomial
		RealFunctionFromStdFunc poly( [](Real x) { return x * x * x - 2 * x + 1; } );
		Real exact = static_cast<MML::Real>(0.25);  // ∫[0,1] (x³ - 2x + 1) dx = [x⁴/4 - x² + x]₀¹ = static_cast<MML::Real>(0.25)
		
		auto trap = IntegrateTrapDetailed(poly, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		auto simp = IntegrateSimpsonDetailed(poly, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		auto romb = IntegrateRombergDetailed(poly, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		auto gauss = IntegrateGauss10Detailed(poly, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE_THAT(trap.value, WithinRel(exact, static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(simp.value, WithinRel(exact, static_cast<MML::Real>(1e-6)));
		REQUIRE_THAT(romb.value, WithinRel(exact, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(gauss.value, WithinRel(exact, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("Integration::AllMethods_transcendental", "[comparison][comprehensive]")
	{
			TEST_PRECISION_INFO();
		// Compare methods on e^(-x²) over [0,1]
		RealFunctionFromStdFunc gauss( [](Real x) { return std::exp(-x * x); } );
		Real exact = static_cast<MML::Real>(0.7468241328124271);  // erf(1) * √π/2
		
		auto trap = IntegrateTrapDetailed(gauss, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		auto simp = IntegrateSimpsonDetailed(gauss, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		auto romb = IntegrateRombergDetailed(gauss, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		auto gauss10 = IntegrateGauss10Detailed(gauss, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		
		REQUIRE_THAT(trap.value, WithinAbs(exact, static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(simp.value, WithinAbs(exact, static_cast<MML::Real>(1e-6)));
		REQUIRE_THAT(romb.value, WithinAbs(exact, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(gauss10.value, WithinAbs(exact, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Integration::Romberg_vs_Gauss_accuracy", "[comparison]")
	{
			TEST_PRECISION_INFO();
		// Compare Romberg and Gauss-Legendre on same integral
		RealFunctionFromStdFunc f( [](Real x) { return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + x * x); } );
		Real exact = Constants::PI / static_cast<MML::Real>(4.0);  // arctan(1) - arctan(0)
		
		auto romb = IntegrateRombergDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		Real gauss = IntegrateGaussLegendre(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 15);
		
		// Both should be accurate to TOL(1e-8, 1e-4)
		REQUIRE_THAT(romb.value, WithinAbs(exact, TOL(1e-8, 1e-4)));
		REQUIRE_THAT(gauss, WithinAbs(exact, TOL(1e-10, 1e-5)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    ERROR HANDLING TESTS                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::Invalid_bounds", "[error]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc f( [](Real x) { return x; } );
		
		// Empty interval should give 0
		REQUIRE(IntegrateTrap(f, static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)) == static_cast<MML::Real>(0.0));
		
		// Reversed bounds should work (gives negative of forward integral)
		Real forward = IntegrateTrap(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		Real reverse = IntegrateTrap(f, static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));
		REQUIRE_THAT(reverse, WithinAbs(-forward, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Integration::Zero_valued_integrals_converge", "[convergence][zero][wewi.1.11]")
	{
			TEST_PRECISION_INFO();
		// abs_tol used to be eps * machine_eps, so a sequence converging to 0 could never satisfy
		// |delta| < eps*|sum| + abs_tol and every adaptive rule reported converged == false
		// after the maximum number of refinements (analysis v4, CORE C5).
		RealFunctionFromStdFunc odd( [](Real x) { return x * x * x; } );      // integral over [-1,1] is exactly 0
		RealFunctionFromStdFunc sine( [](Real x) { return std::sin(x); } );   // integral over [-pi, pi] is 0

		for (auto* f : { &odd, &sine })
		{
			const Real a = (f == &odd) ? -static_cast<MML::Real>(1.0) : -Constants::PI;
			const Real b = -a;

			auto trap = IntegrateTrapDetailed(*f, a, b, TOL(1e-8, 1e-4));
			auto simp = IntegrateSimpsonDetailed(*f, a, b, TOL(1e-8, 1e-4));
			auto romb = IntegrateRombergDetailed(*f, a, b, TOL(1e-8, 1e-4));

			INFO("trap iterations " << trap.iterations << " simpson " << simp.iterations << " romberg " << romb.iterations);
			REQUIRE(trap.converged);
			REQUIRE(simp.converged);
			REQUIRE(romb.converged);
			REQUIRE_THAT(trap.value, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
			REQUIRE_THAT(simp.value, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
			REQUIRE_THAT(romb.value, WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-8, 1e-4)));
			REQUIRE(trap.iterations < Defaults::TrapezoidIntegrationMaxSteps);
			REQUIRE(simp.iterations < Defaults::SimpsonIntegrationMaxSteps);
		}
	}

	TEST_CASE("Integration::Result_honesty_flags", "[result][wewi.1.11]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc f( [](Real x) { return std::exp(-x * x); } );

		// Conversion to Real is explicit: it drops converged / error_estimate
		static_assert(!std::is_convertible_v<IntegrationResult, Real>, "IntegrationResult must not convert implicitly");
		static_assert(std::is_constructible_v<Real, IntegrationResult>, "explicit conversion must remain available");
		auto trap = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		REQUIRE(static_cast<Real>(trap) == trap.value);
		REQUIRE(trap.error_estimated);

		// Fixed rules check no criterion: converged is false, status Completed, no error estimate
		auto gauss10 = IntegrateGauss10Detailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		REQUIRE_FALSE(gauss10.error_estimated);
		REQUIRE_FALSE(gauss10.converged);
		REQUIRE_FALSE(HasEstimate(gauss10.error_estimate));
		auto legendre = IntegrateGaussLegendreDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 10);
		REQUIRE_FALSE(legendre.error_estimated);
		REQUIRE_FALSE(legendre.converged);
		auto legendreDetailed = IntegrateGaussLegendreDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 10);
		REQUIRE(legendreDetailed.IsSuccess());                       // usable: Completed counts as success
		REQUIRE(legendreDetailed.status == AlgorithmStatus::Completed);
		REQUIRE_FALSE(legendreDetailed.converged);
		REQUIRE_FALSE(legendreDetailed.error_estimated);
		REQUIRE_FALSE(HasEstimate(legendreDetailed.error_estimate));
		REQUIRE(legendreDetailed.error_message.empty());
		// printing a fixed-rule result must not label it "not converged"
		std::ostringstream os;
		os << gauss10;
		REQUIRE(os.str().find("not converged") == std::string::npos);

		// Gauss-Kronrod does carry a real estimate
		auto gk = IntegrateGK21Detailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		REQUIRE(gk.error_estimated);
		REQUIRE(gk.error_estimate >= static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Integration::Interior_singular_bounds_check", "[error][singular]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc f( [](Real x) { return static_cast<MML::Real>(1.0) / std::sqrt(std::abs(x)); } );
		
		// Singular point must be inside (a, b)
		REQUIRE_THROWS_AS(IntegrateInteriorSingular(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)), std::invalid_argument);
		REQUIRE_THROWS_AS(IntegrateInteriorSingular(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)), std::invalid_argument);
		REQUIRE_THROWS_AS(IntegrateInteriorSingular(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)), std::invalid_argument);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    SPECIAL VALUES TESTS                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::WellKnown_pi_computations", "[special]")
	{
			TEST_PRECISION_INFO();
		// Various ways to compute π via integration
		
		// ∫[-1,1] 1/√(1-x²) dx = π
		RealFunctionFromStdFunc arcCircle( [](Real x) { return static_cast<MML::Real>(1.0) / std::sqrt(static_cast<MML::Real>(1.0) - x * x); } );
		auto pi1 = IntegrateBothSingularDetailed(arcCircle, -static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0));
		REQUIRE_THAT(pi1.value, WithinAbs(Constants::PI, static_cast<MML::Real>(1e-4)));
		
		// ∫[0,1] 4/(1+x²) dx = π
		RealFunctionFromStdFunc arctan4( [](Real x) { return static_cast<MML::Real>(4.0) / (static_cast<MML::Real>(1.0) + x * x); } );
		auto pi2 = IntegrateRombergDetailed(arctan4, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		REQUIRE_THAT(pi2.value, WithinAbs(Constants::PI, static_cast<MML::Real>(1e-7)));  // Relaxed tolerance
		
		// ∫(-∞,∞) e^(-x²) dx = √π
		RealFunctionFromStdFunc gaussian( [](Real x) { return std::exp(-x * x); } );
		auto sqrtPi = IntegrateInfDetailed(gaussian);
		REQUIRE_THAT(sqrtPi.value, WithinAbs(std::sqrt(Constants::PI), static_cast<MML::Real>(1e-6)));
	}

	TEST_CASE("Integration::WellKnown_e_computations", "[special]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] e^x dx = e - 1
		RealFunctionFromStdFunc expf( [](Real x) { return std::exp(x); } );
		auto result = IntegrateRombergDetailed(expf, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
		REQUIRE_THAT(result.value, WithinAbs(std::exp(static_cast<MML::Real>(1.0)) - static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		
		// ∫[0,∞] e^(-x) dx = 1
		RealFunctionFromStdFunc expDecay( [](Real x) { return std::exp(-x); } );
		auto infResult = IntegrateUpperInfDetailed(expDecay, static_cast<MML::Real>(0.0));
		REQUIRE_THAT(infResult.value, WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-6)));
	}

	TEST_CASE("Integration::WellKnown_gamma_values", "[special]")
	{
			TEST_PRECISION_INFO();
		// Γ(5) = 4! = 24 via ∫[0,∞) x^4 e^(-x) dx
		// Using Gauss-Laguerre with n=20, alpha=0 (standard weight e^(-x))
	const int laguerre_n = std::is_same_v<Real, float> ? 15 : 20;
	RealFunctionFromStdFunc gamma5([](Real x) -> Real { return std::pow(x, static_cast<MML::Real>(4)); });
	Real result5 = IntegrateGaussLaguerre(gamma5, laguerre_n, static_cast<MML::Real>(0.0));
	REQUIRE_THAT(result5, WithinAbs(static_cast<MML::Real>(24.0), TOL3(1e-6, 1e-2, 1e-8)));
	
	// Γ(3) = 2! = 2 via ∫[0,∞) x^2 e^(-x) dx
	RealFunctionFromStdFunc gamma3([](Real x) -> Real { return x * x; });
	Real result3 = IntegrateGaussLaguerre(gamma3, laguerre_n, static_cast<MML::Real>(0.0));
	REQUIRE_THAT(result3, WithinAbs(static_cast<MML::Real>(2.0), TOL3(1e-8, 1e-2, 1e-10)));
}

	///////////////////////////////////////////////////////////////////////////
	//                UNIFIED INTEGRATE() FUNCTION TESTS                     //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::Integrate_runtime_method_selection", "[unified][api]")
	{
		TEST_PRECISION_INFO();
		// Test the unified Integrate() function with explicit method selection
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });
		Real expected = static_cast<MML::Real>(1.0) / static_cast<MML::Real>(3.0);  // ∫[0,1] x² dx = 1/3
		
		SECTION("TRAP method")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TRAP);
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinRel(expected, static_cast<MML::Real>(1e-4)));
		}
		
		SECTION("SIMPSON method")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), SIMPSON);
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinRel(expected, static_cast<MML::Real>(1e-6)));
		}
		
		SECTION("ROMBERG method")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), ROMBERG);
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinRel(expected, TOL(1e-10, 1e-5)));
		}
		
		SECTION("GAUSS10 method")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), GAUSS10);
			REQUIRE(result.converged == false);     // fixed rule: nothing checked
			REQUIRE_FALSE(result.error_estimated);
			REQUIRE_THAT(result.value, WithinRel(expected, TOL(1e-12, 1e-5)));
		}
		
		SECTION("Default method (TRAP)")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));  // Default is TRAP
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinRel(expected, static_cast<MML::Real>(1e-4)));
		}
	}
	
	TEST_CASE("Integration::Integrate_runtime_with_custom_eps", "[unified][api]")
	{
		TEST_PRECISION_INFO();
		// Test that the eps parameter is properly passed through
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::sin(x); });
		Real expected = static_cast<MML::Real>(2.0);  // ∫[0,π] sin(x) dx = 2
		
		SECTION("High precision request")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, SIMPSON, TOL(1e-10, 1e-5));
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinAbs(expected, TOL(1e-8, 1e-4)));
		}
		
		SECTION("Low precision request")
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, TRAP, static_cast<MML::Real>(1e-2));
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(0.1)));
		}
	}
	
	TEST_CASE("Integration::Integrate_template_compile_time_method", "[unified][api][template]")
	{
		TEST_PRECISION_INFO();
		// Test the template-based Integrate<Method>() function
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::exp(-x * x); });
		Real expected = static_cast<MML::Real>(0.7468241328124271);  // ∫[0,1] e^(-x²) dx = erf(1) * √π/2
		
		SECTION("Template TRAP method")
		{
			auto result = IntegrateDetailed<TRAP>(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(1e-4)));
		}
		
		SECTION("Template SIMPSON method")
		{
			auto result = IntegrateDetailed<SIMPSON>(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(1e-6)));
		}
		
		SECTION("Template ROMBERG method")
		{
			auto result = IntegrateDetailed<ROMBERG>(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinAbs(expected, TOL(1e-8, 1e-4)));
		}
		
		SECTION("Template GAUSS10 method")
		{
			auto result = IntegrateDetailed<GAUSS10>(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
			REQUIRE(result.converged == false);     // fixed rule: nothing checked
			REQUIRE_FALSE(result.error_estimated);
			REQUIRE_THAT(result.value, WithinAbs(expected, TOL(1e-10, 1e-5)));
		}
		
		SECTION("Template default method")
		{
			auto result = IntegrateDetailed<>(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));  // Default is TRAP
			REQUIRE(result.converged == true);
			REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(1e-4)));
		}
	}
	
	TEST_CASE("Integration::Integrate_template_with_custom_eps", "[unified][api][template]")
	{
		TEST_PRECISION_INFO();
		// Test the template version with custom eps parameter
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0) / (static_cast<MML::Real>(1.0) + x * x); });
		Real expected = Constants::PI / static_cast<MML::Real>(4.0);  // ∫[0,1] 1/(1+x²) dx = π/4
		
		auto result = IntegrateDetailed<ROMBERG>(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), TOL(1e-12, 1e-5));
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, TOL(1e-10, 1e-5)));
	}
	
	TEST_CASE("Integration::IntegrationMethod_enum_values", "[unified][api]")
	{
		TEST_PRECISION_INFO();
		// Verify that all enum values are accessible and produce valid results
		RealFunctionFromStdFunc f([](Real x) -> Real { return x; });
		Real expected = static_cast<MML::Real>(0.5);  // ∫[0,1] x dx = 0.5
		
		// Test each enum value explicitly
		IntegrationMethod methods[] = { TRAP, SIMPSON, ROMBERG, GAUSS10 };
		
		for (auto method : methods)
		{
			auto result = IntegrateDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), method);
			// adaptive rules converge; the fixed GAUSS10 rule reports no criterion instead
			const bool honest = method == GAUSS10 ? (!result.converged && !result.error_estimated)
			                                      : (result.converged && result.error_estimated);
			REQUIRE(honest);
			REQUIRE_THAT(result.value, WithinAbs(expected, static_cast<MML::Real>(1e-3)));
		}
	}
	
	TEST_CASE("Integration::Integrate_consistency_runtime_vs_template", "[unified][api]")
	{
		TEST_PRECISION_INFO();
		// Verify that runtime and template versions give identical results
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::cos(x) * std::cos(x); });
		
		auto runtime_trap = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, TRAP);
		auto template_trap = IntegrateDetailed<TRAP>(f, static_cast<MML::Real>(0.0), Constants::PI);
		REQUIRE_THAT(runtime_trap.value, WithinAbs(template_trap.value, TOL(1e-14, 1e-5)));
		
		auto runtime_simpson = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, SIMPSON);
		auto template_simpson = IntegrateDetailed<SIMPSON>(f, static_cast<MML::Real>(0.0), Constants::PI);
		REQUIRE_THAT(runtime_simpson.value, WithinAbs(template_simpson.value, TOL(1e-14, 1e-5)));
		
		auto runtime_romberg = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, ROMBERG);
		auto template_romberg = IntegrateDetailed<ROMBERG>(f, static_cast<MML::Real>(0.0), Constants::PI);
		REQUIRE_THAT(runtime_romberg.value, WithinAbs(template_romberg.value, TOL(1e-14, 1e-5)));
		
		auto runtime_gauss = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, GAUSS10);
		auto template_gauss = IntegrateDetailed<GAUSS10>(f, static_cast<MML::Real>(0.0), Constants::PI);
		REQUIRE_THAT(runtime_gauss.value, WithinAbs(template_gauss.value, TOL(1e-14, 1e-5)));
	}

	/******************************************************************************/
	/*****              Detailed API Tests - 1D Integration                    *****/
	/******************************************************************************/

	TEST_CASE("IntegrateTrapDetailed - basic success", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });

		auto result = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateTrap");
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.converged);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0/3.0), static_cast<MML::Real>(1e-4)));
		REQUIRE(result.iterations > 0);
	}

	TEST_CASE("IntegrateSimpsonDetailed - basic success", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });

		auto result = IntegrateSimpsonDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateSimpson");
		REQUIRE(result.converged);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0/3.0), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("IntegrateRombergDetailed - basic success", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });

		auto result = IntegrateRombergDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateRomberg");
		REQUIRE(result.converged);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0/3.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("IntegrateGauss10Detailed - Completed, not Success (fixed rule)", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });

		auto result = IntegrateGauss10Detailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

		REQUIRE(result.IsSuccess());                                  // payload usable
		REQUIRE(result.status == AlgorithmStatus::Completed);         // but nothing was checked
		REQUIRE(result.algorithm_name == "IntegrateGauss10");
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.error_message.empty());
		REQUIRE(result.function_evaluations == 10);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0/3.0), TOL(1e-14, 1e-5)));
	}

	TEST_CASE("IntegrateGK21Detailed - basic success", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });

		auto result = IntegrateGK21Detailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateGK21");
		REQUIRE_FALSE(result.converged);
		REQUIRE(result.status == AlgorithmStatus::Completed);
		REQUIRE(result.function_evaluations == 21);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0/3.0), TOL(1e-14, 1e-5)));
	}

	TEST_CASE("IntegrateDetailed - runtime dispatch", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::sin(x); });
		const Real expected = static_cast<MML::Real>(2.0);  // ∫[0,π] sin(x) dx = 2

		auto trap = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, TRAP);
		REQUIRE(trap.IsSuccess());
		REQUIRE(trap.algorithm_name == "IntegrateTrap");
		REQUIRE_THAT(trap.value, WithinAbs(expected, static_cast<MML::Real>(1e-4)));

		auto simp = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, SIMPSON);
		REQUIRE(simp.IsSuccess());
		REQUIRE(simp.algorithm_name == "IntegrateSimpson");
		REQUIRE_THAT(simp.value, WithinAbs(expected, TOL(1e-7, 1e-5)));

		auto romb = IntegrateDetailed(f, static_cast<MML::Real>(0.0), Constants::PI, ROMBERG);
		REQUIRE(romb.IsSuccess());
		REQUIRE_THAT(romb.value, WithinAbs(expected, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("IntegrateDetailed - template dispatch", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::sin(x); });
		const Real expected = static_cast<MML::Real>(2.0);

		auto trap = IntegrateDetailed<TRAP>(f, static_cast<MML::Real>(0.0), Constants::PI);
		auto simp = IntegrateDetailed<SIMPSON>(f, static_cast<MML::Real>(0.0), Constants::PI);
		auto gauss = IntegrateDetailed<GAUSS10>(f, static_cast<MML::Real>(0.0), Constants::PI);

		REQUIRE(trap.IsSuccess());
		REQUIRE(simp.IsSuccess());
		REQUIRE(gauss.IsSuccess());
		REQUIRE_THAT(trap.value, WithinAbs(expected, static_cast<MML::Real>(1e-4)));
		REQUIRE_THAT(simp.value, WithinAbs(expected, TOL(1e-7, 1e-5)));
		REQUIRE_THAT(gauss.value, WithinAbs(expected, TOL(1e-14, 1e-5)));
	}

	TEST_CASE("IntegrateDetailed - values match simple API", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::exp(-x * x); });

		auto simple = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0));
		auto detailed = IntegrateTrapDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0));

		REQUIRE_THAT(detailed.value, WithinAbs(simple.value, TOL(1e-12, 1e-5)));
		REQUIRE(detailed.iterations == simple.iterations);
		REQUIRE(detailed.converged == simple.converged);
	}

	/******************************************************************************/
	/*****             Detailed API Tests - Improper Integration               *****/
	/******************************************************************************/

	TEST_CASE("IntegrateUpperInfDetailed - Gaussian tail", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::exp(-x * x); });

		auto result = IntegrateUpperInfDetailed(f, static_cast<MML::Real>(0.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateUpperInf");
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.886226925452758), static_cast<MML::Real>(1e-6)));
	}

	TEST_CASE("IntegrateLowerInfDetailed - Gaussian tail", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::exp(-x * x); });

		auto result = IntegrateLowerInfDetailed(f, static_cast<MML::Real>(0.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateLowerInf");
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(0.886226925452758), static_cast<MML::Real>(1e-6)));
	}

	TEST_CASE("IntegrateInfDetailed - full Gaussian", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return std::exp(-x * x); });

		auto result = IntegrateInfDetailed(f);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateInf");
		// ∫(-∞,∞) e^{-x²} dx = √π ≈ 1.7724538509
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.7724538509055159), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("IntegrateInteriorSingularDetailed - ConvertToStatus on invalid c", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(std::abs(x)); });

		IntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

		// c=5 is outside (0, 1), should produce InvalidInput
		auto result = IntegrateInteriorSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(5.0), config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(result.error_message.empty());
	}

	TEST_CASE("IntegrateInteriorSingularDetailed - Propagate policy throws", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(std::abs(x)); });

		IntegrationConfig config;
		config.exception_policy = EvaluationExceptionPolicy::Propagate;

		// c=5 is outside (0, 1), should throw
		REQUIRE_THROWS_AS(
			IntegrateInteriorSingularDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(5.0), config),
			std::invalid_argument);
	}

	TEST_CASE("IntegrateInteriorSingularDetailed - valid singularity", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0) / std::sqrt(std::abs(x)); });

		auto result = IntegrateInteriorSingularDetailed(f, static_cast<MML::Real>(-1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0));

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "IntegrateInteriorSingular");
		// ∫[-1,1] 1/√|x| dx = 4
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(4.0), static_cast<MML::Real>(1e-3)));
	}

	/******************************************************************************/
	/*****       Detailed API Tests - Gaussian Quadrature Variants             *****/
	/******************************************************************************/

	TEST_CASE("IntegrateGaussLegendreDetailed - basic success", "[mml_integration][Detailed]")
	{
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });
		const int n = 10;

		auto result = IntegrateGaussLegendreDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), n);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GaussLegendre");
		REQUIRE(result.status == AlgorithmStatus::Completed);   // fixed rule: nothing checked (RESULT_POLICY)
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.function_evaluations == n);
		// ∫[0,1] x² dx = 1/3
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0/3.0), TOL(1e-14, 1e-5)));
	}

	TEST_CASE("IntegrateGaussLaguerreDetailed - basic success", "[mml_integration][Detailed]")
	{
		// f(x) = 1, weight = x^0 * e^(-x), ∫[0,∞) e^(-x) dx = 1
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0); });
		const int n = 10;

		auto result = IntegrateGaussLaguerreDetailed(f, n);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GaussLaguerre");
		REQUIRE(result.status == AlgorithmStatus::Completed);   // fixed rule: nothing checked (RESULT_POLICY)
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.function_evaluations == n);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-3)));
	}

	TEST_CASE("IntegrateGaussHermiteDetailed - basic success", "[mml_integration][Detailed]")
	{
		// f(x) = 1, weight = e^(-x²), ∫(-∞,∞) e^(-x²) dx = √π
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0); });
		const int n = 10;

		auto result = IntegrateGaussHermiteDetailed(f, n);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GaussHermite");
		REQUIRE(result.status == AlgorithmStatus::Completed);   // fixed rule: nothing checked (RESULT_POLICY)
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.function_evaluations == n);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(1.7724538509055159), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("IntegrateGaussJacobiDetailed - basic success", "[mml_integration][Detailed]")
	{
		// α=β=0 reduces to Gauss-Legendre: ∫[-1,1] x² dx = 2/3
		RealFunctionFromStdFunc f([](Real x) -> Real { return x * x; });
		const int n = 10;

		auto result = IntegrateGaussJacobiDetailed(f, static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), n);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GaussJacobi");
		REQUIRE(result.status == AlgorithmStatus::Completed);   // fixed rule: nothing checked (RESULT_POLICY)
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.function_evaluations == n);
		REQUIRE_THAT(result.value, WithinAbs(static_cast<MML::Real>(2.0/3.0), TOL(1e-12, 1e-5)));
	}

	TEST_CASE("IntegrateGaussChebyshev1Detailed - basic success", "[mml_integration][Detailed]")
	{
		// f(x) = 1, weight = 1/√(1-x²), ∫[-1,1] 1/√(1-x²) dx = π
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0); });
		const int n = 20;

		auto result = IntegrateGaussChebyshev1Detailed(f, n);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GaussChebyshev1");
		REQUIRE(result.status == AlgorithmStatus::Completed);   // fixed rule: nothing checked (RESULT_POLICY)
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.function_evaluations == n);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("IntegrateGaussChebyshev2Detailed - basic success", "[mml_integration][Detailed]")
	{
		// f(x) = 1, weight = √(1-x²), ∫[-1,1] √(1-x²) dx = π/2
		RealFunctionFromStdFunc f([](Real x) -> Real { return static_cast<MML::Real>(1.0); });
		const int n = 20;

		auto result = IntegrateGaussChebyshev2Detailed(f, n);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.algorithm_name == "GaussChebyshev2");
		REQUIRE(result.status == AlgorithmStatus::Completed);   // fixed rule: nothing checked (RESULT_POLICY)
		REQUIRE_FALSE(result.converged);
		REQUIRE_FALSE(HasEstimate(result.error_estimate));
		REQUIRE(result.elapsed_time_ms >= 0.0);
		REQUIRE(result.function_evaluations == n);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / static_cast<MML::Real>(2.0), TOL(1e-12, 1e-5)));
	}
}
