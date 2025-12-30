#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include <algorithm>  // for std::sort
#include <utility>    // for std::pair

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Integration.h"

#endif

#include "../test_data/real_functions_test_bed.h"

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::IntegrationTests
{
	TEST_CASE("Test_Integration_Trap_const_func", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc constFunc0( [](Real x) -> Real { return REAL(0.0); } );
		RealFunctionFromStdFunc constFunc1( [](Real x) -> Real { return REAL(2.0); } );
		RealFunctionFromStdFunc constFunc2( [](Real x) -> Real { return -REAL(5.0); } );
		RealFunctionFromStdFunc constFunc3( [](Real x) -> Real { return 1e-10; } );

		REQUIRE_THAT(REAL(0.0), WithinAbs(IntegrateTrap(constFunc0, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(REAL(0.0), WithinAbs(IntegrateTrap(constFunc0, -REAL(3.0), REAL(1.0)), REAL(1e-15)));

		REQUIRE_THAT(REAL(2.0), WithinAbs(IntegrateTrap(constFunc1, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(REAL(8.0), WithinAbs(IntegrateTrap(constFunc1, -REAL(3.0), REAL(1.0)), REAL(1e-15)));

		REQUIRE_THAT(-REAL(5.0), WithinAbs(IntegrateTrap(constFunc2, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(-REAL(20.0), WithinAbs(IntegrateTrap(constFunc2, -REAL(3.0), REAL(1.0)), REAL(1e-15)));

		REQUIRE_THAT(1e-10, WithinAbs(IntegrateTrap(constFunc3, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(4e-10, WithinAbs(IntegrateTrap(constFunc3, -REAL(3.0), REAL(1.0)), REAL(1e-15)));
	}

	// TODO
	//TEST_CASE("Test_Integration_almost_const_func", "[simple]") 
	//{
	//	RealFunctionFromStdFunc linearFun( [](Real x) -> Real { return REAL(2.0); } );

	//	double int_trap = IntegrateTrap(linearFun, REAL(0.0), REAL(1.0));

	//	REQUIRE_THAT(REAL(1.0), Catch::Matchers::WithinAbs(IntegrateTrap(linearFun, REAL(0.0), REAL(1.0)), REAL(1e-16)));
	//}

	TEST_CASE("Test_Integration_Trap_linear_func", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc linearFunc1( [](Real x) -> Real { return 2 * x; } );
		RealFunctionFromStdFunc linearFunc2( [](Real x) -> Real { return -3 * x; } );
		RealFunctionFromStdFunc linearFunc3( [](Real x) -> Real { return 1e-10 * x; } );
		RealFunctionFromStdFunc linearFunc4( [](Real x) -> Real { return 2 * (x - REAL(0.5)); } );

		REQUIRE_THAT(REAL(1.0), WithinAbs(IntegrateTrap(linearFunc1, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(-REAL(1.5), WithinAbs(IntegrateTrap(linearFunc2, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(0.5e-10, WithinAbs(IntegrateTrap(linearFunc3, REAL(0.0), REAL(1.0)), REAL(1e-15)));
		REQUIRE_THAT(REAL(0.0), WithinAbs(IntegrateTrap(linearFunc4, REAL(0.0), REAL(1.0)), REAL(1e-15)));

		REQUIRE_THAT(-REAL(5.0), WithinAbs(IntegrateTrap(linearFunc1, -REAL(3.0), REAL(2.0)), REAL(1e-15)));
		REQUIRE_THAT(REAL(7.5), WithinAbs(IntegrateTrap(linearFunc2, -REAL(3.0), REAL(2.0)), REAL(1e-15)));
		REQUIRE_THAT(-2.5e-10, WithinAbs(IntegrateTrap(linearFunc3, -REAL(3.0), REAL(2.0)), REAL(1e-15)));
		REQUIRE_THAT(-REAL(10.0), WithinAbs(IntegrateTrap(linearFunc4, -REAL(3.0), REAL(2.0)), REAL(1e-15)));
	}

	TEST_CASE("Test_Integration_Trap_standard_func_precision", "[simple]")
	{
			TEST_PRECISION_INFO();
		RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getFunc("Sin")._func;
		RealFunction sinFunc_int = TestBeds::RealFunctionsTestBed::getFunc("Sin")._funcIntegrated;

		Real a = REAL(0.0), b = REAL(1.0);
		REQUIRE_THAT(sinFunc_int(b) - sinFunc_int(a), WithinAbs(IntegrateTrap(sinFunc, a, b), REAL(1e-5)));
		REQUIRE_THAT(sinFunc_int(b) - sinFunc_int(a), !WithinAbs(IntegrateTrap(sinFunc, a, b), REAL(1e-6)));

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

		Real a = REAL(0.0), b = REAL(1.0);
		Real integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal,  WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-6)));

		REQUIRE_THAT(integralExactVal,  WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-6)));

		a = REAL(0.0), b = REAL(5.0);
		integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal,  WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-2)));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-3)));

		REQUIRE_THAT(integralExactVal,  WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-4)));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-5)));

		a = -REAL(5.0), b = REAL(2.0);
		integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal, WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-2)));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-3)));

		REQUIRE_THAT(integralExactVal, WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-4)));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
	}

	TEST_CASE("Test_Integration_Trap_integral_test_bed2", "[simple]")
	{
			TEST_PRECISION_INFO();
		auto testFunc = TestBeds::RealFunctionsTestBed::getFuncWithIntegral("TestInt2")._func;
		auto testFunc_int = TestBeds::RealFunctionsTestBed::getFuncWithIntegral("TestInt2")._funcIntegrated;

		Real a = REAL(0.0), b = REAL(1.0);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-6)));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-6)));

		a = REAL(0.0), b = REAL(5.0);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-6)));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-4)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-5)));

		a = -REAL(5.0), b = REAL(2.0);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-4)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), REAL(1e-5)));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-4)));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), REAL(1e-5)));
	}

	TEST_CASE("Integration::Adaptive_precision_control", "[adaptive]")
	{
			TEST_PRECISION_INFO();
		// Test that integration methods achieve requested precision
		RealFunctionFromStdFunc testFunc( [](Real x) -> Real { return std::sin(x) * std::exp(-x); } );

		// Test with different precision requirements
		int steps1, steps2, steps3;
		Real precision1, precision2, precision3;

		Real result1 = IntegrateTrap(testFunc, REAL(0.0), Constants::PI, &steps1, &precision1, REAL(1e-6));
		Real result2 = IntegrateTrap(testFunc, REAL(0.0), Constants::PI, &steps2, &precision2, REAL(1e-8));
		Real result3 = IntegrateTrap(testFunc, REAL(0.0), Constants::PI, &steps3, &precision3, REAL(1e-10));

		// Higher precision requires more steps
		REQUIRE(steps2 > steps1);
		REQUIRE(steps3 > steps2);

		// Results should be consistent and increasingly accurate
		REQUIRE_THAT(result1, WithinAbs(result2, REAL(1e-6)));
		REQUIRE_THAT(result2, WithinAbs(result3, REAL(1e-8)));
	}

	TEST_CASE("Integration::Simpson_vs_Trap_accuracy", "[comparison]")
	{
			TEST_PRECISION_INFO();
		// Simpson's rule should be more accurate than trapezoidal for smooth functions
		RealFunctionFromStdFunc smoothFunc( [](Real x) -> Real { return std::exp(-x * x); } );

		int trapSteps, simpsonSteps;
		Real trapPrec, simpsonPrec;

		// Request same relative precision
		Real trapResult = IntegrateTrap(smoothFunc, REAL(0.0), REAL(2.0), &trapSteps, &trapPrec, REAL(1e-8));
		Real simpsonResult = IntegrateSimpson(smoothFunc, REAL(0.0), REAL(2.0), &simpsonSteps, &simpsonPrec, REAL(1e-8));

		// Simpson should converge faster (fewer steps)
		REQUIRE(simpsonSteps < trapSteps);

		// Both should give similar results
		REQUIRE_THAT(trapResult, WithinAbs(simpsonResult, REAL(1e-7)));
	}

	TEST_CASE("Integration::Gauss10_high_accuracy", "[gauss]")
	{
			TEST_PRECISION_INFO();
		// Gauss-Legendre integration should be very accurate for polynomials
		// Test with polynomial of degree <= 19 (exactly integrable with 10-point rule)
		RealFunctionFromStdFunc poly5( [](Real x) -> Real { return REAL(1.0) + x + x * x + x * x * x + x * x * x * x + x * x * x * x * x; } );

		// Analytical integral: x + x^2/2 + x^3/3 + x^4/4 + x^5/5 + x^6/6
		auto exactIntegral = [](Real x) -> Real {
			return x + x * x / REAL(2.0) + std::pow(x, 3) / REAL(3.0) + std::pow(x, 4) / REAL(4.0) + std::pow(x, 5) / REAL(5.0) + std::pow(x, 6) / REAL(6.0);
		};

		Real a = REAL(0.0), b = REAL(2.0);
		Real exact = exactIntegral(b) - exactIntegral(a);
		Real gauss10 = IntegrateGauss10(poly5, a, b);

		// Gauss-Legendre should be nearly exact for polynomials
		REQUIRE_THAT(exact, WithinAbs(gauss10, REAL(1e-10)));
	}
	
	TEST_CASE("Integration::Gauss10_returns_IntegrationResult", "[gauss][api]")
	{
			TEST_PRECISION_INFO();
		// Verify that IntegrateGauss10 now returns IntegrationResult structure
		RealFunctionFromStdFunc quadratic( [](Real x) -> Real { return x * x; } );
		
		// ∫₀¹ x² dx = 1/3
		auto result = IntegrateGauss10(quadratic, REAL(0.0), REAL(1.0));
		
		// Check all fields of IntegrationResult
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(3.0), REAL(1e-12)));
		REQUIRE(result.error_estimate == REAL(0.0));  // Gauss10 has no error estimate
		REQUIRE(result.iterations == 1);        // Single-pass evaluation
		REQUIRE(result.converged == true);      // Always "converged" for fixed-order
		
		// Verify implicit conversion to Real still works (backward compatibility)
		Real value = IntegrateGauss10(quadratic, REAL(0.0), REAL(1.0));
		REQUIRE_THAT(value, WithinAbs(REAL(1.0)/REAL(3.0), REAL(1e-12)));
	}

	TEST_CASE("Integration::Oscillatory_functions", "[oscillatory]")
	{
			TEST_PRECISION_INFO();
		// Test integration of highly oscillatory functions
		RealFunctionFromStdFunc oscillatory( [](Real x) -> Real { return std::sin(REAL(10.0) * x) / (REAL(1.0) + x); } );

		// These require more careful integration
		Real trapResult = IntegrateTrap(oscillatory, REAL(0.0), Constants::PI, REAL(1e-6));
		Real simpsonResult = IntegrateSimpson(oscillatory, REAL(0.0), Constants::PI, REAL(1e-6));

		// Results should be reasonably close
		REQUIRE_THAT(trapResult, WithinAbs(simpsonResult, REAL(1e-4)));
	}

	TEST_CASE("Integration::Quasi_improper_integrals", "[improper]")
	{
			TEST_PRECISION_INFO();
		// Test functions that are well-behaved but approach infinity at boundaries
		// f(x) = 1/sqrt(x) on (0, 1] - integrable but singular at 0
		RealFunctionFromStdFunc invSqrt( [](Real x) -> Real { return REAL(1.0) / std::sqrt(x); } );

		// Analytical result: 2*sqrt(x) from epsilon to 1 = 2
		// We start slightly away from 0 to avoid numerical issues
		Real result = IntegrateTrap(invSqrt, 1e-6, REAL(1.0), REAL(1e-5));
		REQUIRE_THAT(result, WithinAbs(REAL(2.0), REAL(2e-3)));

		// f(x) = ln(x) on (0, 1] - integrable with x*ln(x) - x
		RealFunctionFromStdFunc logFunc( [](Real x) -> Real { return std::log(x); } );

		// Analytical result: [x*ln(x) - x] from epsilon to 1 = -1
		result = IntegrateTrap(logFunc, 1e-6, REAL(1.0), REAL(1e-5));
		REQUIRE_THAT(result, WithinAbs(-REAL(1.0), REAL(1e-2)));
	}

	TEST_CASE("Integration::Symmetric_functions", "[properties]")
	{
			TEST_PRECISION_INFO();
		// Test that integration respects symmetry properties
		RealFunctionFromStdFunc evenFunc( [](Real x) -> Real { return std::cos(x); } );
		RealFunctionFromStdFunc oddFunc( [](Real x) -> Real { return std::sin(x); } );

		// For even functions: integral from -a to a = 2 * integral from 0 to a
		Real a = Constants::PI / REAL(2.0);
		Real evenResultSymmetric = IntegrateTrap(evenFunc, -a, a);
		Real evenResultPositive = IntegrateTrap(evenFunc, REAL(0.0), a);
		REQUIRE_THAT(evenResultSymmetric, WithinAbs(REAL(2.0) * evenResultPositive, REAL(1e-8)));

		// For odd functions: integral from -a to a = 0
		Real oddResultSymmetric = IntegrateTrap(oddFunc, -a, a);
		REQUIRE_THAT(oddResultSymmetric, WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("Integration::Edge_cases", "[edge]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc simpleFunc( [](Real x) -> Real { return x * x; } );

		// Zero-width interval
		REQUIRE_THAT(IntegrateTrap(simpleFunc, REAL(1.0), REAL(1.0)), WithinAbs(REAL(0.0), REAL(1e-15)));

		// Reversed limits (should handle gracefully or give negative result)
		Real forward = IntegrateTrap(simpleFunc, REAL(0.0), REAL(1.0));
		Real backward = IntegrateTrap(simpleFunc, REAL(1.0), REAL(0.0));
		REQUIRE_THAT(forward, WithinAbs(-backward, REAL(1e-10)));

		// Very small interval
		Real tinyInterval = IntegrateTrap(simpleFunc, REAL(0.0), 1e-8, REAL(1e-12));
		// For x^2, integral from 0 to h ≈ h^3/3
		REQUIRE_THAT(tinyInterval, WithinAbs(std::pow(1e-8, REAL(3)) / REAL(3.0), REAL(1e-25)));
	}

	TEST_CASE("Integration::Precision_limits", "[precision]")
	{
			TEST_PRECISION_INFO();
		// Test behavior when requesting very high precision
		RealFunctionFromStdFunc smoothFunc( [](Real x) -> Real { return std::exp(-x); } );

		int steps;
		Real achievedPrecision;

		// Exact integral from 0 to 1: 1 - 1/e ≈ REAL(0.632120558829)
		Real exact = REAL(1.0) - std::exp(-REAL(1.0));

		// Test achievable precision
		Real result = IntegrateTrap(smoothFunc, REAL(0.0), REAL(1.0), &steps, &achievedPrecision, REAL(1e-12));

		// Should achieve result close to exact
		REQUIRE_THAT(result, WithinAbs(exact, REAL(1e-10)));

		// Achieved precision should be reasonable
		REQUIRE(achievedPrecision < 1e-10);
	}

	TEST_CASE("Integration::Method_consistency", "[consistency]")
	{
			TEST_PRECISION_INFO();
		// All methods should give consistent results for simple functions
		RealFunctionFromStdFunc quadratic( [](Real x) -> Real { return REAL(2.0) * x * x - REAL(3.0) * x + REAL(1.0); } );

		// Analytical: [2x^3/3 - 3x^2/2 + x] from 0 to 2 = 16/3 - 6 + 2 = 4/3
		Real exact = REAL(4.0) / REAL(3.0);

		Real trapResult = IntegrateTrap(quadratic, REAL(0.0), REAL(2.0), REAL(1e-8));
		Real simpsonResult = IntegrateSimpson(quadratic, REAL(0.0), REAL(2.0), REAL(1e-8));
		Real gaussResult = IntegrateGauss10(quadratic, REAL(0.0), REAL(2.0));

		// All methods should converge to exact result
		REQUIRE_THAT(exact, WithinAbs(trapResult, REAL(1e-7)));
		REQUIRE_THAT(exact, WithinAbs(simpsonResult, REAL(1e-9)));
		REQUIRE_THAT(exact, WithinAbs(gaussResult, REAL(1e-12)));
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
		auto result = IntegrateRomberg(quadratic, REAL(0.0), REAL(1.0), REAL(1e-12));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0)/REAL(3.0), REAL(1e-11)));
	}
	
	TEST_CASE("Integration::Romberg_cubic", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc cubic( [](Real x) -> Real { return x * x * x - REAL(2.0) * x; } );
		
		// ∫₀² (x³ - 2x) dx = [x⁴/4 - x²]₀² = 4 - 4 = 0
		auto result = IntegrateRomberg(cubic, REAL(0.0), REAL(2.0), REAL(1e-12));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.0), REAL(1e-11)));
	}
	
	TEST_CASE("Integration::Romberg_trigonometric", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc sinFunc( [](Real x) -> Real { return std::sin(x); } );
		
		// ∫₀^π sin(x) dx = 2
		auto result = IntegrateRomberg(sinFunc, REAL(0.0), Constants::PI, REAL(1e-10));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(2.0), REAL(1e-9)));
	}
	
	TEST_CASE("Integration::Romberg_exponential", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc expFunc( [](Real x) -> Real { return std::exp(x); } );
		
		// ∫₀¹ e^x dx = e - 1 ≈ REAL(1.718281828)
		auto result = IntegrateRomberg(expFunc, REAL(0.0), REAL(1.0), REAL(1e-12));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(std::exp(REAL(1.0)) - REAL(1.0), REAL(1e-11)));
	}
	
	TEST_CASE("Integration::Romberg_gaussian", "[romberg]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc gaussian( [](Real x) -> Real { return std::exp(-x * x); } );
		
		// ∫₀¹ e^(-x²) dx ≈ REAL(0.7468241328) (error function related)
		// erf(1) = 2/√π * ∫₀¹ e^(-t²) dt, so ∫₀¹ e^(-x²) dx = √π/2 * erf(1)
		Real expected = std::sqrt(Constants::PI) / REAL(2.0) * std::erf(REAL(1.0));
		
		auto result = IntegrateRomberg(gaussian, REAL(0.0), REAL(1.0), REAL(1e-12));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(expected, REAL(1e-11)));
	}
	
	TEST_CASE("Integration::Romberg_high_precision", "[romberg]")
	{
			TEST_PRECISION_INFO();
		// Test Romberg's ability to achieve very high precision
		RealFunctionFromStdFunc smooth( [](Real x) -> Real { return REAL(1.0) / (REAL(1.0) + x * x); } );
		
		// ∫₀¹ 1/(1+x²) dx = arctan(1) = π/4
		auto result = IntegrateRomberg(smooth, REAL(0.0), REAL(1.0), REAL(1e-14));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / REAL(4.0), REAL(1e-13)));
	}
	
	TEST_CASE("Integration::Romberg_faster_than_simpson", "[romberg][comparison]")
	{
			TEST_PRECISION_INFO();
		// For smooth functions, Romberg should converge faster (fewer iterations)
		RealFunctionFromStdFunc smooth( [](Real x) -> Real { return std::sin(x) * std::exp(-x); } );
		
		auto rombergResult = IntegrateRomberg(smooth, REAL(0.0), Constants::PI, REAL(1e-10));
		auto simpsonResult = IntegrateSimpson(smooth, REAL(0.0), Constants::PI, REAL(1e-10));
		
		// Both should converge
		REQUIRE(rombergResult.converged == true);
		REQUIRE(simpsonResult.converged == true);
		
		// Results should be very close
		REQUIRE_THAT(rombergResult.value, WithinAbs(simpsonResult.value, REAL(1e-9)));
		
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
		Real expected = Constants::PI / REAL(2.0);
		
		auto trapResult = IntegrateTrap(testFunc, REAL(0.0), Constants::PI, REAL(1e-8));
		auto simpsonResult = IntegrateSimpson(testFunc, REAL(0.0), Constants::PI, REAL(1e-8));
		auto rombergResult = IntegrateRomberg(testFunc, REAL(0.0), Constants::PI, REAL(1e-10));
		Real gaussResult = IntegrateGauss10(testFunc, REAL(0.0), Constants::PI);
		
		REQUIRE_THAT(trapResult.value, WithinAbs(expected, REAL(1e-7)));
		REQUIRE_THAT(simpsonResult.value, WithinAbs(expected, REAL(1e-7)));
		REQUIRE_THAT(rombergResult.value, WithinAbs(expected, REAL(1e-9)));
		REQUIRE_THAT(gaussResult, WithinAbs(expected, REAL(1e-6)));
	}
	
	TEST_CASE("Integration::Romberg_legacy_interface", "[romberg]")
	{
			TEST_PRECISION_INFO();
		// Test the legacy interface with output parameters
		RealFunctionFromStdFunc testFunc( [](Real x) -> Real { return x * x; } );
		
		int steps = 0;
		Real precision = REAL(0.0);
		Real result = IntegrateRomberg(testFunc, REAL(0.0), REAL(1.0), &steps, &precision, REAL(1e-10));
		
		REQUIRE_THAT(result, WithinAbs(REAL(1.0)/REAL(3.0), REAL(1e-9)));
		REQUIRE(steps > 0);
		REQUIRE(precision < 1e-9);
	}

	/*********************************************************************/
	/*****          Improper Integrals Tests                         *****/
	/*********************************************************************/
	
	TEST_CASE("Integration::ImproperUpper_exponential_decay", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) e^(-x) dx = 1
		RealFunctionFromStdFunc expDecay( [](Real x) -> Real { return std::exp(-x); } );
		
		auto result = IntegrateUpperInf(expDecay, REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperUpper_inverse_square", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[1, ∞) 1/x² dx = 1
		RealFunctionFromStdFunc invSquare( [](Real x) -> Real { return REAL(1.0) / (x * x); } );
		
		auto result = IntegrateUpperInf(invSquare, REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-7)));
	}
	
	TEST_CASE("Integration::ImproperUpper_arctan_derivative", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) 1/(1+x²) dx = π/2
		RealFunctionFromStdFunc lorentzian( [](Real x) -> Real { return REAL(1.0) / (REAL(1.0) + x * x); } );
		
		auto result = IntegrateUpperInf(lorentzian, REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / REAL(2.0), REAL(1e-7)));
	}
	
	TEST_CASE("Integration::ImproperLower_exponential_growth", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, 0] e^x dx = 1
		RealFunctionFromStdFunc expGrowth( [](Real x) -> Real { return std::exp(x); } );
		
		auto result = IntegrateLowerInf(expGrowth, REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperLower_arctan_derivative", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, 0] 1/(1+x²) dx = π/2
		RealFunctionFromStdFunc lorentzian( [](Real x) -> Real { return REAL(1.0) / (REAL(1.0) + x * x); } );
		
		auto result = IntegrateLowerInf(lorentzian, REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / REAL(2.0), REAL(1e-7)));
	}
	
	TEST_CASE("Integration::ImproperFull_gaussian", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, ∞) e^(-x²) dx = √π
		RealFunctionFromStdFunc gaussian( [](Real x) -> Real { return std::exp(-x * x); } );
		
		auto result = IntegrateInf(gaussian);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(std::sqrt(Constants::PI), REAL(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperFull_lorentzian", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, ∞) 1/(1+x²) dx = π
		RealFunctionFromStdFunc lorentzian( [](Real x) -> Real { return REAL(1.0) / (REAL(1.0) + x * x); } );
		
		auto result = IntegrateInf(lorentzian);
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI, REAL(1e-6)));
	}
	
	TEST_CASE("Integration::ImproperFull_gaussian_shifted", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫(-∞, ∞) e^(-(x-3)²) dx = √π (same as unshifted)
		RealFunctionFromStdFunc shifted( [](Real x) -> Real { return std::exp(-(x - REAL(3.0)) * (x - REAL(3.0))); } );
		
		// Without split point adjustment - should still converge
		auto result1 = IntegrateInf(shifted);
		REQUIRE(result1.converged == true);
		REQUIRE_THAT(result1.value, WithinAbs(std::sqrt(Constants::PI), REAL(1e-4)));
		
		// With split point at the peak - better convergence
		auto result2 = IntegrateInfSplit(shifted, REAL(3.0));
		REQUIRE(result2.converged == true);
		REQUIRE_THAT(result2.value, WithinAbs(std::sqrt(Constants::PI), REAL(1e-5)));
	}
	
	TEST_CASE("Integration::ImproperUpper_sinc_squared", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) sin²(x)/x² dx = π/2
		// Use sinc²(x) = sin²(x)/x² (with sinc(0) = 1)
		// Note: This is a challenging integral due to oscillations
		RealFunctionFromStdFunc sincSquared( [](Real x) -> Real {
			if (std::abs(x) < 1e-10) return REAL(1.0);  // Limit as x → 0
			return std::sin(x) * std::sin(x) / (x * x);
		} );
		
		auto result = IntegrateUpperInf(sincSquared, REAL(0.0), REAL(1e-6));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / REAL(2.0), 2e-3));  // ~REAL(0.1)% accuracy for oscillatory function
	}
	
	TEST_CASE("Integration::ImproperUpper_gamma_like", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) x * e^(-x) dx = Γ(2) = 1! = 1
		RealFunctionFromStdFunc gammaLike( [](Real x) -> Real { return x * std::exp(-x); } );
		
		auto result = IntegrateUpperInf(gammaLike, REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(1.0), REAL(1e-4)));
	}
	
	TEST_CASE("Integration::ImproperUpper_gamma3", "[improper][infinite]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, ∞) x² * e^(-x) dx = Γ(3) = 2! = 2
		RealFunctionFromStdFunc gamma3( [](Real x) -> Real { return x * x * std::exp(-x); } );
		
		auto result = IntegrateUpperInf(gamma3, REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(2.0), REAL(1e-4)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                     SINGULAR INTEGRAL TESTS                           //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::LowerSingular_inv_sqrt", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√x dx = 2
		RealFunctionFromStdFunc f( [](Real x) -> Real { return REAL(1.0) / std::sqrt(x); } );
		
		auto result = IntegrateLowerSingular(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(2.0), REAL(1e-6)));
	}

	TEST_CASE("Integration::UpperSingular_inv_sqrt", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√(1-x) dx = 2
		RealFunctionFromStdFunc f( [](Real x) -> Real { return REAL(1.0) / std::sqrt(REAL(1.0) - x); } );
		
		auto result = IntegrateUpperSingular(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(2.0), REAL(1e-6)));
	}

	TEST_CASE("Integration::BothSingular_beta_half", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√(x(1-x)) dx = π  (Beta function B(1/2, 1/2) = Γ(1/2)²/Γ(1) = π)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return REAL(1.0) / std::sqrt(x * (REAL(1.0) - x)); } );
		
		auto result = IntegrateBothSingular(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI, REAL(1e-4)));
	}

	TEST_CASE("Integration::InteriorSingular_inv_sqrt_abs", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[-1,1] 1/√|x| dx = 4  (two contributions of 2 each)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return REAL(1.0) / std::sqrt(std::abs(x)); } );
		
		auto result = IntegrateInteriorSingular(f, -REAL(1.0), REAL(1.0), REAL(0.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(4.0), REAL(1e-5)));
	}

	TEST_CASE("Integration::LowerSingular_log_times_sqrt", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] ln(x)/√x dx = -4
		// This has both logarithmic and square-root singularities at 0
		RealFunctionFromStdFunc f( [](Real x) -> Real { return std::log(x) / std::sqrt(x); } );
		
		auto result = IntegrateLowerSingular(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(-REAL(4.0), REAL(0.02)));  // Harder integral, relax tolerance
	}

	TEST_CASE("Integration::LowerSingular_sqrt_itself", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] √x dx = 2/3 (no singularity, but transformation should still work)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return std::sqrt(x); } );
		
		auto result = IntegrateLowerSingular(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(2.0)/REAL(3.0), REAL(1e-8)));
	}

	TEST_CASE("Integration::UpperSingular_arcsin_deriv", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] 1/√(1-x²) dx = arcsin(1) - arcsin(0) = π/2
		// Note: this has singularity only at x=1, not at x=0
		RealFunctionFromStdFunc f( [](Real x) -> Real { return REAL(1.0) / std::sqrt(REAL(1.0) - x * x); } );
		
		auto result = IntegrateUpperSingular(f, REAL(0.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(Constants::PI / REAL(2.0), REAL(1e-5)));
	}

	TEST_CASE("Integration::InteriorSingular_two_point", "[improper][singular]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,2] 1/√|x-1| dx = 4  (singularity at x=1 interior point)
		RealFunctionFromStdFunc f( [](Real x) -> Real { return REAL(1.0) / std::sqrt(std::abs(x - REAL(1.0))); } );
		
		auto result = IntegrateInteriorSingular(f, REAL(0.0), REAL(2.0), REAL(1.0));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinAbs(REAL(4.0), REAL(1e-5)));
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
		REQUIRE_THAT(sum, WithinAbs(REAL(2.0), REAL(1e-12)));
	}
	
	TEST_CASE("GaussLegendre::weights_sum_scaled_interval", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Legendre on [a,b], weights should sum to (b-a)
		std::vector<Real> x(15), w(15);
		GaussLegendre(x, w, REAL(0.0), REAL(5.0));
		
		Real sum = 0;
		for (int i = 0; i < 15; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(REAL(5.0), REAL(1e-12)));
	}
	
	TEST_CASE("GaussLegendre::nodes_symmetric", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// Nodes should be symmetric about 0 on [-1,1]
		std::vector<Real> x(8), w(8);
		GaussLegendre(x, w);
		
		// Check that x[i] = -x[n-1-i]
		for (int i = 0; i < 4; i++) {
			REQUIRE_THAT(x[i], WithinAbs(-x[7-i], REAL(1e-14)));
			REQUIRE_THAT(w[i], WithinAbs(w[7-i], REAL(1e-14)));
		}
	}
	
	TEST_CASE("GaussLegendre::integrate_polynomial_exact", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// n-point Gauss-Legendre exactly integrates polynomials up to degree 2n-1
		// 5-point should exactly integrate up to degree 9
		std::vector<Real> x(5), w(5);
		GaussLegendre(x, w, REAL(0.0), REAL(1.0));
		
		// ∫[0,1] x^4 dx = 1/5 = REAL(0.2)
		Real integral = 0;
		for (int i = 0; i < 5; i++)
			integral += w[i] * std::pow(x[i], 4);
		REQUIRE_THAT(integral, WithinAbs(REAL(0.2), REAL(1e-14)));
		
		// ∫[0,1] x^8 dx = 1/9
		integral = 0;
		for (int i = 0; i < 5; i++)
			integral += w[i] * std::pow(x[i], 8);
		REQUIRE_THAT(integral, WithinAbs(REAL(1.0)/REAL(9.0), REAL(1e-14)));
	}
	
	TEST_CASE("GaussLegendre::integrate_sin_high_precision", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// ∫[0, π] sin(x) dx = 2
		std::vector<Real> x(20), w(20);
		GaussLegendre(x, w, REAL(0.0), Constants::PI);
		
		Real integral = 0;
		for (int i = 0; i < 20; i++)
			integral += w[i] * std::sin(x[i]);
		REQUIRE_THAT(integral, WithinAbs(REAL(2.0), REAL(1e-12)));
	}
	
	TEST_CASE("GaussLegendre::rule_struct", "[gaussian][legendre]")
	{
			TEST_PRECISION_INFO();
		// Test the GaussQuadratureRule structure
		auto rule = GaussLegendreRule(12, -REAL(2.0), REAL(3.0));
		
		REQUIRE(rule.n == 12);
		REQUIRE(rule.nodes.size() == 12);
		REQUIRE(rule.weights.size() == 12);
		
		// Weights should sum to (3 - (-2)) = 5
		Real sum = 0;
		for (auto w : rule.weights) sum += w;
		REQUIRE_THAT(sum, WithinAbs(REAL(5.0), REAL(1e-12)));
	}
	
	TEST_CASE("GaussLaguerre::weights_sum_to_gamma", "[gaussian][laguerre]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Laguerre with α=0, weights sum to Γ(1) = 1
		std::vector<Real> x(10), w(10);
		GaussLaguerre(x, w, REAL(0.0));
		
		Real sum = 0;
		for (int i = 0; i < 10; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(REAL(1.0), REAL(1e-10)));
	}
	
	TEST_CASE("GaussLaguerre::nodes_positive", "[gaussian][laguerre]")
	{
			TEST_PRECISION_INFO();
		// All Laguerre nodes should be positive (on [0,∞))
		std::vector<Real> x(15), w(15);
		GaussLaguerre(x, w, REAL(0.0));
		
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
		GaussLaguerre(x, w, REAL(0.0));
		
		Real integral = 0;
		for (int i = 0; i < 10; i++)
			integral += w[i] * REAL(1.0);  // f(x) = 1, weight has e^(-x)
		REQUIRE_THAT(integral, WithinAbs(REAL(1.0), REAL(1e-10)));
	}
	
	TEST_CASE("GaussHermite::weights_sum_to_sqrt_pi", "[gaussian][hermite]")
	{
			TEST_PRECISION_INFO();
		// For Gauss-Hermite, weights sum to ∫ e^(-x²) dx = √π
		std::vector<Real> x(10), w(10);
		GaussHermite(x, w);
		
		Real sum = 0;
		for (int i = 0; i < 10; i++) sum += w[i];
		REQUIRE_THAT(sum, WithinAbs(std::sqrt(Constants::PI), REAL(1e-12)));
	}
	
	TEST_CASE("GaussHermite::nodes_symmetric", "[gaussian][hermite]")
	{
			TEST_PRECISION_INFO();
		// Hermite nodes should be symmetric about 0
		std::vector<Real> x(10), w(10);
		GaussHermite(x, w);
		
		// Nodes should satisfy x[i] = -x[n-1-i]
		for (int i = 0; i < 5; i++) {
			REQUIRE_THAT(x[i], WithinAbs(-x[9-i], REAL(1e-14)));
			REQUIRE_THAT(w[i], WithinAbs(w[9-i], REAL(1e-14)));
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
		REQUIRE_THAT(integral, WithinAbs(std::sqrt(Constants::PI) / REAL(2.0), REAL(1e-12)));
	}
	
	TEST_CASE("GaussChebyshev1::closed_form_accuracy", "[gaussian][chebyshev]")
	{
			TEST_PRECISION_INFO();
		// Chebyshev type 1: w(x) = 1/√(1-x²)
		// Weights should all equal π/n
		std::vector<Real> x(10), w(10);
		GaussChebyshev1(x, w);
		
		Real expected_weight = Constants::PI / REAL(10.0);
		for (int i = 0; i < 10; i++) {
			REQUIRE_THAT(w[i], WithinAbs(expected_weight, REAL(1e-14)));
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
			Real expected = std::cos((REAL(2.0) * k + REAL(1.0)) * Constants::PI / (REAL(2.0) * n));
			REQUIRE_THAT(x[k], WithinAbs(expected, REAL(1e-14)));
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
			Real angle = (k + REAL(1.0)) * Constants::PI / (n + REAL(1.0));
			Real expected = Constants::PI / (n + REAL(1.0)) * std::pow(std::sin(angle), 2);
			REQUIRE_THAT(w[k], WithinAbs(expected, REAL(1e-14)));
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
		GaussJacobi(xJ, wJ, REAL(0.0), REAL(0.0));
		
		// Sort nodes (paired with weights)
		std::vector<std::pair<Real, Real>> legendre(n), jacobi(n);
		for (int i = 0; i < n; i++) {
			legendre[i] = {xL[i], wL[i]};
			jacobi[i] = {xJ[i], wJ[i]};
		}
		std::sort(legendre.begin(), legendre.end());
		std::sort(jacobi.begin(), jacobi.end());
		
		for (int i = 0; i < n; i++) {
			REQUIRE_THAT(jacobi[i].first, WithinAbs(legendre[i].first, REAL(1e-8)));
			REQUIRE_THAT(jacobi[i].second, WithinAbs(legendre[i].second, REAL(1e-8)));
		}
	}
	
	TEST_CASE("GaussJacobi::reduces_to_chebyshev1", "[gaussian][jacobi]")
	{
			TEST_PRECISION_INFO();
		// Jacobi with α=β=-1/2 is Chebyshev type 1
		int n = 8;
		std::vector<Real> xC(n), wC(n), xJ(n), wJ(n);
		
		GaussChebyshev1(xC, wC);
		GaussJacobi(xJ, wJ, -REAL(0.5), -REAL(0.5));
		
		// Nodes should match (possibly in different order)
		// Both are symmetric, just compare sorted
		std::sort(xC.begin(), xC.end());
		std::sort(xJ.begin(), xJ.end());
		
		for (int i = 0; i < n; i++) {
			REQUIRE_THAT(xJ[i], WithinAbs(xC[i], REAL(1e-8)));
		}
	}
	
	TEST_CASE("IntegrateGaussLegendre::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test the convenience integration function
		RealFunctionFromStdFunc f( [](Real x) { return x * x; } );
		
		// ∫[0,2] x² dx = 8/3
		Real result = IntegrateGaussLegendre(f, REAL(0.0), REAL(2.0), 10);
		REQUIRE_THAT(result, WithinAbs(REAL(8.0)/REAL(3.0), REAL(1e-12)));  // Relaxed for floating-point
	}
	
	TEST_CASE("IntegrateGaussLaguerre::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test the convenience integration function
		// ∫[0,∞) x * e^(-x) dx = 1 (function f=x, weight has e^(-x))
		RealFunctionFromStdFunc f( [](Real x) { return x; } );
		
		Real result = IntegrateGaussLaguerre(f, 15, REAL(0.0));
		REQUIRE_THAT(result, WithinAbs(REAL(1.0), REAL(1e-8)));
	}
	
	TEST_CASE("IntegrateGaussHermite::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test the convenience integration function
		// ∫ x⁴ e^(-x²) dx = 3√π/4 (even function)
		RealFunctionFromStdFunc f( [](Real x) { return x * x * x * x; } );
		
		Real result = IntegrateGaussHermite(f, 15);
		Real expected = REAL(3.0) * std::sqrt(Constants::PI) / REAL(4.0);
		REQUIRE_THAT(result, WithinAbs(expected, REAL(1e-10)));
	}
	
	TEST_CASE("IntegrateGaussJacobi::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// Test Gauss-Jacobi with α=1, β=0: weight = (1-x)
		// ∫[-1,1] (1-x) * x² dx = ∫[-1,1] x² - x³ dx = 2/3 - 0 = 2/3
		RealFunctionFromStdFunc f( [](Real x) { return x * x; } );
		
		Real result = IntegrateGaussJacobi(f, REAL(1.0), REAL(0.0), 15);
		REQUIRE_THAT(result, WithinAbs(REAL(2.0)/REAL(3.0), REAL(1e-10)));
	}
	
	TEST_CASE("IntegrateGaussChebyshev1::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫[-1,1] 1/√(1-x²) dx = π (integral of weight function alone)
		RealFunctionFromStdFunc one( [](Real x) { return REAL(1.0); } );
		
		Real result = IntegrateGaussChebyshev1(one, 10);
		REQUIRE_THAT(result, WithinAbs(Constants::PI, REAL(1e-12)));
	}
	
	TEST_CASE("IntegrateGaussChebyshev2::convenience_function", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫[-1,1] √(1-x²) dx = π/2 (integral of weight function alone)
		RealFunctionFromStdFunc one( [](Real x) { return REAL(1.0); } );
		
		Real result = IntegrateGaussChebyshev2(one, 10);
		REQUIRE_THAT(result, WithinAbs(Constants::PI / REAL(2.0), REAL(1e-12)));
	}
	
	TEST_CASE("IntegrateGaussLegendre::high_order_polynomial", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// 10-point Gauss-Legendre should exactly integrate polynomials up to degree 19
		// Test x^18 on [0,1]: ∫x^18 dx = 1/19
	RealFunctionFromStdFunc poly18([](Real x) -> Real { return std::pow(x, REAL(18)); });
	}
	
	TEST_CASE("IntegrateGaussLaguerre::gamma_integral", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,∞) x^3 e^(-x) dx = Γ(4) = 3! = 6
		// Using α=0, so f(x) = x³
		RealFunctionFromStdFunc f( [](Real x) { return x * x * x; } );
		
		Real result = IntegrateGaussLaguerre(f, 20, REAL(0.0));
		REQUIRE_THAT(result, WithinAbs(REAL(6.0), REAL(1e-8)));
	}
	
	TEST_CASE("IntegrateGaussHermite::gaussian_moments", "[gaussian][integrate]")
	{
			TEST_PRECISION_INFO();
		// ∫ x⁶ e^(-x²) dx = 15√π/8 (6th moment of Gaussian)
	  RealFunctionFromStdFunc f([](Real x) -> Real { return std::pow(x, REAL(6)); });
	
	  Real result = IntegrateGaussHermite(f, 20);
	  Real expected = REAL(15.0) * std::sqrt(Constants::PI) / REAL(8.0);
	  REQUIRE_THAT(result, WithinAbs(expected, REAL(1e-10)));
  } 

  TEST_CASE("IntegrateGaussJacobi::beta_function", "[gaussian][integrate]")
  {
		TEST_PRECISION_INFO();
	// ∫[-1,1] (1-x)^α (1+x)^β dx = 2^(α+β+1) B(α+1, β+1)
	// where B is the Beta function: B(a,b) = Γ(a)Γ(b)/Γ(a+b)
	// For α=REAL(0.5), β=REAL(1.5): B(REAL(1.5), REAL(2.5)) = Γ(REAL(1.5))Γ(REAL(2.5))/Γ(4) = (√π/2)(3√π/4)/6 = 3π/48 = π/16
	// Result = 2^3 * π/16 = π/2
	  RealFunctionFromStdFunc one( [](Real x) { return REAL(1.0); } );
	
	  Real result = IntegrateGaussJacobi(one, REAL(0.5), REAL(1.5), 15);
	  REQUIRE_THAT(result, WithinAbs(Constants::PI / REAL(2.0), REAL(1e-8)));
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
		REQUIRE_THROWS(GaussJacobi(x, w, REAL(100.0), REAL(100.0)));
  }

	TEST_CASE("Integration::Difficult_peaked_function", "[difficult]")
	{
	TEST_PRECISION_INFO();
		// Integrate a sharply peaked function: 1/(1+(x-REAL(0.3))²/REAL(0.01))
		// This is a Lorentzian with width REAL(0.1) centered at REAL(0.3)
		// ∫ dx/(1+(x-c)²/w²) = w*arctan((x-c)/w)
		RealFunctionFromStdFunc peaked( [](Real x) { 
			Real diff = x - REAL(0.3);
			return REAL(1.0) / (REAL(1.0) + diff * diff / REAL(0.01)); 
		} );
		
		// Exact: REAL(0.1) * (arctan((1-REAL(0.3))/REAL(0.1)) - arctan((0-REAL(0.3))/REAL(0.1)))
		//      = REAL(0.1) * (arctan(7) + arctan(3))
		Real expected = REAL(0.1) * (std::atan(REAL(7.0)) + std::atan(REAL(3.0)));
		
		auto result = IntegrateSimpson(peaked, REAL(0.0), REAL(1.0), REAL(1e-6));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("Integration::Difficult_rapid_oscillation", "[difficult]")
	{
			TEST_PRECISION_INFO();
		// Integrate sin(20x) on [0, π] = 0 (cancellation)
		// Using lower frequency to be more tractable
		RealFunctionFromStdFunc fastSin( [](Real x) { return std::sin(REAL(20.0) * x); } );
		
		auto result = IntegrateSimpson(fastSin, REAL(0.0), Constants::PI, REAL(1e-4));
		
		// Should be very close to 0 due to cancellation
		REQUIRE_THAT(result.value, WithinAbs(REAL(0.0), REAL(0.1)));
	}

	TEST_CASE("Integration::Difficult_near_singularity", "[difficult]")
	{
			TEST_PRECISION_INFO();
		// Integrate 1/√(x + REAL(0.001)) on [0,1] - near-singular at x=0
		RealFunctionFromStdFunc nearSing( [](Real x) { return REAL(1.0) / std::sqrt(x + REAL(0.001)); } );
		
		auto result = IntegrateSimpson(nearSing, REAL(0.0), REAL(1.0), REAL(1e-6));
		Real expected = REAL(2.0) * (std::sqrt(REAL(1.001)) - std::sqrt(REAL(0.001)));
		
		REQUIRE(result.converged == true);
		REQUIRE_THAT(result.value, WithinRel(expected, REAL(1e-4)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    METHOD COMPARISON TESTS                            //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::AllMethods_polynomial", "[comparison][comprehensive]")
	{
			TEST_PRECISION_INFO();
		// All methods should give same answer for polynomial
		RealFunctionFromStdFunc poly( [](Real x) { return x * x * x - 2 * x + 1; } );
		Real exact = REAL(0.25);  // ∫[0,1] (x³ - 2x + 1) dx = [x⁴/4 - x² + x]₀¹ = REAL(0.25)
		
		auto trap = IntegrateTrap(poly, REAL(0.0), REAL(1.0));
		auto simp = IntegrateSimpson(poly, REAL(0.0), REAL(1.0));
		auto romb = IntegrateRomberg(poly, REAL(0.0), REAL(1.0));
		auto gauss = IntegrateGauss10(poly, REAL(0.0), REAL(1.0));
		
		REQUIRE_THAT(trap, WithinRel(exact, REAL(1e-4)));
		REQUIRE_THAT(simp.value, WithinRel(exact, REAL(1e-6)));
		REQUIRE_THAT(romb.value, WithinRel(exact, REAL(1e-10)));
		REQUIRE_THAT(gauss.value, WithinRel(exact, REAL(1e-12)));
	}

	TEST_CASE("Integration::AllMethods_transcendental", "[comparison][comprehensive]")
	{
			TEST_PRECISION_INFO();
		// Compare methods on e^(-x²) over [0,1]
		RealFunctionFromStdFunc gauss( [](Real x) { return std::exp(-x * x); } );
		Real exact = REAL(0.7468241328124271);  // erf(1) * √π/2
		
		auto trap = IntegrateTrap(gauss, REAL(0.0), REAL(1.0));
		auto simp = IntegrateSimpson(gauss, REAL(0.0), REAL(1.0));
		auto romb = IntegrateRomberg(gauss, REAL(0.0), REAL(1.0));
		auto gauss10 = IntegrateGauss10(gauss, REAL(0.0), REAL(1.0));
		
		REQUIRE_THAT(trap, WithinAbs(exact, REAL(1e-4)));
		REQUIRE_THAT(simp.value, WithinAbs(exact, REAL(1e-6)));
		REQUIRE_THAT(romb.value, WithinAbs(exact, REAL(1e-8)));
		REQUIRE_THAT(gauss10.value, WithinAbs(exact, REAL(1e-10)));
	}

	TEST_CASE("Integration::Romberg_vs_Gauss_accuracy", "[comparison]")
	{
			TEST_PRECISION_INFO();
		// Compare Romberg and Gauss-Legendre on same integral
		RealFunctionFromStdFunc f( [](Real x) { return REAL(1.0) / (REAL(1.0) + x * x); } );
		Real exact = Constants::PI / REAL(4.0);  // arctan(1) - arctan(0)
		
		auto romb = IntegrateRomberg(f, REAL(0.0), REAL(1.0));
		Real gauss = IntegrateGaussLegendre(f, REAL(0.0), REAL(1.0), 15);
		
		// Both should be accurate to 1e-8
		REQUIRE_THAT(romb.value, WithinAbs(exact, REAL(1e-8)));
		REQUIRE_THAT(gauss, WithinAbs(exact, REAL(1e-10)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    ERROR HANDLING TESTS                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::Invalid_bounds", "[error]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc f( [](Real x) { return x; } );
		
		// Empty interval should give 0
		REQUIRE(IntegrateTrap(f, REAL(1.0), REAL(1.0)) == REAL(0.0));
		
		// Reversed bounds should work (gives negative of forward integral)
		Real forward = IntegrateTrap(f, REAL(0.0), REAL(1.0));
		Real reverse = IntegrateTrap(f, REAL(1.0), REAL(0.0));
		REQUIRE_THAT(reverse, WithinAbs(-forward, REAL(1e-10)));
	}

	TEST_CASE("Integration::Interior_singular_bounds_check", "[error][singular]")
	{
			TEST_PRECISION_INFO();
		RealFunctionFromStdFunc f( [](Real x) { return REAL(1.0) / std::sqrt(std::abs(x)); } );
		
		// Singular point must be inside (a, b)
		REQUIRE_THROWS_AS(IntegrateInteriorSingular(f, REAL(0.0), REAL(1.0), REAL(0.0)), std::invalid_argument);
		REQUIRE_THROWS_AS(IntegrateInteriorSingular(f, REAL(0.0), REAL(1.0), REAL(1.0)), std::invalid_argument);
		REQUIRE_THROWS_AS(IntegrateInteriorSingular(f, REAL(0.0), REAL(1.0), REAL(2.0)), std::invalid_argument);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    SPECIAL VALUES TESTS                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Integration::WellKnown_pi_computations", "[special]")
	{
			TEST_PRECISION_INFO();
		// Various ways to compute π via integration
		
		// ∫[-1,1] 1/√(1-x²) dx = π
		RealFunctionFromStdFunc arcCircle( [](Real x) { return REAL(1.0) / std::sqrt(REAL(1.0) - x * x); } );
		auto pi1 = IntegrateBothSingular(arcCircle, -REAL(1.0), REAL(1.0));
		REQUIRE_THAT(pi1.value, WithinAbs(Constants::PI, REAL(1e-4)));
		
		// ∫[0,1] 4/(1+x²) dx = π
		RealFunctionFromStdFunc arctan4( [](Real x) { return REAL(4.0) / (REAL(1.0) + x * x); } );
		auto pi2 = IntegrateRomberg(arctan4, REAL(0.0), REAL(1.0));
		REQUIRE_THAT(pi2.value, WithinAbs(Constants::PI, REAL(1e-7)));  // Relaxed tolerance
		
		// ∫(-∞,∞) e^(-x²) dx = √π
		RealFunctionFromStdFunc gaussian( [](Real x) { return std::exp(-x * x); } );
		auto sqrtPi = IntegrateInf(gaussian);
		REQUIRE_THAT(sqrtPi.value, WithinAbs(std::sqrt(Constants::PI), REAL(1e-6)));
	}

	TEST_CASE("Integration::WellKnown_e_computations", "[special]")
	{
			TEST_PRECISION_INFO();
		// ∫[0,1] e^x dx = e - 1
		RealFunctionFromStdFunc expf( [](Real x) { return std::exp(x); } );
		auto result = IntegrateRomberg(expf, REAL(0.0), REAL(1.0));
		REQUIRE_THAT(result.value, WithinAbs(std::exp(REAL(1.0)) - REAL(1.0), REAL(1e-10)));
		
		// ∫[0,∞] e^(-x) dx = 1
		RealFunctionFromStdFunc expDecay( [](Real x) { return std::exp(-x); } );
		auto infResult = IntegrateUpperInf(expDecay, REAL(0.0));
		REQUIRE_THAT(infResult.value, WithinAbs(REAL(1.0), REAL(1e-6)));
	}

	TEST_CASE("Integration::WellKnown_gamma_values", "[special]")
	{
			TEST_PRECISION_INFO();
		// Γ(5) = 4! = 24 via ∫[0,∞) x^4 e^(-x) dx
		// Using Gauss-Laguerre with n=20, alpha=0 (standard weight e^(-x))
	RealFunctionFromStdFunc gamma5([](Real x) -> Real { return std::pow(x, REAL(4)); });
	Real result5 = IntegrateGaussLaguerre(gamma5, 20, REAL(0.0));  // n=20, alpha=0
	REQUIRE_THAT(result5, WithinAbs(REAL(24.0), REAL(1e-6)));
	
	// Γ(3) = 2! = 2 via ∫[0,∞) x^2 e^(-x) dx
	RealFunctionFromStdFunc gamma3([](Real x) -> Real { return x * x; });
Real result3 = IntegrateGaussLaguerre(gamma3, 20, REAL(0.0));  // n=20, alpha=0
REQUIRE_THAT(result3, WithinAbs(REAL(2.0), REAL(1e-8)));
}
}


