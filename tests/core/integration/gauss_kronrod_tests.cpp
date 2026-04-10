///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        gauss_kronrod_tests.cpp                                             ///
///  Description: Unit tests for Gauss-Kronrod adaptive quadrature                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"
#include <catch2/catch_approx.hpp>

#include "MMLBase.h"
#include "core/Integration/GaussKronrod.h"

#include <cmath>

using namespace MML;
using namespace MML::Integration;
using namespace Catch;

///////////////////////////////////////////////////////////////////////////////////////////
///                        BASIC G7K15 RULE TESTS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateGK15 - basic polynomial integrals", "[gauss-kronrod][integration]")
{
	SECTION("Constant function")
	{
		// ∫₀¹ 5 dx = 5
		auto f = [](Real x) { return 5.0; };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(5.0).epsilon(TOL(1e-14, 1e-5)));
		REQUIRE(result.error_estimate < TOL(1e-14, 1e-5));
	}

	SECTION("Linear function")
	{
		// ∫₀¹ x dx = 0.5
		auto f = [](Real x) { return x; };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(0.5).epsilon(TOL(1e-14, 1e-5)));
	}

	SECTION("Quadratic function")
	{
		// ∫₀¹ x² dx = 1/3
		auto f = [](Real x) { return x * x; };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(1.0 / 3.0).epsilon(TOL(1e-14, 1e-5)));
	}

	SECTION("High-degree polynomial (exact for degree ≤ 29)")
	{
		// ∫₀¹ x^10 dx = 1/11
		auto f = [](Real x) { return std::pow(x, 10); };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(1.0 / 11.0).epsilon(TOL(1e-13, 1e-5)));
	}
}

TEST_CASE("IntegrateGK15 - transcendental functions", "[gauss-kronrod][integration]")
{
	SECTION("sin(x) over [0, π]")
	{
		// ∫₀^π sin(x) dx = 2
		auto f = [](Real x) { return std::sin(x); };
		auto result = IntegrateGK15(f, 0.0, Constants::PI);
		REQUIRE(result.value == Approx(2.0).epsilon(TOL(1e-14, 1e-5)));
	}

	SECTION("exp(x) over [0, 1]")
	{
		// ∫₀¹ e^x dx = e - 1
		auto f = [](Real x) { return std::exp(x); };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(std::exp(1.0) - 1.0).epsilon(TOL(1e-14, 1e-5)));
	}

	SECTION("1/(1+x²) over [0, 1]")
	{
		// ∫₀¹ 1/(1+x²) dx = arctan(1) = π/4
		auto f = [](Real x) { return 1.0 / (1.0 + x * x); };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(Constants::PI / 4.0).epsilon(TOL(1e-14, 1e-5)));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        G10K21 AND G15K31 RULE TESTS                                ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateGK21 - higher accuracy tests", "[gauss-kronrod][integration]")
{
	SECTION("sin(x) over [0, π]")
	{
		auto f = [](Real x) { return std::sin(x); };
		auto result = IntegrateGK21(f, 0.0, Constants::PI);
		REQUIRE(result.value == Approx(2.0).epsilon(TOL(1e-15, 1e-5)));
	}

	SECTION("Gaussian bell curve")
	{
		// ∫₋₃³ exp(-x²) dx ≈ √π * erf(3) ≈ 1.7724538509055159
		auto f = [](Real x) { return std::exp(-x * x); };
		auto result = IntegrateGK21(f, -3.0, 3.0);
		Real expected = std::sqrt(Constants::PI) * std::erf(3.0);
		// Note: TOL(1e-10, 1e-5) relative tolerance due to slight numerical differences
		// in erf() implementations across platforms
		REQUIRE(result.value == Approx(expected).epsilon(TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("IntegrateGK31 - highest accuracy tests", "[gauss-kronrod][integration]")
{
	SECTION("Oscillatory function")
	{
		// ∫₀^π cos(x)² dx = π/2
		auto f = [](Real x) { return std::cos(x) * std::cos(x); };
		auto result = IntegrateGK31(f, 0.0, Constants::PI);
		REQUIRE(result.value == Approx(Constants::PI / 2.0).epsilon(TOL(1e-14, 1e-5)));
	}

	SECTION("High-degree polynomial (exact for degree ≤ 61)")
	{
		// ∫₀¹ x^20 dx = 1/21
		auto f = [](Real x) { return std::pow(x, 20); };
		auto result = IntegrateGK31(f, 0.0, 1.0);
		REQUIRE(result.value == Approx(1.0 / 21.0).epsilon(TOL(1e-14, 1e-5)));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        ADAPTIVE INTEGRATION TESTS                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateGKAdaptive - functions requiring subdivision", "[gauss-kronrod][integration][adaptive]")
{
	SECTION("Peaked function")
	{
		// ∫₀¹ 1/(1+(x-0.5)²/0.001) dx - sharp peak at x=0.5
		auto f = [](Real x) { 
			Real t = x - 0.5;
			return 1.0 / (1.0 + t * t / 0.001); 
		};
		auto result = IntegrateGKAdaptive(f, 0.0, 1.0, TOL(1e-8, 1e-4), TOL(1e-8, 1e-4), 50, GKRule::GK15);
		
		// Exact: √0.001 * (arctan(0.5/√0.001) + arctan(0.5/√0.001))
		Real sqrtEps = std::sqrt(0.001);
		Real expected = sqrtEps * (std::atan(0.5 / sqrtEps) + std::atan(0.5 / sqrtEps));
		
		REQUIRE(result.value == Approx(expected).epsilon(TOL(1e-6, 1e-4)));
		REQUIRE(result.converged);
	}

	SECTION("Oscillatory over long interval")
	{
		// ∫₀^10π sin(x) dx = 0 (by symmetry)
		auto f = [](Real x) { return std::sin(x); };
		auto result = IntegrateGKAdaptive(f, 0.0, 10 * Constants::PI, TOL(1e-8, 1e-4), TOL(1e-8, 1e-4));
		REQUIRE(result.value == Approx(0.0).margin(TOL(1e-6, 1e-4)));
	}

	SECTION("Step-like function")
	{
		// ∫₀¹ tanh(100*(x-0.5)) dx ≈ 0 (antisymmetric about 0.5)
		auto f = [](Real x) { return std::tanh(100.0 * (x - 0.5)); };
		auto result = IntegrateGKAdaptive(f, 0.0, 1.0, 1e-6, 1e-6);
		REQUIRE(result.value == Approx(0.0).margin(1e-4));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        ERROR ESTIMATION TESTS                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Gauss-Kronrod error estimates", "[gauss-kronrod][integration]")
{
	SECTION("Error estimate is meaningful for smooth function")
	{
		auto f = [](Real x) { return std::exp(x); };
		auto result = IntegrateGK15(f, 0.0, 1.0);
		
		Real true_value = std::exp(1.0) - 1.0;
		Real actual_error = std::abs(result.value - true_value);
		
		// Error estimate should bound actual error (conservatively)
		REQUIRE(actual_error < result.error_estimate * 10 + TOL(1e-15, 1e-5));
	}

	SECTION("Higher-order rules give better error")
	{
		// A function where higher-order is beneficial
		auto f = [](Real x) { return std::sin(5.0 * x); };
		
		auto result15 = IntegrateGK15(f, 0.0, Constants::PI);
		auto result21 = IntegrateGK21(f, 0.0, Constants::PI);
		auto result31 = IntegrateGK31(f, 0.0, Constants::PI);
		
		// True value: ∫sin(5x)dx = -cos(5x)/5, from 0 to π: (-cos(5π)+cos(0))/5 = 2/5
		Real true_value = 2.0 / 5.0;
		
		Real err15 = std::abs(result15.value - true_value);
		Real err21 = std::abs(result21.value - true_value);
		Real err31 = std::abs(result31.value - true_value);
		
		// All should be reasonably accurate for this smooth function
		REQUIRE(err15 < TOL(1e-10, 1e-5));
		REQUIRE(err21 < TOL(1e-12, 1e-5));
		REQUIRE(err31 < TOL(1e-14, 1e-5));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        CONVENIENCE FUNCTION TESTS                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Gauss-Kronrod convenience functions", "[gauss-kronrod][integration]")
{
	auto f = [](Real x) { return std::sin(x); };
	Real expected = 2.0;  // ∫₀^π sin(x) dx = 2
	
	SECTION("IntegrateGaussKronrod15")
	{
		Real result = IntegrateGaussKronrod15(f, 0.0, Constants::PI);
		REQUIRE(result == Approx(expected).epsilon(TOL(1e-14, 1e-5)));
	}
	
	SECTION("IntegrateGaussKronrod21")
	{
		Real result = IntegrateGaussKronrod21(f, 0.0, Constants::PI);
		REQUIRE(result == Approx(expected).epsilon(TOL(1e-14, 1e-5)));
	}
	
	SECTION("IntegrateGaussKronrod31")
	{
		Real result = IntegrateGaussKronrod31(f, 0.0, Constants::PI);
		REQUIRE(result == Approx(expected).epsilon(TOL(1e-14, 1e-5)));
	}
	
	SECTION("IntegrateGaussKronrod (adaptive)")
	{
		Real result = IntegrateGaussKronrod(f, 0.0, Constants::PI);
		REQUIRE(result == Approx(expected).epsilon(TOL(1e-10, 1e-5)));
	}
	
	SECTION("IntegrateGaussKronrod with error output")
	{
		Real error = 0.0;
		Real result = IntegrateGaussKronrod(f, 0.0, Constants::PI, &error);
		REQUIRE(result == Approx(expected).epsilon(TOL(1e-10, 1e-5)));
		REQUIRE(error < TOL(1e-10, 1e-5));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        COMPARISON WITH OTHER METHODS                                ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Gauss-Kronrod vs other integration methods", "[gauss-kronrod][integration][comparison]")
{
	SECTION("Comparison with exact polynomial integration")
	{
		// G7K15 should be exact for polynomials up to degree 29
		// G10K21 should be exact for polynomials up to degree 41
		// G15K31 should be exact for polynomials up to degree 61
		
		// Test with degree 25 polynomial (should be exact for all rules)
		auto f = [](Real x) { return std::pow(x, 25); };
		Real expected = 1.0 / 26.0;  // ∫₀¹ x^25 dx = 1/26
		
		auto result15 = IntegrateGK15(f, 0.0, 1.0);
		auto result21 = IntegrateGK21(f, 0.0, 1.0);
		auto result31 = IntegrateGK31(f, 0.0, 1.0);
		
		REQUIRE(result15.value == Approx(expected).epsilon(TOL(1e-12, 1e-5)));
		REQUIRE(result21.value == Approx(expected).epsilon(TOL(1e-14, 1e-5)));
		REQUIRE(result31.value == Approx(expected).epsilon(TOL(1e-14, 1e-5)));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        EDGE CASES                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Gauss-Kronrod edge cases", "[gauss-kronrod][integration]")
{
	SECTION("Zero-width interval")
	{
		auto f = [](Real x) { return x * x; };
		auto result = IntegrateGK15(f, 1.0, 1.0);
		REQUIRE(result.value == Approx(0.0).margin(TOL(1e-15, 1e-5)));
	}
	
	SECTION("Reversed interval (b < a)")
	{
		// ∫₁⁰ x dx = -∫₀¹ x dx = -0.5
		auto f = [](Real x) { return x; };
		auto result = IntegrateGK15(f, 1.0, 0.0);
		REQUIRE(result.value == Approx(-0.5).epsilon(TOL(1e-14, 1e-5)));
	}
	
	SECTION("Large interval")
	{
		// ∫₀^100 sin(x)/x dx (converges slowly)
		auto f = [](Real x) -> Real { 
			if (std::abs(x) < Real(TOL(1e-10, 1e-5))) return Real(1.0);  // sinc(0) = 1
			return std::sin(x) / x; 
		};
		auto result = IntegrateGKAdaptive(f, 0.001, 100.0, 1e-6, 1e-6);
		// Si(100) - Si(0.001) ≈ 1.5622 (Sine integral)
		REQUIRE(result.converged);
		REQUIRE(result.value > 1.0);  // Sanity check
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        DETAILED API TESTS                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IntegrateGK15Detailed - basic success", "[gauss-kronrod][Detailed]")
{
	auto f = [](Real x) { return x * x; };
	auto result = IntegrateGK15Detailed(f, 0.0, 1.0);

	REQUIRE(result.IsSuccess());
	REQUIRE(result.algorithm_name == "IntegrateGK15");
	REQUIRE(result.elapsed_time_ms >= 0.0);
	REQUIRE(result.converged);
	REQUIRE(result.value == Approx(1.0/3.0).epsilon(TOL(1e-14, 1e-5)));
	REQUIRE(result.function_evaluations > 0);
}

TEST_CASE("IntegrateGK21Detailed - basic success", "[gauss-kronrod][Detailed]")
{
	auto f = [](Real x) { return std::sin(x); };
	auto result = IntegrateGK21Detailed(f, 0.0, Constants::PI);

	REQUIRE(result.IsSuccess());
	REQUIRE(result.algorithm_name == "IntegrateGK21");
	REQUIRE(result.value == Approx(2.0).epsilon(TOL(1e-15, 1e-5)));
}

TEST_CASE("IntegrateGK31Detailed - basic success", "[gauss-kronrod][Detailed]")
{
	auto f = [](Real x) { return std::cos(x) * std::cos(x); };
	auto result = IntegrateGK31Detailed(f, 0.0, Constants::PI);

	REQUIRE(result.IsSuccess());
	REQUIRE(result.algorithm_name == "IntegrateGK31");
	REQUIRE(result.value == Approx(Constants::PI / 2.0).epsilon(TOL(1e-14, 1e-5)));
}

TEST_CASE("IntegrateGKAdaptiveDetailed - peaked function", "[gauss-kronrod][Detailed]")
{
	auto f = [](Real x) {
		Real t = x - 0.5;
		return 1.0 / (1.0 + t * t / 0.001);
	};

	auto result = IntegrateGKAdaptiveDetailed(f, 0.0, 1.0);

	REQUIRE(result.IsSuccess());
	REQUIRE(result.algorithm_name == "IntegrateGKAdaptive");
	REQUIRE(result.converged);
	REQUIRE(result.function_evaluations > 15);  // Must subdivide

	Real sqrtEps = std::sqrt(0.001);
	Real expected = sqrtEps * 2.0 * std::atan(0.5 / sqrtEps);
	REQUIRE(result.value == Approx(expected).epsilon(TOL(1e-6, 1e-4)));
}

TEST_CASE("IntegrateGKAdaptiveDetailed - values match simple API", "[gauss-kronrod][Detailed]")
{
	auto f = [](Real x) { return std::exp(-x * x); };

	auto simple = IntegrateGKAdaptive(f, 0.0, 2.0, TOL(1e-10, 1e-5), TOL(1e-10, 1e-5), 50, GKRule::GK15);
	auto detailed = IntegrateGKAdaptiveDetailed(f, 0.0, 2.0, {}, TOL(1e-10, 1e-5), TOL(1e-10, 1e-5), 50, GKRule::GK15);

	REQUIRE(detailed.value == Approx(simple.value).epsilon(TOL(1e-14, 1e-5)));
	REQUIRE(detailed.function_evaluations == simple.function_evals);
}
