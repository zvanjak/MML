///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RichardsonExtrapolationTests.cpp                                    ///
///  Description: Unit tests for Richardson extrapolation algorithms                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "../TestPrecision.h"
#include <catch2/catch_approx.hpp>

#include <core/RichardsonExtrapolation.h>
#include "MMLBase.h"

#include <cmath>

using namespace MML;
using namespace Catch;

namespace MML::Tests::Core::RichardsonExtrapolationTests {

///////////////////////////////////////////////////////////////////////////////////////////
///                            NEVILLE EXTRAPOLATION TESTS                             ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Richardson::NevilleExtrapolate - basic polynomial extrapolation", "[richardson]")
{
	SECTION("Linear extrapolation to zero")
	{
		// f(x) = 2x + 5, want f(0) = 5
		std::vector<Real> y = { 7.0, 9.0, 11.0 };  // f(1)=7, f(2)=9, f(3)=11
		std::vector<Real> x = { 1.0, 2.0, 3.0 };
		
		auto result = Richardson::NevilleExtrapolate(y, x);
		REQUIRE(result.value == Approx(5.0).epsilon(TOL(1e-12, 1e-5)));
	}

	SECTION("Quadratic extrapolation to zero")
	{
		// f(x) = x² + 2x + 3, want f(0) = 3
		std::vector<Real> y = { 6.0, 11.0, 18.0, 27.0 };  // f(1), f(2), f(3), f(4)
		std::vector<Real> x = { 1.0, 2.0, 3.0, 4.0 };
		
		auto result = Richardson::NevilleExtrapolate(y, x);
		REQUIRE(result.value == Approx(3.0).epsilon(TOL(1e-10, 1e-5)));
	}

	SECTION("Romberg-style: h² extrapolation")
	{
		// Simulating trapezoidal integration where error ~ h²
		// If we have values at h, h/2, h/4 with h²-error
		// x values are h² positions
		Real h = 0.5;
		std::vector<Real> x_vals;
		std::vector<Real> y_vals;
		
		// Suppose true value is π, estimates approach it as h→0
		Real true_val = Constants::PI;
		for (int k = 0; k < 5; k++)
		{
			Real hk = h / std::pow(2.0, k);
			x_vals.push_back(hk * hk);  // h² position
			// Simulate: estimate = true + error*h²
			y_vals.push_back(true_val + 0.1 * hk * hk);
		}
		
		auto result = Richardson::NevilleExtrapolate(y_vals, x_vals);
		REQUIRE(result.value == Approx(true_val).epsilon(TOL(1e-10, 1e-5)));
		REQUIRE(result.converged);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         RICHARDSON TABLEAU (DFRIDR-STYLE) TESTS                    ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Richardson::Extrapolate - derivative estimation", "[richardson]")
{
	SECTION("Derivative of sin(x) at x=1")
	{
		// f'(x) = cos(x), true derivative at x=1 is cos(1) ≈ 0.5403023058681398
		Real x0 = 1.0;
		Real true_deriv = std::cos(x0);
		
		// Central difference evaluator: (f(x+h) - f(x-h)) / (2h)
		auto central_diff = [x0](Real h) -> Real {
			return (std::sin(x0 + h) - std::sin(x0 - h)) / (2.0 * h);
		};
		
		auto result = Richardson::Extrapolate(central_diff, 0.5, 1.4, 10, 2.0);
		
		REQUIRE(result.value == Approx(true_deriv).epsilon(TOL(1e-10, 1e-5)));
		REQUIRE(result.converged);
		REQUIRE(result.error_estimate < TOL(1e-10, 1e-5));
	}

	SECTION("Derivative of exp(x) at x=0")
	{
		// f'(0) = exp(0) = 1
		Real x0 = 0.0;
		Real true_deriv = 1.0;
		
		auto central_diff = [x0](Real h) -> Real {
			return (std::exp(x0 + h) - std::exp(x0 - h)) / (2.0 * h);
		};
		
		auto result = Richardson::Extrapolate(central_diff, 0.1, 1.4, 10, 2.0);
		
		REQUIRE(result.value == Approx(true_deriv).epsilon(TOL(1e-12, 1e-5)));
		REQUIRE(result.converged);
	}

	SECTION("Derivative of x³ at x=2")
	{
		// f(x) = x³, f'(x) = 3x², f'(2) = 12
		Real x0 = 2.0;
		Real true_deriv = 12.0;
		
		auto central_diff = [x0](Real h) -> Real {
			auto f = [](Real t) { return t * t * t; };
			return (f(x0 + h) - f(x0 - h)) / (2.0 * h);
		};
		
		auto result = Richardson::Extrapolate(central_diff, 0.5, 1.4, 10, 2.0);
		
		REQUIRE(result.value == Approx(true_deriv).epsilon(TOL(1e-10, 1e-5)));
	}

	SECTION("Second derivative of sin(x) at x=π/4")
	{
		// f''(x) = -sin(x), f''(π/4) = -√2/2 ≈ -0.7071067811865476
		Real x0 = Constants::PI / 4.0;
		Real true_deriv2 = -std::sin(x0);
		
		// Second derivative central difference: (f(x+h) - 2f(x) + f(x-h)) / h²
		auto central_diff2 = [x0](Real h) -> Real {
			return (std::sin(x0 + h) - 2.0 * std::sin(x0) + std::sin(x0 - h)) / (h * h);
		};
		
		auto result = Richardson::Extrapolate(central_diff2, 0.3, 1.4, 10, 2.0);
		
		REQUIRE(result.value == Approx(true_deriv2).epsilon(TOL(1e-8, 1e-4)));
	}
}

TEST_CASE("Richardson::Extrapolate - different shrink factors", "[richardson]")
{
	// Test with con=2.0 (classic halving)
	Real x0 = 1.0;
	Real true_deriv = std::cos(x0);
	
	auto central_diff = [x0](Real h) -> Real {
		return (std::sin(x0 + h) - std::sin(x0 - h)) / (2.0 * h);
	};
	
	SECTION("con = 2.0")
	{
		auto result = Richardson::Extrapolate(central_diff, 0.5, 2.0, 10, 2.0);
		REQUIRE(result.value == Approx(true_deriv).epsilon(TOL(1e-10, 1e-5)));
	}
	
	SECTION("con = 1.5")
	{
		auto result = Richardson::Extrapolate(central_diff, 0.5, 1.5, 10, 2.0);
		REQUIRE(result.value == Approx(true_deriv).epsilon(TOL(1e-10, 1e-5)));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           SIMPLE RICHARDSON TESTS                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Richardson::SimpleExtrapolate - known error order", "[richardson]")
{
	SECTION("Trapezoidal rule error order p=2")
	{
		// True integral of sin(x) from 0 to π is 2
		Real true_val = 2.0;
		
		// Simulated trapezoidal estimates (error ~ h²)
		Real T_h = 1.98;      // estimate at h
		Real T_h2 = 1.995;    // estimate at h/2 (better)
		
		Real improved = Richardson::SimpleExtrapolate(T_h, T_h2, 2.0, 2);
		
		// Improved should be closer to true value than T_h2
		REQUIRE(std::abs(improved - true_val) < std::abs(T_h2 - true_val));
	}

	SECTION("Simpson's rule error order p=4")
	{
		Real true_val = 2.0;
		Real S_h = 1.9998;
		Real S_h2 = 1.99999;
		
		Real improved = Richardson::SimpleExtrapolate(S_h, S_h2, 2.0, 4);
		
		REQUIRE(std::abs(improved - true_val) < std::abs(S_h2 - true_val));
	}
}

TEST_CASE("Richardson::SimpleErrorEstimate", "[richardson]")
{
	// With trapezoidal rule (p=2) and r=2:
	// error estimate = |T_{h/2} - T_h| / (2^2 - 1) = |T_{h/2} - T_h| / 3
	
	Real T_h = 1.9;
	Real T_h2 = 1.975;
	
	Real err = Richardson::SimpleErrorEstimate(T_h, T_h2, 2.0, 2);
	Real expected = std::abs(T_h2 - T_h) / 3.0;
	
	REQUIRE(err == Approx(expected).epsilon(TOL(1e-12, 1e-5)));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                                 EDGE CASES                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Richardson - edge cases", "[richardson]")
{
	SECTION("Neville with empty vectors")
	{
		std::vector<Real> empty;
		auto result = Richardson::NevilleExtrapolate(empty, empty);
		REQUIRE(result.value == 0.0);
		REQUIRE(!result.converged);
	}

	SECTION("Neville with single element")
	{
		std::vector<Real> y = { 5.0 };
		std::vector<Real> x = { 1.0 };
		auto result = Richardson::NevilleExtrapolate(y, x);
		REQUIRE(result.value == 5.0);
		REQUIRE(!result.converged);
	}

	SECTION("Extrapolate with h=0")
	{
		auto eval = [](Real h) { return 1.0; };
		auto result = Richardson::Extrapolate(eval, 0.0);
		REQUIRE(!result.converged);
	}

	SECTION("Extrapolate with constant function")
	{
		// f(x) = 5, derivative = 0
		Real x0 = 1.0;
		auto central_diff = [x0](Real h) -> Real {
			return (5.0 - 5.0) / (2.0 * h);  // Always 0
		};
		
		auto result = Richardson::Extrapolate(central_diff, 0.1);
		REQUIRE(result.value == Approx(0.0).margin(TOL(1e-15, 1e-5)));
	}
}

} // namespace MML::Tests::Core::RichardsonExtrapolationTests
