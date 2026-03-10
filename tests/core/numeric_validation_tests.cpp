///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        numeric_validation_tests.cpp                                        ///
///  Description: Tests for NumericValidation.h - Input validation helpers            ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_all.hpp>

#include "MMLBase.h"
#include "core/NumericValidation.h"

#include <limits>
#include <cmath>

using namespace MML;

namespace MML::Tests::Core::NumericValidationTests {

///////////////////////////////////////////////////////////////////////////////////////////
///                           IsFinite / IsFiniteValue                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsFinite - Normal values", "[NumericValidation][IsFinite]")
{
	REQUIRE(IsFinite(0.0));
	REQUIRE(IsFinite(1.0));
	REQUIRE(IsFinite(-1.0));
	REQUIRE(IsFinite(1e10));
	REQUIRE(IsFinite(-1e10));
	REQUIRE(IsFinite(1e-10));
	REQUIRE(IsFinite(std::numeric_limits<Real>::min()));
	REQUIRE(IsFinite(std::numeric_limits<Real>::max()));
}

TEST_CASE("IsFinite - NaN", "[NumericValidation][IsFinite]")
{
	REQUIRE_FALSE(IsFinite(std::numeric_limits<Real>::quiet_NaN()));
	REQUIRE_FALSE(IsFinite(std::numeric_limits<Real>::signaling_NaN()));
	REQUIRE_FALSE(IsFinite(std::nan("")));
}

TEST_CASE("IsFinite - Infinity", "[NumericValidation][IsFinite]")
{
	REQUIRE_FALSE(IsFinite(std::numeric_limits<Real>::infinity()));
	REQUIRE_FALSE(IsFinite(-std::numeric_limits<Real>::infinity()));
	REQUIRE_FALSE(IsFinite(HUGE_VAL));
}

TEST_CASE("IsFiniteValue - Template version", "[NumericValidation][IsFiniteValue]")
{
	REQUIRE(IsFiniteValue(0.0f));
	REQUIRE(IsFiniteValue(1.0));
	REQUIRE(IsFiniteValue(42));  // int
	
	REQUIRE_FALSE(IsFiniteValue(std::numeric_limits<float>::quiet_NaN()));
	REQUIRE_FALSE(IsFiniteValue(std::numeric_limits<double>::infinity()));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ValidateFinite                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ValidateFinite - Normal values pass", "[NumericValidation][ValidateFinite]")
{
	REQUIRE_NOTHROW(ValidateFinite(0.0, "x"));
	REQUIRE_NOTHROW(ValidateFinite(1.0, "value"));
	REQUIRE_NOTHROW(ValidateFinite(-1e10, "large negative"));
	REQUIRE_NOTHROW(ValidateFinite(std::numeric_limits<Real>::max(), "max"));
}

TEST_CASE("ValidateFinite - NaN throws", "[NumericValidation][ValidateFinite]")
{
	Real nan_val = std::numeric_limits<Real>::quiet_NaN();
	
	REQUIRE_THROWS_AS(ValidateFinite(nan_val, "x"), NumericInputError);
	
	try {
		ValidateFinite(nan_val, "test_param");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("test_param") != std::string::npos);
		REQUIRE(msg.find("NaN") != std::string::npos);
	}
}

TEST_CASE("ValidateFinite - Infinity throws", "[NumericValidation][ValidateFinite]")
{
	Real inf_val = std::numeric_limits<Real>::infinity();
	Real neg_inf = -std::numeric_limits<Real>::infinity();
	
	REQUIRE_THROWS_AS(ValidateFinite(inf_val, "x"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateFinite(neg_inf, "y"), NumericInputError);
	
	try {
		ValidateFinite(inf_val, "infinity_param");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("infinity_param") != std::string::npos);
		REQUIRE(msg.find("Inf") != std::string::npos);
	}
}

TEST_CASE("ValidateFinite - String version", "[NumericValidation][ValidateFinite]")
{
	std::string name = "my_parameter";
	
	REQUIRE_NOTHROW(ValidateFinite(1.0, name));
	REQUIRE_THROWS_AS(ValidateFinite(std::nan(""), name), NumericInputError);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ValidateBounds                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ValidateBounds - Normal bounds pass", "[NumericValidation][ValidateBounds]")
{
	REQUIRE_NOTHROW(ValidateBounds(0.0, 1.0, "test"));
	REQUIRE_NOTHROW(ValidateBounds(-10.0, 10.0, "range"));
	REQUIRE_NOTHROW(ValidateBounds(-1e10, 1e10, "wide"));
	
	// Equal bounds are allowed
	REQUIRE_NOTHROW(ValidateBounds(5.0, 5.0, "equal"));
	
	// Inverted bounds are allowed (caller's responsibility)
	REQUIRE_NOTHROW(ValidateBounds(10.0, 0.0, "inverted"));
}

TEST_CASE("ValidateBounds - NaN lower bound throws", "[NumericValidation][ValidateBounds]")
{
	Real nan_val = std::nan("");
	
	REQUIRE_THROWS_AS(ValidateBounds(nan_val, 1.0, "test"), NumericInputError);
	
	try {
		ValidateBounds(nan_val, 1.0, "bisection");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("bisection") != std::string::npos);
		REQUIRE(msg.find("lower bound") != std::string::npos);
	}
}

TEST_CASE("ValidateBounds - NaN upper bound throws", "[NumericValidation][ValidateBounds]")
{
	Real nan_val = std::nan("");
	
	REQUIRE_THROWS_AS(ValidateBounds(0.0, nan_val, "test"), NumericInputError);
	
	try {
		ValidateBounds(0.0, nan_val, "integration");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("integration") != std::string::npos);
		REQUIRE(msg.find("upper bound") != std::string::npos);
	}
}

TEST_CASE("ValidateBounds - Infinite bounds throw", "[NumericValidation][ValidateBounds]")
{
	Real inf = std::numeric_limits<Real>::infinity();
	
	REQUIRE_THROWS_AS(ValidateBounds(inf, 1.0, "test"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateBounds(-inf, 1.0, "test"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateBounds(0.0, inf, "test"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateBounds(0.0, -inf, "test"), NumericInputError);
}

TEST_CASE("ValidateBounds - String version", "[NumericValidation][ValidateBounds]")
{
	std::string ctx = "my_algorithm";
	
	REQUIRE_NOTHROW(ValidateBounds(0.0, 1.0, ctx));
	REQUIRE_THROWS_AS(ValidateBounds(std::nan(""), 1.0, ctx), NumericInputError);
}

TEST_CASE("ValidateBounds - Default context", "[NumericValidation][ValidateBounds]")
{
	REQUIRE_NOTHROW(ValidateBounds(0.0, 1.0));  // Default "algorithm"
	
	try {
		ValidateBounds(std::nan(""), 1.0);  // Default context
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("algorithm") != std::string::npos);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ValidateTolerance                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ValidateTolerance - Valid tolerances pass", "[NumericValidation][ValidateTolerance]")
{
	REQUIRE_NOTHROW(ValidateTolerance(1e-6, "test"));
	REQUIRE_NOTHROW(ValidateTolerance(1e-12, "high precision"));
	REQUIRE_NOTHROW(ValidateTolerance(0.1, "low precision"));
	REQUIRE_NOTHROW(ValidateTolerance(std::numeric_limits<Real>::min(), "smallest"));
}

TEST_CASE("ValidateTolerance - Zero throws", "[NumericValidation][ValidateTolerance]")
{
	REQUIRE_THROWS_AS(ValidateTolerance(0.0, "test"), NumericInputError);
}

TEST_CASE("ValidateTolerance - Negative throws", "[NumericValidation][ValidateTolerance]")
{
	REQUIRE_THROWS_AS(ValidateTolerance(-1e-6, "test"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateTolerance(-0.1, "test"), NumericInputError);
}

TEST_CASE("ValidateTolerance - NaN throws", "[NumericValidation][ValidateTolerance]")
{
	REQUIRE_THROWS_AS(ValidateTolerance(std::nan(""), "test"), NumericInputError);
}

TEST_CASE("ValidateTolerance - Infinity throws", "[NumericValidation][ValidateTolerance]")
{
	REQUIRE_THROWS_AS(ValidateTolerance(std::numeric_limits<Real>::infinity(), "test"), NumericInputError);
}

TEST_CASE("ValidateTolerance - Error message content", "[NumericValidation][ValidateTolerance]")
{
	try {
		ValidateTolerance(-0.001, "Newton iteration");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("Newton iteration") != std::string::npos);
		REQUIRE(msg.find("tolerance") != std::string::npos);
		REQUIRE(msg.find("positive") != std::string::npos);
	}
}

TEST_CASE("ValidateTolerance - String version", "[NumericValidation][ValidateTolerance]")
{
	std::string ctx = "my_tolerance";
	
	REQUIRE_NOTHROW(ValidateTolerance(1e-6, ctx));
	REQUIRE_THROWS_AS(ValidateTolerance(0.0, ctx), NumericInputError);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ValidateFunctionValue                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ValidateFunctionValue - Normal values pass", "[NumericValidation][ValidateFunctionValue]")
{
	REQUIRE_NOTHROW(ValidateFunctionValue(0.0, "f(x)"));
	REQUIRE_NOTHROW(ValidateFunctionValue(1.0, "f(x)"));
	REQUIRE_NOTHROW(ValidateFunctionValue(-1e10, "f(x)"));
	REQUIRE_NOTHROW(ValidateFunctionValue(1e10, "f(x)"));
}

TEST_CASE("ValidateFunctionValue - NaN throws", "[NumericValidation][ValidateFunctionValue]")
{
	REQUIRE_THROWS_AS(ValidateFunctionValue(std::nan(""), "f(x)"), NumericInputError);
	
	try {
		ValidateFunctionValue(std::nan(""), "function at x=0.5");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("function at x=0.5") != std::string::npos);
		REQUIRE(msg.find("NaN") != std::string::npos);
	}
}

TEST_CASE("ValidateFunctionValue - Infinity throws", "[NumericValidation][ValidateFunctionValue]")
{
	Real inf = std::numeric_limits<Real>::infinity();
	
	REQUIRE_THROWS_AS(ValidateFunctionValue(inf, "f(x)"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateFunctionValue(-inf, "f(x)"), NumericInputError);
	
	try {
		ValidateFunctionValue(inf, "derivative");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("derivative") != std::string::npos);
		REQUIRE(msg.find("Inf") != std::string::npos);
	}
}

TEST_CASE("ValidateFunctionValue - Default context", "[NumericValidation][ValidateFunctionValue]")
{
	REQUIRE_NOTHROW(ValidateFunctionValue(1.0));  // Default "function evaluation"
	
	try {
		ValidateFunctionValue(std::nan(""));
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("function evaluation") != std::string::npos);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           IsFunctionValueValid                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("IsFunctionValueValid - Normal values", "[NumericValidation][IsFunctionValueValid]")
{
	REQUIRE(IsFunctionValueValid(0.0));
	REQUIRE(IsFunctionValueValid(1.0));
	REQUIRE(IsFunctionValueValid(-1e10));
	REQUIRE(IsFunctionValueValid(std::numeric_limits<Real>::max()));
}

TEST_CASE("IsFunctionValueValid - Invalid values", "[NumericValidation][IsFunctionValueValid]")
{
	REQUIRE_FALSE(IsFunctionValueValid(std::nan("")));
	REQUIRE_FALSE(IsFunctionValueValid(std::numeric_limits<Real>::infinity()));
	REQUIRE_FALSE(IsFunctionValueValid(-std::numeric_limits<Real>::infinity()));
}

TEST_CASE("IsFunctionValueValid - No throw on invalid", "[NumericValidation][IsFunctionValueValid]")
{
	// This should never throw, even with bad input
	REQUIRE_NOTHROW(IsFunctionValueValid(std::nan("")));
	REQUIRE_NOTHROW(IsFunctionValueValid(std::numeric_limits<Real>::infinity()));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ValidateMaxIterations                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ValidateMaxIterations - Valid values pass", "[NumericValidation][ValidateMaxIterations]")
{
	REQUIRE_NOTHROW(ValidateMaxIterations(1, "test"));
	REQUIRE_NOTHROW(ValidateMaxIterations(100, "Newton"));
	REQUIRE_NOTHROW(ValidateMaxIterations(10000, "long iteration"));
	REQUIRE_NOTHROW(ValidateMaxIterations(std::numeric_limits<int>::max(), "huge"));
}

TEST_CASE("ValidateMaxIterations - Zero throws", "[NumericValidation][ValidateMaxIterations]")
{
	REQUIRE_THROWS_AS(ValidateMaxIterations(0, "test"), NumericInputError);
}

TEST_CASE("ValidateMaxIterations - Negative throws", "[NumericValidation][ValidateMaxIterations]")
{
	REQUIRE_THROWS_AS(ValidateMaxIterations(-1, "test"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateMaxIterations(-100, "test"), NumericInputError);
}

TEST_CASE("ValidateMaxIterations - Error message", "[NumericValidation][ValidateMaxIterations]")
{
	try {
		ValidateMaxIterations(0, "bisection");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("bisection") != std::string::npos);
		REQUIRE(msg.find("max_iterations") != std::string::npos);
		REQUIRE(msg.find("positive") != std::string::npos);
	}
}

TEST_CASE("ValidateMaxIterations - Default context", "[NumericValidation][ValidateMaxIterations]")
{
	try {
		ValidateMaxIterations(0);
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("algorithm") != std::string::npos);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ValidateStepSize                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ValidateStepSize - Valid values pass", "[NumericValidation][ValidateStepSize]")
{
	REQUIRE_NOTHROW(ValidateStepSize(0.01, "test"));
	REQUIRE_NOTHROW(ValidateStepSize(1e-6, "differentiation"));
	REQUIRE_NOTHROW(ValidateStepSize(1.0, "integration"));
	REQUIRE_NOTHROW(ValidateStepSize(std::numeric_limits<Real>::min(), "tiny"));
}

TEST_CASE("ValidateStepSize - Zero throws", "[NumericValidation][ValidateStepSize]")
{
	REQUIRE_THROWS_AS(ValidateStepSize(0.0, "test"), NumericInputError);
}

TEST_CASE("ValidateStepSize - Negative throws", "[NumericValidation][ValidateStepSize]")
{
	REQUIRE_THROWS_AS(ValidateStepSize(-0.01, "test"), NumericInputError);
	REQUIRE_THROWS_AS(ValidateStepSize(-1e-6, "test"), NumericInputError);
}

TEST_CASE("ValidateStepSize - NaN throws", "[NumericValidation][ValidateStepSize]")
{
	REQUIRE_THROWS_AS(ValidateStepSize(std::nan(""), "test"), NumericInputError);
}

TEST_CASE("ValidateStepSize - Infinity throws", "[NumericValidation][ValidateStepSize]")
{
	REQUIRE_THROWS_AS(ValidateStepSize(std::numeric_limits<Real>::infinity(), "test"), NumericInputError);
}

TEST_CASE("ValidateStepSize - Error message", "[NumericValidation][ValidateStepSize]")
{
	try {
		ValidateStepSize(-0.001, "Romberg");
		FAIL("Should have thrown");
	} catch (const NumericInputError& e) {
		std::string msg = e.what();
		REQUIRE(msg.find("Romberg") != std::string::npos);
		REQUIRE(msg.find("step size") != std::string::npos);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Edge Cases & Combined                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("NumericValidation - Denormalized values are valid", "[NumericValidation][edge]")
{
	Real denorm = std::numeric_limits<Real>::denorm_min();
	
	REQUIRE(IsFinite(denorm));
	REQUIRE_NOTHROW(ValidateFinite(denorm, "denorm"));
	REQUIRE_NOTHROW(ValidateTolerance(denorm, "denorm_tol"));
	REQUIRE(IsFunctionValueValid(denorm));
}

TEST_CASE("NumericValidation - Epsilon values are valid", "[NumericValidation][edge]")
{
	Real eps = std::numeric_limits<Real>::epsilon();
	
	REQUIRE(IsFinite(eps));
	REQUIRE_NOTHROW(ValidateTolerance(eps, "epsilon"));
}

TEST_CASE("NumericValidation - Near-overflow values are valid", "[NumericValidation][edge]")
{
	Real large = std::numeric_limits<Real>::max() / 2;
	
	REQUIRE(IsFinite(large));
	REQUIRE_NOTHROW(ValidateFinite(large, "large"));
	REQUIRE_NOTHROW(ValidateBounds(-large, large, "wide bounds"));
}

} // namespace MML::Tests::Core::NumericValidationTests
