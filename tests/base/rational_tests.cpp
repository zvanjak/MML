///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        rational_tests.cpp                                                  ///
///  Description: Tests for Rational.h (exact int64-backed rational numbers)          ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"
#include "../../mml/base/Rational.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <sstream>

using namespace MML;
using Catch::Approx;

namespace MML::Tests::Base::RationalTests {

TEST_CASE("Rational_construction_and_reduction", "[rational]") {
	Rational a(2, 4);
	REQUIRE(a.Num() == 1);
	REQUIRE(a.Den() == 2);

	Rational b(-3, -4);   // sign pushed to numerator, denominator positive
	REQUIRE(b.Num() == 3);
	REQUIRE(b.Den() == 4);

	Rational c(3, -6);
	REQUIRE(c.Num() == -1);
	REQUIRE(c.Den() == 2);

	Rational zero(0, 5);
	REQUIRE(zero.Num() == 0);
	REQUIRE(zero.Den() == 1);

	Rational whole(7);
	REQUIRE(whole.Num() == 7);
	REQUIRE(whole.Den() == 1);
}

TEST_CASE("Rational_to_real", "[rational]") {
	REQUIRE(Rational(1, 4).ToReal() == Approx(0.25));
	REQUIRE(Rational(-2, 3).ToReal() == Approx(-2.0 / 3.0));
}

TEST_CASE("Rational_arithmetic", "[rational]") {
	REQUIRE(Rational(1, 2) + Rational(1, 3) == Rational(5, 6));
	REQUIRE(Rational(1, 2) - Rational(1, 3) == Rational(1, 6));
	REQUIRE(Rational(2, 3) * Rational(3, 4) == Rational(1, 2));
	REQUIRE(Rational(1, 2) / Rational(3, 4) == Rational(2, 3));
	REQUIRE(-Rational(3, 5) == Rational(-3, 5));

	Rational acc(0);
	acc += Rational(1, 2);
	acc += Rational(1, 3);
	acc += Rational(1, 6);
	REQUIRE(acc == Rational(1));   // 1/2 + 1/3 + 1/6 = 1
}

TEST_CASE("Rational_mixed_integer", "[rational]") {
	REQUIRE(Rational(1, 2) + 1 == Rational(3, 2));
	REQUIRE(Rational(5) == Rational(10, 2));
}

TEST_CASE("Rational_comparisons", "[rational]") {
	REQUIRE(Rational(1, 3) < Rational(1, 2));
	REQUIRE(Rational(2, 4) == Rational(1, 2));
	REQUIRE(Rational(-1, 2) < Rational(0));
	REQUIRE(Rational(3, 4) > Rational(2, 3));
	REQUIRE(Rational(1, 2) <= Rational(1, 2));
	REQUIRE(Rational(1, 2) >= Rational(1, 2));
}

TEST_CASE("Rational_errors", "[rational]") {
	REQUIRE_THROWS_AS(Rational(1, 0), DivisionByZeroError);
	REQUIRE_THROWS_AS(Rational(1, 2) / Rational(0), DivisionByZeroError);
}

TEST_CASE("Rational_stream_output", "[rational]") {
	std::ostringstream os1;
	os1 << Rational(3, 4);
	REQUIRE(os1.str() == "3/4");

	std::ostringstream os2;
	os2 << Rational(5);
	REQUIRE(os2.str() == "5");
}

} // namespace MML::Tests::Base::RationalTests
