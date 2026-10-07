///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        combinatorics_tests.cpp                                             ///
///  Description: Tests for Combinatorics.h (binomials, sequences, partitions, enum)  ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"
#include "../../mml/base/Combinatorics.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

using namespace MML;
using namespace MML::Combinatorics;
using Catch::Approx;

namespace MML::Tests::Base::CombinatoricsTests {

TEST_CASE("Combinatorics_binomial", "[combinatorics]") {
	REQUIRE(BinomialCoefficient(5, 2) == 10);
	REQUIRE(BinomialCoefficient(10, 5) == 252);
	REQUIRE(BinomialCoefficient(52, 5) == 2598960);
	REQUIRE(BinomialCoefficient(0, 0) == 1);
	REQUIRE(BinomialCoefficient(5, 0) == 1);
	REQUIRE(BinomialCoefficient(5, 6) == 0);
	REQUIRE(BinomialCoefficient(5, -1) == 0);

	// exceeds 64 bits -> throws; the Real variant still works
	REQUIRE_THROWS_AS(BinomialCoefficient(70, 35), DomainError);
	REQUIRE(BinomialCoefficientReal(10, 5) == Approx(252.0));
	REQUIRE(BinomialCoefficientReal(70, 35) > 0.0);
}

TEST_CASE("Combinatorics_multinomial_factorials", "[combinatorics]") {
	REQUIRE(Multinomial({ 2, 2, 2 }) == 90);   // 6! / (2!2!2!)
	REQUIRE(Multinomial({ 1, 1, 1 }) == 6);

	REQUIRE(FallingFactorial(5, 3) == 60);     // 5*4*3
	REQUIRE(FallingFactorial(7, 0) == 1);
	REQUIRE(RisingFactorial(5, 3) == 210);     // 5*6*7
}

TEST_CASE("Combinatorics_catalan_fibonacci_lucas", "[combinatorics]") {
	REQUIRE(Catalan(0) == 1);
	REQUIRE(Catalan(4) == 14);
	REQUIRE(Catalan(5) == 42);

	REQUIRE(Fibonacci(0) == 0);
	REQUIRE(Fibonacci(1) == 1);
	REQUIRE(Fibonacci(10) == 55);
	REQUIRE(Fibonacci(20) == 6765);
	REQUIRE(Fibonacci(90) == 2880067194370816120LL);

	REQUIRE(Lucas(0) == 2);
	REQUIRE(Lucas(1) == 1);
	REQUIRE(Lucas(2) == 3);
	REQUIRE(Lucas(10) == 123);
}

TEST_CASE("Combinatorics_bell_stirling", "[combinatorics]") {
	REQUIRE(BellNumber(0) == 1);
	REQUIRE(BellNumber(3) == 5);
	REQUIRE(BellNumber(5) == 52);

	REQUIRE(StirlingFirst(0, 0) == 1);
	REQUIRE(StirlingFirst(4, 1) == 6);     // (4-1)!
	REQUIRE(StirlingFirst(4, 2) == 11);
	REQUIRE(StirlingFirst(5, 2) == 50);
	REQUIRE(StirlingFirst(4, 4) == 1);

	REQUIRE(StirlingSecond(0, 0) == 1);
	REQUIRE(StirlingSecond(4, 2) == 7);
	REQUIRE(StirlingSecond(5, 3) == 25);
	REQUIRE(StirlingSecond(5, 1) == 1);
	REQUIRE(StirlingSecond(4, 4) == 1);
}

TEST_CASE("Combinatorics_zigzag_tangent", "[combinatorics]") {
	REQUIRE(ZigzagNumber(0) == 1);
	REQUIRE(ZigzagNumber(3) == 2);
	REQUIRE(ZigzagNumber(4) == 5);
	REQUIRE(ZigzagNumber(5) == 16);

	REQUIRE(Tangent(1) == 1);
	REQUIRE(Tangent(2) == 2);
	REQUIRE(Tangent(3) == 16);
	REQUIRE(Tangent(4) == 272);
}

TEST_CASE("Combinatorics_bernoulli_exact", "[combinatorics]") {
	REQUIRE(Bernoulli(0) == Rational(1));
	REQUIRE(Bernoulli(1) == Rational(1, 2));    // + convention
	REQUIRE(Bernoulli(2) == Rational(1, 6));
	REQUIRE(Bernoulli(3) == Rational(0));
	REQUIRE(Bernoulli(4) == Rational(-1, 30));
	REQUIRE(Bernoulli(6) == Rational(1, 42));
	REQUIRE(Bernoulli(8) == Rational(-1, 30));
	REQUIRE(Bernoulli(10) == Rational(5, 66));
}

TEST_CASE("Combinatorics_bernoulli_real", "[combinatorics]") {
	REQUIRE(BernoulliReal(0) == Approx(1.0));
	REQUIRE(BernoulliReal(1) == Approx(0.5));
	REQUIRE(BernoulliReal(3) == Approx(0.0));
	REQUIRE(BernoulliReal(2) == Approx(1.0 / 6.0).epsilon(1e-9));
	REQUIRE(BernoulliReal(4) == Approx(-1.0 / 30.0).epsilon(1e-9));
	REQUIRE(BernoulliReal(12) == Approx(-691.0 / 2730.0).epsilon(1e-9));
	REQUIRE(BernoulliReal(20) == Approx(-174611.0 / 330.0).epsilon(1e-9));
	// large n: exact int64 rational overflows, so the zeta-formula fallback is used
	REQUIRE(BernoulliReal(40) == Approx(-261082718496449122051.0 / 13530.0).epsilon(1e-6));
}

TEST_CASE("Combinatorics_partitions", "[combinatorics]") {
	REQUIRE(PartitionCount(0) == 1);
	REQUIRE(PartitionCount(1) == 1);
	REQUIRE(PartitionCount(4) == 5);
	REQUIRE(PartitionCount(5) == 7);
	REQUIRE(PartitionCount(10) == 42);
	REQUIRE(PartitionCount(100) == 190569292LL);

	int count = 0;
	bool allSumTo5 = true;
	ForEachPartition(5, [&](const std::vector<int>& part) {
		++count;
		int s = 0;
		for (int p : part) s += p;
		if (s != 5) allSumTo5 = false;
	});
	REQUIRE(count == 7);
	REQUIRE(allSumTo5);
}

TEST_CASE("Combinatorics_enumeration", "[combinatorics]") {
	int combCount = 0;
	std::vector<int> first, last;
	ForEachCombination(5, 2, [&](const std::vector<int>& c) {
		if (combCount == 0) first = c;
		last = c;
		++combCount;
	});
	REQUIRE(combCount == 10);
	REQUIRE(first == std::vector<int>{ 0, 1 });
	REQUIRE(last == std::vector<int>{ 3, 4 });

	int multiCount = 0;
	ForEachMultiset(3, 2, [&](const std::vector<int>&) { ++multiCount; });
	REQUIRE(multiCount == 6);   // C(3+2-1, 2)

	int permCount = 0;
	ForEachPermutation(4, [&](const std::vector<int>&) { ++permCount; });
	REQUIRE(permCount == 24);
}

} // namespace MML::Tests::Base::CombinatoricsTests
