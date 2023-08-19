#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "utilities/Intervals.h"
#endif

using namespace MML;

TEST_CASE("Test_OpenInterval", "[simple]") {
    OpenInterval int1(0.0, 1.0); 
    // da li je postavljeno sve kako treba
	REQUIRE(0.0 ==  int1.getLowerBound());
    REQUIRE(1.0 ==  int1.getUpperBound());
    REQUIRE(1.0 ==  int1.getLength());
    REQUIRE(true ==  int1.isContinuous());
}

TEST_CASE("Test_OpenInterval_contains", "[simple]") {
    OpenInterval int1(0.0, 1.0); 

    REQUIRE(true ==  int1.contains(0.5));
    REQUIRE(false ==  int1.contains(-0.5));
}

TEST_CASE("Test_OpenInterval_Covering", "[simple]") {
    // kad uzmeš Covering da granice NISU u skupu točaka
}

TEST_CASE("Test_IntervalUnion_simple", "[simple]") {
    // ispravan low i upp bound
    OpenInterval int1(0,1), int2(3,4);
    IntervalUnion union1(int1, int2);

    REQUIRE(0.0 ==  union1.getLowerBound());
    REQUIRE(4.0 ==  union1.getUpperBound());
}

TEST_CASE("Test_IntervalUnion_swap", "[simple]") {
    // ispravan low i upp bound
    OpenInterval int1(0,1), int2(3,4);
    IntervalUnion union1(int2, int1);

    // REQUIRE(0.0 ==  union1.getLowerBound());
    // REQUIRE(4.0 ==  union1.getUpperBound());
}

