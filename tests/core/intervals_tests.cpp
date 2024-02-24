#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Intervals.h"
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

TEST_CASE("Test_CompleR", "[simple]") 
{
    CompleteRInterval int1;

    auto a = int1.getLowerBound();
}

TEST_CASE("Test_Interval", "[simple]") 
{
    // old way definition
    Interval *tangDefInterval = new Interval({  new ClosedOpenInterval(-2.0*Constants::PI, -1.5*Constants::PI), 
                                        new OpenInterval(-1.5*Constants::PI, -0.5*Constants::PI),
                                        new OpenInterval(-0.5*Constants::PI, 0.5*Constants::PI), 
                                        new OpenInterval(0.5*Constants::PI, 1.5*Constants::PI),
                                        new OpenClosedInterval(1.5*Constants::PI, 2.0*Constants::PI) } );

    REQUIRE(tangDefInterval->contains(-2.0*Constants::PI));
}

TEST_CASE("Test_ClosedIntervalWithReccuringPointHoles", "[simple]") 
{
    BaseInterval *tanDef = new ClosedIntervalWithReccuringPointHoles(-5*Constants::PI, 5*Constants::PI, 0.5*Constants::PI, Constants::PI);

    REQUIRE(tanDef->contains(0.0));

    // investigation around half PI
    REQUIRE_FALSE(tanDef->contains(0.5*Constants::PI));
    REQUIRE(tanDef->contains(0.5*Constants::PI + 0.0001));
    REQUIRE(tanDef->contains(0.5*Constants::PI + 1e-8));
    REQUIRE(tanDef->contains(0.5*Constants::PI + 1e-15));
    REQUIRE_FALSE(tanDef->contains(0.5*Constants::PI + 1e-16));
    REQUIRE(tanDef->contains(0.5*Constants::PI - 0.0001));
    REQUIRE(tanDef->contains(0.5*Constants::PI - 1e-8));
    REQUIRE(tanDef->contains(0.5*Constants::PI - 1e-15));
    REQUIRE_FALSE(tanDef->contains(0.5*Constants::PI - 1e-16));

    // verification of all the holes
    REQUIRE_FALSE(tanDef->contains(-0.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains( 1.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains(-1.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains( 2.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains(-2.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains( 3.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains(-3.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains( 4.5*Constants::PI));
    REQUIRE_FALSE(tanDef->contains(-5.5*Constants::PI));

    // precision verification on end holes
    REQUIRE(tanDef->contains(4.5*Constants::PI + 1e-15));
    REQUIRE_FALSE(tanDef->contains(4.5*Constants::PI + 1e-16));
    REQUIRE(tanDef->contains(4.5*Constants::PI - 1e-15));
    REQUIRE_FALSE(tanDef->contains(4.5*Constants::PI - 1e-16));

}

