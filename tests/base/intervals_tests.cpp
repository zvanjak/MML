#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

// Removed: using Catch::Approx; - now using precision-aware RealApprox

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Intervals.h>
#endif

#include <cmath>
#include <algorithm>

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::IntervalsTests {

TEST_CASE("Test_OpenInterval", "[intervals][open]") {
	TEST_PRECISION_INFO();
	OpenInterval int1(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
	// Verify basic properties
	REQUIRE(static_cast<MML::Real>(0.0) == int1.getLowerBound());
	REQUIRE(static_cast<MML::Real>(1.0) == int1.getUpperBound());
	REQUIRE(static_cast<MML::Real>(1.0) == int1.getLength());
	REQUIRE(int1.isContinuous());
}

TEST_CASE("Test_OpenInterval_contains", "[intervals][open]") {
	TEST_PRECISION_INFO();
	OpenInterval int1(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));

	REQUIRE(int1.contains(static_cast<MML::Real>(0.5)));
	REQUIRE_FALSE(int1.contains(-static_cast<MML::Real>(0.5)));
	REQUIRE_FALSE(int1.contains(static_cast<MML::Real>(0.0)));  // Open - excludes lower bound
	REQUIRE_FALSE(int1.contains(static_cast<MML::Real>(1.0)));  // Open - excludes upper bound
}

TEST_CASE("Test_OpenInterval_Covering", "[intervals][open][covering]") {
	TEST_PRECISION_INFO();
	OpenInterval interval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0));
	std::vector<Real> points;
	
	interval.getEquidistantCovering(5, points);
	
	// For open interval, endpoints should NOT be exactly at boundaries
	REQUIRE(points.size() == 5);
	REQUIRE(points[0] > static_cast<MML::Real>(0.0));    // Slightly inside lower bound
	REQUIRE(points[4] < static_cast<MML::Real>(10.0));   // Slightly inside upper bound
	
	// Middle points should be at expected positions (slightly offset at endpoints)
	REQUIRE_THAT(points[2], RealApprox(static_cast<MML::Real>(5.0)).epsilon(static_cast<MML::Real>(0.01)));
}

TEST_CASE("CompleteRInterval - Properties", "[intervals][infinite][complete]")
{
	TEST_PRECISION_INFO();
	CompleteRInterval interval;

	SECTION("Bounds are at numeric limits")
	{
		REQUIRE(interval.getLowerBound() == -std::numeric_limits<Real>::max());
		REQUIRE(interval.getUpperBound() == std::numeric_limits<Real>::max());
	}
	
	SECTION("Contains any real number")
	{
		REQUIRE(interval.contains(static_cast<MML::Real>(0.0)));
		const Real largeFinite = std::sqrt(std::numeric_limits<Real>::max());
		REQUIRE(interval.contains(-largeFinite));
		REQUIRE(interval.contains(largeFinite));
		REQUIRE(interval.contains(std::numeric_limits<Real>::min()));
		REQUIRE(interval.contains(std::numeric_limits<Real>::max()));
	}
	
	SECTION("Is continuous")
	{
		REQUIRE(interval.isContinuous());
	}
}

TEST_CASE("CompleteRWithReccuringPointHoles - Periodic Exclusions", "[intervals][infinite][holes]")
{
	TEST_PRECISION_INFO();
	// Tangent function domain: exclude x = π/2 + nπ
	CompleteRWithReccuringPointHoles tanDomain(static_cast<MML::Real>(0.5) * Constants::PI, Constants::PI);
	
	SECTION("Contains regular points")
	{
		REQUIRE(tanDomain.contains(static_cast<MML::Real>(0.0)));
		REQUIRE(tanDomain.contains(Constants::PI));
		REQUIRE(tanDomain.contains(-Constants::PI));
		REQUIRE(tanDomain.contains(static_cast<MML::Real>(2.0) * Constants::PI));
		REQUIRE(tanDomain.contains(static_cast<MML::Real>(0.25) * Constants::PI));
	}
	
	SECTION("Excludes hole points")
	{
		REQUIRE_FALSE(tanDomain.contains(static_cast<MML::Real>(0.5) * Constants::PI));   // First hole
		REQUIRE_FALSE(tanDomain.contains(-static_cast<MML::Real>(0.5) * Constants::PI));  // Hole at -π/2
		REQUIRE_FALSE(tanDomain.contains(static_cast<MML::Real>(1.5) * Constants::PI));   // Hole at 3π/2
		REQUIRE_FALSE(tanDomain.contains(-static_cast<MML::Real>(1.5) * Constants::PI));  // Hole at -3π/2
		REQUIRE_FALSE(tanDomain.contains(static_cast<MML::Real>(2.5) * Constants::PI));   // Hole at 5π/2
	}
	
	SECTION("Is not continuous")
	{
		REQUIRE_FALSE(tanDomain.isContinuous());
	}

	SECTION("Negative hole indices via std::round")
	{
		// Holes at negative multiples: -π/2, -3π/2, -5π/2, etc.
		// These require correct rounding (not truncation) of negative diff values
		REQUIRE_FALSE(tanDomain.contains(-static_cast<MML::Real>(2.5) * Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(-static_cast<MML::Real>(3.5) * Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(-static_cast<MML::Real>(4.5) * Constants::PI));
		// Points between negative holes should be contained
		REQUIRE(tanDomain.contains(-static_cast<MML::Real>(2.0) * Constants::PI));
		REQUIRE(tanDomain.contains(-static_cast<MML::Real>(3.0) * Constants::PI));
	}
}

TEST_CASE("Test_Interval - Composite with value-based intervals", "[intervals][compound]")
{
	TEST_PRECISION_INFO();
	Interval tangDefInterval;
	tangDefInterval
		.AddInterval(ClosedOpenInterval(-static_cast<MML::Real>(2.0) * Constants::PI, -static_cast<MML::Real>(1.5) * Constants::PI))
		.AddInterval(OpenInterval(-static_cast<MML::Real>(1.5) * Constants::PI, -static_cast<MML::Real>(0.5) * Constants::PI))
		.AddInterval(OpenInterval(-static_cast<MML::Real>(0.5) * Constants::PI, static_cast<MML::Real>(0.5) * Constants::PI))
		.AddInterval(OpenInterval(static_cast<MML::Real>(0.5) * Constants::PI, static_cast<MML::Real>(1.5) * Constants::PI))
		.AddInterval(OpenClosedInterval(static_cast<MML::Real>(1.5) * Constants::PI, static_cast<MML::Real>(2.0) * Constants::PI));

	// Left endpoint of first interval is closed
	REQUIRE(tangDefInterval.contains(-static_cast<MML::Real>(2.0) * Constants::PI));
	// Right endpoint of last interval is closed
	REQUIRE(tangDefInterval.contains(static_cast<MML::Real>(2.0) * Constants::PI));
	// Holes at ±0.5π, ±1.5π are excluded (open endpoints)
	REQUIRE_FALSE(tangDefInterval.contains(-static_cast<MML::Real>(0.5) * Constants::PI));
	REQUIRE_FALSE(tangDefInterval.contains(static_cast<MML::Real>(0.5) * Constants::PI));
	REQUIRE_FALSE(tangDefInterval.contains(-static_cast<MML::Real>(1.5) * Constants::PI));
	REQUIRE_FALSE(tangDefInterval.contains(static_cast<MML::Real>(1.5) * Constants::PI));
	// Interior points are contained
	REQUIRE(tangDefInterval.contains(static_cast<MML::Real>(0.0)));
	REQUIRE(tangDefInterval.contains(Constants::PI));
}

TEST_CASE("Test_ClosedIntervalWithReccuringPointHoles", "[intervals][holes]")
{
		TEST_PRECISION_INFO();
	BaseInterval* tanDef = new ClosedIntervalWithReccuringPointHoles(-5 * Constants::PI, 5 * Constants::PI, static_cast<MML::Real>(0.5) * Constants::PI, Constants::PI);

	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.0)));

	// investigation around half PI
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI + static_cast<MML::Real>(0.0001)));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI + TOL(1e-8, 1e-4)));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI + TOL(1e-13, 5e-5)));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI + TOL3(1e-15, 5e-6, 1e-18)));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI - static_cast<MML::Real>(0.0001)));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI - TOL(1e-8, 1e-4)));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI - TOL(1e-13, 5e-5)));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(0.5) * Constants::PI - TOL3(1e-15, 5e-6, 1e-18)));

	// verification of all the holes
	REQUIRE_FALSE(tanDef->contains(-static_cast<MML::Real>(0.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(1.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-static_cast<MML::Real>(1.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(2.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-static_cast<MML::Real>(2.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(3.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-static_cast<MML::Real>(3.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(4.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-static_cast<MML::Real>(5.5) * Constants::PI));

	// precision verification on end holes
	REQUIRE(tanDef->contains(static_cast<MML::Real>(4.5) * Constants::PI + TOL(1e-13, 5e-5)));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(4.5) * Constants::PI + TOL3(1e-15, 1e-5, 1e-18)));
	REQUIRE(tanDef->contains(static_cast<MML::Real>(4.5) * Constants::PI - TOL(1e-13, 5e-5)));
	REQUIRE_FALSE(tanDef->contains(static_cast<MML::Real>(4.5) * Constants::PI - TOL3(1e-15, 1e-5, 1e-18)));

}

//=============================================================================
// COMPREHENSIVE INTERVAL TESTS - Basic Interval Types
//=============================================================================

TEST_CASE("ClosedInterval - Basic Properties", "[intervals][closed]") {
		TEST_PRECISION_INFO();
	ClosedInterval interval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
	
	SECTION("Bounds and Length") {
		REQUIRE(interval.getLowerBound() == static_cast<MML::Real>(2.0));
		REQUIRE(interval.getUpperBound() == static_cast<MML::Real>(8.0));
		REQUIRE(interval.getLength() == static_cast<MML::Real>(6.0));
		REQUIRE(interval.isContinuous());
	}
	
	SECTION("Contains - Interior Points") {
		REQUIRE(interval.contains(static_cast<MML::Real>(2.0)));   // Lower bound (closed)
		REQUIRE(interval.contains(static_cast<MML::Real>(8.0)));   // Upper bound (closed)
		REQUIRE(interval.contains(static_cast<MML::Real>(5.0)));   // Middle point
		REQUIRE(interval.contains(static_cast<MML::Real>(2.001))); // Just inside lower
		REQUIRE(interval.contains(static_cast<MML::Real>(7.999))); // Just inside upper
	}
	
	SECTION("Contains - Exterior Points") {
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(1.999)));  // Just outside lower
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(8.001)));  // Just outside upper
		REQUIRE_FALSE(interval.contains(-static_cast<MML::Real>(10.0)));  // Far outside
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(100.0)));  // Far outside
	}
}

TEST_CASE("OpenInterval - Basic Properties", "[intervals][open]") {
		TEST_PRECISION_INFO();
	OpenInterval interval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
	
	SECTION("Bounds and Length") {
		REQUIRE(interval.getLowerBound() == static_cast<MML::Real>(2.0));
		REQUIRE(interval.getUpperBound() == static_cast<MML::Real>(8.0));
		REQUIRE(interval.getLength() == static_cast<MML::Real>(6.0));
		REQUIRE(interval.isContinuous());
	}
	
	SECTION("Contains - Excludes Endpoints") {
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(2.0))); // Lower bound excluded
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(8.0))); // Upper bound excluded
		REQUIRE(interval.contains(static_cast<MML::Real>(5.0)));       // Middle point
		REQUIRE(interval.contains(static_cast<MML::Real>(2.001)));     // Just inside lower
		REQUIRE(interval.contains(static_cast<MML::Real>(7.999)));     // Just inside upper
	}
}

TEST_CASE("OpenClosedInterval - Mixed Endpoints", "[intervals][open-closed]") {
		TEST_PRECISION_INFO();
	OpenClosedInterval interval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
	
	SECTION("Contains - Mixed Behavior") {
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(2.0))); // Lower bound excluded (open)
		REQUIRE(interval.contains(static_cast<MML::Real>(8.0)));       // Upper bound included (closed)
		REQUIRE(interval.contains(static_cast<MML::Real>(5.0)));       // Middle point
		REQUIRE(interval.contains(static_cast<MML::Real>(2.001)));
		REQUIRE(interval.contains(static_cast<MML::Real>(7.999)));
	}
}

TEST_CASE("ClosedOpenInterval - Mixed Endpoints", "[intervals][closed-open]") {
		TEST_PRECISION_INFO();
	ClosedOpenInterval interval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
	
	SECTION("Contains - Mixed Behavior") {
		REQUIRE(interval.contains(static_cast<MML::Real>(2.0)));       // Lower bound included (closed)
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(8.0))); // Upper bound excluded (open)
		REQUIRE(interval.contains(static_cast<MML::Real>(5.0)));       // Middle point
		REQUIRE(interval.contains(static_cast<MML::Real>(2.001)));
		REQUIRE(interval.contains(static_cast<MML::Real>(7.999)));
	}
}

TEST_CASE("Infinite Bound Intervals", "[intervals][infinite]") {
		TEST_PRECISION_INFO();
	SECTION("NegInfToClosedInterval") {
		NegInfToClosedInterval interval(static_cast<MML::Real>(5.0));
		
		REQUIRE(interval.contains(-1e100));  // Very negative
		REQUIRE(interval.contains(static_cast<MML::Real>(0.0)));
		REQUIRE(interval.contains(static_cast<MML::Real>(4.999)));
		REQUIRE(interval.contains(static_cast<MML::Real>(5.0)));     // Closed upper bound
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(5.001)));
	}
	
	SECTION("NegInfToOpenInterval") {
		NegInfToOpenInterval interval(static_cast<MML::Real>(5.0));
		
		REQUIRE(interval.contains(-1e100));
		REQUIRE(interval.contains(static_cast<MML::Real>(0.0)));
		REQUIRE(interval.contains(static_cast<MML::Real>(4.999)));
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(5.0))); // Open upper bound
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(5.001)));
	}
	
	SECTION("ClosedToInfInterval") {
		ClosedToInfInterval interval(static_cast<MML::Real>(5.0));
		
		REQUIRE(interval.contains(static_cast<MML::Real>(5.0)));      // Closed lower bound
		REQUIRE(interval.contains(static_cast<MML::Real>(5.001)));
		REQUIRE(interval.contains(static_cast<MML::Real>(100.0)));
		REQUIRE(interval.contains(1e100));    // Very positive
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(4.999)));
	}
	
	SECTION("OpenToInfInterval") {
		OpenToInfInterval interval(static_cast<MML::Real>(5.0));
		
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(5.0))); // Open lower bound
		REQUIRE(interval.contains(static_cast<MML::Real>(5.001)));
		REQUIRE(interval.contains(static_cast<MML::Real>(100.0)));
		REQUIRE(interval.contains(1e100));
		REQUIRE_FALSE(interval.contains(static_cast<MML::Real>(4.999)));
	}
	
	SECTION("CompleteRInterval") {
		CompleteRInterval interval;
		
		REQUIRE(interval.contains(-1e100));
		REQUIRE(interval.contains(static_cast<MML::Real>(0.0)));
		REQUIRE(interval.contains(1e100));
		REQUIRE(interval.contains(-static_cast<MML::Real>(0.000001)));
		REQUIRE(interval.contains(static_cast<MML::Real>(0.000001)));
	}
}

//=============================================================================
// EQUIDISTANT COVERING TESTS
//=============================================================================

TEST_CASE("GetEquidistantCovering - Closed Interval", "[intervals][covering][closed]") {
		TEST_PRECISION_INFO();
	ClosedInterval interval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0));
	std::vector<Real> points;
	
	SECTION("5 Points") {
		interval.getEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE_THAT(points[0], RealApprox(static_cast<MML::Real>(0.0)));
		REQUIRE_THAT(points[1], RealApprox(static_cast<MML::Real>(2.5)));
		REQUIRE_THAT(points[2], RealApprox(static_cast<MML::Real>(5.0)));
		REQUIRE_THAT(points[3], RealApprox(static_cast<MML::Real>(7.5)));
		REQUIRE_THAT(points[4], RealApprox(static_cast<MML::Real>(10.0)));
	}
	
	SECTION("Single Point") {
		interval.getEquidistantCovering(1, points);
		
		REQUIRE(points.size() == 1);
		REQUIRE_THAT(points[0], RealApprox(static_cast<MML::Real>(5.0))); // Midpoint
	}
	
	SECTION("Empty Request") {
		interval.getEquidistantCovering(0, points);
		REQUIRE(points.empty());
	}
	
	SECTION("Many Points") {
		interval.getEquidistantCovering(101, points);
		
		REQUIRE(points.size() == 101);
		REQUIRE_THAT(points[0], RealApprox(static_cast<MML::Real>(0.0)));
		REQUIRE_THAT(points[50], RealApprox(static_cast<MML::Real>(5.0)));
		REQUIRE_THAT(points[100], RealApprox(static_cast<MML::Real>(10.0)));
		
		// Verify equidistant spacing
		Real delta = points[1] - points[0];
		for (size_t i = 1; i < points.size() - 1; i++) {
			REQUIRE_THAT(points[i+1] - points[i], RealApprox(delta).epsilon(TOL(1e-9, 1e-5)));
		}
	}
}

TEST_CASE("GetEquidistantCovering - Open Interval", "[intervals][covering][open]") {
		TEST_PRECISION_INFO();
	OpenInterval interval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0));
	std::vector<Real> points;
	
	SECTION("Points Avoid Endpoints") {
		interval.getEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE(points[0] > static_cast<MML::Real>(0.0));    // Not exactly 0 (open)
		REQUIRE(points[4] < static_cast<MML::Real>(10.0));   // Not exactly 10 (open)
		REQUIRE_THAT(points[2], RealApprox(static_cast<MML::Real>(5.0)).epsilon(static_cast<MML::Real>(0.01))); // Middle point
	}
}

TEST_CASE("GetEquidistantCovering - Infinite Bounds", "[intervals][covering][infinite]") {
		TEST_PRECISION_INFO();
	SECTION("Infinite Lower Bound") {
		NegInfToClosedInterval interval(static_cast<MML::Real>(10.0));
		std::vector<Real> points;
		
		interval.getEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE(points[0] < points[4]); // Ascending order
		REQUIRE_THAT(points[4], RealApprox(static_cast<MML::Real>(10.0))); // Upper bound
		Real expectedLower = std::is_same_v<Real, float> ? Real(-1e6) : Real(-1e10);
		REQUIRE_THAT(points[0], RealApprox(expectedLower).epsilon(static_cast<MML::Real>(0.01))); // Practical lower limit
	}
	
	SECTION("Infinite Upper Bound") {
		ClosedToInfInterval interval(static_cast<MML::Real>(0.0));
		std::vector<Real> points;
		
		interval.getEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE_THAT(points[0], RealApprox(static_cast<MML::Real>(0.0)));
		Real expectedUpper = std::is_same_v<Real, float> ? Real(1e6) : Real(1e10);
		REQUIRE_THAT(points[4], RealApprox(expectedUpper).epsilon(static_cast<MML::Real>(0.01))); // Practical upper limit
	}
}

//=============================================================================
// INTERSECTION TESTS
//=============================================================================

TEST_CASE("Interval::Intersection - Basic Cases", "[intervals][intersection]") {
		TEST_PRECISION_INFO();
	SECTION("Complete Overlap - Same Intervals") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		ClosedInterval b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(8.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(2.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(8.0)));
	}
	
	SECTION("Partial Overlap") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		ClosedInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(10.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(5.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(8.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(6.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(8.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(4.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(9.0)));
	}
	
	SECTION("No Overlap - Gap Between") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
		ClosedInterval b(static_cast<MML::Real>(7.0), static_cast<MML::Real>(10.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(3.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(8.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(6.0))); // In the gap
	}
	
	SECTION("One Contains Other") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0));
		ClosedInterval b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(4.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(6.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
	}
	
	SECTION("Touching at Boundary") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
		ClosedInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE(result.contains(static_cast<MML::Real>(5.0))); // Touching point included (both closed)
		REQUIRE_THAT(result.getLength(), RealApprox(static_cast<MML::Real>(0.0))); // Degenerate interval
	}
}

TEST_CASE("Interval::Intersection - Endpoint Types", "[intervals][intersection][endpoints]") {
		TEST_PRECISION_INFO();
	SECTION("Open meets Closed - Lower Bound") {
		OpenInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		ClosedInterval b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(2.0))); // Open endpoint wins
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(8.0))); // Open endpoint wins (a is open at upper)
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
	}
	
	SECTION("Open meets Closed - Upper Bound") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		OpenInterval b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(2.0))); // Open endpoint wins
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(8.0))); // Open endpoint wins
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
	}
	
	SECTION("Mixed Endpoint Types") {
		OpenClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(8.0));  // (2, 8]
		ClosedOpenInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(10.0)); // [5, 10)
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(5.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(8.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));   // Closed from b
		REQUIRE(result.contains(static_cast<MML::Real>(8.0)));   // Closed from a
		REQUIRE(result.contains(static_cast<MML::Real>(6.5)));
	}
}

//=============================================================================
// UNION TESTS
//=============================================================================

TEST_CASE("Interval::Union - Overlapping Intervals", "[intervals][union]") {
	TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(6.0));
	ClosedInterval b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(8.0));
	
	Interval result = Interval::Union(a, b);
	
	SECTION("Returns Merged Interval") {
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(8.0)));
		REQUIRE_THAT(result.getLength(), RealApprox(static_cast<MML::Real>(6.0)));    // Hull
		REQUIRE_THAT(result.getMeasure(), RealApprox(static_cast<MML::Real>(6.0)));   // Measure (single interval)
		REQUIRE(result.contains(static_cast<MML::Real>(2.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(8.0)));
	}
}

TEST_CASE("Interval::Union - Adjacent Intervals", "[intervals][union]") {
	TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
	ClosedInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(8.0));
	
	Interval result = Interval::Union(a, b);
	
	SECTION("Returns Merged Interval") {
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(8.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));  // Shared endpoint included
	}
}

TEST_CASE("Interval::Union - Disjoint Intervals", "[intervals][union]") {
	TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0));
	ClosedInterval b(static_cast<MML::Real>(6.0), static_cast<MML::Real>(8.0));
	
	Interval result = Interval::Union(a, b);
	
	SECTION("Returns Composite Interval") {
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(8.0)));
		REQUIRE_THAT(result.getLength(), RealApprox(static_cast<MML::Real>(6.0)));    // Hull: 8 - 2 = 6
		REQUIRE_THAT(result.getMeasure(), RealApprox(static_cast<MML::Real>(4.0)));   // Measure: 2 + 2 = 4
		REQUIRE(result.contains(static_cast<MML::Real>(2.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(3.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0)));  // In gap
		REQUIRE(result.contains(static_cast<MML::Real>(7.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(8.0)));
	}
}

TEST_CASE("Interval::Union - One Contains Other", "[intervals][union]") {
	TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0));
	ClosedInterval b(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));
	
	Interval result = Interval::Union(a, b);
	
	SECTION("Returns Larger Interval") {
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(10.0)));
		REQUIRE_THAT(result.getMeasure(), RealApprox(static_cast<MML::Real>(8.0)));
	}
}

TEST_CASE("Interval::Union - Open Endpoints Not Adjacent", "[intervals][union]") {
	TEST_PRECISION_INFO();
	OpenInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));   // (2, 5)
	OpenInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(8.0));   // (5, 8)
	
	Interval result = Interval::Union(a, b);
	
	SECTION("Returns Disjoint (Gap at 5)") {
		// Both are open at 5, so they don't touch
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0)));  // Gap at 5
		REQUIRE(result.contains(static_cast<MML::Real>(4.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(6.0)));
	}
}

//=============================================================================
// DIFFERENCE TESTS
//=============================================================================

TEST_CASE("Interval::Difference - No Intersection", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
	ClosedInterval b(static_cast<MML::Real>(7.0), static_cast<MML::Real>(10.0));
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Original Interval") {
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(5.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(2.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(3.5)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(8.0)));
	}
}

TEST_CASE("Interval::Difference - Complete Containment", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));
	ClosedInterval b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0)); // b completely contains a
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Empty Set") {
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(4.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(6.0)));
	}
}

TEST_CASE("Interval::Difference - Left Cut", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(5.0), static_cast<MML::Real>(10.0));
	ClosedInterval b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(7.0)); // Cuts left part of a
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Right Remainder") {
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(6.0)));  // Removed by b
		REQUIRE(result.contains(static_cast<MML::Real>(7.001)));      // Remaining part
		REQUIRE(result.contains(static_cast<MML::Real>(10.0)));
		// When b=[2,7] removes from a=[5,10], result is (7,10]
		// b's upper is CLOSED at 7, so result starts OPEN at 7
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(7.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(7.0))); // Boundary is open
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(10.0)));
	}
}

TEST_CASE("Interval::Difference - Right Cut", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0));
	ClosedInterval b(static_cast<MML::Real>(7.0), static_cast<MML::Real>(15.0)); // Cuts right part of a
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Left Remainder") {
		REQUIRE(result.contains(static_cast<MML::Real>(2.0)));
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(7.0)));  // Cut here (closed from b)
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(9.0)));  // Removed by b
		REQUIRE_THAT(result.getLowerBound(), RealApprox(static_cast<MML::Real>(2.0)));
		// When b=[7,15] removes from a=[2,10], result is [2,7)
		// b's lower is CLOSED at 7, so result ends OPEN at 7
		REQUIRE_THAT(result.getUpperBound(), RealApprox(static_cast<MML::Real>(7.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(7.0))); // Boundary is open
	}
}

TEST_CASE("Interval::Difference - Middle Cut (Split)", "[intervals][difference]") {
	TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0));
	ClosedInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(7.0)); // Cuts middle, creating two pieces
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Two Intervals") {
		REQUIRE(result.contains(static_cast<MML::Real>(2.0)));      // Left piece
		REQUIRE(result.contains(static_cast<MML::Real>(4.999)));    // Left piece edge
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0))); // Removed (closed from b)
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(6.0))); // Middle removed
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(7.0))); // Removed (closed from b)
		REQUIRE(result.contains(static_cast<MML::Real>(7.001)));    // Right piece edge
		REQUIRE(result.contains(static_cast<MML::Real>(10.0)));     // Right piece
		
		// Hull length is span from 2 to 10
		REQUIRE_THAT(result.getLength(), RealApprox(static_cast<MML::Real>(8.0)));  // Hull: 10 - 2 = 8
		// Measure is sum of segments: (5-2) + (10-7) = 3 + 3 = 6
		REQUIRE_THAT(result.getMeasure(), RealApprox(static_cast<MML::Real>(6.0)));
	}
}

TEST_CASE("Interval::Difference - Endpoint Handling", "[intervals][difference][endpoints]") {
		TEST_PRECISION_INFO();
	SECTION("Open Interval Removes Closed Interval") {
		ClosedInterval a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0));
		OpenInterval b(static_cast<MML::Real>(5.0), static_cast<MML::Real>(7.0)); // Open boundaries
		
		Interval result = Interval::Difference(a, b);
		
		REQUIRE(result.contains(static_cast<MML::Real>(5.0)));  // NOT removed (b is open at 5)
		REQUIRE(result.contains(static_cast<MML::Real>(7.0)));  // NOT removed (b is open at 7)
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(6.0))); // Removed (inside b)
	}
}

//=============================================================================
// COMPLEMENT TESTS
//=============================================================================

TEST_CASE("Interval::Complement - Finite Interval", "[intervals][complement]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(3.0), static_cast<MML::Real>(7.0));
	
	Interval result = Interval::Complement(a);
	
	SECTION("Creates Two Intervals") {
		REQUIRE(result.contains(-static_cast<MML::Real>(1000.0)));  // Left piece (to -inf)
		REQUIRE(result.contains(static_cast<MML::Real>(2.999)));    // Just before lower bound
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(3.0))); // Lower bound excluded
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0))); // Middle excluded
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(7.0))); // Upper bound excluded
		REQUIRE(result.contains(static_cast<MML::Real>(7.001)));    // Just after upper bound
		REQUIRE(result.contains(static_cast<MML::Real>(1000.0)));   // Right piece (to +inf)
	}
}

TEST_CASE("Interval::Complement - Open Interval", "[intervals][complement]") {
		TEST_PRECISION_INFO();
	OpenInterval a(static_cast<MML::Real>(3.0), static_cast<MML::Real>(7.0));
	
	Interval result = Interval::Complement(a);
	
	SECTION("Includes Endpoints") {
		REQUIRE(result.contains(static_cast<MML::Real>(3.0)));  // Included (a is open at 3)
		REQUIRE(result.contains(static_cast<MML::Real>(7.0)));  // Included (a is open at 7)
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0))); // Excluded (inside a)
	}
}

TEST_CASE("Interval::Complement - Infinite Bounds", "[intervals][complement]") {
		TEST_PRECISION_INFO();
	SECTION("Left Infinite") {
		NegInfToClosedInterval a(static_cast<MML::Real>(5.0));
		
		Interval result = Interval::Complement(a);
		
		REQUIRE_FALSE(result.contains(-static_cast<MML::Real>(1000.0))); // Part of a
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0)));     // Part of a (closed)
		REQUIRE(result.contains(static_cast<MML::Real>(5.001)));         // Complement (to +inf)
		REQUIRE(result.contains(static_cast<MML::Real>(1000.0)));
	}
	
	SECTION("Right Infinite") {
		ClosedToInfInterval a(static_cast<MML::Real>(5.0));
		
		Interval result = Interval::Complement(a);
		
		REQUIRE(result.contains(-static_cast<MML::Real>(1000.0)));       // Complement (from -inf)
		REQUIRE(result.contains(static_cast<MML::Real>(4.999)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0)));     // Part of a (closed)
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(1000.0)));  // Part of a
	}
	
	SECTION("Complete Real Line") {
		CompleteRInterval a;
		
		Interval result = Interval::Complement(a);
		
		// Complement of R is empty set
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(0.0)));
		REQUIRE_FALSE(result.contains(-static_cast<MML::Real>(1000.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(1000.0)));
	}
}

//=============================================================================
// COMPOUND INTERVAL TESTS
//=============================================================================

TEST_CASE("Compound Interval - Basic Operations", "[intervals][compound]") {
		TEST_PRECISION_INFO();
	Interval compound;
	compound.AddInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0)))
	        .AddInterval(ClosedInterval(static_cast<MML::Real>(5.0), static_cast<MML::Real>(8.0)))
	        .AddInterval(ClosedInterval(static_cast<MML::Real>(10.0), static_cast<MML::Real>(12.0)));
	
	SECTION("Bounds") {
		REQUIRE_THAT(compound.getLowerBound(), RealApprox(static_cast<MML::Real>(0.0)));
		REQUIRE_THAT(compound.getUpperBound(), RealApprox(static_cast<MML::Real>(12.0)));
		REQUIRE_THAT(compound.getLength(), RealApprox(static_cast<MML::Real>(12.0)));  // Hull: 12 - 0 = 12
		REQUIRE_THAT(compound.getMeasure(), RealApprox(static_cast<MML::Real>(8.0)));  // Measure: 3 + 3 + 2 = 8
	}
	
	SECTION("Contains") {
		REQUIRE(compound.contains(static_cast<MML::Real>(1.0)));    // In first interval
		REQUIRE(compound.contains(static_cast<MML::Real>(6.0)));    // In second interval
		REQUIRE(compound.contains(static_cast<MML::Real>(11.0)));   // In third interval
		REQUIRE_FALSE(compound.contains(static_cast<MML::Real>(4.0)));  // In gap
		REQUIRE_FALSE(compound.contains(static_cast<MML::Real>(9.0)));  // In gap
		REQUIRE_FALSE(compound.contains(-static_cast<MML::Real>(1.0))); // Before all
		REQUIRE_FALSE(compound.contains(static_cast<MML::Real>(13.0))); // After all
	}
	
	SECTION("Interval Contains") {
		REQUIRE(compound.contains(ClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0))));     // In first
		REQUIRE(compound.contains(ClosedInterval(static_cast<MML::Real>(10.5), static_cast<MML::Real>(11.5))));   // In third
		REQUIRE_FALSE(compound.contains(ClosedInterval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(6.0)))); // Spans gap
	}
	
	SECTION("Interval Intersects") {
		REQUIRE(compound.intersects(ClosedInterval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0))));   // Overlaps first
		REQUIRE(compound.intersects(ClosedInterval(static_cast<MML::Real>(7.0), static_cast<MML::Real>(9.0))));   // Overlaps second
		REQUIRE(compound.intersects(ClosedInterval(-static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0))));  // Overlaps first
		REQUIRE_FALSE(compound.intersects(ClosedInterval(static_cast<MML::Real>(3.5), static_cast<MML::Real>(4.5)))); // In gap
		REQUIRE_FALSE(compound.intersects(ClosedInterval(static_cast<MML::Real>(20.0), static_cast<MML::Real>(25.0)))); // After all
	}
}

TEST_CASE("Compound Interval - Equidistant Covering", "[intervals][compound][covering]") {
		TEST_PRECISION_INFO();
	Interval compound;
	compound.AddInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0)))   // Length 3
	        .AddInterval(ClosedInterval(static_cast<MML::Real>(10.0), static_cast<MML::Real>(14.0))); // Length 4
	
	std::vector<Real> points;
	compound.getEquidistantCovering(14, points);
	
	SECTION("Points Distributed Proportionally") {
		REQUIRE(points.size() >= 5); // Should have points
		
		// Count points in each interval
		int count1 = 0, count2 = 0;
		for (Real p : points) {
			if (p >= static_cast<MML::Real>(0.0) && p <= static_cast<MML::Real>(3.0)) count1++;
			if (p >= static_cast<MML::Real>(10.0) && p <= static_cast<MML::Real>(14.0)) count2++;
		}
		
		// Roughly proportional to lengths (3:4 ratio)
		REQUIRE(count1 > 0);
		REQUIRE(count2 > 0);
		// count2/count1 should be roughly 4/3 � static_cast<MML::Real>(1.33)
	}
}

//=============================================================================
// EDGE CASES AND SPECIAL SCENARIOS
//=============================================================================

TEST_CASE("Edge Cases - Degenerate Intervals", "[intervals][edge]") {
		TEST_PRECISION_INFO();
	SECTION("Zero-Length Interval") {
		ClosedInterval point(static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0));
		
		REQUIRE(point.getLowerBound() == static_cast<MML::Real>(5.0));
		REQUIRE(point.getUpperBound() == static_cast<MML::Real>(5.0));
		REQUIRE(point.getLength() == static_cast<MML::Real>(0.0));
		REQUIRE(point.contains(static_cast<MML::Real>(5.0)));
		REQUIRE_FALSE(point.contains(static_cast<MML::Real>(4.999)));
		REQUIRE_FALSE(point.contains(static_cast<MML::Real>(5.001)));
	}
	
	SECTION("Negative Interval (Invalid)") {
		// Reversed bounds are rejected at construction (MinimalMathLibrary-wewi.2.2)
		REQUIRE_THROWS_AS(ClosedInterval(static_cast<MML::Real>(10.0), static_cast<MML::Real>(5.0)), ArgumentError);
	}
}

TEST_CASE("Edge Cases - Precision", "[intervals][edge][precision]") {
		TEST_PRECISION_INFO();
	SECTION("Very Small Interval") {
		ClosedInterval tiny(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5));
		
		REQUIRE_THAT(tiny.getLength(), RealApprox(TOL(1e-10, 1e-5)));
		REQUIRE(tiny.contains(TOL(5e-11, 5e-6)));
		REQUIRE_FALSE(tiny.contains(TOL(2e-10, 2e-5)));
	}
	
	SECTION("Very Large Interval") {
		ClosedInterval huge(-1e100, 1e100);
		
		REQUIRE(huge.contains(static_cast<MML::Real>(0.0)));
		REQUIRE(huge.contains(-5e99));
		REQUIRE(huge.contains(5e99));
	}
}

TEST_CASE("Complex Scenarios - Chained Operations", "[intervals][complex]") {
		TEST_PRECISION_INFO();
	SECTION("Intersection of Differences") {
		ClosedInterval a(static_cast<MML::Real>(0.0), static_cast<MML::Real>(10.0));
		ClosedInterval b(static_cast<MML::Real>(3.0), static_cast<MML::Real>(7.0));
		ClosedInterval c(static_cast<MML::Real>(5.0), static_cast<MML::Real>(12.0));

		// (A \ B) ∩ C, now expressible directly via the closed algebra.
		// A\B = [0,3) ∪ (7,10];  ∩ [5,12] = (7,10]
		CompositeInterval result = (CompositeInterval(a) - b) & c;

		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(1.0)));   // in [0,3) but outside C
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(5.0)));   // removed by B
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(7.0)));   // removed by B (closed at 7)
		REQUIRE(result.contains(static_cast<MML::Real>(8.0)));         // (7,10] ∩ C
		REQUIRE(result.contains(static_cast<MML::Real>(10.0)));
		REQUIRE_FALSE(result.contains(static_cast<MML::Real>(11.0)));  // beyond A
	}
	
	SECTION("Multiple Complements") {
		ClosedInterval a(static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0));
		
		Interval comp1 = Interval::Complement(a);     // (-8, 0) ? (5, 8)
		
		REQUIRE(comp1.contains(-static_cast<MML::Real>(10.0)));
		REQUIRE(comp1.contains(static_cast<MML::Real>(10.0)));
		REQUIRE_FALSE(comp1.contains(static_cast<MML::Real>(2.5)));
	}
}

TEST_CASE("Real-World Use Cases", "[intervals][practical]") {
		TEST_PRECISION_INFO();
	SECTION("Function Domain - Tangent") {
		// tan(x) undefined at x = p/2 + np
		// Domain: R \ {p/2 + np}
		ClosedIntervalWithReccuringPointHoles tanDomain(
			-2 * Constants::PI, 
			2 * Constants::PI,
			static_cast<MML::Real>(0.5) * Constants::PI,  // First hole
			Constants::PI         // Hole spacing
		);
		
		REQUIRE(tanDomain.contains(static_cast<MML::Real>(0.0)));
		REQUIRE(tanDomain.contains(Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(static_cast<MML::Real>(0.5) * Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(-static_cast<MML::Real>(0.5) * Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(static_cast<MML::Real>(1.5) * Constants::PI));
	}
	
	SECTION("Temperature Range with Exclusions") {
		Interval validTemps;
		validTemps.AddInterval(ClosedInterval(-static_cast<MML::Real>(50.0), static_cast<MML::Real>(0.0)))    // Freezing range
		          .AddInterval(OpenInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(100.0)))      // Liquid range (excluding phase boundaries)
		          .AddInterval(ClosedInterval(static_cast<MML::Real>(100.0), static_cast<MML::Real>(200.0))); // Boiling range
		
		REQUIRE(validTemps.contains(-static_cast<MML::Real>(25.0)));   // Freezing
		REQUIRE(validTemps.contains(static_cast<MML::Real>(0.0)));     // Freezing point (in first interval)
		REQUIRE(validTemps.contains(static_cast<MML::Real>(50.0)));    // Liquid
		REQUIRE(validTemps.contains(static_cast<MML::Real>(100.0)));   // Boiling point (in third interval)
		REQUIRE(validTemps.contains(static_cast<MML::Real>(150.0)));   // Above boiling
	}
	
	SECTION("Root Finding Search Interval") {
		ClosedInterval searchSpace(-static_cast<MML::Real>(10.0), static_cast<MML::Real>(10.0));
		std::vector<Real> samplePoints;
		
		searchSpace.getEquidistantCovering(21, samplePoints);
		
		REQUIRE(samplePoints.size() == 21);
		REQUIRE_THAT(samplePoints[0], RealApprox(-static_cast<MML::Real>(10.0)));
		REQUIRE_THAT(samplePoints[10], RealApprox(static_cast<MML::Real>(0.0)));
		REQUIRE_THAT(samplePoints[20], RealApprox(static_cast<MML::Real>(10.0)));
		
		// Verify uniform spacing of static_cast<MML::Real>(1.0)
		for (size_t i = 0; i < samplePoints.size() - 1; i++) {
			REQUIRE_THAT(samplePoints[i+1] - samplePoints[i], RealApprox(static_cast<MML::Real>(1.0)));
		}
	}
}

TEST_CASE("Performance - Large Compound Intervals", "[intervals][performance]") {
		TEST_PRECISION_INFO();
	SECTION("Many Sub-Intervals") {
		Interval large;
		
		// Create 100 disjoint intervals
		for (int i = 0; i < 100; i++) {
			large.AddInterval(ClosedInterval(i * static_cast<MML::Real>(10.0), i * static_cast<MML::Real>(10.0) + static_cast<MML::Real>(5.0)));
		}
		
		REQUIRE_THAT(large.getLowerBound(), RealApprox(static_cast<MML::Real>(0.0)));
		REQUIRE_THAT(large.getUpperBound(), RealApprox(static_cast<MML::Real>(995.0)));
		REQUIRE_THAT(large.getLength(), RealApprox(static_cast<MML::Real>(995.0)));   // Hull: 995 - 0 = 995
		REQUIRE_THAT(large.getMeasure(), RealApprox(static_cast<MML::Real>(500.0)));  // Measure: 100 * 5 = 500
		
		REQUIRE(large.contains(static_cast<MML::Real>(2.5)));     // In first interval
		REQUIRE(large.contains(static_cast<MML::Real>(502.5)));   // In 51st interval
		REQUIRE_FALSE(large.contains(static_cast<MML::Real>(7.0)));   // In gap
		REQUIRE_FALSE(large.contains(static_cast<MML::Real>(1000.0))); // After all
	}
}

//=============================================================================
// CLOSED ALGEBRA (CompositeInterval members / operators)
//=============================================================================

TEST_CASE("CompositeInterval - implicit conversion and union chaining", "[intervals][closed]") {
	TEST_PRECISION_INFO();
	ClosedInterval a(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0));
	ClosedInterval b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
	ClosedInterval c(static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));

	// b and c convert implicitly to CompositeInterval through operator|
	CompositeInterval u = CompositeInterval(a) | b | c;

	REQUIRE(u.contains(static_cast<MML::Real>(0.5)));
	REQUIRE(u.contains(static_cast<MML::Real>(2.5)));
	REQUIRE(u.contains(static_cast<MML::Real>(4.5)));
	REQUIRE_FALSE(u.contains(static_cast<MML::Real>(1.5)));
	REQUIRE_FALSE(u.contains(static_cast<MML::Real>(3.5)));
	REQUIRE_THAT(u.getMeasure(), RealApprox(static_cast<MML::Real>(3.0)));  // 1 + 1 + 1
	REQUIRE_THAT(u.getLength(),  RealApprox(static_cast<MML::Real>(5.0)));  // hull 0..5
}

TEST_CASE("CompositeInterval - Normalize merges overlapping and adjacent", "[intervals][closed][normalize]") {
	TEST_PRECISION_INFO();

	// Overlapping [0,2] u [1,3] = [0,3], measure 3 (not 4)
	auto overlap = CompositeInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0))) | ClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0));
	REQUIRE_THAT(overlap.getMeasure(), RealApprox(static_cast<MML::Real>(3.0)));
	REQUIRE(overlap.contains(static_cast<MML::Real>(2.5)));

	// Adjacent closed [0,1] u [1,2] merges (shared point closed) => [0,2]
	auto adjClosed = CompositeInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0))) | ClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
	REQUIRE_THAT(adjClosed.getMeasure(), RealApprox(static_cast<MML::Real>(2.0)));
	REQUIRE(adjClosed.contains(static_cast<MML::Real>(1.0)));

	// Adjacent open [0,1) u (1,2] do NOT merge: shared point 1 excluded by both
	auto adjOpen = CompositeInterval(ClosedOpenInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0))) | OpenClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
	REQUIRE_FALSE(adjOpen.contains(static_cast<MML::Real>(1.0)));
	REQUIRE(adjOpen.contains(static_cast<MML::Real>(0.5)));
	REQUIRE(adjOpen.contains(static_cast<MML::Real>(1.5)));
}

TEST_CASE("CompositeInterval - chained (a|b|c) & d", "[intervals][closed][intersection]") {
	TEST_PRECISION_INFO();
	auto s = CompositeInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0)))
	       | ClosedInterval(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0))
	       | ClosedInterval(static_cast<MML::Real>(8.0), static_cast<MML::Real>(10.0));

	auto r = s & ClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(9.0));  // [1,2] u [4,6] u [8,9]

	REQUIRE(r.contains(static_cast<MML::Real>(1.5)));
	REQUIRE(r.contains(static_cast<MML::Real>(5.0)));
	REQUIRE(r.contains(static_cast<MML::Real>(8.5)));
	REQUIRE_FALSE(r.contains(static_cast<MML::Real>(0.5)));
	REQUIRE_FALSE(r.contains(static_cast<MML::Real>(3.0)));
	REQUIRE_FALSE(r.contains(static_cast<MML::Real>(9.5)));
	REQUIRE_THAT(r.getMeasure(), RealApprox(static_cast<MML::Real>(4.0)));  // 1 + 2 + 1
}

TEST_CASE("CompositeInterval - complement of a composite and double complement", "[intervals][closed][complement]") {
	TEST_PRECISION_INFO();
	auto s = CompositeInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0))) | ClosedInterval(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));

	auto comp = ~s;  // (-inf,0) u (1,3) u (4,+inf)
	REQUIRE(comp.contains(-static_cast<MML::Real>(5.0)));
	REQUIRE(comp.contains(static_cast<MML::Real>(2.0)));
	REQUIRE(comp.contains(static_cast<MML::Real>(10.0)));
	REQUIRE_FALSE(comp.contains(static_cast<MML::Real>(0.0)));   // 0 belongs to s (closed)
	REQUIRE_FALSE(comp.contains(static_cast<MML::Real>(0.5)));
	REQUIRE_FALSE(comp.contains(static_cast<MML::Real>(3.5)));

	// Double complement recovers the original set (containment-wise)
	auto cc = ~comp;
	REQUIRE(cc.contains(static_cast<MML::Real>(0.0)));
	REQUIRE(cc.contains(static_cast<MML::Real>(0.5)));
	REQUIRE(cc.contains(static_cast<MML::Real>(3.5)));
	REQUIRE_FALSE(cc.contains(static_cast<MML::Real>(2.0)));
}

TEST_CASE("CompositeInterval - difference splits and endpoint types", "[intervals][closed][difference]") {
	TEST_PRECISION_INFO();
	auto r = CompositeInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0))) - ClosedInterval(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)); // [0,1) u (2,3]

	REQUIRE(r.contains(static_cast<MML::Real>(0.0)));
	REQUIRE(r.contains(static_cast<MML::Real>(0.5)));
	REQUIRE_FALSE(r.contains(static_cast<MML::Real>(1.0)));   // removed (b closed at 1)
	REQUIRE_FALSE(r.contains(static_cast<MML::Real>(1.5)));
	REQUIRE_FALSE(r.contains(static_cast<MML::Real>(2.0)));   // removed (b closed at 2)
	REQUIRE(r.contains(static_cast<MML::Real>(2.5)));
	REQUIRE(r.contains(static_cast<MML::Real>(3.0)));
	REQUIRE_THAT(r.getMeasure(), RealApprox(static_cast<MML::Real>(2.0)));  // [0,1) + (2,3]
}

TEST_CASE("CompositeInterval - static overloads accept sets and stay back-compatible", "[intervals][closed]") {
	TEST_PRECISION_INFO();

	// New: static ops accept composite sets
	auto u = Interval::Union(CompositeInterval(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0))),
	                         CompositeInterval(ClosedInterval(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0))));
	REQUIRE(u.contains(static_cast<MML::Real>(0.5)));
	REQUIRE(u.contains(static_cast<MML::Real>(2.5)));
	REQUIRE_FALSE(u.contains(static_cast<MML::Real>(1.5)));

	// Complement of an open set via static overload: (-inf,0] u [1,+inf)
	auto comp = Interval::Complement(CompositeInterval(OpenInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0))));
	REQUIRE(comp.contains(static_cast<MML::Real>(0.0)));
	REQUIRE(comp.contains(static_cast<MML::Real>(1.0)));
	REQUIRE_FALSE(comp.contains(static_cast<MML::Real>(0.5)));

	// Back-compat: two plain BaseIntervals still bind to the original static op
	auto ob = Interval::Union(ClosedInterval(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)), ClosedInterval(static_cast<MML::Real>(0.5), static_cast<MML::Real>(2.0)));
	REQUIRE(ob.contains(static_cast<MML::Real>(1.5)));
	REQUIRE(ob.contains(static_cast<MML::Real>(0.0)));
}

} // namespace MML::Tests::Base::IntervalsTests
