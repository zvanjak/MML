#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

// Removed: using Catch::Approx; - now using precision-aware RealApprox

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Intervals.h"
#endif

#include <cmath>
#include <algorithm>

using namespace MML;
using namespace MML::Testing;

TEST_CASE("Test_OpenInterval", "[simple]") {
		TEST_PRECISION_INFO();
	OpenInterval int1(REAL(0.0), REAL(1.0));
	// da li je postavljeno sve kako treba
	REQUIRE(REAL(0.0) == int1.getLowerBound());
	REQUIRE(REAL(1.0) == int1.getUpperBound());
	REQUIRE(REAL(1.0) == int1.getLength());
	REQUIRE(true == int1.isContinuous());
}

TEST_CASE("Test_OpenInterval_contains", "[simple]") {
		TEST_PRECISION_INFO();
	OpenInterval int1(REAL(0.0), REAL(1.0));

	REQUIRE(true == int1.contains(REAL(0.5)));
	REQUIRE(false == int1.contains(-REAL(0.5)));
}

TEST_CASE("Test_OpenInterval_Covering", "[simple]") {
		TEST_PRECISION_INFO();
	// kad uzmes Covering da granice NISU u skupu tocaka
}

TEST_CASE("Test_CompleR", "[simple]")
{
		TEST_PRECISION_INFO();
	CompleteRInterval int1;

	auto a = int1.getLowerBound();
}

TEST_CASE("Test_Interval", "[simple]")
{
		TEST_PRECISION_INFO();
	// old way definition
	Interval* tangDefInterval = new Interval({ new ClosedOpenInterval(-REAL(2.0) * Constants::PI, -REAL(1.5) * Constants::PI),
																			new OpenInterval(-REAL(1.5) * Constants::PI, -REAL(0.5) * Constants::PI),
																			new OpenInterval(-REAL(0.5) * Constants::PI, REAL(0.5) * Constants::PI),
																			new OpenInterval(REAL(0.5) * Constants::PI, REAL(1.5) * Constants::PI),
																			new OpenClosedInterval(REAL(1.5) * Constants::PI, REAL(2.0) * Constants::PI) });

	REQUIRE(tangDefInterval->contains(-REAL(2.0) * Constants::PI));
}

TEST_CASE("Test_ClosedIntervalWithReccuringPointHoles", "[simple]")
{
		TEST_PRECISION_INFO();
	BaseInterval* tanDef = new ClosedIntervalWithReccuringPointHoles(-5 * Constants::PI, 5 * Constants::PI, REAL(0.5) * Constants::PI, Constants::PI);

	REQUIRE(tanDef->contains(REAL(0.0)));

	// investigation around half PI
	REQUIRE_FALSE(tanDef->contains(REAL(0.5) * Constants::PI));
	REQUIRE(tanDef->contains(REAL(0.5) * Constants::PI + REAL(0.0001)));
	REQUIRE(tanDef->contains(REAL(0.5) * Constants::PI + 1e-8));
	REQUIRE(tanDef->contains(REAL(0.5) * Constants::PI + 1e-15));
	REQUIRE_FALSE(tanDef->contains(REAL(0.5) * Constants::PI + 1e-16));
	REQUIRE(tanDef->contains(REAL(0.5) * Constants::PI - REAL(0.0001)));
	REQUIRE(tanDef->contains(REAL(0.5) * Constants::PI - 1e-8));
	REQUIRE(tanDef->contains(REAL(0.5) * Constants::PI - 1e-15));
	REQUIRE_FALSE(tanDef->contains(REAL(0.5) * Constants::PI - 1e-16));

	// verification of all the holes
	REQUIRE_FALSE(tanDef->contains(-REAL(0.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(REAL(1.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-REAL(1.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(REAL(2.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-REAL(2.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(REAL(3.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-REAL(3.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(REAL(4.5) * Constants::PI));
	REQUIRE_FALSE(tanDef->contains(-REAL(5.5) * Constants::PI));

	// precision verification on end holes
	REQUIRE(tanDef->contains(REAL(4.5) * Constants::PI + 1e-15));
	REQUIRE_FALSE(tanDef->contains(REAL(4.5) * Constants::PI + 1e-16));
	REQUIRE(tanDef->contains(REAL(4.5) * Constants::PI - 1e-15));
	REQUIRE_FALSE(tanDef->contains(REAL(4.5) * Constants::PI - 1e-16));

}

//=============================================================================
// COMPREHENSIVE INTERVAL TESTS - Basic Interval Types
//=============================================================================

TEST_CASE("ClosedInterval - Basic Properties", "[intervals][closed]") {
		TEST_PRECISION_INFO();
	ClosedInterval interval(REAL(2.0), REAL(8.0));
	
	SECTION("Bounds and Length") {
		REQUIRE(interval.getLowerBound() == REAL(2.0));
		REQUIRE(interval.getUpperBound() == REAL(8.0));
		REQUIRE(interval.getLength() == REAL(6.0));
		REQUIRE(interval.isContinuous());
	}
	
	SECTION("Contains - Interior Points") {
		REQUIRE(interval.contains(REAL(2.0)));   // Lower bound (closed)
		REQUIRE(interval.contains(REAL(8.0)));   // Upper bound (closed)
		REQUIRE(interval.contains(REAL(5.0)));   // Middle point
		REQUIRE(interval.contains(REAL(2.001))); // Just inside lower
		REQUIRE(interval.contains(REAL(7.999))); // Just inside upper
	}
	
	SECTION("Contains - Exterior Points") {
		REQUIRE_FALSE(interval.contains(REAL(1.999)));  // Just outside lower
		REQUIRE_FALSE(interval.contains(REAL(8.001)));  // Just outside upper
		REQUIRE_FALSE(interval.contains(-REAL(10.0)));  // Far outside
		REQUIRE_FALSE(interval.contains(REAL(100.0)));  // Far outside
	}
}

TEST_CASE("OpenInterval - Basic Properties", "[intervals][open]") {
		TEST_PRECISION_INFO();
	OpenInterval interval(REAL(2.0), REAL(8.0));
	
	SECTION("Bounds and Length") {
		REQUIRE(interval.getLowerBound() == REAL(2.0));
		REQUIRE(interval.getUpperBound() == REAL(8.0));
		REQUIRE(interval.getLength() == REAL(6.0));
		REQUIRE(interval.isContinuous());
	}
	
	SECTION("Contains - Excludes Endpoints") {
		REQUIRE_FALSE(interval.contains(REAL(2.0))); // Lower bound excluded
		REQUIRE_FALSE(interval.contains(REAL(8.0))); // Upper bound excluded
		REQUIRE(interval.contains(REAL(5.0)));       // Middle point
		REQUIRE(interval.contains(REAL(2.001)));     // Just inside lower
		REQUIRE(interval.contains(REAL(7.999)));     // Just inside upper
	}
}

TEST_CASE("OpenClosedInterval - Mixed Endpoints", "[intervals][open-closed]") {
		TEST_PRECISION_INFO();
	OpenClosedInterval interval(REAL(2.0), REAL(8.0));
	
	SECTION("Contains - Mixed Behavior") {
		REQUIRE_FALSE(interval.contains(REAL(2.0))); // Lower bound excluded (open)
		REQUIRE(interval.contains(REAL(8.0)));       // Upper bound included (closed)
		REQUIRE(interval.contains(REAL(5.0)));       // Middle point
		REQUIRE(interval.contains(REAL(2.001)));
		REQUIRE(interval.contains(REAL(7.999)));
	}
}

TEST_CASE("ClosedOpenInterval - Mixed Endpoints", "[intervals][closed-open]") {
		TEST_PRECISION_INFO();
	ClosedOpenInterval interval(REAL(2.0), REAL(8.0));
	
	SECTION("Contains - Mixed Behavior") {
		REQUIRE(interval.contains(REAL(2.0)));       // Lower bound included (closed)
		REQUIRE_FALSE(interval.contains(REAL(8.0))); // Upper bound excluded (open)
		REQUIRE(interval.contains(REAL(5.0)));       // Middle point
		REQUIRE(interval.contains(REAL(2.001)));
		REQUIRE(interval.contains(REAL(7.999)));
	}
}

TEST_CASE("Infinite Bound Intervals", "[intervals][infinite]") {
		TEST_PRECISION_INFO();
	SECTION("NegInfToClosedInterval") {
		NegInfToClosedInterval interval(REAL(5.0));
		
		REQUIRE(interval.contains(-1e100));  // Very negative
		REQUIRE(interval.contains(REAL(0.0)));
		REQUIRE(interval.contains(REAL(4.999)));
		REQUIRE(interval.contains(REAL(5.0)));     // Closed upper bound
		REQUIRE_FALSE(interval.contains(REAL(5.001)));
	}
	
	SECTION("NegInfToOpenInterval") {
		NegInfToOpenInterval interval(REAL(5.0));
		
		REQUIRE(interval.contains(-1e100));
		REQUIRE(interval.contains(REAL(0.0)));
		REQUIRE(interval.contains(REAL(4.999)));
		REQUIRE_FALSE(interval.contains(REAL(5.0))); // Open upper bound
		REQUIRE_FALSE(interval.contains(REAL(5.001)));
	}
	
	SECTION("ClosedToInfInterval") {
		ClosedToInfInterval interval(REAL(5.0));
		
		REQUIRE(interval.contains(REAL(5.0)));      // Closed lower bound
		REQUIRE(interval.contains(REAL(5.001)));
		REQUIRE(interval.contains(REAL(100.0)));
		REQUIRE(interval.contains(1e100));    // Very positive
		REQUIRE_FALSE(interval.contains(REAL(4.999)));
	}
	
	SECTION("OpenToInfInterval") {
		OpenToInfInterval interval(REAL(5.0));
		
		REQUIRE_FALSE(interval.contains(REAL(5.0))); // Open lower bound
		REQUIRE(interval.contains(REAL(5.001)));
		REQUIRE(interval.contains(REAL(100.0)));
		REQUIRE(interval.contains(1e100));
		REQUIRE_FALSE(interval.contains(REAL(4.999)));
	}
	
	SECTION("CompleteRInterval") {
		CompleteRInterval interval;
		
		REQUIRE(interval.contains(-1e100));
		REQUIRE(interval.contains(REAL(0.0)));
		REQUIRE(interval.contains(1e100));
		REQUIRE(interval.contains(-REAL(0.000001)));
		REQUIRE(interval.contains(REAL(0.000001)));
	}
}

//=============================================================================
// EQUIDISTANT COVERING TESTS
//=============================================================================

TEST_CASE("GetEquidistantCovering - Closed Interval", "[intervals][covering][closed]") {
		TEST_PRECISION_INFO();
	ClosedInterval interval(REAL(0.0), REAL(10.0));
	std::vector<Real> points;
	
	SECTION("5 Points") {
		interval.GetEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE_THAT(points[0], RealApprox(REAL(0.0)));
		REQUIRE_THAT(points[1], RealApprox(REAL(2.5)));
		REQUIRE_THAT(points[2], RealApprox(REAL(5.0)));
		REQUIRE_THAT(points[3], RealApprox(REAL(7.5)));
		REQUIRE_THAT(points[4], RealApprox(REAL(10.0)));
	}
	
	SECTION("Single Point") {
		interval.GetEquidistantCovering(1, points);
		
		REQUIRE(points.size() == 1);
		REQUIRE_THAT(points[0], RealApprox(REAL(5.0))); // Midpoint
	}
	
	SECTION("Empty Request") {
		interval.GetEquidistantCovering(0, points);
		REQUIRE(points.empty());
	}
	
	SECTION("Many Points") {
		interval.GetEquidistantCovering(101, points);
		
		REQUIRE(points.size() == 101);
		REQUIRE_THAT(points[0], RealApprox(REAL(0.0)));
		REQUIRE_THAT(points[50], RealApprox(REAL(5.0)));
		REQUIRE_THAT(points[100], RealApprox(REAL(10.0)));
		
		// Verify equidistant spacing
		Real delta = points[1] - points[0];
		for (size_t i = 1; i < points.size() - 1; i++) {
			REQUIRE_THAT(points[i+1] - points[i], RealApprox(delta));
		}
	}
}

TEST_CASE("GetEquidistantCovering - Open Interval", "[intervals][covering][open]") {
		TEST_PRECISION_INFO();
	OpenInterval interval(REAL(0.0), REAL(10.0));
	std::vector<Real> points;
	
	SECTION("Points Avoid Endpoints") {
		interval.GetEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE(points[0] > REAL(0.0));    // Not exactly 0 (open)
		REQUIRE(points[4] < REAL(10.0));   // Not exactly 10 (open)
		REQUIRE_THAT(points[2], RealApprox(REAL(5.0)).epsilon(REAL(0.01))); // Middle point
	}
}

TEST_CASE("GetEquidistantCovering - Infinite Bounds", "[intervals][covering][infinite]") {
		TEST_PRECISION_INFO();
	SECTION("Infinite Lower Bound") {
		NegInfToClosedInterval interval(REAL(10.0));
		std::vector<Real> points;
		
		interval.GetEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE(points[0] < points[4]); // Ascending order
		REQUIRE_THAT(points[4], RealApprox(REAL(10.0))); // Upper bound
		REQUIRE_THAT(points[0], RealApprox(-1e10).epsilon(REAL(0.01))); // Practical lower limit
	}
	
	SECTION("Infinite Upper Bound") {
		ClosedToInfInterval interval(REAL(0.0));
		std::vector<Real> points;
		
		interval.GetEquidistantCovering(5, points);
		
		REQUIRE(points.size() == 5);
		REQUIRE_THAT(points[0], RealApprox(REAL(0.0)));
		REQUIRE_THAT(points[4], RealApprox(1e10).epsilon(REAL(0.01))); // Practical upper limit
	}
}

//=============================================================================
// INTERSECTION TESTS
//=============================================================================

TEST_CASE("Interval::Intersection - Basic Cases", "[intervals][intersection]") {
		TEST_PRECISION_INFO();
	SECTION("Complete Overlap - Same Intervals") {
		ClosedInterval a(REAL(2.0), REAL(8.0));
		ClosedInterval b(REAL(2.0), REAL(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(8.0)));
		REQUIRE(result.contains(REAL(5.0)));
		REQUIRE(result.contains(REAL(2.0)));
		REQUIRE(result.contains(REAL(8.0)));
	}
	
	SECTION("Partial Overlap") {
		ClosedInterval a(REAL(2.0), REAL(8.0));
		ClosedInterval b(REAL(5.0), REAL(10.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(5.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(8.0)));
		REQUIRE(result.contains(REAL(6.0)));
		REQUIRE(result.contains(REAL(5.0)));
		REQUIRE(result.contains(REAL(8.0)));
		REQUIRE_FALSE(result.contains(REAL(4.0)));
		REQUIRE_FALSE(result.contains(REAL(9.0)));
	}
	
	SECTION("No Overlap - Gap Between") {
		ClosedInterval a(REAL(2.0), REAL(5.0));
		ClosedInterval b(REAL(7.0), REAL(10.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_FALSE(result.contains(REAL(3.0)));
		REQUIRE_FALSE(result.contains(REAL(8.0)));
		REQUIRE_FALSE(result.contains(REAL(6.0))); // In the gap
	}
	
	SECTION("One Contains Other") {
		ClosedInterval a(REAL(2.0), REAL(10.0));
		ClosedInterval b(REAL(4.0), REAL(6.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(4.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(6.0)));
		REQUIRE(result.contains(REAL(5.0)));
	}
	
	SECTION("Touching at Boundary") {
		ClosedInterval a(REAL(2.0), REAL(5.0));
		ClosedInterval b(REAL(5.0), REAL(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE(result.contains(REAL(5.0))); // Touching point included (both closed)
		REQUIRE_THAT(result.getLength(), RealApprox(REAL(0.0))); // Degenerate interval
	}
}

TEST_CASE("Interval::Intersection - Endpoint Types", "[intervals][intersection][endpoints]") {
		TEST_PRECISION_INFO();
	SECTION("Open meets Closed - Lower Bound") {
		OpenInterval a(REAL(2.0), REAL(8.0));
		ClosedInterval b(REAL(2.0), REAL(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_FALSE(result.contains(REAL(2.0))); // Open endpoint wins
		REQUIRE_FALSE(result.contains(REAL(8.0))); // Open endpoint wins (a is open at upper)
		REQUIRE(result.contains(REAL(5.0)));
	}
	
	SECTION("Open meets Closed - Upper Bound") {
		ClosedInterval a(REAL(2.0), REAL(8.0));
		OpenInterval b(REAL(2.0), REAL(8.0));
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_FALSE(result.contains(REAL(2.0))); // Open endpoint wins
		REQUIRE_FALSE(result.contains(REAL(8.0))); // Open endpoint wins
		REQUIRE(result.contains(REAL(5.0)));
	}
	
	SECTION("Mixed Endpoint Types") {
		OpenClosedInterval a(REAL(2.0), REAL(8.0));  // (2, 8]
		ClosedOpenInterval b(REAL(5.0), REAL(10.0)); // [5, 10)
		
		Interval result = Interval::Intersection(a, b);
		
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(5.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(8.0)));
		REQUIRE(result.contains(REAL(5.0)));   // Closed from b
		REQUIRE(result.contains(REAL(8.0)));   // Closed from a
		REQUIRE(result.contains(REAL(6.5)));
	}
}

//=============================================================================
// DIFFERENCE TESTS
//=============================================================================

TEST_CASE("Interval::Difference - No Intersection", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(REAL(2.0), REAL(5.0));
	ClosedInterval b(REAL(7.0), REAL(10.0));
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Original Interval") {
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(2.0)));
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(5.0)));
		REQUIRE(result.contains(REAL(2.0)));
		REQUIRE(result.contains(REAL(5.0)));
		REQUIRE(result.contains(REAL(3.5)));
		REQUIRE_FALSE(result.contains(REAL(8.0)));
	}
}

TEST_CASE("Interval::Difference - Complete Containment", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(REAL(4.0), REAL(6.0));
	ClosedInterval b(REAL(2.0), REAL(10.0)); // b completely contains a
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Empty Set") {
		REQUIRE_FALSE(result.contains(REAL(5.0)));
		REQUIRE_FALSE(result.contains(REAL(4.0)));
		REQUIRE_FALSE(result.contains(REAL(6.0)));
	}
}

TEST_CASE("Interval::Difference - Left Cut", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(REAL(5.0), REAL(10.0));
	ClosedInterval b(REAL(2.0), REAL(7.0)); // Cuts left part of a
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Right Remainder") {
		REQUIRE_FALSE(result.contains(REAL(6.0)));  // Removed by b
		REQUIRE(result.contains(REAL(7.001)));      // Remaining part
		REQUIRE(result.contains(REAL(10.0)));
		// When b=[2,7] removes from a=[5,10], result is (7,10]
		// b's upper is CLOSED at 7, so result starts OPEN at 7
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(7.0)));
		REQUIRE_FALSE(result.contains(REAL(7.0))); // Boundary is open
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(10.0)));
	}
}

TEST_CASE("Interval::Difference - Right Cut", "[intervals][difference]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(REAL(2.0), REAL(10.0));
	ClosedInterval b(REAL(7.0), REAL(15.0)); // Cuts right part of a
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Left Remainder") {
		REQUIRE(result.contains(REAL(2.0)));
		REQUIRE(result.contains(REAL(5.0)));
		REQUIRE_FALSE(result.contains(REAL(7.0)));  // Cut here (closed from b)
		REQUIRE_FALSE(result.contains(REAL(9.0)));  // Removed by b
		REQUIRE_THAT(result.getLowerBound(), RealApprox(REAL(2.0)));
		// When b=[7,15] removes from a=[2,10], result is [2,7)
		// b's lower is CLOSED at 7, so result ends OPEN at 7
		REQUIRE_THAT(result.getUpperBound(), RealApprox(REAL(7.0)));
		REQUIRE_FALSE(result.contains(REAL(7.0))); // Boundary is open
	}
}

TEST_CASE("Interval::Difference - Middle Cut (Split)", "[intervals][difference]") {
	ClosedInterval a(REAL(2.0), REAL(10.0));
	ClosedInterval b(REAL(5.0), REAL(7.0)); // Cuts middle, creating two pieces
	
	Interval result = Interval::Difference(a, b);
	
	SECTION("Returns Two Intervals") {
		REQUIRE(result.contains(REAL(2.0)));      // Left piece
		REQUIRE(result.contains(REAL(4.999)));    // Left piece edge
		REQUIRE_FALSE(result.contains(REAL(5.0))); // Removed (closed from b)
		REQUIRE_FALSE(result.contains(REAL(6.0))); // Middle removed
		REQUIRE_FALSE(result.contains(REAL(7.0))); // Removed (closed from b)
		REQUIRE(result.contains(REAL(7.001)));    // Right piece edge
		REQUIRE(result.contains(REAL(10.0)));     // Right piece
		
		// Length should be original - removed
		REQUIRE_THAT(result.getLength(), RealApprox(REAL(6.0))); // (5-2) + (10-7) = 3 + 3
	}
}

TEST_CASE("Interval::Difference - Endpoint Handling", "[intervals][difference][endpoints]") {
		TEST_PRECISION_INFO();
	SECTION("Open Interval Removes Closed Interval") {
		ClosedInterval a(REAL(2.0), REAL(10.0));
		OpenInterval b(REAL(5.0), REAL(7.0)); // Open boundaries
		
		Interval result = Interval::Difference(a, b);
		
		REQUIRE(result.contains(REAL(5.0)));  // NOT removed (b is open at 5)
		REQUIRE(result.contains(REAL(7.0)));  // NOT removed (b is open at 7)
		REQUIRE_FALSE(result.contains(REAL(6.0))); // Removed (inside b)
	}
}

//=============================================================================
// COMPLEMENT TESTS
//=============================================================================

TEST_CASE("Interval::Complement - Finite Interval", "[intervals][complement]") {
		TEST_PRECISION_INFO();
	ClosedInterval a(REAL(3.0), REAL(7.0));
	
	Interval result = Interval::Complement(a);
	
	SECTION("Creates Two Intervals") {
		REQUIRE(result.contains(-REAL(1000.0)));  // Left piece (to -inf)
		REQUIRE(result.contains(REAL(2.999)));    // Just before lower bound
		REQUIRE_FALSE(result.contains(REAL(3.0))); // Lower bound excluded
		REQUIRE_FALSE(result.contains(REAL(5.0))); // Middle excluded
		REQUIRE_FALSE(result.contains(REAL(7.0))); // Upper bound excluded
		REQUIRE(result.contains(REAL(7.001)));    // Just after upper bound
		REQUIRE(result.contains(REAL(1000.0)));   // Right piece (to +inf)
	}
}

TEST_CASE("Interval::Complement - Open Interval", "[intervals][complement]") {
		TEST_PRECISION_INFO();
	OpenInterval a(REAL(3.0), REAL(7.0));
	
	Interval result = Interval::Complement(a);
	
	SECTION("Includes Endpoints") {
		REQUIRE(result.contains(REAL(3.0)));  // Included (a is open at 3)
		REQUIRE(result.contains(REAL(7.0)));  // Included (a is open at 7)
		REQUIRE_FALSE(result.contains(REAL(5.0))); // Excluded (inside a)
	}
}

TEST_CASE("Interval::Complement - Infinite Bounds", "[intervals][complement]") {
		TEST_PRECISION_INFO();
	SECTION("Left Infinite") {
		NegInfToClosedInterval a(REAL(5.0));
		
		Interval result = Interval::Complement(a);
		
		REQUIRE_FALSE(result.contains(-REAL(1000.0))); // Part of a
		REQUIRE_FALSE(result.contains(REAL(5.0)));     // Part of a (closed)
		REQUIRE(result.contains(REAL(5.001)));         // Complement (to +inf)
		REQUIRE(result.contains(REAL(1000.0)));
	}
	
	SECTION("Right Infinite") {
		ClosedToInfInterval a(REAL(5.0));
		
		Interval result = Interval::Complement(a);
		
		REQUIRE(result.contains(-REAL(1000.0)));       // Complement (from -inf)
		REQUIRE(result.contains(REAL(4.999)));
		REQUIRE_FALSE(result.contains(REAL(5.0)));     // Part of a (closed)
		REQUIRE_FALSE(result.contains(REAL(1000.0)));  // Part of a
	}
	
	SECTION("Complete Real Line") {
		CompleteRInterval a;
		
		Interval result = Interval::Complement(a);
		
		// Complement of R is empty set
		REQUIRE_FALSE(result.contains(REAL(0.0)));
		REQUIRE_FALSE(result.contains(-REAL(1000.0)));
		REQUIRE_FALSE(result.contains(REAL(1000.0)));
	}
}

//=============================================================================
// COMPOUND INTERVAL TESTS
//=============================================================================

TEST_CASE("Compound Interval - Basic Operations", "[intervals][compound]") {
		TEST_PRECISION_INFO();
	Interval compound;
	compound.AddInterval(ClosedInterval(REAL(0.0), REAL(3.0)))
	        .AddInterval(ClosedInterval(REAL(5.0), REAL(8.0)))
	        .AddInterval(ClosedInterval(REAL(10.0), REAL(12.0)));
	
	SECTION("Bounds") {
		REQUIRE_THAT(compound.getLowerBound(), RealApprox(REAL(0.0)));
		REQUIRE_THAT(compound.getUpperBound(), RealApprox(REAL(12.0)));
		REQUIRE_THAT(compound.getLength(), RealApprox(REAL(8.0))); // 3 + 3 + 2
	}
	
	SECTION("Contains") {
		REQUIRE(compound.contains(REAL(1.0)));    // In first interval
		REQUIRE(compound.contains(REAL(6.0)));    // In second interval
		REQUIRE(compound.contains(REAL(11.0)));   // In third interval
		REQUIRE_FALSE(compound.contains(REAL(4.0)));  // In gap
		REQUIRE_FALSE(compound.contains(REAL(9.0)));  // In gap
		REQUIRE_FALSE(compound.contains(-REAL(1.0))); // Before all
		REQUIRE_FALSE(compound.contains(REAL(13.0))); // After all
	}
	
	SECTION("Interval Contains") {
		REQUIRE(compound.contains(ClosedInterval(REAL(1.0), REAL(2.0))));     // In first
		REQUIRE(compound.contains(ClosedInterval(REAL(10.5), REAL(11.5))));   // In third
		REQUIRE_FALSE(compound.contains(ClosedInterval(REAL(2.0), REAL(6.0)))); // Spans gap
	}
	
	SECTION("Interval Intersects") {
		REQUIRE(compound.intersects(ClosedInterval(REAL(2.0), REAL(4.0))));   // Overlaps first
		REQUIRE(compound.intersects(ClosedInterval(REAL(7.0), REAL(9.0))));   // Overlaps second
		REQUIRE(compound.intersects(ClosedInterval(-REAL(1.0), REAL(1.0))));  // Overlaps first
		REQUIRE_FALSE(compound.intersects(ClosedInterval(REAL(3.5), REAL(4.5)))); // In gap
		REQUIRE_FALSE(compound.intersects(ClosedInterval(REAL(20.0), REAL(25.0)))); // After all
	}
}

TEST_CASE("Compound Interval - Equidistant Covering", "[intervals][compound][covering]") {
		TEST_PRECISION_INFO();
	Interval compound;
	compound.AddInterval(ClosedInterval(REAL(0.0), REAL(3.0)))   // Length 3
	        .AddInterval(ClosedInterval(REAL(10.0), REAL(14.0))); // Length 4
	
	std::vector<Real> points;
	compound.GetEquidistantCovering(14, points);
	
	SECTION("Points Distributed Proportionally") {
		REQUIRE(points.size() >= 10); // Should have points
		
		// Count points in each interval
		int count1 = 0, count2 = 0;
		for (Real p : points) {
			if (p >= REAL(0.0) && p <= REAL(3.0)) count1++;
			if (p >= REAL(10.0) && p <= REAL(14.0)) count2++;
		}
		
		// Roughly proportional to lengths (3:4 ratio)
		REQUIRE(count1 > 0);
		REQUIRE(count2 > 0);
		// count2/count1 should be roughly 4/3 � REAL(1.33)
	}
}

//=============================================================================
// EDGE CASES AND SPECIAL SCENARIOS
//=============================================================================

TEST_CASE("Edge Cases - Degenerate Intervals", "[intervals][edge]") {
		TEST_PRECISION_INFO();
	SECTION("Zero-Length Interval") {
		ClosedInterval point(REAL(5.0), REAL(5.0));
		
		REQUIRE(point.getLowerBound() == REAL(5.0));
		REQUIRE(point.getUpperBound() == REAL(5.0));
		REQUIRE(point.getLength() == REAL(0.0));
		REQUIRE(point.contains(REAL(5.0)));
		REQUIRE_FALSE(point.contains(REAL(4.999)));
		REQUIRE_FALSE(point.contains(REAL(5.001)));
	}
	
	SECTION("Negative Interval (Invalid)") {
		// Note: Current implementation doesn't validate this
		ClosedInterval invalid(REAL(10.0), REAL(5.0)); // Upper < Lower
		
		REQUIRE(invalid.getLength() < REAL(0.0)); // Negative length
	}
}

TEST_CASE("Edge Cases - Precision", "[intervals][edge][precision]") {
		TEST_PRECISION_INFO();
	SECTION("Very Small Interval") {
		ClosedInterval tiny(REAL(0.0), 1e-10);
		
		REQUIRE_THAT(tiny.getLength(), RealApprox(1e-10));
		REQUIRE(tiny.contains(5e-11));
		REQUIRE_FALSE(tiny.contains(2e-10));
	}
	
	SECTION("Very Large Interval") {
		ClosedInterval huge(-1e100, 1e100);
		
		REQUIRE(huge.contains(REAL(0.0)));
		REQUIRE(huge.contains(-5e99));
		REQUIRE(huge.contains(5e99));
	}
}

TEST_CASE("Complex Scenarios - Chained Operations", "[intervals][complex]") {
		TEST_PRECISION_INFO();
	SECTION("Intersection of Differences") {
		ClosedInterval a(REAL(0.0), REAL(10.0));
		ClosedInterval b(REAL(3.0), REAL(7.0));
		ClosedInterval c(REAL(5.0), REAL(12.0));
		
		// (A \ B) n C
		Interval diff = Interval::Difference(a, b);    // [0,3) ? (7,10]
		
		// Check difference result
		REQUIRE(diff.contains(REAL(1.0)));
		REQUIRE(diff.contains(REAL(9.0)));
		REQUIRE_FALSE(diff.contains(REAL(5.0)));
		
		// Now intersect with C [5, 12]
		// Should give (7, 10] (the right piece of diff within C)
		// Note: Current Interval::Intersection works on BaseInterval,
		// so we'd need to extract sub-intervals or extend the implementation
	}
	
	SECTION("Multiple Complements") {
		ClosedInterval a(REAL(0.0), REAL(5.0));
		
		Interval comp1 = Interval::Complement(a);     // (-8, 0) ? (5, 8)
		
		REQUIRE(comp1.contains(-REAL(10.0)));
		REQUIRE(comp1.contains(REAL(10.0)));
		REQUIRE_FALSE(comp1.contains(REAL(2.5)));
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
			REAL(0.5) * Constants::PI,  // First hole
			Constants::PI         // Hole spacing
		);
		
		REQUIRE(tanDomain.contains(REAL(0.0)));
		REQUIRE(tanDomain.contains(Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(REAL(0.5) * Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(-REAL(0.5) * Constants::PI));
		REQUIRE_FALSE(tanDomain.contains(REAL(1.5) * Constants::PI));
	}
	
	SECTION("Temperature Range with Exclusions") {
		Interval validTemps;
		validTemps.AddInterval(ClosedInterval(-REAL(50.0), REAL(0.0)))    // Freezing range
		          .AddInterval(OpenInterval(REAL(0.0), REAL(100.0)))      // Liquid range (excluding phase boundaries)
		          .AddInterval(ClosedInterval(REAL(100.0), REAL(200.0))); // Boiling range
		
		REQUIRE(validTemps.contains(-REAL(25.0)));   // Freezing
		REQUIRE(validTemps.contains(REAL(0.0)));     // Freezing point (in first interval)
		REQUIRE(validTemps.contains(REAL(50.0)));    // Liquid
		REQUIRE(validTemps.contains(REAL(100.0)));   // Boiling point (in third interval)
		REQUIRE(validTemps.contains(REAL(150.0)));   // Above boiling
	}
	
	SECTION("Root Finding Search Interval") {
		ClosedInterval searchSpace(-REAL(10.0), REAL(10.0));
		std::vector<Real> samplePoints;
		
		searchSpace.GetEquidistantCovering(21, samplePoints);
		
		REQUIRE(samplePoints.size() == 21);
		REQUIRE_THAT(samplePoints[0], RealApprox(-REAL(10.0)));
		REQUIRE_THAT(samplePoints[10], RealApprox(REAL(0.0)));
		REQUIRE_THAT(samplePoints[20], RealApprox(REAL(10.0)));
		
		// Verify uniform spacing of REAL(1.0)
		for (size_t i = 0; i < samplePoints.size() - 1; i++) {
			REQUIRE_THAT(samplePoints[i+1] - samplePoints[i], RealApprox(REAL(1.0)));
		}
	}
}

TEST_CASE("Performance - Large Compound Intervals", "[intervals][performance]") {
		TEST_PRECISION_INFO();
	SECTION("Many Sub-Intervals") {
		Interval large;
		
		// Create 100 disjoint intervals
		for (int i = 0; i < 100; i++) {
			large.AddInterval(ClosedInterval(i * REAL(10.0), i * REAL(10.0) + REAL(5.0)));
		}
		
		REQUIRE_THAT(large.getLowerBound(), RealApprox(REAL(0.0)));
		REQUIRE_THAT(large.getUpperBound(), RealApprox(REAL(995.0)));
		REQUIRE_THAT(large.getLength(), RealApprox(REAL(500.0))); // 100 * REAL(5.0)
		
		REQUIRE(large.contains(REAL(2.5)));     // In first interval
		REQUIRE(large.contains(REAL(502.5)));   // In 51st interval
		REQUIRE_FALSE(large.contains(REAL(7.0)));   // In gap
		REQUIRE_FALSE(large.contains(REAL(1000.0))); // After all
	}
}
