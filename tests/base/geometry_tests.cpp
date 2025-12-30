#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry.h"
#include "base/VectorN.h"
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MML::Tests::Base::GeometryTests
{
	/*********************************************************************/
	/*****                 Point2Cartesian tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point2Cartesian::Init", "[Point2Cartesian]") {
			TEST_PRECISION_INFO();
		Point2Cartesian a(1, 1);

		REQUIRE(a.X() == REAL(1.0));
		REQUIRE(a.Y() == REAL(1.0));
	}
	TEST_CASE("Point2Cartesian::Dist", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(0, 0);
		Point2Cartesian b(3, 4);
		REQUIRE_THAT(a.Dist(b) , RealWithinRel(REAL(5.0), REAL(1e-5)));
		REQUIRE_THAT(b.Dist(a) , RealWithinRel(REAL(5.0), REAL(1e-5)));
		REQUIRE_THAT(a.Dist(a) , RealWithinRel(REAL(0.0), REAL(1e-5)));
	}

	TEST_CASE("Point2Cartesian::op+", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(1, 2);
		Point2Cartesian b(3, 4);
		Point2Cartesian c = a + b;
		REQUIRE_THAT(c.X() , RealWithinRel(REAL(4.0), REAL(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(REAL(6.0), REAL(1e-5)));
	}

	TEST_CASE("Point2Cartesian::op*", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(2, 3);
		Real scalar = REAL(4.0);
		Point2Cartesian b = a * scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(REAL(8.0), REAL(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(REAL(12.0), REAL(1e-5)));

		Point2Cartesian c = scalar * a;
		REQUIRE_THAT(c.X() , RealWithinRel(REAL(8.0), REAL(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(REAL(12.0), REAL(1e-5)));
	}

	TEST_CASE("Point2Cartesian::op/", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(8, 12);
		Real scalar = REAL(4.0);
		Point2Cartesian b = a / scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(REAL(2.0), REAL(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(REAL(3.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****                 Point2Polar tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point2Polar::Init", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Polar p(REAL(5.0), MML::Constants::PI / 4);
		REQUIRE_THAT(p.R() , RealWithinRel(REAL(5.0), REAL(1e-5)));
		REQUIRE_THAT(p.Phi() , RealWithinRel(MML::Constants::PI / 4, REAL(1e-5)));
	}

	TEST_CASE("Point2Polar::Dist", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Polar a(REAL(2.0), REAL(0.0));
		Point2Polar b(REAL(2.0), MML::Constants::PI);
		// These points are at (2,0) and (-2,0) in Cartesian, so distance is 4
		REQUIRE_THAT(a.Dist(b) , RealWithinRel(REAL(4.0), REAL(1e-5)));
		REQUIRE_THAT(a.Dist(a) , RealWithinRel(REAL(0.0), REAL(1e-5)));
	}

	TEST_CASE("Point2Polar::GetFromCartesian", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian cart(REAL(0.0), REAL(2.0));
		Point2Polar pol(cart);
		REQUIRE_THAT(pol.R() , RealWithinRel(REAL(2.0), REAL(1e-5)));
		REQUIRE_THAT(pol.Phi() , RealWithinRel(MML::Constants::PI / 2, REAL(1e-5)));
	}

	TEST_CASE("Point2Polar::TransfToCartesian", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Polar pol(REAL(3.0), MML::Constants::PI / 6);
		Point2Cartesian cart = pol.TransfToCart();
		REQUIRE_THAT(cart.X(), RealWithinRel(REAL(3.0) * cos(MML::Constants::PI / 6), 1e-14));
		REQUIRE_THAT(cart.Y(), RealWithinRel(REAL(3.0) * sin(MML::Constants::PI / 6), 1e-14));
	}

	/*********************************************************************/
	/*****                 Point3Cartesian tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point3Cartesian::Init", "[Point3Cartesian]") 
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 1, 1);

		REQUIRE(a.X() == REAL(1.0));
		REQUIRE(a.Y() == REAL(1.0));
		REQUIRE(a.Z() == REAL(1.0));
	}
	TEST_CASE("Point3Cartesian::Dist", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(0, 0, 0);
		Point3Cartesian b(1, 2, 2);
		REQUIRE_THAT(a.Dist(b) , RealWithinRel(REAL(3.0), REAL(1e-5)));
		REQUIRE_THAT(b.Dist(a) , RealWithinRel(REAL(3.0), REAL(1e-5)));
		REQUIRE_THAT(a.Dist(a) , RealWithinRel(REAL(0.0), REAL(1e-5)));
	}

	TEST_CASE("Point3Cartesian::op+", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 2, 3);
		Point3Cartesian b(4, 5, 6);
		Point3Cartesian c = a + b;
		REQUIRE_THAT(c.X() , RealWithinRel(REAL(5.0), REAL(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(REAL(7.0), REAL(1e-5)));
		REQUIRE_THAT(c.Z() , RealWithinRel(REAL(9.0), REAL(1e-5)));
	}

	TEST_CASE("Point3Cartesian::op*", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(2, 3, 4);
		Real scalar = REAL(5.0);
		Point3Cartesian b = a * scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(REAL(10.0), REAL(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(REAL(15.0), REAL(1e-5)));
		REQUIRE_THAT(b.Z() , RealWithinRel(REAL(20.0), REAL(1e-5)));

		Point3Cartesian c = scalar * a;
		REQUIRE_THAT(c.X() , RealWithinRel(REAL(10.0), REAL(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(REAL(15.0), REAL(1e-5)));
		REQUIRE_THAT(c.Z() , RealWithinRel(REAL(20.0), REAL(1e-5)));
	}

	TEST_CASE("Point3Cartesian::op/", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(10, 20, 30);
		Real scalar = REAL(10.0);
		Point3Cartesian b = a / scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(REAL(1.0), REAL(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(REAL(2.0), REAL(1e-5)));
		REQUIRE_THAT(b.Z() , RealWithinRel(REAL(3.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****                     Triangle tests                        *****/
	/*********************************************************************/
	TEST_CASE("Triangle::Init", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(REAL(3.0), REAL(4.0), REAL(5.0));
		// Sides should be set correctly
		REQUIRE_THAT(t.A() , RealWithinRel(REAL(3.0), REAL(1e-5)));
		REQUIRE_THAT(t.B() , RealWithinRel(REAL(4.0), REAL(1e-5)));
		REQUIRE_THAT(t.C() , RealWithinRel(REAL(5.0), REAL(1e-5)));
	}

	TEST_CASE("Triangle::Area", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(REAL(3.0), REAL(4.0), REAL(5.0));
		// Heron's formula: s = (a+b+c)/2 = 6, area = sqrt(s*(s-a)*(s-b)*(s-c)) = 6
		Real s = (REAL(3.0) + REAL(4.0) + REAL(5.0)) / REAL(2.0);
		Real expected_area = sqrt(s * (s - REAL(3.0)) * (s - REAL(4.0)) * (s - REAL(5.0)));
		REQUIRE_THAT(t.Area() , RealWithinRel(expected_area, REAL(1e-5)));
	}

	TEST_CASE("Triangle::IsRight", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(REAL(3.0), REAL(4.0), REAL(5.0));
		REQUIRE(t.IsRight());

		Triangle t2(REAL(2.0), REAL(2.0), REAL(3.0));
		REQUIRE_FALSE(t2.IsRight());
	}

	TEST_CASE("Triangle::IsIsosceles", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(REAL(2.0), REAL(2.0), REAL(3.0));
		REQUIRE(t.IsIsosceles());

		Triangle t2(REAL(3.0), REAL(4.0), REAL(5.0));
		REQUIRE_FALSE(t2.IsIsosceles());
	}

	TEST_CASE("Triangle::IsEquilateral", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(REAL(2.0), REAL(2.0), REAL(2.0));
		REQUIRE(t.IsEquilateral());

		Triangle t2(REAL(3.0), REAL(4.0), REAL(5.0));
		REQUIRE_FALSE(t2.IsEquilateral());
	}
} // namespace MML::Tests::Base::GeometryTests