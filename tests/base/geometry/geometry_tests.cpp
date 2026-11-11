#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Geometry/Geometry.h>
#include <mml/base/Vector/VectorN.h>
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

		REQUIRE(a.X() == static_cast<MML::Real>(1.0));
		REQUIRE(a.Y() == static_cast<MML::Real>(1.0));
	}
	TEST_CASE("Point2Cartesian::Dist", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(0, 0);
		Point2Cartesian b(3, 4);
		REQUIRE_THAT(a.Dist(b) , RealWithinRel(static_cast<MML::Real>(5.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Dist(a) , RealWithinRel(static_cast<MML::Real>(5.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(a.Dist(a) , RealWithinRel(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point2Cartesian::op+", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(1, 2);
		Point2Cartesian b(3, 4);
		Point2Cartesian c = a + b;
		REQUIRE_THAT(c.X() , RealWithinRel(static_cast<MML::Real>(4.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(static_cast<MML::Real>(6.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point2Cartesian::op*", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(2, 3);
		Real scalar = static_cast<MML::Real>(4.0);
		Point2Cartesian b = a * scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(static_cast<MML::Real>(8.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(static_cast<MML::Real>(12.0), static_cast<MML::Real>(1e-5)));

		Point2Cartesian c = scalar * a;
		REQUIRE_THAT(c.X() , RealWithinRel(static_cast<MML::Real>(8.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(static_cast<MML::Real>(12.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point2Cartesian::op/", "[Point2Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian a(8, 12);
		Real scalar = static_cast<MML::Real>(4.0);
		Point2Cartesian b = a / scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-5)));
	}

	/*********************************************************************/
	/*****                 Point2Polar tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point2Polar::Init", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Polar p(static_cast<MML::Real>(5.0), MML::Constants::PI / 4);
		REQUIRE_THAT(p.R() , RealWithinRel(static_cast<MML::Real>(5.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(p.Phi() , RealWithinRel(MML::Constants::PI / 4, static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point2Polar::Dist", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Polar a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0));
		Point2Polar b(static_cast<MML::Real>(2.0), MML::Constants::PI);
		// These points are at (2,0) and (-2,0) in Cartesian, so distance is 4
		REQUIRE_THAT(a.Dist(b) , RealWithinRel(static_cast<MML::Real>(4.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(a.Dist(a) , RealWithinRel(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point2Polar::GetFromCartesian", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Cartesian cart(static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0));
		Point2Polar pol(cart);
		REQUIRE_THAT(pol.R() , RealWithinRel(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(pol.Phi() , RealWithinRel(MML::Constants::PI / 2, static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point2Polar::TransfToCartesian", "[Point2Polar]")
	{
			TEST_PRECISION_INFO();
		Point2Polar pol(static_cast<MML::Real>(3.0), MML::Constants::PI / 6);
		Point2Cartesian cart = pol.TransfToCart();
		REQUIRE_THAT(cart.X(), RealWithinRel(static_cast<MML::Real>(3.0) * cos(MML::Constants::PI / 6), TOL(1e-14, 1e-5)));
		REQUIRE_THAT(cart.Y(), RealWithinRel(static_cast<MML::Real>(3.0) * sin(MML::Constants::PI / 6), TOL(1e-14, 1e-5)));
	}

	/*********************************************************************/
	/*****                 Point3Cartesian tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point3Cartesian::Init", "[Point3Cartesian]") 
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 1, 1);

		REQUIRE(a.X() == static_cast<MML::Real>(1.0));
		REQUIRE(a.Y() == static_cast<MML::Real>(1.0));
		REQUIRE(a.Z() == static_cast<MML::Real>(1.0));
	}
	TEST_CASE("Point3Cartesian::Dist", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(0, 0, 0);
		Point3Cartesian b(1, 2, 2);
		REQUIRE_THAT(a.Dist(b) , RealWithinRel(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Dist(a) , RealWithinRel(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(a.Dist(a) , RealWithinRel(static_cast<MML::Real>(0.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point3Cartesian::op+", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 2, 3);
		Point3Cartesian b(4, 5, 6);
		Point3Cartesian c = a + b;
		REQUIRE_THAT(c.X() , RealWithinRel(static_cast<MML::Real>(5.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(static_cast<MML::Real>(7.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(c.Z() , RealWithinRel(static_cast<MML::Real>(9.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point3Cartesian::op*", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(2, 3, 4);
		Real scalar = static_cast<MML::Real>(5.0);
		Point3Cartesian b = a * scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(static_cast<MML::Real>(10.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(static_cast<MML::Real>(15.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Z() , RealWithinRel(static_cast<MML::Real>(20.0), static_cast<MML::Real>(1e-5)));

		Point3Cartesian c = scalar * a;
		REQUIRE_THAT(c.X() , RealWithinRel(static_cast<MML::Real>(10.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(c.Y() , RealWithinRel(static_cast<MML::Real>(15.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(c.Z() , RealWithinRel(static_cast<MML::Real>(20.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Point3Cartesian::op/", "[Point3Cartesian]")
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(10, 20, 30);
		Real scalar = static_cast<MML::Real>(10.0);
		Point3Cartesian b = a / scalar;
		REQUIRE_THAT(b.X() , RealWithinRel(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Y() , RealWithinRel(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(b.Z() , RealWithinRel(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-5)));
	}

	/*********************************************************************/
	/*****                     Triangle tests                        *****/
	/*********************************************************************/
	TEST_CASE("Triangle::Init", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		// Sides should be set correctly
		REQUIRE_THAT(t.A() , RealWithinRel(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(t.B() , RealWithinRel(static_cast<MML::Real>(4.0), static_cast<MML::Real>(1e-5)));
		REQUIRE_THAT(t.C() , RealWithinRel(static_cast<MML::Real>(5.0), static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Triangle::Area", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		// Heron's formula: s = (a+b+c)/2 = 6, area = sqrt(s*(s-a)*(s-b)*(s-c)) = 6
		Real s = (static_cast<MML::Real>(3.0) + static_cast<MML::Real>(4.0) + static_cast<MML::Real>(5.0)) / static_cast<MML::Real>(2.0);
		Real expected_area = sqrt(s * (s - static_cast<MML::Real>(3.0)) * (s - static_cast<MML::Real>(4.0)) * (s - static_cast<MML::Real>(5.0)));
		REQUIRE_THAT(t.Area() , RealWithinRel(expected_area, static_cast<MML::Real>(1e-5)));
	}

	TEST_CASE("Triangle::IsRight", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE(t.IsRight());

		Triangle t2(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		REQUIRE_FALSE(t2.IsRight());
	}

	TEST_CASE("Triangle::IsIsosceles", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		REQUIRE(t.IsIsosceles());

		Triangle t2(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE_FALSE(t2.IsIsosceles());
	}

	TEST_CASE("Triangle::IsEquilateral", "[Triangle]")
	{
			TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0));
		REQUIRE(t.IsEquilateral());

		Triangle t2(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE_FALSE(t2.IsEquilateral());
	}

	/*********************************************************************/
	/*****           Point2Cartesian additional tests                *****/
	/*********************************************************************/
	TEST_CASE("Point2Cartesian::DefaultConstructor", "[Point2Cartesian]")
	{
		Point2Cartesian p;
		REQUIRE(p.X() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Y() == static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Point2Cartesian::MutableAccess", "[Point2Cartesian]")
	{
		Point2Cartesian p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		p.X() = static_cast<MML::Real>(5.0);
		p.Y() = static_cast<MML::Real>(10.0);
		REQUIRE(p.X() == static_cast<MML::Real>(5.0));
		REQUIRE(p.Y() == static_cast<MML::Real>(10.0));
	}

	TEST_CASE("Point2Cartesian::EqualityOperators", "[Point2Cartesian]")
	{
		Point2Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		Point2Cartesian b(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		Point2Cartesian c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0));
		
		REQUIRE(a == b);
		REQUIRE_FALSE(a == c);
		REQUIRE(a != c);
		REQUIRE_FALSE(a != b);
	}

	TEST_CASE("Point2Cartesian::IsEqualTo", "[Point2Cartesian]")
	{
		Point2Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		Point2Cartesian b(static_cast<MML::Real>(1.0) + TOL3(1e-12, 1e-7, 1e-16), static_cast<MML::Real>(2.0) + TOL3(1e-12, 1e-7, 1e-16));
		Point2Cartesian c(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		
		REQUIRE(a.IsEqualTo(b));
		REQUIRE_FALSE(a.IsEqualTo(c));
		REQUIRE(a.IsEqualTo(a));  // same point
	}

	TEST_CASE("Point2Cartesian::op-", "[Point2Cartesian]")
	{
		Point2Cartesian a(static_cast<MML::Real>(5.0), static_cast<MML::Real>(7.0));
		Point2Cartesian b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point2Cartesian c = a - b;
		REQUIRE_THAT(c.X(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Y(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point2Cartesian::CompoundAssignment", "[Point2Cartesian]")
	{
		Point2Cartesian p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		
		p += Point2Cartesian(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
		
		p -= Point2Cartesian(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0));
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		
		p *= static_cast<MML::Real>(2.0);
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(10.0), TOL(1e-10, 1e-5)));
		
		p /= static_cast<MML::Real>(2.0);
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****             Point2Polar additional tests                  *****/
	/*********************************************************************/
	TEST_CASE("Point2Polar::DefaultConstructor", "[Point2Polar]")
	{
		Point2Polar p;
		REQUIRE(p.R() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Phi() == static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Point2Polar::MutableAccess", "[Point2Polar]")
	{
		Point2Polar p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5));
		p.R() = static_cast<MML::Real>(5.0);
		p.Phi() = static_cast<MML::Real>(1.0);
		REQUIRE(p.R() == static_cast<MML::Real>(5.0));
		REQUIRE(p.Phi() == static_cast<MML::Real>(1.0));
	}

	TEST_CASE("Point2Polar::EqualityOperators", "[Point2Polar]")
	{
		Point2Polar a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.5));
		Point2Polar b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.5));
		Point2Polar c(static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0));
		
		REQUIRE(a == b);
		REQUIRE_FALSE(a == c);
		REQUIRE(a != c);
		REQUIRE_FALSE(a != b);
	}

	TEST_CASE("Point2Polar::IsEqualTo", "[Point2Polar]")
	{
		Point2Polar a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.5));
		Point2Polar b(static_cast<MML::Real>(2.0) + TOL(1e-12, 1e-5), static_cast<MML::Real>(0.5));
		
		REQUIRE(a.IsEqualTo(b));
		REQUIRE(a.IsEqualTo(a));
	}

	TEST_CASE("Point2Polar::RoundTrip", "[Point2Polar]")
	{
		TEST_PRECISION_INFO();
		// Convert Cartesian -> Polar -> Cartesian
		Point2Cartesian original(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
		Point2Polar polar(original);
		Point2Cartesian back = polar.TransfToCart();
		
		REQUIRE_THAT(back.X(), RealWithinAbs(original.X(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(back.Y(), RealWithinAbs(original.Y(), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****           Point3Cartesian additional tests                *****/
	/*********************************************************************/
	TEST_CASE("Point3Cartesian::DefaultConstructor", "[Point3Cartesian]")
	{
		Point3Cartesian p;
		REQUIRE(p.X() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Y() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Z() == static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Point3Cartesian::MutableAccess", "[Point3Cartesian]")
	{
		Point3Cartesian p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		p.X() = static_cast<MML::Real>(5.0);
		p.Y() = static_cast<MML::Real>(10.0);
		p.Z() = static_cast<MML::Real>(15.0);
		REQUIRE(p.X() == static_cast<MML::Real>(5.0));
		REQUIRE(p.Y() == static_cast<MML::Real>(10.0));
		REQUIRE(p.Z() == static_cast<MML::Real>(15.0));
	}

	TEST_CASE("Point3Cartesian::EqualityOperators", "[Point3Cartesian]")
	{
		Point3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point3Cartesian b(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point3Cartesian c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0));
		
		REQUIRE(a == b);
		REQUIRE_FALSE(a == c);
		REQUIRE(a != c);
		REQUIRE_FALSE(a != b);
	}

	TEST_CASE("Point3Cartesian::IsEqualTo", "[Point3Cartesian]")
	{
		Point3Cartesian a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point3Cartesian b(static_cast<MML::Real>(1.0) + TOL3(1e-12, 1e-7, 1e-16), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point3Cartesian c(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
		
		REQUIRE(a.IsEqualTo(b));
		REQUIRE_FALSE(a.IsEqualTo(c));
	}

	TEST_CASE("Point3Cartesian::op-", "[Point3Cartesian]")
	{
		Point3Cartesian a(static_cast<MML::Real>(5.0), static_cast<MML::Real>(7.0), static_cast<MML::Real>(9.0));
		Point3Cartesian b(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
		Point3Cartesian c = a - b;
		REQUIRE_THAT(c.X(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Y(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Z(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Cartesian::CompoundAssignment", "[Point3Cartesian]")
	{
		Point3Cartesian p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		
		p += Point3Cartesian(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0));
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Z(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		
		p -= Point3Cartesian(static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.5));
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(1.5), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(2.5), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Z(), RealWithinAbs(static_cast<MML::Real>(3.5), TOL(1e-10, 1e-5)));
		
		p *= static_cast<MML::Real>(2.0);
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Z(), RealWithinAbs(static_cast<MML::Real>(7.0), TOL(1e-10, 1e-5)));
		
		p /= static_cast<MML::Real>(2.0);
		REQUIRE_THAT(p.X(), RealWithinAbs(static_cast<MML::Real>(1.5), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Y(), RealWithinAbs(static_cast<MML::Real>(2.5), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Z(), RealWithinAbs(static_cast<MML::Real>(3.5), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****             Point3Spherical tests                         *****/
	/*********************************************************************/
	TEST_CASE("Point3Spherical::DefaultConstructor", "[Point3Spherical]")
	{
		Point3Spherical p;
		REQUIRE(p.R() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Theta() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Phi() == static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Point3Spherical::ParameterConstructor", "[Point3Spherical]")
	{
		Point3Spherical p(static_cast<MML::Real>(5.0), Constants::PI / 4, Constants::PI / 2);
		REQUIRE_THAT(p.R(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Theta(), RealWithinAbs(Constants::PI / 4, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Phi(), RealWithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Spherical::FromCartesian", "[Point3Spherical]")
	{
		TEST_PRECISION_INFO();
		// Point on z-axis (north pole)
		Point3Cartesian cartZ(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(5.0));
		Point3Spherical sphZ(cartZ);
		REQUIRE_THAT(sphZ.R(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sphZ.Theta(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));  // theta=0 at north pole
		
		// Point on x-axis
		Point3Cartesian cartX(static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));
		Point3Spherical sphX(cartX);
		REQUIRE_THAT(sphX.R(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sphX.Theta(), RealWithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));  // on equator
		REQUIRE_THAT(sphX.Phi(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Spherical::TransfToCart", "[Point3Spherical]")
	{
		TEST_PRECISION_INFO();
		// Unit sphere point at equator, phi=0
		Point3Spherical sph(static_cast<MML::Real>(1.0), Constants::PI / 2, static_cast<MML::Real>(0.0));
		Point3Cartesian cart = sph.TransfToCart();
		REQUIRE_THAT(cart.X(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(cart.Y(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(cart.Z(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Spherical::RoundTrip", "[Point3Spherical]")
	{
		TEST_PRECISION_INFO();
		Point3Cartesian original(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point3Spherical spherical(original);
		Point3Cartesian back = spherical.TransfToCart();
		
		REQUIRE_THAT(back.X(), RealWithinAbs(original.X(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(back.Y(), RealWithinAbs(original.Y(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(back.Z(), RealWithinAbs(original.Z(), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Spherical::EqualityOperators", "[Point3Spherical]")
	{
		Point3Spherical a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(1.0));
		Point3Spherical b(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(1.0));
		Point3Spherical c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(2.0));
		
		REQUIRE(a == b);
		REQUIRE_FALSE(a == c);
		REQUIRE(a != c);
	}

	TEST_CASE("Point3Spherical::Dist", "[Point3Spherical]")
	{
		TEST_PRECISION_INFO();
		Point3Spherical a(static_cast<MML::Real>(1.0), Constants::PI / 2, static_cast<MML::Real>(0.0));  // on x-axis
		Point3Spherical b(static_cast<MML::Real>(1.0), Constants::PI / 2, Constants::PI);  // on negative x-axis
		
		// These are at (1,0,0) and (-1,0,0), distance = 2
		REQUIRE_THAT(a.Dist(a), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(a.Dist(b), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));

		// Points at same theta, different phi (on equatorial plane)
		// a=(1, pi/2, 0) -> (1,0,0), c=(1, pi/2, pi/2) -> (0,1,0), distance = sqrt(2)
		Point3Spherical c(static_cast<MML::Real>(1.0), Constants::PI / 2, Constants::PI / 2);
		REQUIRE_THAT(a.Dist(c), RealWithinAbs(std::sqrt(static_cast<MML::Real>(2.0)), TOL(1e-10, 1e-5)));

		// Same phi, different theta: a=(1, pi/4, 0), d=(1, pi/4, pi)
		// a -> (sin(pi/4), 0, cos(pi/4)), d -> (-sin(pi/4), 0, cos(pi/4))
		// distance = 2*sin(pi/4) = sqrt(2)
		Point3Spherical e(static_cast<MML::Real>(1.0), Constants::PI / 4, static_cast<MML::Real>(0.0));
		Point3Spherical f(static_cast<MML::Real>(1.0), Constants::PI / 4, Constants::PI);
		REQUIRE_THAT(e.Dist(f), RealWithinAbs(std::sqrt(static_cast<MML::Real>(2.0)), TOL(1e-10, 1e-5)));

		// Different radii: (2, pi/2, 0) -> (2,0,0) and (3, pi/2, 0) -> (3,0,0), distance = 1
		Point3Spherical g(static_cast<MML::Real>(2.0), Constants::PI / 2, static_cast<MML::Real>(0.0));
		Point3Spherical h(static_cast<MML::Real>(3.0), Constants::PI / 2, static_cast<MML::Real>(0.0));
		REQUIRE_THAT(g.Dist(h), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));

		// Cross-validate with Cartesian distance
		Point3Spherical s1(static_cast<MML::Real>(2.0), Constants::PI / 3, Constants::PI / 4);
		Point3Spherical s2(static_cast<MML::Real>(3.0), Constants::PI / 6, Constants::PI / 2);
		Point3Cartesian c1 = s1.TransfToCart();
		Point3Cartesian c2 = s2.TransfToCart();
		Real cartDist = c1.Dist(c2);
		REQUIRE_THAT(s1.Dist(s2), RealWithinAbs(cartDist, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****             Point3Cylindrical tests                       *****/
	/*********************************************************************/
	TEST_CASE("Point3Cylindrical::DefaultConstructor", "[Point3Cylindrical]")
	{
		Point3Cylindrical p;
		REQUIRE(p.R() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Phi() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Z() == static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Point3Cylindrical::ParameterConstructor", "[Point3Cylindrical]")
	{
		Point3Cylindrical p(static_cast<MML::Real>(3.0), Constants::PI / 4, static_cast<MML::Real>(5.0));
		REQUIRE_THAT(p.R(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Phi(), RealWithinAbs(Constants::PI / 4, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.Z(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Cylindrical::FromCartesian", "[Point3Cylindrical]")
	{
		TEST_PRECISION_INFO();
		Point3Cartesian cart(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		Point3Cylindrical cyl(cart);
		
		REQUIRE_THAT(cyl.R(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));  // sqrt(9+16)
		REQUIRE_THAT(cyl.Z(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Cylindrical::TransfToCart", "[Point3Cylindrical]")
	{
		TEST_PRECISION_INFO();
		Point3Cylindrical cyl(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(3.0));  // on x-axis at z=3
		Point3Cartesian cart = cyl.TransfToCart();
		
		REQUIRE_THAT(cart.X(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(cart.Y(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(cart.Z(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Cylindrical::RoundTrip", "[Point3Cylindrical]")
	{
		TEST_PRECISION_INFO();
		Point3Cartesian original(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		Point3Cylindrical cylindrical(original);
		Point3Cartesian back = cylindrical.TransfToCart();
		
		REQUIRE_THAT(back.X(), RealWithinAbs(original.X(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(back.Y(), RealWithinAbs(original.Y(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(back.Z(), RealWithinAbs(original.Z(), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Point3Cylindrical::EqualityOperators", "[Point3Cylindrical]")
	{
		Point3Cylindrical a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(2.0));
		Point3Cylindrical b(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(2.0));
		Point3Cylindrical c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(3.0));
		
		REQUIRE(a == b);
		REQUIRE_FALSE(a == c);
		REQUIRE(a != c);
	}

	TEST_CASE("Point3Cylindrical::Dist", "[Point3Cylindrical]")
	{
		TEST_PRECISION_INFO();
		Point3Cylindrical a(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0));  // origin
		Point3Cylindrical b(static_cast<MML::Real>(3.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0));  // (3,0,4) in Cartesian
		
		REQUIRE_THAT(a.Dist(b), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));  // 3-4-5 triangle
	}

	/*********************************************************************/
	/*****             Triangle additional tests                     *****/
	/*********************************************************************/
	TEST_CASE("Triangle::DefaultConstructor", "[Triangle]")
	{
		Triangle t;
		REQUIRE(t.A() == static_cast<MML::Real>(0.0));
		REQUIRE(t.B() == static_cast<MML::Real>(0.0));
		REQUIRE(t.C() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(t.IsValid());
	}

	TEST_CASE("Triangle::IsValid", "[Triangle]")
	{
		Triangle valid(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE(valid.IsValid());
		
		// Degenerate: sum of two sides equals third
		Triangle degenerate(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		REQUIRE_FALSE(degenerate.IsValid());
		
		// Invalid: violates triangle inequality
		Triangle invalid(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(10.0));
		REQUIRE_FALSE(invalid.IsValid());
		
		// Zero side
		Triangle zeroSide(static_cast<MML::Real>(0.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE_FALSE(zeroSide.IsValid());
	}

	TEST_CASE("Triangle::Perimeter", "[Triangle]")
	{
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE_THAT(t.Perimeter(), RealWithinAbs(static_cast<MML::Real>(12.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(t.Semiperimeter(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Triangle::Angles", "[Triangle]")
	{
		TEST_PRECISION_INFO();
		// Right triangle 3-4-5: angle opposite to 5 is 90°
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		
		// AngleC is opposite to side C (the hypotenuse)
		REQUIRE_THAT(t.AngleC(), RealWithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));
		
		// Sum of all angles = π
		Real sum = t.AngleA() + t.AngleB() + t.AngleC();
		REQUIRE_THAT(sum, RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Triangle::Altitudes", "[Triangle]")
	{
		TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		Real area = t.Area();  // = 6 for 3-4-5 triangle
		
		// Altitude = 2*Area/base
		REQUIRE_THAT(t.AltitudeToA(), RealWithinAbs(static_cast<MML::Real>(2.0) * area / static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(t.AltitudeToB(), RealWithinAbs(static_cast<MML::Real>(2.0) * area / static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(t.AltitudeToC(), RealWithinAbs(static_cast<MML::Real>(2.0) * area / static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Triangle::Medians", "[Triangle]")
	{
		TEST_PRECISION_INFO();
		// Equilateral triangle: all medians equal
		Triangle eq(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0));
		REQUIRE_THAT(eq.MedianToA(), RealWithinAbs(eq.MedianToB(), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(eq.MedianToB(), RealWithinAbs(eq.MedianToC(), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Triangle::Inradius_Circumradius", "[Triangle]")
	{
		TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		
		// For 3-4-5 right triangle: inradius = (a+b-c)/2 = (3+4-5)/2 = 1
		REQUIRE_THAT(t.Inradius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		
		// Circumradius = c/2 = 5/2 = 2.5 for right triangle
		REQUIRE_THAT(t.Circumradius(), RealWithinAbs(static_cast<MML::Real>(2.5), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Triangle::Classification", "[Triangle]")
	{
		TEST_PRECISION_INFO();
		
		// Acute triangle: all angles < 90°
		Triangle acute(static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0), static_cast<MML::Real>(7.0));
		REQUIRE(acute.IsAcute());
		REQUIRE_FALSE(acute.IsObtuse());
		REQUIRE_FALSE(acute.IsRight());
		
		// Obtuse triangle: one angle > 90°
		Triangle obtuse(static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.5));
		REQUIRE(obtuse.IsObtuse());
		REQUIRE_FALSE(obtuse.IsAcute());
		
		// Scalene: all sides different
		Triangle scalene(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		REQUIRE(scalene.IsScalene());
	}

	/*********************************************************************/
	/*****                   Circle tests                            *****/
	/*********************************************************************/
	TEST_CASE("Circle::DefaultConstructor", "[Circle]")
	{
		Circle c;
		REQUIRE(c.Radius() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(c.IsValid());
	}

	TEST_CASE("Circle::ParameterConstructor", "[Circle]")
	{
		Circle c(static_cast<MML::Real>(5.0));
		REQUIRE_THAT(c.Radius(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Diameter(), RealWithinAbs(static_cast<MML::Real>(10.0), TOL(1e-10, 1e-5)));
		REQUIRE(c.IsValid());
	}

	TEST_CASE("Circle::Area", "[Circle]")
	{
		TEST_PRECISION_INFO();
		Circle c(static_cast<MML::Real>(2.0));
		REQUIRE_THAT(c.Area(), RealWithinAbs(static_cast<MML::Real>(4.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Circle::Circumference", "[Circle]")
	{
		TEST_PRECISION_INFO();
		Circle c(static_cast<MML::Real>(3.0));
		REQUIRE_THAT(c.Circumference(), RealWithinAbs(static_cast<MML::Real>(6.0) * Constants::PI, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Perimeter(), RealWithinAbs(static_cast<MML::Real>(6.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Circle::ArcAndSector", "[Circle]")
	{
		TEST_PRECISION_INFO();
		Circle c(static_cast<MML::Real>(2.0));
		
		// Quarter circle
		Real angle = Constants::PI / 2;
		REQUIRE_THAT(c.ArcLength(angle), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));  // r*θ = 2*π/2 = π
		REQUIRE_THAT(c.SectorArea(angle), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));  // 0.5*r²*θ = 0.5*4*π/2 = π
	}

	TEST_CASE("Circle::ChordAndSegment", "[Circle]")
	{
		TEST_PRECISION_INFO();
		Circle c(static_cast<MML::Real>(1.0));
		
		// Semicircle: chord = diameter = 2r
		REQUIRE_THAT(c.ChordLength(Constants::PI), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		
		// Segment area for semicircle
		REQUIRE_THAT(c.SegmentArea(Constants::PI), RealWithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Circle::FactoryMethods", "[Circle]")
	{
		TEST_PRECISION_INFO();
		
		Circle fromArea = Circle::FromArea(Constants::PI);
		REQUIRE_THAT(fromArea.Radius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		
		Circle fromCirc = Circle::FromCircumference(static_cast<MML::Real>(2.0) * Constants::PI);
		REQUIRE_THAT(fromCirc.Radius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Circle::InscribedCircumscribed", "[Circle]")
	{
		TEST_PRECISION_INFO();
		Triangle t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		
		Circle inscribed = Circle::Inscribed(t);
		REQUIRE_THAT(inscribed.Radius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		
		Circle circumscribed = Circle::Circumscribed(t);
		REQUIRE_THAT(circumscribed.Radius(), RealWithinAbs(static_cast<MML::Real>(2.5), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                   Ellipse tests                           *****/
	/*********************************************************************/
	TEST_CASE("Ellipse::DefaultConstructor", "[Ellipse]")
	{
		Ellipse e;
		REQUIRE(e.SemiMajor() == static_cast<MML::Real>(0.0));
		REQUIRE(e.SemiMinor() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(e.IsValid());
	}

	TEST_CASE("Ellipse::ParameterConstructor", "[Ellipse]")
	{
		Ellipse e(static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0));
		REQUIRE_THAT(e.SemiMajor(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(e.SemiMinor(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(e.MajorAxis(), RealWithinAbs(static_cast<MML::Real>(10.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(e.MinorAxis(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
		REQUIRE(e.IsValid());
	}

	TEST_CASE("Ellipse::IsCircle", "[Ellipse]")
	{
		Ellipse circle(static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0));
		REQUIRE(circle.IsCircle());
		
		Ellipse notCircle(static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0));
		REQUIRE_FALSE(notCircle.IsCircle());
	}

	TEST_CASE("Ellipse::Area", "[Ellipse]")
	{
		TEST_PRECISION_INFO();
		Ellipse e(static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0));
		REQUIRE_THAT(e.Area(), RealWithinAbs(static_cast<MML::Real>(12.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Ellipse::Eccentricity", "[Ellipse]")
	{
		TEST_PRECISION_INFO();
		// Circle has eccentricity 0
		Ellipse circle(static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0));
		REQUIRE_THAT(circle.Eccentricity(), RealWithinAbs(static_cast<MML::Real>(0.0), TOL(1e-10, 1e-5)));
		
		// e = sqrt(1 - b²/a²) for a=5, b=3: e = sqrt(1 - 9/25) = sqrt(16/25) = 4/5
		Ellipse e(static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0));
		REQUIRE_THAT(e.Eccentricity(), RealWithinAbs(static_cast<MML::Real>(0.8), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Ellipse::FocalProperties", "[Ellipse]")
	{
		TEST_PRECISION_INFO();
		Ellipse e(static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0));
		
		// c = sqrt(a² - b²) = sqrt(25-9) = 4
		REQUIRE_THAT(e.LinearEccentricity(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(e.FocalDistance(), RealWithinAbs(static_cast<MML::Real>(8.0), TOL(1e-10, 1e-5)));
		
		// Semi-latus rectum = b²/a = 9/5 = 1.8
		REQUIRE_THAT(e.SemiLatusRectum(), RealWithinAbs(static_cast<MML::Real>(1.8), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Ellipse::FromFoci", "[Ellipse]")
	{
		TEST_PRECISION_INFO();
		Ellipse e = Ellipse::FromFoci(static_cast<MML::Real>(8.0), static_cast<MML::Real>(10.0));  // focalDist=8, majorAxis=10
		
		REQUIRE_THAT(e.SemiMajor(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(e.SemiMinor(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****              CircularSector tests                         *****/
	/*********************************************************************/
	TEST_CASE("CircularSector::DefaultConstructor", "[CircularSector]")
	{
		CircularSector s;
		REQUIRE(s.Radius() == static_cast<MML::Real>(0.0));
		REQUIRE(s.Angle() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(s.IsValid());
	}

	TEST_CASE("CircularSector::ParameterConstructor", "[CircularSector]")
	{
		CircularSector s(static_cast<MML::Real>(5.0), Constants::PI / 2);
		REQUIRE_THAT(s.Radius(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(s.Angle(), RealWithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(s.AngleDegrees(), RealWithinAbs(static_cast<MML::Real>(90.0), TOL(1e-10, 1e-5)));
		REQUIRE(s.IsValid());
	}

	TEST_CASE("CircularSector::Calculations", "[CircularSector]")
	{
		TEST_PRECISION_INFO();
		CircularSector s(static_cast<MML::Real>(2.0), Constants::PI / 2);  // Quarter of circle with r=2
		
		// Area = 0.5 * r² * θ = 0.5 * 4 * π/2 = π
		REQUIRE_THAT(s.Area(), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
		
		// Arc = r * θ = 2 * π/2 = π
		REQUIRE_THAT(s.ArcLength(), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
		
		// Perimeter = 2r + arc = 4 + π
		REQUIRE_THAT(s.Perimeter(), RealWithinAbs(static_cast<MML::Real>(4.0) + Constants::PI, TOL(1e-10, 1e-5)));
		
		// Chord for 90°: 2r*sin(45°) = 2*2*sqrt(2)/2 = 2*sqrt(2)
		REQUIRE_THAT(s.ChordLength(), RealWithinAbs(static_cast<MML::Real>(2.0) * std::sqrt(static_cast<MML::Real>(2.0)), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****              CircularSegment tests                        *****/
	/*********************************************************************/
	TEST_CASE("CircularSegment::DefaultConstructor", "[CircularSegment]")
	{
		CircularSegment s;
		REQUIRE(s.Radius() == static_cast<MML::Real>(0.0));
		REQUIRE(s.Angle() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(s.IsValid());
	}

	TEST_CASE("CircularSegment::Calculations", "[CircularSegment]")
	{
		TEST_PRECISION_INFO();
		CircularSegment s(static_cast<MML::Real>(2.0), Constants::PI);  // Semicircle
		
		// Area = 0.5 * r² * (θ - sin(θ)) = 0.5 * 4 * (π - 0) = 2π
		REQUIRE_THAT(s.Area(), RealWithinAbs(static_cast<MML::Real>(2.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Arc = r * θ = 2π
		REQUIRE_THAT(s.ArcLength(), RealWithinAbs(static_cast<MML::Real>(2.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Chord = diameter = 4
		REQUIRE_THAT(s.ChordLength(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		
		// Sagitta = r(1 - cos(θ/2)) = 2(1 - cos(π/2)) = 2(1 - 0) = 2
		REQUIRE_THAT(s.Sagitta(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("CircularSegment::FromSagitta", "[CircularSegment]")
	{
		TEST_PRECISION_INFO();
		CircularSegment s = CircularSegment::FromSagitta(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0));  // radius=2, sagitta=2 -> semicircle
		
		REQUIRE_THAT(s.Angle(), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                   Annulus tests                           *****/
	/*********************************************************************/
	TEST_CASE("Annulus::DefaultConstructor", "[Annulus]")
	{
		Annulus a;
		REQUIRE(a.InnerRadius() == static_cast<MML::Real>(0.0));
		REQUIRE(a.OuterRadius() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(a.IsValid());
	}

	TEST_CASE("Annulus::ParameterConstructor", "[Annulus]")
	{
		Annulus a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
		REQUIRE_THAT(a.InnerRadius(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(a.OuterRadius(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(a.Width(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE(a.IsValid());
	}

	TEST_CASE("Annulus::Area", "[Annulus]")
	{
		TEST_PRECISION_INFO();
		Annulus a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0));
		// Area = π(R² - r²) = π(16 - 4) = 12π
		REQUIRE_THAT(a.Area(), RealWithinAbs(static_cast<MML::Real>(12.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Annulus::Circumferences", "[Annulus]")
	{
		TEST_PRECISION_INFO();
		Annulus a(static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0));
		
		REQUIRE_THAT(a.InnerCircumference(), RealWithinAbs(static_cast<MML::Real>(2.0) * Constants::PI, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(a.OuterCircumference(), RealWithinAbs(static_cast<MML::Real>(6.0) * Constants::PI, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(a.MeanCircumference(), RealWithinAbs(static_cast<MML::Real>(4.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Annulus::Circles", "[Annulus]")
	{
		Annulus a(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
		
		Circle inner = a.InnerCircle();
		Circle outer = a.OuterCircle();
		
		REQUIRE_THAT(inner.Radius(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(outer.Radius(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****              RegularPolygon tests                         *****/
	/*********************************************************************/
	TEST_CASE("RegularPolygon::DefaultConstructor", "[RegularPolygon]")
	{
		RegularPolygon p;
		REQUIRE(p.NumSides() == 3);
		REQUIRE(p.SideLength() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(p.IsValid());
	}

	TEST_CASE("RegularPolygon::Square", "[RegularPolygon]")
	{
		TEST_PRECISION_INFO();
		RegularPolygon sq = RegularPolygon::Square(static_cast<MML::Real>(2.0));
		
		REQUIRE(sq.NumSides() == 4);
		REQUIRE_THAT(sq.SideLength(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sq.Perimeter(), RealWithinAbs(static_cast<MML::Real>(8.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sq.Area(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sq.InteriorAngle(), RealWithinAbs(Constants::PI / 2, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("RegularPolygon::Hexagon", "[RegularPolygon]")
	{
		TEST_PRECISION_INFO();
		RegularPolygon hex = RegularPolygon::Hexagon(static_cast<MML::Real>(1.0));
		
		REQUIRE(hex.NumSides() == 6);
		REQUIRE_THAT(hex.Perimeter(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
		
		// Interior angle = (n-2)*π/n = 4π/6 = 2π/3
		REQUIRE_THAT(hex.InteriorAngle(), RealWithinAbs(static_cast<MML::Real>(2.0) * Constants::PI / 3, TOL(1e-10, 1e-5)));
		
		// For regular hexagon with side 1: circumradius = 1, inradius = sqrt(3)/2
		REQUIRE_THAT(hex.Circumradius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(hex.Inradius(), RealWithinAbs(std::sqrt(static_cast<MML::Real>(3.0)) / 2, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("RegularPolygon::EquilateralTriangle", "[RegularPolygon]")
	{
		TEST_PRECISION_INFO();
		RegularPolygon tri = RegularPolygon::EquilateralTriangle(static_cast<MML::Real>(2.0));
		
		REQUIRE(tri.NumSides() == 3);
		REQUIRE_THAT(tri.InteriorAngle(), RealWithinAbs(Constants::PI / 3, TOL(1e-10, 1e-5)));
		
		// Area = (sqrt(3)/4) * s² = sqrt(3)
		REQUIRE_THAT(tri.Area(), RealWithinAbs(std::sqrt(static_cast<MML::Real>(3.0)), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("RegularPolygon::FromCircumradius", "[RegularPolygon]")
	{
		TEST_PRECISION_INFO();
		// Square inscribed in circle of radius sqrt(2) has side 2
		RegularPolygon sq = RegularPolygon::FromCircumradius(4, std::sqrt(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(sq.SideLength(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("RegularPolygon::FromArea", "[RegularPolygon]")
	{
		TEST_PRECISION_INFO();
		// Unit square has area 1
		RegularPolygon sq = RegularPolygon::FromArea(4, static_cast<MML::Real>(1.0));
		REQUIRE_THAT(sq.SideLength(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("RegularPolygon::AllNamedPolygons", "[RegularPolygon]")
	{
		// Just verify they create valid polygons
		REQUIRE(RegularPolygon::Pentagon(static_cast<MML::Real>(1.0)).IsValid());
		REQUIRE(RegularPolygon::Heptagon(static_cast<MML::Real>(1.0)).IsValid());
		REQUIRE(RegularPolygon::Octagon(static_cast<MML::Real>(1.0)).IsValid());
		REQUIRE(RegularPolygon::Nonagon(static_cast<MML::Real>(1.0)).IsValid());
		REQUIRE(RegularPolygon::Decagon(static_cast<MML::Real>(1.0)).IsValid());
		REQUIRE(RegularPolygon::Dodecagon(static_cast<MML::Real>(1.0)).IsValid());
	}

	/*********************************************************************/
	/*****                  Rectangle tests                          *****/
	/*********************************************************************/
	TEST_CASE("Rectangle::DefaultConstructor", "[Rectangle]")
	{
		Rectangle r;
		REQUIRE(r.Width() == static_cast<MML::Real>(0.0));
		REQUIRE(r.Height() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(r.IsValid());
	}

	TEST_CASE("Rectangle::Properties", "[Rectangle]")
	{
		TEST_PRECISION_INFO();
		Rectangle r(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));
		
		REQUIRE(r.IsValid());
		REQUIRE_FALSE(r.IsSquare());
		REQUIRE_THAT(r.Area(), RealWithinAbs(static_cast<MML::Real>(12.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(r.Perimeter(), RealWithinAbs(static_cast<MML::Real>(14.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(r.Diagonal(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(r.AspectRatio(), RealWithinAbs(static_cast<MML::Real>(0.75), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Rectangle::Square", "[Rectangle]")
	{
		Rectangle sq = Rectangle::Square(static_cast<MML::Real>(5.0));
		REQUIRE(sq.IsSquare());
		REQUIRE_THAT(sq.Width(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sq.Height(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Rectangle::GoldenRectangle", "[Rectangle]")
	{
		TEST_PRECISION_INFO();
		Rectangle golden = Rectangle::GoldenRectangle(Constants::GoldenRatio);
		REQUIRE_THAT(golden.AspectRatio(), RealWithinAbs(Constants::GoldenRatio, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Rectangle::CircleRadii", "[Rectangle]")
	{
		TEST_PRECISION_INFO();
		Rectangle r(static_cast<MML::Real>(6.0), static_cast<MML::Real>(8.0));
		
		// Inscribed circle radius = min(w,h)/2 = 3
		REQUIRE_THAT(r.InscribedCircleRadius(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		
		// Circumscribed circle radius = diagonal/2 = 5
		REQUIRE_THAT(r.CircumscribedCircleRadius(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                Parallelogram tests                        *****/
	/*********************************************************************/
	TEST_CASE("Parallelogram::DefaultConstructor", "[Parallelogram]")
	{
		Parallelogram p;
		REQUIRE(p.SideA() == static_cast<MML::Real>(0.0));
		REQUIRE(p.SideB() == static_cast<MML::Real>(0.0));
		REQUIRE(p.Angle() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(p.IsValid());
	}

	TEST_CASE("Parallelogram::Properties", "[Parallelogram]")
	{
		TEST_PRECISION_INFO();
		Parallelogram p(static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0), Constants::PI / 3);  // 60°
		
		REQUIRE(p.IsValid());
		REQUIRE_THAT(p.Perimeter(), RealWithinAbs(static_cast<MML::Real>(14.0), TOL(1e-10, 1e-5)));
		
		// Area = a*b*sin(θ) = 4*3*sin(60°) = 12 * sqrt(3)/2 = 6*sqrt(3)
		REQUIRE_THAT(p.Area(), RealWithinAbs(static_cast<MML::Real>(6.0) * std::sqrt(static_cast<MML::Real>(3.0)), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Parallelogram::IsRectangle", "[Parallelogram]")
	{
		Parallelogram rect(static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0), Constants::PI / 2);  // 90°
		REQUIRE(rect.IsRectangle());
		REQUIRE_FALSE(rect.IsRhombus());
	}

	TEST_CASE("Parallelogram::IsRhombus", "[Parallelogram]")
	{
		Parallelogram rhombus(static_cast<MML::Real>(4.0), static_cast<MML::Real>(4.0), Constants::PI / 3);
		REQUIRE(rhombus.IsRhombus());
		REQUIRE_FALSE(rhombus.IsRectangle());
	}

	TEST_CASE("Parallelogram::IsSquare", "[Parallelogram]")
	{
		Parallelogram sq(static_cast<MML::Real>(4.0), static_cast<MML::Real>(4.0), Constants::PI / 2);
		REQUIRE(sq.IsSquare());
		REQUIRE(sq.IsRectangle());
		REQUIRE(sq.IsRhombus());
	}

	TEST_CASE("Parallelogram::Heights", "[Parallelogram]")
	{
		TEST_PRECISION_INFO();
		Parallelogram p(static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0), Constants::PI / 2);  // Rectangle
		
		REQUIRE_THAT(p.HeightToA(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p.HeightToB(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                    Rhombus tests                          *****/
	/*********************************************************************/
	TEST_CASE("Rhombus::DefaultConstructor", "[Rhombus]")
	{
		Rhombus r;
		REQUIRE(r.Side() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(r.IsValid());
	}

	TEST_CASE("Rhombus::FromDiagonals", "[Rhombus]")
	{
		TEST_PRECISION_INFO();
		Rhombus r(static_cast<MML::Real>(6.0), static_cast<MML::Real>(8.0));  // diagonals 6 and 8
		
		REQUIRE(r.IsValid());
		REQUIRE_THAT(r.Diagonal1(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(r.Diagonal2(), RealWithinAbs(static_cast<MML::Real>(8.0), TOL(1e-10, 1e-5)));
		
		// Side = sqrt((d1/2)² + (d2/2)²) = sqrt(9+16) = 5
		REQUIRE_THAT(r.Side(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		
		// Area = d1*d2/2 = 24
		REQUIRE_THAT(r.Area(), RealWithinAbs(static_cast<MML::Real>(24.0), TOL(1e-10, 1e-5)));
		
		// Perimeter = 4*side = 20
		REQUIRE_THAT(r.Perimeter(), RealWithinAbs(static_cast<MML::Real>(20.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Rhombus::IsSquare", "[Rhombus]")
	{
		Rhombus sq(static_cast<MML::Real>(4.0), static_cast<MML::Real>(4.0));  // Equal diagonals -> square
		REQUIRE(sq.IsSquare());
		
		Rhombus notSq(static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0));
		REQUIRE_FALSE(notSq.IsSquare());
	}

	TEST_CASE("Rhombus::FromSideAndAngle", "[Rhombus]")
	{
		TEST_PRECISION_INFO();
		Rhombus r = Rhombus::FromSideAndAngle(static_cast<MML::Real>(5.0), Constants::PI / 2);  // 90° -> square
		
		REQUIRE(r.IsSquare());
	}

	/*********************************************************************/
	/*****                   Trapezoid tests                         *****/
	/*********************************************************************/
	TEST_CASE("Trapezoid::DefaultConstructor", "[Trapezoid]")
	{
		Trapezoid t;
		REQUIRE(t.ParallelSideA() == static_cast<MML::Real>(0.0));
		REQUIRE(t.ParallelSideB() == static_cast<MML::Real>(0.0));
		REQUIRE(t.Height() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(t.IsValid());
	}

	TEST_CASE("Trapezoid::Properties", "[Trapezoid]")
	{
		TEST_PRECISION_INFO();
		Trapezoid t(static_cast<MML::Real>(6.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0));  // bases 6 and 4, height 3
		
		REQUIRE(t.IsValid());
		REQUIRE_FALSE(t.IsParallelogram());
		
		// Area = (a+b)*h/2 = 10*3/2 = 15
		REQUIRE_THAT(t.Area(), RealWithinAbs(static_cast<MML::Real>(15.0), TOL(1e-10, 1e-5)));
		
		// Median = (a+b)/2 = 5
		REQUIRE_THAT(t.Median(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Trapezoid::IsParallelogram", "[Trapezoid]")
	{
		Trapezoid para(static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(3.0));  // Equal parallel sides
		REQUIRE(para.IsParallelogram());
	}

	TEST_CASE("Trapezoid::Isosceles", "[Trapezoid]")
	{
		TEST_PRECISION_INFO();
		// Isosceles trapezoid with bases 4 and 10, legs of length 5
		Trapezoid t = Trapezoid::Isosceles(static_cast<MML::Real>(10.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0));
		
		// Height = sqrt(5² - 3²) = sqrt(25-9) = 4
		REQUIRE_THAT(t.Height(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                  SphereGeom tests                         *****/
	/*********************************************************************/
	TEST_CASE("SphereGeom::DefaultConstructor", "[SphereGeom]")
	{
		SphereGeom s;
		REQUIRE(s.Radius() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(s.IsValid());
	}

	TEST_CASE("SphereGeom::Properties", "[SphereGeom]")
	{
		TEST_PRECISION_INFO();
		SphereGeom s(static_cast<MML::Real>(3.0));
		
		REQUIRE(s.IsValid());
		REQUIRE_THAT(s.Diameter(), RealWithinAbs(static_cast<MML::Real>(6.0), TOL(1e-10, 1e-5)));
		
		// Volume = (4/3)πr³ = (4/3)π*27 = 36π
		REQUIRE_THAT(s.Volume(), RealWithinAbs(static_cast<MML::Real>(36.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Surface area = 4πr² = 36π
		REQUIRE_THAT(s.SurfaceArea(), RealWithinAbs(static_cast<MML::Real>(36.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("SphereGeom::GreatCircle", "[SphereGeom]")
	{
		SphereGeom s(static_cast<MML::Real>(5.0));
		Circle great = s.GreatCircle();
		REQUIRE_THAT(great.Radius(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("SphereGeom::Cap", "[SphereGeom]")
	{
		TEST_PRECISION_INFO();
		SphereGeom s(static_cast<MML::Real>(2.0));
		
		// Hemisphere: h = r = 2
		REQUIRE_THAT(s.CapArea(static_cast<MML::Real>(2.0)), RealWithinAbs(static_cast<MML::Real>(8.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Cap volume for hemisphere: (2/3)πr³
		Real expectedVol = (static_cast<MML::Real>(2.0) / static_cast<MML::Real>(3.0)) * Constants::PI * static_cast<MML::Real>(8.0);
		REQUIRE_THAT(s.CapVolume(static_cast<MML::Real>(2.0)), RealWithinAbs(expectedVol, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("SphereGeom::FactoryMethods", "[SphereGeom]")
	{
		TEST_PRECISION_INFO();
		
		// From volume = (4/3)π -> r = 1
		SphereGeom fromVol = SphereGeom::FromVolume((static_cast<MML::Real>(4.0) / static_cast<MML::Real>(3.0)) * Constants::PI);
		REQUIRE_THAT(fromVol.Radius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
		
		// From surface area = 4π -> r = 1
		SphereGeom fromArea = SphereGeom::FromSurfaceArea(static_cast<MML::Real>(4.0) * Constants::PI);
		REQUIRE_THAT(fromArea.Radius(), RealWithinAbs(static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                 CylinderGeom tests                        *****/
	/*********************************************************************/
	TEST_CASE("CylinderGeom::DefaultConstructor", "[CylinderGeom]")
	{
		CylinderGeom c;
		REQUIRE(c.Radius() == static_cast<MML::Real>(0.0));
		REQUIRE(c.Height() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(c.IsValid());
	}

	TEST_CASE("CylinderGeom::Properties", "[CylinderGeom]")
	{
		TEST_PRECISION_INFO();
		CylinderGeom c(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
		
		REQUIRE(c.IsValid());
		REQUIRE_THAT(c.Diameter(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		
		// Volume = πr²h = 4π*5 = 20π
		REQUIRE_THAT(c.Volume(), RealWithinAbs(static_cast<MML::Real>(20.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Lateral area = 2πrh = 4π*5 = 20π
		REQUIRE_THAT(c.LateralArea(), RealWithinAbs(static_cast<MML::Real>(20.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Base area = πr² = 4π
		REQUIRE_THAT(c.BaseArea(), RealWithinAbs(static_cast<MML::Real>(4.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Surface area = lateral + 2*base = 20π + 8π = 28π
		REQUIRE_THAT(c.SurfaceArea(), RealWithinAbs(static_cast<MML::Real>(28.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("CylinderGeom::FromVolume", "[CylinderGeom]")
	{
		TEST_PRECISION_INFO();
		CylinderGeom c = CylinderGeom::FromVolume(static_cast<MML::Real>(20.0) * Constants::PI, static_cast<MML::Real>(2.0));
		
		REQUIRE_THAT(c.Radius(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Height(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                   ConeGeom tests                          *****/
	/*********************************************************************/
	TEST_CASE("ConeGeom::DefaultConstructor", "[ConeGeom]")
	{
		ConeGeom c;
		REQUIRE(c.Radius() == static_cast<MML::Real>(0.0));
		REQUIRE(c.Height() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(c.IsValid());
	}

	TEST_CASE("ConeGeom::Properties", "[ConeGeom]")
	{
		TEST_PRECISION_INFO();
		ConeGeom c(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0));  // 3-4-5 right triangle profile
		
		REQUIRE(c.IsValid());
		REQUIRE_THAT(c.SlantHeight(), RealWithinAbs(static_cast<MML::Real>(5.0), TOL(1e-10, 1e-5)));
		
		// Volume = (1/3)πr²h = (1/3)*9π*4 = 12π
		REQUIRE_THAT(c.Volume(), RealWithinAbs(static_cast<MML::Real>(12.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Lateral area = πrs = 3π*5 = 15π
		REQUIRE_THAT(c.LateralArea(), RealWithinAbs(static_cast<MML::Real>(15.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Base area = πr² = 9π
		REQUIRE_THAT(c.BaseArea(), RealWithinAbs(static_cast<MML::Real>(9.0) * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Surface area = 15π + 9π = 24π
		REQUIRE_THAT(c.SurfaceArea(), RealWithinAbs(static_cast<MML::Real>(24.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("ConeGeom::FromSlantHeight", "[ConeGeom]")
	{
		TEST_PRECISION_INFO();
		ConeGeom c = ConeGeom::FromSlantHeight(static_cast<MML::Real>(3.0), static_cast<MML::Real>(5.0));  // r=3, slant=5 -> h=4
		
		REQUIRE_THAT(c.Radius(), RealWithinAbs(static_cast<MML::Real>(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(c.Height(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                    Frustum tests                          *****/
	/*********************************************************************/
	TEST_CASE("Frustum::DefaultConstructor", "[Frustum]")
	{
		Frustum f;
		REQUIRE(f.BottomRadius() == static_cast<MML::Real>(0.0));
		REQUIRE(f.TopRadius() == static_cast<MML::Real>(0.0));
		REQUIRE(f.Height() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(f.IsValid());
	}

	TEST_CASE("Frustum::Properties", "[Frustum]")
	{
		TEST_PRECISION_INFO();
		Frustum f(static_cast<MML::Real>(4.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));  // r1=4, r2=2, h=3
		
		REQUIRE(f.IsValid());
		
		// Slant height = sqrt((r1-r2)² + h²) = sqrt(4+9) = sqrt(13)
		REQUIRE_THAT(f.SlantHeight(), RealWithinAbs(std::sqrt(static_cast<MML::Real>(13.0)), TOL(1e-10, 1e-5)));
		
		// Volume = (π*h/3)(r1² + r1*r2 + r2²) = (π*3/3)(16+8+4) = π*28
		REQUIRE_THAT(f.Volume(), RealWithinAbs(static_cast<MML::Real>(28.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Frustum::Bases", "[Frustum]")
	{
		Frustum f(static_cast<MML::Real>(4.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		
		Circle bottom = f.BottomBase();
		Circle top = f.TopBase();
		
		REQUIRE_THAT(bottom.Radius(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(top.Radius(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                  Tetrahedron tests                        *****/
	/*********************************************************************/
	TEST_CASE("Tetrahedron::DefaultConstructor", "[Tetrahedron]")
	{
		Tetrahedron t;
		REQUIRE(t.AB() == static_cast<MML::Real>(0.0));
	}

	TEST_CASE("Tetrahedron::Regular", "[Tetrahedron]")
	{
		TEST_PRECISION_INFO();
		Tetrahedron t = Tetrahedron::Regular(static_cast<MML::Real>(2.0));
		
		REQUIRE(t.IsRegular());
		REQUIRE_THAT(t.AB(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		
		// Regular tetrahedron with edge a: V = a³/(6*sqrt(2))
		Real expectedVol = static_cast<MML::Real>(8.0) / (static_cast<MML::Real>(6.0) * std::sqrt(static_cast<MML::Real>(2.0)));
		REQUIRE_THAT(t.Volume(), RealWithinAbs(expectedVol, TOL(1e-8, 1e-4)));
	}

	TEST_CASE("Tetrahedron::Faces", "[Tetrahedron]")
	{
		Tetrahedron t = Tetrahedron::Regular(static_cast<MML::Real>(2.0));
		
		Triangle abc = t.FaceABC();
		Triangle abd = t.FaceABD();
		Triangle acd = t.FaceACD();
		Triangle bcd = t.FaceBCD();
		
		// All faces should be equilateral
		REQUIRE(abc.IsEquilateral());
		REQUIRE(abd.IsEquilateral());
		REQUIRE(acd.IsEquilateral());
		REQUIRE(bcd.IsEquilateral());
	}

	/*********************************************************************/
	/*****                   Spheroid tests                          *****/
	/*********************************************************************/
	TEST_CASE("Spheroid::DefaultConstructor", "[Spheroid]")
	{
		Spheroid s;
		REQUIRE(s.EquatorialRadius() == static_cast<MML::Real>(0.0));
		REQUIRE(s.PolarRadius() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(s.IsValid());
	}

	TEST_CASE("Spheroid::Classification", "[Spheroid]")
	{
		Spheroid sphere(static_cast<MML::Real>(5.0), static_cast<MML::Real>(5.0));
		REQUIRE(sphere.IsSphere());
		REQUIRE_FALSE(sphere.IsOblate());
		REQUIRE_FALSE(sphere.IsProlate());
		
		Spheroid oblate(static_cast<MML::Real>(6.0), static_cast<MML::Real>(4.0));  // Earth-like
		REQUIRE(oblate.IsOblate());
		REQUIRE_FALSE(oblate.IsProlate());
		
		Spheroid prolate(static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0));  // Rugby ball-like
		REQUIRE(prolate.IsProlate());
		REQUIRE_FALSE(prolate.IsOblate());
	}

	TEST_CASE("Spheroid::Volume", "[Spheroid]")
	{
		TEST_PRECISION_INFO();
		Spheroid s(static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0));  // a=3, c=2
		
		// Volume = (4/3)π*a²*c = (4/3)π*9*2 = 24π
		REQUIRE_THAT(s.Volume(), RealWithinAbs(static_cast<MML::Real>(24.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Spheroid::Flattening", "[Spheroid]")
	{
		TEST_PRECISION_INFO();
		Spheroid s(static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0));
		
		// Flattening = (a-c)/a = (4-3)/4 = 0.25
		REQUIRE_THAT(s.Flattening(), RealWithinAbs(static_cast<MML::Real>(0.25), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Spheroid::Earth", "[Spheroid]")
	{
		Spheroid earth = Spheroid::Earth();
		
		REQUIRE(earth.IsValid());
		REQUIRE(earth.IsOblate());
		REQUIRE_THAT(earth.EquatorialRadius(), RealWithinAbs(static_cast<MML::Real>(6378137.0), static_cast<MML::Real>(1.0)));
	}

	/*********************************************************************/
	/*****                  TorusGeom tests                          *****/
	/*********************************************************************/
	TEST_CASE("TorusGeom::DefaultConstructor", "[TorusGeom]")
	{
		TorusGeom t;
		REQUIRE(t.MajorRadius() == static_cast<MML::Real>(0.0));
		REQUIRE(t.MinorRadius() == static_cast<MML::Real>(0.0));
		REQUIRE_FALSE(t.IsValid());
	}

	TEST_CASE("TorusGeom::Properties", "[TorusGeom]")
	{
		TEST_PRECISION_INFO();
		TorusGeom t(static_cast<MML::Real>(3.0), static_cast<MML::Real>(1.0));  // R=3, r=1
		
		REQUIRE(t.IsValid());
		REQUIRE(t.IsRingTorus());
		REQUIRE_FALSE(t.IsHornTorus());
		REQUIRE_FALSE(t.IsSpindleTorus());
		
		// Volume = 2π²Rr² = 2π²*3*1 = 6π²
		REQUIRE_THAT(t.Volume(), RealWithinAbs(static_cast<MML::Real>(6.0) * Constants::PI * Constants::PI, TOL(1e-10, 1e-5)));
		
		// Surface area = 4π²Rr = 4π²*3*1 = 12π²
		REQUIRE_THAT(t.SurfaceArea(), RealWithinAbs(static_cast<MML::Real>(12.0) * Constants::PI * Constants::PI, TOL(1e-10, 1e-5)));
		
		REQUIRE_THAT(t.InnerRadius(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(t.OuterRadius(), RealWithinAbs(static_cast<MML::Real>(4.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("TorusGeom::Types", "[TorusGeom]")
	{
		TorusGeom ring(static_cast<MML::Real>(5.0), static_cast<MML::Real>(2.0));
		REQUIRE(ring.IsRingTorus());
		
		TorusGeom horn(static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0));
		REQUIRE(horn.IsHornTorus());
		
		TorusGeom spindle(static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0));
		REQUIRE(spindle.IsSpindleTorus());
	}

	TEST_CASE("TorusGeom::TubeCrossSection", "[TorusGeom]")
	{
		TorusGeom t(static_cast<MML::Real>(5.0), static_cast<MML::Real>(2.0));
		Circle tube = t.TubeCrossSection();
		REQUIRE_THAT(tube.Radius(), RealWithinAbs(static_cast<MML::Real>(2.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                 Type alias tests                          *****/
	/*********************************************************************/
	TEST_CASE("Geometry::TypeAliases", "[Geometry]")
	{
		// Verify type aliases work
		Pnt2Cart p2c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0));
		REQUIRE(p2c.X() == static_cast<MML::Real>(1.0));
		
		Pnt2Pol p2p(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5));
		REQUIRE(p2p.R() == static_cast<MML::Real>(1.0));
		
		Pnt3Cart p3c(static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0));
		REQUIRE(p3c.Z() == static_cast<MML::Real>(3.0));
		
		Pnt3Sph p3s(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(1.0));
		REQUIRE(p3s.Theta() == static_cast<MML::Real>(0.5));
		
		Pnt3Cyl p3cy(static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.5), static_cast<MML::Real>(2.0));
		REQUIRE(p3cy.Z() == static_cast<MML::Real>(2.0));
	}

} // namespace MML::Tests::Base::GeometryTests
