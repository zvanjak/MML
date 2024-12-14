#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry.h"
#include "base/VectorN.h"
#endif

using namespace MML;

namespace MML::Tests::Base::GeometryTests
{
	/*********************************************************************/
	/*****                 Point2Cartesian tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point2Cartesian::Init", "[Point2Cartesian]") {
		Point2Cartesian a(1, 1);

		REQUIRE(a.X() == 1.0);
		REQUIRE(a.Y() == 1.0);
	}
	TEST_CASE("Point2Cartesian::Dist", "[Point2Cartesian]") 
	{
	}
	TEST_CASE("Point2Cartesian::op+", "[Point2Cartesian]") 
	{
	}
	TEST_CASE("Point2Cartesian::op*", "[Point2Cartesian]") 
	{
	}
	TEST_CASE("Point2Cartesian::op/", "[Point2Cartesian]") 
	{
	}

	/*********************************************************************/
	/*****                 Point2Polar tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point2Polar::Init", "[Point2Polar]") 
	{
	}
	TEST_CASE("Point2Polar::Dist", "[Point2Polar]") 
	{
	}
	TEST_CASE("Point2Polar::GetFromCartesian", "[Point2Polar]") 
	{
	}
	TEST_CASE("Point2Polar::TransfToCartesian", "[Point2Polar]") 
	{
	}

	/*********************************************************************/
	/*****                 Point3Cartesian tests                     *****/
	/*********************************************************************/
	TEST_CASE("Point3Cartesian::Init", "[Point3Cartesian]") 
	{
	}
	TEST_CASE("Point3Cartesian::Dist", "[Point3Cartesian]") 
	{
	}
	TEST_CASE("Point3Cartesian::op+", "[Point3Cartesian]") 
	{
	}
	TEST_CASE("Point3Cartesian::op*", "[Point3Cartesian]") 
	{
	}
	TEST_CASE("Point3Cartesian::op/", "[Point3Cartesian]") 
	{
	}

	/*********************************************************************/
	/*****                     Triangle tests                        *****/
	/*********************************************************************/
	TEST_CASE("Triangle::Init", "[Triangle]")
	{
	}
	TEST_CASE("Triangle::Area", "[Triangle]")
	{
	}
	TEST_CASE("Triangle::IsRight", "[Triangle]")
	{
	}
	TEST_CASE("Triangle::IsIsosceles", "[Triangle]")
	{
	}
	TEST_CASE("Triangle::IsEquilateral", "[Triangle]")
	{
	}
} // namespace MML::Tests::Base::GeometryTests