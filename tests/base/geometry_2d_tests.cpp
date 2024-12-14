#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry2D.h"
#endif

using namespace MML;

// TODO 0.9 - BIG, EMPTY!!!
namespace MML::Tests::Base::Geometry2DTests
{
	/*********************************************************************/
	/*****                       Line2D tests                        *****/
	/*********************************************************************/
	TEST_CASE("Line2D::Init", "[Line2D]") {
		Point2Cartesian    start_pnt(1, 1);
		Vector2Cartesian   direction(1, 0);

		Line2D     line(start_pnt, direction);

		REQUIRE(line.StartPoint().X() == 1.0);
		REQUIRE(line.StartPoint().X() == 1.0);
		REQUIRE(line.Direction().X() == 1.0);
		REQUIRE(line.Direction().Y() == 0.0);

		Point2Cartesian    pnt1(1, 1);
		Point2Cartesian    pnt2(2, 2);

		Line2D     line2(pnt1, pnt2);

	}
	TEST_CASE("Line2D::PointOnLine", "[Line2D]")
	{
	}

	/*********************************************************************/
	/*****                   SegmentLine2D tests                     *****/
	/*********************************************************************/
	TEST_CASE("SegmentLine2D::Init", "[SegmentLine2D]")
	{
	}
	TEST_CASE("SegmentLine2D::PointOnSegment", "[SegmentLine2D]")
	{
	}
	TEST_CASE("SegmentLine2D::Length", "[SegmentLine2D]")
	{
	}
	TEST_CASE("SegmentLine2D::Direction", "[SegmentLine2D]")
	{
	}

	/*********************************************************************/
	/*****                     Triangle2D tests                      *****/
	/*********************************************************************/
	TEST_CASE("Triangle2D::Init", "[Triangle2D]")
	{
	}
	TEST_CASE("Triangle2D::Area", "[Triangle2D]")
	{
	}

	/*********************************************************************/
	/*****                     Polygon2D tests                       *****/
	/*********************************************************************/
	TEST_CASE("Polygon2D::Init", "[Polygon2D]")
	{
	}
	TEST_CASE("Polygon2D::Area", "[Polygon2D]")
	{
	}
	TEST_CASE("Polygon2D::IsSimple", "[Polygon2D]")
	{
	}
	TEST_CASE("Polygon2D::IsConvex", "[Polygon2D]")
	{
	}
	TEST_CASE("Polygon2D::IsInside", "[Polygon2D]")
	{
	}
}