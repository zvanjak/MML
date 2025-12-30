#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry2D.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
// using statement removed - RealWithinRel is a function, not a type

namespace MML::Tests::Base::Geometry2DTests
{
	/*********************************************************************/
	/*****                       Line2D tests                        *****/
	/*********************************************************************/
	TEST_CASE("Line2D::Init", "[Line2D]") {
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Point2Cartesian    start_pnt(1, 1);
		Vector2Cartesian   direction(1, 0);

		Line2D     line(start_pnt, direction);

		REQUIRE(line.StartPoint().X() == REAL(REAL(1.0)));
		REQUIRE(line.StartPoint().X() == REAL(REAL(1.0)));
		REQUIRE(line.Direction().X() == REAL(REAL(1.0)));
		REQUIRE(line.Direction().Y() == REAL(REAL(0.0)));

		Point2Cartesian    pnt1(1, 1);
		Point2Cartesian    pnt2(2, 2);

		Line2D     line2(pnt1, pnt2);

		REQUIRE(line2.StartPoint().X() == REAL(REAL(1.0)));
		REQUIRE(line2.StartPoint().Y() == REAL(REAL(1.0)));
		REQUIRE_THAT(line2.Direction().X() , RealWithinRel(REAL(REAL(0.7071067811865475)), 1e-5));
		REQUIRE_THAT(line2.Direction().Y() , RealWithinRel(REAL(REAL(0.7071067811865475)), 1e-5));
	}

	/*********************************************************************/
	/*****                   SegmentLine2D tests                     *****/
	/*********************************************************************/
	TEST_CASE("SegmentLine2D::Init", "[SegmentLine2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(REAL(1.0)), REAL(REAL(2.0)));
		Point2Cartesian p2(REAL(REAL(4.0)), REAL(REAL(6.0)));
		SegmentLine2D seg(p1, p2);
		REQUIRE_THAT(seg.StartPoint().X() , RealWithinRel(REAL(REAL(1.0)), 1e-5));
		REQUIRE_THAT(seg.StartPoint().Y() , RealWithinRel(REAL(REAL(2.0)), 1e-5));
		REQUIRE_THAT(seg.EndPoint().X() , RealWithinRel(REAL(REAL(4.0)), 1e-5));
		REQUIRE_THAT(seg.EndPoint().Y() , RealWithinRel(REAL(REAL(6.0)), 1e-5));
	}

	TEST_CASE("SegmentLine2D::PointOnSegment", "[SegmentLine2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(REAL(0.0)), REAL(REAL(0.0)));
		Point2Cartesian p2(REAL(REAL(2.0)), REAL(REAL(2.0)));
		SegmentLine2D seg(p1, p2);

		// Start and end points should be on the segment
		REQUIRE(seg.PointOnSegment(REAL(REAL(0.0))) == p1);
		REQUIRE(seg.PointOnSegment(REAL(REAL(1.0))) == p2);

		// Midpoint
		Point2Cartesian mid = seg.PointOnSegment(REAL(REAL(0.5)));
		REQUIRE_THAT(mid.X() , RealWithinRel(REAL(REAL(1.0)), 1e-5));
		REQUIRE_THAT(mid.Y() , RealWithinRel(REAL(REAL(1.0)), 1e-5));
	}

	TEST_CASE("SegmentLine2D::Length", "[SegmentLine2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(REAL(1.0)), REAL(REAL(1.0)));
		Point2Cartesian p2(REAL(REAL(4.0)), REAL(REAL(5.0)));
		SegmentLine2D seg(p1, p2);
		Real expected_length = sqrt(pow(REAL(REAL(4.0)) - REAL(REAL(1.0)), 2) + pow(REAL(REAL(5.0)) - REAL(REAL(1.0)), 2));
		REQUIRE_THAT(seg.StartPoint().Dist(seg.EndPoint()) , RealWithinRel(expected_length, REAL(1e-5)));
	}

	TEST_CASE("SegmentLine2D::Direction", "[SegmentLine2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(REAL(0.0)), REAL(REAL(0.0)));
		Point2Cartesian p2(REAL(REAL(3.0)), REAL(REAL(4.0)));
		SegmentLine2D seg(p1, p2);
		Vector2Cartesian dir(p1, p2);
		// Direction should be (3,4) normalized
		Real norm = sqrt(REAL(REAL(3.0)) * REAL(REAL(3.0)) + REAL(REAL(4.0)) * REAL(REAL(4.0)));
		REQUIRE_THAT(dir.X() / norm , RealWithinRel(REAL(REAL(0.6)), 1e-5));
		REQUIRE_THAT(dir.Y() / norm , RealWithinRel(REAL(REAL(0.8)), 1e-5));
	}

	/*********************************************************************/
	/*****                     Triangle2D tests                      *****/
	/*********************************************************************/
	TEST_CASE("Triangle2D::Init", "[Triangle2D]")
	{
			TEST_PRECISION_INFO();
		TEST_PRECISION_INFO();
	}
  TEST_CASE("Triangle2D::Area", "[Triangle2D]")
  {
    
    Point2Cartesian p1(0, 0);
    Point2Cartesian p2(4, 0);
    Point2Cartesian p3(0, 3);
    Triangle2D tri(p1, p2, p3);

    // Assume Triangle2D has Area() method
    Real area = tri.Area();
    REQUIRE_THAT(area , RealWithinRel(REAL(REAL(6.0)), 1e-5));

    // Degenerate triangle (collinear)
    Point2Cartesian p4(1, 1);
    Point2Cartesian p5(2, 2);
    Point2Cartesian p6(3, 3);
    Triangle2D degenerate(p4, p5, p6);
    Real deg_area = degenerate.Area();
    REQUIRE_THAT(deg_area , RealWithinRel(REAL(REAL(0.0)), 1e-5));
  }

	/*********************************************************************/
	/*****                     Polygon2D tests                       *****/
	/*********************************************************************/
	TEST_CASE("Polygon2D::Init", "[Polygon2D]")
	{
			TEST_PRECISION_INFO();
		TEST_PRECISION_INFO();
	}
	TEST_CASE("Polygon2D::Area1", "[Polygon2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Square with side 2, area should be 4
		std::vector<Point2Cartesian> square = {
			{0, 0}, {2, 0}, {2, 2}, {0, 2}
		};
		Polygon2D poly(square);
		REQUIRE_THAT(poly.Area() , RealWithinRel(REAL(REAL(4.0)), 1e-5));
	}

	TEST_CASE("Polygon2D::Area2", "[Polygon2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Triangle as polygon, area should be 6
		std::vector<Point2Cartesian> triangle = {
			{0, 0}, {4, 0}, {0, 3}
		};
		Polygon2D poly(triangle);
		REQUIRE_THAT(poly.Area() , RealWithinRel(REAL(REAL(6.0)), 1e-5));
	}

	TEST_CASE("Polygon2D::Area3", "[Polygon2D]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Degenerate polygon (collinear points), area should be 0
		std::vector<Point2Cartesian> line = {
			{1, 1}, {2, 2}, {3, 3}
		};
		Polygon2D poly(line);
		REQUIRE_THAT(poly.Area() , RealWithinRel(REAL(REAL(0.0)), 1e-5));

		// Regular pentagon, side length 1, area � REAL(REAL(2.3776))
		std::vector<Point2Cartesian> pentagon = {
			{REAL(REAL(0.0)), REAL(REAL(1.0))},
			{REAL(REAL(0.9510565163)), REAL(REAL(0.309017))},
			{REAL(REAL(0.5877852523)), -REAL(REAL(0.809017))},
			{-REAL(REAL(0.5877852523)), -REAL(REAL(0.809017))},
			{-REAL(REAL(0.9510565163)), REAL(REAL(0.309017))}
		};
		Polygon2D poly2(pentagon);
		REQUIRE_THAT(poly2.Area() , RealWithinRel(REAL(REAL(2.3776)), 1e-4));
	}
	TEST_CASE("Polygon2D::IsSimple", "[Polygon2D]")
	{
			TEST_PRECISION_INFO();
		TEST_PRECISION_INFO();
	}
	TEST_CASE("Polygon2D::IsConvex", "[Polygon2D]")
	{
			TEST_PRECISION_INFO();
		TEST_PRECISION_INFO();
	}
	TEST_CASE("Polygon2D::IsInside", "[Polygon2D]")
	{
			TEST_PRECISION_INFO();
		TEST_PRECISION_INFO();
	}

	/*********************************************************************/
	/*****          Line2D Enhanced Functionality Tests              *****/
	/*********************************************************************/
	TEST_CASE("Line2D - Distance to point", "[Line2D][Distance]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Horizontal line at y=2
		Point2Cartesian p1(0, 2);
		Vector2Cartesian dir(1, 0);
		Line2D line(p1, dir);

		// Point above line
		Point2Cartesian p(5, 5);
		REQUIRE_THAT(line.DistanceToPoint(p), WithinAbs(REAL(REAL(3.0)), 1e-10));

		// Point on line
		Point2Cartesian pOn(10, 2);
		REQUIRE_THAT(line.DistanceToPoint(pOn), WithinAbs(REAL(REAL(0.0)), 1e-10));
	}

	TEST_CASE("Line2D - Closest point", "[Line2D][ClosestPoint]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Line through origin at 45 degrees
		Point2Cartesian origin(0, 0);
		Vector2Cartesian dir(1, 1);
		Line2D line(origin, dir);

		Point2Cartesian p(2, 0);
		Point2Cartesian closest = line.ClosestPoint(p);
		
		REQUIRE_THAT(closest.X(), WithinAbs(REAL(REAL(1.0)), 1e-10));
		REQUIRE_THAT(closest.Y(), WithinAbs(REAL(REAL(1.0)), 1e-10));
	}

	TEST_CASE("Line2D - Perpendicular line", "[Line2D][Perpendicular]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Horizontal line
		Line2D horizontal(Point2Cartesian(0, 0), Vector2Cartesian(1, 0));
		
		Point2Cartesian through(5, 3);
		Line2D perp = horizontal.Perpendicular(through);
		
		// Perpendicular should be vertical (direction 0,1 or 0,-1)
		REQUIRE_THAT(std::abs(perp.Direction().X()), WithinAbs(REAL(REAL(0.0)), 1e-10));
		REQUIRE_THAT(std::abs(perp.Direction().Y()), WithinAbs(REAL(REAL(1.0)), 1e-10));
		REQUIRE(perp.Contains(through));
	}

	TEST_CASE("Line2D - Line intersection", "[Line2D][Intersection]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Two intersecting lines
		Line2D line1(Point2Cartesian(0, 0), Vector2Cartesian(1, 0)); // Horizontal
		Line2D line2(Point2Cartesian(0, 0), Vector2Cartesian(0, 1)); // Vertical

		Point2Cartesian intersection;
		REQUIRE(line1.Intersects(line2, &intersection));
		REQUIRE_THAT(intersection.X(), WithinAbs(REAL(REAL(0.0)), 1e-10));
		REQUIRE_THAT(intersection.Y(), WithinAbs(REAL(REAL(0.0)), 1e-10));

		// Parallel lines
		Line2D line3(Point2Cartesian(0, 1), Vector2Cartesian(1, 0));
		REQUIRE_FALSE(line1.Intersects(line3));
	}

	/*********************************************************************/
	/*****        SegmentLine2D Enhanced Functionality Tests         *****/
	/*********************************************************************/
	TEST_CASE("SegmentLine2D - Midpoint", "[SegmentLine2D][Midpoint]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		SegmentLine2D seg(Point2Cartesian(2, 3), Point2Cartesian(8, 7));
		Point2Cartesian mid = seg.Midpoint();
		
		REQUIRE_THAT(mid.X(), WithinAbs(REAL(REAL(5.0)), 1e-10));
		REQUIRE_THAT(mid.Y(), WithinAbs(REAL(REAL(5.0)), 1e-10));
	}

	TEST_CASE("SegmentLine2D - Distance to point", "[SegmentLine2D][Distance]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		SegmentLine2D seg(Point2Cartesian(0, 0), Point2Cartesian(4, 0));

		// Point above segment midpoint
		Point2Cartesian p1(2, 3);
		REQUIRE_THAT(seg.DistanceToPoint(p1), WithinAbs(REAL(REAL(3.0)), 1e-10));

		// Point beyond segment end
		Point2Cartesian p2(6, 0);
		REQUIRE_THAT(seg.DistanceToPoint(p2), WithinAbs(REAL(REAL(2.0)), 1e-10));

		// Point before segment start
		Point2Cartesian p3(-2, 0);
		REQUIRE_THAT(seg.DistanceToPoint(p3), WithinAbs(REAL(REAL(2.0)), 1e-10));
	}

	TEST_CASE("SegmentLine2D - Closest point", "[SegmentLine2D][ClosestPoint]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		SegmentLine2D seg(Point2Cartesian(0, 0), Point2Cartesian(10, 0));

		// Point above segment
		Point2Cartesian p1(5, 5);
		Point2Cartesian closest1 = seg.ClosestPoint(p1);
		REQUIRE_THAT(closest1.X(), WithinAbs(REAL(REAL(5.0)), 1e-10));
		REQUIRE_THAT(closest1.Y(), WithinAbs(REAL(REAL(0.0)), 1e-10));

		// Point beyond end
		Point2Cartesian p2(15, 3);
		Point2Cartesian closest2 = seg.ClosestPoint(p2);
		REQUIRE_THAT(closest2.X(), WithinAbs(REAL(REAL(10.0)), 1e-10));
		REQUIRE_THAT(closest2.Y(), WithinAbs(REAL(REAL(0.0)), 1e-10));
	}

	TEST_CASE("SegmentLine2D - Segment intersection", "[SegmentLine2D][Intersection]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Intersecting segments
		SegmentLine2D seg1(Point2Cartesian(0, 0), Point2Cartesian(4, 4));
		SegmentLine2D seg2(Point2Cartesian(0, 4), Point2Cartesian(4, 0));

		Point2Cartesian intersection;
		REQUIRE(seg1.Intersects(seg2, &intersection));
		REQUIRE_THAT(intersection.X(), WithinAbs(REAL(REAL(2.0)), 1e-10));
		REQUIRE_THAT(intersection.Y(), WithinAbs(REAL(REAL(2.0)), 1e-10));

		// Non-intersecting segments
		SegmentLine2D seg3(Point2Cartesian(0, 0), Point2Cartesian(1, 0));
		SegmentLine2D seg4(Point2Cartesian(2, 0), Point2Cartesian(3, 0));
		REQUIRE_FALSE(seg3.Intersects(seg4));

		// Parallel segments
		SegmentLine2D seg5(Point2Cartesian(0, 0), Point2Cartesian(2, 0));
		SegmentLine2D seg6(Point2Cartesian(0, 1), Point2Cartesian(2, 1));
		REQUIRE_FALSE(seg5.Intersects(seg6));
	}

	/*********************************************************************/
	/*****        Triangle2D Enhanced Functionality Tests            *****/
	/*********************************************************************/
	TEST_CASE("Triangle2D - Perimeter", "[Triangle2D][Perimeter]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// 3-4-5 right triangle
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		REQUIRE_THAT(tri.Perimeter(), WithinAbs(REAL(REAL(12.0)), 1e-10));
	}

	TEST_CASE("Triangle2D - Centroid", "[Triangle2D][Centroid]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(6, 0), Point2Cartesian(0, 6));
		Point2Cartesian centroid = tri.Centroid();
		
		REQUIRE_THAT(centroid.X(), WithinAbs(REAL(REAL(2.0)), 1e-10));
		REQUIRE_THAT(centroid.Y(), WithinAbs(REAL(REAL(2.0)), 1e-10));
	}

	TEST_CASE("Triangle2D - Circumcenter and CircumRadius", "[Triangle2D][Circumcenter]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Right triangle with hypotenuse as diameter of circumcircle
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(4, 0), Point2Cartesian(0, 3));
		
		Point2Cartesian circumcenter = tri.Circumcenter();
		REQUIRE_THAT(circumcenter.X(), WithinAbs(REAL(REAL(2.0)), 1e-10));
		REQUIRE_THAT(circumcenter.Y(), WithinAbs(REAL(REAL(1.5)), 1e-10));

		Real radius = tri.CircumRadius();
		REQUIRE_THAT(radius, WithinAbs(REAL(REAL(2.5)), 1e-10)); // 5/2 (hypotenuse/2)
	}

	TEST_CASE("Triangle2D - Incenter and InRadius", "[Triangle2D][Incenter]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// 3-4-5 right triangle
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		
		Point2Cartesian incenter = tri.Incenter();
		Real inradius = tri.InRadius();
		
		REQUIRE_THAT(inradius, WithinAbs(REAL(REAL(1.0)), 1e-10)); // Area / s = 6 / 6 = 1
		
		// Incenter should be at distance inradius from all sides
		REQUIRE_THAT(incenter.X(), WithinAbs(REAL(REAL(1.0)), 1e-10));
		REQUIRE_THAT(incenter.Y(), WithinAbs(REAL(REAL(1.0)), 1e-10));
	}

	TEST_CASE("Triangle2D - Angles", "[Triangle2D][Angles]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// 3-4-5 right triangle with right angle at origin (0,0)
		// Vertices: pnt1=(0,0), pnt2=(3,0), pnt3=(0,4)
		// Side a: pnt1->pnt2 = 3 (base along x-axis)
		// Side b: pnt2->pnt3 = 5 (hypotenuse)
		// Side c: pnt3->pnt1 = 4 (vertical side)
		// AngleA = angle at pnt3 (opposite to side a)
		// AngleB = angle at pnt1 (opposite to side b) - THIS is 90 degrees!
		// AngleC = angle at pnt2 (opposite to side c)
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		
		Real angleA = tri.AngleA(); // Angle at pnt3 (0,4)
		Real angleB = tri.AngleB(); // Angle at pnt1 (0,0) - the right angle
		Real angleC = tri.AngleC(); // Angle at pnt2 (3,0)
		
		// AngleB is at the origin where perpendicular sides meet - should be 90 degrees
		REQUIRE_THAT(angleB, WithinAbs(Constants::PI / REAL(REAL(2.0)), 1e-8)); // 90 degrees
		
		// Sum of angles should be π
		REQUIRE_THAT(angleA + angleB + angleC, WithinAbs(Constants::PI, REAL(1e-10)));
	}

	TEST_CASE("Triangle2D - Point containment", "[Triangle2D][Contains]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(4, 0), Point2Cartesian(2, 3));
		
		// Inside
		REQUIRE(tri.Contains(Point2Cartesian(2, 1)));
		
		// Outside
		REQUIRE_FALSE(tri.Contains(Point2Cartesian(5, 5)));
		REQUIRE_FALSE(tri.Contains(Point2Cartesian(-1, 0)));
		
		// On vertex
		REQUIRE(tri.Contains(Point2Cartesian(0, 0)));
		
		// On edge
		REQUIRE(tri.Contains(Point2Cartesian(2, 0)));
	}

	TEST_CASE("Triangle2D - Classification with epsilon", "[Triangle2D][Classification]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// Right triangle (fixed floating-point comparison)
		Triangle2D right(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		REQUIRE(right.IsRight());
		REQUIRE_FALSE(right.IsEquilateral());

		// Isosceles triangle
		Triangle2D isosceles(Point2Cartesian(0, 0), Point2Cartesian(2, 0), Point2Cartesian(1, 2));
		REQUIRE(isosceles.IsIsosceles());
		REQUIRE_FALSE(isosceles.IsEquilateral());

		// Equilateral triangle (approximately)
		Real h = std::sqrt(REAL(REAL(3.0))) / REAL(REAL(2.0));
		Triangle2D equilateral(Point2Cartesian(0, 0), Point2Cartesian(1, 0), Point2Cartesian(REAL(REAL(0.5)), h));
		REQUIRE(equilateral.IsEquilateral());
		REQUIRE(equilateral.IsIsosceles());
	}

	TEST_CASE("Triangle2D - Orientation", "[Triangle2D][Orientation]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		// CCW triangle
		Triangle2D ccw(Point2Cartesian(0, 0), Point2Cartesian(2, 0), Point2Cartesian(1, 2));
		REQUIRE(ccw.IsCounterClockwise());
		REQUIRE_FALSE(ccw.IsClockwise());

		// CW triangle
		Triangle2D cw(Point2Cartesian(0, 0), Point2Cartesian(1, 2), Point2Cartesian(2, 0));
		REQUIRE(cw.IsClockwise());
		REQUIRE_FALSE(cw.IsCounterClockwise());
	}

	TEST_CASE("Triangle2D - Cached side lengths", "[Triangle2D][Performance]")
	{
				TEST_PRECISION_INFO();
			TEST_PRECISION_INFO();
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		
		// Multiple calls should use cached values
		Real a1 = tri.A();
		Real a2 = tri.A();
		Real a3 = tri.A();
		
		REQUIRE_THAT(a1, WithinAbs(a2, REAL(1e-15)));
		REQUIRE_THAT(a2, WithinAbs(a3, REAL(1e-15)));
		REQUIRE_THAT(a1, WithinAbs(REAL(REAL(3.0)), 1e-10)); // A() is distance from pnt1 to pnt2
		
		// Verify other sides
		REQUIRE_THAT(tri.B(), WithinAbs(REAL(REAL(5.0)), 1e-10)); // B() is pnt2 to pnt3 (hypotenuse)
		REQUIRE_THAT(tri.C(), WithinAbs(REAL(REAL(4.0)), 1e-10)); // C() is pnt3 to pnt1
	}
}
