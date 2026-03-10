#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/base/Geometry/Geometry2D.h"
#endif

using namespace MML;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;

namespace MML::Tests::Base::Geometry2DTests
{
	/*********************************************************************/
	/*****                       Line2D tests                        *****/
	/*********************************************************************/
	TEST_CASE("Line2D::Init", "[Line2D]") {
		TEST_PRECISION_INFO();
		Point2Cartesian    start_pnt(1, 1);
		Vector2Cartesian   direction(1, 0);

		Line2D     line(start_pnt, direction);

		REQUIRE(line.StartPoint().X() == REAL(1.0));
		REQUIRE(line.StartPoint().X() == REAL(1.0));
		REQUIRE(line.Direction().X() == REAL(1.0));
		REQUIRE(line.Direction().Y() == REAL(0.0));

		Point2Cartesian    pnt1(1, 1);
		Point2Cartesian    pnt2(2, 2);

		Line2D     line2(pnt1, pnt2);

		REQUIRE(line2.StartPoint().X() == REAL(1.0));
		REQUIRE(line2.StartPoint().Y() == REAL(1.0));
		REQUIRE_THAT(line2.Direction().X() , RealWithinRel(REAL(0.7071067811865475), 1e-5));
		REQUIRE_THAT(line2.Direction().Y() , RealWithinRel(REAL(0.7071067811865475), 1e-5));
	}

	/*********************************************************************/
	/*****                   SegmentLine2D tests                     *****/
	/*********************************************************************/
	TEST_CASE("SegmentLine2D::Init", "[SegmentLine2D]")
	{
		TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(1.0), REAL(2.0));
		Point2Cartesian p2(REAL(4.0), REAL(6.0));
		SegmentLine2D seg(p1, p2);
		REQUIRE_THAT(seg.StartPoint().X() , RealWithinRel(REAL(1.0), 1e-5));
		REQUIRE_THAT(seg.StartPoint().Y() , RealWithinRel(REAL(2.0), 1e-5));
		REQUIRE_THAT(seg.EndPoint().X() , RealWithinRel(REAL(4.0), 1e-5));
		REQUIRE_THAT(seg.EndPoint().Y() , RealWithinRel(REAL(6.0), 1e-5));
	}

	TEST_CASE("SegmentLine2D::PointOnSegment", "[SegmentLine2D]")
	{
		TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(0.0), REAL(0.0));
		Point2Cartesian p2(REAL(2.0), REAL(2.0));
		SegmentLine2D seg(p1, p2);

		// Start and end points should be on the segment
		REQUIRE(seg.PointOnSegment(REAL(0.0)) == p1);
		REQUIRE(seg.PointOnSegment(REAL(1.0)) == p2);

		// Midpoint
		Point2Cartesian mid = seg.PointOnSegment(REAL(0.5));
		REQUIRE_THAT(mid.X() , RealWithinRel(REAL(1.0), 1e-5));
		REQUIRE_THAT(mid.Y() , RealWithinRel(REAL(1.0), 1e-5));
	}

	TEST_CASE("SegmentLine2D::Length", "[SegmentLine2D]")
	{
		TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(1.0), REAL(1.0));
		Point2Cartesian p2(REAL(4.0), REAL(5.0));
		SegmentLine2D seg(p1, p2);
		Real expected_length = sqrt(pow(REAL(4.0) - REAL(1.0), 2) + pow(REAL(5.0) - REAL(1.0), 2));
		REQUIRE_THAT(seg.StartPoint().Dist(seg.EndPoint()) , RealWithinRel(expected_length, REAL(1e-5)));
	}

	TEST_CASE("SegmentLine2D::Direction", "[SegmentLine2D]")
	{
		TEST_PRECISION_INFO();
		Point2Cartesian p1(REAL(0.0), REAL(0.0));
		Point2Cartesian p2(REAL(3.0), REAL(4.0));
		SegmentLine2D seg(p1, p2);
		Vector2Cartesian dir(p1, p2);
		// Direction should be (3,4) normalized
		Real norm = sqrt(REAL(3.0) * REAL(3.0) + REAL(4.0) * REAL(4.0));
		REQUIRE_THAT(dir.X() / norm , RealWithinRel(REAL(0.6), 1e-5));
		REQUIRE_THAT(dir.Y() / norm , RealWithinRel(REAL(0.8), 1e-5));
	}

	/*********************************************************************/
	/*****                     Triangle2D tests                      *****/
	/*********************************************************************/
	TEST_CASE("Triangle2D::Init", "[Triangle2D]")
	{
		TEST_PRECISION_INFO();
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		REQUIRE_THAT(tri.Pnt1().X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(tri.Pnt2().X(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THAT(tri.Pnt3().Y(), WithinAbs(REAL(4.0), 1e-10));
	}

	TEST_CASE("Triangle2D::Area", "[Triangle2D]")
	{
		TEST_PRECISION_INFO();
		Point2Cartesian p1(0, 0);
		Point2Cartesian p2(4, 0);
		Point2Cartesian p3(0, 3);
		Triangle2D tri(p1, p2, p3);

		Real area = tri.Area();
		REQUIRE_THAT(area , RealWithinRel(REAL(6.0), 1e-5));

		// Degenerate triangle (collinear)
		Point2Cartesian p4(1, 1);
		Point2Cartesian p5(2, 2);
		Point2Cartesian p6(3, 3);
		Triangle2D degenerate(p4, p5, p6);
		Real deg_area = degenerate.Area();
		REQUIRE_THAT(deg_area , RealWithinRel(REAL(0.0), 1e-5));
	}

	/*********************************************************************/
	/*****                     Polygon2D tests                       *****/
	/*********************************************************************/
	TEST_CASE("Polygon2D::Init", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		Polygon2D poly({Point2Cartesian(0, 0), Point2Cartesian(1, 0), Point2Cartesian(0, 1)});
		REQUIRE(poly.NumVertices() == 3);
	}

	TEST_CASE("Polygon2D::Area1", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		// Square with side 2, area should be 4
		std::vector<Point2Cartesian> square = {
			{0, 0}, {2, 0}, {2, 2}, {0, 2}
		};
		Polygon2D poly(square);
		REQUIRE_THAT(poly.Area() , RealWithinRel(REAL(4.0), 1e-5));
	}

	TEST_CASE("Polygon2D::Area2", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		// Triangle as polygon, area should be 6
		std::vector<Point2Cartesian> triangle = {
			{0, 0}, {4, 0}, {0, 3}
		};
		Polygon2D poly(triangle);
		REQUIRE_THAT(poly.Area() , RealWithinRel(REAL(6.0), 1e-5));
	}

	TEST_CASE("Polygon2D::Area3", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		// Degenerate polygon (collinear points), area should be 0
		std::vector<Point2Cartesian> line = {
			{1, 1}, {2, 2}, {3, 3}
		};
		Polygon2D poly(line);
		REQUIRE_THAT(poly.Area() , RealWithinRel(REAL(0.0), 1e-5));

		// Regular pentagon, side length 1, area ≈ 2.3776
		std::vector<Point2Cartesian> pentagon = {
			{REAL(0.0), REAL(1.0)},
			{REAL(0.9510565163), REAL(0.309017)},
			{REAL(0.5877852523), -REAL(0.809017)},
			{-REAL(0.5877852523), -REAL(0.809017)},
			{-REAL(0.9510565163), REAL(0.309017)}
		};
		Polygon2D poly2(pentagon);
		REQUIRE_THAT(poly2.Area() , RealWithinRel(REAL(2.3776), 1e-4));
	}

	TEST_CASE("Polygon2D::IsSimple", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		// Simple square
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		REQUIRE(square.IsSimple());
	}

	TEST_CASE("Polygon2D::IsConvex", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		// Convex square
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		REQUIRE(square.IsConvex());
	}
	TEST_CASE("Polygon2D::IsInside", "[Polygon2D]")
	{
		TEST_PRECISION_INFO();
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		REQUIRE(square.Contains(Point2Cartesian(1, 1)));
		REQUIRE_FALSE(square.Contains(Point2Cartesian(5, 5)));
	}

	/*********************************************************************/
	/*****          Line2D Enhanced Functionality Tests              *****/
	/*********************************************************************/
	TEST_CASE("Line2D - Distance to point", "[Line2D][Distance]")
	{
		TEST_PRECISION_INFO();
		// Horizontal line at y=2
		Point2Cartesian p1(0, 2);
		Vector2Cartesian dir(1, 0);
		Line2D line(p1, dir);

		// Point above line
		Point2Cartesian p(5, 5);
		REQUIRE_THAT(line.DistanceToPoint(p), WithinAbs(REAL(3.0), 1e-10));

		// Point on line
		Point2Cartesian pOn(10, 2);
		REQUIRE_THAT(line.DistanceToPoint(pOn), WithinAbs(REAL(0.0), 1e-10));
	}

	TEST_CASE("Line2D - Closest point", "[Line2D][ClosestPoint]")
	{
		TEST_PRECISION_INFO();
		// Line through origin at 45 degrees
		Point2Cartesian origin(0, 0);
		Vector2Cartesian dir(1, 1);
		Line2D line(origin, dir);

		Point2Cartesian p(2, 0);
		Point2Cartesian closest = line.ClosestPoint(p);
		
		REQUIRE_THAT(closest.X(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(closest.Y(), WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Line2D - Perpendicular line", "[Line2D][Perpendicular]")
	{
		TEST_PRECISION_INFO();
		// Horizontal line
		Line2D horizontal(Point2Cartesian(0, 0), Vector2Cartesian(1, 0));
		
		Point2Cartesian through(5, 3);
		Line2D perp = horizontal.Perpendicular(through);
		
		// Perpendicular should be vertical (direction 0,1 or 0,-1)
		REQUIRE_THAT(std::abs(perp.Direction().X()), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(std::abs(perp.Direction().Y()), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE(perp.Contains(through));
	}

	TEST_CASE("Line2D - Line intersection", "[Line2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		// Two intersecting lines
		Line2D line1(Point2Cartesian(0, 0), Vector2Cartesian(1, 0)); // Horizontal
		Line2D line2(Point2Cartesian(0, 0), Vector2Cartesian(0, 1)); // Vertical

		Point2Cartesian intersection;
		REQUIRE(line1.Intersects(line2, &intersection));
		REQUIRE_THAT(intersection.X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(intersection.Y(), WithinAbs(REAL(0.0), 1e-10));

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
		SegmentLine2D seg(Point2Cartesian(2, 3), Point2Cartesian(8, 7));
		Point2Cartesian mid = seg.Midpoint();
		
		REQUIRE_THAT(mid.X(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(mid.Y(), WithinAbs(REAL(5.0), 1e-10));
	}

	TEST_CASE("SegmentLine2D - Distance to point", "[SegmentLine2D][Distance]")
	{
		TEST_PRECISION_INFO();
		SegmentLine2D seg(Point2Cartesian(0, 0), Point2Cartesian(4, 0));

		// Point above segment midpoint
		Point2Cartesian p1(2, 3);
		REQUIRE_THAT(seg.DistanceToPoint(p1), WithinAbs(REAL(3.0), 1e-10));

		// Point beyond segment end
		Point2Cartesian p2(6, 0);
		REQUIRE_THAT(seg.DistanceToPoint(p2), WithinAbs(REAL(2.0), 1e-10));

		// Point before segment start
		Point2Cartesian p3(-2, 0);
		REQUIRE_THAT(seg.DistanceToPoint(p3), WithinAbs(REAL(2.0), 1e-10));
	}

	TEST_CASE("SegmentLine2D - Closest point", "[SegmentLine2D][ClosestPoint]")
	{
		TEST_PRECISION_INFO();
		SegmentLine2D seg(Point2Cartesian(0, 0), Point2Cartesian(10, 0));

		// Point above segment
		Point2Cartesian p1(5, 5);
		Point2Cartesian closest1 = seg.ClosestPoint(p1);
		REQUIRE_THAT(closest1.X(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(closest1.Y(), WithinAbs(REAL(0.0), 1e-10));

		// Point beyond end
		Point2Cartesian p2(15, 3);
		Point2Cartesian closest2 = seg.ClosestPoint(p2);
		REQUIRE_THAT(closest2.X(), WithinAbs(REAL(10.0), 1e-10));
		REQUIRE_THAT(closest2.Y(), WithinAbs(REAL(0.0), 1e-10));
	}

	TEST_CASE("SegmentLine2D - Segment intersection", "[SegmentLine2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		// Intersecting segments
		SegmentLine2D seg1(Point2Cartesian(0, 0), Point2Cartesian(4, 4));
		SegmentLine2D seg2(Point2Cartesian(0, 4), Point2Cartesian(4, 0));

		Point2Cartesian intersection;
		REQUIRE(seg1.Intersects(seg2, &intersection));
		REQUIRE_THAT(intersection.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(intersection.Y(), WithinAbs(REAL(2.0), 1e-10));

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
		// 3-4-5 right triangle
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		REQUIRE_THAT(tri.Perimeter(), WithinAbs(REAL(12.0), 1e-10));
	}

	TEST_CASE("Triangle2D - Centroid", "[Triangle2D][Centroid]")
	{
		TEST_PRECISION_INFO();
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(6, 0), Point2Cartesian(0, 6));
		Point2Cartesian centroid = tri.Centroid();
		
		REQUIRE_THAT(centroid.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(centroid.Y(), WithinAbs(REAL(2.0), 1e-10));
	}

	TEST_CASE("Triangle2D - Circumcenter and CircumRadius", "[Triangle2D][Circumcenter]")
	{
		TEST_PRECISION_INFO();
		// Right triangle with hypotenuse as diameter of circumcircle
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(4, 0), Point2Cartesian(0, 3));
		
		Point2Cartesian circumcenter = tri.Circumcenter();
		REQUIRE_THAT(circumcenter.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(circumcenter.Y(), WithinAbs(REAL(1.5), 1e-10));

		Real radius = tri.CircumRadius();
		REQUIRE_THAT(radius, WithinAbs(REAL(2.5), 1e-10)); // 5/2 (hypotenuse/2)
	}

	TEST_CASE("Triangle2D - Incenter and InRadius", "[Triangle2D][Incenter]")
	{
		TEST_PRECISION_INFO();
		// 3-4-5 right triangle
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		
		Point2Cartesian incenter = tri.Incenter();
		Real inradius = tri.InRadius();
		
		REQUIRE_THAT(inradius, WithinAbs(REAL(1.0), 1e-10)); // Area / s = 6 / 6 = 1
		
		// Incenter should be at distance inradius from all sides
		REQUIRE_THAT(incenter.X(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(incenter.Y(), WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Triangle2D - Angles", "[Triangle2D][Angles]")
	{
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
		REQUIRE_THAT(angleB, WithinAbs(Constants::PI / REAL(2.0), 1e-8)); // 90 degrees
		
		// Sum of angles should be π
		REQUIRE_THAT(angleA + angleB + angleC, WithinAbs(Constants::PI, REAL(1e-10)));
	}

	TEST_CASE("Triangle2D - Point containment", "[Triangle2D][Contains]")
	{
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
		// Right triangle (fixed floating-point comparison)
		Triangle2D right(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		REQUIRE(right.IsRight());
		REQUIRE_FALSE(right.IsEquilateral());

		// Isosceles triangle
		Triangle2D isosceles(Point2Cartesian(0, 0), Point2Cartesian(2, 0), Point2Cartesian(1, 2));
		REQUIRE(isosceles.IsIsosceles());
		REQUIRE_FALSE(isosceles.IsEquilateral());

		// Equilateral triangle (approximately)
		Real h = std::sqrt(REAL(3.0)) / REAL(2.0);
		Triangle2D equilateral(Point2Cartesian(0, 0), Point2Cartesian(1, 0), Point2Cartesian(REAL(0.5), h));
		REQUIRE(equilateral.IsEquilateral());
		REQUIRE(equilateral.IsIsosceles());
	}

	TEST_CASE("Triangle2D - Orientation", "[Triangle2D][Orientation]")
	{
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
		Triangle2D tri(Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4));
		
		// Multiple calls should use cached values
		Real a1 = tri.A();
		Real a2 = tri.A();
		Real a3 = tri.A();
		
		REQUIRE_THAT(a1, WithinAbs(a2, REAL(1e-15)));
		REQUIRE_THAT(a2, WithinAbs(a3, REAL(1e-15)));
		REQUIRE_THAT(a1, WithinAbs(REAL(3.0), 1e-10)); // A() is distance from pnt1 to pnt2
		
		// Verify other sides
		REQUIRE_THAT(tri.B(), WithinAbs(REAL(5.0), 1e-10)); // B() is pnt2 to pnt3 (hypotenuse)
		REQUIRE_THAT(tri.C(), WithinAbs(REAL(4.0), 1e-10)); // C() is pnt3 to pnt1
	}

	TEST_CASE("Triangle2D - Signed area and vertex accessors", "[Triangle2D][SignedArea]")
	{
		TEST_PRECISION_INFO();
		// CCW triangle - positive signed area
		Triangle2D ccw(Point2Cartesian(0, 0), Point2Cartesian(4, 0), Point2Cartesian(2, 3));
		Real signedAreaCCW = ccw.SignedArea();
		REQUIRE(signedAreaCCW > 0);
		REQUIRE_THAT(std::abs(signedAreaCCW), WithinAbs(REAL(6.0), 1e-10)); // Area = 0.5 * base * height = 0.5 * 4 * 3 = 6

		// CW triangle - negative signed area
		Triangle2D cw(Point2Cartesian(0, 0), Point2Cartesian(2, 3), Point2Cartesian(4, 0));
		Real signedAreaCW = cw.SignedArea();
		REQUIRE(signedAreaCW < 0);
		REQUIRE_THAT(std::abs(signedAreaCW), WithinAbs(REAL(6.0), 1e-10));

		// Vertex accessors
		REQUIRE_THAT(ccw.Pnt1().X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(ccw.Pnt2().X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(ccw.Pnt3().Y(), WithinAbs(REAL(3.0), 1e-10));
		
		// Mutable accessor invalidates cache
		ccw.Pnt1() = Point2Cartesian(1, 0);
		REQUIRE_THAT(ccw.A(), WithinAbs(REAL(3.0), 1e-10)); // Distance (1,0) to (4,0) = 3
	}

	/*********************************************************************/
	/*****        Polygon2D Extended Functionality Tests             *****/
	/*********************************************************************/
	TEST_CASE("Polygon2D - Perimeter calculation", "[Polygon2D][Perimeter]")
	{
		TEST_PRECISION_INFO();
		// Square 2x2
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		REQUIRE_THAT(square.Perimeter(), WithinAbs(REAL(8.0), 1e-10));
		
		// Regular triangle
		Polygon2D triangle({
			Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(0, 4)
		});
		REQUIRE_THAT(triangle.Perimeter(), WithinAbs(REAL(12.0), 1e-10)); // 3 + 4 + 5 = 12
	}

	TEST_CASE("Polygon2D - Centroid calculation", "[Polygon2D][Centroid]")
	{
		TEST_PRECISION_INFO();
		// Square centered at (1,1)
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		Point2Cartesian centroid = square.Centroid();
		REQUIRE_THAT(centroid.X(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(centroid.Y(), WithinAbs(REAL(1.0), 1e-10));
		
		// Single point
		Polygon2D singlePoint({Point2Cartesian(5, 7)});
		Point2Cartesian c1 = singlePoint.Centroid();
		REQUIRE_THAT(c1.X(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(c1.Y(), WithinAbs(REAL(7.0), 1e-10));
		
		// Two points (line segment)
		Polygon2D twoPoints({Point2Cartesian(0, 0), Point2Cartesian(4, 6)});
		Point2Cartesian c2 = twoPoints.Centroid();
		REQUIRE_THAT(c2.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(c2.Y(), WithinAbs(REAL(3.0), 1e-10));
	}

	TEST_CASE("Polygon2D - Bounding box", "[Polygon2D][BoundingBox]")
	{
		TEST_PRECISION_INFO();
		Polygon2D poly({
			Point2Cartesian(1, 2), Point2Cartesian(5, 3),
			Point2Cartesian(3, 7), Point2Cartesian(0, 4)
		});
		auto bbox = poly.GetBoundingBox();
		REQUIRE_THAT(bbox.minX, WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(bbox.maxX, WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(bbox.minY, WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(bbox.maxY, WithinAbs(REAL(7.0), 1e-10));
		REQUIRE_THAT(bbox.Width(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(bbox.Height(), WithinAbs(REAL(5.0), 1e-10));
		
		Point2Cartesian center = bbox.Center();
		REQUIRE_THAT(center.X(), WithinAbs(REAL(2.5), 1e-10));
		REQUIRE_THAT(center.Y(), WithinAbs(REAL(4.5), 1e-10));
	}

	TEST_CASE("Polygon2D - Winding number", "[Polygon2D][WindingNumber]")
	{
		TEST_PRECISION_INFO();
		// CCW square
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(4, 0),
			Point2Cartesian(4, 4), Point2Cartesian(0, 4)
		});
		
		// Point inside - should have winding number 1 (or -1 for CW)
		int wn_inside = square.WindingNumber(Point2Cartesian(2, 2));
		REQUIRE(wn_inside != 0);
		
		// Point outside
		int wn_outside = square.WindingNumber(Point2Cartesian(10, 10));
		REQUIRE(wn_outside == 0);
	}

	TEST_CASE("Polygon2D - Edge and vertex access", "[Polygon2D][Access]")
	{
		TEST_PRECISION_INFO();
		Polygon2D poly({
			Point2Cartesian(0, 0), Point2Cartesian(3, 0),
			Point2Cartesian(3, 4), Point2Cartesian(0, 4)
		});
		
		// NumVertices
		REQUIRE(poly.NumVertices() == 4);
		
		// Vertex access with bounds checking
		REQUIRE_THAT(poly.Vertex(0).X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(poly.Vertex(2).X(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THROWS(poly.Vertex(5)); // Out of range
		REQUIRE_THROWS(poly.Vertex(-1));
		
		// Edge access
		SegmentLine2D edge0 = poly.Edge(0); // (0,0) to (3,0)
		REQUIRE_THAT(edge0.Length(), WithinAbs(REAL(3.0), 1e-10));
		
		SegmentLine2D edge1 = poly.Edge(1); // (3,0) to (3,4)
		REQUIRE_THAT(edge1.Length(), WithinAbs(REAL(4.0), 1e-10));
		
		REQUIRE_THROWS(poly.Edge(10)); // Out of range
		
		// Index operator
		REQUIRE_THAT(poly[0].X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(poly[3].Y(), WithinAbs(REAL(4.0), 1e-10));
	}

	TEST_CASE("Polygon2D - Vertex manipulation", "[Polygon2D][Manipulation]")
	{
		TEST_PRECISION_INFO();
		Polygon2D poly;
		REQUIRE(poly.NumVertices() == 0);
		
		// AddVertex
		poly.AddVertex(Point2Cartesian(0, 0));
		poly.AddVertex(Point2Cartesian(1, 0));
		poly.AddVertex(Point2Cartesian(1, 1));
		REQUIRE(poly.NumVertices() == 3);
		
		// Clear
		poly.Clear();
		REQUIRE(poly.NumVertices() == 0);
		
		// Rebuild and reverse
		Polygon2D ccw({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0), Point2Cartesian(1, 2)
		});
		REQUIRE(ccw.IsCounterClockwise());
		
		ccw.Reverse();
		REQUIRE(ccw.IsClockwise());
	}

	TEST_CASE("Polygon2D - Triangle decomposition", "[Polygon2D][Triangularization]")
	{
		TEST_PRECISION_INFO();
		// Square
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		auto triangles = square.Triangularization();
		REQUIRE(triangles.size() == 2);
		
		// Total area should match
		Real totalArea = 0;
		for (const auto& tri : triangles) {
			totalArea += tri.Area();
		}
		REQUIRE_THAT(totalArea, WithinAbs(REAL(4.0), 1e-10)); // 2x2 square
		
		// Pentagon
		Polygon2D pentagon({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0), Point2Cartesian(3, 1),
			Point2Cartesian(1.5, 3), Point2Cartesian(-0.5, 1)
		});
		auto pentTriangles = pentagon.Triangularization();
		REQUIRE(pentTriangles.size() == 3); // n-2 triangles for convex polygon
	}

	TEST_CASE("Polygon2D - Orientation and simplicity", "[Polygon2D][Properties]")
	{
		TEST_PRECISION_INFO();
		// Simple convex polygon (CCW square)
		Polygon2D square({
			Point2Cartesian(0, 0), Point2Cartesian(4, 0),
			Point2Cartesian(4, 4), Point2Cartesian(0, 4)
		});
		REQUIRE(square.IsSimple());
		REQUIRE(square.IsConvex());
		REQUIRE(square.IsCounterClockwise());
		REQUIRE_FALSE(square.IsClockwise());
		
		// Self-intersecting polygon (bowtie)
		Polygon2D bowtie({
			Point2Cartesian(0, 0), Point2Cartesian(2, 2),
			Point2Cartesian(2, 0), Point2Cartesian(0, 2)
		});
		REQUIRE_FALSE(bowtie.IsSimple());
		
		// Concave polygon (L-shape)
		Polygon2D lShape({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0), Point2Cartesian(2, 1),
			Point2Cartesian(1, 1), Point2Cartesian(1, 2), Point2Cartesian(0, 2)
		});
		REQUIRE(lShape.IsSimple());
		REQUIRE_FALSE(lShape.IsConvex());
	}

	TEST_CASE("Polygon2D - Contains with ray casting", "[Polygon2D][Contains]")
	{
		TEST_PRECISION_INFO();
		// L-shaped concave polygon
		Polygon2D lShape({
			Point2Cartesian(0, 0), Point2Cartesian(3, 0), Point2Cartesian(3, 1),
			Point2Cartesian(1, 1), Point2Cartesian(1, 3), Point2Cartesian(0, 3)
		});
		
		// Inside L
		REQUIRE(lShape.Contains(Point2Cartesian(0.5, 0.5)));
		REQUIRE(lShape.Contains(Point2Cartesian(0.5, 2)));
		
		// Outside L (in the cut-out corner)
		REQUIRE_FALSE(lShape.Contains(Point2Cartesian(2, 2)));
		
		// Completely outside
		REQUIRE_FALSE(lShape.Contains(Point2Cartesian(-1, -1)));
		REQUIRE_FALSE(lShape.Contains(Point2Cartesian(5, 5)));
		
		// Legacy alias
		REQUIRE(lShape.IsInside(Point2Cartesian(0.5, 0.5)));
	}

	TEST_CASE("Polygon2D - Signed area orientation", "[Polygon2D][SignedArea]")
	{
		TEST_PRECISION_INFO();
		// CCW square - positive area
		Polygon2D ccw({
			Point2Cartesian(0, 0), Point2Cartesian(2, 0),
			Point2Cartesian(2, 2), Point2Cartesian(0, 2)
		});
		REQUIRE(ccw.SignedArea() > 0);
		REQUIRE_THAT(ccw.Area(), WithinAbs(REAL(4.0), 1e-10));
		
		// CW square - negative signed area
		Polygon2D cw({
			Point2Cartesian(0, 0), Point2Cartesian(0, 2),
			Point2Cartesian(2, 2), Point2Cartesian(2, 0)
		});
		REQUIRE(cw.SignedArea() < 0);
		REQUIRE_THAT(cw.Area(), WithinAbs(REAL(4.0), 1e-10)); // Unsigned area is positive
	}

	/*********************************************************************/
	/*****                 Circle2D Tests                            *****/
	/*********************************************************************/
	TEST_CASE("Circle2D - Default constructor", "[Circle2D][Constructor]")
	{
		TEST_PRECISION_INFO();
		Circle2D c; // Unit circle at origin
		REQUIRE_THAT(c.Center().X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(c.Center().Y(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(c.Radius(), WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Circle2D - Point and radius constructor", "[Circle2D][Constructor]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(Point2Cartesian(3, 4), 5);
		REQUIRE_THAT(c.Center().X(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THAT(c.Center().Y(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(c.Radius(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(c.Diameter(), WithinAbs(REAL(10.0), 1e-10));
		
		// Negative radius throws
		REQUIRE_THROWS(Circle2D(Point2Cartesian(0, 0), -1));
	}

	TEST_CASE("Circle2D - Coordinate constructor", "[Circle2D][Constructor]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(2, 3, 4);
		REQUIRE_THAT(c.Center().X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(c.Center().Y(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THAT(c.Radius(), WithinAbs(REAL(4.0), 1e-10));
		
		REQUIRE_THROWS(Circle2D(0, 0, -1)); // Negative radius
	}

	TEST_CASE("Circle2D - Factory methods", "[Circle2D][Factory]")
	{
		TEST_PRECISION_INFO();
		// FromDiameter
		Circle2D fromDiam = Circle2D::FromDiameter(Point2Cartesian(0, 0), Point2Cartesian(4, 0));
		REQUIRE_THAT(fromDiam.Center().X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(fromDiam.Center().Y(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(fromDiam.Radius(), WithinAbs(REAL(2.0), 1e-10));
		
		// FromThreePoints - circumcircle of right triangle at origin
		Circle2D fromThree = Circle2D::FromThreePoints(
			Point2Cartesian(0, 0), Point2Cartesian(4, 0), Point2Cartesian(0, 3)
		);
		REQUIRE_THAT(fromThree.Radius(), WithinAbs(REAL(2.5), 1e-10)); // hypotenuse/2 = 5/2
		
		// UnitCircle
		Circle2D unit = Circle2D::UnitCircle();
		REQUIRE_THAT(unit.Center().X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(unit.Center().Y(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(unit.Radius(), WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Circle2D - Area and circumference", "[Circle2D][Measurements]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 3);
		REQUIRE_THAT(c.Area(), WithinAbs(REAL(9.0) * Constants::PI, 1e-10));
		REQUIRE_THAT(c.Circumference(), WithinAbs(REAL(6.0) * Constants::PI, 1e-10));
		REQUIRE_THAT(c.Perimeter(), WithinAbs(REAL(6.0) * Constants::PI, 1e-10));
	}

	TEST_CASE("Circle2D - Point on circle at angle", "[Circle2D][PointAt]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 2);
		
		// At 0 radians (right)
		Point2Cartesian p0 = c.PointAt(0);
		REQUIRE_THAT(p0.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(p0.Y(), WithinAbs(REAL(0.0), 1e-10));
		
		// At π/2 radians (top)
		Point2Cartesian p90 = c.PointAt(Constants::PI / 2);
		REQUIRE_THAT(p90.X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(p90.Y(), WithinAbs(REAL(2.0), 1e-10));
		
		// At π radians (left)
		Point2Cartesian p180 = c.PointAt(Constants::PI);
		REQUIRE_THAT(p180.X(), WithinAbs(REAL(-2.0), 1e-10));
		REQUIRE_THAT(p180.Y(), WithinAbs(REAL(0.0), 1e-10));
	}

	TEST_CASE("Circle2D - Arc, sector, and chord", "[Circle2D][Geometry]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 4);
		
		// Arc length for 90 degrees
		Real arcLen = c.ArcLength(Constants::PI / 2);
		REQUIRE_THAT(arcLen, WithinAbs(REAL(2.0) * Constants::PI, 1e-10)); // r * theta = 4 * π/2 = 2π
		
		// Sector area for 90 degrees
		Real sectorArea = c.SectorArea(Constants::PI / 2);
		REQUIRE_THAT(sectorArea, WithinAbs(REAL(4.0) * Constants::PI, 1e-10)); // 0.5 * r² * theta = 0.5 * 16 * π/2 = 4π
		
		// Chord length for 90 degrees
		Real chordLen = c.ChordLength(Constants::PI / 2);
		REQUIRE_THAT(chordLen, WithinAbs(REAL(4.0) * std::sqrt(REAL(2.0)), 1e-10)); // 2r*sin(theta/2) = 8*sin(π/4) = 8*√2/2 = 4√2
	}

	TEST_CASE("Circle2D - Distance calculations", "[Circle2D][Distance]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 5);
		
		// Distance from center
		REQUIRE_THAT(c.DistanceFromCenter(Point2Cartesian(3, 4)), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(c.DistanceFromCenter(Point2Cartesian(0, 0)), WithinAbs(REAL(0.0), 1e-10));
		
		// Signed distance (negative inside, positive outside)
		REQUIRE(c.SignedDistance(Point2Cartesian(0, 0)) < 0); // Inside
		REQUIRE_THAT(c.SignedDistance(Point2Cartesian(0, 0)), WithinAbs(REAL(-5.0), 1e-10));
		REQUIRE(c.SignedDistance(Point2Cartesian(10, 0)) > 0); // Outside
		REQUIRE_THAT(c.SignedDistance(Point2Cartesian(10, 0)), WithinAbs(REAL(5.0), 1e-10));
		
		// Distance to boundary
		REQUIRE_THAT(c.DistanceToPoint(Point2Cartesian(0, 0)), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(c.DistanceToPoint(Point2Cartesian(10, 0)), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(c.DistanceToPoint(Point2Cartesian(5, 0)), WithinAbs(REAL(0.0), 1e-10)); // On boundary
	}

	TEST_CASE("Circle2D - Point containment", "[Circle2D][Contains]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 3);
		
		// Inside
		REQUIRE(c.Contains(Point2Cartesian(0, 0)));
		REQUIRE(c.Contains(Point2Cartesian(1, 1)));
		
		// On boundary
		REQUIRE(c.Contains(Point2Cartesian(3, 0)));
		REQUIRE(c.IsOnBoundary(Point2Cartesian(3, 0)));
		REQUIRE(c.IsOnBoundary(Point2Cartesian(0, 3)));
		
		// Outside
		REQUIRE_FALSE(c.Contains(Point2Cartesian(4, 0)));
		REQUIRE_FALSE(c.IsOnBoundary(Point2Cartesian(0, 0)));
	}

	TEST_CASE("Circle2D - Closest point", "[Circle2D][ClosestPoint]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 4);
		
		// Point outside
		Point2Cartesian closest = c.ClosestPoint(Point2Cartesian(8, 0));
		REQUIRE_THAT(closest.X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(closest.Y(), WithinAbs(REAL(0.0), 1e-10));
		
		// Point inside
		closest = c.ClosestPoint(Point2Cartesian(1, 0));
		REQUIRE_THAT(closest.X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(closest.Y(), WithinAbs(REAL(0.0), 1e-10));
		
		// Point at center - returns arbitrary point (at angle 0)
		closest = c.ClosestPoint(Point2Cartesian(0, 0));
		REQUIRE_THAT(closest.Dist(c.Center()), WithinAbs(REAL(4.0), 1e-10));
	}

	TEST_CASE("Circle2D - Tangent vector", "[Circle2D][Tangent]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 5);
		
		// Tangent at 0 radians (should be perpendicular to radius, pointing "up")
		Vector2Cartesian t0 = c.TangentAt(0);
		REQUIRE_THAT(t0.X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(t0.Y(), WithinAbs(REAL(1.0), 1e-10));
		
		// Tangent at π/2 radians
		Vector2Cartesian t90 = c.TangentAt(Constants::PI / 2);
		REQUIRE_THAT(t90.X(), WithinAbs(REAL(-1.0), 1e-10));
		REQUIRE_THAT(t90.Y(), WithinAbs(REAL(0.0), 1e-10));
	}

	TEST_CASE("Circle2D - Circle-circle intersection", "[Circle2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		Circle2D c1(0, 0, 3);
		Circle2D c2(4, 0, 3);
		
		// Overlapping circles
		REQUIRE(c1.Intersects(c2));
		
		// Touching circles (externally tangent)
		Circle2D c3(6, 0, 3);
		REQUIRE(c1.Intersects(c3)); // Just touching
		
		// Separated circles
		Circle2D c4(10, 0, 2);
		REQUIRE_FALSE(c1.Intersects(c4));
		
		// One inside the other (not intersecting in geometric sense)
		Circle2D c5(0, 0, 1);
		REQUIRE_FALSE(c1.Intersects(c5)); // c5 is fully inside c1
	}

	TEST_CASE("Circle2D - Circle-line intersection", "[Circle2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 5);
		
		// Line through center
		Line2D lineThrough(Point2Cartesian(0, 0), Vector2Cartesian(1, 0));
		REQUIRE(c.Intersects(lineThrough));
		
		// Tangent line
		Line2D tangent(Point2Cartesian(5, 0), Vector2Cartesian(0, 1));
		REQUIRE(c.Intersects(tangent));
		
		// Line outside circle
		Line2D outside(Point2Cartesian(10, 0), Vector2Cartesian(0, 1));
		REQUIRE_FALSE(c.Intersects(outside));
	}

	TEST_CASE("Circle2D - Bounding box", "[Circle2D][BoundingBox]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(3, 4, 2);
		auto bbox = c.GetBoundingBox();
		
		REQUIRE_THAT(bbox.minX, WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(bbox.maxX, WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(bbox.minY, WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(bbox.maxY, WithinAbs(REAL(6.0), 1e-10));
	}

	TEST_CASE("Circle2D - Translation and scaling", "[Circle2D][Transform]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 3);
		
		// Translate
		Circle2D translated = c.Translated(Vector2Cartesian(5, 7));
		REQUIRE_THAT(translated.Center().X(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(translated.Center().Y(), WithinAbs(REAL(7.0), 1e-10));
		REQUIRE_THAT(translated.Radius(), WithinAbs(REAL(3.0), 1e-10));
		
		// Scale
		Circle2D scaled = c.Scaled(2);
		REQUIRE_THAT(scaled.Radius(), WithinAbs(REAL(6.0), 1e-10));
		
		// Scale with non-positive factor throws
		REQUIRE_THROWS(c.Scaled(0));
		REQUIRE_THROWS(c.Scaled(-1));
	}

	TEST_CASE("Circle2D - Mutable accessors", "[Circle2D][Accessors]")
	{
		TEST_PRECISION_INFO();
		Circle2D c(0, 0, 1);
		
		c.Center() = Point2Cartesian(2, 3);
		REQUIRE_THAT(c.Center().X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(c.Center().Y(), WithinAbs(REAL(3.0), 1e-10));
		
		c.Radius() = 5;
		REQUIRE_THAT(c.Radius(), WithinAbs(REAL(5.0), 1e-10));
	}

	/*********************************************************************/
	/*****                   Box2D Tests                             *****/
	/*********************************************************************/
	TEST_CASE("Box2D - Default constructor", "[Box2D][Constructor]")
	{
		TEST_PRECISION_INFO();
		Box2D box; // Unit box from (0,0) to (1,1)
		REQUIRE_THAT(box.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box.MinY(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box.MaxX(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(box.MaxY(), WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Box2D - Point constructor with normalization", "[Box2D][Constructor]")
	{
		TEST_PRECISION_INFO();
		// Normal order
		Box2D box1(Point2Cartesian(0, 0), Point2Cartesian(4, 3));
		REQUIRE_THAT(box1.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box1.MaxX(), WithinAbs(REAL(4.0), 1e-10));
		
		// Reversed order - should normalize
		Box2D box2(Point2Cartesian(4, 3), Point2Cartesian(0, 0));
		REQUIRE_THAT(box2.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box2.MaxX(), WithinAbs(REAL(4.0), 1e-10));
	}

	TEST_CASE("Box2D - Coordinate constructor", "[Box2D][Constructor]")
	{
		TEST_PRECISION_INFO();
		Box2D box(1, 2, 5, 6);
		REQUIRE_THAT(box.MinX(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(box.MinY(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(box.MaxX(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(box.MaxY(), WithinAbs(REAL(6.0), 1e-10));
	}

	TEST_CASE("Box2D - Factory methods", "[Box2D][Factory]")
	{
		TEST_PRECISION_INFO();
		// FromCenterAndSize
		Box2D box1 = Box2D::FromCenterAndSize(Point2Cartesian(5, 5), 4, 2);
		REQUIRE_THAT(box1.MinX(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THAT(box1.MaxX(), WithinAbs(REAL(7.0), 1e-10));
		REQUIRE_THAT(box1.MinY(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(box1.MaxY(), WithinAbs(REAL(6.0), 1e-10));
		
		// FromCenterAndHalfExtents
		Box2D box2 = Box2D::FromCenterAndHalfExtents(Point2Cartesian(0, 0), 2, 3);
		REQUIRE_THAT(box2.MinX(), WithinAbs(REAL(-2.0), 1e-10));
		REQUIRE_THAT(box2.MaxX(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(box2.MinY(), WithinAbs(REAL(-3.0), 1e-10));
		REQUIRE_THAT(box2.MaxY(), WithinAbs(REAL(3.0), 1e-10));
		
		// FromPoints
		std::vector<Point2Cartesian> points = {
			Point2Cartesian(1, 2), Point2Cartesian(5, 3),
			Point2Cartesian(3, 7), Point2Cartesian(0, 4)
		};
		Box2D box3 = Box2D::FromPoints(points);
		REQUIRE_THAT(box3.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box3.MaxX(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(box3.MinY(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(box3.MaxY(), WithinAbs(REAL(7.0), 1e-10));
		
		// FromPoints with empty vector throws
		REQUIRE_THROWS(Box2D::FromPoints({}));
		
		// UnitBox
		Box2D unit = Box2D::UnitBox();
		REQUIRE_THAT(unit.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(unit.MaxX(), WithinAbs(REAL(1.0), 1e-10));
	}

	TEST_CASE("Box2D - Dimensions and measurements", "[Box2D][Measurements]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 3);
		
		REQUIRE_THAT(box.Width(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(box.Height(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THAT(box.Diagonal(), WithinAbs(REAL(5.0), 1e-10)); // 3-4-5 triangle
		REQUIRE_THAT(box.AspectRatio(), WithinAbs(REAL(4.0) / REAL(3.0), 1e-10));
		REQUIRE_THAT(box.Area(), WithinAbs(REAL(12.0), 1e-10));
		REQUIRE_THAT(box.Perimeter(), WithinAbs(REAL(14.0), 1e-10));
		
		Vector2Cartesian size = box.Size();
		REQUIRE_THAT(size.X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(size.Y(), WithinAbs(REAL(3.0), 1e-10));
		
		Vector2Cartesian halfExt = box.HalfExtents();
		REQUIRE_THAT(halfExt.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(halfExt.Y(), WithinAbs(REAL(1.5), 1e-10));
	}

	TEST_CASE("Box2D - Center and corners", "[Box2D][Corners]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 6);
		
		Point2Cartesian center = box.Center();
		REQUIRE_THAT(center.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(center.Y(), WithinAbs(REAL(3.0), 1e-10));
		
		REQUIRE_THAT(box.BottomLeft().X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box.BottomLeft().Y(), WithinAbs(REAL(0.0), 1e-10));
		
		REQUIRE_THAT(box.BottomRight().X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(box.BottomRight().Y(), WithinAbs(REAL(0.0), 1e-10));
		
		REQUIRE_THAT(box.TopLeft().X(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(box.TopLeft().Y(), WithinAbs(REAL(6.0), 1e-10));
		
		REQUIRE_THAT(box.TopRight().X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(box.TopRight().Y(), WithinAbs(REAL(6.0), 1e-10));
		
		auto corners = box.Corners();
		REQUIRE(corners.size() == 4);
	}

	TEST_CASE("Box2D - ToPolygon", "[Box2D][Polygon]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 3, 2);
		Polygon2D poly = box.ToPolygon();
		
		REQUIRE(poly.NumVertices() == 4);
		REQUIRE_THAT(poly.Area(), WithinAbs(REAL(6.0), 1e-10));
	}

	TEST_CASE("Box2D - Point containment", "[Box2D][Contains]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 3);
		
		// Inside
		REQUIRE(box.Contains(Point2Cartesian(2, 1.5)));
		
		// On boundary
		REQUIRE(box.Contains(Point2Cartesian(0, 0)));
		REQUIRE(box.Contains(Point2Cartesian(4, 3)));
		REQUIRE(box.Contains(Point2Cartesian(2, 0)));
		
		// Outside
		REQUIRE_FALSE(box.Contains(Point2Cartesian(-1, 1)));
		REQUIRE_FALSE(box.Contains(Point2Cartesian(5, 1)));
		REQUIRE_FALSE(box.Contains(Point2Cartesian(2, 4)));
	}

	TEST_CASE("Box2D - Box containment", "[Box2D][Contains]")
	{
		TEST_PRECISION_INFO();
		Box2D outer(0, 0, 10, 10);
		Box2D inner(2, 2, 6, 6);
		Box2D partial(5, 5, 15, 15);
		
		REQUIRE(outer.Contains(inner));
		REQUIRE_FALSE(outer.Contains(partial));
		REQUIRE_FALSE(inner.Contains(outer));
	}

	TEST_CASE("Box2D - Distance to point", "[Box2D][Distance]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 3);
		
		// Point inside - distance is 0
		REQUIRE_THAT(box.DistanceToPoint(Point2Cartesian(2, 1)), WithinAbs(REAL(0.0), 1e-10));
		
		// Point directly outside edge
		REQUIRE_THAT(box.DistanceToPoint(Point2Cartesian(6, 1)), WithinAbs(REAL(2.0), 1e-10));
		
		// Point outside corner
		REQUIRE_THAT(box.DistanceToPoint(Point2Cartesian(7, 7)), WithinAbs(REAL(5.0), 1e-10)); // 3-4-5 triangle
	}

	TEST_CASE("Box2D - Closest point", "[Box2D][ClosestPoint]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 3);
		
		// Point inside - returns the point itself (clamped)
		Point2Cartesian inside = box.ClosestPoint(Point2Cartesian(2, 1));
		REQUIRE_THAT(inside.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(inside.Y(), WithinAbs(REAL(1.0), 1e-10));
		
		// Point outside - projects to boundary
		Point2Cartesian outside = box.ClosestPoint(Point2Cartesian(6, 5));
		REQUIRE_THAT(outside.X(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(outside.Y(), WithinAbs(REAL(3.0), 1e-10));
	}

	TEST_CASE("Box2D - Box-box intersection test", "[Box2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		Box2D box1(0, 0, 4, 3);
		Box2D box2(2, 1, 6, 5);
		Box2D box3(10, 10, 12, 12);
		
		// Overlapping
		REQUIRE(box1.Intersects(box2));
		REQUIRE(box2.Intersects(box1));
		
		// Non-overlapping
		REQUIRE_FALSE(box1.Intersects(box3));
		
		// Touching (sharing an edge)
		Box2D box4(4, 0, 8, 3);
		REQUIRE(box1.Intersects(box4));
	}

	TEST_CASE("Box2D - Intersection calculation", "[Box2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		Box2D box1(0, 0, 4, 3);
		Box2D box2(2, 1, 6, 5);
		
		Box2D isect = box1.Intersection(box2);
		REQUIRE_THAT(isect.MinX(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(isect.MinY(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(isect.MaxX(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(isect.MaxY(), WithinAbs(REAL(3.0), 1e-10));
		
		// Non-overlapping returns empty box
		Box2D box3(10, 10, 12, 12);
		Box2D empty = box1.Intersection(box3);
		REQUIRE_THAT(empty.Area(), WithinAbs(REAL(0.0), 1e-10));
	}

	TEST_CASE("Box2D - Union calculation", "[Box2D][Union]")
	{
		TEST_PRECISION_INFO();
		Box2D box1(0, 0, 4, 3);
		Box2D box2(2, 1, 6, 5);
		
		Box2D u = box1.Union(box2);
		REQUIRE_THAT(u.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(u.MinY(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(u.MaxX(), WithinAbs(REAL(6.0), 1e-10));
		REQUIRE_THAT(u.MaxY(), WithinAbs(REAL(5.0), 1e-10));
	}

	TEST_CASE("Box2D - ExpandToInclude", "[Box2D][Expand]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 2, 2);
		
		box.ExpandToInclude(Point2Cartesian(5, 3));
		REQUIRE_THAT(box.MaxX(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(box.MaxY(), WithinAbs(REAL(3.0), 1e-10));
		
		box.ExpandToInclude(Point2Cartesian(-1, -1));
		REQUIRE_THAT(box.MinX(), WithinAbs(REAL(-1.0), 1e-10));
		REQUIRE_THAT(box.MinY(), WithinAbs(REAL(-1.0), 1e-10));
	}

	TEST_CASE("Box2D - Expanded margin", "[Box2D][Expand]")
	{
		TEST_PRECISION_INFO();
		Box2D box(1, 2, 3, 4);
		Box2D expanded = box.Expanded(1);
		
		REQUIRE_THAT(expanded.MinX(), WithinAbs(REAL(0.0), 1e-10));
		REQUIRE_THAT(expanded.MinY(), WithinAbs(REAL(1.0), 1e-10));
		REQUIRE_THAT(expanded.MaxX(), WithinAbs(REAL(4.0), 1e-10));
		REQUIRE_THAT(expanded.MaxY(), WithinAbs(REAL(5.0), 1e-10));
	}

	TEST_CASE("Box2D - IsDegenerate", "[Box2D][Degenerate]")
	{
		TEST_PRECISION_INFO();
		Box2D normal(0, 0, 2, 3);
		REQUIRE_FALSE(normal.IsDegenerate());
		
		Box2D flat(0, 0, 2, 0);
		REQUIRE(flat.IsDegenerate());
		
		Box2D point(1, 1, 1, 1);
		REQUIRE(point.IsDegenerate());
	}

	TEST_CASE("Box2D - Translation", "[Box2D][Transform]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 3);
		Box2D translated = box.Translated(Vector2Cartesian(5, 7));
		
		REQUIRE_THAT(translated.MinX(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(translated.MinY(), WithinAbs(REAL(7.0), 1e-10));
		REQUIRE_THAT(translated.MaxX(), WithinAbs(REAL(9.0), 1e-10));
		REQUIRE_THAT(translated.MaxY(), WithinAbs(REAL(10.0), 1e-10));
	}

	TEST_CASE("Box2D - Scaling", "[Box2D][Transform]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 2);
		Box2D scaled = box.Scaled(2);
		
		// Scales from center, so center stays the same
		Point2Cartesian center = scaled.Center();
		REQUIRE_THAT(center.X(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(center.Y(), WithinAbs(REAL(1.0), 1e-10));
		
		REQUIRE_THAT(scaled.Width(), WithinAbs(REAL(8.0), 1e-10));
		REQUIRE_THAT(scaled.Height(), WithinAbs(REAL(4.0), 1e-10));
	}

	TEST_CASE("Box2D - Circle intersection", "[Box2D][Intersection]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 4, 3);
		
		// Circle overlapping box
		Circle2D c1(2, 1.5, 1);
		REQUIRE(box.Intersects(c1));
		
		// Circle touching box edge
		Circle2D c2(5, 1.5, 1);
		REQUIRE(box.Intersects(c2));
		
		// Circle outside box
		Circle2D c3(10, 10, 1);
		REQUIRE_FALSE(box.Intersects(c3));
		
		// Circle containing box
		Circle2D c4(2, 1.5, 10);
		REQUIRE(box.Intersects(c4));
	}

	TEST_CASE("Box2D - Mutable accessors", "[Box2D][Accessors]")
	{
		TEST_PRECISION_INFO();
		Box2D box(0, 0, 1, 1);
		
		box.Min() = Point2Cartesian(2, 3);
		box.Max() = Point2Cartesian(5, 6);
		
		REQUIRE_THAT(box.MinX(), WithinAbs(REAL(2.0), 1e-10));
		REQUIRE_THAT(box.MinY(), WithinAbs(REAL(3.0), 1e-10));
		REQUIRE_THAT(box.MaxX(), WithinAbs(REAL(5.0), 1e-10));
		REQUIRE_THAT(box.MaxY(), WithinAbs(REAL(6.0), 1e-10));
	}

	/*********************************************************************/
	/*****           Line2D Additional Tests                         *****/
	/*********************************************************************/
	TEST_CASE("Line2D - Contains point", "[Line2D][Contains]")
	{
		TEST_PRECISION_INFO();
		Line2D line(Point2Cartesian(0, 0), Vector2Cartesian(1, 1));
		
		// Points on line
		REQUIRE(line.Contains(Point2Cartesian(0, 0)));
		REQUIRE(line.Contains(Point2Cartesian(5, 5)));
		REQUIRE(line.Contains(Point2Cartesian(-3, -3)));
		
		// Points off line
		REQUIRE_FALSE(line.Contains(Point2Cartesian(1, 0)));
		REQUIRE_FALSE(line.Contains(Point2Cartesian(0, 1)));
	}

	/*********************************************************************/
	/*****           Error handling tests                            *****/
	/*********************************************************************/
	TEST_CASE("Geometry2D - Degenerate triangle error handling", "[Triangle2D][Error]")
	{
		TEST_PRECISION_INFO();
		// Collinear points form degenerate triangle
		Triangle2D degen(Point2Cartesian(0, 0), Point2Cartesian(1, 0), Point2Cartesian(2, 0));
		
		// Circumcenter should throw for degenerate
		REQUIRE_THROWS(degen.Circumcenter());
		REQUIRE_THROWS(degen.CircumRadius());
	}

	TEST_CASE("Geometry2D - Empty polygon error handling", "[Polygon2D][Error]")
	{
		TEST_PRECISION_INFO();
		Polygon2D empty;
		
		REQUIRE_THROWS(empty.Centroid());
		REQUIRE_THROWS(empty.GetBoundingBox());
	}
}
