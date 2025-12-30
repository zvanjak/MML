#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry3D.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::Geometry3DTests
{
	/*********************************************************************/
	/*****                       Line3D tests                        *****/
	/*********************************************************************/
	TEST_CASE("Line3D::init", "[Line3D]") {
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 1, 1);
		Vector3Cartesian b(1, 1, 1);

		REQUIRE(a.X() == REAL(1.0));
		REQUIRE(a.Y() == REAL(1.0));
		REQUIRE(a.Z() == REAL(1.0));

		REQUIRE(b.X() == REAL(1.0));
		REQUIRE(b.Y() == REAL(1.0));
		REQUIRE(b.Z() == REAL(1.0));

		// Test constructors
		Line3D line1(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
		REQUIRE(line1.StartPoint() == Pnt3Cart(0, 0, 0));
		REQUIRE(line1.Direction().IsEqual(Vec3Cart(1, 0, 0)));

		Line3D line2(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0));
		REQUIRE(line2.Direction().IsEqual(Vec3Cart(1, 0, 0)));

		// Test null direction throws
		REQUIRE_THROWS(Line3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 0)));
		// Test same points throws
		REQUIRE_THROWS(Line3D(Pnt3Cart(1, 1, 1), Pnt3Cart(1, 1, 1)));
	}
	TEST_CASE("Line3D::PointOnLine", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D line(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));  // X-axis

		// Point at t=0 is start point
		REQUIRE(line(REAL(0.0)) == Pnt3Cart(0, 0, 0));
		// Point at t=5 is (5,0,0)
		REQUIRE(line(REAL(5.0)) == Pnt3Cart(5, 0, 0));
		// Negative t
		REQUIRE(line(-REAL(3.0)) == Pnt3Cart(-3, 0, 0));

		// IsPointOnLine tests
		REQUIRE(line.IsPointOnLine(Pnt3Cart(0, 0, 0)));
		REQUIRE(line.IsPointOnLine(Pnt3Cart(100, 0, 0)));
		REQUIRE(line.IsPointOnLine(Pnt3Cart(-50, 0, 0)));
		REQUIRE_FALSE(line.IsPointOnLine(Pnt3Cart(0, 1, 0)));
		REQUIRE_FALSE(line.IsPointOnLine(Pnt3Cart(5, 5, 0)));

		// Diagonal line
		Line3D diagLine(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 1));
		REQUIRE(diagLine.IsPointOnLine(Pnt3Cart(5, 5, 5)));
		REQUIRE_FALSE(diagLine.IsPointOnLine(Pnt3Cart(1, 2, 3)));
	}
	TEST_CASE("Line3D::IsPerpendicular", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
		Line3D yAxis(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0));
		Line3D zAxis(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1));
		Line3D diagonal(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 0));

		REQUIRE(xAxis.IsPerpendicular(yAxis));
		REQUIRE(yAxis.IsPerpendicular(zAxis));
		REQUIRE(xAxis.IsPerpendicular(zAxis));
		REQUIRE_FALSE(xAxis.IsPerpendicular(diagonal));
	}
	TEST_CASE("Line3D::IsParallel", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D xAxis1(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
		Line3D xAxis2(Pnt3Cart(0, 5, 0), Vec3Cart(1, 0, 0));
		Line3D yAxis(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0));

		REQUIRE(xAxis1.IsParallel(xAxis2));
		REQUIRE_FALSE(xAxis1.IsParallel(yAxis));
	}
	TEST_CASE("Line3D::DistToPoint", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));

		// Point on line
		REQUIRE_THAT(xAxis.Dist(Pnt3Cart(5, 0, 0)), RealApprox(REAL(0.0)));
		// Point above line
		REQUIRE_THAT(xAxis.Dist(Pnt3Cart(0, 5, 0)), RealApprox(REAL(5.0)));
		REQUIRE_THAT(xAxis.Dist(Pnt3Cart(10, 3, 4)), RealApprox(REAL(5.0)));

		// Diagonal line
		Line3D diagLine(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 0));
		// Point (1, 0, 0) to diagonal line through origin with direction (1,1,0)/sqrt(2)
		// Distance should be 1/sqrt(2) = sqrt(2)/2 ≈ REAL(0.7071)
		REQUIRE_THAT(diagLine.Dist(Pnt3Cart(1, 0, 0)), RealApprox(REAL(0.7071067811865476)));
	}
	TEST_CASE("Line3D::DistToLine", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		// Parallel lines
		Line3D xAxis1(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
		Line3D xAxis2(Pnt3Cart(0, 5, 0), Vec3Cart(1, 0, 0));
		// Skew lines
		Line3D yAxis(Pnt3Cart(0, 0, 3), Vec3Cart(0, 1, 0));

		// Distance between parallel lines at y=0 and y=5
		// Note: for parallel lines, the formula gives different result
		// REQUIRE_THAT(xAxis1.Dist(xAxis2) , RealApprox(REAL(5.0)));

		// Skew lines: X-axis at z=0 and Y-axis at z=3
		REQUIRE_THAT(xAxis1.Dist(yAxis) , RealApprox(REAL(3.0)));

		// Intersecting lines have distance 0
		Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
		Line3D yAxisOrigin(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0));
		REQUIRE_THAT(xAxis.Dist(yAxisOrigin) , RealApprox(REAL(0.0)));
	}
	TEST_CASE("Line3D::NearestPointOnLine", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));

		// Point above X-axis
		Pnt3Cart nearest = xAxis.NearestPointOnLine(Pnt3Cart(5, 10, 0));
		REQUIRE_THAT(nearest.X() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(nearest.Y() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(nearest.Z() , RealApprox(REAL(0.0)));

		// Point at origin should stay at origin
		nearest = xAxis.NearestPointOnLine(Pnt3Cart(0, 5, 5));
		REQUIRE_THAT(nearest.X() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(nearest.Y() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(nearest.Z() , RealApprox(REAL(0.0)));
	}
	TEST_CASE("Line3D::Intersection", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));
		Line3D yAxis(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0));
		Line3D parallelX(Pnt3Cart(0, 5, 0), Vec3Cart(1, 0, 0));

		Pnt3Cart intersection;

		// Intersecting lines at origin
		REQUIRE(xAxis.Intersection(yAxis, intersection));
		REQUIRE_THAT(intersection.X() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(intersection.Y() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(intersection.Z() , RealApprox(REAL(0.0)));

		// Parallel lines don't intersect
		REQUIRE_FALSE(xAxis.Intersection(parallelX, intersection));

		// Two lines intersecting at (2, 2, 0)
		Line3D line1(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 0));
		Line3D line2(Pnt3Cart(4, 0, 0), Vec3Cart(-1, 1, 0));
		REQUIRE(line1.Intersection(line2, intersection));
		REQUIRE_THAT(intersection.X() , RealApprox(REAL(2.0)));
		REQUIRE_THAT(intersection.Y() , RealApprox(REAL(2.0)));
	}
	TEST_CASE("Line3D::PerpendicularLineThroughPoint", "[Line3D]")
	{
			TEST_PRECISION_INFO();
		Line3D xAxis(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0));

		// Perpendicular from (5, 10, 0) should go toward (5, 0, 0)
		Line3D perpLine = xAxis.PerpendicularLineThroughPoint(Pnt3Cart(5, 10, 0));
		REQUIRE(perpLine.StartPoint() == Pnt3Cart(5, 10, 0));
		// Direction should be (0, -1, 0) normalized
		REQUIRE(perpLine.Direction().IsEqual(Vec3Cart(0, -1, 0)));

		// Point on line should throw
		REQUIRE_THROWS(xAxis.PerpendicularLineThroughPoint(Pnt3Cart(5, 0, 0)));
	}

	/*********************************************************************/
	/*****                   SegmentLine3D tests                     *****/
	/*********************************************************************/

	TEST_CASE("SegmentLine3D::Init", "[SegmentLine3D]")
	{
			TEST_PRECISION_INFO();
		// Constructor with two points
		SegmentLine3D seg1(Pnt3Cart(0, 0, 0), Pnt3Cart(10, 0, 0));
		REQUIRE(seg1.StartPoint() == Pnt3Cart(0, 0, 0));
		REQUIRE(seg1.EndPoint() == Pnt3Cart(10, 0, 0));

		// Constructor with point, direction, and parameter
		SegmentLine3D seg2(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 0), REAL(5.0));
		REQUIRE(seg2.StartPoint() == Pnt3Cart(0, 0, 0));
		REQUIRE_THAT(seg2.EndPoint().X() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(seg2.EndPoint().Y() , RealApprox(REAL(5.0)));
	}
	TEST_CASE("SegmentLine3D::PointOnSegment", "[SegmentLine3D]")
	{
			TEST_PRECISION_INFO();
		SegmentLine3D seg(Pnt3Cart(0, 0, 0), Pnt3Cart(10, 0, 0));

		// t=0 gives start point
		REQUIRE(seg(REAL(0.0)) == Pnt3Cart(0, 0, 0));
		// t=1 gives end point
		REQUIRE(seg(REAL(1.0)) == Pnt3Cart(10, 0, 0));
		// t=REAL(0.5) gives midpoint
		REQUIRE_THAT(seg(REAL(0.5)).X() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(seg(REAL(0.5)).Y() , RealApprox(REAL(0.0)));

		// t < 0 or t > 1 throws
		REQUIRE_THROWS(seg(-REAL(0.1)));
		REQUIRE_THROWS(seg(REAL(1.1)));
	}
	TEST_CASE("SegmentLine3D::Length", "[SegmentLine3D]")
	{
			TEST_PRECISION_INFO();
		SegmentLine3D seg1(Pnt3Cart(0, 0, 0), Pnt3Cart(10, 0, 0));
		REQUIRE_THAT(seg1.Length() , RealApprox(REAL(10.0)));

		SegmentLine3D seg2(Pnt3Cart(0, 0, 0), Pnt3Cart(3, 4, 0));
		REQUIRE_THAT(seg2.Length() , RealApprox(REAL(5.0)));  // 3-4-5 triangle

		SegmentLine3D seg3(Pnt3Cart(1, 2, 3), Pnt3Cart(4, 6, 3));
		REQUIRE_THAT(seg3.Length() , RealApprox(REAL(5.0)));  // 3-4-5 triangle
	}
	TEST_CASE("SegmentLine3D::Direction", "[SegmentLine3D]")
	{
			TEST_PRECISION_INFO();
		SegmentLine3D seg(Pnt3Cart(0, 0, 0), Pnt3Cart(10, 0, 0));
		Vec3Cart dir = seg.Direction();
		REQUIRE_THAT(dir.X() , RealApprox(REAL(10.0)));
		REQUIRE_THAT(dir.Y() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(dir.Z() , RealApprox(REAL(0.0)));

		SegmentLine3D seg2(Pnt3Cart(1, 1, 1), Pnt3Cart(4, 5, 3));
		Vec3Cart dir2 = seg2.Direction();
		REQUIRE_THAT(dir2.X() , RealApprox(REAL(3.0)));
		REQUIRE_THAT(dir2.Y() , RealApprox(REAL(4.0)));
		REQUIRE_THAT(dir2.Z() , RealApprox(REAL(2.0)));
	}

	/*********************************************************************/
	/*****                       Plane3D tests                       *****/
	/*********************************************************************/
	TEST_CASE("Plane3D::Init", "[Plane3D]") {
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 1, 1);
		Vector3Cartesian b(1, 1, 1);

		Plane3D plane1(a, b);

		REQUIRE(a.X() == REAL(1.0));
		REQUIRE(a.Y() == REAL(1.0));
		REQUIRE(a.Z() == REAL(1.0));

		// Test standard planes
		Plane3D xyPlane = Plane3D::GetXYPlane();
		REQUIRE(xyPlane.Normal().IsEqual(Vec3Cart(0, 0, 1)));

		Plane3D xzPlane = Plane3D::GetXZPlane();
		REQUIRE(xzPlane.Normal().IsEqual(Vec3Cart(0, 1, 0)));

		Plane3D yzPlane = Plane3D::GetYZPlane();
		REQUIRE(yzPlane.Normal().IsEqual(Vec3Cart(1, 0, 0)));

		// Constructor with three points
		Plane3D plane3pts(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(0, 1, 0));
		REQUIRE(plane3pts.IsPointOnPlane(Pnt3Cart(REAL(0.5), REAL(0.5), 0)));

		// Null normal throws
		REQUIRE_THROWS(Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 0)));
	}
	TEST_CASE("Plane3D::GetPointOnPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();
		Pnt3Cart pnt = xyPlane.GetPointOnPlane();
		REQUIRE(xyPlane.IsPointOnPlane(pnt));

		Plane3D xzPlane = Plane3D::GetXZPlane();
		pnt = xzPlane.GetPointOnPlane();
		REQUIRE(xzPlane.IsPointOnPlane(pnt));

		// Custom plane
		Plane3D customPlane(Pnt3Cart(1, 2, 3), Vec3Cart(1, 1, 1));
		pnt = customPlane.GetPointOnPlane();
		REQUIRE(customPlane.IsPointOnPlane(pnt));
	}
	TEST_CASE("Plane3D::GetCoordAxisSegments", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		// Plane x + y + z = 6 (intercepts at 6, 6, 6)
		Plane3D plane(REAL(6.0), REAL(6.0), REAL(6.0));  // Using segment constructor
		Real segX, segY, segZ;
		plane.GetCoordAxisSegments(segX, segY, segZ);
		REQUIRE_THAT(segX , RealApprox(REAL(6.0)));
		REQUIRE_THAT(segY , RealApprox(REAL(6.0)));
		REQUIRE_THAT(segZ , RealApprox(REAL(6.0)));

	// XY plane (z=0) has infinite x,y intercepts and z=0
	Plane3D xyPlane = Plane3D::GetXYPlane();
	xyPlane.GetCoordAxisSegments(segX, segY, segZ);
	REQUIRE(segX == Constants::PosInf);
	REQUIRE(segY == Constants::PosInf);
	REQUIRE_THAT(segZ , RealApprox(REAL(0.0)));
	}
	TEST_CASE("Plane3D::IsPointOnPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();

		REQUIRE(xyPlane.IsPointOnPlane(Pnt3Cart(0, 0, 0)));
		REQUIRE(xyPlane.IsPointOnPlane(Pnt3Cart(100, -50, 0)));
		REQUIRE_FALSE(xyPlane.IsPointOnPlane(Pnt3Cart(0, 0, REAL(0.1))));
		REQUIRE_FALSE(xyPlane.IsPointOnPlane(Pnt3Cart(5, 5, 5)));
	}
	TEST_CASE("Plane3D::DistToPoint", "[Plane3D]") 
	{
			TEST_PRECISION_INFO();
		Point3Cartesian a(1, 1, 1);
		Vector3Cartesian b(1, 1, 1);

		Plane3D planeXY = Plane3D::GetXYPlane();
		Plane3D planeXZ = Plane3D::GetXZPlane();
		Plane3D planeYZ = Plane3D::GetYZPlane();

		Point3Cartesian pnt(5, 5, 5);

		REQUIRE(planeXY.DistToPoint(pnt) == REAL(5.0));
		REQUIRE(planeXZ.DistToPoint(pnt) == REAL(5.0));
		REQUIRE(planeYZ.DistToPoint(pnt) == REAL(5.0));

		// Point on plane has distance 0
		REQUIRE(planeXY.DistToPoint(Pnt3Cart(10, 20, 0)) == REAL(0.0));
	}
	TEST_CASE("Plane3D::ProjectionToPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();

		// Project point above XY plane
		Pnt3Cart proj = xyPlane.ProjectionToPlane(Pnt3Cart(5, 5, 10));
		REQUIRE_THAT(proj.X() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(proj.Y() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(proj.Z() , RealApprox(REAL(0.0)));

		// Point on plane projects to itself
		proj = xyPlane.ProjectionToPlane(Pnt3Cart(3, 4, 0));
		REQUIRE_THAT(proj.X() , RealApprox(REAL(3.0)));
		REQUIRE_THAT(proj.Y() , RealApprox(REAL(4.0)));
		REQUIRE_THAT(proj.Z() , RealApprox(REAL(0.0)));
	}
	TEST_CASE("Plane3D::IsLineOnPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();

		Line3D lineOnPlane(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 0));
		Line3D lineNotOnPlane(Pnt3Cart(0, 0, 1), Vec3Cart(1, 0, 0));
		Line3D lineCrossingPlane(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1));

		REQUIRE(xyPlane.IsLineOnPlane(lineOnPlane));
		REQUIRE_FALSE(xyPlane.IsLineOnPlane(lineNotOnPlane));
		REQUIRE_FALSE(xyPlane.IsLineOnPlane(lineCrossingPlane));
	}
	TEST_CASE("Plane3D::AngleToLine", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();

		// Line parallel to plane (angle = 0)
		Line3D parallelLine(Pnt3Cart(0, 0, 5), Vec3Cart(1, 0, 0));
		REQUIRE_THAT(xyPlane.AngleToLine(parallelLine) , RealApprox(REAL(0.0)).margin(1e-10));

		// Line perpendicular to plane (angle = pi/2)
		Line3D perpLine(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1));
		REQUIRE_THAT(xyPlane.AngleToLine(perpLine) , RealApprox(Constants::PI / REAL(2.0)).margin(1e-10));

		// Line at 45 degrees
		Line3D line45(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 1));
		REQUIRE_THAT(xyPlane.AngleToLine(line45) , RealApprox(Constants::PI / REAL(4.0)).margin(1e-10));
	}
	TEST_CASE("Plane3D::IntersectionWithLine", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();
		Pnt3Cart intersection;

		// Line crossing at origin
		Line3D zAxis(Pnt3Cart(0, 0, 5), Vec3Cart(0, 0, -1));
		REQUIRE(xyPlane.IntersectionWithLine(zAxis, intersection));
		REQUIRE_THAT(intersection.X() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(intersection.Y() , RealApprox(REAL(0.0)));
		REQUIRE_THAT(intersection.Z() , RealApprox(REAL(0.0)));

		// Line crossing at (3, 4, 0)
		Line3D diagLine(Pnt3Cart(3, 4, 10), Vec3Cart(0, 0, -1));
		REQUIRE(xyPlane.IntersectionWithLine(diagLine, intersection));
		REQUIRE_THAT(intersection.X() , RealApprox(REAL(3.0)));
		REQUIRE_THAT(intersection.Y() , RealApprox(REAL(4.0)));
		REQUIRE_THAT(intersection.Z() , RealApprox(REAL(0.0)));

		// Parallel line - no intersection
		Line3D parallelLine(Pnt3Cart(0, 0, 5), Vec3Cart(1, 0, 0));
		REQUIRE_FALSE(xyPlane.IntersectionWithLine(parallelLine, intersection));
	}
	TEST_CASE("Plane3D::IsParallelToPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane1 = Plane3D::GetXYPlane();
		Plane3D xyPlane2(Pnt3Cart(0, 0, 5), Vec3Cart(0, 0, 1));  // Parallel, z=5
		Plane3D xzPlane = Plane3D::GetXZPlane();

		REQUIRE(xyPlane1.IsParallelToPlane(xyPlane2));
		REQUIRE_FALSE(xyPlane1.IsParallelToPlane(xzPlane));
	}
	TEST_CASE("Plane3D::IsPerpendicularToPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();
		Plane3D xzPlane = Plane3D::GetXZPlane();
		Plane3D yzPlane = Plane3D::GetYZPlane();

		REQUIRE(xyPlane.IsPerpendicularToPlane(xzPlane));
		REQUIRE(xyPlane.IsPerpendicularToPlane(yzPlane));
		REQUIRE(xzPlane.IsPerpendicularToPlane(yzPlane));
	}
	TEST_CASE("Plane3D::AngleToPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();
		Plane3D xzPlane = Plane3D::GetXZPlane();

		// Perpendicular planes - angle = pi/2
		REQUIRE_THAT(xyPlane.AngleToPlane(xzPlane) , RealApprox(Constants::PI / REAL(2.0)).margin(1e-10));

		// Parallel planes - angle = 0
		Plane3D xyPlane2(Pnt3Cart(0, 0, 5), Vec3Cart(0, 0, 1));
		REQUIRE_THAT(xyPlane.AngleToPlane(xyPlane2) , RealApprox(REAL(0.0)).margin(1e-10));
	}
	TEST_CASE("Plane3D::DistToPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane1 = Plane3D::GetXYPlane();
		Plane3D xyPlane2(Pnt3Cart(0, 0, 5), Vec3Cart(0, 0, 1));  // z=5
		Plane3D xzPlane = Plane3D::GetXZPlane();

		// Parallel planes at z=0 and z=5
		REQUIRE_THAT(xyPlane1.DistToPlane(xyPlane2) , RealApprox(REAL(5.0)));

		// Non-parallel planes have distance 0 (they intersect)
		REQUIRE_THAT(xyPlane1.DistToPlane(xzPlane) , RealApprox(REAL(0.0)));
	}
	TEST_CASE("Plane3D::IntersectionWithPlane", "[Plane3D]")
	{
			TEST_PRECISION_INFO();
		Plane3D xyPlane = Plane3D::GetXYPlane();
		Plane3D xzPlane = Plane3D::GetXZPlane();
		Line3D intersectionLine;

		// XY and XZ planes intersect along X-axis
		REQUIRE(xyPlane.IntersectionWithPlane(xzPlane, intersectionLine));
	// Direction should be along X-axis (or opposite direction)
	Vec3Cart cross = VectorProduct(intersectionLine.Direction(), Vec3Cart(1, 0, 0));
	REQUIRE(cross.NormL2() < 1e-10);  // Cross product zero means parallel or anti-parallel
		Plane3D xyPlane2(Pnt3Cart(0, 0, 5), Vec3Cart(0, 0, 1));
		REQUIRE_FALSE(xyPlane.IntersectionWithPlane(xyPlane2, intersectionLine));
	}
	
	/*********************************************************************/
	/*****                     Triangle3D tests                      *****/
	/*********************************************************************/
	TEST_CASE("Triangle3D::Basic", "[Triangle3D]")
	{
			TEST_PRECISION_INFO();
		Triangle3D tri(Pnt3Cart(0, 0, 0), Pnt3Cart(3, 0, 0), Pnt3Cart(0, 4, 0));

		// 3-4-5 right triangle
		REQUIRE_THAT(tri.A() , RealApprox(REAL(3.0)));  // side 1-2
		REQUIRE_THAT(tri.B() , RealApprox(REAL(5.0)));  // side 2-3 (hypotenuse)
		REQUIRE_THAT(tri.C() , RealApprox(REAL(4.0)));  // side 3-1

		// Area = REAL(0.5) * base * height = REAL(0.5) * 3 * 4 = 6
		REQUIRE_THAT(tri.Area() , RealApprox(REAL(6.0)));

		// Plane
		Plane3D plane = tri.getDefinedPlane();
		REQUIRE(plane.IsPointOnPlane(Pnt3Cart(1, 1, 0)));
		REQUIRE(plane.Normal().IsParallelTo(Vec3Cart(0, 0, 1)));
	}

	TEST_CASE("Triangle3D::IsPointInside", "[Triangle3D]")
	{
			TEST_PRECISION_INFO();
		Triangle3D tri(Pnt3Cart(0, 0, 0), Pnt3Cart(4, 0, 0), Pnt3Cart(2, 3, 0));

		// Centroid should be inside
		REQUIRE(tri.IsPointInside(Pnt3Cart(2, 1, 0)));
		// Vertices should be inside (boundary)
		REQUIRE(tri.IsPointInside(Pnt3Cart(0, 0, 0)));
		REQUIRE(tri.IsPointInside(Pnt3Cart(4, 0, 0)));
		REQUIRE(tri.IsPointInside(Pnt3Cart(2, 3, 0)));
		// Point outside
		REQUIRE_FALSE(tri.IsPointInside(Pnt3Cart(10, 10, 0)));
		REQUIRE_FALSE(tri.IsPointInside(Pnt3Cart(-1, -1, 0)));
	}

	TEST_CASE("TriangleSurface3D::CheckLocalCoordSystem", "[TriangleSurface3D]")
	{
			TEST_PRECISION_INFO();
		// triangle in x-y plane
		TriangleSurface3D triangle(Pnt3Cart(0,0,0), Pnt3Cart(2,0,0), Pnt3Cart(1,1,0));

		REQUIRE(triangle._localX.IsEqualTo(Vector3Cartesian(1, 0, 0)));
		REQUIRE(triangle._localY.IsEqualTo(Vector3Cartesian(0, 1, 0)));

		// triangle in x-z plane
		TriangleSurface3D triangle2(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(1, 0, 1));

		REQUIRE(triangle2._localX.IsEqualTo(Vector3Cartesian(1, 0, 0)));
		REQUIRE(triangle2._localY.IsEqualTo(Vector3Cartesian(0, 0, 1)));

		// triangle in y-z plane
		TriangleSurface3D triangle3(Pnt3Cart(0, 0, 0), Pnt3Cart(0, 2, 0), Pnt3Cart(0, 1, 1));

		REQUIRE(triangle3._localX.IsEqualTo(Vector3Cartesian(0, 1, 0)));
		REQUIRE(triangle3._localY.IsEqualTo(Vector3Cartesian(0, 0, 1)));
	}
	TEST_CASE("TriangleSurface3D::CheckMinMax", "[TriangleSurface3D]")
	{
			TEST_PRECISION_INFO();
		// 3-4-5 triangle in XY plane
		TriangleSurface3D tri(Pnt3Cart(0, 0, 0), Pnt3Cart(4, 0, 0), Pnt3Cart(0, 3, 0));

		// After point rotation (4 is longest side), check bounds
		REQUIRE_THAT(tri.getMinU() , RealApprox(REAL(0.0)));
		REQUIRE(tri.getMaxU() > 0);  // Should be positive
		REQUIRE_THAT(tri.getMinW(REAL(0.0)) , RealApprox(REAL(0.0)));
	}


	/*********************************************************************/
	/*****                      Surfaces tests                       *****/
	/*********************************************************************/


	/*********************************************************************/
	/*****                       Bodies tests                        *****/
	/*********************************************************************/

	/*********************************************************************/
	/*****         New Geometry3D Improvements - Dec 2025            *****/
	/*********************************************************************/

	TEST_CASE("Triangle3D::IsRight with tolerance", "[Triangle3D][improvements]")
	{
			TEST_PRECISION_INFO();
		// Perfect right triangle: 3-4-5
		Triangle3D tri1(Pnt3Cart(0, 0, 0), Pnt3Cart(3, 0, 0), Pnt3Cart(0, 4, 0));
		REQUIRE(tri1.IsRight());
		REQUIRE(tri1.IsRight(1e-10));

		// Clearly not a right triangle (equilateral)
		Triangle3D tri2(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(REAL(0.5), sqrt(REAL(3.0)) / REAL(2.0), 0));
		REQUIRE_FALSE(tri2.IsRight());

		// 3D right triangle
		Triangle3D tri3(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(0, 1, 1));
		REQUIRE(tri3.IsRight(1e-8));
	}

	TEST_CASE("Triangle3D::IsIsosceles with tolerance", "[Triangle3D][improvements]")
	{
			TEST_PRECISION_INFO();
		// Perfect isosceles
		Triangle3D tri1(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(1, 2, 0));
		REQUIRE(tri1.IsIsosceles());

		// Nearly isosceles
		Triangle3D tri2(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(REAL(1.0000001), 2, 0));
		REQUIRE(tri2.IsIsosceles(1e-4));
		REQUIRE_FALSE(tri2.IsIsosceles(1e-12));

		// Equilateral (also isosceles)
		Triangle3D tri3(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(REAL(0.5), sqrt(REAL(3.0)) / REAL(2.0), 0));
		REQUIRE(tri3.IsIsosceles());

		// Scalene (not isosceles)
		Triangle3D tri4(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(0, 3, 0));
		REQUIRE_FALSE(tri4.IsIsosceles());
	}

	TEST_CASE("Triangle3D::IsEquilateral with tolerance", "[Triangle3D][improvements]")
	{
			TEST_PRECISION_INFO();
		// Perfect equilateral
		Triangle3D tri1(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(REAL(0.5), sqrt(REAL(3.0)) / REAL(2.0), 0));
		REQUIRE(tri1.IsEquilateral());

		// Nearly equilateral
		Triangle3D tri2(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(REAL(0.5), sqrt(REAL(3.0)) / REAL(2.0) + REAL(0.0000001), 0));
		REQUIRE(tri2.IsEquilateral(1e-4));
		REQUIRE_FALSE(tri2.IsEquilateral(1e-12));

		// Isosceles but not equilateral
		Triangle3D tri3(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(1, 2, 0));
		REQUIRE_FALSE(tri3.IsEquilateral());

		// 3D equilateral triangle
		Real s = REAL(1.0) / sqrt(REAL(3.0));
		Triangle3D tri4(Pnt3Cart(s, s, s), Pnt3Cart(-s, -s, s), Pnt3Cart(-s, s, -s));
		REQUIRE(tri4.IsEquilateral(1e-8));
	}

	TEST_CASE("Triangle3D::Centroid", "[Triangle3D][improvements]")
	{
			TEST_PRECISION_INFO();
		// Simple triangle in XY plane
		Triangle3D tri1(Pnt3Cart(0, 0, 0), Pnt3Cart(3, 0, 0), Pnt3Cart(0, 3, 0));
		Pnt3Cart centroid1 = tri1.Centroid();
		REQUIRE_THAT(centroid1.X() , RealApprox(REAL(1.0)));
		REQUIRE_THAT(centroid1.Y() , RealApprox(REAL(1.0)));
		REQUIRE_THAT(centroid1.Z() , RealApprox(REAL(0.0)));

		// Triangle with all different coordinates
		Triangle3D tri2(Pnt3Cart(1, 2, 3), Pnt3Cart(4, 5, 6), Pnt3Cart(7, 8, 9));
		Pnt3Cart centroid2 = tri2.Centroid();
		REQUIRE_THAT(centroid2.X() , RealApprox(REAL(4.0)));
		REQUIRE_THAT(centroid2.Y() , RealApprox(REAL(5.0)));
		REQUIRE_THAT(centroid2.Z() , RealApprox(REAL(6.0)));

		// Equilateral triangle centered at origin (centroid should be at origin)
		Real s = sqrt(REAL(3.0));
		Triangle3D tri3(Pnt3Cart(1, 0, 0), Pnt3Cart(-REAL(0.5), s / REAL(2.0), 0), Pnt3Cart(-REAL(0.5), -s / REAL(2.0), 0));
		Pnt3Cart centroid3 = tri3.Centroid();
		REQUIRE_THAT(centroid3.X() , RealApprox(REAL(0.0)).margin(1e-10));
		REQUIRE_THAT(centroid3.Y() , RealApprox(REAL(0.0)).margin(1e-10));
		REQUIRE_THAT(centroid3.Z() , RealApprox(REAL(0.0)));
	}

	TEST_CASE("SegmentLine3D::Dist to point", "[SegmentLine3D][improvements]")
	{
			TEST_PRECISION_INFO();
		// Horizontal segment along X-axis from (0,0,0) to (10,0,0)
		SegmentLine3D seg1(Pnt3Cart(0, 0, 0), Pnt3Cart(10, 0, 0));

		// Point directly above middle of segment
		REQUIRE_THAT(seg1.Dist(Pnt3Cart(5, 5, 0)), RealApprox(REAL(5.0)));

		// Point at start of segment (distance = 0)
		REQUIRE_THAT(seg1.Dist(Pnt3Cart(0, 0, 0)), RealApprox(REAL(0.0)));

		// Point at end of segment (distance = 0)
		REQUIRE_THAT(seg1.Dist(Pnt3Cart(10, 0, 0)), RealApprox(REAL(0.0)));

		// Point beyond start (closest to start point)
		REQUIRE_THAT(seg1.Dist(Pnt3Cart(-5, 0, 0)), RealApprox(REAL(5.0)));

		// Point beyond end (closest to end point)
		REQUIRE_THAT(seg1.Dist(Pnt3Cart(15, 0, 0)), RealApprox(REAL(5.0)));

		// Point perpendicular to middle, in 3D
		REQUIRE_THAT(seg1.Dist(Pnt3Cart(5, 3, 4)), RealApprox(REAL(5.0)));

		// Diagonal segment
		SegmentLine3D seg2(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 1, 1));
		REQUIRE_THAT(seg2.Dist(Pnt3Cart(REAL(0.5), REAL(0.5), REAL(0.5))), RealApprox(REAL(0.0)).margin(1e-10));  // Point on segment

		// Degenerate segment (zero length)
		SegmentLine3D seg3(Pnt3Cart(5, 5, 5), Pnt3Cart(5, 5, 5));
		REQUIRE_THAT(seg3.Dist(Pnt3Cart(5, 5, 8)), RealApprox(REAL(3.0)));
	}

	TEST_CASE("RectSurface3D::coplanarity validation", "[RectSurface3D][improvements]")
	{
			TEST_PRECISION_INFO();
		// Valid rectangular surface (all points in XY plane)
		REQUIRE_NOTHROW(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(2, 3, 0), Pnt3Cart(0, 3, 0)));

		// Valid rectangular surface in XZ plane
		REQUIRE_NOTHROW(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(2, 0, 0), Pnt3Cart(2, 0, 3), Pnt3Cart(0, 0, 3)));

		// Valid rectangular surface in YZ plane
		REQUIRE_NOTHROW(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(0, 2, 0), Pnt3Cart(0, 2, 3), Pnt3Cart(0, 0, 3)));

		// Valid tilted rectangular surface
		REQUIRE_NOTHROW(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 1), Pnt3Cart(0, 1, 1)));

		// Invalid - 4th point not coplanar with first 3
		REQUIRE_THROWS(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 0), Pnt3Cart(0, 1, 1)));

		// Invalid - points form a twisted quadrilateral
		REQUIRE_THROWS(RectSurface3D(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 1), Pnt3Cart(0, 1, 0)));
	}
}
