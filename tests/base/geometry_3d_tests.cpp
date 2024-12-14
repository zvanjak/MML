#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry3D.h"
#endif

using namespace MML;

namespace MML::Tests::Base::Geometry3DTests
{
	/*********************************************************************/
	/*****                       Line3D tests                        *****/
	/*********************************************************************/
	TEST_CASE("Line3D::init", "[Line3D]") {
		Point3Cartesian a(1, 1, 1);
		Vector3Cartesian b(1, 1, 1);

		REQUIRE(a.X() == 1.0);
		REQUIRE(a.Y() == 1.0);
		REQUIRE(a.Z() == 1.0);

		REQUIRE(b.X() == 1.0);
		REQUIRE(b.Y() == 1.0);
		REQUIRE(b.Z() == 1.0);
	}
	TEST_CASE("Line3D::PointOnLine", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::IsPerpendicular", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::IsParallel", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::DistToPoint", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::DistToLine", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::NearestPointOnLine", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::Intersection", "[Line3D]")
	{
	}
	TEST_CASE("Line3D::PerpendicularLineThroughPoint", "[Line3D]")
	{
	}

	/*********************************************************************/
	/*****                   SegmentLine3D tests                     *****/
	/*********************************************************************/

	TEST_CASE("SegmentLine3D::Init", "[SegmentLine3D]")
	{
	}
	TEST_CASE("SegmentLine3D::PointOnSegment", "[SegmentLine3D]")
	{
	}
	TEST_CASE("SegmentLine3D::Length", "[SegmentLine3D]")
	{
	}
	TEST_CASE("SegmentLine3D::Direction", "[SegmentLine3D]")
	{
	}

	/*********************************************************************/
	/*****                       Plane3D tests                       *****/
	/*********************************************************************/
	TEST_CASE("Plane3D::Init", "[Plane3D]") {
		Point3Cartesian a(1, 1, 1);
		Vector3Cartesian b(1, 1, 1);

		Plane3D plane1(a, b);

		REQUIRE(a.X() == 1.0);
		REQUIRE(a.Y() == 1.0);
		REQUIRE(a.Z() == 1.0);
	}
	TEST_CASE("Plane3D::GetPointOnPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::GetCoordAxisSegments", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::IsPointOnPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::DistToPoint", "[Plane3D]") 
	{
		Point3Cartesian a(1, 1, 1);
		Vector3Cartesian b(1, 1, 1);

		Plane3D planeXY = Plane3D::GetXYPlane();
		Plane3D planeXZ = Plane3D::GetXZPlane();
		Plane3D planeYZ = Plane3D::GetYZPlane();

		Point3Cartesian pnt(5, 5, 5);

		REQUIRE(planeXY.DistToPoint(pnt) == 5.0);
		REQUIRE(planeXZ.DistToPoint(pnt) == 5.0);
		REQUIRE(planeYZ.DistToPoint(pnt) == 5.0);
	}
	TEST_CASE("Plane3D::ProjectionToPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::IsLineOnPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::AngleToLine", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::IntersectionWithLine", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::IsParallelToPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::IsPerpendicularToPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::AngleToPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::DistToPlane", "[Plane3D]")
	{
	}
	TEST_CASE("Plane3D::IntersectionWithPlane", "[Plane3D]")
	{
	}
	
	/*********************************************************************/
	/*****                     Triangle3D tests                      *****/
	/*********************************************************************/


	/*********************************************************************/
	/*****                      Surfaces tests                       *****/
	/*********************************************************************/


	/*********************************************************************/
	/*****                       Bodies tests                        *****/
	/*********************************************************************/
}