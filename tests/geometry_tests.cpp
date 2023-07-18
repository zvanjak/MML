#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Geometry.h"
#include "basic_types/Geometry2D.h"
#include "basic_types/Geometry3D.h"
#endif


TEST_CASE("Test_Point2Cartesian", "[simple]") {
    MML::Point2Cartesian a(1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
}
// Dist
// oper + Vector
// oper - Point

TEST_CASE("Test_Vector2Cartesian_GetUnitVector", "[simple]") {
    MML::Vector2Cartesian a(4,0);
    
    auto unit = a.GetUnitVector();

    REQUIRE(unit.X() == 1.0);
    REQUIRE(unit.Y() == 0.0);
}

TEST_CASE("Test_Point3Cartesian", "[simple]") {
    MML::Point3Cartesian a(1,1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);
}
// Dist
// oper + Vector
// oper - Point

TEST_CASE("Test_Vector3Cartesian_init", "[simple]") {
    MML::Vector3Cartesian a(1,1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);

    MML::Vector3Cartesian b{2,2,2};
    
    REQUIRE(b.X() == 2.0);
    REQUIRE(b.Y() == 2.0);
    REQUIRE(b.Z() == 2.0);

    MML::VectorN<Real,3> c{3,3,3};
    MML::Vector3Cartesian d(c);

    REQUIRE(d.X() == 3.0);
    REQUIRE(d.Y() == 3.0);
    REQUIRE(d.Z() == 3.0);
}

TEST_CASE("Test_Vector3Cartesian_GetUnitVector", "[simple]") {
    MML::Vector3Cartesian a(4,0,0);
    
    auto unit = a.GetUnitVector();

    REQUIRE(unit.X() == 1.0);
    REQUIRE(unit.Y() == 0.0);
    REQUIRE(unit.Z() == 0.0);
}

// IsParallelTo
// IsPerpendicularTo
// ScalarProd
// VectorProd
// VectorsAngle

TEST_CASE("Test_Line2D_init", "[simple]") {
    MML::Point2Cartesian    start_pnt(1,1);
    MML::Vector2Cartesian   direction(1,0);
    
    MML::Line2D     line(start_pnt, direction);
    
    REQUIRE(line.StartPoint().X() == 1.0);
    REQUIRE(line.StartPoint().X() == 1.0);
    REQUIRE(line.Direction().X() == 1.0);
    REQUIRE(line.Direction().Y() == 0.0);

    MML::Point2Cartesian    pnt1(1,1);
    MML::Point2Cartesian    pnt2(2,2);

    MML::Line2D     line2(pnt1, pnt2);

}

// PointOnLine
// IsPerpendicular

TEST_CASE("Test_Line3D_init", "[simple]") {
    MML::Point3Cartesian a(1,1,1);
    MML::Vector3Cartesian b(1,1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);
}

// PointOnLine  
// IsPerpendicular

// SegmentLine2D
// SegmentLine3D

// Plane3D

TEST_CASE("Test_Plane3D_init", "[simple]") {
    MML::Point3Cartesian a(1,1,1);
    MML::Vector3Cartesian b(1,1,1);

    MML::Plane3D plane1(a, b);
    
    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);
}

TEST_CASE("Test_Plane3D_DistToPoint", "[simple]") {
    MML::Point3Cartesian a(1,1,1);
    MML::Vector3Cartesian b(1,1,1);

    MML::Plane3D planeXY = MML::Plane3D::GetXYPlane();
    MML::Plane3D planeXZ = MML::Plane3D::GetXZPlane();
    MML::Plane3D planeYZ = MML::Plane3D::GetYZPlane();
    
    MML::Point3Cartesian pnt(5, 5, 5);

    REQUIRE(planeXY.DistToPoint(pnt) == 5.0);
    REQUIRE(planeXZ.DistToPoint(pnt) == 5.0);
    REQUIRE(planeYZ.DistToPoint(pnt) == 5.0);
}