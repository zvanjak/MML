#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry.h"
#include "base/Geometry2D.h"
#include "base/Geometry3D.h"
#endif

using namespace MML;

// TODO - BIG!!! - add tests for all functions
TEST_CASE("Test_Point2Cartesian", "[simple]") {
    Point2Cartesian a(1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
}
// Dist
// oper + Vector
// oper - Point

TEST_CASE("Test_Vector2Cartesian_GetUnitVector", "[simple]") {
    Vector2Cartesian a(4,0);
    
    auto unit = a.GetUnitVector();

    REQUIRE(unit.X() == 1.0);
    REQUIRE(unit.Y() == 0.0);
}

TEST_CASE("Test_Point3Cartesian", "[simple]") {
    Point3Cartesian a(1,1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);
}
// Dist
// oper + Vector
// oper - Point

TEST_CASE("Test_Vector3Cartesian_init", "[simple]") {
    Vector3Cartesian a(1,1,1);

    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);

    Vector3Cartesian b{2,2,2};
    
    REQUIRE(b.X() == 2.0);
    REQUIRE(b.Y() == 2.0);
    REQUIRE(b.Z() == 2.0);

    VectorN<Real,3> c{3,3,3};
    Vector3Cartesian d(c);

    REQUIRE(d.X() == 3.0);
    REQUIRE(d.Y() == 3.0);
    REQUIRE(d.Z() == 3.0);
}

TEST_CASE("Test_Vector3Cartesian_GetUnitVector", "[simple]") {
    Vector3Cartesian a(4,0,0);
    
    auto unit = a.GetAsUnitVector();

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
    Point2Cartesian    start_pnt(1,1);
    Vector2Cartesian   direction(1,0);
    
    Line2D     line(start_pnt, direction);
    
    REQUIRE(line.StartPoint().X() == 1.0);
    REQUIRE(line.StartPoint().X() == 1.0);
    REQUIRE(line.Direction().X() == 1.0);
    REQUIRE(line.Direction().Y() == 0.0);

    Point2Cartesian    pnt1(1,1);
    Point2Cartesian    pnt2(2,2);

    Line2D     line2(pnt1, pnt2);

}

// PointOnLine
// IsPerpendicular

TEST_CASE("Test_Line3D_init", "[simple]") {
    Point3Cartesian a(1,1,1);
    Vector3Cartesian b(1,1,1);

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
    Point3Cartesian a(1,1,1);
    Vector3Cartesian b(1,1,1);

    Plane3D plane1(a, b);
    
    REQUIRE(a.X() == 1.0);
    REQUIRE(a.Y() == 1.0);
    REQUIRE(a.Z() == 1.0);
}

TEST_CASE("Test_Plane3D_DistToPoint", "[simple]") {
    Point3Cartesian a(1,1,1);
    Vector3Cartesian b(1,1,1);

    Plane3D planeXY = Plane3D::GetXYPlane();
    Plane3D planeXZ = Plane3D::GetXZPlane();
    Plane3D planeYZ = Plane3D::GetYZPlane();
    
    Point3Cartesian pnt(5, 5, 5);

    REQUIRE(planeXY.DistToPoint(pnt) == 5.0);
    REQUIRE(planeXZ.DistToPoint(pnt) == 5.0);
    REQUIRE(planeYZ.DistToPoint(pnt) == 5.0);
}