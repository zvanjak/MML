#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Geometry3D.h"
#endif

using namespace MML;

void Demo_Geometry_3D_basic()
{
    Point3Cartesian pnt1;
    Point3Cartesian pnt2(1.0, 2.0, 3.0);
    
    Vector3Cartesian vec1;
    Vector3Cartesian vec2(1.0, 2.0, 3.0);
    Vector3Cartesian vec3(pnt1, pnt2);

    Line3D line1(Point3Cartesian(0.0, 0.0, 1.0), Point3Cartesian(1.0, 0.0, 0.0));
    Line3D line2(Point3Cartesian(0.0, 0.0, 1.0), Vector3Cartesian(0.0, 1.0, 2.0));

    // pravci
    Line3D  x_axis(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian{1.0, 0.0, 0.0});
    Line3D  y_axis(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian{0.0, 1.0, 0.0});
    Line3D  z_axis(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian{0.0, 0.0, 1.0});

    SegmentLine3D seg_line1(Point3Cartesian(0.0, 0.0, 1.0), Point3Cartesian(1.0, 0.0, 0.0));
    SegmentLine3D seg_line2(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(0.0, 1.0, 0.0), 2.0);

    Plane3D plane1(Point3Cartesian(0.0, 0.0, 0.0), Vector3Cartesian(0.0, 0.0, 1.0));
    Plane3D plane2(Point3Cartesian(0.0, 0.0, 0.0), Point3Cartesian(0.0, 1.0, 0.0), Point3Cartesian(1.0, 0.0, 0.0));
    Plane3D plane3(1.0, 2.0, 3.0, 4.0);
    Plane3D plane4(1.0, 2.0, 3.0);

    // Triangle3D triangle1(Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0), Point2Cartesian(0.0, 1.0));

    // Polygon3D polygon1({Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0), Point2Cartesian(0.0, 1.0)});
}

void Demo_Geometry_3D_operations()
{
    Point3Cartesian  pnt1(1.0, 2.0, 1.0);
    Vector3Cartesian vec1(1.0, 2.0, 3.0);
    Vector3Cartesian vec2{3.0, -3.0, 1.0};

    // Point3Cartesian operations
    auto pnt2 = pnt1 + vec1;
    auto pnt3 = pnt1 - vec1;

    // Vector3Cartesian operations
    auto vec3 = vec1 + vec1;
    auto vec4 = vec1 - vec1;
    auto vec5 = 3.0 * vec1;
    auto vec6 = vec1 * 3.0;
    auto vec7 = vec1 / 3.0;

    auto norm = vec1.NormL2();
    auto scal_prod = ScalarProd(vec1, vec2);
    auto vec_prod  = VectorProd(vec1, vec2);

    bool b = vec1.IsParallelTo(vec2);
    bool c = vec1.IsPerpendicularTo(vec2);
}

void Demo_Plane3D_operations()
{
    Plane3D plane1(1.0, 3.0, -2.0, 4.0);
    Plane3D plane2(-1.0, 2.0, 3.0);
    Plane3D plane3(4.0, -1.0, 2.0);

    auto dist = plane1.DistToPoint(Point3Cartesian(1.0, 2.0, 3.0));
    bool ison = plane1.IsPointOnPlane(Point3Cartesian(1.0, 2.0, 3.0));
    auto proj = plane1.ProjectionToPlane(Point3Cartesian(1.0, 2.0, 3.0));
    auto isline = plane1.IsLineOnPlane(Line3D(Point3Cartesian(1.0, 2.0, 3.0), Vector3Cartesian(1.0, 2.0, 3.0)));
    double angle = plane1.AngleToLine(Line3D(Point3Cartesian(1.0, 2.0, 3.0), Vector3Cartesian(1.0, 2.0, 3.0)));
    
    Point3Cartesian outPnt;
    auto pnt_int = plane1.IntersectionWithLine(Line3D(Point3Cartesian(1.0, 2.0, 3.0), Vector3Cartesian(1.0, 2.0, 3.0)), outPnt);
    bool b2 = plane1.IsParallelToPlane(plane2);
    bool b3 = plane1.IsPerpendicularToPlane(plane2);
    auto angle2 = plane1.AngleToPlane(plane2);

    Line3D  outLine;
    bool b4 = plane1.IntersectionWithPlane(plane2, outLine);
}

void Demo_Geometry_3D()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         GEOMETRY 3D                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_Geometry_3D_basic();
    Demo_Geometry_3D_operations();

    Demo_Plane3D_operations();
}