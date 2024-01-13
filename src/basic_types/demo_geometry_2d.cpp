#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Geometry.h"
#include "base/Geometry2D.h"
#include "base/Geometry3D.h"
#endif

using namespace MML;

void Demo_Geometry_2D_basic()
{
    Point2Cartesian pnt1;
    Point2Cartesian pnt2(1.0, 2.0);
    
    Vector2Cartesian vec1;
    Vector2Cartesian vec2(1.0, 2.0);
    Vector2Cartesian vec3(pnt1, pnt2);

    Point2Polar pnt_polar_1;
    Point2Polar pnt_polar_2(1.0, 2.0);
    
    Vector2Polar vec_polar_1;
    Vector2Polar vec_polar_2(1.0, 2.0);

    Line2D line1(Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0));
    Line2D line2(Point2Cartesian(0.0, 0.0), Vector2Cartesian(0.0, 1.0));

    SegmentLine2D seg_line1(Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0));
    SegmentLine2D seg_line2(Point2Cartesian(0.0, 0.0), Vector2Cartesian(0.0, 1.0), 1.0);

    Triangle2D triangle1(Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0), Point2Cartesian(0.0, 1.0));

    Polygon2D polygon1({Point2Cartesian(0.0, 0.0), Point2Cartesian(1.0, 0.0), Point2Cartesian(0.0, 1.0)});
}

void Demo_Geometry_2D_operations()
{
    Point2Cartesian pnt1(1.0, 2.0);
    Vector2Cartesian vec1(1.0, 2.0);

    // Point2Cartesian operations
    auto pnt2 = pnt1 + vec1;
    auto pnt3 = pnt1 - vec1;

    // Vector2Cartesian operations
    Vector2Cartesian vec2(3.0, -3.0);
    auto vec3 = vec1 + vec1;
    auto vec4 = vec1 - vec1;
    auto vec5 = 3.0 * vec1;
    auto vec6 = vec1 * 3.0;
    auto vec7 = vec1 / 3.0;
    auto norm = vec1.NormL2();
    auto scal_prod = vec1.ScalarProductCartesian(vec2);
}

void Demo_Geometry_2D()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         GEOMETRY 2D                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_Geometry_2D_basic();
    Demo_Geometry_2D_operations();
}