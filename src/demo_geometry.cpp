#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/Geometry.h"
#include "basic_types/Geometry2D.h"
#include "basic_types/Geometry3D.h"
#endif

using namespace MML;

void Demo_Geometry_2D()
{
    MML::Point2Cartesian pnt1;
    MML::Point2Cartesian pnt2(1.0, 2.0);
    
    MML::Vector2Cartesian vec1;
    MML::Vector2Cartesian vec2(1.0, 2.0);

    MML::Point2Polar pnt_polar_1;
    MML::Point2Polar pnt_polar_2(1.0, 2.0);
    
    MML::Vector2Polar vec_polar_1;
    MML::Vector2Polar vec_polar_2(1.0, 2.0);

    MML::Line2D line1(MML::Point2Cartesian(0.0, 0.0), MML::Point2Cartesian(1.0, 0.0));
    MML::Line2D line2(MML::Point2Cartesian(0.0, 0.0), MML::Vector2Cartesian(0.0, 1.0));

    MML::SegmentLine2D seg_line1(MML::Point2Cartesian(0.0, 0.0), MML::Point2Cartesian(1.0, 0.0));
    MML::SegmentLine2D seg_line2(MML::Point2Cartesian(0.0, 0.0), MML::Vector2Cartesian(0.0, 1.0), 1.0);

    MML::Triangle2D triangle1(MML::Point2Cartesian(0.0, 0.0), MML::Point2Cartesian(1.0, 0.0), MML::Point2Cartesian(0.0, 1.0));

    MML::Polygon2D polygon1({MML::Point2Cartesian(0.0, 0.0), MML::Point2Cartesian(1.0, 0.0), MML::Point2Cartesian(0.0, 1.0)});
}

void Demo_Geometry()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          GEOMETRY                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_Geometry_2D();

    // deklaracije - različite vrste točaka, 2D, 3D, Cartesian, Spherical
    MML::Point3Cartesian b{1.0, 2.0, 3.0};

    // pravci
    MML::Line3D  x_axis(MML::Point3Cartesian(0.0, 0.0, 0.0), MML::Vector3Cartesian{1.0, 0.0, 0.0});
    MML::Line3D  y_axis(MML::Point3Cartesian(0.0, 0.0, 0.0), MML::Vector3Cartesian{0.0, 1.0, 0.0});
    MML::Line3D  z_axis(MML::Point3Cartesian(0.0, 0.0, 0.0), MML::Vector3Cartesian{0.0, 0.0, 1.0});

    // duzine
    
    // poligoni

    // ravnine

    // tijela


}