#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/Geometry.h"
#endif

void Demo_Geometry()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          GEOMETRY                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::Point3Cartesian a({1.0, 2.0, 3.0});

    MML::Line3D  line1(a, MML::Vector3Cartesian({1.0, 0.0, 0.0}));
}