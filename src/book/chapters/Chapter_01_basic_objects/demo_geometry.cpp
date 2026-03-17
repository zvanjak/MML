///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_geometry.cpp
///  Description: Geometry overview header for Chapter 01
///               Acts as introduction before 2D and 3D geometry demos
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include <iostream>

#include "mml/base/Geometry/Geometry.h"
#endif

using namespace MML;

void Demo_Geometry()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                          GEOMETRY                             ****\n";
    std::cout << "***********************************************************************\n";

    // This is just a header - the actual demos are in:
    //   Demo_Geometry_2D() - Points, lines, shapes in 2D
    //   Demo_Geometry_3D() - Points, lines, planes in 3D
}
