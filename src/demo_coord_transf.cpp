#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/CoordTransf.h"
#endif


void Demo_CoordTransf()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         COORD TRANSF                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::CoordTransfSphericalToCartesian transf;

    MML::VectorN<Real, 3> p1{1.0, 1.0, 1.0};
    auto p1Spher = transf.transfInverse(p1);

    std::cout << "Cartesian: " << p1 << std::endl << "Spherical: " << p1Spher << std::endl;

    MML::Vector3Cartesian p2{1.0, 1.0, 1.0};
    auto p2Spher = transf.transfInverse(p2);

    std::cout << "Cartesian: " << p2 << std::endl << "Spherical: " << p2Spher << std::endl;

}