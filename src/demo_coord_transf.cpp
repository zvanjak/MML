#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/CoordTransf.h"
#endif

using namespace MML;

void Demo_CoordTransf_Spherical()
{
    MML::CoordTransfSphericalToCartesian transf;

    MML::VectorN<Real, 3> p1{1.0, 1.0, 1.0};
    auto p1Spher = transf.transfInverse(p1);

    std::cout << "Cartesian: " << p1 << std::endl << "Spherical: " << p1Spher << std::endl;

    MML::Vector3Cartesian p2{2.0, -1.0, 3.0};
    auto p2Spher = MML::CoordTransfSpherToCart.transfInverse(p2);

    std::cout << "Cartesian: " << p2 << std::endl << "Spherical: " << p2Spher << std::endl;

    auto p3Spher = MML::CoordTransfSpherToCart.transfInverse(MML::Vector3Cartesian(2.0, -1.0, 3.0));

    std::cout << "Cartesian: " << MML::Vector3Cartesian(2.0, -1.0, 3.0) << std::endl << "Spherical: " << p3Spher << std::endl;    
}

void Demo_CoordTransf_Rectilinear()
{

}

void Demo_CoordTransf()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         COORD TRANSF                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_CoordTransf_Spherical();
}