#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "basic_types/CoordTransf.h"
#endif

using namespace MML;

// TODO - demo unit vektora
void Demo_CoordTransf_Spherical()
{
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "Direct spherical transformation\n";

    Vector3Cartesian p1{5.0, 2.0, -3.0};
    auto p1Spher      = CoordTransfCartToSpher.transf(p1);
    auto p1BackTransf = CoordTransfSpherToCart.transf(p1Spher);

    std::cout << "Cartesian   : " << p1 << std::endl;
    std::cout << "Spherical   : " << p1Spher << std::endl;
    std::cout << "Back transf : " << p1BackTransf << std::endl;

    std::cout << "Inverse spherical transformation\n";
    
    Vector3Cartesian p2{5.0, 2.0, -3.0};
    auto p2Spher      = CoordTransfSpherToCart.transfInverse(p2);
    auto p2BackTransf = CoordTransfCartToSpher.transfInverse(p2Spher);

    std::cout << "Cartesian   : " << p2 << std::endl;
    std::cout << "Spherical   : " << p2Spher << std::endl;
    std::cout << "Back transf : " << p2BackTransf << std::endl;
}

void Demo_CoordTransf_Cylindrical()
{
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "Direct cylindrical transformation\n";

    Vector3Cartesian p1{5.0, 2.0, -3.0};
    auto p1Cyl        = CoordTransfCartToCyl.transf(p1);
    auto p1BackTransf = CoordTransfCylToCart.transf(p1Cyl);

    std::cout << "Cartesian   : " << p1 << std::endl;
    std::cout << "Cylindrical : " << p1Cyl << std::endl;
    std::cout << "Back transf : " << p1BackTransf << std::endl;

    std::cout << "Inverse cylindrical transformation\n";
    
    Vector3Cartesian p2{5.0, 2.0, -3.0};
    auto p2Cyl        = CoordTransfCylToCart.transfInverse(p2);
    auto p2BackTransf = CoordTransfCartToCyl.transfInverse(p2Cyl);

    std::cout << "Cartesian   : " << p2 << std::endl;
    std::cout << "Spherical   : " << p2Cyl << std::endl;
    std::cout << "Back transf : " << p2BackTransf << std::endl;
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
    Demo_CoordTransf_Cylindrical();
}