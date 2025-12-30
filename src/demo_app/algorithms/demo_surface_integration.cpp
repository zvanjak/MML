#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "base/Geometry3DBodies.h"

#include "core/Integration/SurfaceIntegration.h"
#include "core/Fields.h"
#endif


using namespace std;
using namespace MML;

void Demo_Surface_integration()
{
    std::cout << endl;
    std::cout << "***********************************************************************" << endl;
    std::cout << "****                  SURFACE INTEGRATION                          ****" << endl;
    std::cout << "***********************************************************************" << endl;

    Cube3D cube(7);
    static VectorFunction<3> field([](const VectorN<Real, 3> &x_cart)   { return Fields::InverseRadialPotentialForceFieldCart(1/(4*Constants::PI), x_cart); });

    Real integral = SurfaceIntegration::SurfaceIntegral(field, cube);
    
    std::cout << "Surface integral of the vector field over the cube: " << integral << std::endl;
}
 