#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Fields.h"
#include "core/Function.h"

#include "algorithms/SurfaceIntegration.h"
#endif

using namespace std;
using namespace MML;

void Demo_Surface_integration()
{
    std::cout << endl;
    std::cout << "***********************************************************************" << endl;
    std::cout << "****                  SURFACE INTEGRATION                          ****" << endl;
    std::cout << "***********************************************************************" << endl;

    Cube3D cube(10.0);
    static VectorFunction<3> field([](const VectorN<Real, 3> &x_cart)   { return InverseRadialPotentialForceFieldCart(x_cart); });

    Real integral = SurfaceIntegration::SurfaceIntegral(field, cube);
    
    std::cout << "Surface integral of the vector field over the cube: " << integral << std::endl;
}
 