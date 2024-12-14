#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"

#include "algorithms/SurfaceIntegration.h"
#endif

#include "../test_data/Fields.h"

using namespace MML;

void Docs_Demo_Integration_surface()
{
  std::cout << std::endl;
  std::cout << "***********************************************************************" << std::endl;
  std::cout << "****                  SURFACE INTEGRATION                          ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  Cube3D cube(7);
  static VectorFunction<3> field([](const VectorN<Real, 3>& x_cart) { return Fields::InverseRadialPotentialForceFieldCart(1 / (4 * Constants::PI), x_cart); });

  Real integral = SurfaceIntegration::SurfaceIntegral(field, cube);

  std::cout << "Surface integral of the vector field over the cube: " << integral << std::endl;

}
