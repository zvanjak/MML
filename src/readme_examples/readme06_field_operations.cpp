///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme06_field_operations.cpp                                       ///
///  Description: README Field Operations section demo                                ///
///               Demonstrates scalar and vector field operations                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/VectorN.h"
#include "base/Function.h"
#include "core/FieldOperations.h"
#endif

using namespace MML;

void Readme_FieldOperations()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            README - Field Operations                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Scalar field: gravitational potential φ(x,y,z) = -1/r
    ScalarFunction<3> potential([](const VectorN<Real, 3>& x) {
        return -1.0 / x.NormL2();
    });
    VectorN<Real, 3> pos{1.0, 2.0, 2.0};

    // Gradient ∇φ (gives force direction)
    auto grad = ScalarFieldOperations::GradientCart<3>(potential, pos);
    std::cout << "Gradient at (1,2,2): " << grad << std::endl;

    // Laplacian ∇²φ (zero outside mass for gravity!)
    Real laplacian = ScalarFieldOperations::LaplacianCart<3>(potential, pos);
    std::cout << "Laplacian at (1,2,2): " << laplacian << std::endl;

    // Vector field: rotating velocity field
    VectorFunction<3> velocity([](const VectorN<Real, 3>& x) -> VectorN<Real, 3> {
        return {x[1], -x[0], x[2]};
    });

    // Divergence ∇·v (compression/expansion rate)
    Real div = VectorFieldOperations::DivCart<3>(velocity, pos);
    std::cout << "Divergence at (1,2,2): " << div << std::endl;

    // Curl ∇×v (rotation/vorticity)
    auto curl = VectorFieldOperations::CurlCart(velocity, pos);
    std::cout << "Curl at (1,2,2): " << curl << std::endl;

/* Expected OUTPUT:
    Gradient at (1,2,2): [0.0370370370, 0.0740740741, 0.0740740741]  (proportional to r/r³)
    Laplacian at (1,2,2): 0.0 (point mass, outside source)
    Divergence at (1,2,2): 1.0 (∂z/∂z = 1)
    Curl at (1,2,2): [0, 0, -2] (constant rotation about z-axis)
*/
}
