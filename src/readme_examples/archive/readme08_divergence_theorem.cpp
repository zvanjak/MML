///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme08_divergence_theorem.cpp                                     ///
///  Description: README Gauss's Divergence Theorem verification demo                 ///
///               Demonstrates volume vs surface integration                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "base/Function.h"
#include "base/Geometry/Geometry3D.h"
#include "base/Geometry/Geometry3DBodies.h"

#include "core/Integration.h"
#include "core/Integration/SurfaceIntegration.h"
#include "core/FieldOperations.h"
#endif

#include <iomanip>

using namespace MML;

void Readme_DivergenceTheorem()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            README - Gauss's Divergence Theorem                ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Verify ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS over a unit cube
    // Define vector field F(x,y,z) = (sin(xy), cos(yz), exp(xz))
    VectorFunction<3> F([](const VectorN<Real, 3>& p) {
        return VectorN<Real, 3>{ sin(p[0]*p[1]), cos(p[1]*p[2]), exp(p[0]*p[2]) };
    });

    // Compute divergence NUMERICALLY - no analytical formula needed!
    ScalarFunctionFromStdFunc<3> divF([&F](const VectorN<Real, 3>& p) {
        return VectorFieldOperations::DivCart<3>(F, p);
    });

    // Integration limits for unit cube [0,1]³
    auto y_lo = [](Real) { return 0.0; };  auto y_hi = [](Real) { return 1.0; };
    auto z_lo = [](Real,Real) { return 0.0; };  auto z_hi = [](Real,Real) { return 1.0; };

    // Volume integral with different methods
    Real volGauss = Integrate3D(divF, GAUSS10, 0.0, 1.0, y_lo, y_hi, z_lo, z_hi).value;
    Real volTrap  = Integrate3D(divF, TRAP,    0.0, 1.0, y_lo, y_hi, z_lo, z_hi).value;

    // Surface integral (flux) through all 6 faces of the cube
    Cube3D unitCube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
    Real flux = SurfaceIntegration::SurfaceIntegral(F, unitCube, 1e-8);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Vector field F(x,y,z) = (sin(xy), cos(yz), exp(xz))" << std::endl;
    std::cout << std::endl;
    std::cout << "Volume integral (GAUSS10): " << volGauss << std::endl;
    std::cout << "Volume integral (TRAP):    " << volTrap << std::endl;
    std::cout << "Surface integral (flux):   " << flux << std::endl;
    std::cout << std::endl;
    std::cout << "GAUSS10 vs Surface error:  " << std::scientific << std::abs(volGauss - flux) << std::endl;

/* Expected OUTPUT:
    Vector field F(x,y,z) = (sin(xy), cos(yz), exp(xz))
    
    Volume integral (GAUSS10): 1.01945051
    Volume integral (TRAP):    1.01946656
    Surface integral (flux):   1.01944989
    
    GAUSS10 vs Surface error:  6.14417000e-07  (Divergence theorem verified!)
*/
}
