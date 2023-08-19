#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/CoordTransf.h"
#include "basic_types/Fields.h"
#include "algorithms/FieldOperations.h"

#endif

using namespace MML;


VectorN<Real, 3> PotentialCartesianGradient(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotCart(InverseRadialFieldFuncCart);    

    return ScalarFieldOperations::GradientCart<3>(fPotCart, x);
}

VectorN<Real, 3> PotentialSphericalGradient(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotSpher(InverseRadialFieldFuncSpher);    
    Vector3Spherical pos = x;

    return ScalarFieldOperations::GradientSpher(fPotSpher, pos);
}

VectorN<Real, 3> PotentialCylindricalGradient(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotCart(InverseRadialFieldFuncCyl);    
    Vector3Cylindrical pos {x};

    return ScalarFieldOperations::GradientCyl(fPotCart, pos);
}

VectorN<Real, 3> SimpleVectorFunc(const VectorN<Real, 3> &x )
{
    VectorN<Real, 3> ret;

    ret[0] = x[0]*x[0] * x[2]*x[2];
    ret[1] = -2 * x[1]*x[1] * x[2]*x[2];
    ret[2] = x[0] * x[1]*x[1] * x[2];

    return ret;
}

void Demo_gradient()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           GRADIENT                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
}

void Demo_divergence()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          DIVERGENCE                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    VectorFunction<3> fSimpleVecFunc(SimpleVectorFunc);

    Vector3Cartesian y1_cart{1.0, -1.0, 1.0};
    auto div_simple = VectorFieldOperations::DivCart<3>(fSimpleVecFunc, y1_cart);
    std::cout << "Divergence of simple vec. field at " << y1_cart << " = " << div_simple << std::endl;

    VectorFunction<3> fCartGrad(PotentialCartesianGradient);
    VectorFunction<3> fSpherGrad(PotentialSphericalGradient);
    VectorFunction<3> fCylGrad(PotentialCylindricalGradient);

    Vector3Cartesian y2_cart{0.2, -0.2, 0.2};
    Vector3Spherical y2_spher{ CoordTransfSpherToCart.transfInverse(y2_cart) };
    Vector3Cylindrical y2_cyl{ CoordTransfCartToCyl.transf(y2_cart) };

    auto div_grad       = VectorFieldOperations::DivCart<3>(fCartGrad, y2_cart);
    auto div_grad_spher = VectorFieldOperations::DivSpher(fSpherGrad, y2_spher);
    auto div_grad_cyl   = VectorFieldOperations::DivCyl(fCylGrad, y2_cyl);

    std::cout << "Div of cartesian gradient field at   : " << y2_cart << " = " << div_grad << std::endl;
    std::cout << "Div of spherical gradient field at   : " << y2_spher << " = " << div_grad_spher << std::endl;
    std::cout << "Div of cylindrical gradient field at : " << y2_cyl << " = " << div_grad_cyl << std::endl;       
}

void Demo_curl()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                             CURL                              ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
}

void Demo_Laplacian()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           LAPLACIAN                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
}

void Demo_Field_operations()
{
    Demo_gradient();
    Demo_divergence();
    Demo_curl();
    Demo_Laplacian();
}