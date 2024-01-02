#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordTransf.h"
#include "core/FieldOperations.h"

#include "basic_types/Fields.h"
#include "basic_types/Curves.h"
#endif

using namespace MML;


void Demo_gradient()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           GRADIENT                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    ScalarFunction<3> fPotCart([](const VectorN<Real, 3> &x) -> Real { return InverseRadialPotentialFieldCart(x); });
    ScalarFunction<3> fPotSpher([](const VectorN<Real, 3> &x) -> Real { return InverseRadialPotentialFieldSpher(x); });
    ScalarFunction<3> fPotCyl([](const VectorN<Real, 3> &x) -> Real { return InverseRadialPotentialFieldCyl(x); });

    // calculating field gradient around circle
    Curves::Circle3DXZ circle(1.0);
    std::cout << "            Position                   Cartesian gradient              Spherical gradient         Spher.grad.covar.transf. to Cart        Cylindrical gradient         Cyl.grad.covar.transf. to Cart" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=0.3)
    {
        Vector3Cartesian   pos = circle(t);
        Vector3Cartesian   grad_cart  = ScalarFieldOperations::GradientCart<3>(fPotCart, pos);
        Vector3Spherical   grad_spher = ScalarFieldOperations::GradientSpher(fPotSpher, CoordTransfCartToSpher.transf(pos));
        Vector3Cylindrical grad_cyl   = ScalarFieldOperations::GradientCyl(fPotCyl, CoordTransfCartToCyl.transf(pos));

        Vector3Cartesian spher_grad_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, pos);
        Vector3Cartesian cyl_grad_transf_to_cart   = CoordTransfCylToCart.transfVecCovariant(grad_cyl, pos);

        std::cout << pos.to_string(8,3) << " = " << grad_cart.to_string(8,4) 
                                        << "  "  << grad_spher.to_string(8,4) 
                                        << "  "  << spher_grad_transf_to_cart.to_string(9,4) 
                                        << "  "  << grad_cyl.to_string(8,4) 
                                        << "  "  << cyl_grad_transf_to_cart.to_string(10,4) 
                                        << std::endl;
    }
}

void Demo_Laplacian()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           LAPLACIAN                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    ScalarFunction<3> fPotCart([](const VectorN<Real, 3> &x) -> Real { return InverseRadialPotentialFieldCart(x); });
    ScalarFunction<3> fPotSpher([](const VectorN<Real, 3> &x) -> Real { return InverseRadialPotentialFieldSpher(x); });
    ScalarFunction<3> fPotCyl([](const VectorN<Real, 3> &x) -> Real { return InverseRadialPotentialFieldCyl(x); });

    Curves::Circle3DXZ circle(1.0);
    std::cout << "            Position                 Cart. laplacian         Spher. laplacian         Cylin. laplacian" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=0.3)
    {
        Vector3Cartesian pos = circle(t);
        Real lapl_cart  = ScalarFieldOperations::LaplacianCart(fPotCart, pos);
        Real lapl_spher = ScalarFieldOperations::LaplacianSpher(fPotSpher, CoordTransfCartToSpher.transf(pos));
        Real lapl_cyl   = ScalarFieldOperations::LaplacianCyl(fPotCyl, CoordTransfCartToCyl.transf(pos));

        std::cout << pos.to_string(8,3) << " = " << std::setw(20) << lapl_cart << "    " << std::setw(20) << lapl_spher << "    " << std::setw(20) << lapl_cyl << std::endl;
    }    
}


VectorN<Real, 3> GradientOfCartesianPotential(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotCart(InverseRadialPotentialFieldCart);    

    return ScalarFieldOperations::GradientCart<3>(fPotCart, x);
}

VectorN<Real, 3> GradientOfSphericalPotential(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotSpher(InverseRadialPotentialFieldSpher);    
    Vector3Spherical pos = x;

    return ScalarFieldOperations::GradientSpher(fPotSpher, pos);
}

VectorN<Real, 3> GradientOfCylindricalPotential(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotCart(InverseRadialPotentialFieldCyl);    
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

void Demo_divergence()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          DIVERGENCE                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    VectorFunction<3> fSimpleVecFunc(SimpleVectorFunc);

    Vector3Cartesian x1_cart{1.0, -1.0, 1.0};
    auto div_simple = VectorFieldOperations::DivCart<3>(fSimpleVecFunc, x1_cart);
    std::cout << "Divergence of simple vec. field at " << x1_cart << " = " << div_simple << std::endl;

    VectorFunction<3> fCartGrad(GradientOfCartesianPotential);
    VectorFunction<3> fSpherGrad(GradientOfSphericalPotential);
    VectorFunction<3> fCylGrad(GradientOfCylindricalPotential);

    Curves::Circle3DXZ circle(1.0);
    std::cout << "            Position                 Cart. divergence    Spher. divergence  Cylin. divergence" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=0.3)
    {
        Vector3Cartesian   pos_cart = circle(t);
        Vector3Spherical   pos_spher{ CoordTransfCartToSpher.transf(pos_cart) };
        Vector3Cylindrical pos_cyl{ CoordTransfCartToCyl.transf(pos_cart) };

        auto div_grad       = VectorFieldOperations::DivCart<3>(fCartGrad, pos_cart);
        auto div_grad_spher = VectorFieldOperations::DivSpher(fSpherGrad, pos_spher);
        auto div_grad_cyl   = VectorFieldOperations::DivCyl(fCylGrad, pos_cyl);

        std::cout << pos_cart.to_string(8,3) << " = " << std::setw(20) << div_grad << "    " << div_grad_spher << "    " << div_grad_cyl << std::endl;
    }
}

void Demo_curl()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                             CURL                              ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    VectorFunction<3> fCartGrad(GradientOfCartesianPotential);
    VectorFunction<3> fSpherGrad(GradientOfSphericalPotential);
    VectorFunction<3> fCylGrad(GradientOfCylindricalPotential);

    Curves::Circle3DXZ circle(1.0);
    std::cout << "            Position                                Cartesian curl                                       Spherical curl                                     Cylindrical curl" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=0.3)
    {
        Vector3Cartesian   pos_cart = circle(t);
        Vector3Spherical   pos_spher{ CoordTransfCartToSpher.transf(pos_cart) };
        Vector3Cylindrical pos_cyl{ CoordTransfCartToCyl.transf(pos_cart) };

        auto curl_cart  = VectorFieldOperations::CurlCart(fCartGrad, pos_cart);
        auto curl_spher = VectorFieldOperations::CurlSpher(fSpherGrad, pos_spher);
        auto curl_cyl   = VectorFieldOperations::CurlCyl(fCylGrad, pos_cyl);

        std::cout << pos_cart.to_string(8,3) << " = " << curl_cart.to_string(15,8) 
                                             << "  " << curl_spher.to_string(15,8) 
                                             << "  " << curl_cyl.to_string(15,8) 
                                             << std::endl;
    }    
}


void Demo_Field_operations()
{
    Demo_gradient();
    Demo_divergence();
    Demo_curl();
    Demo_Laplacian();
}