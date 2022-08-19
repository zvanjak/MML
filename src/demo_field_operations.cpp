#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/CoordTransf.h"
#include "algorithms/FieldOperations.h"

#endif

double PotentialCart(const MML::VectorN<3> &x )
{
    double r = x.NormL2();

    return 10.0 / r;
}

MML::VectorN<3> PotentialCartesianGradient(const MML::VectorN<3> &x )
{
    MML::ScalarFunctionFromFuncPtr<3> fPotCart(PotentialCart);    

    return MML::FieldOperations<3>::GradientCart(fPotCart, x);
}

double PotentialSpher(const MML::VectorN<3> &x )
{
    return 10.0 / x[0];
}

MML::VectorN<3> PotentialSphericalGradient(const MML::VectorN<3> &x )
{
    MML::ScalarFunctionFromFuncPtr<3> fPotSpher(PotentialSpher);    
    MML::Vector3Spherical pos = x;

    return MML::FieldOperations<3>::GradientSpher(fPotSpher, pos);
}

MML::VectorN<3> SimpleVectorFunc(const MML::VectorN<3> &x )
{
    MML::VectorN<3> ret;

    ret[0] = x[0]*x[0] * x[2]*x[2];
    ret[1] = -2 * x[1]*x[1] * x[2]*x[2];
    ret[2] = x[0] * x[1]*x[1] * x[2];

    return ret;
}

void Demo_Fields_old()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                       FIELD OPERATIONS                        ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::CoordTransfCartesianToSpherical transfCartToSpher;
    MML::CoordTransfSphericalToCartesian transfSpherToCart;

    MML::Vector3Cartesian p_cart{1.0, 1.0, 1.0};
    MML::Vector3Spherical p_spher(transfCartToSpher.transf(p_cart));
    std::cout << "Cartesian : " << p_cart << std::endl << "Spherical: " << p_spher << std::endl;


    std::cout << "\n******  Working with potential in cartesian coordinates:  *******\n";
    MML::ScalarFunctionFromFuncPtr<3> fPotCart(PotentialCart);    
    std::cout << "Field at " << p_cart << " = " << PotentialCart(p_cart) << std::endl;

    auto grad_cart = MML::FieldOperations<3>::GradientCart(fPotCart, p_cart);
    std::cout << "Gradient (Cartesian) at " << p_cart << " = " << grad_cart << std::endl;

    auto grad_transf_to_spher = transfSpherToCart.covariantTransf(grad_cart, p_spher);
    std::cout << "Transf.grad. (spher) at " << p_spher << " = " << grad_transf_to_spher << std::endl;


    std::cout << "\n******   Working with potential in spherical coordinates  *******\n";
    MML::ScalarFunctionFromFuncPtr<3> fPotSpher(PotentialSpher);    
    std::cout << "Field at  : " << p_spher << " = " << PotentialSpher(p_spher) << std::endl;

    auto grad_spher = MML::FieldOperations<3>::GradientSpher(fPotSpher, p_spher);
    std::cout << "Gradient (spherical) at " << p_spher << " = " << grad_spher << std::endl;

    auto grad_transf_to_cart = transfCartToSpher.covariantTransf(grad_spher, p_cart);
    std::cout << "Transf.grad. (Cart.) at " << p_cart << " = " << grad_transf_to_cart << std::endl;


    std::cout << "\n******  Transforming contravariant vector of velocity *******\n";
    MML::Vector3Cartesian v_cart{1.0, 0.0, 0.0};
    std::cout << "v_cart    : " << v_cart << std::endl;

    std::cout << "AT POINT:\n";
    MML::Vector3Cartesian x1_cart{-1.0, 0.0, 0.0};
    MML::Vector3Spherical x1_spher{ transfCartToSpher.transf(x1_cart) };

    std::cout << "Cartesian : " << x1_cart << std::endl << "Spherical : " << x1_spher << std::endl;
     
    auto v_transf_to_spher = transfSpherToCart.contravariantTransf(v_cart, x1_cart);
    std::cout << "contravar transf. to spher at cart.pnt : " << x1_cart << " = " << v_transf_to_spher << std::endl;

    auto v_back_transf_to_cart = transfCartToSpher.contravariantTransf(v_transf_to_spher, x1_spher);
    std::cout << "back transf. to cartesian at spher.pnt : " << x1_spher << " = " << v_back_transf_to_cart << std::endl;


    std::cout << "\n******  Divergence  *******\n";
    MML::VectorFunctionFromFuncPtr<3> fSimpleVecFunc(SimpleVectorFunc);

    MML::Vector3Cartesian y1_cart{1.0, -1.0, 1.0};
    auto div_simple = MML::FieldOperations<3>::DivCart(fSimpleVecFunc, y1_cart);
    std::cout << "Div simple vec. field at " << y1_cart << " = " << div_simple << std::endl;

    MML::VectorFunctionFromFuncPtr<3> fCartGrad(PotentialCartesianGradient);

    MML::Vector3Cartesian y2_cart{0.2, -0.2, 0.2};
    auto div_grad = MML::FieldOperations<3>::DivCart(fCartGrad, y2_cart);
    std::cout << "Div of cartesian gradient field at " << y2_cart << " = " << div_grad << std::endl;

    MML::VectorFunctionFromFuncPtr<3> fSpherGrad(PotentialSphericalGradient);

    MML::Vector3Cartesian y3_cart{0.2, 0.2, 0.2};
    auto div_grad_spher = MML::FieldOperations<3>::DivSpher(fSpherGrad, y3_cart);
    std::cout << "Div of spherical gradient field at " << y3_cart << " = " << div_grad_spher << std::endl;
}

void Demo_Fields_Test()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****             COVARIANT TRANSFORMATION OF GRADIENT              ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::CoordTransfSphericalToCartesian transfSpherToCart;

    MML::Vector3Cartesian p_cart{1.0, 1.0, 1.0};
    MML::Vector3Spherical p_spher(transfSpherToCart.transfInverse(p_cart));
    std::cout << "Spherical: " << p_spher << std::endl << "Cartesian: " << p_cart << std::endl;

    std::cout << "\n******   Working with SPHERICAL TO CARTESIAN transformation   *******\n";
    std::cout << "\n******   Potential in spherical coordinates  *******\n";
    MML::ScalarFunctionFromFuncPtr<3> fPotSpher(PotentialSpher);    
    std::cout << "Field at  : " << p_spher << " = " << PotentialSpher(p_spher) << std::endl;

    MML::Vector3Spherical grad_spher = MML::FieldOperations<3>::GradientSpher(fPotSpher, p_spher);
    std::cout << "Gradient (spherical) at " << p_spher << " = " << grad_spher << std::endl;

    MML::Vector3Cartesian grad_transf_to_cart = transfSpherToCart.covariantTransf(grad_spher, p_cart);
    std::cout << "Grad.transf. (Cart.) at " << p_cart << " = " << grad_transf_to_cart << std::endl;

    MML::Vector3Spherical back_transf_to_spher = transfSpherToCart.covariantTransf(grad_transf_to_cart, p_spher);
    std::cout << "Back transf. (Spher) at " << p_spher << " = " << back_transf_to_spher << std::endl;

    std::cout << "\n******   Working with  CARTESIAN TO SPHERICAL transformation   *******\n";
    MML::CoordTransfCartesianToSpherical transfCartToSpher;

    std::cout << "\n******  Potential in cartesian coordinates:  *******\n";
    MML::ScalarFunctionFromFuncPtr<3> fPotCart(PotentialCart);    
    std::cout << "Field at " << p_cart << " = " << PotentialCart(p_cart) << std::endl;

    MML::Vector3Cartesian grad_cart = MML::FieldOperations<3>::GradientCart(fPotCart, p_cart);
    std::cout << "Gradient (Cartesian) at " << p_cart << " = " << grad_cart << std::endl;

    MML::Vector3Spherical grad_transf_to_spher = transfCartToSpher.covariantTransf(grad_cart, p_spher);
    std::cout << "Grad.transf. (spher) at " << p_spher << " = " << grad_transf_to_spher << std::endl;

    MML::Vector3Cartesian back_transf_to_cart = transfCartToSpher.covariantTransf(grad_transf_to_spher, p_cart);
    std::cout << "Grad.transf. (spher) at " << p_cart << " = " << back_transf_to_cart << std::endl;

    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****          CONTRAVARIANT TRANSFORMATION OF VELOCITY             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::Vector3Cartesian v_cart{1.0, 0.0, 0.0};
    std::cout << "v_cart    : " << v_cart << std::endl;

    std::cout << "AT POINT:\n";
    MML::Vector3Cartesian x1_cart{-1.0, 0.0, 0.0};
    MML::Vector3Spherical x1_spher{ transfSpherToCart.transfInverse(x1_cart) };

    std::cout << "Cartesian : " << x1_cart << std::endl << "Spherical : " << x1_spher << std::endl;
     
    MML::Vector3Spherical v_transf_to_spher = transfSpherToCart.contravariantTransf(v_cart, x1_cart);
    std::cout << "contravar transf. to spher at cart.pnt : " << x1_cart << " = " << v_transf_to_spher << std::endl;

    MML::Vector3Cartesian v_back_transf_to_cart = transfSpherToCart.contravariantTransf(v_transf_to_spher, x1_spher);
    std::cout << "back transf. to cartesian at spher.pnt : " << x1_spher << " = " << v_back_transf_to_cart << std::endl;


    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          DIVERGENCE                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::VectorFunctionFromFuncPtr<3> fSimpleVecFunc(SimpleVectorFunc);

    MML::Vector3Cartesian y1_cart{1.0, -1.0, 1.0};
    auto div_simple = MML::FieldOperations<3>::DivCart(fSimpleVecFunc, y1_cart);
    std::cout << "Div simple vec. field at " << y1_cart << " = " << div_simple << std::endl;

    MML::VectorFunctionFromFuncPtr<3> fCartGrad(PotentialCartesianGradient);

    MML::Vector3Cartesian y2_cart{0.2, -0.2, 0.2};
    auto div_grad = MML::FieldOperations<3>::DivCart(fCartGrad, y2_cart);
    std::cout << "Div of cartesian gradient field at " << y2_cart << " = " << div_grad << std::endl;

    MML::VectorFunctionFromFuncPtr<3> fSpherGrad(PotentialSphericalGradient);

    MML::Vector3Spherical y2_spher{ transfSpherToCart.transfInverse(y2_cart) };

    auto div_grad_spher = MML::FieldOperations<3>::DivSpher(fSpherGrad, y2_spher);
    std::cout << "Div of spherical gradient field at " << y2_spher << " = " << div_grad_spher << std::endl;

}

void Demo_Field_operations()
{
    Demo_Fields_old();

    Demo_Fields_Test();
}