#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordTransf.h"
#include "core/Fields.h"
#include "core/FieldOperations.h"
#include "core/Curves.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

// TODO 0.9 - odraditi line work integral
// 	- kako to odraditi u sfernim koord?
void Readme_vector_field_operations()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****               README - vector field operations                ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
/*  INTRO
**Fields and field operations - grad, div, curl, Laplacian**

Using as example inverse radial field, with its potential and force field, demonstrate calculation of field gradient, divergence, curl and Laplacian
Calculations are performed in Cartesian and spherical coordinates, along circle in XZ-plane, and covariant vector transformation is also demonstrated
*/
    // Setting up fields and creating scalar potential and vector force field from predefined functions
    static ScalarFunction<3> pot_cart_exact([](const VectorN<Real, 3> &x_cart) -> Real   { return -Fields::InverseRadialPotentialFieldCart(x_cart); });
    static ScalarFunction<3> pot_spher_exact([](const VectorN<Real, 3> &x_spher) -> Real { return -Fields::InverseRadialPotentialFieldSpher(x_spher); });
    static VectorFunction<3> force_field_cart_exact([](const VectorN<Real, 3> &x_cart)   { return Fields::InverseRadialPotentialForceFieldCart(x_cart); });
    static VectorFunction<3> force_field_spher_exact([](const VectorN<Real, 3> &x_spher) { return Fields::InverseRadialPotentialForceFieldSph(x_spher); });

    // if we have only potential, we can numerical calculate force field from it
    static VectorFunction<3> force_field_cart_from_grad{ [](const VectorN<Real, 3> &x_cart)  { return -ScalarFieldOperations::GradientCart<3>(pot_cart_exact, x_cart); } };  
    static VectorFunction<3> force_field_spher_from_grad([](const VectorN<Real, 3> &x_spher) { return -ScalarFieldOperations::GradientSpher(pot_spher_exact,x_spher); });

    // if we have potential in one coord. system, and we need force field in another, we can calculate it 
    // by first calculating force field in the same coord. system as potential, and then transforming it covariantly to desired coordinates
    VectorFunction<3> force_field_cart_from_spher_pot{ [](const VectorN<Real, 3> &x_cart)
    {
        auto force_spher     = -ScalarFieldOperations::GradientSpher(pot_spher_exact, CoordTransfCartToSpher.transf(x_cart));
        VectorN<Real, 3> ret = CoordTransfSpherToCart.transfVecCovariant(force_spher, x_cart);
        return ret;        // returned force field vector is now in Cartesian coordinates
    } };
    VectorFunction<3> force_field_spher_from_cart_pot{ [](const VectorN<Real, 3> &x_spher)
    {
        auto force_cart      = -ScalarFieldOperations::GradientCart<3>(pot_cart_exact, CoordTransfSpherToCart.transf(x_spher));
        VectorN<Real, 3> ret = CoordTransfCartToSpher.transfVecCovariant(force_cart, x_spher);
        return ret;        // returned force field vector is now in Spherical coordinates
    } };

    // calculating potential and force around circle
    Curves3D::Circle3DXZ circle(10.0);

    // calculating field Gradient in Cart. and Spherical coordinates around circle
    std::cout << "         Position          Pot.(Cart)    Pot.(Spher)        Force exact (Cart. vector)              Force exact (Sph. vector)              Force num.grad (Cart.vector)        Force num.grad (Sph.vector)" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=1)
    {
        Vector3Cartesian   pos_cart = circle(t);
        Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos_cart);

        double pot_cart  = pot_cart_exact(pos_cart);
        double pot_spher = pot_spher_exact(pos_spher);

        Vector3Cartesian   force_cart_exact    = force_field_cart_exact(pos_cart);
        Vector3Spherical   force_spher_exact   = force_field_spher_exact(pos_spher);
        Vector3Cartesian   force_cart_numgrad  = force_field_cart_from_grad(pos_cart);
        Vector3Spherical   force_spher_numgrad = force_field_spher_from_grad(pos_spher);

        std::cout << pos_cart.to_string(6,3) << "       " << pot_cart << "           "  << pot_spher << "    "
                << "  "  << force_cart_exact.to_string(10,4) 
                << "  "  << force_spher_exact.to_string(10,4) 
                << "  "  << force_cart_numgrad.to_string(10,4) 
                << "  "  << force_spher_numgrad.to_string(10,4) 
                << std::endl;
    }

    // calculating field Divergence of four defined vector fields around circle
    std::cout << "            Position               Cart.Div of               Cart.Div. of         Sph.Div of              Sph.Div of   " << std::endl;
    std::cout << "                                  exact Cart.field          num.grad.field     exact Spher. field    num.grad.sph.field" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=1.0)
    {
        Vector3Cartesian   pos_cart  = circle(t);
        Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos_cart);

        auto div_cart_exact  = VectorFieldOperations::DivCart<3>(force_field_cart_exact, pos_cart);
        auto div_spher_exact = VectorFieldOperations::DivSpher(force_field_spher_exact, pos_spher);
        
        auto div_grad_cart_numgrad  = VectorFieldOperations::DivCart<3>(force_field_cart_from_grad, pos_cart);
        auto div_grad_spher_numgrad = VectorFieldOperations::DivSpher(force_field_spher_from_grad, pos_spher);
        
        std::cout << pos_cart.to_string(7,3) << " " << std::setw(20) << div_cart_exact << "      " << std::setw(20) << div_grad_cart_numgrad << "      " 
                                             << div_spher_exact << "      " << div_grad_spher_numgrad << std::endl;
    }

    // calculating vector field Curl and field Laplacian
    std::cout << "            Position                    Cartesian curl                        Spherical curl                   Cartesian Laplacian     Spherical Laplacian" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=1)
    {
        Vector3Cartesian   pos = circle(t);
        Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos);

        Vector3Cartesian   curl_cart  = VectorFieldOperations::CurlCart(force_field_cart_exact, pos);
        Vector3Spherical   curl_spher = VectorFieldOperations::CurlSpher(force_field_spher_exact, pos_spher);

        double   lapl_cart  = ScalarFieldOperations::LaplacianCart(pot_cart_exact, pos);
        double   lapl_spher = ScalarFieldOperations::LaplacianSpher(pot_spher_exact, pos_spher);

        std::cout << pos.to_string(7,3) << " = " << curl_cart.to_string(10,4) 
                                        << "  "  << curl_spher.to_string(10,4)
                                        << "        " << std::setw(16) << lapl_cart
                                        << "        " << std::setw(16) << lapl_spher
                                        << std::endl;
    }    

/* OUTPUT
         Position          Pot.(Cart)    Pot.(Spher)        Force exact (Cart. vector)              Force exact (Sph. vector)              Force num.grad (Cart.vector)        Force num.grad (Sph.vector)
[    10,      0,      0]       -0.1           -0.1      [     -0.01,         -0,         -0]  [     -0.01,          0,          0]  [     -0.01,         -0,         -0]  [     -0.01,         -0,         -0]
[   5.4,      0,   8.41]       -0.1           -0.1      [ -0.005403,         -0,  -0.008415]  [     -0.01,          0,          0]  [ -0.005403,         -0,  -0.008415]  [     -0.01,         -0,         -0]
[ -4.16,      0,   9.09]       -0.1           -0.1      [  0.004161,         -0,  -0.009093]  [     -0.01,          0,          0]  [  0.004161,         -0,  -0.009093]  [     -0.01,         -0,         -0]
[  -9.9,      0,   1.41]       -0.1           -0.1      [    0.0099,         -0,  -0.001411]  [     -0.01,          0,          0]  [    0.0099,         -0,  -0.001411]  [     -0.01,         -0,         -0]
[ -6.54,      0,  -7.57]       -0.1           -0.1      [  0.006536,         -0,   0.007568]  [     -0.01,          0,          0]  [  0.006536,         -0,   0.007568]  [     -0.01,         -0,         -0]
[  2.84,      0,  -9.59]       -0.1           -0.1      [ -0.002837,         -0,   0.009589]  [     -0.01,          0,          0]  [ -0.002837,         -0,   0.009589]  [     -0.01,         -0,         -0]
[   9.6,      0,  -2.79]       -0.1           -0.1      [ -0.009602,         -0,   0.002794]  [     -0.01,          0,          0]  [ -0.009602,         -0,   0.002794]  [     -0.01,         -0,         -0]
            Position               Cart.Div of               Cart.Div. of         Sph.Div of              Sph.Div of
                                  exact Cart.field          num.grad.field     exact Spher. field    num.grad.sph.field
[     10,       0,       0]     -6.505213035e-18          -9.889788207e-12      -9.728329253e-16      -2.036693314e-12
[    5.4,       0,    8.41]      1.140363845e-15          -1.403275347e-11      -9.728329253e-16      -2.036693314e-12
[  -4.16,       0,    9.09]     -2.109857428e-15          -2.632835248e-11      -9.728329253e-16      -2.036693314e-12
[   -9.9,       0,    1.41]       1.86168355e-15           7.353272722e-13      -9.728329253e-16      -2.036693314e-12
[  -6.54,       0,   -7.57]     -7.254396736e-16          -1.703953005e-11      -9.728329253e-16      -2.036693314e-12
[   2.84,       0,   -9.59]     -3.102769777e-15          -1.202754346e-11      -9.728329253e-16      -2.036693314e-12
[    9.6,       0,   -2.79]      7.181755191e-16           5.078564305e-12      -9.728329253e-16      -2.036693314e-12
            Position                    Cartesian curl                        Spherical curl                   Cartesian Laplacian     Spherical Laplacian
[     10,       0,       0] = [         0,          0,          0]  [         0,          0,          0]        -2.865156029e-14         7.464116131e-13
[    5.4,       0,    8.41] = [         0,  2.408e-15,          0]  [         0,          0,          0]        -3.885927604e-12         7.464116131e-13
[  -4.16,       0,    9.09] = [         0,  1.806e-16,          0]  [         0,          0,          0]         5.927392041e-12         7.464116131e-13
[   -9.9,       0,    1.41] = [         0, -1.715e-15,          0]  [         0,          0,          0]         5.756142091e-13         7.464116131e-13
[  -6.54,       0,   -7.57] = [         0,  6.622e-16,          0]  [         0,          0,          0]        -3.369075635e-12         7.464116131e-13
[   2.84,       0,   -9.59] = [         0,  1.625e-15,          0]  [         0,          0,          0]        -2.460620032e-12         7.464116131e-13
[    9.6,       0,   -2.79] = [         0, -2.769e-15,          0]  [         0,          0,          0]         6.892316579e-13         7.464116131e-13
*/
}
