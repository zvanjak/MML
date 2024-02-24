#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordTransf.h"
#include "core/FieldOperations.h"

#include "core/Fields.h"
#include "core/Curves.h"
#endif

using namespace MML;

void Investigating_vector_field_operations()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            INVESTIGATING - vector field operations            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    static const VectorN<Real, 3> x1{ 0.5, 0.0, 0.0 }, x2{ -0.5, 0.0, 0.0 };
    static const Real m1 = 0.5, m2 = 0.5;
    static const Real G = 1;

    ////////////////////////////////////////////////////////////////////////////////////////
    // two pot . func( cart, spher), wtih exact values
    static ScalarFunction<3> gr_potfield_cart_exact{ [](const VectorN<Real, 3> &x) 
    {
        return G * m1 / (x - x1).NormL2() + G * m2 / (x - x2).NormL2();
    } };
    static ScalarFunction<3> gr_potfield_spher_exact{ [](const VectorN<Real, 3> &xSpher) 
    {
        // input vector is Vector3Spherical position
        return gr_potfield_cart_exact(CoordTransfSpherToCart.transf(xSpher));
    } };

    ////////////////////////////////////////////////////////////////////////////////////////
    // defining force field of two masses as VectorFunction in Cartesian coord. - EXACT
    static VectorFunction<3> gr_forcefield_cart_exact{ [](const VectorN<Real, 3> &x) 
    {
        return -G * m1 * (x - x1) / std::pow((x - x1).NormL2(), 3) - G * m2 * (x - x2) / std::pow((x - x2).NormL2(), 3);
    } };
    // spher_exact ... nemoguce

    ////////////////////////////////////////////////////////////////////////////////////////
    // defining force field of two masses as VectorFunction - BASED ON NUM.GRADIENT CALCULATION OF POTENTIAL
    static VectorFunction<3> gr_forcefield_cart_num_from_grad{ [](const VectorN<Real, 3> &x) 
    {
        return (-1) * ScalarFieldOperations::GradientCart<3>(gr_potfield_cart_exact, x);
    } };    
    // defining force field of two masses as VectorFunction - BASED ON NUM.GRADIENT CALCULATION OF SPHERICAL POTENTIAL
    static VectorFunction<3> gr_forcefield_spher_num_from_grad{ [](const VectorN<Real, 3> &x) 
    {
        // input vector is Vector3Spherical position
        return (-1) * ScalarFieldOperations::GradientSpher(gr_potfield_spher_exact, x);
    } };    

    ////////////////////////////////////////////////////////////////////////////////////////
    // AKO NEMAM POTENCIJAL U SFERNIM KOORDINATAMA! a treba mi sferni gradijent
    static VectorFunction<3> grforcefield_spher_from_cart_pot{ [](const VectorN<Real, 3> &x)
    {
        // ulazna pozicija - spherical, 
        // racunamo numerical Cartesian gradijent, za koji nam treba Cartesian position
        Vector3Cartesian x_cart = CoordTransfSpherToCart.transf(x);
        auto force_cart = (-1) * ScalarFieldOperations::GradientCart<3>(gr_potfield_cart_exact, x_cart);

        // i onda ga KOVARIJANTNO TRANSFORMIRAMO u spher vektor gradijenta
        VectorN<Real, 3> ret = CoordTransfCartToSpher.transfVecCovariant(force_cart, x);

        // vraceni vektor je sada u force field spher koordinatama
        return ret;
    } };    



    // auto const &pot_cart_exact = gr_potfield_cart_exact;
    // auto const &pot_spher_exact = gr_potfield_spher_exact;
    // auto const &force_field_cart_from_grad = gr_forcefield_cart_exact;
    // auto const &force_field_spher_from_grad = gr_forcefield_spher_num_from_grad;
    // auto const &force_field_spher_from_cart_pot = grforcefield_spher_from_cart_pot;

/*  INTRO
**Fields and field operations - grad, div, curl, Laplacian**

Using as example inverse radial field, with its potential and force field, demonstrate calculation of field gradient, divergence, curl and Laplacian
Calculations are performed in Cartesian and spherical coordinates, along circle in XZ-plane, and covariant vector transformation is also demonstrated
*/
    // creating scalar potential and vector force field, from predefined functions
    static ScalarFunction<3> fPotCart([](const VectorN<Real, 3> &x_cart) -> Real { return Fields::InverseRadialPotentialFieldCart(x_cart); });
    static ScalarFunction<3> fPotSpher([](const VectorN<Real, 3> &x_spher) -> Real { return Fields::InverseRadialPotentialFieldSpher(x_spher); });
    static VectorFunction<3> fForceFieldCartExact([](const VectorN<Real, 3> &x_cart)  { return Fields::InverseRadialPotentialForceFieldCart(x_cart); });
    static VectorFunction<3> fForceFieldSpherExact([](const VectorN<Real, 3> &x_spher)  { return Fields::InverseRadialPotentialForceFieldSph(x_spher); });

    // if we have only potential, we can numerical calculate force field from it
    static VectorFunction<3> fForceFieldCart_numgrad{ [](const VectorN<Real, 3> &x_cart)  { return (-1) * ScalarFieldOperations::GradientCart<3>(fPotCart, x_cart); } };  
    static VectorFunction<3> fForceFieldSpher_numgrad([](const VectorN<Real, 3> &x_spher) { return -1 * ScalarFieldOperations::GradientSpher(fPotSpher,x_spher); });

    // if we have potential in one coord. system, and we need force field in another, we can calculate it 
    // by first calculating force field in the same coord. system as potential, and then transforming it covariantly
    VectorFunction<3> fForceFieldCart_from_spher_pot{ [](const VectorN<Real, 3> &x_cart)
    {
        auto force_spher = (-1) * ScalarFieldOperations::GradientSpher(fPotSpher, CoordTransfCartToSpher.transf(x_cart));

        VectorN<Real, 3> ret = CoordTransfSpherToCart.transfVecCovariant(force_spher, x_cart);

        return ret;        // vraceni vektor je sada u force field Cartesian koordinatama
    } };
    VectorFunction<3> fForceFieldSpher_from_cart_pot{ [](const VectorN<Real, 3> &x_spher)
    {
        // ulazna pozicija - spherical, racunamo numerical Cartesian gradijent, za koji nam treba Cartesian position
        Vector3Cartesian x_cart = CoordTransfSpherToCart.transf(x_spher);
        auto force_cart = (-1) * ScalarFieldOperations::GradientCart<3>(fPotCart, x_cart);

        // i onda ga KOVARIJANTNO TRANSFORMIRAMO u spher vektor gradijenta
        VectorN<Real, 3> ret = CoordTransfCartToSpher.transfVecCovariant(force_cart, x_spher);

        return ret;        // vraceni vektor je sada u force field spher koordinatama
    } };

    auto const &pot_cart_exact = fPotCart;
    auto const &pot_spher_exact = fPotSpher;
    
    auto const &force_field_cart_exact = fForceFieldCartExact;
    auto const &force_field_spher_exact = fForceFieldSpherExact;
    
    auto const &force_field_cart_from_grad = fForceFieldCart_numgrad;
    auto const &force_field_spher_from_grad = fForceFieldSpher_numgrad;
    
    auto const &force_field_cart_from_spher_pot = fForceFieldCart_from_spher_pot;
    auto const &force_field_spher_from_cart_pot = fForceFieldSpher_from_cart_pot;


    // calculating potential and force around circle
    Curves3D::Circle3DXZ circle(10.0);

    // calculating field Cart. and Sph. gradient around circle (and demonstrating covariant vector transformation of gradient)
    std::cout << "         Position          Pot.(Cart)    Pot.(Spher)        Force exact (Cart. vector)              Force exact (Sph. vector)          Force num.grad (Cart.vector)        Force num.grad (Sph.vector)" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=1)
    {
        Vector3Cartesian   pos_cart = circle(t);
        Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos_cart);

        double pot_cart  = pot_cart_exact(pos_cart);
        double pot_spher = pot_spher_exact(pos_spher);

        Vector3Cartesian   force_cart_exact  = force_field_cart_exact(pos_cart);
        Vector3Spherical   force_spher_exact = force_field_spher_exact(pos_spher);

        Vector3Cartesian   force_cart  = force_field_cart_from_grad(pos_cart);
        Vector3Spherical   force_spher = force_field_spher_from_grad(pos_spher);

        Vector3Spherical   force_spher_from_cart_pot = force_field_spher_from_cart_pot(pos_spher);

        std::cout << pos_cart.to_string(6,3) << "       " << pot_cart << "           "  << pot_spher
                << "      "  << force_cart_exact.to_string(10,4) 
                << "  "  << force_spher_exact.to_string(10,4) 
                << "  "  << force_cart.to_string(10,4) 
                << "  "  << force_spher.to_string(10,4) 
                << std::endl;
    }

    // calculating Divergence of four defined vector fields around circle
    std::cout << "            Position               Cart.Div of               Cart.Div. of         Sph.Div of               Sph.Div of    " << std::endl;
    std::cout << "                                  exact Cart.field          num.grad.field     exact Spher. field      num.grad.sph.field" << std::endl;
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
    // calculating vector field curl
    std::cout << "            Position                    Cartesian curl                        Spherical curl" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=1)
    {
        Vector3Cartesian   pos = circle(t);
        Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos);

        Vector3Cartesian   curl_cart  = VectorFieldOperations::CurlCart(force_field_cart_exact, pos);
        Vector3Spherical   curl_spher = VectorFieldOperations::CurlSpher(force_field_spher_exact, pos_spher);

        std::cout << pos.to_string(7,3) << " = " << curl_cart.to_string(10,4) 
                                        << "  "  << curl_spher.to_string(10,4) 
                                        << std::endl;
    }    
    // calculating scalar field Laplacian
    std::cout << "            Position              Cartesian Laplacian     Spherical Laplacian" << std::endl;
    for(double t=0.0; t<2*Constants::PI; t+=1)
    {
        Vector3Cartesian   pos = circle(t);
        Vector3Spherical   pos_spher = CoordTransfCartToSpher.transf(pos);

        double   lapl_cart  = ScalarFieldOperations::LaplacianCart(pot_cart_exact, pos);
        double   lapl_spher = ScalarFieldOperations::LaplacianSpher(pot_spher_exact, pos_spher);

        std::cout << pos.to_string(7,3) << "        " << std::setw(16) << lapl_cart
                                        << "        " << std::setw(16) << lapl_spher
                                        << std::endl;
    }       
}

void Demo_gradient()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           GRADIENT                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    ScalarFunction<3> fPotCart([](const VectorN<Real, 3> &x) -> Real { return Fields::InverseRadialPotentialFieldCart(x); });
    ScalarFunction<3> fPotSpher([](const VectorN<Real, 3> &x) -> Real { return Fields::InverseRadialPotentialFieldSpher(x); });
    ScalarFunction<3> fPotCyl([](const VectorN<Real, 3> &x) -> Real { return Fields::InverseRadialPotentialFieldCyl(x); });

    // calculating field gradient around circle
    Curves3D::Circle3DXZ circle(1.0);
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

    ScalarFunction<3> fPotCart([](const VectorN<Real, 3> &x) -> Real { return Fields::InverseRadialPotentialFieldCart(x); });
    ScalarFunction<3> fPotSpher([](const VectorN<Real, 3> &x) -> Real { return Fields::InverseRadialPotentialFieldSpher(x); });
    ScalarFunction<3> fPotCyl([](const VectorN<Real, 3> &x) -> Real { return Fields::InverseRadialPotentialFieldCyl(x); });

    Curves3D::Circle3DXZ circle(1.0);
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
    ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);    

    return ScalarFieldOperations::GradientCart<3>(fPotCart, x);
}

VectorN<Real, 3> GradientOfSphericalPotential(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotSpher(Fields::InverseRadialPotentialFieldSpher);    
    Vector3Spherical pos = x;

    return ScalarFieldOperations::GradientSpher(fPotSpher, pos);
}

VectorN<Real, 3> GradientOfCylindricalPotential(const VectorN<Real, 3> &x )
{
    ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCyl);    
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

    Curves3D::Circle3DXZ circle(1.0);
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

    Curves3D::Circle3DXZ circle(1.0);
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