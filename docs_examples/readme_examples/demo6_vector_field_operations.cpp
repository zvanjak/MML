#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordTransf.h"
#include "core/FieldOperations.h"

#include "core/Fields.h"
#include "core/Curves.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"
#include "../test_data/diff_eq_systems_test_bed.h"
#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

// TODO 0.7 - bolji primjer - skalarno polje potencijala, iz njega generirati vektorsko polje, i iz njega grad, div, curl, Laplacian
void Readme_vector_field_operations()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****               README - vector field operations                ****" << std::endl;
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
