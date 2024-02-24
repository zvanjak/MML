#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordTransf.h"
#include "core/Fields.h"
#include "core/FieldOperations.h"

#endif

using namespace MML;

using std::cout;
using std::endl;


void Demo_velocity_contravar_transf()
{
    cout << "***********************************************************************" << std::endl;
    cout << "****          CONTRAVARIANT TRANSFORMATION OF VELOCITY             ****" << std::endl;
    cout << "***********************************************************************" << std::endl;

    Vec3Cart v_cart{1.0, 1.0, 0.0};
    cout << "v_cart    : " << v_cart << std::endl;

    cout << "AT POINT:\n";
    Vector3Cartesian x1_cart{1.0, 1.0, 0.0};
    Vector3Spherical x1_spher{ CoordTransfCartToSpher.transf(x1_cart) };

    cout << "Cartesian : " << x1_cart << std::endl << "Spherical : " << x1_spher << std::endl;
     
    Vector3Spherical v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, x1_cart);
    cout << "contravar transf. to spher at cart.pnt : " << x1_cart << " = " << v_transf_to_spher << std::endl;

    Vector3Cartesian v_back_transf_to_cart = CoordTransfSpherToCart.transfVecContravariant(v_transf_to_spher, x1_spher);
    cout << "back transf. to cartesian at spher.pnt : " << x1_spher << " = " << v_back_transf_to_cart << std::endl;
}

void Demo_gradient_covariant_transf()
{
    cout << std::endl;
    cout << "***********************************************************************" << std::endl;
    cout << "****             COVARIANT TRANSFORMATION OF GRADIENT              ****" << std::endl;
    cout << "***********************************************************************" << std::endl;

    Vector3Cartesian p_cart{1.0, 1.0, 1.0};
    Vector3Spherical p_spher(CoordTransfSpherToCart.transfInverse(p_cart));
    
    cout << "Spherical: " << p_spher << std::endl;
    cout << "Cartesian: " << p_cart << std::endl;

    cout << "\n******   Working with SPHERICAL TO CARTESIAN transformation   *******\n";
    cout << "\n******   Potential in spherical coordinates  *******\n";
    
    ScalarFunction<3> fPotSpher(Fields::InverseRadialPotentialFieldSpher);    
    Vector3Spherical  grad_spher = ScalarFieldOperations::GradientSpher(fPotSpher, p_spher);

    cout << "Field at     : " << p_spher << " = " << Fields::InverseRadialPotentialFieldSpher(p_spher) << std::endl;
    cout << "Grad(sph) at : " << p_spher << " = " << grad_spher << std::endl;

    Vector3Cartesian grad_transf_to_cart  = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p_cart);
    Vector3Spherical back_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_transf_to_cart, p_spher);

    cout << "Grad.transf. (Cart.) at " << p_cart  << " = " << grad_transf_to_cart << std::endl;
    cout << "Back transf. (Spher) at " << p_spher << " = " << back_transf_to_spher << std::endl;

    cout << "\n******   Working with CARTESIAN TO SPHERICAL transformation   *******\n";
    cout << "\n******   Potential in cartesian coordinates:  *******\n";
    
    ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);    
    Vector3Cartesian  grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p_cart);

    cout << "Field at      : " << p_cart << " = " << Fields::InverseRadialPotentialFieldCart(p_cart) << std::endl;
    cout << "Grad.Cart. at : " << p_cart << " = " << grad_cart << std::endl;

    Vector3Spherical grad_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
    Vector3Cartesian back_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_transf_to_spher, p_cart);

    cout << "Grad.transf. (Spher) at " << p_spher << " = " << grad_transf_to_spher << std::endl;
    cout << "Back transf. (Cart.) at " << p_cart << " = " << back_transf_to_cart << std::endl;
}

void Demo_Covar_Contravar_transformations()
{
    Demo_velocity_contravar_transf();

    Demo_gradient_covariant_transf();
}