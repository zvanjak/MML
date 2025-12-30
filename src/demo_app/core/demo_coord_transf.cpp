#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordSystem.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransf3D.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"

#include "core/FieldOperations.h"
#include "core/Fields.h"
#endif


using namespace MML;

// Recilinear - orthogonal
// Rectilinear - oblique
// 2D polar - orthogonal!
// 2D eliptic - oblique
// 3D spherical
// 3D cylindrical
// general

// active vs passive transformation

// transformacija tocke
// covar i contravar bazni vektori u tocki
// kako dobiti unit vektore
// metricki tenzor u tocki
// transformacija vektora

void Demo_2DPolar()
{
	std::cout << "***************    OBLIQUE RECTILINEAR   ****************" << std::endl;

	// new basis vectors (keeping z coordinate out of it ;)
	Vector3Cartesian e1_base{ 1, 3, 0 }, e2_base{ 4, 0, 0 }, e3_base{ 0, 0, 1 };

	CoordTransfCartesianToOblique3D transf(e1_base, e2_base, e3_base);

	std::cout << "Is right-handed? - " << transf.IsRightHanded() << std::endl;

	std::cout << "New basis vectors:\n";
	std::cout << "e1 = "; e1_base.Print(std::cout, 10, 5) << std::endl;
	std::cout << "e2 = "; e2_base.Print(std::cout, 10, 5) << std::endl;
	std::cout << "e3 = "; e3_base.Print(std::cout, 10, 5) << std::endl;

	std::cout << "Calculated dual basis vectors:\n";
	std::cout << "e1 dual = "; transf.Dual(0).Print(std::cout, 10, 5) << std::endl;
	std::cout << "e2 dual = "; transf.Dual(1).Print(std::cout, 10, 5) << std::endl;
	std::cout << "e3 dual = "; transf.Dual(2).Print(std::cout, 10, 5) << std::endl;

	//std::cout << "\nTransf. matrix:\n" << transf.getAlpha() << std::endl;
	//std::cout << "Inverse transf. matrix:\n" << transf.getTransf() << std::endl;

	Vector3Cartesian vec_A{ 7, 2, 0 };
	std::cout << "Vector A (orig)      = "; vec_A.Print(std::cout, 10, 5) << std::endl;

	// std::cout << "alpha * A            = " << transf._alpha * vec_A << std::endl;
	// std::cout << "transf * (alpha * A) = " << transf._transf * (transf._alpha * vec_A) << std::endl;

	std::cout << "\nTransformed to new basis:\n";
	auto contravar_coef = transf.transf(vec_A);
	auto back_contravar = transf.transfInverse(contravar_coef);
	std::cout << "Contravar. coeff.   = "; contravar_coef.Print(std::cout, 10, 5) << std::endl;
	std::cout << "Back transf coeff   = "; back_contravar.Print(std::cout, 10, 5) << std::endl;

	std::cout << "\nCalculating contravariant and covariant components in new basis:\n";
	Vector3Cartesian x_dummy{ 1.0, 1.0, 1.0 };          // point of transf. application is irrelevant (it is linear)

	auto contravar_vec = transf.transfVecContravariant(vec_A, x_dummy);
	auto contravar_vec2 = transf.transfInverseVecContravariant(contravar_vec, x_dummy);
	std::cout << "Contravar. comp. of A = "; contravar_vec.Print(std::cout, 10, 5) << std::endl;
	std::cout << "Back transf.contravar = "; contravar_vec2.Print(std::cout, 10, 5) << std::endl;

	auto covar_vec = transf.transfVecCovariant(vec_A, x_dummy);
	auto covar_vec2 = transf.transfInverseVecCovariant(covar_vec, x_dummy);
	std::cout << "Covariant  comp. of A = "; covar_vec.Print(std::cout, 10, 5) << std::endl;
	std::cout << "Back transf.covariant = "; covar_vec2.Print(std::cout, 10, 5) << std::endl;
}

void Demo_CoordTransf_RectilinearOblique()
{
	std::cout << "\n****    RECTILINEAR OBLIQUE COORD SYSTEM   ****" << std::endl;

	//    Vector3Cartesian e1_base{1, 0, 0}, e2_base{0, 1, 0}, e3_base{0, 0.1, 1};      // ovo definira dot koord system
	Vector3Cartesian e1_base{ 1, 3, 0 }, e2_base{ 4, 0, 0 }, e3_base{ 0, 0, 1 };

	CoordTransfCartesianToOblique3D transf(e1_base, e2_base, e3_base);

	std::cout << "e1 = " << e1_base << std::endl;
	std::cout << "e2 = " << e2_base << std::endl;
	std::cout << "e3 = " << e3_base << std::endl;

	std::cout << "e1 dual = " << transf.Dual(0) << std::endl;
	std::cout << "e2 dual = " << transf.Dual(1) << std::endl;
	std::cout << "e3 dual = " << transf.Dual(2) << std::endl;

	Vector3Cartesian vec_A{ 7.0, 2.0, 0.0 };
	std::cout << "\nVector A: " << vec_A << std::endl;

	auto contravar_coef = transf.transf(vec_A);

	std::cout << "Contravar. coeff.: " << contravar_coef << std::endl;
	VectorN<Real, 3> contra_expanded = contravar_coef[0] * e1_base +
		contravar_coef[1] * e2_base +
		contravar_coef[2] * e3_base;
	std::cout << "Expanded to orig.: " << contra_expanded << std::endl;

	auto contravar_trans = transf.transfVecContravariant(vec_A, Vector3Cartesian(1, 0, -1));
	std::cout << "Contravar transf.: " << contravar_trans << std::endl;

	Vector3Cartesian covar_coef{ ScalarProduct(vec_A, Vector3Cartesian(e1_base)),
															 ScalarProduct(vec_A, Vector3Cartesian(e2_base)),
															 ScalarProduct(vec_A, Vector3Cartesian(e3_base)) };
	std::cout << "Covar. coeff.    : " << covar_coef << std::endl;

	VectorN<Real, 3> covar_expanded = covar_coef[0] * transf.Dual(0) +
		covar_coef[1] * transf.Dual(1) +
		covar_coef[2] * transf.Dual(2);
	std::cout << "Expanded to orig.: " << covar_expanded << std::endl;

	auto covar_trans = transf.transfVecCovariant(vec_A, Vector3Cartesian(1, 0, -1));
	std::cout << "Covar transf.: " << covar_trans << std::endl;

}

void Demo_CoordTransf_Jacobian()
{
	Vector3Cartesian p1{ 2.0, 2.0, 5 };
	auto   p1Spher = CoordTransfCartToSpher.transf(p1);

	MatrixNM<Real, 3, 3> jac = CoordTransfCartToSpher.jacobian(p1Spher);
}

void Demo_CoordTrans_RectilinearOrthogonal()
{
	// vectors of new base
	Vector3Cartesian b1{ 0, 1, 0 }, b2{ -1, 0, 0 }, b3{ 0, 0, 1 };

	std::cout << "e1 = " << b1 << std::endl;
	std::cout << "e2 = " << b2 << std::endl;
	std::cout << "e3 = " << b3 << std::endl;

	CoordTransf3DCartOrthogonal transf(b1, b2, b3);

	Vector3Cartesian origPnt{ 2, -3, 4 };
	auto newPnt = transf.transf(origPnt);
	auto backPnt = transf.transfInverse(newPnt);

	std::cout << "Original point: " << origPnt << std::endl;
	std::cout << "New point: " << newPnt << std::endl;
	std::cout << "Back transformed point: " << backPnt << std::endl;

	std::cout << "Is right-handed " << transf.IsRightHanded() << std::endl;
}

void Demo_CoordTransf_Spherical()
{
	std::cout << "***********************************************************************\n";
	std::cout << "****              SPHERICAL COORDINATE TRANSFORMATION              ****\n";
	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Point (position) transformation - direct \n";

	//	Vector3Cartesian p1{ 5.0, 2.0, -3.0 };
	Vector3Cartesian p1{ 2.0, 2.0, 1.0 };
	auto p1Spher = CoordTransfCartToSpher.transf(p1);
	auto p1BackTransf = CoordTransfSpherToCart.transf(p1Spher);

	std::cout << "Cartesian   : " << p1 << std::endl;
	std::cout << "Spherical   : " << p1Spher << std::endl;
	std::cout << "Back transf : " << p1BackTransf << std::endl;

	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Point (position) transformation - inverse \n";

	Vector3Cartesian p2{ 5.0, 2.0, -3.0 };
	auto p2Spher = CoordTransfSpherToCart.transfInverse(p2);
	auto p2BackTransf = CoordTransfCartToSpher.transfInverse(p2Spher);

	std::cout << "Cartesian   : " << p2 << std::endl;
	std::cout << "Spherical   : " << p2Spher << std::endl;
	std::cout << "Back transf : " << p2BackTransf << std::endl;

	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Unit-vectors:\n";

	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Vector transformation - contravariant\n";

	Vector3Cartesian v1{ 2.0, 2.0, 0.0 };
	std::cout << "Vector (contravariant)       : " << v1 << std::endl;
	auto v1Spher = CoordTransfCartToSpher.transfVecContravariant(v1, p1);
	std::cout << "contravar transf. to spher at cart.pnt : " << p1 << " = " << v1Spher << std::endl;

	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Vector transformation - covariant\n";

	auto v2Spher = CoordTransfCartToSpher.transfVecCovariant(v1, p1Spher);
	std::cout << "covariant transf. to spher at cart.pnt : " << p1 << " = " << v2Spher << std::endl;

	std::cout << "Point (cartesian): " << p1 << std::endl;
	std::cout << "Point (spherical): " << p1Spher << std::endl;

	// calculate gradient in Cartesian system at given point
	ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);
	Vector3Cartesian grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p1);
	std::cout << "Gradient in Cartesian coords:" << grad_cart << std::endl;

	// calculate gradient in Spherical system at given point
	ScalarFunction<3> fPotSpher(Fields::InverseRadialPotentialFieldSpher);
	Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(fPotSpher, p1Spher);
	std::cout << "Gradient in Spherical coords:" << grad_spher << std::endl;

	auto grad_cart_transf = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p1);
	std::cout << "Cartesian gradient from spherical (covariant transf.):" << grad_cart_transf << std::endl;

	auto grad_spher_transf = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p1Spher);
	std::cout << "Spherical gradient from Cartesian (covariant transf.):" << grad_spher_transf << std::endl;

	auto grad_cart_transf_inv = CoordTransfCartToSpher.transfInverseVecCovariant(grad_spher, p1);
	std::cout << "Cartesian gradient from spherical (inverse covariant transf.):" << grad_cart_transf_inv << std::endl;

	auto grad_spher_transf_inv = CoordTransfSpherToCart.transfInverseVecCovariant(grad_cart, p1Spher);
	std::cout << "Spherical gradient from spherical (inverse covariant transf.):" << grad_spher_transf_inv << std::endl;
}

void Demo_CoordTransf_Cylindrical()
{
	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "Direct cylindrical transformation\n";

	Vector3Cartesian p1{ 5.0, 2.0, -3.0 };
	auto p1Cyl = CoordTransfCartToCyl.transf(p1);
	auto p1BackTransf = CoordTransfCylToCart.transf(p1Cyl);

	std::cout << "Cartesian   : " << p1 << std::endl;
	std::cout << "Cylindrical : " << p1Cyl << std::endl;
	std::cout << "Back transf : " << p1BackTransf << std::endl;

	std::cout << "Inverse cylindrical transformation\n";

	Vector3Cartesian p2{ 5.0, 2.0, -3.0 };
	auto p2Cyl = CoordTransfCylToCart.transfInverse(p2);
	auto p2BackTransf = CoordTransfCartToCyl.transfInverse(p2Cyl);

	std::cout << "Cartesian   : " << p2 << std::endl;
	std::cout << "Spherical   : " << p2Cyl << std::endl;
	std::cout << "Back transf : " << p2BackTransf << std::endl;
}

void Demo_GetUnitVector()
{
	std::cout << "-----------          DEMO GET UNIT VECTOR          ----------\n";

	Vector3Cartesian p1{ 1.0, 2.0, 3 };
	auto   p1Spher = CoordTransfCartToSpher.transf(p1);

	Real r = p1Spher[0];
	Real theta = p1Spher[1];
	Real phi = p1Spher[2];

	std::cout << "Point (cartesian)       : " << p1 << std::endl;
	std::cout << "Point (spherical - rad) : " << p1Spher << std::endl;
	std::cout << "Point (spherical - deg) : "; p1Spher.PrintDeg(std::cout, 15, 10);

	std::cout << "\nCOVARIANT basis vectors - using getBasisVec(i)\n";
	Vector3Cartesian vec_r = CoordTransfSpherToCart.getBasisVec(0, p1Spher);
	Vector3Cartesian vec_theta = CoordTransfSpherToCart.getBasisVec(1, p1Spher);
	Vector3Cartesian vec_phi = CoordTransfSpherToCart.getBasisVec(2, p1Spher);
	std::cout << "basis(0), vec_r (Cart)      = " << vec_r << std::endl;
	std::cout << "basic(1), vec_theta (Cart)  = " << vec_theta << std::endl;
	std::cout << "basis(2), vec_phi (Cart)    = " << vec_phi << std::endl << std::endl;

	std::cout << "vec_r.GetAsUnitVector()     = " << vec_r.GetAsUnitVector() << std::endl;
	std::cout << "vec_theta.GetAsUnitVector() = " << vec_theta.GetAsUnitVector() << std::endl;
	std::cout << "vec_phi.GetAsUnitVector()   = " << vec_phi.GetAsUnitVector() << std::endl << std::endl;

	Vector3Cartesian vec_r2    {      sin(theta) * cos(phi),     sin(theta) * sin(phi),      cos(theta) };
	Vector3Cartesian vec_theta2{  r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), -r * sin(theta) };
	Vector3Cartesian vec_phi2  { -r * sin(theta) * sin(phi), r * sin(theta) * cos(phi)             ,0.0 };
	std::cout << "Explicit form for Sph.base vectors (r, theta, phi), expressed in terms of (i, j, k)\n";
	std::cout << "Unit vector r (Cart)        = " << vec_r2 << std::endl;
	std::cout << "Unit vector theta (Cart)    = " << vec_theta2 << std::endl;
	std::cout << "Unit vector phi (Cart)      = " << vec_phi2 << std::endl << std::endl;

	Vector3Cartesian vec_r1    {  sin(theta) * cos(phi), sin(theta) * sin(phi),  cos(theta) };
	Vector3Cartesian vec_theta1{  cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
	Vector3Cartesian vec_phi1  { -sin(phi)             , cos(phi)             ,         0.0 };
	std::cout << "Explicit form for UNIT Sph.base vectors (r, theta, phi), expressed in terms of (i, j, k)\n";
	std::cout << "Unit vector r (Cart)        = " << vec_r1 << std::endl;
	std::cout << "Unit vector theta (Cart)    = " << vec_theta1 << std::endl;
	std::cout << "Unit vector phi (Cart)      = " << vec_phi1 << std::endl << std::endl;

	std::cout << "\nCONTRAVARIANT basis vectors - using getContravarBasisVec(i)\n";
	Vector3Spherical vec_i = CoordTransfCartToSpher.getContravarBasisVec(0, p1Spher);
	Vector3Spherical vec_j = CoordTransfCartToSpher.getContravarBasisVec(1, p1Spher);
	Vector3Spherical vec_k = CoordTransfCartToSpher.getContravarBasisVec(2, p1Spher);
	std::cout << "basis(0), vec_i (Spher)           = " << vec_i << std::endl;
	std::cout << "basis(0), vec_j (Spher)           = " << vec_j << std::endl;
	std::cout << "basis(0), vec_k (Spher)           = " << vec_k << std::endl << std::endl;

	std::cout << "vec_i.GetAsUnitVectorAtPos(p1Sph) = " << vec_i.GetAsUnitVectorAtPos(p1Spher) << std::endl;
	std::cout << "vec_j.GetAsUnitVectorAtPos(p1Sph) = " << vec_j.GetAsUnitVectorAtPos(p1Spher) << std::endl;
	std::cout << "vec_k.GetAsUnitVectorAtPos(p1Sph) = " << vec_k.GetAsUnitVectorAtPos(p1Spher) << std::endl << std::endl;

	Vector3Spherical vec_i2{ sin(theta) * cos(phi), r * cos(theta) * cos(phi), -r * sin(theta) * sin(phi) };
	Vector3Spherical vec_j2{ sin(theta) * sin(phi), r * cos(theta) * sin(phi),  r * sin(theta) * cos(phi) };
	Vector3Spherical vec_k2{ cos(theta)           ,-r * sin(theta)           ,                 0.0 };
	std::cout << "Explicit form for Cart base vectors (i, j, k), expressed in terms of (r, theta, phi)\n";
	std::cout << "Unit vector i (Spher)  = " << vec_i2 << std::endl;
	std::cout << "Unit vector j (Spher)  = " << vec_j2 << std::endl;
	std::cout << "Unit vector k (Spher)  = " << vec_k2 << std::endl << std::endl;

	Vector3Spherical vec_i1{ sin(theta) * cos(phi),  cos(theta) * cos(phi), -sin(phi) };
	Vector3Spherical vec_j1{ sin(theta) * sin(phi),  cos(theta) * sin(phi),  cos(phi) };
	Vector3Spherical vec_k1{ cos(theta)           , -sin(theta)           ,  0.0 };
	std::cout << "Explicit form for UNIT Cart base vectors (i, j, k), expressed in terms of (r, theta, phi)\n";
	std::cout << "Unit vector i (Spher)  = " << vec_i1 << std::endl;
	std::cout << "Unit vector j (Spher)  = " << vec_j1 << std::endl;
	std::cout << "Unit vector k (Spher)  = " << vec_k1 << std::endl << std::endl;
}

void Demo_SpherToCart_BasisVectors()
{
	std::cout << "-----------------------------------------------------------------------\n";
	std::cout << "--            SPHERICAL BASIS (UNIT) VECTORS                         --\n";
	std::cout << "-----------------------------------------------------------------------\n";
	
	// https://physics.stackexchange.com/questions/546479/conversion-of-cartesian-position-and-velocity-to-spherical-velocity

	Vector3Cartesian pos_cart{  1.0, 2.0, 3.0 };
	Vector3Cartesian posCart2{ -1.0,-2.0,-3.0 };
	Vector3Cartesian posCart3{ -3.0, 2.0, 1.0 };

	Vector3Spherical posSpher1 = CoordTransfCartToSpher.transf(pos_cart);
	Real r     = posSpher1[0];
	Real theta = posSpher1[1];
	Real phi   = posSpher1[2];
	
	pos_cart.PrintLine(std::cout, "posCart1        = ", 10,3);
	posSpher1.PrintLine(std::cout, "posSpher1       = ", 10,3);
	std::cout << "posSpher1 (deg) = "; posSpher1.PrintDeg(std::cout, 10, 3);

	std::cout << "\nCOVARIANT BASIS VECTORS - at posSpher1, (e_r, e_theta, e_phi) expressed in terms of (i, j, k) = (e_x, e_y, e_z)\n";
	std::cout << "getBasisVec(i) - NOT normalized!!!\n";
	Vector3Cartesian e_r_cart     = CoordTransfSpherToCart.getBasisVec(0, posSpher1);
	Vector3Cartesian e_theta_cart = CoordTransfSpherToCart.getBasisVec(1, posSpher1);
	Vector3Cartesian e_phi_cart   = CoordTransfSpherToCart.getBasisVec(2, posSpher1);
	e_r_cart.PrintLine(std::cout, "e_r_cart           = ", 10, 3);
	e_theta_cart.PrintLine(std::cout, "e_theta_cart       = ", 10, 3);
	e_phi_cart.PrintLine(std::cout, "e_phi_cart         = ", 10, 3);

	//std::cout << "\nExplicit form for Spher.base vectors (e_r, e_theta, e_phi) \n";
	//Vector3Cartesian e_r_cart1 = CoordTransfSpherToCart.getBasisVectorExplicit(0, posSpher1);
	//Vector3Cartesian e_theta_cart1 = CoordTransfSpherToCart.getBasisVectorExplicit(1, posSpher1);
	//Vector3Cartesian e_phi_cart1 = CoordTransfSpherToCart.getBasisVectorExplicit(2, posSpher1);
	//e_r_cart1.PrintLine(std::cout, "e_r_cart           = ", 10, 3);
	//e_theta_cart1.PrintLine(std::cout, "e_theta_cart       = ", 10, 3);
	//e_phi_cart1.PrintLine(std::cout, "e_phi_cart         = ", 10, 3);

	std::cout << "\nNormalizing base vectors with GetAsUnitVector()\n";
	e_r_cart.GetAsUnitVector().PrintLine(std::cout, "e_r_cart.GetAs ... = ", 10, 3);
	e_theta_cart.GetAsUnitVector().PrintLine(std::cout, "e_theta_cart.GetAs = ", 10, 3);
	e_phi_cart.GetAsUnitVector().PrintLine(std::cout, "e_phi_cart.GetAs . = ", 10, 3);

	std::cout << "\nExplicit form for UNIT Spher.base vectors (e_r, e_theta, e_phi) \n";
	Vector3Cartesian e_r_cart_unit = CoordTransfSpherToCart.getUnitBasisVec(0, posSpher1);
	Vector3Cartesian e_theta_cart_unit = CoordTransfSpherToCart.getUnitBasisVec(1, posSpher1);
	Vector3Cartesian e_phi_cart_unit = CoordTransfSpherToCart.getUnitBasisVec(2, posSpher1);
	e_r_cart_unit.PrintLine(std::cout, "e_r_cart_unit1     = ", 10, 3);
	e_theta_cart_unit.PrintLine(std::cout, "e_theta_cart_unit1 = ", 10, 3);
	e_phi_cart_unit.PrintLine(std::cout, "e_phi_cart_unit1   = ", 10, 3);

	std::cout << "\nINVERSE COVARIANT BASIS VECTORS - at posSpher1, (i, j, k) = (e_x, e_y, e_z) expressed in terms of (e_r, e_theta, e_phi)\n";
	std::cout << "getInverseBasisVec(i) - NOT normalized\n";
	Vector3Spherical e_i_sph = CoordTransfSpherToCart.getInverseBasisVec(0, posSpher1);
	Vector3Spherical e_j_sph = CoordTransfSpherToCart.getInverseBasisVec(1, posSpher1);
	Vector3Spherical e_k_sph = CoordTransfSpherToCart.getInverseBasisVec(2, posSpher1);
	e_i_sph.PrintLine(std::cout, "e_i_sph            = ", 10, 3);
	e_j_sph.PrintLine(std::cout, "e_j_sph            = ", 10, 3);
	e_k_sph.PrintLine(std::cout, "e_k_sph            = ", 10, 3);

	//std::cout << "\nExplicit form for Cart.base vectors (i, j, k) = (e_x, e_y, e_z)\n";
	//Vector3Spherical e_i_sph1 = CoordTransfSpherToCart.getInverseBasisVectorExplicit(0, posSpher1);
	//Vector3Spherical e_j_sph1 = CoordTransfSpherToCart.getInverseBasisVectorExplicit(1, posSpher1);
	//Vector3Spherical e_k_sph1 = CoordTransfSpherToCart.getInverseBasisVectorExplicit(2, posSpher1);
	//e_i_sph1.PrintLine(std::cout, "e_i_sph            = ", 10, 3);
	//e_j_sph1.PrintLine(std::cout, "e_j_sph            = ", 10, 3);
	//e_k_sph1.PrintLine(std::cout, "e_k_sph            = ", 10, 3);

	//// We can obtain same set of vectors by using getContravarBasisVec, BUT WITH THE INVERSE TRANSFORMATION - CoordTransfCartToSpher!!!
	//std::cout << "\nInverse covariant basis vectors, at posSpher1, expressed in terms of (e_r, e_theta, e_phi)  - NOT normalized - using getContravarBasisVec(i)\n";
	//Vector3Spherical e_i_sph = CoordTransfCartToSpher.getContravarBasisVec(0, posSpher1);
	//Vector3Spherical e_j_sph = CoordTransfCartToSpher.getContravarBasisVec(1, posSpher1);
	//Vector3Spherical e_k_sph = CoordTransfCartToSpher.getContravarBasisVec(2, posSpher1);
	//e_i_sph.PrintLine(std::cout, "e_i_sph           = ", 10, 3);
	//e_j_sph.PrintLine(std::cout, "e_j_sph           = ", 10, 3);
	//e_k_sph.PrintLine(std::cout, "e_k_sph           = ", 10, 3);

	std::cout << "\nNormalizing inverse basis vectors with GetAsUnitVectorAtPos(posSpher1)\n";
	e_i_sph.GetAsUnitVectorAtPos(posSpher1).PrintLine(std::cout, "e_i_sph.GetAs .... = ", 10, 3);
	e_j_sph.GetAsUnitVectorAtPos(posSpher1).PrintLine(std::cout, "e_j_sph.GetAs .... = ", 10, 3);
	e_k_sph.GetAsUnitVectorAtPos(posSpher1).PrintLine(std::cout, "e_k_sph.GetAs .... = ", 10, 3);

	std::cout << "\nExplicit form for UNIT Cart.base vectors (i, j, k) = (e_x, e_y, e_z)\n";
	Vector3Spherical e_i_sph_unit = CoordTransfSpherToCart.getInverseUnitBasisVec(0, posSpher1);
	Vector3Spherical e_j_sph_unit = CoordTransfSpherToCart.getInverseUnitBasisVec(1, posSpher1);
	Vector3Spherical e_k_sph_unit = CoordTransfSpherToCart.getInverseUnitBasisVec(2, posSpher1);
	e_i_sph_unit.PrintLine(std::cout, "e_i_sph_unit       = ", 10, 3);
	e_j_sph_unit.PrintLine(std::cout, "e_j_sph_unit       = ", 10, 3);
	e_k_sph_unit.PrintLine(std::cout, "e_k_sph_unit       = ", 10, 3);

	std::cout << "\nExpanding unit cart. basis vectors in terms of covariant spher. basis vectors\n";
	auto vec_i = e_i_sph_unit[0] * e_r_cart_unit + e_i_sph_unit[1] * e_theta_cart_unit + e_i_sph_unit[2] * e_phi_cart_unit;
	auto vec_j = e_j_sph_unit[0] * e_r_cart_unit + e_j_sph_unit[1] * e_theta_cart_unit + e_j_sph_unit[2] * e_phi_cart_unit;
	auto vec_k = e_k_sph_unit[0] * e_r_cart_unit + e_k_sph_unit[1] * e_theta_cart_unit + e_k_sph_unit[2] * e_phi_cart_unit;

	vec_i.PrintLine(std::cout, "vec_i     = e_i_sph_unit[0]      * e_r_cart_unit + e_i_sph_unit[1]      * e_theta_cart_unit + e_i_sph_unit[2]      * e_phi_cart_unit = ", 10, 3);
	vec_j.PrintLine(std::cout, "vec_j     = e_j_sph_unit[0]      * e_r_cart_unit + e_j_sph_unit[1]      * e_theta_cart_unit + e_j_sph_unit[2]      * e_phi_cart_unit = ", 10, 3);
	vec_k.PrintLine(std::cout, "vec_k     = e_k_sph_unit[0]      * e_r_cart_unit + e_k_sph_unit[1]      * e_theta_cart_unit + e_k_sph_unit[2]      * e_phi_cart_unit = ", 10, 3);

	std::cout << "\nExpanding unit spher. basis vectors in terms of covariant Cart. basis vectors\n";
	auto vec_r	   = e_r_cart_unit[0] * e_i_sph_unit     + e_r_cart_unit[1] * e_j_sph_unit     + e_r_cart_unit[2] * e_k_sph_unit;
	auto vec_theta = e_theta_cart_unit[0] * e_i_sph_unit + e_theta_cart_unit[1] * e_j_sph_unit + e_theta_cart_unit[2] * e_k_sph_unit;
	auto vec_phi	 = e_phi_cart_unit[0] * e_i_sph_unit   + e_phi_cart_unit[1] * e_j_sph_unit   + e_phi_cart_unit[2] * e_k_sph_unit;

	vec_r.PrintLine    (std::cout, "vec_r     = e_r_cart_unit[0]     * e_i_sph_unit  + e_r_cart_unit[1]     * e_j_sph_unit      + e_r_cart_unit[2]     * e_k_sph_unit    = ", 10, 3);
	vec_theta.PrintLine(std::cout, "vec_theta = e_theta_cart_unit[0] * e_i_sph_unit  + e_theta_cart_unit[1] * e_j_sph_unit      + e_theta_cart_unit[2] * e_k_sph_unit    = ", 10, 3);
	vec_phi.PrintLine  (std::cout, "vec_phi   = e_phi_cart_unit[0]   * e_i_sph_unit  + e_phi_cart_unit[1]   * e_j_sph_unit      + e_phi_cart_unit[2]   * e_k_sph_unit    = ", 10, 3);

	MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3>   metricSpherToCart(CoordTransfSpherToCart);

	auto metricAtPos2 = metricSpherToCart(posSpher1);
	std::cout << "Metric at pos:" << std::endl;
	metricAtPos2.Print(std::cout, 10, 4);

	MetricTensorSpherical metricSpher;
	auto spher_metric = metricSpher(posSpher1);
	std::cout << "Spherical Metric Tensor:" << spher_metric << std::endl;

	Vector3Cartesian v_cart{ 3.0, 4.0, 5.0 };
	auto speedUnitCart = v_cart.GetAsUnitVector();

	v_cart.PrintLine(std::cout, "\nv_cart          = ", 10, 4);

	Vector3Spherical v_spher_contra = CoordTransfCartToSpher.transfVecContravariant(v_cart, pos_cart);
	std::cout << "transfVecContravariant(v_cart, pos_cart):" << std::endl;
	v_spher_contra.PrintLine(std::cout, "v_spher_contra  = ", 10, 4);

	std::cout << "\nExpanding spherical speed vector in terms of covariant spher. basis vectors\n";
	auto vec_r1 = v_spher_contra[0] * e_r_cart + v_spher_contra[1] * e_theta_cart + v_spher_contra[2] * e_phi_cart;
	vec_r1.PrintLine(std::cout, "vec_r = v_spher_contra[0] * e_r_cart + v_spher_contra[1] * e_theta_cart + v_spher_contra[2] * e_phi_cart = ", 10, 4);

	std::cout << "\nProjecting v_cart on spher. unit basis vectors:\n";
	Vector3Spherical vec_covar{e_r_cart_unit * v_cart, e_theta_cart_unit * v_cart, e_phi_cart_unit * v_cart};
	vec_covar.PrintLine(std::cout, "vec_covar = Vector3Spherical(e_r_cart_unit * v_cart, e_theta_cart_unit * v_cart, e_phi_cart_unit * v_cart)     = ", 10, 4);
}

void Demo_CoordTransf()
{
	std::cout << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                         COORD TRANSF                          ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;


	//Demo_CoordTrans_RectilinearOrthogonal();
	//Demo_CoordTransf_RectilinearOblique();

	//Demo_2DPolar();

	//Demo_CoordTransf_Spherical();
	//Demo_CoordTransf_Cylindrical();
	//Demo_GetUnitVector ();
	Demo_SpherToCart_BasisVectors();

	//Demo_CoordTransf_Jacobian();
}