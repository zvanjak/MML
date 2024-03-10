#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordSystem.h"
#include "core/CoordTransf.h"
#endif

using namespace MML;

void Demo_GetUnitVector()
{
    std::cout << "-----------          DEMO GET UNIT VECTOR          ----------\n";

    Vector3Cartesian p1{2.0, 2.0, 5};
    auto   p1Spher = CoordTransfCartToSpher.transf(p1);
    
    Real r     = p1Spher[0];
    Real theta = p1Spher[1];
    Real phi   = p1Spher[2];

    std::cout << "Point (cartesian)       : " << p1 << std::endl;
    std::cout << "Point (spherical - rad) : " << p1Spher << std::endl;
    std::cout << "Point (spherical - deg) : "; p1Spher.PrintDeg(std::cout, 15, 10); 

    std::cout << "\nCOVARIANT basis vectors\n";
    Vector3Cartesian vec_r     = CoordTransfSpherToCart.getCovariantBasisVec(0, p1Spher);
    Vector3Cartesian vec_theta = CoordTransfSpherToCart.getCovariantBasisVec(1, p1Spher);
    Vector3Cartesian vec_phi   = CoordTransfSpherToCart.getCovariantBasisVec(2, p1Spher);
    std::cout << "Unit vector r     (calc)      : " << vec_r << std::endl;
    std::cout << "Unit vector theta (calc)      : " << vec_theta << std::endl;
    std::cout << "Unit vector phi   (calc)      : " << vec_phi << std::endl<< std::endl;
    
    Vector3Cartesian vec_r2     {                sin(theta) * cos(phi),                sin(theta) * sin(phi),      cos(theta) };
    Vector3Cartesian vec_theta2 { r *            cos(theta) * cos(phi), r *            cos(theta) * sin(phi), -r * sin(theta) };
    Vector3Cartesian vec_phi2   {-r*sin(theta) * sin(phi)             , r*sin(theta) * cos(phi)             ,      0.0 };
    std::cout << "Unit vector r                 : " << vec_r2 << std::endl;
    std::cout << "Unit vector theta             : " << vec_theta2 << std::endl;
    std::cout << "Unit vector phi               : " << vec_phi2 << std::endl << std::endl;

    std::cout << "Unit vector r     unit (calc) : " << vec_r.GetAsUnitVector() << std::endl;
    std::cout << "Unit vector theta unit (calc) : " << vec_theta.GetAsUnitVector() << std::endl;
    std::cout << "Unit vector phi   unit (calc) : " << vec_phi.GetAsUnitVector() << std::endl<< std::endl;

    Vector3Cartesian vec_r1     {  sin(theta) * cos(phi), sin(theta) * sin(phi),  cos(theta) };
    Vector3Cartesian vec_theta1 {  cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
    Vector3Cartesian vec_phi1   { -sin(phi)             , cos(phi)             ,  0.0 };
    std::cout << "Explicit form for UNIT Sph.base vectors (r, theta, phi), expressed in terms of (i, j, k)\n";
    std::cout << "Unit vector r                 : " << vec_r1 << std::endl;
    std::cout << "Unit vector theta             : " << vec_theta1 << std::endl;
    std::cout << "Unit vector phi               : " << vec_phi1 << std::endl << std::endl;

    std::cout << "\nCONTRAVARIANT basis vectors\n";
    Vector3Spherical vec_i = CoordTransfCartToSpher.getContravariantBasisVec(0, p1Spher);
    Vector3Spherical vec_j = CoordTransfCartToSpher.getContravariantBasisVec(1, p1Spher);
    Vector3Spherical vec_k = CoordTransfCartToSpher.getContravariantBasisVec(2, p1Spher);
    std::cout << "Unit vector i (calc)      : " << vec_i << std::endl;
    std::cout << "Unit vector j (calc)      : " << vec_j << std::endl;
    std::cout << "Unit vector k (calc)      : " << vec_k << std::endl << std::endl;
    
    Vector3Spherical vec_i2     { sin(theta) * cos(phi), r * cos(theta) * cos(phi), -r*sin(theta) * sin(phi) };
    Vector3Spherical vec_j2     { sin(theta) * sin(phi), r * cos(theta) * sin(phi),  r*sin(theta) * cos(phi) };
    Vector3Spherical vec_k2     { cos(theta)           ,-r * sin(theta)           ,                 0.0 };
    std::cout << "Unit vector i             : " << vec_i2 << std::endl;
    std::cout << "Unit vector j             : " << vec_j2 << std::endl;
    std::cout << "Unit vector k             : " << vec_k2 << std::endl << std::endl;

    std::cout << "Unit vector i unit (calc) : " << vec_i.GetAsUnitVectorAtPos(p1Spher) << std::endl;
    std::cout << "Unit vector j unit (calc) : " << vec_j.GetAsUnitVectorAtPos(p1Spher) << std::endl;
    std::cout << "Unit vector k unit (calc) : " << vec_k.GetAsUnitVectorAtPos(p1Spher) << std::endl << std::endl;

    Vector3Spherical vec_i1     { sin(theta) * cos(phi),  cos(theta) * cos(phi), -sin(phi) };
    Vector3Spherical vec_j1     { sin(theta) * sin(phi),  cos(theta) * sin(phi),  cos(phi) };
    Vector3Spherical vec_k1     { cos(theta)           , -sin(theta)           ,  0.0 };
    std::cout << "Explicit form for UNIT Cart base vectors (i, j, k), expressed in terms of (r, theta, phi)\n";
    std::cout << "Unit vector i             : " << vec_i1 << std::endl;
    std::cout << "Unit vector j             : " << vec_j1 << std::endl;
    std::cout << "Unit vector k             : " << vec_k1 << std::endl << std::endl;
}

void Demo_CoordTransf_Spherical()
{
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "Direct spherical transformation\n";

    Vector3Cartesian p1{5.0, 2.0, -3.0};
    auto p1Spher      = CoordTransfCartToSpher.transf(p1);
    auto p1BackTransf = CoordTransfSpherToCart.transf(p1Spher);

    std::cout << "Cartesian   : " << p1 << std::endl;
    std::cout << "Spherical   : " << p1Spher << std::endl;
    std::cout << "Back transf : " << p1BackTransf << std::endl;

    std::cout << "Inverse spherical transformation\n";
    
    Vector3Cartesian p2{5.0, 2.0, -3.0};
    auto p2Spher      = CoordTransfSpherToCart.transfInverse(p2);
    auto p2BackTransf = CoordTransfCartToSpher.transfInverse(p2Spher);

    std::cout << "Cartesian   : " << p2 << std::endl;
    std::cout << "Spherical   : " << p2Spher << std::endl;
    std::cout << "Back transf : " << p2BackTransf << std::endl;
}

void Demo_CoordTransf_Cylindrical()
{
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "Direct cylindrical transformation\n";

    Vector3Cartesian p1{5.0, 2.0, -3.0};
    auto p1Cyl        = CoordTransfCartToCyl.transf(p1);
    auto p1BackTransf = CoordTransfCylToCart.transf(p1Cyl);

    std::cout << "Cartesian   : " << p1 << std::endl;
    std::cout << "Cylindrical : " << p1Cyl << std::endl;
    std::cout << "Back transf : " << p1BackTransf << std::endl;

    std::cout << "Inverse cylindrical transformation\n";
    
    Vector3Cartesian p2{5.0, 2.0, -3.0};
    auto p2Cyl        = CoordTransfCylToCart.transfInverse(p2);
    auto p2BackTransf = CoordTransfCartToCyl.transfInverse(p2Cyl);

    std::cout << "Cartesian   : " << p2 << std::endl;
    std::cout << "Spherical   : " << p2Cyl << std::endl;
    std::cout << "Back transf : " << p2BackTransf << std::endl;
}

void Demo_2DPolar()
{
  std::cout << "***************    OBLIQUE RECTILINEAR   ****************" << std::endl;

  // new basis vectors (keeping z coordinate out of it ;)
  Vector3Cartesian e1_base{ 1, 3, 0 }, e2_base{ 4, 0, 0 }, e3_base{ 0, 0, 1 };

  CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

  std::cout << "Is right-handed? - " << transf.IsRightHanded() << std::endl;

  std::cout << "New basis vectors:\n";
  std::cout << "e1 = "; e1_base.Print(std::cout, 10, 5) << std::endl;
  std::cout << "e2 = "; e2_base.Print(std::cout, 10, 5) << std::endl;
  std::cout << "e3 = "; e3_base.Print(std::cout, 10, 5) << std::endl;

  std::cout << "Calculated dual basis vectors:\n";
  std::cout << "e1 dual = "; transf.Dual(0).Print(std::cout, 10, 5) << std::endl;
  std::cout << "e2 dual = "; transf.Dual(1).Print(std::cout, 10, 5) << std::endl;
  std::cout << "e3 dual = "; transf.Dual(2).Print(std::cout, 10, 5) << std::endl;

  std::cout << "\nTransf. matrix:\n" << transf.getAlpha() << std::endl;
  std::cout << "Inverse transf. matrix:\n" << transf.getTransf() << std::endl;

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

void Demo_Coord_Rectilinear()
{
  std::cout << "\n****    RECTILINEAR COORD SYSTEM   ****" << std::endl;

  //    Vector3Cartesian e1_base{1, 0, 0}, e2_base{0, 1, 0}, e3_base{0, 0.1, 1};      // ovo definira dot koord system
  Vector3Cartesian e1_base{ 1, 3, 0 }, e2_base{ 4, 0, 0 }, e3_base{ 0, 0, 1 };

  CoordTransfRectilinear transf(e1_base, e2_base, e3_base);

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

  Vector3Cartesian covar_coef{ ScalarProd(vec_A, Vector3Cartesian(e1_base)),
                                   ScalarProd(vec_A, Vector3Cartesian(e2_base)),
                                   ScalarProd(vec_A, Vector3Cartesian(e3_base)) };
  std::cout << "Covar. coeff.    : " << covar_coef << std::endl;

  VectorN<Real, 3> covar_expanded = covar_coef[0] * transf.Dual(0) +
    covar_coef[1] * transf.Dual(1) +
    covar_coef[2] * transf.Dual(2);
  std::cout << "Expanded to orig.: " << covar_expanded << std::endl;
}

void Demo_Coord_OrthogonalCartesian()
{
  std::cout << "\n****    ORTHOGONAL CARTESIAN COORD SYSTEM   ****" << std::endl;

  Vector3Cartesian b1{ 1, 0, 0 }, b2{ 0, 1, 0 }, b3{ 0, 0, 1 };

  CoordSystemOrthogonalCartesian csys(b1, b2, b3);

  std::cout << "Is orthogonal = " << csys.isOrthogonal() << std::endl;
}

void Demo_CoordTransf_Jacobian()
{
    Vector3Cartesian p1{2.0, 2.0, 5};
    auto   p1Spher = CoordTransfCartToSpher.transf(p1);
    
    MatrixNM<Real, 3, 3> jac = CoordTransfCartToSpher.jacobian(p1Spher);    
}

void Demo_CoordTransf()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         COORD TRANSF                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_CoordTransf_Spherical();
    Demo_CoordTransf_Cylindrical();
    Demo_GetUnitVector ();

    Demo_2DPolar();
    Demo_Coord_Rectilinear();
    Demo_Coord_OrthogonalCartesian();

    Demo_CoordTransf_Jacobian();
}