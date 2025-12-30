#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/Quaternions.h"
#include "base/Tensor.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/CoordTransf/CoordTransf3D.h"
#include "core/CoordTransf/CoordTransfLorentz.h"
#include "core/Fields.h"

#endif


using namespace MML;
using namespace MML::Utils;

void Docs_Demo_CoordTransf_Transf()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****           DEMO COORDINATE TRANSFORMATIONS                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;
  
  std::cout << "******      Direct spherical transformation     *****\n";

  Vector3Cartesian p1{ 5.0, 2.0, -3.0 };
  auto p1Spher = CoordTransfCartToSpher.transf(p1);
  auto p1BackTransf = CoordTransfSpherToCart.transf(p1Spher);

  std::cout << "Cartesian   : " << p1 << std::endl;
  std::cout << "Spherical   : " << p1Spher << std::endl;
  std::cout << "Back transf : " << p1BackTransf << std::endl;

  std::cout << "*****       Inverse spherical transformation   *****\n";

  Vector3Cartesian p2{ 5.0, 2.0, -3.0 };
  auto p2Spher = CoordTransfSpherToCart.transfInverse(p2);
  auto p2BackTransf = CoordTransfCartToSpher.transfInverse(p2Spher);

  std::cout << "Cartesian   : " << p2 << std::endl;
  std::cout << "Spherical   : " << p2Spher << std::endl;
  std::cout << "Back transf : " << p2BackTransf << std::endl;

  std::cout << "*****       Direct cylindrical transformation      *****\n";

  Vector3Cartesian p3{ 5.0, 2.0, -3.0 };
  auto p3Cyl = CoordTransfCartToCyl.transf(p3);
  auto p3BackTransf = CoordTransfCylToCart.transf(p3Cyl);

  std::cout << "Cartesian   : " << p3 << std::endl;
  std::cout << "Cylindrical : " << p3Cyl << std::endl;
  std::cout << "Back transf : " << p3BackTransf << std::endl;

  std::cout << "*****      Inverse cylindrical transformation      *****\n";

  Vector3Cartesian p4{ 5.0, 2.0, -3.0 };
  auto p4Cyl = CoordTransfCylToCart.transfInverse(p4);
  auto p4BackTransf = CoordTransfCartToCyl.transfInverse(p4Cyl);

  std::cout << "Cartesian   : " << p4 << std::endl;
  std::cout << "Spherical   : " << p4Cyl << std::endl;
  std::cout << "Back transf : " << p4BackTransf << std::endl;

/* OUTPUT
******      Direct spherical transformation     *****
Cartesian   : [              5,               2,              -3]
Spherical   : [    6.164414003,     2.079063573,    0.3805063771]
Back transf : [              5,               2,              -3]
*****       Inverse spherical transformation   *****
Cartesian   : [              5,               2,              -3]
Spherical   : [    6.164414003,     2.079063573,    0.3805063771]
Back transf : [              5,               2,              -3]
*****       Direct cylindrical transformation      *****
Cartesian   : [              5,               2,              -3]
Cylindrical : [    5.385164807,    0.3805063771,              -3]
Back transf : [              5,               2,              -3]
*****      Inverse cylindrical transformation      *****
Cartesian   : [              5,               2,              -3]
Spherical   : [    5.385164807,    0.3805063771,              -3]
Back transf : [              5,               2,              -3]
*/
}

void Docs_Demo_CoordTransf_CovarContravar_Vec_transf()
{
  std::cout << "\n***********************************************************************" << std::endl;
  std::cout << "****          CONTRAVARIANT TRANSFORMATION OF VELOCITY             ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  Vec3Cart v_cart{ 1.0, 1.0, 0.0 };
  std::cout << "v_cart    : " << v_cart << std::endl;

  std::cout << "AT POINT:\n";
  Vector3Cartesian x1_cart{ 1.0, -2.0, 1.0 };
  Vector3Spherical x1_spher{ CoordTransfCartToSpher.transf(x1_cart) };

  std::cout << "Cartesian : " << x1_cart << std::endl << "Spherical : " << x1_spher << std::endl;

  Vector3Spherical v_transf_to_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, x1_cart);
  std::cout << "contravar transf. to spher at cart.pnt = " << v_transf_to_spher << std::endl;

  Vector3Cartesian v_back_transf_to_cart = CoordTransfSpherToCart.transfVecContravariant(v_transf_to_spher, x1_spher);
  std::cout << "back transf. to cartesian at spher.pnt = " << v_back_transf_to_cart << std::endl;

  std::cout << "\n***********************************************************************" << std::endl;
  std::cout << "****             COVARIANT TRANSFORMATION OF GRADIENT              ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  Vector3Cartesian p_cart{ 1.0, 1.0, 1.0 };
  Vector3Spherical p_spher(CoordTransfSpherToCart.transfInverse(p_cart));

  std::cout << "Spherical: " << p_spher << std::endl;
  std::cout << "Cartesian: " << p_cart << std::endl;

  std::cout << "\n******   Working with SPHERICAL TO CARTESIAN transformation   *******\n";
  std::cout << "\n******   Potential in spherical coordinates  *******\n";

  ScalarFunction<3> fPotSpher(Fields::InverseRadialPotentialFieldSpher);
  Vector3Spherical  grad_spher = ScalarFieldOperations::GradientSpher(fPotSpher, p_spher);

  std::cout << "Field at     : " << p_spher << " = " << Fields::InverseRadialPotentialFieldSpher(p_spher) << std::endl;
  std::cout << "Grad(sph) at : " << p_spher << " = " << grad_spher << std::endl;

  Vector3Cartesian grad_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p_cart);
  Vector3Spherical back_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_transf_to_cart, p_spher);

  std::cout << "Grad.transf. (Cart.) at " << p_cart << " = " << grad_transf_to_cart << std::endl;
  std::cout << "Back transf. (Spher) at " << p_spher << " = " << back_transf_to_spher << std::endl;

  std::cout << "\n******   Working with CARTESIAN TO SPHERICAL transformation   *******\n";
  std::cout << "\n******   Potential in cartesian coordinates:  *******\n";

  ScalarFunction<3> fPotCart(Fields::InverseRadialPotentialFieldCart);
  Vector3Cartesian  grad_cart = ScalarFieldOperations::GradientCart<3>(fPotCart, p_cart);

  std::cout << "Field at      : " << p_cart << " = " << Fields::InverseRadialPotentialFieldCart(p_cart) << std::endl;
  std::cout << "Grad.Cart. at : " << p_cart << " = " << grad_cart << std::endl;

  Vector3Spherical grad_transf_to_spher = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
  Vector3Cartesian back_transf_to_cart = CoordTransfSpherToCart.transfVecCovariant(grad_transf_to_spher, p_cart);

  std::cout << "Grad.transf. (Spher) at " << p_spher << " = " << grad_transf_to_spher << std::endl;
  std::cout << "Back transf. (Cart.) at " << p_cart << " = " << back_transf_to_cart << std::endl;

/* OUTPUT
***********************************************************************
****          CONTRAVARIANT TRANSFORMATION OF VELOCITY             ****
***********************************************************************
v_cart    : [              1,               1,               0]
AT POINT:
Cartesian : [              1,              -2,               1]
Spherical : [    2.449489743,     1.150261992,    -1.107148718]
contravar transf. to spher at cart.pnt = [  -0.4082482796,  -0.07453559991,    0.6000000052]
back transf. to cartesian at spher.pnt = [    1.000000006,    0.9999999863, 7.274731734e-09]

***********************************************************************
****             COVARIANT TRANSFORMATION OF GRADIENT              ****
***********************************************************************
Spherical: [    1.732050808,    0.9553166181,    0.7853981634]
Cartesian: [              1,               1,               1]

******   Working with SPHERICAL TO CARTESIAN transformation   *******

******   Potential in spherical coordinates  *******
Field at     : [    1.732050808,    0.9553166181,    0.7853981634] = 0.5773502692
Grad(sph) at : [    1.732050808,    0.9553166181,    0.7853981634] = [  -0.3333333333,               0,               0]
Grad.transf. (Cart.) at [              1,               1,               1] = [  -0.1924500897,   -0.1924500897,   -0.1924500897]
Back transf. (Spher) at [    1.732050808,    0.9553166181,    0.7853981634] = [  -0.3333333333, 2.942091015e-15, 2.370326158e-14]

******   Working with CARTESIAN TO SPHERICAL transformation   *******

******   Potential in cartesian coordinates:  *******
Field at      : [              1,               1,               1] = 0.5773502692
Grad.Cart. at : [              1,               1,               1] = [  -0.1924500897,   -0.1924500897,   -0.1924500897]
Grad.transf. (Spher) at [    1.732050808,    0.9553166181,    0.7853981634] = [  -0.3333333333, 2.942091015e-15, 2.370326158e-14]
Back transf. (Cart.) at [              1,               1,               1] = [  -0.1924500897,   -0.1924500897,   -0.1924500897]
*/
}

void Docs_Demo_CoordTransf_CovarContravarBasis()
{
  std::cout << "\n***********************************************************************" << std::endl;
  std::cout << "****         DEMO COVARIANT & CONTRAVARIANT BASIS VECTORS          ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  Vector3Cartesian p1{ 2.0, 2.0, 5 };
  auto   p1Spher = CoordTransfCartToSpher.transf(p1);

  Real r = p1Spher[0];
  Real theta = p1Spher[1];
  Real phi = p1Spher[2];

  std::cout << "Point (cartesian)       : " << p1 << std::endl;
  std::cout << "Point (spherical - rad) : " << p1Spher << std::endl;
  std::cout << "Point (spherical - deg) : "; p1Spher.PrintDeg(std::cout, 15, 10);

  std::cout << "\n*****       Using Spherical to Cartesian transformation      *****\n";

  std::cout << "\nCOVARIANT basis vectors\n";
  Vector3Cartesian vec_co_1 = CoordTransfSpherToCart.getBasisVec(0, p1Spher);
  Vector3Cartesian vec_co_2 = CoordTransfSpherToCart.getBasisVec(1, p1Spher);
  Vector3Cartesian vec_co_3 = CoordTransfSpherToCart.getBasisVec(2, p1Spher);
  std::cout << "getBasisVec 1 : " << vec_co_1 << std::endl;
  std::cout << "getBasisVec 2 : " << vec_co_2 << std::endl;
  std::cout << "getBasisVec 3 : " << vec_co_3 << std::endl << std::endl;

  Vector3Cartesian vec_r2    {      sin(theta) * cos(phi),     sin(theta) * sin(phi),      cos(theta) };
  Vector3Cartesian vec_theta2{  r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), -r * sin(theta) };
  Vector3Cartesian vec_phi2  { -r * sin(theta) * sin(phi), r * sin(theta) * cos(phi),      0.0 };
  
  //std::cout << "Explicit formula 1     : " << vec_r2 << std::endl;
  //std::cout << "Explicit formula 2     : " << vec_theta2 << std::endl;
  //std::cout << "Explicit formula 3     : " << vec_phi2 << std::endl << std::endl;

  std::cout << "GetAsUnitVector 1      : " << vec_co_1.GetAsUnitVector() << std::endl;
  std::cout << "GetAsUnitVector 2      : " << vec_co_2.GetAsUnitVector() << std::endl;
  std::cout << "GetAsUnitVector 3      : " << vec_co_3.GetAsUnitVector() << std::endl << std::endl;

  std::cout << "CONTRAVARIANT basis vectors\n";
  Vector3Cartesian vec_contra_1 = CoordTransfSpherToCart.getContravarBasisVec(0, p1);
  Vector3Cartesian vec_contra_2 = CoordTransfSpherToCart.getContravarBasisVec(1, p1);
  Vector3Cartesian vec_contra_3 = CoordTransfSpherToCart.getContravarBasisVec(2, p1);
  std::cout << "getContravarBasisVec 1 : " << vec_contra_1 << std::endl;
  std::cout << "getContravarBasisVec 2 : " << vec_contra_2 << std::endl;
  std::cout << "getContravarBasisVec 3 : " << vec_contra_3 << std::endl << std::endl;

  Vector3Cartesian vec_i2{        sin(theta) * cos(phi),       sin(theta) * sin(phi),        cos(theta) };
  Vector3Cartesian vec_j2{  1/r * cos(theta) * cos(phi), 1/r * cos(theta) * sin(phi), -1/r * sin(theta) };
  Vector3Cartesian vec_k2{ -1/r * sin(phi) / sin(theta), 1/r * cos(phi) / sin(theta),        0.0 };

  //std::cout << "Explicit formula 1         : " << vec_i2 << std::endl;
  //std::cout << "Explicit formula 2         : " << vec_j2 << std::endl;
  //std::cout << "Explicit formula 3         : " << vec_k2 << std::endl << std::endl;

  std::cout << "GetAsUnitVectorAtPos 1     : " << vec_contra_1.GetAsUnitVectorAtPos(p1Spher) << std::endl;
  std::cout << "GetAsUnitVectorAtPos 2     : " << vec_contra_2.GetAsUnitVectorAtPos(p1Spher) << std::endl;
  std::cout << "GetAsUnitVectorAtPos 3     : " << vec_contra_3.GetAsUnitVectorAtPos(p1Spher) << std::endl << std::endl;


  std::cout << "Veifying orthogonality of covariant and contravariant basis vectors:\n";
  std::cout << ScalarProduct(vec_co_3, vec_contra_1) << std::endl;
  std::cout << ScalarProduct(vec_co_3, vec_contra_2) << std::endl;
  std::cout << ScalarProduct(vec_co_3, vec_contra_3) << std::endl;

  std::cout << "\n*****       Using Cartesian to Spherical transformation      *****\n";

  std::cout << "\nCOVARIANT basis vectors\n";
  Vector3Spherical vec2_co_1 = CoordTransfCartToSpher.getBasisVec(0, p1);
  Vector3Spherical vec2_co_2 = CoordTransfCartToSpher.getBasisVec(1, p1);
  Vector3Spherical vec2_co_3 = CoordTransfCartToSpher.getBasisVec(2, p1);
  std::cout << "getBasisVec 1 : " << vec2_co_1 << std::endl;
  std::cout << "getBasisVec 2 : " << vec2_co_2 << std::endl;
  std::cout << "getBasisVec 3 : " << vec2_co_3 << std::endl << std::endl;

  std::cout << "GetAsUnitVector 1      : " << vec2_co_1.GetAsUnitVectorAtPos(p1Spher) << std::endl;
  std::cout << "GetAsUnitVector 2      : " << vec2_co_2.GetAsUnitVectorAtPos(p1Spher) << std::endl;
  std::cout << "GetAsUnitVector 3      : " << vec2_co_3.GetAsUnitVectorAtPos(p1Spher) << std::endl << std::endl;

  std::cout << "CONTRAVARIANT basis vectors\n";
  Vector3Spherical vec2_contra_1 = CoordTransfCartToSpher.getContravarBasisVec(0, p1Spher);
  Vector3Spherical vec2_contra_2 = CoordTransfCartToSpher.getContravarBasisVec(1, p1Spher);
  Vector3Spherical vec2_contra_3 = CoordTransfCartToSpher.getContravarBasisVec(2, p1Spher);
  std::cout << "getContravarBasisVec 1 : " << vec2_contra_1 << std::endl;
  std::cout << "getContravarBasisVec 2 : " << vec2_contra_2 << std::endl;
  std::cout << "getContravarBasisVec 3 : " << vec2_contra_3 << std::endl << std::endl;

  std::cout << "GetAsUnitVectorAtPos 1     : " << vec2_contra_1.GetAsUnitVectorAtPos(p1) << std::endl;
  std::cout << "GetAsUnitVectorAtPos 2     : " << vec2_contra_2.GetAsUnitVectorAtPos(p1) << std::endl;
  std::cout << "GetAsUnitVectorAtPos 3     : " << vec2_contra_3.GetAsUnitVectorAtPos(p1) << std::endl << std::endl;

/* OUTPUT
Point (cartesian)       : [              2,               2,               5]
Point (spherical - rad) : [    5.744562647,    0.5148059551,    0.7853981634]
Point (spherical - deg) : [    5.7445626465, 29.4962084966, 45.0000000000 ]

*****       Using Spherical to Cartesian transformation      *****

COVARIANT basis vectors
getBasisVec 1 : [   0.3481553119,    0.3481553119,    0.8703882798]
getBasisVec 2 : [   3.5355339059,    3.5355339059,   -2.8284271247]
getBasisVec 3 : [  -2.0000000000,    2.0000000000,    0.0000000000]

GetAsUnitVector 1      : [   0.3481553119,    0.3481553119,    0.8703882798]
GetAsUnitVector 2      : [   0.6154574549,    0.6154574549,   -0.4923659639]
GetAsUnitVector 3      : [  -0.7071067812,    0.7071067812,    0.0000000000]

CONTRAVARIANT basis vectors
getContravarBasisVec 1 : [   0.3481553119,    0.3481553119,    0.8703882798]
getContravarBasisVec 2 : [   0.1071373911,    0.1071373911,   -0.0857099129]
getContravarBasisVec 3 : [  -0.2500000000,    0.2500000000,    0.0000000000]

GetAsUnitVectorAtPos 1     : [   0.3481553119,    0.3481553119,    0.8703882798]
GetAsUnitVectorAtPos 2     : [   0.6154574549,    0.6154574549,   -0.4923659639]
GetAsUnitVectorAtPos 3     : [  -0.7071067812,    0.7071067812,    0.0000000000]

Veifying orthogonality of covariant and contravariant basis vectors:
-0.0000000000
-0.0000000000
1.0000000000

*****       Using Cartesian to Spherical transformation      *****

COVARIANT basis vectors
getBasisVec 1 : [   0.3481553119,    0.1071373911,   -0.2500000000]
getBasisVec 2 : [   0.3481553119,    0.1071373911,    0.2500000000]
getBasisVec 3 : [   0.8703882798,   -0.0857099129,    0.0000000000]

GetAsUnitVector 1      : [   0.3481553119,    0.0186502259,   -0.0883883476]
GetAsUnitVector 2      : [   0.3481553119,    0.0186502259,    0.0883883476]
GetAsUnitVector 3      : [   0.8703882798,   -0.0149201807,    0.0000000000]

CONTRAVARIANT basis vectors
getContravarBasisVec 1 : [   0.3481553119,    3.5355339059,   -2.0000000000]
getContravarBasisVec 2 : [   0.3481553119,    3.5355339059,    2.0000000000]
getContravarBasisVec 3 : [   0.8703882798,   -2.8284271247,    0.0000000000]

GetAsUnitVectorAtPos 1     : [   0.3481553119,    1.7677669530,   -1.0997501703]
GetAsUnitVectorAtPos 2     : [   0.3481553119,    1.7677669530,    1.0997501703]
GetAsUnitVectorAtPos 3     : [   0.8703882798,   -1.4142135624,    0.0000000000]
*/
}

void Docs_Demo_CoordTransf_Tensor_transf()
{
    std::cout << "\n***********************************************************************" << std::endl;
    std::cout << "****                   TENSOR TRANSFORMATIONS                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Create a simple rank-2 tensor (stress-like) in Cartesian coordinates
    Tensor2<3> tensorCart(2, 0);  // Contravariant rank (2,0)
    tensorCart(0, 0) = 10.0;  // xx component
    tensorCart(1, 1) = 20.0;  // yy component  
    tensorCart(2, 2) = 30.0;  // zz component
    tensorCart(0, 1) = 5.0;   // xy component
    tensorCart(1, 0) = 5.0;   // yx component (symmetric)

    std::cout << "\nOriginal tensor (Cartesian, contravariant):" << std::endl;
    std::cout << "T^xx = " << tensorCart(0, 0) << ", T^yy = " << tensorCart(1, 1) 
              << ", T^zz = " << tensorCart(2, 2) << std::endl;
    std::cout << "T^xy = " << tensorCart(0, 1) << " (off-diagonal)" << std::endl;

    // Transform point
    Vector3Cartesian posCart{1.0, 1.0, 1.0};
    Vector3Spherical posSpher = CoordTransfCartToSpher.transf(posCart);

    std::cout << "\nPosition (Cart): " << posCart << std::endl;
    std::cout << "Position (Spher): " << posSpher << std::endl;

    // Transform tensor to spherical coordinates
    Tensor2<3> tensorSpher = CoordTransfCartToSpher.transfTensor2(tensorCart, posSpher);

    std::cout << "\nTransformed tensor (Spherical):" << std::endl;
    std::cout << "T^rr = " << tensorSpher(0, 0) << std::endl;
    std::cout << "T^θθ = " << tensorSpher(1, 1) << std::endl;
    std::cout << "T^φφ = " << tensorSpher(2, 2) << std::endl;
    std::cout << "T^rθ = " << tensorSpher(0, 1) << ", T^rφ = " << tensorSpher(0, 2) << std::endl;

    // Transform back to verify
    Tensor2<3> tensorBackToCart = CoordTransfSpherToCart.transfTensor2(tensorSpher, posCart);

    std::cout << "\nBack-transformed to Cartesian (verification):" << std::endl;
    std::cout << "T^xx = " << tensorBackToCart(0, 0) << " (expected: 10)" << std::endl;
    std::cout << "T^yy = " << tensorBackToCart(1, 1) << " (expected: 20)" << std::endl;
    std::cout << "T^zz = " << tensorBackToCart(2, 2) << " (expected: 30)" << std::endl;
    std::cout << "T^xy = " << tensorBackToCart(0, 1) << " (expected: 5)" << std::endl;
}

void Docs_Demo_CoordTransf_3D_Rotations()
{
    std::cout << "\n***********************************************************************" << std::endl;
    std::cout << "****                   3D ROTATION TRANSFORMATIONS                  ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Vec3Cart v{1.0, 0.0, 0.0};  // Unit vector along x-axis

    std::cout << "\nOriginal vector: " << v << std::endl;

    // Rotation about Z-axis by 90 degrees
    std::cout << "\n=== Rotation about Z-axis (90°) ===" << std::endl;
    CoordTransfCart3DRotationZAxis rotZ(Constants::PI / 2);
    Vec3Cart v_rotZ = rotZ.transf(v);
    std::cout << "After rotation: " << v_rotZ << std::endl;
    std::cout << "Expected: (0, 1, 0)" << std::endl;

    // Rotation about Y-axis by 90 degrees
    std::cout << "\n=== Rotation about Y-axis (90°) ===" << std::endl;
    CoordTransfCart3DRotationYAxis rotY(Constants::PI / 2);
    Vec3Cart v_rotY = rotY.transf(v);
    std::cout << "After rotation: " << v_rotY << std::endl;
    std::cout << "Expected: (0, 0, -1)" << std::endl;

    // Rotation about X-axis (vector along y)
    std::cout << "\n=== Rotation about X-axis (90°) ===" << std::endl;
    Vec3Cart vy{0.0, 1.0, 0.0};  // Vector along y-axis
    CoordTransfCart3DRotationXAxis rotX(Constants::PI / 2);
    Vec3Cart vy_rotX = rotX.transf(vy);
    std::cout << "Vector (0,1,0) rotated: " << vy_rotX << std::endl;
    std::cout << "Expected: (0, 0, 1)" << std::endl;

    // Composed rotations: first X, then Y, then Z
    std::cout << "\n=== Composed Rotations ===" << std::endl;
    CoordTransfCart3DRotationXAxis rx(Constants::PI / 6);  // 30°
    CoordTransfCart3DRotationYAxis ry(Constants::PI / 4);  // 45°
    CoordTransfCart3DRotationZAxis rz(Constants::PI / 3);  // 60°

    Vec3Cart p{1.0, 0.0, 0.0};
    Vec3Cart p1 = rx.transf(p);
    Vec3Cart p2 = ry.transf(p1);
    Vec3Cart p3 = rz.transf(p2);
    std::cout << "Original: " << p << std::endl;
    std::cout << "After Rx(30°): " << p1 << std::endl;
    std::cout << "After Ry(45°): " << p2 << std::endl;
    std::cout << "After Rz(60°): " << p3 << std::endl;

    // Inverse rotation verification
    std::cout << "\n=== Inverse Rotation (verification) ===" << std::endl;
    Vec3Cart p3_inv = rz.transfInverse(p3);
    Vec3Cart p2_inv = ry.transfInverse(p3_inv);
    Vec3Cart p1_inv = rx.transfInverse(p2_inv);
    std::cout << "Back through inverse: " << p1_inv << " (expected: 1,0,0)" << std::endl;
}

void Docs_Demo_CoordTransf_Quaternion()
{
    std::cout << "\n***********************************************************************" << std::endl;
    std::cout << "****                   QUATERNION ROTATIONS                         ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Create rotation: 90° about Z-axis using quaternion
    Vec3Cart axis{0.0, 0.0, 1.0};  // Z-axis
    Real angle = Constants::PI / 2;  // 90 degrees
    Quaternion q = Quaternion::FromAxisAngle(axis, angle);

    std::cout << "\nQuaternion from axis-angle:" << std::endl;
    std::cout << "Axis: " << axis << ", Angle: " << angle << " rad" << std::endl;
    std::cout << "q = " << q << std::endl;

    // Use quaternion rotation transform
    CoordTransfCart3DRotationQuaternion quatRot(q);

    Vec3Cart v{1.0, 0.0, 0.0};
    Vec3Cart rotated = quatRot.transf(v);
    std::cout << "\nRotate (1,0,0) by 90° about Z:" << std::endl;
    std::cout << "Result: " << rotated << " (expected: 0,1,0)" << std::endl;

    // Compose two quaternion rotations
    std::cout << "\n=== Quaternion Composition ===" << std::endl;
    Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart{1.0, 0.0, 0.0}, Constants::PI / 4);  // 45° about X
    Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart{0.0, 1.0, 0.0}, Constants::PI / 4);  // 45° about Y

    Quaternion composed = q2 * q1;  // q2 applied after q1
    CoordTransfCart3DRotationQuaternion composedRot(composed);

    Vec3Cart p{1.0, 0.0, 0.0};
    Vec3Cart p_composed = composedRot.transf(p);
    std::cout << "Original: " << p << std::endl;
    std::cout << "After composed rotation (45° X, then 45° Y): " << p_composed << std::endl;

    // Get rotation axis and angle from quaternion
    std::cout << "\n=== Extract Axis-Angle from Quaternion ===" << std::endl;
    Vec3Cart extractedAxis = quatRot.GetRotationAxis();
    Real extractedAngle = quatRot.GetRotationAngle();
    std::cout << "Extracted axis: " << extractedAxis << std::endl;
    std::cout << "Extracted angle: " << extractedAngle << " rad (" 
              << extractedAngle * 180 / Constants::PI << "°)" << std::endl;
}

void Docs_Demo_CoordTransf_Lorentz()
{
    std::cout << "\n***********************************************************************" << std::endl;
    std::cout << "****                   LORENTZ TRANSFORMATION                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Minkowski coordinates: (ct, x, y, z)
    // Velocity in units of c
    Real v = 0.6;  // 60% speed of light
    Real gamma = 1.0 / std::sqrt(1.0 - v * v);

    std::cout << "\nLorentz boost along x-axis:" << std::endl;
    std::cout << "Velocity v = " << v << "c" << std::endl;
    std::cout << "Lorentz factor γ = " << gamma << " (expected: " << 1.0/std::sqrt(1-0.36) << ")" << std::endl;

    CoordTransfLorentzXAxis lorentz(v);

    // Event in rest frame: (ct=0, x=1, y=0, z=0) - particle at x=1 at t=0
    Vector4Minkowski event{0.0, 1.0, 0.0, 0.0};
    std::cout << "\nEvent in rest frame S: (ct=" << event[0] << ", x=" << event[1] 
              << ", y=" << event[2] << ", z=" << event[3] << ")" << std::endl;

    // Transform to moving frame S'
    Vector4Minkowski eventPrime = lorentz.transf(event);
    std::cout << "Event in frame S' (moving at 0.6c):" << std::endl;
    std::cout << "(ct'=" << eventPrime[0] << ", x'=" << eventPrime[1] 
              << ", y'=" << eventPrime[2] << ", z'=" << eventPrime[3] << ")" << std::endl;

    // Length contraction demo
    std::cout << "\n=== Length Contraction ===" << std::endl;
    std::cout << "Object at rest in S: length = 1 meter" << std::endl;
    std::cout << "Length in S': L' = L/γ = " << 1.0 / gamma << " meters" << std::endl;
    std::cout << "x' coordinate: " << eventPrime[1] << " (should equal L/γ = " << 1.0/gamma << ")" << std::endl;

    // Time dilation demo
    std::cout << "\n=== Time Dilation ===" << std::endl;
    // Clock at x=0 in S, ticking from t=0 to t=1 (in units where c=1, so ct=1)
    Vector4Minkowski tick1{0.0, 0.0, 0.0, 0.0};  // Event 1: clock at origin, t=0
    Vector4Minkowski tick2{1.0, 0.0, 0.0, 0.0};  // Event 2: clock at origin, t=1

    Vector4Minkowski tick1Prime = lorentz.transf(tick1);
    Vector4Minkowski tick2Prime = lorentz.transf(tick2);

    Real dt = tick2[0] - tick1[0];  // Time in S
    Real dtPrime = tick2Prime[0] - tick1Prime[0];  // Time in S'

    std::cout << "Clock in S: ticks from t=0 to t=1 (Δt = " << dt << ")" << std::endl;
    std::cout << "Same clock in S': Δt' = " << dtPrime << " (expected: γ = " << gamma << ")" << std::endl;

    // Inverse transformation verification
    std::cout << "\n=== Inverse Transformation (verification) ===" << std::endl;
    Vector4Minkowski eventBack = lorentz.transfInverse(eventPrime);
    std::cout << "Original event: (ct=" << event[0] << ", x=" << event[1] << ")" << std::endl;
    std::cout << "Back-transformed: (ct=" << eventBack[0] << ", x=" << eventBack[1] << ")" << std::endl;
}

void Docs_Demo_CoordTransf_Jacobian()
{
    std::cout << "\n***********************************************************************" << std::endl;
    std::cout << "****                   JACOBIAN COMPUTATION                         ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Compute Jacobian for spherical to Cartesian transformation
    Vector3Spherical pos{2.0, Constants::PI / 4, Constants::PI / 6};  // r=2, θ=45°, φ=30°

    std::cout << "\nPosition (spherical): r=" << pos[0] << ", θ=" << pos[1] 
              << " (" << pos[1] * 180 / Constants::PI << "°)"
              << ", φ=" << pos[2] << " (" << pos[2] * 180 / Constants::PI << "°)" << std::endl;

    MatrixNM<Real, 3, 3> J = CoordTransfSpherToCart.jacobian(pos);

    std::cout << "\nJacobian matrix J^i_j = ∂x^i/∂q^j:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "  [";
        for (int j = 0; j < 3; j++) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << J(i, j);
        }
        std::cout << " ]" << std::endl;
    }

    // Jacobian determinant (volume element) - manual 3x3 determinant
    Real detJ = J(0,0) * (J(1,1)*J(2,2) - J(1,2)*J(2,1))
              - J(0,1) * (J(1,0)*J(2,2) - J(1,2)*J(2,0))
              + J(0,2) * (J(1,0)*J(2,1) - J(1,1)*J(2,0));
    Real r = pos[0];
    Real theta = pos[1];
    Real expectedDetJ = r * r * std::sin(theta);

    std::cout << "\nJacobian determinant:" << std::endl;
    std::cout << "det(J) = " << detJ << std::endl;
    std::cout << "Expected: r²sin(θ) = " << expectedDetJ << std::endl;

    std::cout << "\nVolume element transformation:" << std::endl;
    std::cout << "dV_cart = |det(J)| dV_spher = r²sin(θ) dr dθ dφ" << std::endl;
}

void Docs_Demo_Coord_Transf()
{
  Docs_Demo_CoordTransf_Transf();
  Docs_Demo_CoordTransf_CovarContravar_Vec_transf();
  Docs_Demo_CoordTransf_CovarContravarBasis();
  Docs_Demo_CoordTransf_Tensor_transf();
  Docs_Demo_CoordTransf_3D_Rotations();
  Docs_Demo_CoordTransf_Quaternion();
  Docs_Demo_CoordTransf_Lorentz();
  Docs_Demo_CoordTransf_Jacobian();
}