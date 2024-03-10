# Coordinate transformations

Defined interfaces:

~~~C++
template<typename VectorFrom, typename VectorTo, int N>
class ICoordTransf
{
public:
	virtual       VectorTo            transf(const VectorFrom& in) const = 0;
	virtual const IScalarFunction<N>& coordTransfFunc(int i) const = 0;
};

template<typename VectorFrom, typename VectorTo, int N>
class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
{
public:
	virtual       VectorFrom          transfInverse(const VectorTo& in) const = 0;
	virtual const IScalarFunction<N>& inverseCoordTransfFunc(int i) const = 0;
};
~~~

Abstract realizations of those interfaces:

~~~C++
template<typename VectorFrom, typename VectorTo, int N>
class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
{
public:
    virtual VectorTo getCovariantBasisVec(int ind, const VectorFrom &pos);

    virtual VectorTo getUnitVector(int ind, const VectorFrom &pos);

    MatrixNM<Real, N, N> jacobian(const VectorN<Real, N> &x);

    VectorTo   transfVecContravariant(const VectorFrom &vec, const VectorFrom &pos);
    VectorFrom transfInverseVecCovariant(const VectorTo &vec, const VectorFrom &pos);
}; 

template<typename VectorFrom, typename VectorTo, int N>
class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
                               public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
{
public:
    virtual VectorFrom getContravariantBasisVec(int ind, const VectorTo &pos);

    virtual VectorFrom getUnitVectorInverse(int ind, const VectorTo &pos);

    VectorTo   transfVecCovariant(const VectorFrom &vec, const VectorTo &pos);
    VectorFrom transfInverseVecContravariant(const VectorTo &vec, const VectorTo &pos);

    Tensor2<N> transfTensor2(const Tensor2<N> &tensor, const VectorFrom &pos);
    Tensor3<N> transfTensor3(const Tensor3<N> &tensor, const VectorFrom &pos);  
    Tensor4<N> transfTensor4(const Tensor4<N> &tensor, const VectorFrom &pos);
    Tensor5<N> transfTensor5(const Tensor5<N> &tensor, const VectorFrom &pos);
}; 
~~~

## Defined coordinate transformations
  - polar <-> cartesian 2d (CoordTransfPolarToCartesian2D)
  - rotation around axis 2d (CoordTransfCart2DRotation)
  - rotation around axis 3d (CoordTransfCart3DRotationXAxis, CoordTransfCart3DRotationYAxis, CoordTransfCart3DRotationZAxis)
  - spherical <-> cartesian 3d
  - cylindrical <-> cartesian 3d

## Defined static members
~~~C++
	static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
	static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;
	static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
	static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
~~~

## Examples

DEMO COORDINATE TRANSFORMATIONS

~~~C++
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
~~~

Covariant and contravariant vector transformations

~~~C++
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
~~~

Covariant and contravariant basis vectors

~~~C++
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
Vector3Cartesian vec_co_1 = CoordTransfSpherToCart.getCovariantBasisVec(0, p1Spher);
Vector3Cartesian vec_co_2 = CoordTransfSpherToCart.getCovariantBasisVec(1, p1Spher);
Vector3Cartesian vec_co_3 = CoordTransfSpherToCart.getCovariantBasisVec(2, p1Spher);
std::cout << "getCovariantBasisVec 1 : " << vec_co_1 << std::endl;
std::cout << "getCovariantBasisVec 2 : " << vec_co_2 << std::endl;
std::cout << "getCovariantBasisVec 3 : " << vec_co_3 << std::endl << std::endl;

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
Vector3Cartesian vec_contra_1 = CoordTransfSpherToCart.getContravariantBasisVec(0, p1);
Vector3Cartesian vec_contra_2 = CoordTransfSpherToCart.getContravariantBasisVec(1, p1);
Vector3Cartesian vec_contra_3 = CoordTransfSpherToCart.getContravariantBasisVec(2, p1);
std::cout << "getContravariantBasisVec 1 : " << vec_contra_1 << std::endl;
std::cout << "getContravariantBasisVec 2 : " << vec_contra_2 << std::endl;
std::cout << "getContravariantBasisVec 3 : " << vec_contra_3 << std::endl << std::endl;

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
std::cout << vec_co_3.ScalarProductCartesian(vec_contra_1) << std::endl;
std::cout << vec_co_3.ScalarProductCartesian(vec_contra_2) << std::endl;
std::cout << vec_co_3.ScalarProductCartesian(vec_contra_3) << std::endl;

std::cout << "\n*****       Using Cartesian to Spherical transformation      *****\n";

std::cout << "\nCOVARIANT basis vectors\n";
Vector3Spherical vec2_co_1 = CoordTransfCartToSpher.getCovariantBasisVec(0, p1);
Vector3Spherical vec2_co_2 = CoordTransfCartToSpher.getCovariantBasisVec(1, p1);
Vector3Spherical vec2_co_3 = CoordTransfCartToSpher.getCovariantBasisVec(2, p1);
std::cout << "getCovariantBasisVec 1 : " << vec2_co_1 << std::endl;
std::cout << "getCovariantBasisVec 2 : " << vec2_co_2 << std::endl;
std::cout << "getCovariantBasisVec 3 : " << vec2_co_3 << std::endl << std::endl;

std::cout << "GetAsUnitVector 1      : " << vec2_co_1.GetAsUnitVectorAtPos(p1Spher) << std::endl;
std::cout << "GetAsUnitVector 2      : " << vec2_co_2.GetAsUnitVectorAtPos(p1Spher) << std::endl;
std::cout << "GetAsUnitVector 3      : " << vec2_co_3.GetAsUnitVectorAtPos(p1Spher) << std::endl << std::endl;

std::cout << "CONTRAVARIANT basis vectors\n";
Vector3Spherical vec2_contra_1 = CoordTransfCartToSpher.getContravariantBasisVec(0, p1Spher);
Vector3Spherical vec2_contra_2 = CoordTransfCartToSpher.getContravariantBasisVec(1, p1Spher);
Vector3Spherical vec2_contra_3 = CoordTransfCartToSpher.getContravariantBasisVec(2, p1Spher);
std::cout << "getContravariantBasisVec 1 : " << vec2_contra_1 << std::endl;
std::cout << "getContravariantBasisVec 2 : " << vec2_contra_2 << std::endl;
std::cout << "getContravariantBasisVec 3 : " << vec2_contra_3 << std::endl << std::endl;

std::cout << "GetAsUnitVectorAtPos 1     : " << vec2_contra_1.GetAsUnitVectorAtPos(p1) << std::endl;
std::cout << "GetAsUnitVectorAtPos 2     : " << vec2_contra_2.GetAsUnitVectorAtPos(p1) << std::endl;
std::cout << "GetAsUnitVectorAtPos 3     : " << vec2_contra_3.GetAsUnitVectorAtPos(p1) << std::endl << std::endl;

/* OUTPUT
Point (cartesian)       : [              2,               2,               5]
Point (spherical - rad) : [    5.744562647,    0.5148059551,    0.7853981634]
Point (spherical - deg) : [    5.7445626465, 29.4962084966, 45.0000000000 ]

*****       Using Spherical to Cartesian transformation      *****

COVARIANT basis vectors
getCovariantBasisVec 1 : [   0.3481553119,    0.3481553119,    0.8703882798]
getCovariantBasisVec 2 : [   3.5355339059,    3.5355339059,   -2.8284271247]
getCovariantBasisVec 3 : [  -2.0000000000,    2.0000000000,    0.0000000000]

GetAsUnitVector 1      : [   0.3481553119,    0.3481553119,    0.8703882798]
GetAsUnitVector 2      : [   0.6154574549,    0.6154574549,   -0.4923659639]
GetAsUnitVector 3      : [  -0.7071067812,    0.7071067812,    0.0000000000]

CONTRAVARIANT basis vectors
getContravariantBasisVec 1 : [   0.3481553119,    0.3481553119,    0.8703882798]
getContravariantBasisVec 2 : [   0.1071373911,    0.1071373911,   -0.0857099129]
getContravariantBasisVec 3 : [  -0.2500000000,    0.2500000000,    0.0000000000]

GetAsUnitVectorAtPos 1     : [   0.3481553119,    0.3481553119,    0.8703882798]
GetAsUnitVectorAtPos 2     : [   0.6154574549,    0.6154574549,   -0.4923659639]
GetAsUnitVectorAtPos 3     : [  -0.7071067812,    0.7071067812,    0.0000000000]

Veifying orthogonality of covariant and contravariant basis vectors:
-0.0000000000
-0.0000000000
1.0000000000

*****       Using Cartesian to Spherical transformation      *****

COVARIANT basis vectors
getCovariantBasisVec 1 : [   0.3481553119,    0.1071373911,   -0.2500000000]
getCovariantBasisVec 2 : [   0.3481553119,    0.1071373911,    0.2500000000]
getCovariantBasisVec 3 : [   0.8703882798,   -0.0857099129,    0.0000000000]

GetAsUnitVector 1      : [   0.3481553119,    0.0186502259,   -0.0883883476]
GetAsUnitVector 2      : [   0.3481553119,    0.0186502259,    0.0883883476]
GetAsUnitVector 3      : [   0.8703882798,   -0.0149201807,    0.0000000000]

CONTRAVARIANT basis vectors
getContravariantBasisVec 1 : [   0.3481553119,    3.5355339059,   -2.0000000000]
getContravariantBasisVec 2 : [   0.3481553119,    3.5355339059,    2.0000000000]
getContravariantBasisVec 3 : [   0.8703882798,   -2.8284271247,    0.0000000000]

GetAsUnitVectorAtPos 1     : [   0.3481553119,    1.7677669530,   -1.0997501703]
GetAsUnitVectorAtPos 2     : [   0.3481553119,    1.7677669530,    1.0997501703]
GetAsUnitVectorAtPos 3     : [   0.8703882798,   -1.4142135624,    0.0000000000]
*/
~~~

