#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;


void Example1_Vector()
{
  std::cout << std::endl;
  std::cout << "***********************************************************************" << std::endl;
  std::cout << "****                            VECTOR                             ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

  // REAL vectors
  Vector<Real>    vec1(5);                    // init vector with 5 elements
  Vector<Real>    vec2(3, 3.14159);           // init with constant value
  Vector<Real>    vec3({ 1.5, -2.1, 0.48 });  // init with list of values
  
  // using defined aliases (typedefs)
  Vector<Real>   vec4(vec3);               // init with copy ctor
  Vector<Real> vec5 = vec2;              // init with assignment

  // initializing from std::vector and C/C++ array 
  std::vector<Real> std_vec{ -1.0, 5.0, -2.0 };
  float  arr[5] = { -1.0, 5.0, -2.0, 10.0, 4.0 };

  Vector<Real>    vec6(std_vec);             // init with std::vector<>
  VectorFlt       vec7(5, arr);              // init with C/C++ array

  // COMPLEX vectors
  Vector<Complex>  vec_c1({ 1.0, 2.0, 3.0 });     // init with list of real values
  VecC             vec_c2({ Complex(1,1),
                            Complex(-1,2),
                            Complex(2, -0.5) });  // init with list of complex values

  // vector operations
	Vector<Real> v1 = vec2 / 2.0 - 1.5 * vec6;
  v1 /= vec4.NormL2();

	Vector<Complex> v2 = vec_c2 * Complex(0.5, -1.5) + Utils::AddVec(v1, vec_c1) + Utils::MulVec(Complex(0, 1), v1);
  
  Real    scalar_prod = Utils::ScalarProduct(vec2, vec3);
  Complex scalar_prod_cmplx = Utils::ScalarProduct(vec_c1, vec_c2);

  // I/O
	std::cout << "Real vectors:" << std::endl;
  std::cout << "v1 = " << v1 << std::endl;
  std::cout << "v1 = " << v1.to_string(10, 5) << std::endl;
	std::cout << "v1 = "; v1.Print(std::cout, 7, 5); std::cout << std::endl;

	std::cout << "\nComplex vectors:" << std::endl;
  std::cout << "v2 = " << v2 << std::endl;
  std::cout << "v2 = " << v2.to_string(10, 6) << std::endl;
  std::cout << "v2 = "; v2.Print(std::cout, 10, 3); std::cout << std::endl;
}