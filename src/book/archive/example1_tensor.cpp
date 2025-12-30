#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Tensor.h"
#endif

using namespace MML;

void Example1_Tensor()
{
  std::cout << "***********************************************************************" << std::endl;
  std::cout << "****                           TENSORS                             ****" << std::endl;
  std::cout << "***********************************************************************" << std::endl;

	Tensor2<2> t2(2, 0, { 1.0, 2.0, 
                        3.0, 4.0 });

	Tensor2<3> t2_3(1, 1, { 1.0, 2.0, 2.0,
													3.0, 4.0, 0.5,
													5.0, 6.0, 1.5 });
	
	Tensor2<4> t2_4(2, 0, {-1.0, 0.0, 0.0, 0.0,
													0.0, 1.0, 0.0, 0.0,
													0.0, 0.0, 1.0, 0.0,
													0.0, 0.0, 0.0, 1.0 });

	Tensor3<3> t3(1, 2, { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
												10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
												19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0 });

	std::cout << "t2   = " << t2 << std::endl;
	std::cout << "t2_3 = " << t2_3 << std::endl;
	std::cout << "t2_4 = " << t2_4 << std::endl;

	// arithmetic operations
	Tensor2<3> t2_3_2 = t2_3 + t2_3;
	Tensor2<3> t2_3_3 = t2_3 - t2_3;
	Tensor2<3> t2_3_4 = t2_3 * 2.0;
	Tensor2<3> t2_3_5 = t2_3 / 2.0;
	Tensor2<3> t2_3_6 =  2.0 * t2_3;

	// contraction
	Real val1 = t2_3.Contract();
	std::cout << "t2_3.Contract() = " << val1 << std::endl;

	VectorN<Real, 3> vec_c1 = t3.Contract(0, 1);
	std::cout << "t3.Contract(0, 1) = " << vec_c1 << std::endl;
	VectorN<Real, 3> vec_c2 = t3.Contract(0, 2);
	std::cout << "t3.Contract(0, 2) = " << vec_c2 << std::endl;

	// second rank tensor operating on two vectors
	VectorN<Real, 3> vec1{ 1.0, -2.0, 3.0 }, vec2{ 0.5, 1.0, -1.5 };
	Real val2 = t2_3(vec1, vec2);
	std::cout << "t2_3(vec1, vec2) = " << val2 << std::endl;

	VectorN<Real, 4> vec3{ 1.0, -2.0, 3.0, 1.0 }, vec4{ 0.5, 1.0, -1.5, 1.0 };
	Real val3 = t2_4(vec3, vec4);
	std::cout << "t2_4(vec3, vec4) = " << val3 << std::endl;
}