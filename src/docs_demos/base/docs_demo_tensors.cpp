#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Tensor.h"
#include "base/VectorN.h"
#endif

using namespace MML;

void Docs_Demo_Tensors()
{
	std::cout << "***********************************************************" << std::endl;
	std::cout << "*****                       Tensors                 *******" << std::endl;

	// Tensor2 construction with different index variance
	std::cout << "\n*** Tensor2 Construction ***" << std::endl;
	
	Tensor2<3> T_mixed(1, 1);   // 1 covariant, 1 contravariant (mixed tensor)
	Tensor2<3> T_contra(0, 2);  // 2 contravariant indices
	Tensor2<3> T_covar(2, 0);   // 2 covariant indices
	
	std::cout << "Mixed tensor (1 covar, 1 contra): NumCovar=" << T_mixed.NumCovar() 
	          << ", NumContravar=" << T_mixed.NumContravar() << std::endl;

	// With initial values - 2D metric tensor
	Tensor2<2> metric(2, 0, {
		1, 0,
		0, 1
	});
	std::cout << "\n2D Euclidean metric tensor:" << std::endl;
	std::cout << metric << std::endl;

	// Index access
	std::cout << "*** Tensor2 Index Access ***" << std::endl;
	Tensor2<3> T(1, 1);
	T(0, 0) = 1.0;
	T(1, 1) = 2.0;
	T(2, 2) = 3.0;
	T(0, 1) = 0.5;
	
	std::cout << "T(0,0) = " << T(0, 0) << std::endl;
	std::cout << "T(1,1) = " << T(1, 1) << std::endl;
	std::cout << "T.at(0,1) = " << T.at(0, 1) << std::endl;
	
	// Index variance queries
	std::cout << "Is index 0 contravariant? " << (T.IsContravar(0) ? "yes" : "no") << std::endl;
	std::cout << "Is index 1 covariant? " << (T.IsCovar(1) ? "yes" : "no") << std::endl;

	// Tensor operations
	std::cout << "\n*** Tensor2 Operations ***" << std::endl;
	Tensor2<2> A(1, 1, {2, 0, 0, 3});
	Tensor2<2> B(1, 1, {1, 0, 0, 1});
	
	auto sum = A + B;
	auto scaled = A * 2.0;
	
	std::cout << "A = " << A << std::endl;
	std::cout << "B = " << B << std::endl;
	std::cout << "A + B = " << sum << std::endl;
	std::cout << "A * 2.0 = " << scaled << std::endl;

	// Contraction (trace for mixed tensors)
	std::cout << "\n*** Contraction ***" << std::endl;
	Real trace_A = A.Contract();
	std::cout << "A.Contract() = " << trace_A << " (should be 2+3=5)" << std::endl;

	// Tensor evaluation on vectors
	std::cout << "\n*** Tensor Evaluation ***" << std::endl;
	VectorN<Real, 2> v1{1, 0};
	VectorN<Real, 2> v2{0, 1};
	Real result = A(v1, v2);
	std::cout << "A(v1, v2) where v1={1,0}, v2={0,1} = " << result << std::endl;

	// Tensor3 - Christoffel symbols example
	std::cout << "\n*** Tensor3 (Rank-3) ***" << std::endl;
	Tensor3<3> christoffel(2, 1);  // 2 covariant, 1 contravariant
	christoffel(0, 1, 1) = 0.5;
	christoffel(1, 0, 1) = 0.5;
	std::cout << "Christoffel symbol Gamma^0_11 = " << christoffel(0, 1, 1) << std::endl;
	std::cout << "NumCovar=" << christoffel.NumCovar() << ", NumContravar=" << christoffel.NumContravar() << std::endl;
}