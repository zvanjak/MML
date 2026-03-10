///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme02_vectors_matrices.cpp                                       ///
///  Description: README Vectors & Matrices section demo                              ///
///               Demonstrates vector/matrix creation, arithmetic, properties         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;

void Readme_VectorsMatrices()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                README - Vectors & Matrices                    ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Real vectors and matrices
    Vector<Real> vec1{1.5, -2.0, 0.5}, vec2{1.0, 1.0, -3.0};
    Matrix<Real> mat_3x3{3, 3, {1.0,  2.0, -1.0,
                               -1.0,  5.0,  6.0,
                                3.0, -2.0,  1.0}};

    // Vector and matrix arithmetic
    Vector<Real> result = Real{2.0} * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
    std::cout << "Result: " << result << std::endl;

    // Complex vectors and matrices
    VectorComplex vec_cmplx{{Complex(1,1), Complex(-1,2)}};
    MatrixComplex mat_cmplx{2, 2, {Complex(0.5,1),  Complex(-1,2),
                                   Complex(-1,-2), Complex(-2,2)}};

    // Matrix properties
    std::cout << "IsOrthogonal:  " << Utils::IsOrthogonal(mat_3x3) << std::endl;
    std::cout << "IsHermitian:   " << Utils::IsHermitian(mat_cmplx) << std::endl;

/* Expected OUTPUT:
    Result: [   -3.137858162,     3.922322703,    -8.629109946]
    IsOrthogonal:  0
    IsHermitian:   1
*/
}
