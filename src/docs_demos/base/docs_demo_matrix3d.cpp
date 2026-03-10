///////////////////////////////////////////////////////////////////////////////////////////
///  File:        docs_demo_matrix3d.cpp                                              ///
///  Description: Brief demonstration of Matrix3D.h - 3D tensor class                 ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Matrix/Matrix3D.h"
#endif

#include <iostream>

using namespace MML;

void Docs_Demo_Matrix3D()
{
    std::cout << "=== Matrix3D - 3D Tensor Demo ===" << std::endl;
    
    // Create 3x4x5 tensor initialized to zero
    Matrix3D<Real> tensor(3, 4, 5, 0.0);
    std::cout << "Created 3x4x5 tensor: " << tensor.Dim1() << "x" << tensor.Dim2() << "x" << tensor.Dim3() << std::endl;
    std::cout << "Total elements: " << tensor.Size() << std::endl;
    
    // Element access - both notations
    tensor[0][0][0] = 1.0;     // Bracket notation
    tensor(1, 2, 3) = 2.5;     // Parenthesis notation
    std::cout << "tensor[0][0][0] = " << tensor[0][0][0] << std::endl;
    std::cout << "tensor(1,2,3) = " << tensor(1, 2, 3) << std::endl;
    
    // Fill and scalar operations
    tensor.Fill(1.0);
    tensor *= 2.0;
    std::cout << "After Fill(1.0) and *=2.0: tensor[1][1][1] = " << tensor[1][1][1] << std::endl;
    
    // Element-wise operations
    Matrix3D<Real> other(3, 4, 5, 1.0);
    Matrix3D<Real> sum = tensor + other;
    std::cout << "After adding tensor with ones: sum[1][1][1] = " << sum[1][1][1] << std::endl;
    
    // Apply function
    sum.Apply([](Real x) { return x * x; });
    std::cout << "After Apply(x*x): sum[1][1][1] = " << sum[1][1][1] << std::endl;
    
    std::cout << "=== Demo Complete ===" << std::endl;
}
