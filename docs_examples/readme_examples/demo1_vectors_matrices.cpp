#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "core/MatrixUtils.h"

#endif

using namespace MML;

void Readme_vectors_matrices()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - vectors, matrices                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Vector<double>  vec1({ 1.5, -2.0, 0.5 }), vec2({ 1.0, 1.0, -3.0 }); 
    VectorComplex   vec_cmplx({ Complex(1,1), Complex(-1,2) });
    VectorComplex   vec_cmplx2({ Complex(1,1), Complex(-1,2), Complex(2.5, -1.5) });

    Matrix<double>  mat_3x3{ 3, 3, { 1.0, 2.0, -1.0, 
                                   -1.0, 5.0, 6.0, 
                                    3.0, 1.0, 1.0 }};  
    MatrixComplex   mat_cmplx(2,2, { Complex(1,1),  Complex(-1,2), 
                                    Complex(2, -0.5), Complex(1,1) });
    MatrixComplex   mat_cmplx2(2,3, { Complex(1,1),  Complex(-1,2), Complex(1.5, -2), 
                                    Complex(2, -0.5), Complex(1,1), Complex(-1, 1) });
    Matrix<double>  unit_mat3 = MML::Matrix<Real>::GetUnitMatrix(3);

    Vector<double> v = 2.0 * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
    VectorComplex  vc = vec_cmplx * mat_cmplx / Complex(1.5, -1.5) / 2.0;

    // combining real and complex vectors and matrices requires special functions
    VectorComplex vc2 = MatrixUtils::MulVecMat(vec_cmplx2, mat_3x3);
    MatrixComplex mc = MatrixUtils::MulMat(mat_cmplx2, mat_3x3);
}
