#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"
#endif

using namespace MML;

void Readme_vectors_matrices()
{
    Vector<double>  vec_dbl_3({ 1.0, 2.0, 3.0 }); 
    VectorComplex   vec_cmplx_2({ Complex(1,1), Complex(-1,2), Complex(2, -0.5) });

    Matrix<Real>   mat_3x3{ 3, 3, { 1.0, 2.0, -1.0, 
                                   -1.0, 5.0, 6.0, 
                                    3.0, 1.0, 1.0 }};  
    MatrixComplex  mat_cmplx(2,2, { Complex(1,1),  Complex(-1,2), 
                                    Complex(2, -0.5), Complex(1,1) });
    Matrix<Real>   unit_mat3 = MML::Matrix<Real>::GetUnitMatrix(3);

}
