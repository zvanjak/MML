#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/CoreUtils.h"

#endif

using namespace MML;

void Readme_vectors_matrices()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - vectors, matrices                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Vector<Real>    vec1({ 1.5, -2.0, 0.5 }), vec2({ 1.0, 1.0, -3.0 }); 
    VectorComplex   vec_cmplx1({ Complex(1,1), Complex(-1,2) });
    VectorComplex   vec_cmplx2({ Complex(1,1), Complex(-1,2), Complex(2.5, -1.5) });

    Matrix<Real>    mat_3x3{ 3, 3, { 1.0,  2.0, -1.0, 
                                    -1.0,  5.0,  6.0, 
                                     3.0, -2.0,  1.0 }};  
    MatrixComplex   mat_cmplx(2,2, { Complex(0.5,1), Complex(-1,2), 
                                     Complex(-1,-2), Complex(-2,2) });
    MatrixComplex   mat_cmplx2(2,3, { Complex(1,2),    Complex(-1,1), Complex(1.5,-2), 
                                      Complex(2,-0.5), Complex(3,-2), Complex(-1,1) });
    Matrix<Real>  unit_mat3 = MML::Matrix<Real>::GetUnitMatrix(3);

    Vector<Real> v_real  = Real{2.0} * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
    VectorComplex  v_cmplx = vec_cmplx1 * mat_cmplx / Complex(1.5, -1.5) / 2.0;

    std::cout << "v_real  = " << v_real << std::endl;
    std::cout << "v_cmplx = " << v_cmplx << std::endl;

    // combining real and complex vectors and matrices requires special functions
    VectorComplex cvec2 = MatrixUtils::MulVecMat(vec_cmplx2, mat_3x3);
    MatrixComplex cmat2 = MatrixUtils::MulMat(mat_cmplx2, mat_3x3);

    std::cout << "Matrix mat_3x3   = " << mat_3x3 << std::endl;
    std::cout << "Matrix mat_cmplx = " << mat_cmplx << std::endl;

    std::cout << "IsOrthogonal(mat_3x3)  = " << MatrixUtils::IsOrthogonal(mat_3x3) << std::endl;
    std::cout << "IsHermitian(mat_cmplx) = " << MatrixUtils::IsHermitian(mat_cmplx) << std::endl;
    std::cout << "IsUnitary(mat_cmplx)   = " << MatrixUtils::IsUnitary(mat_cmplx) << std::endl;

/* OUTPUT
    v_real  = [   -3.137858162,     3.922322703,    -8.629109946]
    v_cmplx = [        (0.5,1), (-0,-1.666666667)]
    Matrix mat_3x3   = Rows: 3 Cols: 3
    [          1,          2,         -1,  ]
    [         -1,          5,          6,  ]
    [          3,         -2,          1,  ]

    Matrix mat_cmplx = Rows: 2 Cols: 2
    [    (0.5,1),     (-1,2),  ]
    [    (-1,-2),     (-2,2),  ]

    IsOrthogonal(mat_3x3)  = 0
    IsHermitian(mat_cmplx) = 1
    IsUnitary(mat_cmplx)   = 0
*/
}
