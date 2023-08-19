#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#include "core/Matrix.h"
#include "core/Matrix3D.h"
#endif

using namespace MML;

void Matrix_initializations() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******            Matrix Initializations           *******" << std::endl;

    MML::Matrix<Real> a;                                // empty matrix
    MML::Matrix<Real> b(2,2);                           // ampty matrix 2 x 2
    MML::Matrix<Real> c(2,2, {1.0, 0.0, 0.0, 1.0});     // matrix with initialized values
    MML::Matrix<Real> d(c);
    MML::Matrix<Real> e = c;
    MML::Matrix<Real> f = MML::Matrix<Real>::GetUnitMatrix(3);

    MML::MatrixInt      mat_int;
    MML::MatrixDbl      mat_dbl(3,3);
    MML::MatrixComplex  mat_cmplx(2,2, { Complex(1,1),  Complex(-1,2), 
                                         Complex(2, -0.5), Complex(1,1) });

    MML::MatI          mat_i(3,3);
    MML::MatF          mat_d(3,3);
    MML::MatC          mat_c(3,3);

    MML::Matrix3D<float>  mat_3d(3,3,3);

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "f = " << f << std::endl;
}

void Matrix_vector_init_operations() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******        Matrix-vector init operations        *******" << std::endl;

    Vector a({1.0, 1.0, 1.0});
    Matrix<Real> matA = Matrix<Real>::RowMatrixFromVector(a);
    Matrix<Real> matB = Matrix<Real>::ColumnMatrixFromVector(a);

    std::cout << "Vector a = " << a << std::endl;
    std::cout << "Matrix matA = Matrix::RowMatrixFromVector(a);\nmatA = " << matA << std::endl;
    std::cout << "Matrix matB = Matrix::ColMatrixFromVector(a);\nmatB = " << matB << std::endl;

    Matrix<Real> m1(2,2, {1.0, -1.0, 1.5, 3.0});
    Vector vecRow  = Matrix<Real>::VectorFromRow(m1, 0);
    Vector vecCol  = Matrix<Real>::VectorFromColumn(m1, 0);
    Vector vecDiag = Matrix<Real>::VectorFromDiagonal(m1);

    std::cout << "Matrix m1 = " << m1 << std::endl;
    std::cout << "Vector vecRow = Matrix::VectorFromRow(a, 0)    = " << vecRow << std::endl;
    std::cout << "Vector vecCol = Matrix::VectorFromColumn(a, 0) = " << vecCol << std::endl;
    std::cout << "Vector vecCol = Matrix::VectorFromDiagonal(a)  = " << vecDiag << std::endl;
}

void Matrix_accessing_elements() 
{
    // ElemAt + operator[]
}

void Basic_mat_operations() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******         Matrix Basic operations             *******" << std::endl;

    Matrix<Real> m1(2,2, {1.0, -1.0, 1.5, 3.0}), m2(2,2, {-2.5, 3.0, 1, -1.0});
    m2.MakeUnitMatrix();

    std::cout << "m1 = " << m1 << std::endl;
    std::cout << "m2 = " << m2 << std::endl;

    std::cout << "m1 + m2 = " << m1 + m2  << std::endl;
    std::cout << "m1 - m2 = " << m1 - m2 << std::endl;
    std::cout << "m1 * m2 = " << m1 * m2 << std::endl;
    std::cout << "2.0 * m1 = " << 2.0 * m1 << std::endl;
    std::cout << "m1 * 2.0 = " << m1 * 2.0 << std::endl;
    std::cout << "m1 / 2.0 = " << m1 / 2.0 << std::endl;
}

void Matrix_Vector_mul() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******          Matrix Vector multiplication       *******" << std::endl;

    Vector v1({1.0, 2.0});
    Matrix<Real> m1(2,2, {1.0, -1.0, 1.5, 3.0});
    
    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "m1 = " << m1 << std::endl;

    std::cout << "v1 * m1 = " << v1 * m1 << std::endl;
    std::cout << "m1 * v1 = " << m1 * v1 << std::endl;
}

void Matrix_Matrix_mul() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******          Matrix Matrix multiplication       *******" << std::endl;

    MML::Matrix<Real> m3(1, 3, {1.0, 1.0, 1.0});
    MML::Matrix<Real> m4(3, 4, {1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0, 
                            0.0, 0.0, 1.0, 1.0});
    MML::Matrix<Real> m5 = m3 * m4;

    std::cout << "m3 = " << m3 << std::endl;
    std::cout << "m4 = " << m4 << std::endl;
    std::cout << "m3 * m4 = " << m5  << std::endl;
}

void Matrix_invert()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******          Matrix Invert                      *******" << std::endl;

    MML::Matrix<Real> m1(2, 2, {1.0, -1.0, 1.5, 3.0});

    std::cout << "m1       = " << m1 << std::endl;
    auto m2 = m1.GetInverse();
    
    std::cout << "m2 (inv) = " << m2 << std::endl;

    auto munit = m1 * m2;
    std::cout << "m1 * m2 = " << munit << std::endl;
}

void Matrix_transpose()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******          Matrix Transpose                   *******" << std::endl;

    MML::Matrix<Real> m1(2, 2, {1.0, -1.0, 1.5, 3.0});
    std::cout << "m1          = " << m1 << std::endl;

    auto m2 = m1.GetTranspose();
    std::cout << "m2 (transp) = " << m2 << std::endl;

    m1.Transpose();
    std::cout << "Transposing m1 (in place) = " << m1 << std::endl;
}

void Demo_Matrix() 
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                              MATRIX                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Matrix_initializations();
    Matrix_vector_init_operations();
    Matrix_accessing_elements();

    Basic_mat_operations();
    Matrix_Vector_mul();
    Matrix_Matrix_mul();

    Matrix_invert();
    Matrix_transpose();
}