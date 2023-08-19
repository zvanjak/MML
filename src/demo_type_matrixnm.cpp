#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/MatrixNM.h"
#endif

using namespace MML;

void MatrixNM_initializations() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******             MatrixNM initialization         *******" << std::endl;

    MML::MatrixNM<Real,2,2> a;
    MML::MatrixNM<Real,2,2> b({1.0, 0.0, 0.0, 1.0});
    MML::MatrixNM<Real,2,2> c(b);
    MML::MatrixNM<Real,2,2> d = c;
    auto e = MML::MatrixNM<Real,3,3>::GetUnitMatrix();

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "e = " << e << std::endl;
}

void MatrixNM_vector_init_operations() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******        MatrixNM - VectorN init operations    ******" << std::endl;

    MML::VectorN<Real, 3> a({1.0, 1.0, 1.0});
    MML::MatrixNM<Real,1,3> matA = MML::MatrixNM<Real,1,3>::RowMatrixFromVector(a);
    auto matAauto = MML::MatrixNM<Real,1,3>::RowMatrixFromVector(a);
    MML::MatrixNM<Real,3,1> matB = MML::MatrixNM<Real,3,1>::ColumnMatrixFromVector(a);
    auto matBauto = MML::MatrixNM<Real,3,1>::ColumnMatrixFromVector(a);

    std::cout << "Vector a = " << a << std::endl;
    std::cout << "Matrix matA = Matrix::RowMatrixFromVector(a);\nmatA = " << matA << std::endl;
    std::cout << "Matrix matB = Matrix::ColMatrixFromVector(a);\nmatB = " << matB << std::endl;

    MML::MatrixNM<Real,2,2> m1({1.0, -1.0, 1.5, 3.0});
    MML::VectorN<Real, 2> vecRow = MML::MatrixNM<Real,2,2>::VectorFromRow(m1, 0);
    MML::VectorN<Real, 2> vecCol = MML::MatrixNM<Real,2,2>::VectorFromColumn(m1, 0);
    MML::VectorN<Real, 2> vecDiag = MML::MatrixNM<Real,2,2>::VectorFromDiagonal(m1);

    std::cout << "Matrix m1 = " << m1 << std::endl;
    std::cout << "Vector vecRow = Matrix::VectorFromRow(a,0)     = " << vecRow << std::endl;
    std::cout << "Vector vecCol = Matrix::VectorFromColumn(a, 0) = " << vecCol << std::endl;
    std::cout << "Vector vecCol = Matrix::VectorFromDiagonal(a)  = " << vecDiag << std::endl;
}

void MatrixNM_accessing_elements() 
{
    // ElemAt + operator[]
}

void Basic_MatrixNM_operations() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******             MatrixNM Basic operation        *******" << std::endl;

    MML::MatrixNM<Real,2,2> m1({1.0, -1.0, 1.5, 3.0}), m2;
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

void MatrixNM_VectorN_mul() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******       MatrixNM VectorN multiplication       *******" << std::endl;

    MML::VectorN<Real, 2> v1({1.0, 2.0});
    MML::MatrixNM<Real,2,2> m1({1.0, -1.0, 1.5, 3.0});
    
    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "m1 = " << m1 << std::endl;

    std::cout << "v1 * m1 = " << v1 * m1  << std::endl;
    std::cout << "m1 * v1 = " << m1 * v1 << std::endl;
}

void MatrixNM_MatrixNM_mul() 
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******        MatrixNM MatrixNM multiplication     *******" << std::endl;

    MML::MatrixNM<Real,1,3> m3{1.0, 1.0, 1.0};
    MML::MatrixNM<Real,3,4> m4{1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0, 
                            0.0, 0.0, 1.0, 1.0};
    MML::MatrixNM<Real,1,4> m5  = m3 * m4;

    std::cout << "m3 = " << m3 << std::endl;
    std::cout << "m4 = " << m4 << std::endl;
    std::cout << "m3 * m4 = " << m5  << std::endl;
}

void MatrixNM_Invert()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******             MatrixNM Invert                 *******" << std::endl;

    MML::MatrixNM<Real,2,2> m1({1.0, -1.0, 1.5, 3.0});

    std::cout << "m1 = " << m1 << std::endl;
    auto m2 = m1.GetInverse();
    
    std::cout << "m2 (inv) = " << m2 << std::endl;

    auto munit = m1 * m2;
    std::cout << "m1 * m2 = " << munit << std::endl;
}

void MatrixNM_transpose()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*******             MatrixNM Transpose              *******" << std::endl;

    MML::MatrixNM<Real,2,2> m1({1.0, -1.0, 1.5, 3.0});
    std::cout << "m1 = " << m1 << std::endl;

    auto m2 = m1.GetTranspose();
    std::cout << "m2 (transp) = " << m2 << std::endl;

    m1.Transpose();
    std::cout << "Transposing m1 (in place) = " << m1 << std::endl;
}

void Demo_MatrixNM()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         MATRIX N_M                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
 
    MatrixNM_initializations();
    MatrixNM_vector_init_operations();
    MatrixNM_accessing_elements();
 
    Basic_MatrixNM_operations();
    MatrixNM_VectorN_mul();
    MatrixNM_MatrixNM_mul();

    MatrixNM_Invert();
    MatrixNM_transpose();
}