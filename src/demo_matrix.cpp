#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#endif

void Basic_mat_operations() 
{
    MML::Matrix m1(2,2, {1.0, -1.0, 1.5, 3.0}), m2(2,2);
    m2.MakeUnitMatrix();

    std::cout << "m1 = " << m1 << std::endl;
    std::cout << "m2 = " << m2 << std::endl;

    std::cout << "m1 + m2 = " << m1 + m2  << std::endl;
    std::cout << "m1 - m2 = " << m1 - m2 << std::endl;
    std::cout << "m1 * m2 = " << m1 * m2 << std::endl;
}

void Matrix_vector_operations() 
{
    MML::Vector a({1.0, 1.0, 1.0});
    MML::Matrix matA = MML::Matrix::RowMatrixFromVector(a);
    MML::Matrix matB = MML::Matrix::ColumnMatrixFromVector(a);

    std::cout << "Vector a = " << a << std::endl;
    std::cout << "Matrix matA = Matrix::RowMatrixFromVector(a) = " << matA << std::endl;
    std::cout << "Matrix matB = Matrix::ColMatrixFromVector(a) = " << matB << std::endl;

    MML::Matrix m1(2,2, {1.0, -1.0, 1.5, 3.0});
    MML::Vector vecRow = MML::Matrix::VectorFromRow(m1, 0);
    MML::Vector vecCol = MML::Matrix::VectorFromColumn(m1, 0);

    std::cout << "Matrix m1 = " << m1 << std::endl;
    std::cout << "Vector vecRow = Matrix::VectorFromRow(a,0) = " << vecRow << std::endl;
    std::cout << "Vector vecCol = Matrix::VectorFromColumn(a, 0) = " << vecCol << std::endl;
}

void Demo_Matrix() 
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                              MATRIX                           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::Matrix a;
    MML::Matrix b(2,2);
    MML::Matrix c(2,2, {1.0, 0.0, 0.0, 1.0});
    MML::Matrix d(c);
    MML::Matrix e = c;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "e = " << e << std::endl;

    Basic_mat_operations();
    Matrix_vector_operations();
}