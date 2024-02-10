#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/BaseUtils.h"

#include "core/CoreUtils.h"

#include "algorithms/EigenSystemSolvers.h"

#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

//////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_Symmetric_Matrix_Eigen_solver 3 x 3", "[simple]")
{
	MatrixSym<Real>       mat = TestBeds::symm_mat_3x3;

    SymmMatEigenSolver eigen_solver(mat);

    for(int i=0; i<3; i++)
    {
        Real            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Real>    eigen_vector = eigen_solver.getEigenvector(i);
        
        REQUIRE(eigen_value == Approx(TestBeds::symm_mat_3x3_eigen_val[i]).epsilon(1e-14));
        REQUIRE_THAT(eigen_value, Catch::Matchers::WithinAbs(TestBeds::symm_mat_3x3_eigen_val[i],1e-14));

        REQUIRE(true == eigen_vector.IsEqual(TestBeds::symm_mat_3x3_eigen_vecs[i], 1e-10));

        REQUIRE(true == Vector<Real>::AreEqual(mat * eigen_vector, eigen_value * eigen_vector, 1e-14));
        REQUIRE(false == Vector<Real>::AreEqual(mat * eigen_vector, eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Symmetric_Matrix_Eigen_solver 5 x 5", "[simple]")
{
	MatrixSym<Real>       mat = TestBeds::symm_mat_5x5;

    SymmMatEigenSolver eigen_solver(mat);

    for(int i=0; i<5; i++)
    {
        Real            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Real>    eigen_vector = eigen_solver.getEigenvector(i);
        
        REQUIRE(eigen_value == Approx(TestBeds::symm_mat_5x5_eigen_val[i]).epsilon(1e-14));
        REQUIRE(true == eigen_vector.IsEqual(TestBeds::symm_mat_5x5_eigen_vecs[i], 1e-10));

        REQUIRE(true == Vector<Real>::AreEqual(mat * eigen_vector, eigen_value * eigen_vector, 1e-14));
        REQUIRE(false == Vector<Real>::AreEqual(mat * eigen_vector, eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Symmetric_Matrix_Eigen_solver 10 x 10", "[simple]")
{
	MatrixSym<Real>       mat = TestBeds::symm_mat_10x10;

    SymmMatEigenSolver eigen_solver(mat);

    for(int i=0; i<10; i++)
    {
        Real            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Real>    eigen_vector = eigen_solver.getEigenvector(i);
        
        REQUIRE(eigen_value == Approx(TestBeds::symm_mat_10x10_eigen_val[i]).epsilon(1e-14));
        REQUIRE(true == eigen_vector.IsEqual(TestBeds::symm_mat_10x10_eigen_vecs[i], 1e-10));

        REQUIRE(true == Vector<Real>::AreEqual(mat * eigen_vector, eigen_value * eigen_vector, 1e-14));
        REQUIRE(false == Vector<Real>::AreEqual(mat * eigen_vector, eigen_value * eigen_vector, 1e-16));
    }  
}

//////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 3 x 3", "[simple]")
{
	Matrix<Real>       mat = TestBeds::mat_3x3;

    EigenSolver eigen_solver(mat);

    for(int i=0; i<3; i++) 
    {
        Complex            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Complex>    eigen_vector = eigen_solver.getEigenvector(i);
       
        REQUIRE(true == Utils::AreEqual(eigen_value, TestBeds::mat_3x3_eigen_val[i], 1e-9));
        REQUIRE(true == Utils::AreEqual(eigen_vector, TestBeds::mat_3x3_eigen_vecs[i], 1e-9));

        REQUIRE(true == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-14));
        REQUIRE(false  == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 5 x 5", "[simple]")
{
	Matrix<Real>       mat = TestBeds::mat_5x5;

    EigenSolver eigen_solver(mat);

    for(int i=0; i<5; i++)               
    {
        Complex            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Complex>    eigen_vector = eigen_solver.getEigenvector(i);
       
        REQUIRE(true == Utils::AreEqual(eigen_value, TestBeds::mat_5x5_eigen_val[i], 1e-9));
        REQUIRE(true == Utils::AreEqual(eigen_vector, TestBeds::mat_5x5_eigen_vecs[i], 1e-9));

        REQUIRE(true == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-13));
        REQUIRE(false  == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 8 x 8", "[simple]")
{
	Matrix<Real>       mat = TestBeds::mat_8x8;

    EigenSolver eigen_solver(mat);

    for(int i=0; i<5; i++)               
    {
        Complex            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Complex>    eigen_vector = eigen_solver.getEigenvector(i);
       
        REQUIRE(true == Utils::AreEqual(eigen_value, TestBeds::mat_8x8_eigen_val[i], 1e-9));
        REQUIRE(true == Utils::AreEqual(eigen_vector, TestBeds::mat_8x8_eigen_vecs[i], 1e-9));

        REQUIRE(true == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-13));
        REQUIRE(false  == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 10 x 10", "[simple]")
{
	Matrix<Real>       mat = TestBeds::mat_10x10;

    EigenSolver eigen_solver(mat);

    for(int i=0; i<10; i++)               
    {
        Complex            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Complex>    eigen_vector = eigen_solver.getEigenvector(i);
       
        REQUIRE(true == Utils::AreEqual(eigen_value, TestBeds::mat_10x10_eigen_val[i], 1e-9));
        REQUIRE(true == Utils::AreEqual(eigen_vector, TestBeds::mat_10x10_eigen_vecs[i], 1e-9));

        REQUIRE(true == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-13));
        REQUIRE(false  == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 20 x 20", "[simple]")
{
	Matrix<Real>       mat = TestBeds::mat_20x20;

    EigenSolver eigen_solver(mat);

    for(int i=0; i<20; i++)
    {
        Complex            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Complex>    eigen_vector = eigen_solver.getEigenvector(i);
       
        REQUIRE(true == Utils::AreEqual(eigen_value, TestBeds::mat_20x20_eigen_val[i], 1e-9));
        REQUIRE(true == Utils::AreEqual(eigen_vector, TestBeds::mat_20x20_eigen_vecs[i], 1e-9));

        REQUIRE(true == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-13));
        REQUIRE(false  == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-15));
    }  
}

TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 50 x 50", "[simple]")
{
	Matrix<Real>       mat = TestBeds::mat_50x50;

    EigenSolver eigen_solver(mat);

    for(int i=0; i<50; i++)
    {
        Complex            eigen_value  = eigen_solver.getEigenvalues()[i];
        Vector<Complex>    eigen_vector = eigen_solver.getEigenvector(i);
       
        REQUIRE(true == Utils::AreEqual(eigen_value, TestBeds::mat_50x50_eigen_val[i], 1e-9));
        REQUIRE(true == Utils::AreEqual(eigen_vector, TestBeds::mat_50x50_eigen_vecs[i], 1e-9));

        REQUIRE(true == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-12));
        REQUIRE(false  == Utils::AreEqual(MatrixUtils::MulMatVec(mat, eigen_vector), eigen_value * eigen_vector, 1e-15));
    }  
}