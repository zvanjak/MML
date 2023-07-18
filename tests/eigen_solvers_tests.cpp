#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#include "algorithms/EigenSystemSolvers.h"
#endif

#include "../test_data/linear_alg_eq_test_bed.h"


TEST_CASE("Test_Symmetric_Matrix_Eigen_solver 4 x 4", "[simple]")
{
	MML::Matrix<Real>       mat = MML::Tests::LinearAlgEqTestBed::symm_mat4;

    MML::SymmMatEigenSolver eigen_solver(mat);

    MML::Vector<Real>       eigen_values = eigen_solver.d;
    MML::Vector<Real>       eigen_vector1 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 0);
    MML::Vector<Real>       eigen_vector2 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 1);
    MML::Vector<Real>       eigen_vector3 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 2);
    MML::Vector<Real>       eigen_vector4 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 3);

	REQUIRE(true == MML::Vector<Real>::AreEqual(mat * eigen_vector1, eigen_values[0] * eigen_vector1, 1e-14));  
	REQUIRE(false == MML::Vector<Real>::AreEqual(mat * eigen_vector1, eigen_values[0] * eigen_vector1, 1e-15));  
}


TEST_CASE("Test_Unsymmetric_Matrix_Eigen_solver 5 x 5", "[simple]")
{
	MML::Matrix<Real>       mat = MML::Tests::LinearAlgEqTestBed::mat4;

    MML::Unsymmeig eigen_solver(mat);

    MML::Vector<Complex>    eigen_values = eigen_solver.wri;
    MML::Vector<Real>       eigen_vectors1 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 0);
    MML::Vector<Real>       eigen_vectors2 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 1);
    MML::Vector<Real>       eigen_vectors3 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 2);
    MML::Vector<Real>       eigen_vectors4 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 3);
    MML::Vector<Real>       eigen_vectors5 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 4);

	REQUIRE(true == MML::Vector<Real>::AreEqual(mat * eigen_vectors1, eigen_values[0].real() * eigen_vectors1, 1e-13));  
	REQUIRE(false == MML::Vector<Real>::AreEqual(mat * eigen_vectors1, eigen_values[0].real() * eigen_vectors1, 1e-14));  

	REQUIRE(true == MML::Vector<Real>::AreEqual(mat * eigen_vectors2, eigen_values[1].real() * eigen_vectors2, 1e-14));  
	REQUIRE(false == MML::Vector<Real>::AreEqual(mat * eigen_vectors2, eigen_values[1].real() * eigen_vectors2, 1e-15));  

	REQUIRE(true == MML::Vector<Real>::AreEqual(mat * eigen_vectors3, eigen_values[2].real() * eigen_vectors3, 1e-14));  
	REQUIRE(false == MML::Vector<Real>::AreEqual(mat * eigen_vectors3, eigen_values[2].real() * eigen_vectors3, 1e-15));  

}