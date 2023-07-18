#include <iostream>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#include "algorithms/EigenSystemSolvers.h"
#endif
#include "../test_data/linear_alg_eq_test_bed.h"

void Test_Symmetric_Eigen_Solver()
{
    std::cout << "CALCULATING EIGENVALUES OF SYMMETRIC MATRIX:\n";

    MML::Matrix<Real>     origMat = MML::Tests::LinearAlgEqTestBed::symm_mat4;
    MML::Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);

    MML::SymmMatEigenSolver eigen_solver(matcopy, true);

    std::cout << "Eigenvalues " << eigen_solver.d << std::endl;
    std::cout << "Eigenvectors " << eigen_solver.z << std::endl;

    MML::Vector<Real>       eigen_values = eigen_solver.d;
    MML::Vector<Real>       eigen_vectors1 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 0);
    MML::Vector<Real>       eigen_vectors2 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 1);
    MML::Vector<Real>       eigen_vectors3 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 2);
    MML::Vector<Real>       eigen_vectors4 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 3);
    MML::Vector<Real>       eigen_vectors5 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.z, 4);

    std::cout << "Mat * eig_vec1      = " << matcopy * eigen_vectors1 << std::endl;
    std::cout << "eig_val1 * eig_vec1 = " << eigen_values[0] * eigen_vectors1 << std::endl;

    std::cout << "Mat * eig_vec2      = " << matcopy * eigen_vectors2 << std::endl;
    std::cout << "eig_val2 * eig_vec2 = " << eigen_values[1] * eigen_vectors2 << std::endl;

    std::cout << "Mat * eig_vec3      = " << matcopy * eigen_vectors3 << std::endl;
    std::cout << "eig_val3 * eig_vec3 = " << eigen_values[2] * eigen_vectors3 << std::endl;

	std::cout << "Mat * eig_vec4      = " << matcopy * eigen_vectors4 << std::endl;
    std::cout << "eig_val4 * eig_vec4 = " << eigen_values[3] * eigen_vectors4 << std::endl;

	std::cout << "Mat * eig_vec5      = " << matcopy * eigen_vectors5 << std::endl;
    std::cout << "eig_val5 * eig_vec5 = " << eigen_values[4] * eigen_vectors5 << std::endl;
}

void Test_Unsymmetric_Eigen_Solver()
{
    std::cout << "CALCULATING EIGENVALUES OF UNSYMMETRIC MATRIX:\n";

    MML::Matrix<Real>     origMat = MML::Tests::LinearAlgEqTestBed::mat4;
    MML::Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);

    MML::Unsymmeig eigen_solver(matcopy, true, false);

    std::cout << "Eigenvalues C " << eigen_solver.wri << std::endl;
    std::cout << "Exact eigenvalues = 23.8652 + 0. I, -3.42447 + 1.92848 I, -3.42447 - 1.92848 I,  3.8746 + 0. I, -0.290852 + 0. I\n";
    std::cout << "Eigenvectors " << eigen_solver.zz << std::endl;

    MML::Vector<Complex>       eigen_values = eigen_solver.wri;
    MML::Vector<Real>          eigen_vectors1 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 0);
    MML::Vector<Real>       eigen_vectors2 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 1);
    MML::Vector<Real>       eigen_vectors3 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 2);
    //MML::Vector<Real>       eigen_vectors4 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 3);
    //MML::Vector<Real>       eigen_vectors5 = MML::Matrix<Real>::VectorFromColumn(eigen_solver.zz, 4);

	std::cout << "Mat * eig_vec1      = " << matcopy * eigen_vectors1 << std::endl;
    std::cout << "eig_val1 * eig_vec1 = " << eigen_values[0].real() * eigen_vectors1 << std::endl;

	std::cout << "Mat * eig_vec2      = " << matcopy * eigen_vectors2 << std::endl;
    std::cout << "eig_val4 * eig_vec2 = " << eigen_values[1].real() * eigen_vectors2 << std::endl;

	std::cout << "Mat * eig_vec3      = " << matcopy * eigen_vectors3 << std::endl;
    std::cout << "eig_val5 * eig_vec3 = " << eigen_values[2].real() * eigen_vectors3 << std::endl;
}

void Demo_EigenSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     EIGENVALUE  SOLVERS                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Test_Symmetric_Eigen_Solver();    
    Test_Unsymmetric_Eigen_Solver();    
}