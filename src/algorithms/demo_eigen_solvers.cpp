#include <iostream>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"

#include "core/CoreUtils.h"

#include "algorithms/EigenSystemSolvers.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

void Test_Symmetric_Eigen_Solver()
{
    std::cout << "CALCULATING EIGENVALUES OF SYMMETRIC MATRIX:\n";

    MatrixSym<Real>     origMat = TestBeds::symm_mat_5x5;
    Matrix<Real>        matcopy = origMat.GetAsMatrix();

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);

    SymmMatEigenSolver eigen_solver(origMat, true);

    std::cout << "Eigenvalues " << eigen_solver.d << std::endl;
    std::cout << "Eigenvectors " << eigen_solver.z << std::endl;

    Vector<Real>       eigen_values = eigen_solver.d;
    Vector<Real>       eigen_vectors1 = eigen_solver.z.VectorFromColumn(0);
    Vector<Real>       eigen_vectors2 = eigen_solver.z.VectorFromColumn(1);
    Vector<Real>       eigen_vectors3 = eigen_solver.z.VectorFromColumn(2);
    Vector<Real>       eigen_vectors4 = eigen_solver.z.VectorFromColumn(3);
    Vector<Real>       eigen_vectors5 = eigen_solver.z.VectorFromColumn(4);

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

void Test_Unsymmetric_Eigen_Solver_Single_Mat(Matrix<Real> origMat)
{
    std::cout << "*******************************************************************\n";
    std::cout << "*****   CALCULATING EIGENVALUES OF UNSYMMETRIC MATRIX   " << origMat.RowNum() << "x" << origMat.ColNum() << "   *****\n";

    int width = 15;
    int prec  = 10;

    int dim = origMat.RowNum();
    Matrix<Real>  matcopy(origMat);

    EigenSolver eigen_solver(matcopy, true, false);

    std::cout << "Initial matrix:\n";  matcopy.Print(std::cout,10,3);  std::cout << std::endl;

    std::cout << "Num real eigenvalues    : " << eigen_solver.getNumReal() << std::endl;
    std::cout << "Num complex eigenvalues : " << eigen_solver.getNumComplex() << "\n\n";
    
    std::cout << "Eigenvalues         : "; eigen_solver.getEigenvalues().Print(std::cout,width,prec); std::cout << std::endl;
    std::cout << "Real eigenvalues    : "; eigen_solver.getRealEigenvalues().Print(std::cout,width,prec); std::cout << std::endl;
    std::cout << "Complex eigenvalues : "; eigen_solver.getComplexEigenvalues().Print(std::cout,width,prec); std::cout << "\n\n";
     
    std::cout << "Eigenvectors matrix:\n"; eigen_solver.zz.Print(std::cout, width, 7); std::cout << std::endl;

    Matrix<Complex> mat_cmplx = MatrixUtils::CmplxMatFromRealMat(matcopy);

    for(int i=0; i<dim; i++ )
    {
        if( eigen_solver.isRealEigenvalue(i) )
        {
            double        eigen_val = eigen_solver.getEigenvalues()[i].real();
            Vector<Real>  eigen_vec = eigen_solver.getRealPartEigenvector(i);

            std::cout << "Eigen value  " << i+1 << " (real)    = " << eigen_val << std::endl;
            std::cout << "Eigen vector " << i+1 << " (real)    = "; eigen_vec.Print(std::cout,width,prec); std::cout << std::endl;

            std::cout << "Mat * eig_vec            = "; (matcopy * eigen_vec).Print(std::cout,width,prec); std::cout << std::endl;
            std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";            
        }
        else
        {
            Complex          eigen_val = eigen_solver.getEigenvalues()[i];
            Vector<Complex>  eigen_vec = eigen_solver.getEigenvector(i);
            
            std::cout << "Eigen value  " << i+1 << " (complex) = " << eigen_val << std::endl;
            std::cout << "Eigen vector " << i+1 << " (complex) = "; eigen_solver.getEigenvector(i).Print(std::cout, width, prec); std::cout << std::endl;

            std::cout << "Mat * eig_vec            = "; (mat_cmplx * eigen_vec).Print(std::cout, width, prec); std::cout << std::endl;
            std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";            
        }
    }
}

void Demo_EigenSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     EIGENVALUE  SOLVERS                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    //Test_Symmetric_Eigen_Solver();    

    Test_Unsymmetric_Eigen_Solver_Single_Mat(TestBeds::mat_50x50_1); 
    // Test_Unsymmetric_Eigen_Solver_Single_Mat(TestBeds::mat_6x6_test1); 
    // Test_Unsymmetric_Eigen_Solver_Single_Mat(TestBeds::mat_8x8_test1); 
}