#include "MMLBase.h"


#include "../include/base/Vector.h"
#include "../include/base/Matrix.h"
#include "../include/base/BaseUtils.h"

#include "../include/core/CoreUtils.h"
#include "../include/core/LinAlgEqSolvers.h"

#include "../include/algorithms/EigenSystemSolvers.h"

#include "../test_data/linear_alg_eq_systems_real_defs.h"
#include "../test_data/linear_alg_eq_systems_complex_defs.h"

using namespace MML;

// TODO - LOW, EASY - generate random complex matrix, i ispis da se lako copy paste u complex_defs.h
Matrix<Real> GenerateRandomMatrix(int n, int m, Real min, Real max)
{
    Matrix<Real> ret(n, m);

    srand((unsigned int)time(NULL));

    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
        {
            Real rnd = min + rand() / (Real)RAND_MAX * (max - min);
            ret(i,j) = rnd;
        }

    return ret;
}

void MatrixAnalyzer(std::string matName, const Matrix<Real> &origMat, const Vector<Real> &rhs)
{
    std::cout << "*******************************************************************\n";
    std::cout << "*****   MATRIX ANALYZER FOR MATRIX   " << origMat.RowNum() << "x" << origMat.ColNum() << "   *****\n";

    int width = 15;
    int prec  = 10;

    int dim = origMat.RowNum();
    Matrix<Real>  matcopy(origMat);
    Vector<Real>  rhscopy(rhs);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    Vector<Real>     vecSol(rhs.size());

    LUDecompositionSolver<Real> luSolver(matcopy);
    luSolver.Solve(rhscopy, vecSol);

    std::cout << "\nSolution:\n" << vecSol << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * vecSol << std::endl;

    Matrix<Real> inv(rhs.size(),rhs.size());
    luSolver.inverse(inv);

    std::cout << "Inverse " << inv << std::endl;
    std::cout << "Orig * inv = " << origMat * inv << std::endl;    

    EigenSolver eigen_solver(matcopy, true, false);

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
    
    std::cout << "    const static inline MML::Vector     " << matName << "_rhs0_sol{";
    for(int i=0; i<dim; i++ )
    {
        std::cout << std::setprecision(15) << vecSol[i];
        if( i < dim-1 )
            std::cout << ", ";
        else
            std::cout << "};\n\n";
    }

    std::cout << "    const static inline VectorComplex   " << matName << "_eigen_val{ Complex(" << eigen_solver.getEigenvalues()[0].real() << ", " << eigen_solver.getEigenvalues()[0].imag() << "),\n";
    for(int i=1; i<dim; i++ )
    {
        std::cout << "                                                             Complex(" << eigen_solver.getEigenvalues()[i].real() << ", " << eigen_solver.getEigenvalues()[i].imag() << ")";
        if( i < dim-1 )
            std::cout << ",\n";
        else
            std::cout << "};\n";
    }

    std::cout << "    const static inline std::vector<VectorComplex> " << matName << "_eigen_vecs\n";
    std::cout << std::setprecision(10) << "    {\n";
    for(int i=0; i<dim; i++ )
    {
        std::cout << "        VectorComplex{ Complex(" << eigen_solver.getEigenvector(i)[0].real() << ", " << eigen_solver.getEigenvector(i)[0].imag() << "),";
        for(int j=1; j<dim; j++ )
        {
            std::cout << "Complex(" << eigen_solver.getEigenvector(i)[j].real() << ", " << eigen_solver.getEigenvector(i)[j].imag() << ")";
            if( j < dim-1 )
                std::cout << ", ";
            else
                std::cout << " },\n";
        }
        // if( i < dim-1 )
        //     std::cout << ",\n";
        // else
        //     std::cout << "\n";
    }
    std::cout << "    };\n";
}


void SymmMatrixAnalyzer(std::string matName, const MatrixSym<Real> &origMat, const Vector<Real> &rhs)
{
    std::cout << "*******************************************************************\n";
    std::cout << "*****   SYMMETRIC MATRIX ANALYZER FOR MATRIX   " << origMat.RowNum() << "x" << origMat.ColNum() << "   *****\n";

    int width = 15;
    int prec  = 10;

    int dim = origMat.RowNum();
    Vector<Real>  rhscopy(rhs);

    std::cout << "Initial matrix:\n";    origMat.Print(std::cout,10,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    auto matcopy = origMat.GetAsMatrix();
    Vector<Real>   vecSol(rhs.size());
    
    LUDecompositionSolver<Real> luSolver(matcopy);
    luSolver.Solve(rhscopy, vecSol);

    std::cout << "\nSolution:\n" << vecSol << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * vecSol << std::endl;

    Matrix<Real> inv(rhs.size(),rhs.size());
    luSolver.inverse(inv);

    std::cout << "Inverse " << inv << std::endl;
    std::cout << "Orig * inv = " << origMat * inv << std::endl;    

    SymmMatEigenSolver eigen_solver(origMat, true);

    std::cout << "Eigenvalues         : "; eigen_solver.getEigenvalues().Print(std::cout,width,prec); std::cout << std::endl;
    std::cout << "Eigenvectors matrix:\n"; eigen_solver.z.Print(std::cout, width, 7); std::cout << std::endl;

    auto matcopy2 = origMat.GetAsMatrix();
    for(int i=0; i<dim; i++ )
    {
        double        eigen_val = eigen_solver.getEigenvalues()[i];
        Vector<Real>  eigen_vec = eigen_solver.getEigenvector(i);

        std::cout << "Eigen value  " << i+1 << " (real)    = " << eigen_val << std::endl;
        std::cout << "Eigen vector " << i+1 << " (real)    = "; eigen_vec.Print(std::cout,width,prec); std::cout << std::endl;

        std::cout << "Mat * eig_vec            = "; (matcopy2 * eigen_vec).Print(std::cout,width,prec); std::cout << std::endl;
        std::cout << "eig_val * eig_vec        = "; (eigen_val * eigen_vec).Print(std::cout, width, prec); std::cout << "\n\n";            
    }
    
    std::cout << "    const static inline MML::Vector     " << matName << "_rhs0_sol{";
    for(int i=0; i<dim; i++ )
    {
        std::cout << std::setprecision(15) << vecSol[i];
        if( i < dim-1 )
            std::cout << ", ";
        else
            std::cout << "};\n\n";
    }

    std::cout << "    const static inline Vector<Real> " << matName << "_eigen_val{ " << eigen_solver.getEigenvalues()[0] << ", ";
    for(int i=1; i<dim; i++ )
    {
        std::cout << eigen_solver.getEigenvalues()[i];
        if( i < dim-1 )
            std::cout << ", ";
        else
            std::cout << " };\n";
    }

    std::cout << "    const static inline std::vector<Vector<Real>> " << matName << "_eigen_vecs\n";
    std::cout << std::setprecision(10) << "    {\n";
    for(int i=0; i<dim; i++ )
    {
        std::cout << "        Vector<Real>{ " << eigen_solver.getEigenvector(i)[0] << ", ";
        for(int j=1; j<dim; j++ )
        {
            std::cout  << eigen_solver.getEigenvector(i)[j];
            if( j < dim-1 )
                std::cout << ", ";
            else
                std::cout << " },\n";
        }
    }
    std::cout << "    };\n";
}


void ComplexMatrixAnalyzer(std::string matName, const Matrix<Complex> &origMat, const Vector<Complex> &rhs)
{
    std::cout << "*******************************************************************\n";
    std::cout << "*****   COMPLEX MATRIX ANALYZER FOR MATRIX   " << origMat.RowNum() << "x" << origMat.ColNum() << "   *****\n";

    int width = 15;
    int prec  = 10;

    int dim = origMat.RowNum();
    Matrix<Complex>  matcopy(origMat);
    Vector<Complex>  rhscopy(rhs);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,15,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,15,3);

    Vector<Complex>     vecSol(rhs.size());

    LUDecompositionSolver<Complex> luSolver(matcopy);
    luSolver.Solve(rhscopy, vecSol);

    std::cout << "\nSolution:\n" << vecSol << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * vecSol << std::endl;

    Matrix<Complex> inv(rhs.size(),rhs.size());
    luSolver.inverse(inv);

    std::cout << "Inverse:\n";
    inv.Print(std::cout, 20, 5);
    std::cout << "\nOrig * inv:\n";
    (origMat * inv).Print(std::cout, 20, 5, 1e-10);
    
    std::cout << "\n    const static inline MML::Vector<Complex>  " << matName << "_rhs0_sol{ ";
    for(int i=0; i<dim; i++ )
    {
        std::cout << "Complex";
        std::cout << std::setprecision(15) << vecSol[i];
        if( i < dim-1 )
            std::cout << ", \n";
        else
            std::cout << "};\n\n";
    }

    std::cout << "    };\n";
}


int main()
{
    // auto mat = GenerateRandomMatrix(50, 50, -100.0, 100.0);

    // // mat.Print(std::cout, 7, 2);
    // std::cout << "{ ";
    // for(int i=0; i<mat.RowNum(); i++)
    // {
    //     for(int j=0; j<mat.ColNum(); j++)
    //         std::cout << std::setw(6) << std::setprecision(4) << mat(i,j) << ", ";
        
    //     std::cout  << std::endl;
    // }

//    MatrixAnalyzer("mat_50x50", TestBeds::mat_50x50, TestBeds::mat_50x50_rhs0);
//    SymmMatrixAnalyzer("symm_mat_10x10", TestBeds::symm_mat_10x10, TestBeds::symm_mat_10x10_rhs0);
    ComplexMatrixAnalyzer("mat_cmplx_1_5x5", TestBeds::mat_cmplx_1_5x5, TestBeds::mat_cmplx_1_5x5_rhs0);

}