#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/LinAlgEqSolvers.h"

#include "algorithms/MatrixAlg.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

void Matrix_Simple_demo()
{
    // initialize a matrix with 5 rows and 3 columns
    Matrix<Real> mat(3, 3, {1, 2, 3,
                                 4, 5, 6,
                                 7, 8, 9});
                                 
    // right side vector
    Vector<Real> rhs({1, 2, 3});

    GaussJordanSolver<Real>::Solve(mat, rhs);
}

void Test_GaussJordan_solver()
{
    std::cout << "SOLVING VIA GAUSS-JORDAN ELIMINATION:\n";

    Matrix<Real>     origMat = TestBeds::mat_5x5;
    Matrix<Real>     matcopy(origMat);
    Matrix<Real>     rhscopy(TestBeds::mat_5x5_rhs_multi);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    GaussJordanSolver<Real>::Solve(matcopy, rhscopy);

    Vector<Real> solVec = rhscopy.VectorFromColumn(0);
    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying matrix with solution vector: " << origMat * solVec << std::endl;

    std::cout << "Inverse " << matcopy << std::endl;
    std::cout << "Orig * inv = " << origMat * matcopy << std::endl;
}

void Test_LU_decomposition_solver()
{
    std::cout << "\nSOLVING VIA LU DECOMPOSITION:\n";
    
    Matrix<Real>     origMat = TestBeds::mat_5x5;
    Matrix<Real>     matcopy(origMat);
    Matrix<Real>     rhscopy(TestBeds::mat_5x5_rhs_multi);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    Matrix<Real>     matSol(5, 2);

    LUDecompositionSolver<Real> luSolver(matcopy);
    luSolver.Solve(rhscopy, matSol);

    Vector<Real> solVec = rhscopy.VectorFromColumn(0);
    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;

    Matrix<Real> inv(5,5);
    luSolver.inverse(inv);

    std::cout << "Inverse " << inv << std::endl;
    std::cout << "Orig * inve = " << origMat * inv << std::endl;    
}

void Test_QR_decomposition_solver()
{
    std::cout << "\nSOLVING VIA QR DECOMPOSITION:\n";

    Matrix<Real>   mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                    1.6, 1.5, 1.1, 0.7, 5.0,
                                    3.8, 8.0, 9.6, 5.4, 8.8,
                                    4.6, 8.2, 8.4, 0.4, 8.0,
                                    2.6, 2.9, 0.1, 9.6, 7.7});

    Vector<Real>     solVec(5), vecrhs4{1.1, 4.7, 0.1, 9.3, 0.4 };

    Matrix<Real>     origMat = mat4;
    Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n" << vecrhs4 << std::endl;

    QRDecompositionSolver qrSolver(matcopy);
    qrSolver.Solve(vecrhs4, solVec);

    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;  
}

void Test_Cholesky_decomposition_solver()
{
    std::cout << "\nSOLVING VIA Cholesky DECOMPOSITION:\n";

    // example positive definite matrix
    Matrix<Real>  mat4{3, 3, { 2.0, -1.0, 0.0, 
                              -1.0,  3.0, -2.0, 
                               0.0, -2.0,  4.0 }};

    bool b = MatrixUtils::IsPositiveDefinite(mat4);

    Vector<Real>     solVec(3), vecrhs4{1.1, 4.7, 0.1 };

    Matrix<Real>     origMat = mat4;
    Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n" << vecrhs4 << std::endl;

    CholeskyDecompositionSolver choleskySolver(matcopy);
    choleskySolver.Solve(vecrhs4, solVec);

    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;  
}

void Test_SVD_decomposition_solver()
{
    std::cout << "\nSOLVING VIA SVD DECOMPOSITION:\n";

    Matrix<Real>     mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7});


    Vector<Real>     vecrhs4{1.1, 4.7, 0.1, 9.3, 0.4 }, solVec(5);

    Matrix<Real>     origMat = mat4;
    Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10, 3);
    std::cout << "Right side:\n";        std::cout << vecrhs4 << std::endl;

    SVDecompositionSolver svdSolver(matcopy);
    svdSolver.Solve(vecrhs4, solVec);

    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;
}

void Demo_LinearAlgEqSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     LINEAR ALG. EQ. SOLVERS                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Test_GaussJordan_solver();    
    // Test_LU_decomposition_solver();
    // Test_QR_decomposition_solver();
    Test_Cholesky_decomposition_solver();
//    Test_SVD_decomposition_solver();
}