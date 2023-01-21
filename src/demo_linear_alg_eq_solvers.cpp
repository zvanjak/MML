#include <iostream>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/LinAlgEqSolvers.h"
#include "../test_data/linear_alg_eq_test_bed.h"
#endif

void Test_GaussJordan_solver()
{
    std::cout << "SOLVING VIA GAUSS-JORDAN ELIMINATION:\n";

    MML::Matrix<Real>     origMat = MML::Tests::LinearAlgEqTestBed::mat4;
    MML::Matrix<Real>     matcopy(origMat);
    MML::Matrix<Real>     rhscopy(MML::Tests::LinearAlgEqTestBed::mat4_rhs4);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    MML::GaussJordanSolver::Solve(matcopy, rhscopy);

    MML::Vector<Real> solVec = MML::Matrix<Real>::VectorFromColumn(rhscopy, 0);
    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying matrix with solution vector: " << origMat * solVec << std::endl;

    std::cout << "Inverse " << matcopy << std::endl;
    std::cout << "Orig * inv = " << origMat * matcopy << std::endl;
}

void Test_LU_decomposition_solver()
{
    std::cout << "\nSOLVING VIA LU DECOMPOSITION:\n";
    
    MML::Matrix<Real>     origMat = MML::Tests::LinearAlgEqTestBed::mat4;
    MML::Matrix<Real>     matcopy(origMat);
    MML::Matrix<Real>     rhscopy(MML::Tests::LinearAlgEqTestBed::mat4_rhs4);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n";        rhscopy.Print(std::cout,10,3);

    MML::Matrix<Real>     matSol(5, 2);

    MML::LUDecompositionSolver luSolver(matcopy);
    luSolver.Solve(rhscopy, matSol);

    MML::Vector<Real> solVec = MML::Matrix<Real>::VectorFromColumn(rhscopy, 0);
    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;

    MML::Matrix<Real> inv(5,5);
    luSolver.inverse(inv);

    std::cout << "Inverse " << inv << std::endl;
    std::cout << "Orig * inve = " << origMat * inv << std::endl;    
}

void Test_QR_decomposition_solver()
{
    std::cout << "\nSOLVING VIA QR DECOMPOSITION:\n";

    MML::Matrix<Real>   mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                    1.6, 1.5, 1.1, 0.7, 5.0,
                                    3.8, 8.0, 9.6, 5.4, 8.8,
                                    4.6, 8.2, 8.4, 0.4, 8.0,
                                    2.6, 2.9, 0.1, 9.6, 7.7});

    MML::Vector<Real>     solVec(5), vecrhs4{1.1, 4.7, 0.1, 9.3, 0.4 };

    MML::Matrix<Real>     origMat = mat4;
    MML::Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10,3);
    std::cout << "Right side:\n" << vecrhs4 << std::endl;

    MML::QRDecompositionSolver qrSolver(matcopy);
    qrSolver.Solve(vecrhs4, solVec);

    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;  
}

void Test_SVD_decomposition_solver()
{
    std::cout << "\nSOLVING VIA SV DECOMPOSITION:\n";

    MML::Matrix<Real>     mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7});

    MML::Vector<Real>     vecrhs4{1.1, 4.7, 0.1, 9.3, 0.4 }, solVec(5);

    MML::Matrix<Real>     origMat = mat4;
    MML::Matrix<Real>     matcopy(origMat);

    std::cout << "Initial matrix:\n";    matcopy.Print(std::cout,10, 3);
    std::cout << "Right side:\n";        std::cout << vecrhs4 << std::endl;

    MML::SVDecompositionSolver svdSolver(matcopy);
    svdSolver.solve(vecrhs4, solVec);

    std::cout << "Solution:\n" << solVec << std::endl;
    std::cout << "Multiplying solution with matrix: " << origMat * solVec << std::endl;
}

void Demo_LinearAlgEqSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     LINEAR ALG. EQ. SOLVERS                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Test_GaussJordan_solver();    
    Test_LU_decomposition_solver();
    Test_QR_decomposition_solver();
    Test_SVD_decomposition_solver();
}