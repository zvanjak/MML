#include <iostream>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/LinAlgEqSolvers.h"
#endif

void Test_gaussj()
{
    std::cout << "SOLVING VIA GAUSS-JORDAN ELIMINATION:\n";

    MML::Matrix     m1(2, 2, {1.0, -1.0, 
                              1.5, 3.0});
    MML::Matrix     m1_rhs(2, 1, {1.0, 2.0});

    MML::Matrix     mat0(3, 3, {1.0, 2.0, -1.0, 
                                -1.0, 5.0, 6.0, 
                                3.0, 1.0, 1.0 });
    MML::Matrix     rhs0(3, 1, {1.0, 2.0, 1.0});

    MML::Matrix     mat1(3, 3, {1.0, 0.0, 0.0,
                                0.0, 2.0, 0.0,
                                0.0, 0.0, 3.0});
    MML::Matrix     rhs1(3, 2, {1.0, 0.0, 0.0,
                                1.0, 1.0, 1.0});

    MML::Matrix     mat2(3, 3, {1.0, 2.0, 3.0,
                                2.0, 2.0, 3.0,
                                3.0, 3.0, 3.0});
    MML::Matrix     rhs2(3, 2, {1.0, 1.0, 1.0,
                                1.0, 2.0, 3.0});

    MML::Matrix     mat3(5, 5, {1.0, 2.0, 3.0, 4.0, 5.0,
                                2.0, 3.0, 4.0, 5.0, 1.0,
                                3.0, 4.0, 5.0, 1.0, 2.0,
                                4.0, 5.0, 1.0, 2.0, 3.0,
                                5.0, 1.0, 2.0, 3.0, 4.0});
    MML::Matrix     rhs3(5, 2, {1.0, 1.0, 1.0, 1.0, 1.0,
                                1.0, 2.0, 3.0, 4.0, 5.0});

    MML::Matrix     mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7});
    MML::Matrix     rhs4(5, 2, {1.1, 1.6, 
                                4.7, 9.1, 
                                0.1, 4.0, 
                                9.3, 8.4, 
                                0.4, 4.1});    
 
    
    // MML::Matrix&    origMat = mat4;
    // MML::Matrix     matcopy(origMat);
    // MML::Matrix     rhscopy(rhs4);

    MML::Matrix&    origMat = m1;
    MML::Matrix     matcopy(origMat);
    MML::Matrix     rhscopy(m1_rhs);

    std::cout << "Initial matrix:\n";
    matcopy.Print(std::cout,10,3);

    std::cout << "Right side:\n";
    rhscopy.Print(std::cout,10,3);

    MML::GaussJordanSolver::Solve(matcopy, rhscopy);

    std::cout << "Solution:\n";
    MML::Vector solVec = MML::Matrix::VectorFromColumn(rhscopy, 0);
    //MML::Vector    solVec({rhscopy[0][0], rhscopy[1][0], rhscopy[2][0], rhscopy[3][0], rhscopy[4][0]});
    std::cout << solVec << std::endl;

    // inicijaliziramo rhs1Vec sa rješenjem
    MML::Vector    res = origMat * solVec;

    std::cout << "Multiplying solution with matrix: " << res << std::endl;

    std::cout << "Inverse " << matcopy << std::endl;
    std::cout << "Orig * inve = " << origMat * matcopy << std::endl;
}

void Test_LU_decomposition_solver()
{
    std::cout << "\nSOLVING VIA LU DECOMPOSITION:\n";

    MML::Matrix     mat0(3, 3, {1.0, 2.0, -1.0, 
                                -1.0, 5.0, 6.0, 
                                3.0, 1.0, 1.0 });
    MML::Matrix     rhs0(3, 1, {1.0, 2.0, 1.0});

    MML::Matrix     mat1(3, 3, {1.0, 0.0, 0.0,
                                0.0, 2.0, 0.0,
                                0.0, 0.0, 3.0});
    MML::Matrix     rhs1(3, 2, {1.0, 0.0, 0.0,
                                1.0, 1.0, 1.0});

    MML::Matrix     mat2(3, 3, {1.0, 2.0, 3.0,
                                2.0, 2.0, 3.0,
                                3.0, 3.0, 3.0});
    MML::Matrix     rhs2(3, 2, {1.0, 1.0, 1.0,
                                1.0, 2.0, 3.0});

    MML::Matrix     mat3(5, 5, {1.0, 2.0, 3.0, 4.0, 5.0,
                                2.0, 3.0, 4.0, 5.0, 1.0,
                                3.0, 4.0, 5.0, 1.0, 2.0,
                                4.0, 5.0, 1.0, 2.0, 3.0,
                                5.0, 1.0, 2.0, 3.0, 4.0});
    MML::Matrix     rhs3(5, 2, {1.0, 1.0, 1.0, 1.0, 1.0,
                                1.0, 2.0, 3.0, 4.0, 5.0});

    MML::Matrix     mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7});
    MML::Matrix     rhs4(5, 2, {1.1, 1.6, 
                                4.7, 9.1, 
                                0.1, 4.0, 
                                9.3, 8.4, 
                                0.4, 4.1});    
    
    MML::Matrix&    origMat = mat4;
    MML::Matrix     matcopy(origMat);
    MML::Matrix     rhscopy(rhs4);

    std::cout << "Initial matrix:\n";
    matcopy.Print(std::cout,10,3);

    std::cout << "Right side:\n";
    rhscopy.Print(std::cout,10,3);

    MML::LUDecompositionSolver luSolver(matcopy);

    MML::Matrix     matSol(5, 2);
    luSolver.Solve(rhscopy, matSol);

    std::cout << "Solution:\n";
    MML::Vector    solVec({matSol[0][0], matSol[1][0], matSol[2][0], matSol[3][0], matSol[4][0]});
    std::cout << solVec << std::endl;

    // inicijaliziramo rhs1Vec sa rješenjem
    MML::Vector    res = origMat * solVec;

    std::cout << "Multiplying solution with matrix: " << res << std::endl;
}

void Test_QR_decomposition_solver()
{
    std::cout << "\nSOLVING VIA QR DECOMPOSITION:\n";

    MML::Matrix     mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7});

    MML::Vector     vecrhs4{1.1, 4.7, 0.1, 9.3, 0.4 }, vecSol(5);

    MML::Matrix&    origMat = mat4;
    MML::Matrix     matcopy(origMat);

    std::cout << "Initial matrix:\n";
    matcopy.Print(std::cout,10,3);

    std::cout << "Right side:\n";
    std::cout << vecrhs4 << std::endl;

    MML::QRDecompositionSolver qrSolver(matcopy);

    qrSolver.Solve(vecrhs4, vecSol);

    std::cout << "Solution:\n";
    std::cout << vecSol << std::endl;

    // inicijaliziramo rhs1Vec sa rješenjem
    MML::Vector    res = origMat * vecSol;

    std::cout << "Multiplying solution with matrix: " << res << std::endl;
}


void Test_SVD_decomposition_solver()
{
    std::cout << "\nSOLVING VIA SV DECOMPOSITION:\n";

    MML::Matrix     mat4(5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                1.6, 1.5, 1.1, 0.7, 5.0,
                                3.8, 8.0, 9.6, 5.4, 8.8,
                                4.6, 8.2, 8.4, 0.4, 8.0,
                                2.6, 2.9, 0.1, 9.6, 7.7});

    MML::Vector     vecrhs4{1.1, 4.7, 0.1, 9.3, 0.4 }, vecSol(5);

    MML::Matrix&    origMat = mat4;
    MML::Matrix     matcopy(origMat);

    std::cout << "Initial matrix:\n";
    matcopy.Print(std::cout,10, 3);

    std::cout << "Right side:\n";
    std::cout << vecrhs4 << std::endl;

    MML::SVDecompositionSolver svdSolver(matcopy);

    svdSolver.solve(vecrhs4, vecSol);

    std::cout << "Solution:\n";
    std::cout << vecSol << std::endl;

    // inicijaliziramo rhs1Vec sa rješenjem
    MML::Vector    res = origMat * vecSol;

    std::cout << "Multiplying solution with matrix: " << res << std::endl;
}

void Demo_LinearAlgEqSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     LINEAR ALG. EQ. SOLVERS                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Test_gaussj();    
    Test_LU_decomposition_solver();
    Test_QR_decomposition_solver();
    Test_SVD_decomposition_solver();
}