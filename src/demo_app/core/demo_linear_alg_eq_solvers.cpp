#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/LinAlgEqSolvers.h"

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


void Demo_LinearAlgEqSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     LINEAR ALG. EQ. SOLVERS                   ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Test_GaussJordan_solver();    
}