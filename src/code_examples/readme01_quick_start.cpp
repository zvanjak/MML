///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme01_quick_start.cpp                                            ///
///  Description: README Quick Start / First Program section demo                     ///
///               Demonstrates solving, determinant, and eigenvalue calculations      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/core/LinAlgEqSolvers.h>
#include <mml/algorithms/EigenSystemSolvers.h>
#endif

using namespace MML;

void Readme_QuickStart()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - Quick Start                         ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Create a matrix and solve a linear system
    Matrix<Real> A{3, 3, {4,  1,  2,
                          1, -5, -3,
                         -1,  1,  6}};
    Vector<Real> b{1, 4, -3};
    
    LUSolver<Real> solver(A);
    Vector<Real> x = solver.Solve(b);
    
    std::cout << "Solution: " << x << std::endl;
    std::cout << "Residual: " << (A * x - b).NormL2() << std::endl;
    std::cout << "Determinant: " << solver.det() << std::endl;

    auto eigenResult = EigenSolver::Solve(A);
    std::cout << "Eigenvalues:" << std::endl;
    for (const auto& eigenvalue : eigenResult.eigenvalues)
        std::cout << "  " << eigenvalue << std::endl;

/* Expected OUTPUT:
Solution: [   0.5378151261,   -0.4957983193,   -0.3277310924]
Residual: 0.0000000000
Determinant: -119.0000000000
Eigenvalues:
   4.8937181442 + 0.9530217116i
   4.8937181442 - 0.9530217116i
  -4.7874362884 + 0.0000000000i
*/
}


