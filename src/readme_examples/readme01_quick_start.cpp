///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme01_quick_start.cpp                                            ///
///  Description: README Quick Start / First Program section demo                     ///
///               Demonstrates basic matrix creation and linear system solving        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "core/LinAlgEqSolvers.h"
#endif

using namespace MML;

void Readme_QuickStart()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - Quick Start                         ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Create a matrix and solve a linear system
    Matrix<Real> A{3, 3, {4, 1, 2, 
                          1, 5, 1, 
                          2, 1, 6}};
    Vector<Real> b{1, 2, 3};
    
    LUSolver<Real> solver(A);
    Vector<Real> x = solver.Solve(b);
    
    std::cout << "Solution: " << x << std::endl;
    std::cout << "Residual: " << (A * x - b).NormL2() << std::endl;

/* Expected OUTPUT:
    Solution: [ -0.06382978723,    0.3191489362,    0.4680851064]
    Residual: 0
*/
}
