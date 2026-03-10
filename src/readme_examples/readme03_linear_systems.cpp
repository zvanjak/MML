///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme03_linear_systems.cpp                                         ///
///  Description: README Linear Systems & Eigenvalues section demo                    ///
///               Demonstrates LUSolver and EigenSolver                               ///
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
#include "algorithms/EigenSystemSolvers.h"
#endif

using namespace MML;

void Readme_LinearSystems()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            README - Linear Systems & Eigenvalues              ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // Define a linear system Ax = b
    Matrix<Real> A{5, 5, {0.2,  4.1, -2.1, -7.4,  1.6,
                          1.6,  1.5, -1.1,  0.7,  5.0,
                         -3.8, -8.0,  9.6, -5.4, -7.8,
                          4.6, -8.2,  8.4,  0.4,  8.0,
                         -2.6,  2.9,  0.1, -9.6, -2.7}};
    Vector<Real> b{1.1, 4.7, 0.1, 9.3, 0.4};

    // Solve using LU decomposition
    LUSolver<Real> luSolver(A);
    Vector<Real> x = luSolver.Solve(b);

    std::cout << "Solution: " << x << std::endl;
    std::cout << "Verification A*x: " << (A * x) << std::endl;
    std::cout << "Residual norm: " << (A*x - b).NormL2() << std::endl;

    // Compute eigenvalues
    auto eigenResult = EigenSolver::Solve(A);

    std::cout << "Eigenvalues:" << std::endl;
    for (const auto& ev : eigenResult.eigenvalues)
        std::cout << "  lambda = " << ev.real << " + " << ev.imag << "i" << std::endl;

    // QR decomposition
    QRSolver<Real> qr(A);
    Matrix<Real> Q = qr.GetQ();
    Matrix<Real> R = qr.GetR();
    std::cout << "QR valid: " << Matrix<Real>::AreEqual(A, Q * R, 1e-10) << std::endl;

/* Expected OUTPUT:
    Solution: [   -5.568500786,    -5.944693206,    -5.007620645,    -1.393638021,     3.598760994]
    Verification A*x: [            1.1,             4.7,             0.1,             9.3,             0.4]
    Residual norm: 1.24445e-14
    Eigenvalues:
      lambda = -2.470087376 + 12.99433106i
      ...
    QR valid: 1
*/
}
