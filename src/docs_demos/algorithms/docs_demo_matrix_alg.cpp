///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_matrix_alg.cpp                                            ///
///  Description: Demonstration of MatrixAlg.h functionality                          ///
///               - Hessenberg reduction                                              ///
///               - Matrix definiteness classification                                ///
///               - Positive/negative definite checks                                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#include "algorithms/MatrixAlg.h"
#include "core/MatrixUtils.h"
#endif

#include <iostream>
#include <iomanip>
#include <string>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                            HELPER: DEFINITENESS TO STRING                           ///
///////////////////////////////////////////////////////////////////////////////////////////

std::string DefinitenessToString(MatrixAlg::Definiteness def)
{
    switch (def) {
        case MatrixAlg::Definiteness::PositiveDefinite:      return "Positive Definite";
        case MatrixAlg::Definiteness::PositiveSemiDefinite:  return "Positive Semi-Definite";
        case MatrixAlg::Definiteness::NegativeDefinite:      return "Negative Definite";
        case MatrixAlg::Definiteness::NegativeSemiDefinite:  return "Negative Semi-Definite";
        case MatrixAlg::Definiteness::Indefinite:            return "Indefinite";
        default:                                              return "Unknown";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           HESSENBERG REDUCTION                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_HessenbergReduction()
{
    std::cout << "\n=== HESSENBERG REDUCTION ===\n" << std::endl;
    
    std::cout << "Hessenberg form H = Q^T * A * Q has zeros below the first subdiagonal." << std::endl;
    std::cout << "This is an important preprocessing step for eigenvalue algorithms.\n" << std::endl;
    
    // Create a general 4x4 matrix
    Matrix<Real> A(4, 4, {4.0,  1.0, -2.0,  2.0,
                          1.0,  2.0,  0.0,  1.0,
                         -2.0,  0.0,  3.0, -2.0,
                          2.0,  1.0, -2.0, -1.0});
    
    std::cout << "Original matrix A:" << std::endl;
    A.Print(std::cout, 8, 3);
    
    // Reduce to Hessenberg form
    MatrixAlg::HessenbergResult result = MatrixAlg::ReduceToHessenberg(A);
    
    std::cout << "\nHessenberg form H:" << std::endl;
    result.H.Print(std::cout, 10, 5);
    
    std::cout << "\nOrthogonal transformation matrix Q:" << std::endl;
    result.Q.Print(std::cout, 10, 5);
    
    // Verify the transformation
    std::cout << "\n--- Verification ---" << std::endl;
    std::cout << "Checking H = Q^T * A * Q:" << std::endl;
    
    Matrix<Real> QtAQ = result.Q.transpose() * A * result.Q;
    Real maxError = Utils::MaxAbsDiff(result.H, QtAQ);
    std::cout << "Max difference ||H - Q^T*A*Q||_∞ = " << std::scientific << maxError << std::fixed << std::endl;
    
    std::cout << "\nChecking Q is orthogonal (Q^T * Q = I):" << std::endl;
    Matrix<Real> QtQ = result.Q.transpose() * result.Q;
    Matrix<Real> I = Matrix<Real>::Identity(4);
    maxError = Utils::MaxAbsDiff(QtQ, I);
    std::cout << "Max difference ||Q^T*Q - I||_∞ = " << std::scientific << maxError << std::fixed << std::endl;
    
    std::cout << "\nIs H upper Hessenberg? " << (Utils::IsUpperHessenberg(result.H) ? "Yes" : "No") << std::endl;
    
    std::cout << "\nTrace(A) = " << Utils::Trace(A) << ", Trace(H) = " << Utils::Trace(result.H) << std::endl;
    std::cout << "(Trace is preserved under similarity transformation)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DEFINITENESS CLASSIFICATION                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_DefinitenessClassification()
{
    std::cout << "\n=== DEFINITENESS CLASSIFICATION ===\n" << std::endl;
    
    std::cout << "A symmetric matrix A is:" << std::endl;
    std::cout << "  - Positive definite if x^T*A*x > 0 for all x ≠ 0" << std::endl;
    std::cout << "  - Positive semi-definite if x^T*A*x >= 0 for all x" << std::endl;
    std::cout << "  - Negative definite if x^T*A*x < 0 for all x ≠ 0" << std::endl;
    std::cout << "  - Indefinite if it has both positive and negative eigenvalues\n" << std::endl;
    
    // Positive definite: identity matrix
    MatrixSym<Real> identity(3);
    for (int i = 0; i < 3; i++) identity(i, i) = 1.0;
    
    std::cout << "--- Identity Matrix (3x3) ---" << std::endl;
    std::cout << "Eigenvalues: 1, 1, 1 (all positive)" << std::endl;
    std::cout << "Classification: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(identity)) << std::endl;
    std::cout << "IsPositiveDefinite: " << (MatrixAlg::IsPositiveDefinite(identity) ? "Yes" : "No") << std::endl;
    
    // Positive definite: diagonally dominant
    MatrixSym<Real> posdef(3);
    posdef(0, 0) = 4.0; posdef(0, 1) = 1.0; posdef(0, 2) = 0.0;
                        posdef(1, 1) = 5.0; posdef(1, 2) = 1.0;
                                            posdef(2, 2) = 6.0;
    
    std::cout << "\n--- Positive Definite Matrix ---" << std::endl;
    std::cout << "Matrix:" << std::endl;
    // Print symmetric matrix
    Matrix<Real> pdfull = posdef.GetAsMatrix();
    pdfull.Print(std::cout, 6, 2);
    std::cout << "Classification: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(posdef)) << std::endl;
    std::cout << "IsPositiveDefinite: " << (MatrixAlg::IsPositiveDefinite(posdef) ? "Yes" : "No") << std::endl;
    
    // Positive semi-definite (rank deficient)
    MatrixSym<Real> possemidef(3);
    possemidef(0, 0) = 1.0; possemidef(0, 1) = 1.0; possemidef(0, 2) = 1.0;
                            possemidef(1, 1) = 1.0; possemidef(1, 2) = 1.0;
                                                    possemidef(2, 2) = 1.0;
    
    std::cout << "\n--- Positive Semi-Definite Matrix ---" << std::endl;
    std::cout << "Matrix (all 1's, rank 1):" << std::endl;
    Matrix<Real> psdfull = possemidef.GetAsMatrix();
    psdfull.Print(std::cout, 6, 2);
    std::cout << "Classification: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(possemidef)) << std::endl;
    std::cout << "IsPositiveSemiDefinite: " << (MatrixAlg::IsPositiveSemiDefinite(possemidef) ? "Yes" : "No") << std::endl;
    std::cout << "IsPositiveDefinite: " << (MatrixAlg::IsPositiveDefinite(possemidef) ? "Yes" : "No") << std::endl;
    
    // Negative definite
    MatrixSym<Real> negdef(3);
    negdef(0, 0) = -4.0; negdef(0, 1) = 1.0; negdef(0, 2) = 0.0;
                         negdef(1, 1) = -5.0; negdef(1, 2) = 1.0;
                                              negdef(2, 2) = -6.0;
    
    std::cout << "\n--- Negative Definite Matrix ---" << std::endl;
    std::cout << "Matrix:" << std::endl;
    Matrix<Real> ndfull = negdef.GetAsMatrix();
    ndfull.Print(std::cout, 6, 2);
    std::cout << "Classification: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(negdef)) << std::endl;
    std::cout << "IsNegativeDefinite: " << (MatrixAlg::IsNegativeDefinite(negdef) ? "Yes" : "No") << std::endl;
    
    // Indefinite
    MatrixSym<Real> indef(3);
    indef(0, 0) = 1.0;  indef(0, 1) = 0.0; indef(0, 2) = 0.0;
                        indef(1, 1) = -2.0; indef(1, 2) = 0.0;
                                            indef(2, 2) = 3.0;
    
    std::cout << "\n--- Indefinite Matrix ---" << std::endl;
    std::cout << "Matrix (diagonal with eigenvalues 1, -2, 3):" << std::endl;
    Matrix<Real> indefull = indef.GetAsMatrix();
    indefull.Print(std::cout, 6, 2);
    std::cout << "Classification: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(indef)) << std::endl;
    std::cout << "IsIndefinite: " << (MatrixAlg::IsIndefinite(indef) ? "Yes" : "No") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         DEFINITENESS FOR GENERAL MATRICES                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_GeneralMatrixDefiniteness()
{
    std::cout << "\n=== DEFINITENESS FOR GENERAL MATRICES ===\n" << std::endl;
    
    std::cout << "For non-symmetric matrices, definiteness uses the symmetric part:" << std::endl;
    std::cout << "  A_sym = (A + A^T) / 2\n" << std::endl;
    
    // Asymmetric but positive definite symmetric part
    Matrix<Real> A(3, 3, {5.0,  2.0, 1.0,
                          0.0,  4.0, 2.0,
                          0.0,  0.0, 3.0});
    
    std::cout << "Upper triangular matrix A:" << std::endl;
    A.Print(std::cout, 6, 2);
    
    std::cout << "\nSymmetric part (A + A^T) / 2:" << std::endl;
    Matrix<Real> Asym(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Asym(i, j) = 0.5 * (A(i, j) + A(j, i));
    Asym.Print(std::cout, 6, 2);
    
    std::cout << "\nClassification of A: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(A)) << std::endl;
    std::cout << "IsPositiveDefinite(A): " << (MatrixAlg::IsPositiveDefinite(A) ? "Yes" : "No") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                              SYMMETRY CHECKS                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_SymmetryChecks()
{
    std::cout << "\n=== SYMMETRY CHECKS ===\n" << std::endl;
    
    // Symmetric matrix
    Matrix<Real> sym(3, 3, {1.0, 2.0, 3.0,
                            2.0, 4.0, 5.0,
                            3.0, 5.0, 6.0});
    
    std::cout << "Symmetric matrix:" << std::endl;
    sym.Print(std::cout, 6, 2);
    std::cout << "IsSymmetric: " << (MatrixAlg::IsSymmetric(sym) ? "Yes" : "No") << std::endl;
    std::cout << "IsSkewSymmetric: " << (MatrixAlg::IsSkewSymmetric(sym) ? "Yes" : "No") << std::endl;
    
    // Skew-symmetric matrix
    Matrix<Real> skew(3, 3, { 0.0,  2.0, -3.0,
                             -2.0,  0.0,  4.0,
                              3.0, -4.0,  0.0});
    
    std::cout << "\nSkew-symmetric matrix:" << std::endl;
    skew.Print(std::cout, 6, 2);
    std::cout << "IsSymmetric: " << (MatrixAlg::IsSymmetric(skew) ? "Yes" : "No") << std::endl;
    std::cout << "IsSkewSymmetric: " << (MatrixAlg::IsSkewSymmetric(skew) ? "Yes" : "No") << std::endl;
    
    // General matrix
    Matrix<Real> gen(3, 3, {1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 9.0});
    
    std::cout << "\nGeneral matrix (neither symmetric nor skew-symmetric):" << std::endl;
    gen.Print(std::cout, 6, 2);
    std::cout << "IsSymmetric: " << (MatrixAlg::IsSymmetric(gen) ? "Yes" : "No") << std::endl;
    std::cout << "IsSkewSymmetric: " << (MatrixAlg::IsSkewSymmetric(gen) ? "Yes" : "No") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                       APPLICATIONS OF DEFINITENESS                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_DefinitenessApplications()
{
    std::cout << "\n=== APPLICATIONS OF DEFINITENESS ===\n" << std::endl;
    
    std::cout << "Matrix definiteness is important for:\n" << std::endl;
    
    std::cout << "1. OPTIMIZATION" << std::endl;
    std::cout << "   - Hessian positive definite at x* ⟹ x* is a local minimum" << std::endl;
    std::cout << "   - Hessian negative definite at x* ⟹ x* is a local maximum" << std::endl;
    std::cout << "   - Hessian indefinite at x* ⟹ x* is a saddle point\n" << std::endl;
    
    // Example: Hessian of f(x,y) = x² + 2y²
    MatrixSym<Real> hessianMin(2);
    hessianMin(0, 0) = 2.0;   // ∂²f/∂x²
    hessianMin(0, 1) = 0.0;   // ∂²f/∂x∂y
    hessianMin(1, 1) = 4.0;   // ∂²f/∂y²
    
    std::cout << "   Example: f(x,y) = x² + 2y²" << std::endl;
    std::cout << "   Hessian at (0,0):" << std::endl;
    Matrix<Real> hm = hessianMin.GetAsMatrix();
    std::cout << "   "; hm.Print(std::cout, 6, 2);
    std::cout << "   Definiteness: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(hessianMin)) << std::endl;
    std::cout << "   ⟹ (0,0) is a minimum\n" << std::endl;
    
    // Example: Saddle point
    MatrixSym<Real> hessianSaddle(2);
    hessianSaddle(0, 0) = 2.0;   // ∂²f/∂x²
    hessianSaddle(0, 1) = 0.0;   // ∂²f/∂x∂y  
    hessianSaddle(1, 1) = -4.0;  // ∂²f/∂y²
    
    std::cout << "   Example: f(x,y) = x² - 2y²" << std::endl;
    std::cout << "   Hessian at (0,0):" << std::endl;
    Matrix<Real> hs = hessianSaddle.GetAsMatrix();
    std::cout << "   "; hs.Print(std::cout, 6, 2);
    std::cout << "   Definiteness: " << DefinitenessToString(MatrixAlg::ClassifyDefiniteness(hessianSaddle)) << std::endl;
    std::cout << "   ⟹ (0,0) is a saddle point\n" << std::endl;
    
    std::cout << "2. LINEAR SYSTEMS" << std::endl;
    std::cout << "   - Positive definite A allows Cholesky factorization: A = LL^T" << std::endl;
    std::cout << "   - Conjugate gradient method converges for positive definite A\n" << std::endl;
    
    std::cout << "3. COVARIANCE MATRICES" << std::endl;
    std::cout << "   - Covariance matrices are always positive semi-definite" << std::endl;
    std::cout << "   - Positive definite ⟹ no redundant variables\n" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                               MAIN ENTRY POINT                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MatrixAlg()
{
    std::cout << "###################################################################" << std::endl;
    std::cout << "###                MatrixAlg.h - DEMONSTRATION                  ###" << std::endl;
    std::cout << "###################################################################" << std::endl;
    
    Demo_HessenbergReduction();
    Demo_DefinitenessClassification();
    Demo_GeneralMatrixDefiniteness();
    Demo_SymmetryChecks();
    Demo_DefinitenessApplications();
    
    std::cout << "\n###################################################################" << std::endl;
    std::cout << "###                   DEMONSTRATION COMPLETE                    ###" << std::endl;
    std::cout << "###################################################################\n" << std::endl;
}
