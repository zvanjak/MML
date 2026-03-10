///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_matrix_utils.cpp                                          ///
///  Description: Demonstration of MatrixUtils.h functionality                        ///
///               - Matrix property checks (triangular, symmetric, etc.)              ///
///               - Matrix norms and scalar quantities                                ///
///               - SVD-based utilities (rank, null space, fundamental subspaces)     ///
///               - Pseudoinverse and condition number                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Matrix/Matrix.h"
#include "core/MatrixUtils.h"
#endif

#include <iostream>
#include <iomanip>
#include <string>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                              MATRIX PROPERTY CHECKS                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_PropertyChecks()
{
    std::cout << "\n=== MATRIX PROPERTY CHECKS ===\n" << std::endl;
    
    // Create various test matrices
    Matrix<Real> upper(3, 3, {1.0, 2.0, 3.0,
                              0.0, 4.0, 5.0,
                              0.0, 0.0, 6.0});
    
    Matrix<Real> lower(3, 3, {1.0, 0.0, 0.0,
                              2.0, 3.0, 0.0,
                              4.0, 5.0, 6.0});
    
    Matrix<Real> diag(3, 3, {2.0, 0.0, 0.0,
                             0.0, 3.0, 0.0,
                             0.0, 0.0, 4.0});
    
    Matrix<Real> symmetric(3, 3, {1.0, 2.0, 3.0,
                                   2.0, 4.0, 5.0,
                                   3.0, 5.0, 6.0});
    
    Matrix<Real> skewSymmetric(3, 3, {0.0, 2.0, -3.0,
                                      -2.0, 0.0, 4.0,
                                       3.0, -4.0, 0.0});
    
    Matrix<Real> hessenberg(4, 4, {1.0, 2.0, 3.0, 4.0,
                                    5.0, 6.0, 7.0, 8.0,
                                    0.0, 9.0, 10.0, 11.0,
                                    0.0, 0.0, 12.0, 13.0});
    
    std::cout << "Upper triangular matrix:" << std::endl;
    upper.Print(std::cout, 6, 2);
    std::cout << "  IsUpperTriangular: " << (Utils::IsUpperTriangular(upper) ? "true" : "false") << std::endl;
    std::cout << "  IsLowerTriangular: " << (Utils::IsLowerTriangular(upper) ? "true" : "false") << std::endl;
    std::cout << "  IsDiagonal: " << (Utils::IsDiagonal(upper) ? "true" : "false") << std::endl;
    
    std::cout << "\nLower triangular matrix:" << std::endl;
    lower.Print(std::cout, 6, 2);
    std::cout << "  IsUpperTriangular: " << (Utils::IsUpperTriangular(lower) ? "true" : "false") << std::endl;
    std::cout << "  IsLowerTriangular: " << (Utils::IsLowerTriangular(lower) ? "true" : "false") << std::endl;
    
    std::cout << "\nDiagonal matrix:" << std::endl;
    diag.Print(std::cout, 6, 2);
    std::cout << "  IsDiagonal: " << (Utils::IsDiagonal(diag) ? "true" : "false") << std::endl;
    std::cout << "  IsUpperTriangular: " << (Utils::IsUpperTriangular(diag) ? "true" : "false") << std::endl;
    std::cout << "  IsLowerTriangular: " << (Utils::IsLowerTriangular(diag) ? "true" : "false") << std::endl;
    
    std::cout << "\nSymmetric matrix:" << std::endl;
    symmetric.Print(std::cout, 6, 2);
    std::cout << "  IsSymmetric: " << (Utils::IsSymmetric(symmetric) ? "true" : "false") << std::endl;
    std::cout << "  IsSkewSymmetric: " << (Utils::IsSkewSymmetric(symmetric) ? "true" : "false") << std::endl;
    
    std::cout << "\nSkew-symmetric matrix:" << std::endl;
    skewSymmetric.Print(std::cout, 6, 2);
    std::cout << "  IsSymmetric: " << (Utils::IsSymmetric(skewSymmetric) ? "true" : "false") << std::endl;
    std::cout << "  IsSkewSymmetric: " << (Utils::IsSkewSymmetric(skewSymmetric) ? "true" : "false") << std::endl;
    
    std::cout << "\nUpper Hessenberg matrix:" << std::endl;
    hessenberg.Print(std::cout, 6, 2);
    std::cout << "  IsUpperHessenberg: " << (Utils::IsUpperHessenberg(hessenberg) ? "true" : "false") << std::endl;
    std::cout << "  IsUpperTriangular: " << (Utils::IsUpperTriangular(hessenberg) ? "true" : "false") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MATRIX NORMS AND SCALAR QUANTITIES                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_NormsAndScalars()
{
    std::cout << "\n=== MATRIX NORMS AND SCALAR QUANTITIES ===\n" << std::endl;
    
    Matrix<Real> A(3, 3, {1.0, -2.0, 3.0,
                          4.0,  5.0, -6.0,
                         -7.0,  8.0,  9.0});
    
    std::cout << "Matrix A:" << std::endl;
    A.Print(std::cout, 6, 2);
    
    // Trace
    std::cout << "\nTrace(A) = " << Utils::Trace(A) << " (sum of diagonal: 1 + 5 + 9 = 15)" << std::endl;
    
    // Matrix norms
    std::cout << "\nMatrix Norms:" << std::endl;
    std::cout << "  Frobenius norm ||A||_F = " << Utils::FrobeniusNorm(A) << std::endl;
    std::cout << "  1-norm ||A||_1 = " << Utils::OneNorm(A) << " (max column sum)" << std::endl;
    std::cout << "  Infinity norm ||A||_inf = " << Utils::InfinityNorm(A) << " (max row sum)" << std::endl;
    
    // Determinant
    std::cout << "\nDeterminant: det(A) = " << Utils::Det(A) << std::endl;
    
    // Compare two matrices
    Matrix<Real> B = A;
    B(1, 1) += 0.5;
    std::cout << "\nMatrix B (A with B[1,1] += 0.5):" << std::endl;
    B.Print(std::cout, 6, 2);
    std::cout << "Max absolute difference |A - B| = " << Utils::MaxAbsDiff(A, B) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        RANK AND SPECIAL MATRIX PROPERTIES                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_RankAndSpecialProperties()
{
    std::cout << "\n=== RANK AND SPECIAL MATRIX PROPERTIES ===\n" << std::endl;
    
    // Full rank matrix
    Matrix<Real> fullRank(3, 3, {1.0, 2.0, 3.0,
                                  0.0, 1.0, 4.0,
                                  5.0, 6.0, 0.0});
    
    // Rank-deficient matrix (row 3 = row 1 + row 2)
    Matrix<Real> rankDeficient(3, 3, {1.0, 2.0, 3.0,
                                       4.0, 5.0, 6.0,
                                       5.0, 7.0, 9.0});
    
    std::cout << "Full rank matrix:" << std::endl;
    fullRank.Print(std::cout, 6, 2);
    std::cout << "  Rank (Gaussian): " << Utils::RankGaussian(fullRank) << std::endl;
    std::cout << "  Rank (SVD): " << Utils::RankSVD(fullRank) << std::endl;
    std::cout << "  Determinant: " << Utils::Det(fullRank) << std::endl;
    
    std::cout << "\nRank-deficient matrix (row 3 = row 1 + row 2):" << std::endl;
    rankDeficient.Print(std::cout, 6, 2);
    std::cout << "  Rank (Gaussian): " << Utils::RankGaussian(rankDeficient) << std::endl;
    std::cout << "  Rank (SVD): " << Utils::RankSVD(rankDeficient) << std::endl;
    std::cout << "  Nullity: " << Utils::Nullity(rankDeficient) << std::endl;
    std::cout << "  Determinant: " << Utils::Det(rankDeficient) << " (≈0 since rank < n)" << std::endl;
    
    // Nilpotent matrix: N^k = 0 for some k
    Matrix<Real> nilpotent(3, 3, {0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   0.0, 0.0, 0.0});
    
    std::cout << "\nNilpotent matrix (N^3 = 0):" << std::endl;
    nilpotent.Print(std::cout, 6, 2);
    std::cout << "  IsNilpotent: " << (Utils::IsNilpotent(nilpotent) ? "true" : "false") << std::endl;
    
    // Unipotent matrix: (U-I)^k = 0 for some k
    Matrix<Real> unipotent(3, 3, {1.0, 1.0, 0.0,
                                   0.0, 1.0, 1.0,
                                   0.0, 0.0, 1.0});
    
    std::cout << "\nUnipotent matrix ((U-I)^k = 0):" << std::endl;
    unipotent.Print(std::cout, 6, 2);
    std::cout << "  IsUnipotent: " << (Utils::IsUnipotent(unipotent) ? "true" : "false") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           SVD-BASED LINEAR ALGEBRA                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_SVDUtilities()
{
    std::cout << "\n=== SVD-BASED LINEAR ALGEBRA ===\n" << std::endl;
    
    // Well-conditioned matrix
    Matrix<Real> wellCond(3, 3, {10.0, 2.0, 3.0,
                                  4.0, 15.0, 6.0,
                                  7.0, 8.0, 20.0});
    
    std::cout << "Well-conditioned matrix:" << std::endl;
    wellCond.Print(std::cout, 7, 2);
    
    // Full SVD result
    Utils::SVDResult svdResult = Utils::ComputeSVD(wellCond);
    
    std::cout << "\nSingular values: ";
    svdResult.singularValues.Print(std::cout, 8, 4);
    std::cout << "Rank: " << svdResult.rank << std::endl;
    std::cout << "Condition number: " << svdResult.conditionNumber << std::endl;
    
    // Just get condition number directly
    std::cout << "\nCondition number (direct): " << Utils::ConditionNumber(wellCond) << std::endl;
    
    // Ill-conditioned matrix
    Matrix<Real> illCond(3, 3, {1.0, 1.0, 1.0,
                                 1.0, 1.0, 1.00001,
                                 1.0, 1.00001, 1.00002});
    
    std::cout << "\nIll-conditioned matrix (nearly singular):" << std::endl;
    illCond.Print(std::cout, 10, 5);
    std::cout << "Condition number: " << std::scientific << Utils::ConditionNumber(illCond) << std::fixed << std::endl;
    std::cout << "(Large condition number → digits of precision lost when solving Ax = b)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                          FUNDAMENTAL SUBSPACES                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_FundamentalSubspaces()
{
    std::cout << "\n=== FUNDAMENTAL SUBSPACES ===\n" << std::endl;
    
    // 3x4 rank-2 matrix to demonstrate non-trivial null space
    Matrix<Real> A(3, 4, {1.0, 2.0, 3.0, 4.0,
                           2.0, 4.0, 6.0, 8.0,  // Row 2 = 2 * Row 1
                           1.0, 3.0, 2.0, 5.0});
    
    std::cout << "Matrix A (3x4, rank 2):" << std::endl;
    A.Print(std::cout, 6, 2);
    
    Utils::FundamentalSubspaces fs = Utils::ComputeFundamentalSubspaces(A);
    
    std::cout << "\nRank: " << fs.rank << std::endl;
    std::cout << "Expected: rank(A) + nullity(A) = " << fs.rank << " + " 
              << (A.cols() - fs.rank) << " = " << A.cols() << " = number of columns" << std::endl;
    
    std::cout << "\nColumn space (range) - " << fs.columnSpace.cols() << " basis vectors in R^3:" << std::endl;
    if (fs.columnSpace.cols() > 0)
        fs.columnSpace.Print(std::cout, 8, 4);
    else
        std::cout << "  (empty - matrix is zero)" << std::endl;
    
    std::cout << "\nRow space - " << fs.rowSpace.cols() << " basis vectors in R^4:" << std::endl;
    if (fs.rowSpace.cols() > 0)
        fs.rowSpace.Print(std::cout, 8, 4);
    
    std::cout << "\nNull space (kernel) - " << fs.nullSpace.cols() << " basis vectors in R^4:" << std::endl;
    if (fs.nullSpace.cols() > 0) {
        fs.nullSpace.Print(std::cout, 8, 4);
        
        // Verify: A * null_vector should be zero
        std::cout << "\nVerification: A * (first null space vector) =" << std::endl;
        Vector<Real> nullVec(fs.nullSpace.rows());
        for (int i = 0; i < fs.nullSpace.rows(); i++)
            nullVec[i] = fs.nullSpace(i, 0);
        Vector<Real> result = A * nullVec;
        result.Print(std::cout, 10, 6);
        std::cout << "(Should be ≈ 0)" << std::endl;
    }
    
    std::cout << "\nLeft null space - " << fs.leftNullSpace.cols() << " basis vectors in R^3:" << std::endl;
    if (fs.leftNullSpace.cols() > 0)
        fs.leftNullSpace.Print(std::cout, 8, 4);
    else
        std::cout << "  (trivial - rank equals number of rows)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                            PSEUDOINVERSE                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Pseudoinverse()
{
    std::cout << "\n=== PSEUDOINVERSE (MOORE-PENROSE) ===\n" << std::endl;
    
    // Rectangular matrix (overdetermined case)
    Matrix<Real> A(4, 2, {1.0, 2.0,
                           3.0, 4.0,
                           5.0, 6.0,
                           7.0, 8.0});
    
    std::cout << "Overdetermined system A (4x2):" << std::endl;
    A.Print(std::cout, 6, 2);
    
    Matrix<Real> Aplus = Utils::PseudoInverse(A);
    std::cout << "\nPseudoinverse A+ (2x4):" << std::endl;
    Aplus.Print(std::cout, 10, 6);
    
    // Verify pseudoinverse properties
    Matrix<Real> AAplus = A * Aplus;
    Matrix<Real> AplusA = Aplus * A;
    
    std::cout << "\nVerification of pseudoinverse properties:" << std::endl;
    std::cout << "A * A+ (should be 4x4 projection matrix):" << std::endl;
    AAplus.Print(std::cout, 10, 6);
    
    std::cout << "\nA+ * A (should be 2x2 identity for full column rank):" << std::endl;
    AplusA.Print(std::cout, 10, 6);
    
    // Least squares solution
    Vector<Real> b({1.0, 2.0, 3.0, 4.0});
    Vector<Real> x = Aplus * b;
    
    std::cout << "\nLeast squares solution x = A+ * b:" << std::endl;
    std::cout << "b = "; b.Print(std::cout, 6, 2);
    std::cout << "x = A+ * b = "; x.Print(std::cout, 10, 6);
    
    Vector<Real> Ax = A * x;
    std::cout << "A * x = "; Ax.Print(std::cout, 10, 6);
    std::cout << "(Best approximation to b in the column space of A)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                      FADDEEV-LEVERRIER ALGORITHM                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_FaddeevLeverrier()
{
    std::cout << "\n=== FADDEEV-LEVERRIER ALGORITHM ===\n" << std::endl;
    
    // 3x3 matrix with known eigenvalues
    Matrix<Real> A(3, 3, {1.0, 2.0, 3.0,
                           0.0, 4.0, 5.0,
                           0.0, 0.0, 6.0});  // Triangular: eigenvalues are 1, 4, 6
    
    std::cout << "Upper triangular matrix (eigenvalues: 1, 4, 6):" << std::endl;
    A.Print(std::cout, 6, 2);
    
    PolynomRealFunc charPoly;
    Real det;
    Matrix<Real> inv(3, 3);
    
    Utils::FaddeevAlg(A, charPoly, det, inv);
    
    std::cout << "\nCharacteristic polynomial coefficients:" << std::endl;
    std::cout << "det(λI - A) = ";
    for (int i = charPoly.degree(); i >= 0; i--) {
        if (i < charPoly.degree() && charPoly[i] >= 0) std::cout << "+";
        std::cout << charPoly[i];
        if (i > 1) std::cout << "λ^" << i << " ";
        else if (i == 1) std::cout << "λ ";
    }
    std::cout << std::endl;
    
    std::cout << "\nDeterminant: " << det << " (expected: 1*4*6 = 24)" << std::endl;
    
    std::cout << "\nInverse matrix:" << std::endl;
    inv.Print(std::cout, 10, 6);
    
    // Verify: A * A^(-1) = I
    Matrix<Real> product = A * inv;
    std::cout << "\nVerification A * A^(-1):" << std::endl;
    product.Print(std::cout, 10, 6);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         SIMILARITY TRANSFORMATION                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_SimilarityTransform()
{
    std::cout << "\n=== SIMILARITY TRANSFORMATION ===\n" << std::endl;
    
    // Diagonal matrix with eigenvalues 1, 2, 3
    Matrix<Real> D(3, 3, {1.0, 0.0, 0.0,
                           0.0, 2.0, 0.0,
                           0.0, 0.0, 3.0});
    
    // Orthogonal transformation matrix (rotation)
    Real angle = Constants::PI / 4;  // 45 degrees
    Matrix<Real> Q(3, 3, {std::cos(angle), -std::sin(angle), 0.0,
                           std::sin(angle),  std::cos(angle), 0.0,
                           0.0,              0.0,             1.0});
    
    std::cout << "Diagonal matrix D:" << std::endl;
    D.Print(std::cout, 6, 2);
    
    std::cout << "\nOrthogonal matrix Q (rotation by 45°):" << std::endl;
    Q.Print(std::cout, 8, 4);
    
    Matrix<Real> A = Utils::SimilarityTransform(Q, D);
    std::cout << "\nSimilar matrix A = Q^T * D * Q:" << std::endl;
    A.Print(std::cout, 10, 6);
    
    std::cout << "\nTrace(D) = " << Utils::Trace(D) << std::endl;
    std::cout << "Trace(A) = " << Utils::Trace(A) << " (preserved under similarity)" << std::endl;
    
    std::cout << "\nDet(D) = " << Utils::Det(D) << std::endl;
    std::cout << "Det(A) = " << Utils::Det(A) << " (preserved under similarity)" << std::endl;
    
    std::cout << "\n(Trace and determinant are invariant under similarity transformation)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                               MAIN ENTRY POINT                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MatrixUtils()
{
    std::cout << "###################################################################" << std::endl;
    std::cout << "###               MatrixUtils.h - DEMONSTRATION                 ###" << std::endl;
    std::cout << "###################################################################" << std::endl;
    
    Demo_PropertyChecks();
    Demo_NormsAndScalars();
    Demo_RankAndSpecialProperties();
    Demo_SVDUtilities();
    Demo_FundamentalSubspaces();
    Demo_Pseudoinverse();
    Demo_FaddeevLeverrier();
    Demo_SimilarityTransform();
    
    std::cout << "\n###################################################################" << std::endl;
    std::cout << "###                   DEMONSTRATION COMPLETE                    ###" << std::endl;
    std::cout << "###################################################################\n" << std::endl;
}
