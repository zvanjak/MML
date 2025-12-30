///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        test_precision_linalg.cpp                                           ///
///  Description: Comprehensive linear algebra precision tests                        ///
///               Compares LU, Cholesky, QR, SVD solvers and decompositions           ///
///                                                                                   ///
///  Tests:       - Linear system solving accuracy (Ax = b)                           ///
///               - Matrix decomposition verification                                 ///
///               - Condition number analysis                                         ///
///               - Ill-conditioned system behavior                                   ///
///               - Eigenvalue computation accuracy                                   ///
///               - Overdetermined system (least squares)                             ///
///               - Matrix inversion accuracy                                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "PrecisionTestFramework.h"

#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/MatrixSym.h"

#include "core/LinAlgEqSolvers.h"
#include "algorithms/MatrixAlg.h"
#include "algorithms/EigenSystemSolvers.h"

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace MML;
using namespace MML::PrecisionTesting;

///////////////////////////////////////////////////////////////////////////////////////////
//                        HELPER FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Create a Hilbert matrix of size n x n
 * 
 * Hilbert matrix: H[i][j] = 1 / (i + j + 1)
 * Known to be very ill-conditioned - condition number grows exponentially with n
 */
Matrix<Real> CreateHilbertMatrix(int n)
{
    Matrix<Real> H(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            H[i][j] = 1.0 / (i + j + 1);
    return H;
}

/**
 * @brief Create a Vandermonde matrix for given nodes
 * V[i][j] = x[i]^j
 */
Matrix<Real> CreateVandermondeMatrix(const Vector<Real>& x)
{
    int n = static_cast<int>(x.size());
    Matrix<Real> V(n, n);
    for (int i = 0; i < n; i++) {
        Real xi_power = 1.0;
        for (int j = 0; j < n; j++) {
            V[i][j] = xi_power;
            xi_power *= x[i];
        }
    }
    return V;
}

/**
 * @brief Create a symmetric positive definite matrix
 * Uses A = B * B^T where B is random to guarantee SPD
 */
Matrix<Real> CreateSPDMatrix(int n, unsigned seed = 42)
{
    Matrix<Real> B(n, n);
    std::srand(seed);
    
    // Create random matrix with diagonal dominance for stability
    for (int i = 0; i < n; i++) {
        Real row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                B[i][j] = (std::rand() % 100 - 50) / 50.0;  // [-1, 1]
                row_sum += std::abs(B[i][j]);
            }
        }
        // Make diagonally dominant
        B[i][i] = row_sum + 1.0 + (std::rand() % 100) / 100.0;
    }
    
    // A = B * B^T is SPD
    Matrix<Real> A(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0.0;
            for (int k = 0; k < n; k++)
                A[i][j] += B[i][k] * B[j][k];
        }
    }
    return A;
}

/**
 * @brief Create a symmetric positive definite matrix (MatrixSym)
 */
MatrixSym<Real> CreateSPDMatrixSym(int n, unsigned seed = 42)
{
    Matrix<Real> A = CreateSPDMatrix(n, seed);
    MatrixSym<Real> S(n);
    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
            S(i, j) = A[i][j];
    return S;
}

/**
 * @brief Create a known-eigenvalue matrix for testing
 * Diagonal matrix with eigenvalues 1, 2, ..., n
 */
Matrix<Real> CreateDiagonalEigenMatrix(int n)
{
    Matrix<Real> D(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            D[i][j] = (i == j) ? static_cast<Real>(i + 1) : 0.0;
        }
    }
    return D;
}

/**
 * @brief Create a rotation matrix in 2D (embedded in nxn)
 * Has eigenvalues e^{i*theta} - complex for non-trivial rotations
 */
Matrix<Real> CreateRotationMatrix(int n, Real theta)
{
    Matrix<Real> R = Matrix<Real>::GetUnitMatrix(n);
    Real c = std::cos(theta);
    Real s = std::sin(theta);
    R[0][0] = c;  R[0][1] = -s;
    R[1][0] = s;  R[1][1] = c;
    return R;
}

/**
 * @brief Compute matrix-vector product A * x
 */
Vector<Real> MatVecMul(const Matrix<Real>& A, const Vector<Real>& x)
{
    int m = A.RowNum();
    int n = A.ColNum();
    Vector<Real> y(m);
    for (int i = 0; i < m; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++)
            y[i] += A[i][j] * x[j];
    }
    return y;
}

/**
 * @brief Compute vector difference norm: ||x - y||
 */
Real VectorDiffNorm(const Vector<Real>& x, const Vector<Real>& y)
{
    Real sum = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    return std::sqrt(sum);
}

/**
 * @brief Compute vector 2-norm: ||x||
 */
Real VectorNorm(const Vector<Real>& x)
{
    Real sum = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        sum += x[i] * x[i];
    return std::sqrt(sum);
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 1: BASIC LINEAR SYSTEM SOLVING
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test basic linear system solving with well-conditioned matrix
 * 
 * Tests solving Ax = b where A has known solution x_exact = [1, 2, 3, ...]^T
 * and b = A * x_exact
 */
void Test_BasicLinearSystem(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 1: Basic Linear System Solving ===\n";
    
    const int n = 5;
    
    // Create a well-conditioned matrix (diagonally dominant)
    Matrix<Real> A(n, n);
    A[0][0] = 10; A[0][1] = -1; A[0][2] =  2; A[0][3] =  0; A[0][4] =  0;
    A[1][0] = -1; A[1][1] = 11; A[1][2] = -1; A[1][3] =  3; A[1][4] =  0;
    A[2][0] =  2; A[2][1] = -1; A[2][2] = 10; A[2][3] = -1; A[2][4] =  0;
    A[3][0] =  0; A[3][1] =  3; A[3][2] = -1; A[3][3] =  8; A[3][4] =  0;
    A[4][0] =  0; A[4][1] =  0; A[4][2] =  0; A[4][3] = -1; A[4][4] =  5;
    
    // Known exact solution
    Vector<Real> x_exact(n);
    for (int i = 0; i < n; i++)
        x_exact[i] = static_cast<Real>(i + 1);  // [1, 2, 3, 4, 5]
    
    // Compute RHS: b = A * x_exact
    Vector<Real> b = MatVecMul(A, x_exact);
    
    // Test Gauss-Jordan
    try {
        Matrix<Real> A_copy = A;
        Vector<Real> x_gj = GaussJordanSolver<Real>::SolveConst(A, b);
        Real error = VectorDiffNorm(x_gj, x_exact);
        suite.addResult("GaussJordan", "Basic5x5", 0.0, error);
        std::cout << "  GaussJordan: error = " << std::scientific << error << "\n";
    } catch (const std::exception& e) {
        std::cout << "  GaussJordan: FAILED - " << e.what() << "\n";
    }
    
    // Test LU Solver
    try {
        LUSolver<Real> lu(A);
        Vector<Real> x_lu = lu.Solve(b);
        Real error = VectorDiffNorm(x_lu, x_exact);
        suite.addResult("LU", "Basic5x5", 0.0, error);
        std::cout << "  LU:          error = " << std::scientific << error << "\n";
    } catch (const std::exception& e) {
        std::cout << "  LU:          FAILED - " << e.what() << "\n";
    }
    
    // Test QR Solver
    try {
        QRSolver<Real> qr(A);
        Vector<Real> x_qr = qr.Solve(b);
        Real error = VectorDiffNorm(x_qr, x_exact);
        suite.addResult("QR", "Basic5x5", 0.0, error);
        std::cout << "  QR:          error = " << std::scientific << error << "\n";
    } catch (const std::exception& e) {
        std::cout << "  QR:          FAILED - " << e.what() << "\n";
    }
    
    // Test SVD Solver
    try {
        SVDecompositionSolver svd(A);
        Vector<Real> x_svd = svd.Solve(b);
        Real error = VectorDiffNorm(x_svd, x_exact);
        suite.addResult("SVD", "Basic5x5", 0.0, error);
        std::cout << "  SVD:         error = " << std::scientific << error << "\n";
    } catch (const std::exception& e) {
        std::cout << "  SVD:         FAILED - " << e.what() << "\n";
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 2: SYMMETRIC POSITIVE DEFINITE SYSTEMS
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test Cholesky decomposition for SPD matrices
 */
void Test_CholeskySPD(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 2: Symmetric Positive Definite Systems (Cholesky) ===\n";
    
    // Test various SPD matrix sizes
    for (int n : {3, 5, 8}) {
        Matrix<Real> A = CreateSPDMatrix(n);
        
        // Known exact solution
        Vector<Real> x_exact(n);
        for (int i = 0; i < n; i++)
            x_exact[i] = static_cast<Real>(i + 1);
        
        Vector<Real> b = MatVecMul(A, x_exact);
        
        std::string func_name = "SPD" + std::to_string(n) + "x" + std::to_string(n);
        
        // Test Cholesky
        try {
            CholeskySolver<Real> chol(A);
            Vector<Real> x_chol = chol.Solve(b);
            Real error = VectorDiffNorm(x_chol, x_exact);
            suite.addResult("Cholesky", func_name, 0.0, error);
            std::cout << "  Cholesky " << n << "x" << n << ":  error = " << std::scientific << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  Cholesky " << n << "x" << n << ":  FAILED - " << e.what() << "\n";
        }
        
        // Compare with LU
        try {
            LUSolver<Real> lu(A);
            Vector<Real> x_lu = lu.Solve(b);
            Real error = VectorDiffNorm(x_lu, x_exact);
            suite.addResult("LU", func_name, 0.0, error);
            std::cout << "  LU       " << n << "x" << n << ":  error = " << std::scientific << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  LU       " << n << "x" << n << ":  FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 3: ILL-CONDITIONED SYSTEMS (HILBERT MATRIX)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test solver accuracy on Hilbert matrices (increasing ill-conditioning)
 * 
 * Hilbert matrices are notorious for being extremely ill-conditioned.
 * Condition number of H_n ≈ O((1+√2)^(4n) / √n)
 */
void Test_HilbertMatrix(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 3: Ill-Conditioned Systems (Hilbert Matrix) ===\n";
    std::cout << "  (Testing numerical stability as condition number increases)\n";
    
    for (int n : {3, 4, 5, 6}) {
        Matrix<Real> H = CreateHilbertMatrix(n);
        
        // Known exact solution
        Vector<Real> x_exact(n);
        for (int i = 0; i < n; i++)
            x_exact[i] = 1.0;  // All ones
        
        Vector<Real> b = MatVecMul(H, x_exact);
        
        std::string func_name = "Hilbert" + std::to_string(n);
        
        std::cout << "\n  Hilbert " << n << "x" << n << " matrix:\n";
        
        // Test LU
        try {
            LUSolver<Real> lu(H);
            Vector<Real> x_lu = lu.Solve(b);
            Real error = VectorDiffNorm(x_lu, x_exact);
            suite.addResult("LU", func_name, 0.0, error);
            std::cout << "    LU:  error = " << std::scientific << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "    LU:  FAILED - " << e.what() << "\n";
        }
        
        // Test QR
        try {
            QRSolver<Real> qr(H);
            Vector<Real> x_qr = qr.Solve(b);
            Real error = VectorDiffNorm(x_qr, x_exact);
            suite.addResult("QR", func_name, 0.0, error);
            std::cout << "    QR:  error = " << std::scientific << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "    QR:  FAILED - " << e.what() << "\n";
        }
        
        // Test SVD
        try {
            SVDecompositionSolver svd(H);
            Vector<Real> x_svd = svd.Solve(b);
            Real error = VectorDiffNorm(x_svd, x_exact);
            suite.addResult("SVD", func_name, 0.0, error);
            std::cout << "    SVD: error = " << std::scientific << error << "\n";
            
            // Also report condition number approximation
            Real inv_cond = svd.inv_condition();
            std::cout << "    (1/condition ≈ " << inv_cond << ")\n";
        } catch (const std::exception& e) {
            std::cout << "    SVD: FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 4: EIGENVALUE COMPUTATION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test eigenvalue computation accuracy for symmetric matrices
 */
void Test_EigenvalueAccuracy(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 4: Eigenvalue Computation (Jacobi) ===\n";
    
    // Test 1: Diagonal matrix - eigenvalues are diagonal elements
    {
        const int n = 5;
        MatrixSym<Real> D(n);
        Vector<Real> exact_eigenvalues(n);
        for (int i = 0; i < n; i++) {
            D(i, i) = static_cast<Real>(i + 1);  // Eigenvalues: 1, 2, 3, 4, 5
            exact_eigenvalues[i] = static_cast<Real>(i + 1);
        }
        
        auto result = SymmMatEigenSolverJacobi::Solve(D);
        
        // Sort exact eigenvalues (already sorted)
        Real max_error = 0.0;
        for (int i = 0; i < n; i++) {
            Real error = std::abs(result.eigenvalues[i] - exact_eigenvalues[i]);
            max_error = std::max(max_error, error);
        }
        
        suite.addResult("Jacobi", "Diag5x5", 0.0, max_error);
        std::cout << "  Diagonal 5x5:  max eigenvalue error = " << std::scientific << max_error << "\n";
        std::cout << "    Iterations: " << result.iterations << ", Converged: " << (result.converged ? "Yes" : "No") << "\n";
    }
    
    // Test 2: 2x2 symmetric matrix with known eigenvalues
    {
        MatrixSym<Real> A(2);
        A(0, 0) = 3.0;
        A(0, 1) = 1.0;
        A(1, 1) = 3.0;
        // Eigenvalues: 3±1 = {2, 4}
        Vector<Real> exact_eigenvalues(2);
        exact_eigenvalues[0] = 2.0;
        exact_eigenvalues[1] = 4.0;
        
        auto result = SymmMatEigenSolverJacobi::Solve(A);
        
        Real max_error = 0.0;
        for (int i = 0; i < 2; i++) {
            Real error = std::abs(result.eigenvalues[i] - exact_eigenvalues[i]);
            max_error = std::max(max_error, error);
        }
        
        suite.addResult("Jacobi", "Symm2x2", 0.0, max_error);
        std::cout << "  Symmetric 2x2: max eigenvalue error = " << std::scientific << max_error << "\n";
    }
    
    // Test 3: Random SPD matrix - verify A*v = λ*v
    {
        const int n = 4;
        MatrixSym<Real> S = CreateSPDMatrixSym(n, 123);
        
        auto result = SymmMatEigenSolverJacobi::Solve(S);
        
        // Verify A*v = λ*v for each eigenpair
        Real max_residual = 0.0;
        for (int i = 0; i < n; i++) {
            Vector<Real> v(n);
            for (int j = 0; j < n; j++)
                v[j] = result.eigenvectors[j][i];  // Column i
            
            // Compute A*v
            Vector<Real> Av(n);
            for (int j = 0; j < n; j++) {
                Av[j] = 0.0;
                for (int k = 0; k < n; k++)
                    Av[j] += (j <= k ? S(j, k) : S(k, j)) * v[k];
            }
            
            // Compute λ*v
            Vector<Real> lambda_v(n);
            for (int j = 0; j < n; j++)
                lambda_v[j] = result.eigenvalues[i] * v[j];
            
            Real residual = VectorDiffNorm(Av, lambda_v);
            max_residual = std::max(max_residual, residual);
        }
        
        suite.addResult("Jacobi", "SPD4x4_residual", 0.0, max_residual);
        std::cout << "  SPD 4x4:       max eigenpair residual = " << std::scientific << max_residual << "\n";
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 5: QR DECOMPOSITION VERIFICATION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Verify QR decomposition: A = Q*R, Q^T*Q = I, R upper triangular
 */
void Test_QRDecomposition(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 5: QR Decomposition Verification ===\n";
    
    for (int n : {3, 5, 8}) {
        // Create test matrix
        Matrix<Real> A(n, n);
        std::srand(42 + n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = (std::rand() % 200 - 100) / 10.0;
        
        std::string func_name = "QR" + std::to_string(n) + "x" + std::to_string(n);
        
        try {
            QRSolver<Real> qr(A);
            Matrix<Real> Q = qr.GetQ();
            Matrix<Real> R = qr.GetR();
            
            // Verify Q is orthogonal: Q^T * Q = I
            Real ortho_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real dot = 0.0;
                    for (int k = 0; k < n; k++)
                        dot += Q[k][i] * Q[k][j];
                    Real expected = (i == j) ? 1.0 : 0.0;
                    ortho_error = std::max(ortho_error, std::abs(dot - expected));
                }
            }
            
            // Verify R is upper triangular
            Real lower_error = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < i; j++)
                    lower_error = std::max(lower_error, std::abs(R[i][j]));
            
            // Verify A = Q*R
            Real recon_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real sum = 0.0;
                    for (int k = 0; k < n; k++)
                        sum += Q[i][k] * R[k][j];
                    recon_error = std::max(recon_error, std::abs(A[i][j] - sum));
                }
            }
            
            suite.addResult("QR_Ortho", func_name, 0.0, ortho_error);
            suite.addResult("QR_Upper", func_name, 0.0, lower_error);
            suite.addResult("QR_Recon", func_name, 0.0, recon_error);
            
            std::cout << "  " << n << "x" << n << ": Ortho=" << std::scientific << ortho_error
                      << ", Upper=" << lower_error
                      << ", Recon=" << recon_error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  " << n << "x" << n << ": FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 6: SVD DECOMPOSITION VERIFICATION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Verify SVD decomposition: A = U*S*V^T
 */
void Test_SVDDecomposition(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 6: SVD Decomposition Verification ===\n";
    
    for (int n : {3, 4, 5}) {
        // Create test matrix
        Matrix<Real> A(n, n);
        std::srand(123 + n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = (std::rand() % 200 - 100) / 10.0;
        
        std::string func_name = "SVD" + std::to_string(n) + "x" + std::to_string(n);
        
        try {
            SVDecompositionSolver svd(A);
            Matrix<Real> U = svd.getU();
            Matrix<Real> V = svd.getV();
            Vector<Real> W = svd.getW();
            
            // Verify U is orthogonal: U^T * U = I
            Real U_ortho_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real dot = 0.0;
                    for (int k = 0; k < n; k++)
                        dot += U[k][i] * U[k][j];
                    Real expected = (i == j) ? 1.0 : 0.0;
                    U_ortho_error = std::max(U_ortho_error, std::abs(dot - expected));
                }
            }
            
            // Verify V is orthogonal: V^T * V = I
            Real V_ortho_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real dot = 0.0;
                    for (int k = 0; k < n; k++)
                        dot += V[k][i] * V[k][j];
                    Real expected = (i == j) ? 1.0 : 0.0;
                    V_ortho_error = std::max(V_ortho_error, std::abs(dot - expected));
                }
            }
            
            // Verify A = U * diag(W) * V^T
            Real recon_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real sum = 0.0;
                    for (int k = 0; k < n; k++)
                        sum += U[i][k] * W[k] * V[j][k];  // V^T means V[j][k]
                    recon_error = std::max(recon_error, std::abs(A[i][j] - sum));
                }
            }
            
            // Check singular values are non-negative and descending
            bool sv_valid = true;
            for (int i = 0; i < n; i++) {
                if (W[i] < 0) sv_valid = false;
                if (i > 0 && W[i] > W[i-1]) sv_valid = false;
            }
            
            suite.addResult("SVD_U_Ortho", func_name, 0.0, U_ortho_error);
            suite.addResult("SVD_V_Ortho", func_name, 0.0, V_ortho_error);
            suite.addResult("SVD_Recon", func_name, 0.0, recon_error);
            
            std::cout << "  " << n << "x" << n << ": U_Ortho=" << std::scientific << U_ortho_error
                      << ", V_Ortho=" << V_ortho_error
                      << ", Recon=" << recon_error
                      << ", SV_Valid=" << (sv_valid ? "Yes" : "No") << "\n";
        } catch (const std::exception& e) {
            std::cout << "  " << n << "x" << n << ": FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 7: LU DECOMPOSITION VERIFICATION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Verify LU decomposition produces correct determinant
 */
void Test_LUDeterminant(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 7: LU Decomposition - Determinant Verification ===\n";
    
    // Test 1: 2x2 matrix with known determinant
    {
        Matrix<Real> A(2, 2);
        A[0][0] = 4; A[0][1] = 3;
        A[1][0] = 6; A[1][1] = 3;
        // det = 4*3 - 3*6 = 12 - 18 = -6
        Real exact_det = -6.0;
        
        try {
            LUSolver<Real> lu(A);
            Real computed_det = lu.det();
            Real error = std::abs(computed_det - exact_det);
            suite.addResult("LU_Det", "2x2", 0.0, error);
            std::cout << "  2x2: det=" << computed_det << " (exact=" << exact_det << "), error=" << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  2x2: FAILED - " << e.what() << "\n";
        }
    }
    
    // Test 2: 3x3 matrix with known determinant
    {
        Matrix<Real> A(3, 3);
        A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
        A[1][0] = 4; A[1][1] = 5; A[1][2] = 6;
        A[2][0] = 7; A[2][1] = 8; A[2][2] = 10;
        // det = 1*(50-48) - 2*(40-42) + 3*(32-35) = 2 + 4 - 9 = -3
        Real exact_det = -3.0;
        
        try {
            LUSolver<Real> lu(A);
            Real computed_det = lu.det();
            Real error = std::abs(computed_det - exact_det);
            suite.addResult("LU_Det", "3x3", 0.0, error);
            std::cout << "  3x3: det=" << computed_det << " (exact=" << exact_det << "), error=" << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  3x3: FAILED - " << e.what() << "\n";
        }
    }
    
    // Test 3: Diagonal matrix - det = product of diagonal
    {
        const int n = 5;
        Matrix<Real> D(n, n);
        Real exact_det = 1.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                D[i][j] = (i == j) ? static_cast<Real>(i + 1) : 0.0;
            exact_det *= (i + 1);  // 1 * 2 * 3 * 4 * 5 = 120
        }
        
        try {
            LUSolver<Real> lu(D);
            Real computed_det = lu.det();
            Real error = std::abs(computed_det - exact_det);
            suite.addResult("LU_Det", "Diag5x5", 0.0, error);
            std::cout << "  Diag5x5: det=" << computed_det << " (exact=" << exact_det << "), error=" << error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  Diag5x5: FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 8: MATRIX INVERSION ACCURACY
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test matrix inversion: A * A^{-1} = I
 */
void Test_MatrixInversion(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 8: Matrix Inversion Accuracy ===\n";
    
    for (int n : {3, 4, 5}) {
        // Create a non-singular matrix
        Matrix<Real> A(n, n);
        std::srand(456 + n);
        for (int i = 0; i < n; i++) {
            Real row_sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    A[i][j] = (std::rand() % 100 - 50) / 50.0;
                    row_sum += std::abs(A[i][j]);
                }
            }
            A[i][i] = row_sum + 1.0;  // Diagonally dominant
        }
        
        std::string func_name = "Inv" + std::to_string(n) + "x" + std::to_string(n);
        
        // Test LU inversion
        try {
            LUSolver<Real> lu(A);
            Matrix<Real> A_inv(n, n);
            lu.inverse(A_inv);
            
            // Compute A * A^{-1}
            Real identity_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real sum = 0.0;
                    for (int k = 0; k < n; k++)
                        sum += A[i][k] * A_inv[k][j];
                    Real expected = (i == j) ? 1.0 : 0.0;
                    identity_error = std::max(identity_error, std::abs(sum - expected));
                }
            }
            
            suite.addResult("LU_Inv", func_name, 0.0, identity_error);
            std::cout << "  LU " << n << "x" << n << ":      ||A*A^{-1} - I|| = " << std::scientific << identity_error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  LU " << n << "x" << n << ":      FAILED - " << e.what() << "\n";
        }
        
        // Test Cholesky inversion (for SPD)
        Matrix<Real> S = CreateSPDMatrix(n, 789 + n);
        try {
            CholeskySolver<Real> chol(S);
            Matrix<Real> S_inv(n, n);
            chol.inverse(S_inv);
            
            Real identity_error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Real sum = 0.0;
                    for (int k = 0; k < n; k++)
                        sum += S[i][k] * S_inv[k][j];
                    Real expected = (i == j) ? 1.0 : 0.0;
                    identity_error = std::max(identity_error, std::abs(sum - expected));
                }
            }
            
            suite.addResult("Chol_Inv", func_name, 0.0, identity_error);
            std::cout << "  Cholesky " << n << "x" << n << ": ||S*S^{-1} - I|| = " << std::scientific << identity_error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  Cholesky " << n << "x" << n << ": FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 9: LEAST SQUARES (OVERDETERMINED SYSTEMS)
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test least squares solution for overdetermined systems
 */
void Test_LeastSquares(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 9: Least Squares (Overdetermined Systems) ===\n";
    
    // Linear regression: fit y = mx + c to points
    // Points: (0, 0.1), (1, 2.1), (2, 3.9), (3, 6.1), (4, 8.0)
    // Expected: m ≈ 2, c ≈ 0
    {
        const int m = 5;  // equations
        const int n = 2;  // unknowns (m and c)
        
        Matrix<Real> A(m, n);
        Vector<Real> b(m);
        
        // Setup: A[i] = [x_i, 1], b[i] = y_i
        Real x_vals[] = {0, 1, 2, 3, 4};
        Real y_vals[] = {0.1, 2.1, 3.9, 6.1, 8.0};
        
        for (int i = 0; i < m; i++) {
            A[i][0] = x_vals[i];  // coefficient for m
            A[i][1] = 1.0;        // coefficient for c
            b[i] = y_vals[i];
        }
        
        // QR least squares
        try {
            QRSolver<Real> qr(A);
            Vector<Real> x = qr.LeastSquaresSolve(b);
            
            // Expected: m ≈ 2.0, c ≈ 0.0
            Real slope_error = std::abs(x[0] - 2.0);
            Real intercept_error = std::abs(x[1] - 0.0);
            
            suite.addResult("QR_LS", "LinReg", 0.0, slope_error);
            std::cout << "  QR: slope=" << x[0] << " (err=" << slope_error << ")"
                      << ", intercept=" << x[1] << " (err=" << intercept_error << ")\n";
        } catch (const std::exception& e) {
            std::cout << "  QR: FAILED - " << e.what() << "\n";
        }
        
        // SVD least squares
        try {
            SVDecompositionSolver svd(A);
            Vector<Real> x = svd.Solve(b);
            
            Real slope_error = std::abs(x[0] - 2.0);
            Real intercept_error = std::abs(x[1] - 0.0);
            
            suite.addResult("SVD_LS", "LinReg", 0.0, slope_error);
            std::cout << "  SVD: slope=" << x[0] << " (err=" << slope_error << ")"
                      << ", intercept=" << x[1] << " (err=" << intercept_error << ")\n";
        } catch (const std::exception& e) {
            std::cout << "  SVD: FAILED - " << e.what() << "\n";
        }
    }
    
    // Polynomial fit: fit y = ax^2 + bx + c to 10 noisy points
    {
        const int m = 10;  // equations
        const int n = 3;   // unknowns (a, b, c)
        
        Matrix<Real> A(m, n);
        Vector<Real> b(m);
        
        // Generate points from y = x^2 - 2x + 1 with small noise
        std::srand(999);
        for (int i = 0; i < m; i++) {
            Real x = static_cast<Real>(i) * 0.5;
            A[i][0] = x * x;  // x^2
            A[i][1] = x;      // x
            A[i][2] = 1.0;    // 1
            Real noise = (std::rand() % 100 - 50) / 500.0;  // [-0.1, 0.1]
            b[i] = x * x - 2 * x + 1 + noise;
        }
        
        // QR least squares
        try {
            QRSolver<Real> qr(A);
            Vector<Real> x = qr.LeastSquaresSolve(b);
            
            // Expected: a ≈ 1.0, b ≈ -2.0, c ≈ 1.0
            Real a_error = std::abs(x[0] - 1.0);
            Real b_error = std::abs(x[1] - (-2.0));
            Real c_error = std::abs(x[2] - 1.0);
            Real max_error = std::max({a_error, b_error, c_error});
            
            suite.addResult("QR_LS", "PolyFit", 0.0, max_error);
            std::cout << "  QR Poly: a=" << x[0] << ", b=" << x[1] << ", c=" << x[2] << ", max_err=" << max_error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  QR Poly: FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        TEST 10: HESSENBERG REDUCTION
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Test Hessenberg reduction: verify H = Q^T * A * Q
 */
void Test_HessenbergReduction(PrecisionTestSuite& suite)
{
    std::cout << "\n=== Test 10: Hessenberg Reduction ===\n";
    
    for (int n : {3, 4, 5}) {
        Matrix<Real> A(n, n);
        std::srand(555 + n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = (std::rand() % 200 - 100) / 10.0;
        
        std::string func_name = "Hess" + std::to_string(n) + "x" + std::to_string(n);
        
        try {
            auto result = MatrixAlg::ReduceToHessenberg(A);
            
            // Check H is upper Hessenberg
            bool is_hess = MatrixAlg::IsUpperHessenberg(result.H, 1e-10);
            
            // Check Q is orthogonal
            bool is_ortho = MatrixAlg::IsOrthogonal(result.Q, 1e-10);
            
            // Verify H = Q^T * A * Q
            Matrix<Real> QTAQ = MatrixAlg::SimilarityTransform(result.Q, A);
            Real recon_error = MatrixAlg::MaxAbsDiff(result.H, QTAQ);
            
            suite.addResult("Hessenberg", func_name, 0.0, recon_error);
            std::cout << "  " << n << "x" << n << ": Hess=" << (is_hess ? "Yes" : "No")
                      << ", Ortho=" << (is_ortho ? "Yes" : "No")
                      << ", Recon=" << std::scientific << recon_error << "\n";
        } catch (const std::exception& e) {
            std::cout << "  " << n << "x" << n << ": FAILED - " << e.what() << "\n";
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//                        MAIN TEST RUNNER
///////////////////////////////////////////////////////////////////////////////////////////

void Test_Precision_LinearAlgebra()
{
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                         LINEAR ALGEBRA PRECISION TESTS                                             ║\n";
    std::cout << "║  Comprehensive testing of LU, Cholesky, QR, SVD solvers and eigenvalue computation                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════════════════════════════════════════╝\n";
    
    PrecisionTestSuite suite("Linear Algebra", "Testing linear system solvers and matrix decompositions");
    
    // Run all tests
    Test_BasicLinearSystem(suite);
    Test_CholeskySPD(suite);
    Test_HilbertMatrix(suite);
    Test_EigenvalueAccuracy(suite);
    Test_QRDecomposition(suite);
    Test_SVDDecomposition(suite);
    Test_LUDeterminant(suite);
    Test_MatrixInversion(suite);
    Test_LeastSquares(suite);
    Test_HessenbergReduction(suite);
    
    // Print summary
    std::cout << "\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "                                       SUMMARY\n";
    std::cout << "════════════════════════════════════════════════════════════════════════════════════════════════════\n";
    
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Export CSV
    suite.exportCSV("precision_linalg_results.csv");
    std::cout << "\n  Results exported to: precision_linalg_results.csv\n";
}
