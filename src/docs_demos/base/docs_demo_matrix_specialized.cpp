///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_matrix_specialized.cpp                                    ///
///  Description: Documentation examples for specialized matrix types                 ///
///               MatrixSym, TridiagonalMatrix, BandDiagonalMatrix                   ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#include "base/Matrix/MatrixTriDiag.h"
#include "base/Matrix/MatrixBandDiag.h"
#include "base/Vector/Vector.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                           Symmetric Matrix Construction                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixSym_Construction()
{
    std::cout << "\n=== Symmetric Matrix Construction ===\n\n";

    // Empty matrix
    MatrixSym<double> empty;
    std::cout << "Empty matrix: dim = " << empty.Dim() << std::endl;

    // Size-based construction (initialized to zero)
    MatrixSym<double> zeros(3);
    std::cout << "\nMatrixSym(3) - zeros:" << std::endl;
    zeros.Print(std::cout, 8, 3);

    // Lower triangular elements (row-wise)
    // For a 3x3 symmetric matrix, we store:
    // [a00, a10, a11, a20, a21, a22]
    MatrixSym<double> S(3, {1.0,           // row 0: a00
                            2.0, 3.0,      // row 1: a10, a11
                            4.0, 5.0, 6.0  // row 2: a20, a21, a22
                           });
    std::cout << "\nMatrixSym from lower triangular {1, 2, 3, 4, 5, 6}:" << std::endl;
    S.Print(std::cout, 8, 3);
    std::cout << "  (Notice: S(0,2) = S(2,0) = 4.0 due to symmetry)" << std::endl;

    // Uniform value construction
    MatrixSym<double> uniform(3, 7.5);
    std::cout << "\nMatrixSym(3, 7.5) - uniform value:" << std::endl;
    uniform.Print(std::cout, 8, 3);

    // Unit (identity) matrix
    auto identity = MatrixSym<double>::Identity(4);
    std::cout << "\nMatrixSym::Identity(4):" << std::endl;
    identity.Print(std::cout, 8, 3);

    // From full matrix (symmetrizes via (A + A^T)/2)
    Matrix<double> M(3, 3, {1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 9.0});
    std::cout << "\nOriginal non-symmetric matrix M:" << std::endl;
    M.Print(std::cout, 8, 3);
    
    auto fromFull = MatrixSym<double>::FromFullMatrix(M);
    std::cout << "MatrixSym::FromFullMatrix(M) - symmetrized:" << std::endl;
    fromFull.Print(std::cout, 8, 3);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Symmetric Matrix Properties                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixSym_Properties()
{
    std::cout << "\n=== Symmetric Matrix Properties ===\n\n";

    MatrixSym<double> S(3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    std::cout << "Symmetric matrix S:" << std::endl;
    S.Print(std::cout, 8, 3);

    // Dimensions
    std::cout << "\nDimensions:" << std::endl;
    std::cout << "  Dim():    " << S.Dim() << std::endl;
    std::cout << "  RowNum(): " << S.rows() << std::endl;
    std::cout << "  ColNum(): " << S.cols() << std::endl;
    std::cout << "  empty():  " << (S.empty() ? "true" : "false") << std::endl;

    // Storage efficiency
    size_t fullStorage = S.Dim() * S.Dim();
    size_t actualStorage = S.size();  // Uses compact storage
    std::cout << "\nStorage efficiency:" << std::endl;
    std::cout << "  Full matrix storage:    " << fullStorage << " elements" << std::endl;
    std::cout << "  MatrixSym storage:      " << actualStorage << " elements" << std::endl;
    std::cout << "  Savings:                " << (100.0 * (fullStorage - actualStorage) / fullStorage) << "%" << std::endl;

    // Trace
    std::cout << "\nTrace (sum of diagonal): " << S.Trace() << " (1 + 3 + 6 = 10)" << std::endl;

    // Norms
    std::cout << "\nMatrix norms:" << std::endl;
    std::cout << "  Frobenius norm: " << std::fixed << std::setprecision(4) << S.NormFrobenius() << std::endl;
    std::cout << "  Infinity norm:  " << S.NormInf() << std::endl;
    std::cout << "  1-norm:         " << S.Norm1() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Symmetric Matrix Element Access                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixSym_Access()
{
    std::cout << "\n=== Symmetric Matrix Element Access ===\n\n";

    MatrixSym<double> S(3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    std::cout << "Symmetric matrix S:" << std::endl;
    S.Print(std::cout, 8, 3);

    // Element access - symmetry is automatic
    std::cout << "\nElement access (symmetry is automatic):" << std::endl;
    std::cout << "  S(0,0) = " << S(0,0) << std::endl;
    std::cout << "  S(1,0) = " << S(1,0) << " = S(0,1) = " << S(0,1) << std::endl;
    std::cout << "  S(2,0) = " << S(2,0) << " = S(0,2) = " << S(0,2) << std::endl;
    std::cout << "  S(2,1) = " << S(2,1) << " = S(1,2) = " << S(1,2) << std::endl;

    // Modifying elements - symmetry is maintained
    std::cout << "\nModifying S(0,2) = 99.0:" << std::endl;
    S(0,2) = 99.0;
    S.Print(std::cout, 8, 3);
    std::cout << "  Notice: Both S(0,2) and S(2,0) are now 99.0" << std::endl;

    // Convert to full matrix
    std::cout << "\nConvert to full Matrix:" << std::endl;
    Matrix<double> full = S.GetAsMatrix();
    std::cout << "  S.GetAsMatrix() returned " << full.rows() << "x" << full.cols() << " matrix" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Symmetric Matrix Arithmetic                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixSym_Arithmetic()
{
    std::cout << "\n=== Symmetric Matrix Arithmetic ===\n\n";

    MatrixSym<double> A(3, {1.0, 0.0, 2.0, 0.0, 0.0, 3.0});
    MatrixSym<double> B(3, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0});

    std::cout << "Matrix A:" << std::endl;
    A.Print(std::cout, 8, 3);
    std::cout << "Matrix B:" << std::endl;
    B.Print(std::cout, 8, 3);

    // Addition (result is symmetric)
    MatrixSym<double> sum = A + B;
    std::cout << "A + B (symmetric result):" << std::endl;
    sum.Print(std::cout, 8, 3);

    // Subtraction
    MatrixSym<double> diff = A - B;
    std::cout << "A - B (symmetric result):" << std::endl;
    diff.Print(std::cout, 8, 3);

    // Scalar multiplication
    MatrixSym<double> scaled = A * 2.0;
    std::cout << "A * 2.0:" << std::endl;
    scaled.Print(std::cout, 8, 3);

    // Matrix-vector multiplication
    Vector<double> v({1.0, 2.0, 3.0});
    Vector<double> Av = A * v;
    std::cout << "A * [1, 2, 3] = " << Av << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Tridiagonal Matrix                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixTriDiag()
{
    std::cout << "\n=== Tridiagonal Matrix ===\n\n";

    std::cout << "Tridiagonal matrices have non-zero elements only on:\n";
    std::cout << "  - Main diagonal\n";
    std::cout << "  - One diagonal above (super-diagonal)\n";
    std::cout << "  - One diagonal below (sub-diagonal)\n\n";

    // Create 4x4 tridiagonal matrix using initializer list
    // Format: d0, c0, a1, d1, c1, a2, d2, c2, a3, d3
    // where: a = below, d = main, c = above
    TridiagonalMatrix<double> T(4, {
        2.0, -1.0,          // row 0: d0, c0
        -1.0, 2.0, -1.0,    // row 1: a1, d1, c1
        -1.0, 2.0, -1.0,    // row 2: a2, d2, c2
        -1.0, 2.0           // row 3: a3, d3
    });

    std::cout << "Tridiagonal matrix T (4x4, typical Poisson-like):" << std::endl;
    T.Print(std::cout, 8, 3);

    // Storage efficiency
    int n = T.RowNum();
    size_t fullStorage = n * n;
    size_t triStorage = 3 * n - 2;  // Main + above + below
    std::cout << "\nStorage comparison:" << std::endl;
    std::cout << "  Full matrix:       " << fullStorage << " elements" << std::endl;
    std::cout << "  Tridiagonal:       " << triStorage << " elements" << std::endl;
    std::cout << "  Savings:           " << (100.0 * (fullStorage - triStorage) / fullStorage) << "%" << std::endl;

    // Element access
    std::cout << "\nElement access:" << std::endl;
    std::cout << "  T(0,0) = " << T(0,0) << " (main diagonal)" << std::endl;
    std::cout << "  T(0,1) = " << T(0,1) << " (super-diagonal)" << std::endl;
    std::cout << "  T(1,0) = " << T(1,0) << " (sub-diagonal)" << std::endl;
    std::cout << "  T(0,2) = " << T(0,2) << " (zero - outside band)" << std::endl;

    // Solve tridiagonal system using Thomas algorithm
    Vector<double> rhs({1.0, 2.0, 3.0, 4.0});
    std::cout << "\nSolving T * x = rhs where rhs = " << rhs << std::endl;
    
    Vector<double> solution = T.Solve(rhs);
    std::cout << "Solution x = " << solution << std::endl;

    std::cout << "\nNote: Tridiagonal solve uses Thomas algorithm (O(n) complexity)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Band Diagonal Matrix                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixBandDiag()
{
    std::cout << "\n=== Band Diagonal Matrix ===\n\n";

    std::cout << "Band diagonal matrices have bandwidth parameters:" << std::endl;
    std::cout << "  - m1: number of sub-diagonals (below main)" << std::endl;
    std::cout << "  - m2: number of super-diagonals (above main)" << std::endl;
    std::cout << "  - Total bandwidth: m1 + m2 + 1" << std::endl;
    std::cout << std::endl;

    // Create 5x5 matrix with 1 sub-diagonal and 2 super-diagonals
    int n = 5, m1 = 1, m2 = 2;
    int bandwidth = m1 + m2 + 1;  // = 4
    
    // Create the compact data storage matrix (n x bandwidth)
    Matrix<Real> data(n, bandwidth);
    
    // Fill the band data
    // Column m1 = main diagonal, columns 0..m1-1 = subdiagonals, columns m1+1..m1+m2 = superdiagonals
    for (int i = 0; i < n; ++i) {
        data[i][m1] = 4.0;  // Main diagonal
    }
    for (int i = 1; i < n; ++i) {
        data[i][m1 - 1 + (i-1 - (i-1))] = 1.0;  // Sub-diagonal (simplified: element at data[i][0])
    }
    // Simpler approach: set values directly in compact storage
    // data[i][j-i+m1] corresponds to matrix element (i,j)
    data[0][m1] = 4.0; data[0][m1+1] = 2.0; data[0][m1+2] = 0.5;  // row 0
    data[1][m1-1] = 1.0; data[1][m1] = 4.0; data[1][m1+1] = 2.0; data[1][m1+2] = 0.5;  // row 1
    data[2][m1-1] = 1.0; data[2][m1] = 4.0; data[2][m1+1] = 2.0; data[2][m1+2] = 0.5;  // row 2
    data[3][m1-1] = 1.0; data[3][m1] = 4.0; data[3][m1+1] = 2.0;  // row 3
    data[4][m1-1] = 1.0; data[4][m1] = 4.0;  // row 4
    
    // Create the band diagonal matrix
    BandDiagonalMatrix B(n, m1, m2, data);

    std::cout << "Band matrix B (5x5, m1=1, m2=2):" << std::endl;
    std::cout << "Full matrix representation:" << std::endl;
    Matrix<Real> fullB = B.ToFullMatrix();
    fullB.Print(std::cout, 8, 3);

    // Storage efficiency
    size_t fullStorage = n * n;
    size_t bandStorage = n * (m1 + m2 + 1);
    std::cout << "\nStorage comparison:" << std::endl;
    std::cout << "  Full matrix:  " << fullStorage << " elements" << std::endl;
    std::cout << "  Band storage: " << bandStorage << " elements" << std::endl;
    std::cout << "  Savings:      " << (100.0 * (fullStorage - bandStorage) / fullStorage) << "%" << std::endl;

    // Matrix-vector multiplication
    Vector<Real> x({1.0, 2.0, 3.0, 4.0, 5.0});
    std::cout << "\nMatrix-vector multiplication: B * x" << std::endl;
    std::cout << "x = " << x << std::endl;
    
    Vector<Real> result = B * x;
    std::cout << "B * x = " << result << std::endl;

    // Demonstrate band queries
    std::cout << "\nBand queries:" << std::endl;
    std::cout << "  Dimension: " << B.GetDimension() << std::endl;
    std::cout << "  Lower bandwidth (m1): " << B.GetLowerBandwidth() << std::endl;
    std::cout << "  Upper bandwidth (m2): " << B.GetUpperBandwidth() << std::endl;
    std::cout << "  IsInBand(0,4): " << (B.IsInBand(0, 4) ? "yes" : "no") << std::endl;
    std::cout << "  IsInBand(4,0): " << (B.IsInBand(4, 0) ? "yes" : "no") << std::endl;

    std::cout << "\nNote: Band matrix operations exploit sparsity for efficiency." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Use Cases and Applications                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_SpecializedMatrix_Applications()
{
    std::cout << "\n=== Applications of Specialized Matrices ===\n\n";

    std::cout << "SYMMETRIC MATRICES (MatrixSym):" << std::endl;
    std::cout << "  - Covariance matrices in statistics" << std::endl;
    std::cout << "  - Hessian matrices in optimization" << std::endl;
    std::cout << "  - Mass/stiffness matrices in FEM" << std::endl;
    std::cout << "  - Graph adjacency matrices (undirected graphs)" << std::endl;

    std::cout << "\nTRIDIAGONAL MATRICES (TridiagonalMatrix):" << std::endl;
    std::cout << "  - 1D finite difference discretization" << std::endl;
    std::cout << "  - Cubic spline interpolation" << std::endl;
    std::cout << "  - Thomas algorithm for O(n) solving" << std::endl;
    std::cout << "  - Heat equation implicit schemes" << std::endl;

    std::cout << "\nBAND DIAGONAL MATRICES (BandDiagonalMatrix):" << std::endl;
    std::cout << "  - Higher-order finite differences" << std::endl;
    std::cout << "  - 2D/3D problems with structured grids" << std::endl;
    std::cout << "  - B-spline collocation methods" << std::endl;
    std::cout << "  - General sparse systems with banded structure" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Main Demo Runner                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MatrixSpecialized()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "       MML Documentation Demo: Specialized Matrices" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    Demo_MatrixSym_Construction();
    Demo_MatrixSym_Properties();
    Demo_MatrixSym_Access();
    Demo_MatrixSym_Arithmetic();
    Demo_MatrixTriDiag();
    Demo_MatrixBandDiag();
    Demo_SpecializedMatrix_Applications();

    std::cout << "\n" << std::string(70, '=') << std::endl;
}
