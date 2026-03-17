///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_matrix.cpp
///  Description: Comprehensive Matrix demonstrations for Chapter 01
///               Shows construction, operations, and practical applications
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/Matrix/MatrixNM.h"
#include "mml/base/Polynom.h"

#include "mml/core/MatrixUtils.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        Matrix Construction Methods                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix_Construction()
{
    std::cout << "\n=== Matrix Construction ===\n\n";

    // Basic construction
    std::cout << "--- Basic Construction ---\n";
    Matrix<Real> a;                                    // Default (empty)
    Matrix<Real> b(2, 2);                              // 2x2, uninitialized
    Matrix<Real> c(2, 2, { 1.0, -2.0, 3.75, -1.0 });   // From initializer list
    Matrix<Real> d(c);                                 // Copy constructor
    Matrix<Real> e = c;                                // Assignment

    std::cout << "c (2x2 from list):\n" << c << "\n";

    // Different element types
    std::cout << "--- Different Types ---\n";
    MatrixInt mat_int(2, 2, { 1, 2, 3, 4 });
    std::cout << "Integer matrix:\n" << mat_int << "\n";

    MatrixComplex mat_cmplx(2, 2, { Complex(1,1), Complex(-1,2),
                                    Complex(2,-0.5), Complex(1,1) });
    std::cout << "Complex matrix:\n" << mat_cmplx << "\n";

    // Unit/Identity matrix
    std::cout << "--- Special Matrices ---\n";
    Matrix<Real> I3 = Matrix<Real>::Identity(3);
    std::cout << "3x3 Identity:\n" << I3 << "\n";

    // From vectors
    std::cout << "--- From Vectors ---\n";
    Vector<Real> vec_a({ 1, 2, 3, 4, 5 });
    Matrix<Real> diag = Utils::DiagonalMatrixFromVector(vec_a);
    Matrix<Real> row_mat = Utils::RowMatrixFromVector<Real>(vec_a);
    Matrix<Real> col_mat = Utils::ColumnMatrixFromVector<Real>(vec_a);

    std::cout << "Diagonal from vector:\n" << diag << "\n";
    std::cout << "Row matrix: " << row_mat << "\n";
    std::cout << "Column matrix:\n" << col_mat << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Submatrix Operations                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix_Submatrix()
{
    std::cout << "\n=== Submatrix Operations ===\n\n";

    Matrix<Real> M(5, 5, { 1,  2,  3,  4,  5,
                           6,  7,  8,  9,  10,
                           11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20,
                           21, 22, 23, 24, 25 });

    std::cout << "Original 5x5 matrix:\n" << M << "\n";

    // Extract submatrix (row, col, numRows, numCols)
    Matrix<Real> sub1 = M.GetSubmatrix(1, 1, 3, 3);  // 3x3 from center
    std::cout << "3x3 center submatrix (starting at row 1, col 1):\n" << sub1 << "\n";

    Matrix<Real> sub2 = M.GetSubmatrix(0, 0, 2, 5);  // First 2 rows
    std::cout << "First 2 rows:\n" << sub2 << "\n";

    // Transpose
    Matrix<Real> A(2, 3, { 1, 2, 3, 4, 5, 6 });
    std::cout << "A (2x3):\n" << A << "\n";
    std::cout << "A^T (3x2):\n" << A.transpose() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Matrix Arithmetic                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix_Arithmetic()
{
    std::cout << "\n=== Matrix Arithmetic ===\n\n";

    Matrix<Real> A(2, 2, { 1, 2, 3, 4 });
    Matrix<Real> B(2, 2, { 5, 6, 7, 8 });

    std::cout << "A =\n" << A;
    std::cout << "B =\n" << B << "\n";

    // Addition and subtraction
    std::cout << "--- Addition/Subtraction ---\n";
    std::cout << "A + B =\n" << (A + B);
    std::cout << "A - B =\n" << (A - B);

    // Scalar multiplication
    std::cout << "--- Scalar Operations ---\n";
    std::cout << "2.0 * A =\n" << (2.0 * A);
    std::cout << "A / 2.0 =\n" << (A / 2.0);

    // Matrix multiplication
    std::cout << "--- Matrix Multiplication ---\n";
    Matrix<Real> M1(3, 2, { 1.0, -1.0, 
                            2.0, 0.5,
                            1.5, 3.0 });
    Matrix<Real> M2(2, 4, {-2.5, 3.0, 1.0, -1.0, 
                           0.5,-2.0, 5.0, -3.0 });
    
    std::cout << "M1 (3x2):\n" << M1;
    std::cout << "M2 (2x4):\n" << M2;
    std::cout << "M1 * M2 (3x4):\n" << (M1 * M2);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Matrix-Vector Multiplication                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix_Vector_Mul()
{
    std::cout << "\n=== Matrix-Vector Multiplication ===\n\n";

    Vector<Real> v({ 1.0, 2.0 });
    Matrix<Real> M(2, 2, { 1.0, -1.0, 1.5, 3.0 });

    std::cout << "v = " << v << "\n";
    std::cout << "M =\n" << M;

    // v * M: row vector times matrix
    std::cout << "\nv * M = " << (v * M) << "  (row vector × matrix)\n";
    
    // M * v: matrix times column vector
    std::cout << "M * v = " << (M * v) << "  (matrix × column vector)\n";

    // Physics: Linear transformation
    std::cout << "\n--- Physics: 2D Rotation Matrix ---\n";
    double theta = Constants::PI / 4;  // 45 degrees
    Matrix<Real> R(2, 2, { std::cos(theta), -std::sin(theta),
                           std::sin(theta),  std::cos(theta) });
    Vector<Real> p({ 1.0, 0.0 });  // Point on X-axis

    std::cout << "Rotation matrix (45°):\n" << R;
    std::cout << "Original point: " << p << "\n";
    std::cout << "Rotated point:  " << (R * p) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Matrix Properties                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix_Properties()
{
    std::cout << "\n=== Matrix Properties ===\n\n";

    Matrix<Real> A(3, 3, { 3, 1, 5, 
                           3, 3, 1, 
                           4, 6, 4 });

    std::cout << "A =\n" << A;

    // Dimensions
    std::cout << "Rows: " << A.rows() << ", Cols: " << A.cols() << "\n";

    // Trace (sum of diagonal)
    std::cout << "Trace: " << A.trace() << " (3 + 3 + 4 = 10)\n";

    // Norms
    std::cout << "NormL1: " << A.NormL1() << "\n";

    // Check properties
    std::cout << "\n--- Matrix Type Checks ---\n";
    Matrix<Real> I = Matrix<Real>::Identity(3);
    Matrix<Real> sym(2, 2, { 1, 2, 2, 4 });

    std::cout << "I.isDiagonal(): " << (I.isDiagonal() ? "true" : "false") << "\n";
    std::cout << "sym is symmetric: [1,2; 2,4]\n" << sym;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Faddeev Algorithm - Advanced Example                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Faddeev_Algorithm()
{
    std::cout << "\n=== Faddeev Algorithm ===\n";
    std::cout << "(Computes characteristic polynomial, determinant, and inverse)\n\n";

    Matrix<Real> A(3, 3, { 3, 1, 5, 
                           3, 3, 1, 
                           4, 6, 4 });

    std::cout << "A =\n" << A;

    // Use the library function
    Real detValue;
    PolynomRealFunc charPoly;
    Matrix<Real> Ainv;

    Utils::FaddeevAlg(A, charPoly, detValue, Ainv);
    
    std::cout << "Characteristic polynomial: " << charPoly << "\n";
    std::cout << "Determinant: " << detValue << "\n";
    std::cout << "\nA^(-1) =\n" << Ainv;
    std::cout << "\nVerification A * A^(-1) =\n" << (A * Ainv);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Complex Matrices                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Complex_Matrices()
{
    std::cout << "\n=== Complex Matrices ===\n\n";

    // Quantum mechanics: Pauli matrices
    std::cout << "--- Pauli Spin Matrices ---\n";
    
    MatrixComplex sigma_x(2, 2, { Complex(0,0), Complex(1,0),
                                   Complex(1,0), Complex(0,0) });
    MatrixComplex sigma_y(2, 2, { Complex(0,0), Complex(0,-1),
                                   Complex(0,1), Complex(0,0) });
    MatrixComplex sigma_z(2, 2, { Complex(1,0), Complex(0,0),
                                   Complex(0,0), Complex(-1,0) });

    std::cout << "σx =\n" << sigma_x;
    std::cout << "σy =\n" << sigma_y;
    std::cout << "σz =\n" << sigma_z;

    // Verify: σx² = I
    std::cout << "σx² (should be I):\n" << (sigma_x * sigma_x);

    // Quantum gate: Hadamard
    std::cout << "\n--- Hadamard Gate ---\n";
    double s = 1.0 / std::sqrt(2.0);
    MatrixComplex H(2, 2, { Complex(s,0), Complex(s,0),
                            Complex(s,0), Complex(-s,0) });
    std::cout << "H = (1/√2) *\n" << H;

    // Apply to |0⟩ state
    Vector<Complex> ket_0({ Complex(1,0), Complex(0,0) });
    std::cout << "|0⟩ = " << ket_0 << "\n";
    std::cout << "H|0⟩ = " << (H * ket_0) << " (superposition)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // Computer Graphics: Scaling and rotation combined
    std::cout << "--- Graphics: Combined Transformation ---\n";
    
    double scale = 2.0;
    double angle = Constants::PI / 6;  // 30 degrees
    
    // Scale matrix
    Matrix<Real> S(2, 2, { scale, 0, 0, scale });
    
    // Rotation matrix
    Matrix<Real> R(2, 2, { std::cos(angle), -std::sin(angle),
                           std::sin(angle),  std::cos(angle) });
    
    // Combined: first scale, then rotate
    Matrix<Real> T = R * S;
    
    Vector<Real> point({ 1.0, 0.0 });
    std::cout << "Original point: " << point << "\n";
    std::cout << "After scale (2x): " << (S * point) << "\n";
    std::cout << "After scale + rotate (30°): " << (T * point) << "\n";

    // Linear system representation
    std::cout << "\n--- Linear Systems: Ax = b ---\n";
    Matrix<Real> A(2, 2, { 2, 1, 1, 3 });
    Vector<Real> b({ 5, 8 });
    
    std::cout << "System: A =\n" << A;
    std::cout << "b = " << b << "\n";
    std::cout << "(Solution methods covered in linear systems chapter)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Matrix()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                      MATRICES IN MML                          ****\n";
    std::cout << "****         Dynamic (Matrix<T>) and Fixed (MatrixNM<T,N,M>)       ****\n";
    std::cout << "***********************************************************************\n";

    Demo_Matrix_Construction();
    Demo_Matrix_Submatrix();
    Demo_Matrix_Arithmetic();
    Demo_Matrix_Vector_Mul();
    Demo_Matrix_Properties();
    Demo_Faddeev_Algorithm();
    Demo_Complex_Matrices();
    Demo_Matrix_Applications();
}