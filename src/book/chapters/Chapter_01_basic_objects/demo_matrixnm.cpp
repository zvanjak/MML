///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_matrixnm.cpp
///  Description: Comprehensive MatrixNM demonstrations for Chapter 01
///               Shows compile-time fixed-size matrices with full type safety
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Matrix/MatrixNM.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/BaseUtils.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        MatrixNM Overview                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_Overview()
{
    std::cout << "\n=== MatrixNM<T, N, M> Overview ===\n\n";
    
    std::cout << "MatrixNM provides compile-time fixed-size matrices:\n";
    std::cout << "  - Size known at compile time (N rows, M columns)\n";
    std::cout << "  - Stack allocation (no heap, faster)\n";
    std::cout << "  - Dimension mismatches caught at compile time\n";
    std::cout << "  - Ideal for small matrices (3x3, 4x4 transforms)\n";
    std::cout << "\nType aliases for common sizes:\n";
    std::cout << "  Mat22 = MatrixNM<Real, 2, 2>\n";
    std::cout << "  Mat33 = MatrixNM<Real, 3, 3>\n";
    std::cout << "  Mat44 = MatrixNM<Real, 4, 4>\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Construction Methods                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_Construction()
{
    std::cout << "\n=== MatrixNM Construction ===\n\n";

    // Default construction (zeros)
    std::cout << "--- Default and Value Initialization ---\n";
    MatrixNM<Real, 2, 2> zero;  // All zeros
    MatrixNM<Real, 2, 2> ones(1.0);  // All ones
    
    std::cout << "Default (zeros):\n" << zero << "\n";
    std::cout << "All 1.0:\n" << ones << "\n";

    // Initializer list (row-major)
    std::cout << "--- Initializer List (Row-Major) ---\n";
    MatrixNM<Real, 2, 3> m1({ 1, 2, 3,    // Row 0
                              4, 5, 6 }); // Row 1
    std::cout << "2x3 from initializer list:\n" << m1 << "\n";

    // Identity matrix
    std::cout << "--- Identity Matrices ---\n";
    auto I2 = MatrixNM<Real, 2, 2>::Identity();
    auto I3 = MatrixNM<Real, 3, 3>::Identity();
    auto I4 = MatrixNM<Real, 4, 4>::Identity();
    
    std::cout << "2x2 Identity:\n" << I2 << "\n";
    std::cout << "3x3 Identity:\n" << I3 << "\n";

    // Copy and assignment
    std::cout << "--- Copy Construction ---\n";
    MatrixNM<Real, 2, 2> original({1, 2, 3, 4});
    MatrixNM<Real, 2, 2> copy1(original);  // Copy constructor
    MatrixNM<Real, 2, 2> copy2 = original; // Copy assignment
    
    std::cout << "Original:\n" << original << "\n";
    std::cout << "Copy:\n" << copy1 << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector-Matrix Conversions                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_VectorConversions()
{
    std::cout << "\n=== Vector-Matrix Conversions ===\n\n";

    VectorN<Real, 3> v{1.0, 2.0, 3.0};
    std::cout << "Vector v = " << v << "\n\n";

    // Row and column matrices from vector
    std::cout << "--- Row and Column Matrices ---\n";
    auto rowMat = Utils::RowMatrixFromVector<Real, 3>(v);
    auto colMat = Utils::ColumnMatrixFromVector<Real, 3>(v);
    
    std::cout << "Row matrix (1x3):\n" << rowMat << "\n";
    std::cout << "Column matrix (3x1):\n" << colMat << "\n";

    // Extract vectors from matrix
    std::cout << "--- Extract Vectors from Matrix ---\n";
    MatrixNM<Real, 3, 3> M({
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    });
    
    std::cout << "Matrix M:\n" << M << "\n";
    
    auto row0 = M.VectorFromRow(0);
    auto col1 = M.VectorFromColumn(1);
    auto diag = M.VectorFromDiagonal();
    
    std::cout << "Row 0: " << row0 << "\n";
    std::cout << "Column 1: " << col1 << "\n";
    std::cout << "Diagonal: " << diag << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Arithmetic Operations                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_Arithmetic()
{
    std::cout << "\n=== MatrixNM Arithmetic ===\n\n";

    MatrixNM<Real, 2, 2> A({1, 2, 3, 4});
    MatrixNM<Real, 2, 2> B({5, 6, 7, 8});
    
    std::cout << "A =\n" << A << "\n";
    std::cout << "B =\n" << B << "\n";

    // Addition and subtraction
    std::cout << "--- Addition/Subtraction ---\n";
    std::cout << "A + B =\n" << (A + B) << "\n";
    std::cout << "A - B =\n" << (A - B) << "\n";

    // Scalar operations
    std::cout << "--- Scalar Multiplication ---\n";
    std::cout << "2.0 * A =\n" << (2.0 * A) << "\n";
    std::cout << "A / 2.0 =\n" << (A / 2.0) << "\n";

    // Matrix-matrix multiplication
    std::cout << "--- Matrix Multiplication ---\n";
    std::cout << "A * B =\n" << (A * B) << "\n";
    
    // Non-square multiplication
    MatrixNM<Real, 2, 3> M23({1, 2, 3, 4, 5, 6});
    MatrixNM<Real, 3, 2> M32({1, 2, 3, 4, 5, 6});
    
    std::cout << "M23 (2x3):\n" << M23 << "\n";
    std::cout << "M32 (3x2):\n" << M32 << "\n";
    std::cout << "M23 * M32 (2x2):\n" << (M23 * M32) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Matrix-Vector Multiplication                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_VectorMul()
{
    std::cout << "\n=== Matrix-Vector Multiplication ===\n\n";

    MatrixNM<Real, 2, 2> M({1, 2, 3, 4});
    VectorN<Real, 2> v{1, 2};
    
    std::cout << "M =\n" << M << "\n";
    std::cout << "v = " << v << "\n\n";

    // M * v (matrix times column vector)
    auto Mv = M * v;
    std::cout << "M * v = " << Mv << "  (column vector result)\n";

    // v * M (row vector times matrix)
    auto vM = v * M;
    std::cout << "v * M = " << vM << "  (row vector result)\n";

    // Practical: 2D transformation
    std::cout << "\n--- 2D Transformation Example ---\n";
    Real theta = Constants::PI / 4;  // 45 degrees
    MatrixNM<Real, 2, 2> Rot({
        std::cos(theta), -std::sin(theta),
        std::sin(theta),  std::cos(theta)
    });
    
    VectorN<Real, 2> point{1, 0};
    VectorN<Real, 2> rotated = Rot * point;
    
    std::cout << "Rotation matrix (45 deg):\n" << Rot << "\n";
    std::cout << "Point: " << point << "\n";
    std::cout << "Rotated: " << rotated << " (approx 45 deg CCW)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Transpose and Inverse                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_TransposeInverse()
{
    std::cout << "\n=== Transpose and Inverse ===\n\n";

    // Transpose
    std::cout << "--- Transpose ---\n";
    MatrixNM<Real, 2, 3> M({1, 2, 3, 4, 5, 6});
    auto MT = M.transpose();
    
    std::cout << "M (2x3):\n" << M << "\n";
    std::cout << "M^T (3x2):\n" << MT << "\n";

    // In-place transpose (only for square)
    MatrixNM<Real, 3, 3> S({1, 2, 3, 4, 5, 6, 7, 8, 9});
    std::cout << "Square matrix S:\n" << S << "\n";
    S.Transpose();
    std::cout << "After S.Transpose():\n" << S << "\n";

    // Inverse
    std::cout << "--- Inverse ---\n";
    MatrixNM<Real, 2, 2> A({4, 7, 2, 6});
    auto Ainv = A.GetInverse();
    
    std::cout << "A:\n" << A << "\n";
    std::cout << "A^(-1):\n" << Ainv << "\n";
    std::cout << "A * A^(-1) (should be I):\n" << (A * Ainv) << "\n";

    // 3x3 inverse
    MatrixNM<Real, 3, 3> B({
        1, 2, 3,
        0, 1, 4,
        5, 6, 0
    });
    auto Binv = B.GetInverse();
    
    std::cout << "B (3x3):\n" << B << "\n";
    std::cout << "B^(-1):\n" << Binv << "\n";
    std::cout << "B * B^(-1):\n" << (B * Binv) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Special Matrices and Properties                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_Special()
{
    std::cout << "\n=== Special Matrices and Properties ===\n\n";

    // Symmetric matrix
    std::cout << "--- Symmetric Matrix ---\n";
    MatrixNM<Real, 3, 3> sym({
        1, 2, 3,
        2, 4, 5,
        3, 5, 6
    });
    std::cout << "Symmetric matrix:\n" << sym << "\n";

    // Diagonal matrix
    std::cout << "--- Diagonal Matrix ---\n";
    MatrixNM<Real, 3, 3> diag;
    diag(0, 0) = 1;
    diag(1, 1) = 2;
    diag(2, 2) = 3;
    std::cout << "Diagonal matrix:\n" << diag << "\n";

    // Rotation matrices
    std::cout << "--- 3D Rotation Matrices ---\n";
    Real angle = Constants::PI / 6;  // 30 degrees
    
    // Rotation around Z-axis
    MatrixNM<Real, 3, 3> Rz({
        std::cos(angle), -std::sin(angle), 0,
        std::sin(angle),  std::cos(angle), 0,
        0,                0,               1
    });
    std::cout << "Rz (30 deg around Z):\n" << Rz << "\n";

    // Verify orthogonality: R^T * R = I
    auto RzT = Rz.transpose();
    std::cout << "Rz^T * Rz (should be I):\n" << (RzT * Rz) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // 4x4 homogeneous transform (graphics/robotics)
    std::cout << "--- 4x4 Homogeneous Transform ---\n";
    Real tx = 5, ty = 3, tz = 1;
    MatrixNM<Real, 4, 4> translate({
        1, 0, 0, tx,
        0, 1, 0, ty,
        0, 0, 1, tz,
        0, 0, 0, 1
    });
    
    std::cout << "Translation matrix (tx=5, ty=3, tz=1):\n" << translate << "\n";

    // Scale matrix
    Real sx = 2, sy = 2, sz = 2;
    MatrixNM<Real, 4, 4> scale({
        sx, 0,  0,  0,
        0,  sy, 0,  0,
        0,  0,  sz, 0,
        0,  0,  0,  1
    });
    std::cout << "Scale matrix (2x):\n" << scale << "\n";

    // Combined transform
    auto combined = translate * scale;
    std::cout << "Combined (translate * scale):\n" << combined << "\n";

    // Physics: Moment of inertia tensor
    std::cout << "--- Physics: Moment of Inertia Tensor ---\n";
    // For a uniform rectangular block
    Real m = 10.0;  // mass (kg)
    Real a = 2.0, b = 1.0, c = 0.5;  // dimensions (m)
    
    MatrixNM<Real, 3, 3> I({
        m*(b*b + c*c)/12, 0, 0,
        0, m*(a*a + c*c)/12, 0,
        0, 0, m*(a*a + b*b)/12
    });
    
    std::cout << "Inertia tensor for rectangular block:\n";
    std::cout << "  m=" << m << " kg, dimensions: " << a << "x" << b << "x" << c << " m\n";
    std::cout << I << "\n";

    // Stress tensor (engineering)
    std::cout << "--- Engineering: Stress Tensor ---\n";
    MatrixNM<Real, 3, 3> stress({
        100e6, 20e6,  0,      // sigma_xx, tau_xy, tau_xz
        20e6,  50e6,  0,      // tau_yx, sigma_yy, tau_yz  
        0,     0,     30e6    // tau_zx, tau_zy, sigma_zz
    });
    
    std::cout << "Stress tensor (Pa):\n" << stress << "\n";
    std::cout << "(Note: symmetric for equilibrium)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixNM()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                      MATRIXNM IN MML                          ****\n";
    std::cout << "****           Fixed-Size Matrices with Compile-Time Safety        ****\n";
    std::cout << "***********************************************************************\n";

    Demo_MatrixNM_Overview();
    Demo_MatrixNM_Construction();
    Demo_MatrixNM_VectorConversions();
    Demo_MatrixNM_Arithmetic();
    Demo_MatrixNM_VectorMul();
    Demo_MatrixNM_TransposeInverse();
    Demo_MatrixNM_Special();
    Demo_MatrixNM_Applications();
}


