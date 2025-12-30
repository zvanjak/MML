# Matrix\<T\> - Dynamic-Size Matrix

**File**: `mml/base/Matrix.h`

Dynamic-size matrix with runtime dimensions. Template parameter `T` can be any numeric type.

## Runnable Examples

All code examples in this document are available as runnable code in `src/docs_demos/docs_demo_matrix.cpp`. Build and run `MML_DocsApp` to see them in action.

| Demo Function | Description |
|---------------|-------------|
| `Docs_Demo_Matrix_initializations()` | Creation patterns, factory methods, type aliases |
| `Docs_Demo_Basic_mat_operations()` | +, -, *, /, unary minus, compound assignment |
| `Docs_Demo_Matrix_Row_Col_Access()` | VectorFromRow, VectorFromColumn, Init methods |
| `Docs_Demo_Matrix_Submatrix()` | GetSubmatrix, GetLower, GetUpper |
| `Docs_Demo_Matrix_Properties()` | IsUnit, IsDiagonal, IsSymmetric, Trace |
| `Docs_Demo_Matrix_Norms()` | NormL1, NormL2, NormLInf, Utils::Det |
| `Docs_Demo_Matrix_Equality()` | ==, !=, IsEqualTo with tolerance |
| `Docs_Demo_Matrix_Vector_mul()` | Matrix-vector multiplication |
| `Docs_Demo_Matrix_Matrix_mul()` | Matrix-matrix multiplication |
| `Docs_Demo_Matrix_invert()` | GetInverse() |
| `Docs_Demo_Matrix_transpose()` | Transpose(), GetTranspose() |

## Quick Reference

```cpp
Matrix<double> A(3, 3);              // 3×3 zero matrix
Matrix<double> B(3, 3, 2.5);         // 3×3 filled with 2.5
Matrix<double> C{3, 3, {1,2,3,       // From initializer list
                        4,5,6,
                        7,8,9}};
Matrix<double> I = Matrix<double>::GetUnitMatrix(3);  // Identity matrix
```

## Construction

```cpp
// Empty and sized matrices
Matrix<double> empty;
Matrix<double> zeros(5, 3);          // 5×3 zero matrix
Matrix<double> filled(4, 4, 1.5);    // 4×4 filled with 1.5

// From initializer list (row-wise)
Matrix<double> A{3, 3, {1, 2, 3,
                        4, 5, 6,
                        7, 8, 9}};

// From std::vector<std::vector>
std::vector<std::vector<double>> data = {{1,2},{3,4}};
Matrix<double> B(data);

// Special matrices
auto I = Matrix<double>::GetUnitMatrix(5);     // 5×5 identity
auto diag = Matrix<double>::GetDiagonalMatrix(Vector<double>({1,2,3})); // Diagonal matrix
```

## Access & Properties

```cpp
Matrix<double> A{3, 3, {1,2,3,4,5,6,7,8,9}};

// Dimensions
int rows = A.RowNum();    // 3
int cols = A.ColNum();    // 3

// Element access: (row, col) - 0-based indexing
double val = A(1, 2);     // Get element at row 1, col 2
A(0, 0) = 10.0;           // Set element

// Row/column access
Vector<double> row = A.VectorFromRow(1);    // Get row 1
Vector<double> col = A.VectorFromColumn(2); // Get column 2
A.InitRowWithVector(0, Vector<double>({9,8,7}));   // Set row
A.InitColWithVector(1, Vector<double>({1,1,1}));   // Set column

// Diagonal
Vector<double> diag = A.GetDiagonal();  // Main diagonal
```

## Arithmetic Operations

```cpp
Matrix<double> A{2,2,{1,2,3,4}}, B{2,2,{5,6,7,8}};

// Addition/Subtraction
Matrix<double> sum = A + B;
Matrix<double> diff = A - B;
A += B;  // In-place addition

// Scalar multiplication
Matrix<double> scaled = A * 2.0;
Matrix<double> scaled2 = 3.0 * A;

// Matrix multiplication
Matrix<double> C = A * B;

// Matrix-Vector multiplication
Vector<double> v({1, 2});
Vector<double> result = A * v;
```

## Matrix Operations

```cpp
#include "core/MatrixUtils.h"  // For Utils::Det

Matrix<double> A{3,3,{1,2,3,4,5,6,7,8,9}};

// Transpose
Matrix<double> At = A.GetTranspose();  // Returns copy
A.Transpose();                          // In-place (square only)

// Determinant (free function)
double det = Utils::Det(A);

// Inverse (if exists)
Matrix<double> Ainv = A.GetInverse();

// Trace
double tr = A.Trace();

// Norms
double normL2 = A.NormL2();      // Frobenius norm
double normLinf = A.NormLInf();  // Max norm
double normL1 = A.NormL1();      // Sum of absolute values
```

## Submatrices & Reshaping

```cpp
Matrix<double> A{4,4, /* data */};

// Extract submatrix (startRow, startCol, numRows, numCols)
Matrix<double> sub = A.GetSubmatrix(0, 1, 2, 2);  // 2×2 block starting at (0,1)

// Extract triangular parts
Matrix<double> L = A.GetLower();  // Lower triangular (with diagonal)
Matrix<double> U = A.GetUpper();  // Upper triangular (with diagonal)

// Resize
A.Resize(5, 5);                    // Discard old data
A.Resize(5, 5, true);              // Preserve existing elements
```

## I/O and Printing

```cpp
Matrix<double> A{2,2,{1.5, 2.7, 3.14, 4.9}};

// Print to console
A.Print(std::cout, 8, 3);  // width=8, precision=3
// Output:
// Rows: 2 Cols: 2
// [     1.5,     2.7,  ]
// [    3.14,     4.9,  ]

// Stream operator
std::cout << A << std::endl;
```

## Examples

### Example 1: Solving Linear System
```cpp
#include "MML.h"
using namespace MML;

Matrix<double> A{3,3,{2,-1,0, -1,2,-1, 0,-1,2}};
Vector<double> b{1, 0, 1};

LUDecompositionSolver<double> solver(A);
Vector<double> x = solver.Solve(b);

std::cout << "Solution: " << x << std::endl;
std::cout << "Verification: " << (A * x) << std::endl;
```

### Example 2: Matrix Properties
```cpp
#include "core/MatrixUtils.h"

Matrix<double> A{3,3,{1,0,0, 0,2,0, 0,0,3}};

std::cout << "Determinant: " << Utils::Det(A) << std::endl; // 6
std::cout << "Trace: " << A.Trace() << std::endl;           // 6
std::cout << "L2 norm: " << A.NormL2() << std::endl;        // sqrt(14)
```

### Example 3: Matrix Decomposition
```cpp
#include "core/LinearAlgEqSolvers.h"

Matrix<double> A{3,3,{4,12,-16, 12,37,-43, -16,-43,98}};

// LU Decomposition
LUSolver<double> lu(A);
double det = lu.det();  // Determinant via LU

// QR Decomposition  
QRSolver<double> qr(A);
// Use qr.Solve(b) for solving linear systems
```

## Specialized Matrix Types

- **MatrixSym** - Symmetric matrices (efficient storage)
- **MatrixTridiag** - Tridiagonal matrices
- **MatrixBandDiag** - Band-diagonal matrices

See [Matrices.md](Matrices.md) for complete overview.

## Type Aliases

```cpp
typedef Matrix<Real>    MatrixReal;
typedef Matrix<Complex> MatrixComplex;
```

## See Also
- [Matrices.md](Matrices.md) - Matrix types overview
- [MatrixNM.md](MatrixNM.md) - Fixed-size matrices
- [Linear_equations_solvers.md](../core/Linear_equations_solvers.md) - Solving linear systems
- [Eigen_solvers.md](../algorithms/Eigen_solvers.md) - Eigenvalue computation
