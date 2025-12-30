# Matrices — Overview

Together with [vectors](Vectors.md), matrices are the basic building blocks of numerical computation. 

MML provides two complementary matrix types, depending on whether size is known at compile time:
- **`Matrix<T>`** — Dynamic-size matrix with runtime dimensions
- **`MatrixNM<T, N, M>`** — Fixed-size matrix with compile-time dimensions

Both matrix types have similar operations for initialization, element access, row/column manipulation, arithmetic, equality comparison, and I/O. Core linear algebra utilities include decompositions, inversions, and solvers.

## Highlights
- ✅ **Dual matrix families**: Dynamic (`Matrix<T>`) and fixed-size (`MatrixNM<T,N,M>`)
- ✅ **Performance optimized**: Stack allocation for fixed-size matrices
- ✅ **Rich operations**: All standard matrix operations and decompositions
- ✅ **Specialized types**: Symmetric, tridiagonal, band-diagonal matrices
- ✅ **Solver integration**: Works seamlessly with linear system and eigenvalue solvers

## Documentation
- **[Matrix.md](Matrix.md)** - Dynamic-size matrices (`Matrix<T>`)
- **[MatrixNM.md](MatrixNM.md)** - Fixed-size matrices (`MatrixNM<T, N, M>`)

### Runnable Examples

See individual documentation files for runnable demos:
- [Matrix.md](Matrix.md) - `Docs_Demo_Matrix()` functions
- [MatrixNM.md](MatrixNM.md) - `Docs_Demo_MatrixNM()` functions

## Choosing a Matrix Type

- **`Matrix<T>`**: Use for I/O operations, data pipelines, large matrices, or runtime-determined sizes
- **`MatrixNM<T,N,M>`**: Use for small matrices (2×2, 3×3, 4×4) with known dimensions for better performance

## Quick Example

```cpp
#include "core/MatrixUtils.h"
using namespace MML;

// Dynamic-size matrix
Matrix<Real> A{3, 3, {2, -1,  0,
                     -1,  2, -1,
                      0, -1,  2}};

// Fixed-size matrix (identity)
MatrixNM<Real, 3, 3> I = MatrixNM<Real, 3, 3>::GetUnitMatrix();

// Matrix operations
double det = Utils::Det(A);         // Determinant via Utils
Matrix<Real> At = A.GetTranspose(); // Returns transposed copy
Matrix<Real> Ainv = A.GetInverse(); // Returns inverse copy
```

## Specialized Matrix Types

MML provides efficient storage for special matrix structures:
- **MatrixSym** - Symmetric matrices (stores only upper triangle)
- **MatrixTridiag** - Tridiagonal matrices (three diagonals)
- **MatrixBandDiag** - Band-diagonal matrices (configurable bandwidth)

## Related Documentation
- [Linear_equations_solvers.md](../core/Linear_equations_solvers.md) - Solving linear systems
- [Eigen_solvers.md](../algorithms/Eigen_solvers.md) - Computing eigenvalues and eigenvectors
- [Vectors.md](Vectors.md) - Vector types overview

