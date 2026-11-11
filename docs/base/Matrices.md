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
#include <mml/algorithms/MatrixAlg.h>
using namespace MML;

// Dynamic-size matrix
Matrix<Real> A{3, 3, {2, -1,  0,
                     -1,  2, -1,
                      0, -1,  2}};

// Fixed-size matrix (identity)
MatrixNM<Real, 3, 3> I = MatrixNM<Real, 3, 3>::Identity();

// Matrix operations
double det = MatrixAlg::Determinant(A);
Matrix<Real> At = A.Transposed(); // Returns transposed copy
Matrix<Real> Ainv = MatrixAlg::Inverse(A); // Returns inverse copy
```

## Specialized Matrix Types

MML provides efficient storage for special matrix structures:
- **MatrixSym** - Symmetric matrices (stores only lower triangle)
- **MatrixTridiag** - Tridiagonal matrices (three diagonals)
- **MatrixBandDiag** - Band-diagonal matrices (configurable bandwidth)

### Scale and failure contracts

`TridiagonalMatrix::IsEqualTo` uses an inclusive **absolute** tolerance on each
stored band entry, matching the comparison convention for dynamic, fixed-size
and symmetric matrices. Complex differences use their magnitude, not separate
real/imaginary component thresholds. Dimensions must match; negative or
nonfinite tolerances, or any nonfinite stored entry in either operand, return
false (including self-comparison). A zero tolerance requests exact equality
of finite entries. Integral comparisons avoid overflow, unsigned wraparound,
and rounding a large difference down to the tolerance.

The three-vector constructor's unused boundary slots (`below[0]` and
`above[dim-1]`) remain part of the stored-entry comparison, preserving existing
behavior; they are also checked for nonfinite values.

```cpp
const Real nan = std::numeric_limits<Real>::quiet_NaN();
const TridiagonalMatrix<Real> zero(1, {0}), invalid(1, {nan}), one(1, {1});
const bool finiteOnly = invalid.IsEqualTo(zero); // false
const bool invalidTolerance = zero.IsEqualTo(zero, nan); // false
const bool inclusiveBoundary = zero.IsEqualTo(one, 1); // true
```

This example runs in `Docs_Demo_Matrix_ScaleContracts()`, selected with
`MML_DocsApp --arithmetic-scale-contracts`.

`MatrixSym::NormFrobenius` and `BandDiagonalMatrix::NormFrobenius` use scaled
sum-of-squares accumulation. They do not square unscaled entries, so a finite,
representable norm is retained even when individual squared entries would
overflow or underflow. Symmetric off-diagonal entries are counted twice;
complex norms include both real and imaginary components.

`TridiagonalMatrix::Solve` normalizes the Thomas recurrence by the largest
active coefficient magnitude (largest component magnitude for complex entries).
Pivot rejection is relative to that scale, not
an absolute epsilon floor, and normalized elimination avoids products
underflowing solely because the whole system was scaled down.
Uniformly scaling a well-conditioned matrix and its right-hand side therefore
preserves the solution, including for subnormal coefficients.

This is still a nonpivoting solver: exact zero or numerically negligible
elimination pivots raise `SingularMatrixError`, even if a different row ordering
could solve the system. Incorrect vector dimensions raise `VectorDimensionError`;
nonfinite active coefficients or right-hand-side entries raise
`NumericInputError`. The caller-provided solution may be partially modified
if a later elimination pivot fails.

The executable counterpart is `Docs_Demo_Matrix_ScaleContracts()` in
[docs_demo_matrix_scale_contracts.cpp](../../src/docs_demos/base/docs_demo_matrix_scale_contracts.cpp).

## Related Documentation
- [MatrixAnalysisContracts.md](../algorithms/MatrixAnalysisContracts.md) - Canonical matrix-analysis architecture, numerical contracts, caching, and migration
- [Linear_equations_solvers.md](../core/Linear_equations_solvers.md) - Solving linear systems
- [Eigen_solvers.md](../algorithms/Eigen_solvers.md) - Computing eigenvalues and eigenvectors
- [Vectors.md](Vectors.md) - Vector types overview


Matrix inversion now belongs to algorithms, not base containers; see
[the inversion and formatting migration](../algorithms/MatrixInversionAndFormatting.md).
