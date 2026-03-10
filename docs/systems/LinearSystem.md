# LinearSystem - Unified Linear Algebra Facade

`LinearSystem<T>` is MML's crown jewel facade for linear algebra, providing unified access to ALL linear system solving capabilities through a single, intuitive interface.

## Overview

```cpp
#include "systems/LinearSystem.h"
using namespace MML;

// Create and solve - it's that simple!
Matrix<Real> A(3, 3, {4, 1, 0, 1, 4, 1, 0, 1, 4});
Vector<Real> b({5, 6, 5});
LinearSystem<Real> sys(A, b);

Vector<Real> x = sys.Solve();  // Auto-selects best solver!
```

## Design Philosophy

1. **One Interface, All Analysis** - No hunting through multiple solver classes
2. **Smart Defaults** - Automatically selects optimal solver based on matrix properties
3. **Expert Control** - Force specific methods when needed for comparison/teaching
4. **Rich Diagnostics** - Condition numbers, residuals, stability assessment
5. **Efficient** - Lazy decomposition caching, factorization reuse
6. **Educational** - Clear API teaches linear algebra concepts

## Features

### Construction Options

```cpp
// System with single RHS: Ax = b
LinearSystem<Real> sys1(A, b);

// System with multiple RHS: AX = B (each column is a RHS)
LinearSystem<Real> sys2(A, B);

// Analysis only (no RHS)
LinearSystem<Real> sys3(A);
```

### Automatic Solver Selection

The `Solve()` method intelligently selects the best algorithm:

```cpp
auto x = sys.Solve();  // Returns solution vector
```

**Selection Logic:**
1. **Symmetric Positive Definite** → Cholesky (fastest, most stable)
2. **Symmetric** → LU with pivoting
3. **Triangular** → Back/Forward substitution (O(n²))
4. **Diagonally Dominant** → Gauss-Jordan
5. **General** → LU decomposition

### Specific Solvers

For comparison or when you know best:

```cpp
auto x1 = sys.SolveByGaussJordan();  // Classic elimination
auto x2 = sys.SolveByLU();           // LU decomposition with pivoting
auto x3 = sys.SolveByCholesky();     // For SPD matrices (fastest)
auto x4 = sys.SolveByQR();           // QR decomposition
auto x5 = sys.SolveBySVD();          // Most robust, handles rank deficiency
```

### Iterative Methods (For Large Sparse Systems)

```cpp
auto x = sys.SolveIterative(IterativeMethod::Auto);      // Auto-select
auto x = sys.SolveIterative(IterativeMethod::Jacobi);    // Parallelizable
auto x = sys.SolveIterative(IterativeMethod::GaussSeidel); // Faster convergence
auto x = sys.SolveIterative(IterativeMethod::SOR);       // Fastest, needs tuning
```

### Least Squares (Overdetermined Systems)

```cpp
// System with more equations than unknowns
Matrix<Real> A(100, 3, data);  // 100 observations, 3 unknowns
Vector<Real> b(100);

LinearSystem<Real> sys(A, b);
auto x = sys.SolveLeastSquares();  // Minimizes ||Ax - b||₂
```

## Solution Verification

Always verify your solutions!

```cpp
Vector<Real> x = sys.Solve();

// Quick check
Real residual = sys.ResidualNorm(x);        // ||Ax - b||
Real relResidual = sys.RelativeResidual(x); // ||Ax - b|| / ||b||

// Full verification
auto verify = sys.Verify(x);
std::cout << "Absolute residual: " << verify.absoluteResidual << "\n";
std::cout << "Relative residual: " << verify.relativeResidual << "\n";
std::cout << "Backward error: " << verify.backwardError << "\n";
std::cout << "Estimated forward error: " << verify.estimatedForwardError << "\n";
std::cout << "Is accurate: " << (verify.isAccurate ? "Yes" : "No") << "\n";
```

## Matrix Analysis

### Dimension Properties

```cpp
int m = sys.Rows();
int n = sys.Cols();

bool square = sys.IsSquare();
bool overdetermined = sys.IsOverdetermined();  // m > n
bool underdetermined = sys.IsUnderdetermined(); // m < n
```

### Structure Detection

```cpp
bool symmetric = sys.IsSymmetric();
bool spd = sys.IsPositiveDefinite();
bool diagDom = sys.IsDiagonallyDominant();
bool upper = sys.IsUpperTriangular();
bool lower = sys.IsLowerTriangular();
bool diag = sys.IsDiagonal();
```

### Numerical Properties

```cpp
Real det = sys.Determinant();
int rank = sys.Rank();
int nullity = sys.Nullity();  // = n - rank
Real cond = sys.ConditionNumber();

// Stability assessment
auto stability = sys.AssessStability();  // WellConditioned, Moderate, IllConditioned, Singular
int digitsLost = sys.ExpectedDigitsLost();  // log₁₀(cond)
```

## Cached Decompositions

Decompositions are computed lazily and cached:

```cpp
auto lu = sys.GetLU();        // LU decomposition
auto qr = sys.GetQR();        // QR decomposition
auto svd = sys.GetSVD();      // SVD decomposition
auto chol = sys.GetCholesky(); // Cholesky (throws if not SPD)
```

Each decomposition struct provides:
- Factor matrices (L, U, Q, R, U, S, V, etc.)
- Permutation information
- Numerical properties

## Eigenanalysis

```cpp
// Full eigensystem (handles complex eigenvalues)
auto eigen = sys.GetEigen();
for (const auto& ev : eigen.eigenvalues)
    std::cout << ev.real << " + " << ev.imag << "i\n";

// Real eigenvalues only
Vector<Real> eigs = sys.Eigenvalues();

// Symmetric matrices (faster, guaranteed real)
Vector<Real> eigs = sys.EigenvaluesSymmetric();

// Spectral radius (largest |eigenvalue|)
Real rho = sys.SpectralRadius();

// Check for complex eigenvalues
if (sys.HasComplexEigenvalues())
    std::cout << "Matrix has complex eigenvalues\n";
```

## Fundamental Subspaces

```cpp
Matrix<Real> null = sys.NullSpace();    // Orthonormal basis for null(A)
Matrix<Real> range = sys.ColumnSpace(); // Orthonormal basis for col(A)

// Matrix operations
Matrix<Real> inv = sys.Inverse();        // A⁻¹ (throws if singular)
Matrix<Real> pinv = sys.PseudoInverse(); // A⁺ (Moore-Penrose)
```

## Comprehensive Analysis

Get everything at once:

```cpp
auto analysis = sys.Analyze();

std::cout << "=== System Analysis ===" << "\n";
std::cout << "Dimensions: " << analysis.rows << " x " << analysis.cols << "\n";
std::cout << "Rank: " << analysis.rank << "\n";
std::cout << "Condition number: " << analysis.conditionNumber << "\n";
std::cout << "Stability: " << analysis.stabilityDescription << "\n";
std::cout << "Recommended solver: " << analysis.recommendedSolver << "\n";

// Solution exists?
if (analysis.hasUniqueSolution)
    std::cout << "System has unique solution\n";
else if (analysis.hasInfiniteSolutions)
    std::cout << "System has infinitely many solutions\n";
else
    std::cout << "System has no solution\n";
```

## Convenience Functions

For one-off operations without creating a LinearSystem object:

```cpp
// Quick solve
auto x = SolveLinearSystem(A, b);

// Quick least squares
auto x = SolveLeastSquares(A, b);

// Quick analysis
auto analysis = AnalyzeMatrix(A);
```

## Multiple Right-Hand Sides

Efficiently solve Ax₁ = b₁, Ax₂ = b₂, ... using the same factorization:

```cpp
Matrix<Real> B(n, numRHS);  // Each column is a RHS
LinearSystem<Real> sys(A, B);
Matrix<Real> X = sys.SolveMultiple();  // Each column is a solution
```

## Error Handling

LinearSystem throws meaningful exceptions:

```cpp
try {
    auto x = sys.Solve();
} catch (const MatrixDimensionError& e) {
    // Dimension mismatch
} catch (const SingularMatrixError& e) {
    // Matrix is singular
} catch (const ConvergenceError& e) {
    // Iterative method didn't converge
}
```

## Example: Comparing Solvers

```cpp
Matrix<Real> A(100, 100);  // Fill with your data
Vector<Real> b(100);

LinearSystem<Real> sys(A, b);

// Test different solvers
std::vector<std::pair<std::string, Vector<Real>>> solutions = {
    {"Auto", sys.Solve()},
    {"GaussJordan", sys.SolveByGaussJordan()},
    {"LU", sys.SolveByLU()},
    {"QR", sys.SolveByQR()},
    {"SVD", sys.SolveBySVD()}
};

// Compare residuals
for (const auto& [name, x] : solutions) {
    std::cout << name << ": residual = " << sys.ResidualNorm(x) << "\n";
}
```

## Performance Tips

1. **Use Auto** - `Solve()` picks optimal algorithm for your matrix
2. **SPD Matrices** - If you know your matrix is SPD, Cholesky is 2x faster than LU
3. **Multiple RHS** - Use `SolveMultiple()` instead of multiple `Solve()` calls
4. **Caching** - Decompositions are cached; reuse the LinearSystem object
5. **Large Sparse** - Use iterative methods for matrices > 1000×1000

## Implementation Details

LinearSystem wraps these existing MML components:

| Feature | Underlying Implementation |
|---------|---------------------------|
| Gauss-Jordan | `GaussJordanSolver` |
| LU | `LUSolver` |
| Cholesky | `CholeskySolver` |
| QR | `QRSolver` |
| SVD | `SVDecompositionSolver` |
| Jacobi/GS/SOR | `JacobiSolver`, `GaussSeidelSolver`, `SORSolver` |
| Eigenvalues | `EigenSolver`, `SymmMatEigenSolverJacobi` |
| Matrix utilities | `MatrixAlg` |

## See Also

- [LinAlgEqSolvers.h](../core/LinAlgEqSolvers.md) - Individual direct solvers
- [LinAlgEqSolvers_iterative.h](../core/LinAlgEqSolvers_iterative.md) - Iterative solvers
- [MatrixAlg.h](../algorithms/MatrixAlg.md) - Matrix utilities
- [EigenSystemSolvers.h](../algorithms/EigenSystemSolvers.md) - Eigenvalue solvers
