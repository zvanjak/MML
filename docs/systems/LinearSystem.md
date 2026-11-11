# LinearSystem - Unified Linear Algebra Facade

`LinearSystem<Real>` combines coefficient-matrix analysis with right-hand-side-aware solving and diagnostics. Use `MatrixAlg` for stateless matrix queries and `MatrixAnalyzer<Scalar>` for cached matrix-only analysis.

## Overview

```cpp
#include <mml/systems/LinearSystem.h>
using namespace MML;
using namespace MML::Systems;

// Create and solve - it's that simple!
Matrix<Real> A(3, 3, {4, 1, 0, 1, 4, 1, 0, 1, 4});
Vector<Real> b({5, 6, 5});
LinearSystem<Real> sys(A, b);

Vector<Real> x = sys.Solve();  // Auto-selects best solver!
```

## Design Philosophy

1. **RHS-aware facade** - Solving, consistency, recommendations, and verification in one interface
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
1. **Rank-deficient, wide, singular, or ill-conditioned** → SVD
2. **Exactly triangular square** → Back/Forward substitution (O(n²))
3. **Tall full-column-rank** → QR least squares
4. **Exactly symmetric positive definite** → Cholesky
5. **General square** → LU with pivoting

`Analyze(threshold)` applies the same SVD threshold to matrix rank, per-RHS
consistency, stability, and the typed solver recommendation.

Structural solver selection uses zero absolute and relative comparison tolerances:
triangular substitution requires actual zeros in the unused triangle, and Cholesky
selection requires exactly matching off-diagonal pairs. Nonzero stored entries
are not dropped merely because they are small in the current units. Upper versus
lower substitution uses the same exact structural rule.

Public property queries such as `IsUpperTriangular()`, `IsLowerTriangular()` and
`IsSymmetric()` retain their approximate default tolerances. Consequently, an
approximate structural flag in `Analyze().matrix` can be true while the solver
recommendation deliberately chooses a general method for the full stored matrix.
Rank, conditioning and positive-definiteness decisions retain their existing
numerical thresholds; this is not an exact-arithmetic guarantee for the solve.

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

Verification requires finite matrix, RHS and solution entries, a matching solution
dimension and a finite positive tolerance. Nonfinite inputs raise `NumericInputError`;
a matrix-RHS instance cannot be verified through the single-vector interface.

Residual and denominator magnitudes are scaled before ratios are formed. Thus an
absolute residual or norm may exceed the scalar range while a relative residual
remains finite. A nonzero RHS always uses its actual norm, however small. For an
exactly zero RHS, the relative residual is zero when the residual is zero and
infinity otherwise. Residual terms spanning an unresolvable accumulation range
raise `MatrixNumericalError`, rather than silently certifying a lost contribution.
These guarantees apply to the stored `Real` entries: a scale that underflows to
zero changes the problem, while overflowed inputs are rejected rather than
reconstructed. An exactly zero matrix and RHS give zero residual for every finite
candidate; `isAccurate` is not a uniqueness test.

`backwardError` retains the existing diagnostic convention
`||Ax-b||_2 / (||A||_infinity * ||x||_2)`: an A-only mixed-norm proxy, **not**
the standard joint A/b normwise backward error. A zero denominator gives zero
for an exact residual and infinity otherwise. `estimatedForwardError` uses the
condition number times relative residual; an unavailable/nonfinite condition
number leaves this estimate unknown (NaN). `isAccurate` means only that the
finite relative residual is below the requested tolerance, not that forward
solution error was independently bounded.

### Dimension Properties

```cpp
int m = sys.Rows();
int n = sys.Cols();

bool square = sys.IsSquare();
bool overdetermined = sys.IsTall(); // m > n
bool underdetermined = sys.IsWide(); // m < n
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
auto stability = sys.AssessStability();  // WellConditioned, ModeratelyConditioned, IllConditioned, Singular
std::optional<int> digitsLost = sys.ExpectedDigitsLost(); // nullopt when singular
```

## Cached Decompositions

Decompositions are computed lazily and cached:

```cpp
const auto& lu = sys.LUDecompose();
const auto& qr = sys.QRDecompose();
const auto& svd = sys.SVDDecompose();
const auto& chol = sys.CholeskyDecompose(); // throws if not SPD
```

Each decomposition struct provides:
- Factor matrices (L, U, Q, R, U, S, V, etc.)
- Permutation information
- Numerical properties

## Eigenanalysis

```cpp
// Full eigensystem (handles complex eigenvalues)
const auto& eigen = sys.Eigensystem();
for (const auto& ev : eigen.eigenvalues)
    std::cout << ev.real() << " + " << ev.imag() << "i\n";

// General real matrices can have complex eigenvalues
Vector<Complex> eigs = sys.Eigenvalues();

// Symmetric matrices (faster, guaranteed real)
Vector<Real> eigs = sys.SymmetricEigenvalues();

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
std::cout << "Dimensions: " << analysis.matrix.rows << " x " << analysis.matrix.cols << "\n";
std::cout << "Rank: " << analysis.matrix.rank << "\n";
std::cout << "Condition number: " << analysis.matrix.conditionNumber << "\n";
std::cout << analysis.report;

for (SolutionStatus status : analysis.solutionStatuses) {
    if (status == SolutionStatus::Unique)
        std::cout << "Unique solution\n";
    else if (status == SolutionStatus::Infinite)
        std::cout << "Infinitely many solutions\n";
    else if (status == SolutionStatus::Inconsistent)
        std::cout << "No solution\n";
}
```

## Convenience Functions

For one-off operations without creating a LinearSystem object:

```cpp
// Quick solve
auto x = SolveLinearSystem(A, b);

// Quick least squares
auto x = SolveLeastSquares(A, b);

// Quick analysis
auto analysis = MatrixAnalyzer<Real>(A).Analyze();
```

## Multiple Right-Hand Sides

`SolveMultiple()` follows the same automatic selection as vector `Solve()`: exact
triangular substitution, Cholesky or LU for suitable square matrices, QR least
squares for tall full-column-rank matrices, and SVD for wide, rank-deficient or
ill-conditioned matrices. SVD returns minimum-norm least-squares solutions under
the existing singular-value threshold policy; inconsistent columns do not become
exact solutions merely because a result is returned.

Each column uses the same solver instance within a call, rather than refactorizing
per RHS. The general-square LU solver also remains cached across calls. Other
solving factorizations are constructed once per call; the separately exposed
matrix-analysis caches retain their existing behavior.
Nonfinite matrix-RHS entries raise `NumericInputError`. A nonfinite computed
solution raises `MatrixNumericalError`, rather than returning a usable-looking
solution matrix; other coefficient/factorization errors propagate normally.

```cpp
Matrix<Real> A(3, 2, {1, 0, 0, 1, 1, 1});
Matrix<Real> B(3, 2, {1, 1, 2, 2, 3, 4});
LinearSystem<Real> sys(A, B);
const Matrix<Real>& storedB = sys.getRHSMatrix();
Matrix<Real> X = sys.SolveMultiple();
// X columns: [1, 2] (exact) and [4/3, 7/3] (least squares).
```

For A with shape m-by-n and B with shape m-by-p, the result is n-by-p.
A matrix RHS with zero columns returns n-by-zero after the usual coefficient
matrix validation/solver selection; an empty coefficient matrix remains invalid.

RHS mode is determined by the constructor, not by the number of columns:
even a one-column B is a matrix RHS. `getRHS()` returns the owned vector RHS;
`getRHSMatrix()` returns the owned matrix RHS. Both return const references and
throw `InvalidStateError` for the wrong mode or when no RHS was supplied.
`GetRHS()` remains as a deprecated forwarding alias to `getRHS()`.

Every vector-solving entry point (including the specific and iterative solvers)
and vector verification (`ResidualNorm`, `RelativeResidual`, `Verify`) requires
a vector RHS. `SolveMultiple()` requires a matrix RHS. Wrong-mode calls throw
`InvalidStateError` before accessing an absent RHS or invoking a solver.
Matrix-only queries and RHS-aware `Analyze()` remain available in every mode.

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
