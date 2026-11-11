# MML 2.0 Matrix Analysis Contract

This document defines the public matrix-analysis architecture and behavior in MML 2.0.

## Architecture

The dependency direction is:

```text
base Matrix types
    -> core decomposition and solver kernels
    -> MatrixAlg stateless analysis algorithms
    -> MatrixAnalyzer<Scalar> cached matrix facade
    -> LinearSystem<Real> RHS-aware solving facade
```

`MatrixAlg` owns public stateless matrix-analysis semantics. Numerical kernels remain in the core solver layer. `MatrixAnalyzer<Scalar>` owns an immutable matrix snapshot and adds lazy, configuration-aware caching. `LinearSystem<Real>` composes a real matrix analyzer with right-hand-side state, solving, consistency classification, recommendations, and verification.

Algorithms-layer types do not depend on `Systems::LinearSystem`.

`MatrixAnalyzer::Analyze()` returns numerical data only. Presentation is optional and lives in
`mml/tools/analyzers/MatrixAnalysisReports.h`; algorithms never include reporting tools.
`MatrixAnalysis::report` has been removed. Render a snapshot explicitly:

```cpp
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#include <mml/tools/analyzers/MatrixAnalysisReports.h>

Matrix<Real> matrix(2, 2, {2, 0, 0, 3});
MatrixAnalyzer<Real> matrixAnalyzer(matrix);
const auto matrixAnalysis = matrixAnalyzer.Analyze();
AnalysisReports::PrintMatrixAnalysis(std::cout, matrixAnalysis, "diagonal");
```

Reports do not request new factorizations. They accept real or complex analysis snapshots,
print optional definiteness classifications only when present, and preserve the caller's stream
formatting. `LinearSystem::Analyze().report` remains available and composes its matrix section
through this tools-side formatter.

Runnable counterparts: `Docs_Demo_AnalyzerReports` and `Demo_CachedAndComplexAnalysis` in
[function demos](../../src/docs_demos/algorithms/docs_demo_function_analyzer.cpp) and
[matrix demos](../../src/docs_demos/algorithms/docs_demo_matrix_alg.cpp).
Run `MML_DocsApp --analyzer-reports`.

## Choosing An Interface

Use `MatrixAlg` for independent stateless operations:

```cpp
#include <mml/algorithms/MatrixAlg.h>

Matrix<Real> A(3, 2, {1, 2, 2, 4, 3, 6});
int rank = MatrixAlg::Rank(A);
Matrix<Real> nullSpace = MatrixAlg::NullSpace(A);
Real norm = MatrixAlg::FrobeniusNorm(A);
```

Use `MatrixAnalyzer` when multiple requests should share factorizations:

```cpp
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>

MatrixAnalyzer<Real> analyzer(A);
const auto& svd = analyzer.SVDDecompose();
int rank = analyzer.Rank();
Real condition = analyzer.ConditionNumber();
const auto& spaces = analyzer.FundamentalSubspacesOf();
```

Use `LinearSystem<Real>` when a right-hand side affects the result:

```cpp
#include <mml/systems/LinearSystem.h>

Matrix<Real> coefficients(2, 2, {3, 1, 1, 2});
Vector<Real> rhs({9, 8});
Systems::LinearSystem<Real> system(coefficients, rhs);

Vector<Real> solution = system.Solve();
auto verification = system.Verify(solution);
auto analysis = system.Analyze();
```

## Naming And Ownership

Public matrix-analysis operations use PascalCase:

```cpp
MatrixAlg::Rank(matrix);
analyzer.Rank();
system.Rank();

analyzer.Rows();
analyzer.Cols();
analyzer.IsSquare();
```

Lowercase analysis aliases are not compatibility APIs. Container operations such as `Matrix::rows()` and distinct sparse-matrix or vector APIs are unaffected.

`MatrixOps` provides base-level comparison, similarity transformation, and Gram-Schmidt helpers. It does not own matrix-analysis queries or forwarding wrappers. Norms, structure predicates, determinant, rank, conditioning, subspaces, definiteness, and eigenanalysis belong to `MatrixAlg`.

`LinearSystem<Real>` forwards matrix-only operations to its `MatrixAnalyzer<Real>`. Its only separate factorization cache is the LU solver used to reuse a factorization across repeated right-hand-side solves.

## Scalar Model

```cpp
template<MMLScalar Scalar>
using MatrixMagnitude = decltype(std::abs(std::declval<Scalar>()));

template<MMLScalar Scalar>
using MatrixComplexScalar = std::conditional_t<
    MMLComplex<Scalar>,
    std::remove_cvref_t<Scalar>,
    std::complex<std::remove_cvref_t<Scalar>>>;
```

| Quantity | Return scalar |
|---|---|
| Matrix entries, trace, determinant | `Scalar` |
| Inverse, pseudoinverse, decomposition factors | `Matrix<Scalar>` |
| Norms, tolerances, singular values, condition numbers | `MatrixMagnitude<Scalar>` |
| Sparsity, residual magnitudes, spectral radius | `MatrixMagnitude<Scalar>` |
| Hermitian eigenvalues and definiteness pivots | `MatrixMagnitude<Scalar>` |
| General eigensystem eigenvalues and eigenvectors | `MatrixComplexScalar<Scalar>` |
| Rank, nullity, dimensions, expected digits lost | `int` |

For MML's configured `Complex`, `MatrixMagnitude<Complex>` is `Real`.

Complex support includes SVD, Hermitian Cholesky and definiteness, general complex eigenanalysis, pseudoinverse, fundamental subspaces, and cached `MatrixAnalyzer<Complex>` operations.

## Ownership And Caching

`MatrixAnalyzer<Scalar>` stores its matrix snapshot by value. Mutating the source matrix after construction cannot invalidate a cached result.

Expensive results are lazy. A cache key includes every argument that changes the result, including thresholds, tolerances, and maximum iterations. Requests with different configurations never share an entry. Successful results and failures are cached. Comprehensive analysis reuses cached decompositions instead of recomputing them.

## Tolerances And SVD Thresholds

Approximate structure predicates accept `MatrixComparisonTolerance<Magnitude>{absolute, relative}`. With

$$
s = \max(1, \max_{i,j}|a_{ij}|),
$$

a difference $d$ is treated as zero when

$$
|d| \le absolute + relative\,s.
$$

Symmetry applies this rule to $a_{ij}-a_{ji}$, Hermitian symmetry to $a_{ij}-\overline{a_{ji}}$, skew predicates to the corresponding sum, and triangular or diagonal predicates to entries required to be zero.

SVD-derived operations accept `std::optional<MatrixMagnitude<Scalar>> threshold = std::nullopt`. The automatic threshold is

$$
\tau = \max(m,n)\,\epsilon\,\sigma_{max}.
$$

The resolved threshold is reused for rank, nullity, condition number, pseudoinverse, and all four fundamental subspaces. Explicit thresholds are absolute, finite, and nonnegative.

## Dimensions And Failures

Shape queries are defined for every matrix, including empty matrices. A $0\times0$ matrix is square as a shape, but it is not positive definite.

| Operation | Singular matrix | Wrong shape or empty input |
|---|---|---|
| `Determinant` | returns `Scalar{0}` | throws `MatrixDimensionError` |
| `Rank`, `Nullity`, subspaces | returns numerical result | empty input throws |
| `ConditionNumber` | returns positive infinity | empty input throws |
| `Inverse` | throws `SingularMatrixError` | throws |
| `PseudoInverse` | returns numerical result | empty input throws |
| Requested decomposition | throws its specific MML error | throws |

A returned decomposition satisfies its invariants; result types do not use a `valid` flag produced by catch-all exception suppression.

`MatrixStability` values are `WellConditioned`, `ModeratelyConditioned`, `IllConditioned`, and `Singular`. Finite condition numbers below $10^4$ are well-conditioned, values below $10^8$ are moderately conditioned, and larger finite values are ill-conditioned. Infinite values are singular.

`ExpectedDigitsLost` returns $\lfloor\log_{10}(\max(\kappa,1))\rfloor$ for finite condition numbers and `std::nullopt` for singular/infinite values.

## Structure And Definiteness

- Symmetric means $A=A^T$.
- Skew-symmetric means $A=-A^T$.
- Hermitian means $A=A^*$.
- Skew-Hermitian means $A=-A^*$.
- Orthogonal uses transpose; unitary uses conjugate transpose.
- Diagonal dominance compares entry magnitudes and supports real and complex matrices.

`ClassifyDefiniteness` accepts only a symmetric real matrix or Hermitian complex matrix. A nonsymmetric or non-Hermitian input is rejected and is never silently replaced by its symmetric or Hermitian part.

The classifications are `PositiveDefinite`, `PositiveSemidefinite`, `NegativeDefinite`, `NegativeSemidefinite`, `Indefinite`, and `ZeroSemidefinite`. For `ZeroSemidefinite`, both semidefinite predicates return true while strict definiteness predicates return false.

## Decomposition Results

Canonical result types are:

- `LUDecomposition<Scalar>`: $PA=LU$, with permutation indices and determinant.
- `QRDecomposition<Scalar>`: economy QR for $m\ge n$.
- `SVDDecomposition<Scalar>`: full $U$ and $V$, magnitude-valued singular values, rank, and resolved threshold.
- `CholeskyDecomposition<Scalar>`: $A=LL^*$.
- `FundamentalSubspaces<Scalar>`: column, row, null, and left-null bases with one rank/threshold policy.
- `FaddeevLeVerrierResult<Scalar>`: characteristic coefficients, determinant, and optional inverse.

### Faddeev-LeVerrier numerical contract

`FaddeevLeVerrier` requires a nonempty square matrix with finite real or complex entries.
The coefficients are ordered from the monic leading coefficient to the constant term of
`det(lambda I - A)`. The determinant is the constant term with the appropriate parity sign.

The recurrence uses power-of-two normalization and rescales each coefficient according to
its degree. This avoids intermediate trace overflow for cases such as `1e154 * I` in double
precision, whose 2-by-2 determinant is the representable value `1e308`. Normalization that
loses an input component, nonfinite recurrence arithmetic, and coefficient overflow or
nonzero-to-zero rescaling underflow throw `MatrixNumericalError`; no partially corrupted
result is returned.

The optional inverse is computed independently using pivoted LU, not by dividing recurrence
matrices by a potentially cancelled constant coefficient. An absent inverse means LU
classified the matrix as singular (including its numerical pivot threshold). A present
inverse has finite entries and both `A * inverse` and `inverse * A` are finite and agree
entrywise with identity within `n * PrecisionValues<MatrixMagnitude<Scalar>>::MatrixIsEqualTolerance`.
An unrepresentable inverse or failed identity check throws `MatrixNumericalError`, rather
than masquerading as singularity or publishing an unusable matrix.

Power-of-two scaling does not remove the inherent cancellation sensitivity of the
Faddeev-LeVerrier coefficients. In particular, a rounded zero constant term is not a
singularity certificate and does not determine optional-inverse presence. For a determinant
without the full characteristic polynomial, prefer the LU-based `Determinant` API. These
checks do not promise coefficient accuracy for every conditioning or dynamic range.

Runnable counterpart: `Demo_FaddeevLeverrier` in
`src/docs_demos/core/docs_demo_matrix_analysis.cpp`, selected with
`MML_DocsApp --faddeev-leverrier`.

```cpp
const Real scale = 1e154; // Double-precision example.
const Matrix<Real> scaledIdentity{2, 2, {scale, 0, 0, scale}};
const auto scaled = MatrixAlg::FaddeevLeVerrier(scaledIdentity);
std::cout << "Determinant: " << scaled.determinant << '\n'; // Approximately 1e308.
if (scaled.inverse) {
    std::cout << "A * inverse:\n";
    (scaledIdentity * *scaled.inverse).Print(std::cout, 10, 6);
}
```

For complex matrices, orthogonality and subspaces use the Hermitian inner product:

$$
A=U\Sigma V^*, \qquad A^+=V\Sigma^+U^*.
$$

## Eigenanalysis

No API named `Eigenvalues` discards imaginary parts.

`EigensystemResult<Scalar>` contains complex-counterpart eigenvalues and eigenvectors plus convergence, iteration, residual, status, algorithm-name, and error-message metadata.

`SelfAdjointEigensystemResult<Scalar>` contains real magnitude-valued eigenvalues and `Matrix<Scalar>` unitary/orthogonal eigenvectors with the same diagnostics.

- General real matrices return `Vector<Complex>` and `Matrix<Complex>`.
- Symmetric real matrices use the symmetric Jacobi solver.
- Hermitian complex matrices use the Hermitian Jacobi solver and return real eigenvalues.
- General complex matrices use the complex shifted QR solver.
- Spectral radius is magnitude-valued.

Invalid dimensions and invalid self-adjoint preconditions throw. Structured eigensystem APIs
retain partial eigenvalues and eigenvectors on numerical nonconvergence, with
`converged == false` and failure diagnostics. The general real facade preserves the kernel's
status and error message rather than relabeling every failure as iteration exhaustion.
Imaginary components are never truncated.

Value-only `Eigenvalues`, `SymmetricEigenvalues`, `HermitianEigenvalues`, and `SpectralRadius`
helpers require both convergence and a usable status; otherwise they throw
`MatrixNumericalError` including the status and available solver message. The same rule applies
to cached `MatrixAnalyzer` helpers, including `HasComplexEigenvalues`. Use the structured
eigensystem APIs to inspect partial results instead. A failed cached request remains available
for inspection and does not prevent a separate request with a larger iteration budget.

Definiteness classification and the definiteness component of `Analyze()` also reject
nonconverged eigensolves rather than classifying partial eigenvalues.

## Matrix And System Reports

`MatrixAnalysis<Scalar>` contains dimensions, structure flags, sparsity, rank, nullity, optional determinant, condition number, stability, optional digits lost, optional definiteness, and formatted report text. It does not contain decomposition factors or eigenvectors. `Analyze()` performs at most one SVD and one definiteness eigensolve.

`Systems::SystemAnalysis<Real>` contains one `MatrixAnalysis<Real>` plus:

- one `SolutionStatus` per right-hand side,
- an optional `MultipleRHSAggregate`,
- an optional typed `LinearSolverRecommendation`,
- formatted report text.

Consistency uses $rank(A)=rank([A|b])$. A consistent system is unique when $rank(A)$ equals the number of unknowns and otherwise has infinitely many solutions. Without a right-hand side, solution statuses are empty and no solver is recommended.

Solver selection currently prioritizes rank-deficient, singular, and ill-conditioned systems for SVD, then square triangular substitution, tall QR, symmetric positive-definite Cholesky, and otherwise LU.

## Complexity

| Operation | Typical leading complexity |
|---|---|
| Structure predicates and norms | $O(mn)$ |
| LU or Cholesky for $n\times n$ | $O(n^3)$ |
| Economy QR with $m\ge n$ | $O(mn^2)$ |
| SVD | $O(mn\min(m,n))$ |
| Dense eigensystem | $O(n^3)$ |

An identical cached analyzer request avoids repeating the factorization. APIs returning matrices or vectors by value still incur normal copy costs.

## MML 2.0 Removed APIs

The migration removes rather than wraps:

- high-level analysis functions formerly exposed through `Utils`, including `IsOrthonormalColumns`,
- dense `Matrix::isDiagonal`, `Matrix::isDiagonallyDominant`, and `Matrix::isSymmetric`,
- dense `Matrix::NormL1`, `Matrix::NormL2`, and `Matrix::NormLInf`,
- the free `AnalyzeMatrix` wrapper,
- `MatrixAnalyzer::GetEigen` and `MatrixAnalyzer::EigenvaluesSymmetric`,
- `LinearSystem::GetEigen`, `LinearSystem::EigenvaluesSymmetric`, and lowercase analysis aliases,
- `Kernel` and `Range` matrix-analysis aliases,
- packed or imaginary-discarding public eigen APIs,
- duplicate analysis/decomposition result types and matrix-analysis caches.

Use these replacements:

| Removed surface | Canonical replacement |
|---|---|
| Legacy `Utils` analysis | `MatrixAlg` or `MatrixAnalyzer` |
| Dense matrix structure members | `MatrixAlg::IsDiagonal`, `IsSymmetric`, `IsDiagonallyDominant` |
| Dense matrix norm members | `MatrixAlg::OneNorm`, `FrobeniusNorm`, `InfinityNorm` |
| `Kernel`, `Range` | `NullSpace`, `ColumnSpace` |
| Legacy eigen APIs | `Eigensystem`, `Eigenvalues`, `SymmetricEigensystem`, `HermitianEigensystem` |
| `AnalyzeMatrix(matrix)` | `MatrixAnalyzer<Scalar>(matrix).Analyze()` |

The matrix $1$-norm is the maximum absolute column sum, and the matrix infinity norm is the maximum absolute row sum. These standard induced norms differ from the removed entrywise sum and maximum helpers.

## API Inventory

| Category | Canonical operations |
|---|---|
| Shape | `Rows`, `Cols`, `IsSquare`, `IsTall`, `IsWide` |
| Structure | `IsUpperHessenberg`, `IsUpperTriangular`, `IsLowerTriangular`, `IsDiagonal`, `IsSymmetric`, `IsSkewSymmetric`, `IsHermitian`, `IsSkewHermitian`, `IsDiagonallyDominant`, `Sparsity` |
| Scalar properties | `Trace`, `FrobeniusNorm`, `OneNorm`, `InfinityNorm`, `Determinant`, `ConditionNumber`, `ConditionNumber1`, `ConditionNumberInfinity`, `ExpectedDigitsLost`, `AssessStability` |
| Algebraic properties | `IsOrthogonal`, `IsUnitary`, `IsNilpotent`, `IsUnipotent`, definiteness classification and predicates |
| Factorizations | `LUDecompose`, `QRDecompose`, `SVDDecompose`, `CholeskyDecompose` |
| Rank and spaces | `Rank`, `Nullity`, `NullSpace`, `ColumnSpace`, `RowSpace`, `LeftNullSpace`, `FundamentalSubspacesOf` |
| Derived matrices | `Inverse`, `PseudoInverse`, `ReduceToHessenberg` |
| Eigenanalysis | `Eigensystem`, `SymmetricEigensystem`, `HermitianEigensystem`, `Eigenvalues`, `SymmetricEigenvalues`, `HermitianEigenvalues`, `SpectralRadius` |
| Comprehensive | `Analyze` returning `MatrixAnalysis<Scalar>` |

`RankGaussian` remains an explicitly named educational alternative. `FaddeevLeVerrier` remains a single result-returning algorithm. Similarity transformation, entry comparison, and Gram-Schmidt remain outside `MatrixAnalyzer` because they transform matrices or construct bases rather than report properties.
