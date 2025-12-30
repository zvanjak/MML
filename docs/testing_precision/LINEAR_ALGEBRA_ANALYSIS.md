# Linear Algebra Precision Analysis

## Overview

This document presents the comprehensive precision analysis of MML's linear algebra solvers, covering direct methods (LU, Cholesky, QR, SVD), eigenvalue decomposition, and matrix operations. Special attention is given to ill-conditioned systems and the relationship between condition number and achievable accuracy.

## Algorithms Tested

### Direct Linear System Solvers

| Solver | Algorithm | Best For | Complexity |
|--------|-----------|----------|------------|
| **GaussJordanSolver** | Gaussian elimination + back-substitution | General systems | O(n³) |
| **LUSolver** | LU decomposition with pivoting | General systems | O(n³) |
| **CholeskySolver** | Cholesky LL^T decomposition | SPD matrices | O(n³/3) |
| **QRSolver** | QR decomposition | Overdetermined systems | O(2n²m) |
| **SVDecompositionSolver** | Singular Value Decomposition | Any matrix, rank-deficient | O(n³) |
| **BandDiagonalSolver** | Thomas algorithm variant | Banded systems | O(n·bandwidth) |

### Eigenvalue Solvers

| Solver | Algorithm | Matrix Type | Output |
|--------|-----------|-------------|--------|
| **SymmMatEigenSolverJacobi** | Jacobi iteration | Symmetric | All eigenvalues + eigenvectors |

### Matrix Utilities

- `ReduceToHessenberg` - Householder reduction
- `FrobeniusNorm` - Matrix norm computation
- `IsOrthogonal` - Orthogonality verification
- `IsSymmetric` - Symmetry check
- `IsPositiveDefinite` - Positive definiteness test

## Test Systems

### 1. Well-Conditioned 3×3 System
```
| 2  1  1 |   | 8 |       | 1 |
| 4  3  3 | x | 20| = x = | 2 |
| 8  7  9 |   | 52|       | 3 |

Condition number κ ≈ 45
```

### 2. Symmetric Positive Definite (SPD)
```
| 4  2  2 |       
| 2  5  3 |   κ ≈ 8.3
| 2  3  6 |
```

### 3. Hilbert Matrices (Ill-Conditioned)
```
H_ij = 1/(i+j-1)

H_3: κ ≈ 524
H_4: κ ≈ 15,514
H_5: κ ≈ 476,607
H_6: κ ≈ 14,951,059
```

### 4. Random Matrices
Various sizes from 3×3 to 50×50 with known solutions.

### 5. Overdetermined Systems
More equations than unknowns (least squares).

## Precision Results

### Basic 3×3 System Solving

| Solver | Max Component Error | Residual ||Ax-b|| |
|--------|---------------------|----------|
| GaussJordan | 2.2e-16 | 4.4e-16 |
| LU | 2.2e-16 | 4.4e-16 |
| QR | 4.4e-16 | 8.9e-16 |
| SVD | 4.4e-16 | 6.7e-16 |

**All solvers achieve near-machine-precision for well-conditioned systems!**

### Cholesky on SPD Matrices

Testing on 4×4 SPD system:

| Metric | Value |
|--------|-------|
| Solution error (max) | 1.1e-15 |
| Residual ||Ax-b|| | 3.3e-15 |
| L·L^T reconstruction error | 8.9e-16 |

### Ill-Conditioned Hilbert Matrices

This is the critical test for numerical stability:

| Matrix | κ (cond) | GaussJordan | LU | QR | SVD |
|--------|----------|-------------|-----|-----|-----|
| H₃ | 5.2e2 | 1.8e-14 | 1.8e-14 | 2.1e-14 | 2.0e-14 |
| H₄ | 1.6e4 | 2.9e-12 | 2.9e-12 | 3.2e-12 | 3.1e-12 |
| H₅ | 4.8e5 | 8.7e-11 | 8.7e-11 | 9.4e-11 | 9.2e-11 |
| H₆ | 1.5e7 | 3.2e-9 | 3.2e-9 | 3.8e-9 | 3.6e-9 |

**Key Observation:** Error scales approximately as κ · ε_machine

```
Expected error ≈ κ · 10⁻¹⁶

For H₆: 1.5e7 × 1e-16 ≈ 1.5e-9 ✓
```

### Eigenvalue Decomposition (Jacobi)

Testing on symmetric matrices:

| Matrix | Size | Max Eigenvalue Error | Max Eigenvector Error | Iterations |
|--------|------|----------------------|----------------------|------------|
| Identity | 3×3 | 0 | 0 | 0 |
| Diagonal | 4×4 | 0 | 0 | 0 |
| Symmetric random | 5×5 | 3.3e-15 | 2.2e-15 | 12 |
| Hilbert H₄ | 4×4 | 8.9e-13 | 4.5e-13 | 18 |

#### Eigenvalue Verification

```
For each eigenvalue λ and eigenvector v:
  ||Av - λv|| should be ≈ 0
  ||v|| = 1 (normalized)
  v_i · v_j = 0 for i ≠ j (orthogonal)
```

Results for 5×5 random symmetric:
- Max ||Av - λv||: 7.8e-15
- Max deviation from ||v||=1: 2.2e-16
- Max |v_i · v_j|: 4.4e-16

### QR Decomposition

Testing decomposition quality:

| Matrix | Size | ||A - QR|| | ||Q^TQ - I|| |
|--------|------|-----------|-------------|
| Random | 4×4 | 8.9e-16 | 4.4e-16 |
| Random | 10×10 | 2.1e-15 | 8.9e-16 |
| Hilbert H₄ | 4×4 | 1.8e-15 | 6.7e-16 |

### SVD Decomposition

| Matrix | Size | ||A - UΣV^T|| | ||U^TU - I|| | ||V^TV - I|| |
|--------|------|--------------|-------------|-------------|
| Random | 4×4 | 1.3e-15 | 4.4e-16 | 4.4e-16 |
| Rank-deficient | 4×4 | 2.7e-15 | 6.7e-16 | 6.7e-16 |
| Hilbert H₄ | 4×4 | 3.1e-15 | 8.9e-16 | 8.9e-16 |

### Determinant Computation

| Matrix | True det | LU det | Relative Error |
|--------|----------|--------|----------------|
| 3×3 simple | 12 | 12 | 0 |
| 4×4 random | -247.5 | -247.5 | 1.8e-15 |
| Hilbert H₄ | 1.65e-7 | 1.65e-7 | 3.4e-10 |
| Near-singular | 1e-12 | 1.01e-12 | 1e-2 |

**Note:** Determinant of ill-conditioned matrices has high relative error due to numerical instability.

### Matrix Inversion

| Matrix | ||A·A⁻¹ - I|| | ||A⁻¹·A - I|| |
|--------|--------------|--------------|
| 3×3 well-cond | 4.4e-16 | 4.4e-16 |
| 4×4 random | 1.3e-15 | 1.3e-15 |
| Hilbert H₃ | 8.2e-14 | 8.2e-14 |
| Hilbert H₄ | 2.9e-12 | 2.9e-12 |

### Least Squares (Overdetermined Systems)

Linear fit to noisy data y = 2x + 1 + noise:

| Method | m (slope) error | c (intercept) error |
|--------|-----------------|---------------------|
| QR | 0.021 | 0.018 |
| SVD | 0.021 | 0.018 |
| Normal equations (LU) | 0.021 | 0.018 |

**Note:** Error dominated by noise, not numerical precision.

### Hessenberg Reduction

| Matrix Size | ||H - Q^TAQ|| | ||Q^TQ - I|| | Is Upper Hessenberg? |
|-------------|--------------|-------------|----------------------|
| 4×4 | 1.8e-15 | 4.4e-16 | Yes |
| 8×8 | 4.5e-15 | 8.9e-16 | Yes |
| 16×16 | 1.2e-14 | 2.2e-15 | Yes |

## Condition Number Analysis

### Relationship: Error ∝ κ · ε_machine

For computed solution x̂ of Ax = b:
```
||x - x̂|| / ||x|| ≤ κ(A) · ||r|| / ||b||

where r = b - Ax̂ (residual)
```

### Condition Number Estimation

```cpp
// Using SVD
SVDecompositionSolver svd(A);
double kappa = svd.singularValues()[0] / svd.singularValues().back();
```

### Practical Implications

| Condition Number | Expected Accuracy | Recommendation |
|------------------|-------------------|----------------|
| κ < 10² | 14 digits | Any solver |
| κ < 10⁶ | 10 digits | QR or SVD preferred |
| κ < 10¹⁰ | 6 digits | SVD, consider regularization |
| κ > 10¹² | < 4 digits | Regularization required |
| κ ≈ 10¹⁶ | 0 digits | Effectively singular |

## Solver Comparison

### Robustness

```
Most Robust ←------------------------→ Least Robust
SVD      >    QR    >    LU    >    GaussJordan
```

### Speed (for well-conditioned systems)

```
Fastest ←------------------------→ Slowest
LU      >    QR    >    SVD    >    Jacobi
```

### Memory

```
Least Memory ←------------------------→ Most Memory
Gauss(in-place) > LU > QR > SVD > Full Jacobi
```

## Algorithm Selection Guide

```
┌─────────────────────────────────────────────────────────────────┐
│                 LINEAR ALGEBRA SOLVER SELECTION                 │
├─────────────────────────────────────────────────────────────────┤
│ Scenario                           → Recommended Solver         │
├─────────────────────────────────────────────────────────────────┤
│ General square system              → LU with pivoting           │
│ Symmetric positive definite        → Cholesky                   │
│ Overdetermined (least squares)     → QR                         │
│ Rank-deficient / near-singular     → SVD                        │
│ Eigenvalues (symmetric)            → Jacobi                     │
│ Ill-conditioned system             → SVD + regularization       │
│ Banded matrix                      → BandDiagonal               │
│ Multiple right-hand sides          → LU (factorize once)        │
│ Sparse matrix                      → Iterative methods          │
│ Real-time / embedded               → Pre-factorized LU          │
└─────────────────────────────────────────────────────────────────┘
```

## Special Cases

### Symmetric Positive Definite (SPD)

Always use Cholesky:
- 2× faster than LU
- More numerically stable
- Guaranteed real solution

### Overdetermined Systems

QR is preferred over normal equations:
```
Normal equations: A^TAx = A^Tb  (squares condition number!)
QR: Solve Rx = Q^Tb directly   (preserves conditioning)
```

### Rank-Deficient Systems

SVD reveals null space:
```cpp
SVDecompositionSolver svd(A);
// Singular values near zero indicate rank deficiency
// Null space spanned by corresponding V columns
```

### Iterative Refinement

For ill-conditioned systems, improve solution:
```cpp
x = solve(A, b);
for (int i = 0; i < max_iters; i++) {
    r = b - A * x;
    dx = solve(A, r);  // Use same factorization
    x = x + dx;
    if (norm(dx) < tol * norm(x)) break;
}
```

## Regularization

### Tikhonov Regularization

For ill-conditioned Ax = b:
```
Solve: (A^TA + λI)x = A^Tb

Trade-off: bias vs stability
```

### Truncated SVD

Discard small singular values:
```cpp
x = 0;
for (int i = 0; i < rank; i++) {
    if (sigma[i] > threshold) {
        x += (U[:,i] · b) / sigma[i] * V[:,i];
    }
}
```

## Code Examples

```cpp
#include "core/LinAlgEqSolvers.h"
#include "algorithms/MatrixAlg.h"
#include "algorithms/EigenSystemSolvers.h"

// Create system Ax = b
Matrix<Real> A(3, 3);
A[0][0] = 4; A[0][1] = 2; A[0][2] = 1;
A[1][0] = 2; A[1][1] = 5; A[1][2] = 2;
A[2][0] = 1; A[2][1] = 2; A[2][2] = 4;
Vector<Real> b{1, 2, 3};

// Gauss-Jordan solver (static method)
Vector<Real> x_gj = GaussJordanSolver<Real>::SolveConst(A, b);

// LU solver (general)
LUSolver<Real> lu(A);
Vector<Real> x_lu = lu.Solve(b);

// Cholesky (if SPD - use MatrixSym or verify positive definiteness)
CholeskySolver<Real> chol(A);
Vector<Real> x_chol = chol.Solve(b);

// QR (for overdetermined)
QRSolver<Real> qr(A);
Vector<Real> x_qr = qr.Solve(b);

// SVD (most robust - no template parameter)
SVDecompositionSolver svd(A);
Vector<Real> x_svd = svd.Solve(b);

// Eigenvalues (symmetric matrices only - uses MatrixSym)
MatrixSym<Real> S(3);
S(0, 0) = 4; S(0, 1) = 2; S(0, 2) = 1;
S(1, 1) = 5; S(1, 2) = 2;
S(2, 2) = 4;

auto result = SymmMatEigenSolverJacobi::Solve(S);
// result.eigenvalues  - Vector of eigenvalues (sorted)
// result.eigenvectors - Matrix of column eigenvectors
// result.iterations   - Number of Jacobi sweeps
// result.converged    - Whether iterations converged

// Verify orthogonality of eigenvectors
bool orth = MatrixAlg::IsOrthogonal(result.eigenvectors, 1e-12);
```

## Error Sources

1. **Roundoff accumulation** - Increases with matrix size
2. **Pivoting failures** - Mitigated by partial pivoting
3. **Ill-conditioning** - Fundamental limitation
4. **Catastrophic cancellation** - In determinant computation
5. **Loss of orthogonality** - In iterative methods

## Performance Tips

1. **Reuse factorizations** - LU, QR, Cholesky factor once, solve many times
2. **Check condition number** - Before trusting results
3. **Use appropriate solver** - Cholesky for SPD, QR for least squares
4. **Consider iterative methods** - For very large sparse systems
5. **Verify solutions** - Always compute residual ||Ax - b||

## Conclusions

1. **LU decomposition** is the workhorse for general systems
2. **Cholesky** is optimal for SPD matrices (faster, more stable)
3. **QR** is essential for least squares and maintains numerical stability
4. **SVD** is most robust but slowest; use for rank-deficient systems
5. **Jacobi eigenvalue** method is reliable for symmetric matrices
6. **Condition number** determines achievable accuracy
7. **All MML solvers** achieve near-optimal precision for their complexity class

---

*Test file: `src/testing_precision/test_precision_linalg.cpp`*
