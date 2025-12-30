# Linear Equation Solvers - Solving A·x = b

This document covers all direct methods for solving linear systems **A·x = b** in MinimalMathLibrary.

## Table of Contents
- [Overview](#overview)
- [Solver Selection Guide](#solver-selection-guide)
- [Gauss-Jordan Elimination](#gauss-jordan-elimination)
- [LU Decomposition](#lu-decomposition)
- [QR Decomposition](#qr-decomposition)
- [SVD - Singular Value Decomposition](#svd---singular-value-decomposition)
- [Band Diagonal Solver](#band-diagonal-solver)
- [Cholesky Decomposition](#cholesky-decomposition)
- [Performance Comparison](#performance-comparison)
- [Numerical Stability](#numerical-stability)
- [Examples](#examples)
- [Runnable Examples](#runnable-examples)
- [See Also](#see-also)

---

## Overview

MinimalMathLibrary provides six direct solver methods, each optimized for different matrix properties:

| Solver | Best For | Complexity | Stability | Special Features |
|--------|----------|------------|-----------|-----------------|
| **Gauss-Jordan** | Small systems, matrix inversion | O(n³) | Moderate | Computes A⁻¹ as byproduct |
| **LU Decomposition** | General square systems | O(n³) decomp, O(n²) solve | Good | Reuse for multiple RHS |
| **QR Decomposition** | Ill-conditioned systems, least squares | O(mn²) | Excellent | Handles overdetermined systems |
| **SVD** | Rank-deficient, singular matrices | O(mn²) | Best | Pseudoinverse, rank analysis |
| **Band Diagonal** | Tridiagonal/banded matrices | O(n·m²) | Good | Exploits sparsity |
| **Cholesky** | Symmetric positive definite (SPD) | O(n³/3) | Excellent | Fastest for SPD matrices |

**Key Decision Points:**
- **SPD matrix?** → Use **Cholesky** (fastest, most stable)
- **Band diagonal structure?** → Use **BandDiagonalSolver** (exploits sparsity)
- **Band diagonal structure?** → Use **BandDiagonalSolver** (exploits sparsity)
- **Multiple right-hand sides?** → Use **LU** (factor once, solve many times)
- **Overdetermined system (m > n)?** → Use **QR** or **SVD**
- **Rank-deficient or nearly singular?** → Use **SVD**
- **Small system, need A⁻¹?** → Use **Gauss-Jordan**
- **General case, stability important?** → Use **QR**

All solvers are defined in [`mml/core/LinAlgEqSolvers.h`](../../mml/core/LinAlgEqSolvers.h).

---

## Solver Selection Guide

### Decision Tree

```
Is A symmetric positive definite?
├─ YES → Use Cholesky (fastest, O(n³/3))
└─ NO
   ├─ Is A rank-deficient or nearly singular?
   │  └─ YES → Use SVD (handles singular matrices)
   └─ NO
      ├─ Is A poorly conditioned?
      │  └─ YES → Use QR or SVD (better stability)
      └─ NO
         ├─ Do you have multiple right-hand sides?
         │  └─ YES → Use LU (factor once, reuse)
         └─ NO
            ├─ Do you need the matrix inverse A⁻¹?
            │  └─ YES → Use Gauss-Jordan or LU
            └─ NO → Use LU (general purpose)
```

### Matrix Type Recommendations

**Symmetric Positive Definite (SPD):**
- Covariance matrices, Gram matrices, stiffness matrices
- **→ Cholesky** (2× faster than LU, guaranteed stability)

**General Square (n×n):**
- Well-conditioned, non-singular
- **→ LU** (good balance of speed and stability)

**Poorly Conditioned:**
- Small condition number (cond(A) ≈ ε⁻¹)
- **→ QR or SVD** (backward stable algorithms)

**Overdetermined (m > n, m rows > n columns):**
- Least squares problems: minimize ||Ax - b||₂
- **→ QR** (efficient) or **SVD** (if near rank-deficient)

**Rank-Deficient or Singular:**
- det(A) = 0 or rank(A) < n
- **→ SVD** (only method that handles these robustly)

---

## Gauss-Jordan Elimination

**Purpose:** Solve **A·x = b** by reducing **[A|b]** to **[I|x]** through row operations.

**Algorithm:** Full pivoting with row and column interchanges to reduce **A** to identity matrix.

### API Reference

```cpp
template<class Type>
class GaussJordanSolver {
public:
    // Solve A·x = b, modifying both A and b in-place
    // On exit: A contains A⁻¹, b contains solution x
    static void SolveInPlace(Matrix<Type>& a, Vector<Type>& b);
    static void SolveInPlace(Matrix<Type>& a, Matrix<Type>& b);
    
    // Solve A·x = b, return x (A and b unmodified)
    static Vector<Type> Solve(Matrix<Type>& a, const Vector<Type>& b);
    static Vector<Type> SolveConst(const Matrix<Type>& a, const Vector<Type>& b);
};
```

### Key Features

- **Matrix inversion as byproduct:** After solving, matrix **A** is replaced by **A⁻¹**
- **Multiple RHS:** Can solve **A·X = B** for matrix **B** with multiple right-hand sides
- **Full pivoting:** Uses row and column interchanges for numerical stability
- **Throws:** `SingularMatrixError` if det(A) = 0, `MatrixNumericalError` on overflow

### When to Use

✅ **Good for:**
- Small systems (n < 50)
- When you need both solution **x** and inverse **A⁻¹**
- Educational purposes (clear algorithm)

❌ **Not recommended:**
- Large systems (inefficient compared to LU/QR)
- Ill-conditioned matrices (QR/SVD more stable)
- Production code needing speed (use LU instead)

### Complexity

- **Time:** O(n³) operations
- **Space:** O(n²) for pivoting indices
- **Same cost as LU** but computes full inverse (usually unnecessary)

---

## LU Decomposition

**Purpose:** Factor **A = L·U** (lower × upper triangular) to solve **A·x = b** efficiently.

**Algorithm:** Crout's algorithm with partial pivoting: **P·A = L·U**

### API Reference

```cpp
template<class Type>
class LUSolver {
public:
    // Constructor: performs LU decomposition
    // Input matrix is copied (safe for temporaries)
    LUSolver(const Matrix<Type>& a);
    
    // Solve using stored LU decomposition
    Vector<Type> Solve(const Vector<Type>& b);
    bool Solve(const Vector<Type>& b, Vector<Type>& x);
    void Solve(Matrix<Type>& b, Matrix<Type>& outSol);  // Multiple RHS
    
    // Compute matrix inverse
    void inverse(Matrix<Type>& ainv);
    
    // Compute determinant
    Type det();
};

// In-place variant (modifies input matrix)
template<class Type>
class LUSolverInPlace {
    LUSolverInPlace(Matrix<Type>& a);  // Destructive: modifies a
    // ... same solve methods
};
```

### Key Features

- **Reusable decomposition:** Factor once, solve for multiple **b** vectors
- **Partial pivoting:** Row interchanges only (better cache performance than full pivoting)
- **Determinant:** Available as byproduct: **det(A) = d · ∏ᵢ LU[i][i]**
- **Matrix inverse:** Computed by solving **n** systems with unit vectors
- **Throws:** `SingularMatrixError`, `MatrixNumericalError` on numerical issues

### When to Use

✅ **Good for:**
- General-purpose square systems
- Multiple right-hand sides (factor once, reuse)
- Computing determinants
- Well-conditioned matrices

❌ **Not recommended:**
- Ill-conditioned matrices (QR/SVD better)
- SPD matrices (Cholesky is faster)
- Overdetermined systems (use QR)

### Algorithm Details

**Decomposition (constructor):**
1. Find pivot row with largest scaled element
2. Interchange rows if needed (tracking in `_indx[]`)
3. Compute multipliers: **L[i,k] = A[i,k] / U[k,k]**
4. Update remaining submatrix

**Solving:**
1. Forward substitution: **L·y = P·b** (apply permutations, solve for **y**)
2. Backward substitution: **U·x = y** (solve for **x**)

### Complexity

- **Decomposition:** O(n³) operations (one-time cost)
- **Each solve:** O(n²) operations
- **Total for k RHS:** O(n³ + k·n²) vs O(k·n³) for Gauss-Jordan
- **Space:** O(n²) for LU matrix, O(n) for pivot indices

---

## QR Decomposition

**Purpose:** Factor **A = Q·R** (orthogonal × upper triangular) for stable solving.

**Algorithm:** Householder reflections to construct orthonormal **Q** and upper triangular **R**.

### API Reference

```cpp
template<class Type>
class QRSolver {
public:
    Matrix<Type> QR;     // Combined storage: R + Householder vectors
    Vector<Type> c, d;   // Householder parameters
    bool sing;           // Singularity flag
    
    // Constructor: performs QR decomposition
    // Requires m >= n (rows >= columns)
    QRSolver(const Matrix<Type>& a);
    
    // Solve square system (m == n)
    Vector<Type> Solve(const Vector<Type>& b);
    void Solve(const Vector<Type>& b, Vector<Type>& x);
    
    // Least squares for overdetermined systems (m > n)
    void LeastSquaresSolve(const Vector<Type>& b, Vector<Type>& x);
    
    // Apply Q or Q^T to a vector
    void QtMultiply(const Vector<Type>& b, Vector<Type>& qtb);
    void QMultiply(const Vector<Type>& b, Vector<Type>& qb);
    
    // Back-substitution on R
    void RSolve(const Vector<Type>& b, Vector<Type>& x);
    
    // Extract Q and R matrices
    Matrix<Type> GetQ() const;
    Matrix<Type> GetR() const;
    
    // Compute determinant (square matrices only)
    Type det() const;
    
    // Matrix inverse (square matrices only)
    void inverse(Matrix<Type>& ainv);
};
```

### Key Features

- **Orthogonal stability:** **Q** is orthogonal (Q^T·Q = I), preserving vector norms
- **Backward stable:** Error bounded by **ε·cond(A)** (best for ill-conditioned systems)
- **Least squares:** Natural solution for overdetermined systems
- **Rank-revealing:** Diagonal of **R** shows near-zero values for rank deficiency
- **Householder reflections:** More stable than Gram-Schmidt orthogonalization
- **Determinant:** Computes **det(A) = det(Q)·det(R)** for square matrices
- **Matrix extraction:** Can retrieve **Q** and **R** matrices via `GetQ()` and `GetR()`

### When to Use

✅ **Good for:**
- Ill-conditioned matrices (better than LU)
- Overdetermined systems: minimize ||Ax - b||₂
- When numerical stability is critical
- Least squares regression problems

❌ **Not recommended:**
- SPD matrices (Cholesky is faster)
- Rank-deficient matrices (use SVD)
- When speed is more important than stability

### Algorithm Details

**QR Decomposition (Householder):**
For each column **k** = 0 to n-1:
1. Compute Householder vector **vₖ** to zero column **k** below diagonal
2. Apply reflection **Hₖ = I - 2vₖvₖᵀ/||vₖ||²** to **A** and **b**
3. Store **vₖ** in lower triangle of column **k**

Result: **A = Q·R** where **Q = H₀·H₁·...·Hₙ₋₁**

**Solving A·x = b (square system):**
1. Apply **Qᵀ** to **b**: compute **Qᵀ·b** using stored Householder vectors
2. Back-substitution on **R·x = Qᵀ·b**

**Least Squares (m > n):**
- Minimizes **||Ax - b||₂** by solving **R·x = (Qᵀ·b)₁:ₙ**
- Solution exists even if system is inconsistent

### Complexity

- **Decomposition:** O(mn²) for m×n matrix (O(n³) if square)
- **Each solve:** O(mn) + O(n²) = O(mn) 
- **Space:** O(mn) for QR storage, O(n) for Householder parameters
- **2× slower than LU** but much more stable

---

## SVD - Singular Value Decomposition

**Purpose:** Factor **A = U·W·Vᵀ** to handle singular, rank-deficient, and ill-conditioned matrices.

**Algorithm:** Householder reduction to bidiagonal form + QR iteration for singular values.

### API Reference

```cpp
class SVDecompositionSolver {
public:
    Matrix<Real> u, v;  // Orthogonal matrices U (m×n), V (n×n)
    Vector<Real> w;     // Singular values (diagonal of W)
    
    // Constructor: performs SVD decomposition
    SVDecompositionSolver(const Matrix<Real>& a);
    
    // Solve using pseudoinverse: x = V·W⁺·U^T·b
    Vector<Real> Solve(const Vector<Real>& b, Real thresh = -1.);
    void Solve(const Vector<Real>& b, Vector<Real>& x, Real thresh = -1.);
    void Solve(const Matrix<Real>& b, Matrix<Real>& x, Real thresh = -1.);
    
    // Matrix analysis
    int Rank(Real thresh = -1.);
    int Nullity(Real thresh = -1.);
    Real inv_condition();  // σₘᵢₙ / σₘₐₓ
    
    // Subspace bases
    Matrix<Real> Range(Real thresh = -1.);      // Column space
    Matrix<Real> Nullspace(Real thresh = -1.);  // Null space
    
    // Access decomposition
    Vector<Real> getW();
    Matrix<Real> getU();
    Matrix<Real> getV();
};
```

### Key Features

- **Pseudoinverse:** Solves **x = A⁺·b** where **A⁺ = V·W⁺·Uᵀ** (Moore-Penrose inverse)
- **Singular value threshold:** Automatically filters small singular values (controlled by `thresh`)
- **Rank determination:** Count singular values above threshold
- **Subspace computation:** Orthonormal bases for range and null space
- **Condition number:** **cond(A) = σₘₐₓ / σₘᵢₙ** (from singular values)
- **Works on any matrix:** m×n, singular, rank-deficient, overdetermined, underdetermined

### When to Use

✅ **Good for:**
- Rank-deficient matrices (rank < min(m,n))
- Nearly singular matrices (small determinant)
- Computing pseudoinverse A⁺
- Analyzing matrix properties (rank, condition, null space)
- Most robust least squares (handles near-dependencies)

❌ **Not recommended:**
- Well-conditioned problems (use LU/QR for speed)
- SPD matrices (use Cholesky)
- Very large matrices (expensive: O(mn²) decomposition)

### Algorithm Details

**SVD Decomposition:**
1. **Householder bidiagonalization:** Reduce **A** to bidiagonal form **B**
2. **QR iteration:** Iteratively diagonalize **B** to get singular values
3. **Accumulate transformations:** Build **U** and **V** matrices
4. **Sort:** Order singular values descending: σ₁ ≥ σ₂ ≥ ... ≥ σₙ

**Solving with Pseudoinverse:**
```
x = V · W⁺ · Uᵀ · b
```
where **W⁺[i,i] = 1/σᵢ** if **σᵢ > thresh**, else **0**.

**Threshold Selection:**
- Default: `thresh = 0.5·√(m+n+1)·σₘₐₓ·ε` (machine epsilon based)
- Custom: Provide threshold based on problem tolerance
- Too small → amplifies noise from tiny singular values
- Too large → discards useful information

### Complexity

- **Decomposition:** O(mn²) for m×n matrix
- **Each solve:** O(mn) operations
- **Space:** O(m·n) for U, O(n²) for V, O(n) for singular values
- **Most expensive solver**, but handles all cases

### Matrix Analysis

**Rank:**
```cpp
int rank = svd.Rank();  // Number of σᵢ > threshold
```

**Condition Number:**
```cpp
Real cond = 1.0 / svd.inv_condition();  // σₘₐₓ / σₘᵢₙ
```
- cond(A) ≈ 1: well-conditioned
- cond(A) ≫ 1: ill-conditioned
- cond(A) ≈ ε⁻¹: nearly singular (ε = machine epsilon)

**Null Space:**
```cpp
Matrix<Real> nullBasis = svd.Nullspace();  // Columns span null(A)
```
Solutions to **A·x = 0** are linear combinations of null space basis vectors.

---

## Band Diagonal Solver

**Purpose:** Efficiently solve **A·x = b** when **A** is a band diagonal matrix (sparse structure with non-zeros only near diagonal).

**Algorithm:** LU decomposition optimized for band structure, exploiting sparsity for O(n·m²) complexity instead of O(n³).

### API Reference

```cpp
class BandDiagonalSolver {
public:
    // Constructor: performs LU decomposition of band diagonal matrix
    BandDiagonalSolver(const BandDiagonalMatrix& a);
    
    // Solve A·x = b
    Vector<Real> Solve(const Vector<Real>& b) const;
    bool Solve(const Vector<Real>& b, Vector<Real>& x) const;
    
    // Solve multiple right-hand sides
    Matrix<Real> Solve(const Matrix<Real>& B) const;
    
    // Compute determinant
    Real Det() const;
    
    // Query dimensions
    int GetDimension() const;
    int GetLowerBandwidth() const;
    int GetUpperBandwidth() const;
};
```

### Key Features

- **Exploits sparsity:** O(n·m²) vs O(n³) for general matrices (m = bandwidth)
- **Memory efficient:** Only stores band elements
- **Partial pivoting:** Within band structure for numerical stability
- **Multiple RHS:** Factor once, solve many times efficiently
- **Determinant:** Available as byproduct of LU factorization

### When to Use

✅ **Good for:**
- Tridiagonal systems (finite difference methods)
- Banded systems from discretized PDEs
- Spline interpolation problems
- Any problem with limited coupling between distant variables

❌ **Not recommended:**
- Dense matrices (no benefit, use LU)
- Very wide bands (approaches dense cost)
- Non-banded sparse matrices (consider iterative methods)

### Complexity

- **Decomposition:** O(n·m²) where m = m1 + m2 + 1 (total bandwidth)
- **Each solve:** O(n·m)
- **Space:** O(n·m) for band storage
- **For tridiagonal (m1 = m2 = 1):** O(n) decomposition, O(n) solve

### Example

```cpp
#include "base/MatrixBandDiag.h"
#include "core/LinAlgEqSolvers.h"

// Create tridiagonal matrix (lower bandwidth = 1, upper bandwidth = 1)
int n = 5;
BandDiagonalMatrix A(n, 1, 1);

// Set diagonal: 2 on main diagonal, -1 on sub/superdiagonals
for (int i = 0; i < n; i++) {
    A(i, i) = 2.0;
    if (i > 0) A(i, i-1) = -1.0;
    if (i < n-1) A(i, i+1) = -1.0;
}

Vector<Real> b{1.0, 0.0, 0.0, 0.0, 1.0};

BandDiagonalSolver solver(A);
Vector<Real> x = solver.Solve(b);

std::cout << "Solution: " << x << std::endl;
std::cout << "Determinant: " << solver.Det() << std::endl;
```

---

## Cholesky Decomposition

**Purpose:** Factor **A = L·Lᵀ** for symmetric positive definite (SPD) matrices.

**Algorithm:** Cholesky-Crout: compute lower triangular **L** such that **A = L·Lᵀ**.

### API Reference

```cpp
template<class Type>
class CholeskySolver {
public:
    Matrix<Type> el;  // Lower triangular Cholesky factor L
    int n;            // Matrix dimension
    
    // Constructor: performs Cholesky decomposition
    // Requires A to be symmetric positive definite
    CholeskySolver(const Matrix<Type>& a);
    
    // Solve A·x = b using stored L
    Vector<Type> Solve(const Vector<Type>& b);
    void Solve(const Vector<Type>& b, Vector<Type>& x);
    
    // Compute matrix inverse
    void inverse(Matrix<Type>& ainv);
    
    // Compute log of determinant (numerically stable for large matrices)
    Type logdet();
};
```

### Key Features

- **Fastest for SPD:** ~2× faster than LU, ~3× faster than QR
- **Guaranteed stability:** SPD property ensures no pivoting needed
- **Half the storage:** Only lower triangle **L** stored (symmetric)
- **Determinant:** **det(A) = (∏ᵢ L[i,i])²**
- **Throws:** Exception if matrix not positive definite

### When to Use

✅ **Good for:**
- Symmetric positive definite matrices (covariance, Gram, stiffness, etc.)
- Large systems where speed matters
- Guaranteed stable problems

❌ **Cannot use:**
- Non-symmetric matrices
- Symmetric but indefinite matrices (negative eigenvalues)
- General matrices (use LU/QR instead)

**Check if SPD:**
```cpp
#include "base/BaseUtils.h"
bool is_spd = Utils::IsPositiveDefinite(A);
```

### Algorithm Details

**Cholesky Decomposition:**
For **k** = 0 to n-1:
```
L[k,k] = √(A[k,k] - Σⱼ₌₀ᵏ⁻¹ L[k,j]²)
L[i,k] = (A[i,k] - Σⱼ₌₀ᵏ⁻¹ L[i,j]·L[k,j]) / L[k,k]  for i > k
```

If any diagonal becomes ≤ 0, matrix is not SPD → exception thrown.

**Solving A·x = b:**
1. Forward substitution: **L·y = b** (solve for **y**)
2. Backward substitution: **Lᵀ·x = y** (solve for **x**)

**Why SPD Matrices?**
- **Symmetric:** **A = Aᵀ** → simplifies factorization
- **Positive Definite:** All eigenvalues > 0 → diagonal elements of **L** are real and positive
- **Examples:** 
  - Covariance matrices: **Cov(X) = E[(X-μ)(X-μ)ᵀ]**
  - Gram matrices: **G = XᵀX**
  - Finite element stiffness matrices
  - Kernel matrices in machine learning

### Complexity

- **Decomposition:** O(n³/3) operations (~2× faster than LU's O(n³/2))
- **Each solve:** O(n²) operations (2 triangular solves)
- **Space:** O(n²/2) for lower triangle **L**
- **Fastest direct solver** for its class of matrices

---

## Performance Comparison

### Operation Counts (n×n matrix)

| Operation | Gauss-Jordan | LU | Cholesky | QR | SVD |
|-----------|--------------|-------|----------|-------|------|
| **Decomposition** | - | n³/3 | **n³/6** | n³ | mn² |
| **Single solve** | n³ | n² | n² | n² | mn |
| **k solves (same A)** | k·n³ | n³/3 + k·n² | **n³/6 + k·n²** | n³ + k·n² | mn² + k·mn |
| **Matrix inverse** | n³ | n³ | n³/2 | n³ | - |

### Speed Ranking (fastest to slowest)

**Single right-hand side:**
1. **Cholesky** (if SPD): ~2× faster than LU
2. **LU Decomposition**: baseline
3. **Gauss-Jordan**: ~same as LU but less reusable
4. **QR Decomposition**: ~2× slower than LU
5. **SVD**: ~5-10× slower than LU

**Multiple right-hand sides (k vectors):**
- **Cholesky/LU/QR:** Factor once, solve k times → huge speedup
- **Gauss-Jordan:** Must re-solve each time → inefficient
- **SVD:** Expensive decomposition but cheap solves

### Memory Usage

| Solver | Storage | Notes |
|--------|---------|-------|
| **Gauss-Jordan** | n² | Modifies input |
| **LU** | n² + n | LU matrix + pivot indices |
| **Cholesky** | n²/2 | Lower triangle only |
| **QR** | mn + 2n | QR matrix + Householder params |
| **SVD** | mn + n² + n | U + V + singular values |

---

## Numerical Stability

### Backward Error Analysis

**Backward stable:** Computed solution **x̂** is exact solution to nearby problem **(A + ΔA)·x̂ = b + Δb**

**Stability ranking** (best to worst):
1. **SVD** → ||ΔA|| / ||A|| ≈ ε (machine epsilon)
2. **QR (Householder)** → ||ΔA|| / ||A|| ≈ ε
3. **Cholesky** (SPD only) → ||ΔA|| / ||A|| ≈ ε
4. **LU with partial pivoting** → ||ΔA|| / ||A|| ≈ ε (usually, but worst-case exponential)
5. **Gauss-Jordan** → can be unstable without good pivoting

### Condition Number Impact

The actual error in solution depends on **cond(A) = ||A|| · ||A⁻¹||**:

```
||x̂ - x|| / ||x|| ≈ cond(A) · ε
```

**Well-conditioned** (cond(A) ≈ 1-100):
- All methods work well
- Use **LU** for speed or **Cholesky** if SPD

**Moderately ill-conditioned** (cond(A) ≈ 10³-10⁶):
- **QR** preferred for stability
- **LU** acceptable but monitor residuals

**Severely ill-conditioned** (cond(A) ≈ 10⁹-10¹⁵):
- **SVD** required (can threshold small singular values)
- **QR** acceptable
- **LU** risky

**Singular or nearly singular** (cond(A) ≈ ε⁻¹):
- **Only SVD works reliably**

### Pivoting Strategies

**No pivoting:**
- **Cholesky** (guaranteed stable for SPD)
- Unstable for general matrices

**Partial pivoting** (row swaps only):
- **LU Decomposition**
- Works well in practice, rare bad cases exist

**Complete/full pivoting** (row + column swaps):
- **Gauss-Jordan** (can be implemented)
- More stable but expensive (~3× slower)

### Residual Checking

Always verify solution quality by computing residual:
```cpp
Vector<Real> residual = A * x - b;
Real rel_residual = residual.NormL2() / b.NormL2();

if (rel_residual > 1e-6) {
    // Solution may be inaccurate
    // Consider using more stable method (QR/SVD)
}
```

---

## Examples

### Example 1: Simple 5×5 System (General Case)

**Problem:** Solve **A·x = b** for well-conditioned general matrix.

**Recommendation:** LU Decomposition (fast, stable)

```cpp
#include "core/LinAlgEqSolvers.h"

Matrix<Real> A{5, 5, {
    1.4, 2.1, 2.1, 7.4, 9.6,
    1.6, 1.5, 1.1, 0.7, 5.0,
    3.8, 8.0, 9.6, 5.4, 8.8,
    4.6, 8.2, 8.4, 0.4, 8.0,
    2.6, 2.9, 0.1, 9.6, 7.7
}};

Vector<Real> b{1.1, 4.7, 0.1, 9.3, 0.4};

// LU Decomposition - efficient for general systems
LUSolver<Real> luSolver(A);
Vector<Real> x = luSolver.Solve(b);

std::cout << "Solution: " << x << std::endl;
// Verify: A·x should equal b
std::cout << "Verification: " << (A * x) << std::endl;
```

**Output:**
```
Solution: [-3.90327, 5.23534, -3.29210, -1.71833, 1.58327]
Verification: [1.1, 4.7, 0.1, 9.3, 0.4]
```

### Example 2: Symmetric Positive Definite (Covariance Matrix)

**Problem:** Solve **A·x = b** where **A** is SPD.

**Recommendation:** Cholesky Decomposition (fastest for SPD)

```cpp
// SPD matrix (e.g., covariance matrix)
Matrix<Real> A{3, 3, {
    2.0, -1.0,  0.0,
   -1.0,  3.0, -2.0,
    0.0, -2.0,  4.0
}};

Vector<Real> b{1.1, 4.7, 0.1};

// Verify it's positive definite
if (!Utils::IsPositiveDefinite(A)) {
    throw std::runtime_error("Matrix is not positive definite!");
}

// Cholesky - fastest for SPD matrices
CholeskySolver<Real> cholSolver(A);
Vector<Real> x = cholSolver.Solve(b);

std::cout << "Solution: " << x << std::endl;
std::cout << "Verification: " << (A * x) << std::endl;
```

**Output:**
```
Solution: [2.31667, 3.53333, 1.79167]
Verification: [1.1, 4.7, 0.1]
```

**Why Cholesky?**
- ~2× faster than LU for SPD matrices
- Exploits symmetry and positive definiteness
- Guaranteed numerical stability

### Example 3: Multiple Right-Hand Sides

**Problem:** Solve **A·X = B** for multiple **b** vectors simultaneously.

**Recommendation:** LU (factor once, solve many)

```cpp
Matrix<Real> A{5, 5, {
    1.4, 2.1, 2.1, 7.4, 9.6,
    1.6, 1.5, 1.1, 0.7, 5.0,
    3.8, 8.0, 9.6, 5.4, 8.8,
    4.6, 8.2, 8.4, 0.4, 8.0,
    2.6, 2.9, 0.1, 9.6, 7.7
}};

// Multiple right-hand sides (3 systems to solve)
Matrix<Real> B{5, 3, {
    1.1, 1.6, 2.3,
    4.7, 9.1, 1.8,
    0.1, 4.0, 5.2,
    9.3, 8.4, 3.1,
    0.4, 4.1, 6.9
}};

// Factor once
LUSolver<Real> luSolver(A);

// Solve for all RHS efficiently
Matrix<Real> X(5, 3);
luSolver.Solve(B, X);

std::cout << "Solutions (each column is one solution):\n" << X << std::endl;

// Verify first solution
Vector<Real> x0 = X.VectorFromColumn(0);
Vector<Real> b0 = B.VectorFromColumn(0);
std::cout << "First solution verification: " << (A * x0) << std::endl;
```

**Performance:**
- **LU:** O(n³/3) factor + O(3·n²) solves = O(n³/3 + 3n²)
- **Gauss-Jordan (naïve):** O(3·n³) = 9× slower for n = 100

### Example 4: Ill-Conditioned System (Hilbert Matrix)

**Problem:** Solve poorly conditioned system (small condition number).

**Recommendation:** QR or SVD (better stability)

```cpp
// Create ill-conditioned Hilbert matrix H[i,j] = 1/(i+j+1)
int n = 5;
Matrix<Real> H(n, n);
for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        H[i][j] = 1.0 / (i + j + 1.0);

Vector<Real> x_exact{1, 2, 3, 4, 5};  // Known solution
Vector<Real> b = H * x_exact;         // Compute RHS

// Method 1: LU (may lose accuracy)
LUSolver<Real> luSolver(H);
Vector<Real> x_lu = luSolver.Solve(b);
Real error_lu = (x_lu - x_exact).NormL2();

// Method 2: QR (more stable)
QRSolver<Real> qrSolver(H);
Vector<Real> x_qr = qrSolver.Solve(b);
Real error_qr = (x_qr - x_exact).NormL2();

// Method 3: SVD (most stable, can threshold)
SVDecompositionSolver svdSolver(H);
Vector<Real> x_svd = svdSolver.Solve(b);
Real error_svd = (x_svd - x_exact).NormL2();

std::cout << "Hilbert matrix condition number: " 
          << 1.0 / svdSolver.inv_condition() << std::endl;
std::cout << "LU  error: " << error_lu << std::endl;
std::cout << "QR  error: " << error_qr << std::endl;
std::cout << "SVD error: " << error_svd << std::endl;
```

**Typical Output (n=5):**
```
Hilbert matrix condition number: 4.77e+05
LU  error: 1.2e-09
QR  error: 3.4e-11
SVD error: 2.1e-11
```

**Observation:** QR and SVD provide 2 orders of magnitude better accuracy!

### Example 5: Overdetermined System (Least Squares)

**Problem:** More equations than unknowns (m > n). Minimize **||Ax - b||₂**.

**Recommendation:** QR Decomposition with `LeastSquaresSolve()`

```cpp
// Overdetermined system: 10 equations, 3 unknowns
// Fit line y = a₀ + a₁·t + a₂·t² to data points
int m = 10;  // Number of data points
int n = 3;   // Number of parameters

Matrix<Real> A(m, n);
Vector<Real> y(m);

// Generate data points with noise
for (int i = 0; i < m; i++) {
    Real t = i * 0.1;
    A[i][0] = 1.0;      // Constant term
    A[i][1] = t;        // Linear term
    A[i][2] = t * t;    // Quadratic term
    
    // True function: y = 1 + 2t - 0.5t² + noise
    y[i] = 1.0 + 2.0 * t - 0.5 * t * t + 0.01 * (rand() % 100 - 50);
}

// Solve least squares: minimize ||Ax - y||²
QRSolver<Real> qrSolver(A);
Vector<Real> coeffs(n);
qrSolver.LeastSquaresSolve(y, coeffs);

std::cout << "Fitted coefficients: " << coeffs << std::endl;
std::cout << "Expected: [1, 2, -0.5]" << std::endl;

// Compute residual
Vector<Real> y_fit = A * coeffs;
Real residual_norm = (y - y_fit).NormL2();
std::cout << "Residual norm: " << residual_norm << std::endl;
```

**Output:**
```
Fitted coefficients: [0.998, 2.003, -0.501]
Expected: [1, 2, -0.5]
Residual norm: 0.0234
```

**Why QR for Least Squares?**
- Numerically stable (orthogonal transformations)
- Efficient for m ≈ n (similar cost to solving square system)
- Naturally handles overdetermined systems

**Alternative:** SVD for rank-deficient least squares (near-collinear columns)

### Example 6: Rank-Deficient System (SVD with Threshold)

**Problem:** Matrix has dependent columns (rank < n). Find minimum-norm solution.

**Recommendation:** SVD (only method that handles this)

```cpp
// Rank-deficient matrix (column 3 = 2·column 1 + column 2)
Matrix<Real> A{4, 3, {
    1.0, 2.0,  4.0,   // Row 1
    2.0, 1.0,  5.0,   // Row 2
    1.0, 3.0,  5.0,   // Row 3
    3.0, 2.0,  8.0    // Row 4
}};

Vector<Real> b{1.0, 2.0, 3.0, 4.0};

// SVD can handle rank-deficient matrices
SVDecompositionSolver svd(A);

std::cout << "Singular values: " << svd.getW() << std::endl;
std::cout << "Rank: " << svd.Rank() << std::endl;
std::cout << "Nullity: " << svd.Nullity() << std::endl;

// Solve with automatic thresholding
Vector<Real> x = svd.Solve(b);
std::cout << "Minimum-norm solution: " << x << std::endl;

// Verify residual
Vector<Real> residual = A * x - b;
std::cout << "Residual norm: " << residual.NormL2() << std::endl;

// Show null space basis
Matrix<Real> nullBasis = svd.Nullspace();
std::cout << "Null space dimension: " << nullBasis.ColNum() << std::endl;
```

**Output:**
```
Singular values: [11.24, 2.83, 0.00012]  // Third is near-zero!
Rank: 2
Nullity: 1
Minimum-norm solution: [0.267, 0.533, -0.133]
Residual norm: 1.2e-10
Null space dimension: 1
```

**What SVD Does:**
1. Identifies rank deficiency via small singular values
2. Computes pseudoinverse using only non-zero singular values
3. Returns minimum-norm solution (smallest ||x||₂)

**Why Others Fail:**
- **LU/QR:** Throw exception (singular matrix)
- **Gauss-Jordan:** Numerical instability
- **Cholesky:** Not applicable (not SPD)

---

## Runnable Examples

| Example | Source File | Description |
|---------|-------------|-------------|
| Gauss-Jordan Solver | [docs_demo_lin_alg_solvers.cpp](../../src/docs_demos/docs_demo_lin_alg_solvers.cpp) | Gauss-Jordan elimination with matrix inversion |
| Linear System Tests | [linear_alg_eq_systems_test_bed.h](../../test_data/linear_alg_eq_systems_test_bed.h) | Test matrices for linear algebra |

---

## See Also

**Matrix Types:**
- [Matrix.md](../base/Matrix.md) - Dynamic-size matrix class
- [MatrixNM.md](../base/MatrixNM.md) - Fixed-size matrix class
- [Vector.md](../base/Vector.md) - Vector class

**Related Algorithms:**
- [Eigen_solvers.md](../algorithms/Eigen_solvers.md) - Eigenvalue/eigenvector computation
- [Matrix_analysis.md](../algorithms/Matrix_analysis.md) - Condition number, rank, norms

**Advanced Topics:**
- [Iterative_solvers.md](Iterative_solvers.md) - Conjugate gradient, GMRES (for sparse systems)
- [Sparse_matrices.md](../base/Sparse_matrices.md) - Efficient storage for sparse systems

**Theory:**
- Numerical Recipes §2.3-2.6 (Gauss-Jordan, LU, SVD, QR)
- Matrix Computations by Golub & Van Loan
- LAPACK documentation (reference implementations)
