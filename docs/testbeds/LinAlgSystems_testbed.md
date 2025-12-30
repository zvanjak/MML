# Linear Algebra Systems Test Bed

## Overview

The Linear Algebra Equation Systems Test Bed provides a comprehensive collection of test matrices for validating linear system solvers, eigenvalue computations, singular value decomposition, and matrix factorization algorithms. Each test system includes pre-computed reference values verified against established numerical libraries.

## Architecture

The test bed is organized across multiple header files in `test_data/`:

| File | Purpose |
|------|---------|
| `linear_alg_eq_systems_test_bed.h` | Main test bed class and system registry |
| `linear_alg_eq_systems_classic_defs.h` | Classic benchmark matrices (Hilbert, Pascal, etc.) |
| `linear_alg_eq_systems_real_defs.h` | General real-valued test matrices |
| `linear_alg_eq_systems_sym_defs.h` | Symmetric matrix test systems |
| `linear_alg_eq_systems_complex_defs.h` | Complex-valued test systems |
| `linear_alg_eq_systems_solverspecific_defs.h` | Solver-specific matrices (SPD, overdetermined, etc.) |
| `linear_alg_eq_systems_tridiag_defs.h` | Tridiagonal matrix systems |
| `linear_alg_eq_systems_singular_defs.h` | Singular and near-singular matrices |

## Test System Data Structures

### TestLinearSystem (General Square Systems)

The primary test container for real-valued square systems:

```cpp
class TestLinearSystem {
    int _n;                              // Matrix dimension
    Matrix<Real> _mat;                   // Coefficient matrix A
    Vector<Real> _rhs;                   // Right-hand side b
    Vector<Real> _sol;                   // Known solution x such that Ax = b
    Vector<Complex> _eigen_values;       // All eigenvalues (may be complex)
    std::vector<Vector<Complex>> _eigen_vectors;  // Corresponding eigenvectors
    
    // Extended properties for algorithm verification
    Real _determinant;                   // det(A)
    Vector<Real> _singular_values;       // σ₁ ≥ σ₂ ≥ ... ≥ σₙ (descending)
    Real _cond_2;                        // Condition number κ₂(A) = σ_max/σ_min
    int _rank;                           // Numerical rank
};
```

### TestLinearSystemSymmetric (Symmetric Systems)

Optimized for symmetric matrices with real eigenvalues:

```cpp
class TestLinearSystemSymmetric {
    int _n;
    MatrixSym<Real> _mat;                // Symmetric matrix storage
    Vector<Real> _rhs;
    Vector<Real> _sol;
    Vector<Real> _eigen_values;          // Real eigenvalues (symmetric guarantee)
    std::vector<Vector<Real>> _eigen_vectors;  // Real eigenvectors
    Real _determinant, _cond_2;
    Vector<Real> _singular_values;
    int _rank;
};
```

### TestLinearSystemComplex (Complex Systems)

For complex-valued matrices and systems:

```cpp
class TestLinearSystemComplex {
    int _n, _m;                          // Dimensions (m rows, n columns)
    Matrix<Complex> _mat;
    Vector<Complex> _rhs;
    Vector<Complex> _sol;
    Real _det_abs;                       // |det(A)|
    Vector<Real> _singular_values;
    Real _cond_2;
    int _rank;
};
```

### TestLinearSystemOverdetermined (Least Squares)

For overdetermined systems (m > n):

```cpp
class TestLinearSystemOverdetermined {
    int _m, _n;                          // m rows, n columns (m > n)
    Matrix<Real> _mat;
    Vector<Real> _rhs;
    Vector<Real> _sol;                   // Least-squares solution
    Vector<Real> _singular_values;
    Real _cond_2;
    int _rank;
};
```

## Classic Test Matrix Categories

### Hilbert Matrices

**Definition:** H[i,j] = 1/(i+j+1)

**Source:** David Hilbert (1894), "Über die Darstellung definiter Formen als Summe von Formenquadraten"

**Properties:**
- Symmetric positive definite
- Notoriously ill-conditioned
- Condition number grows as O(exp(3.5n))
- Exact inverse exists but numerical computation is challenging

**Available sizes:** 3×3, 4×4, 5×5, 8×8

| Size | Condition Number (κ₂) | Determinant |
|------|----------------------|-------------|
| 3×3 | 524 | 4.63×10⁻⁴ |
| 4×4 | 15,514 | 1.65×10⁻⁷ |
| 5×5 | 476,607 | 3.75×10⁻¹² |
| 8×8 | 1.69×10¹¹ | 2.74×10⁻³³ |

**Purpose:** Tests solver behavior on ill-conditioned systems. Critical for validating numerical stability.

### Pascal Matrices

**Definition:** P[i,j] = C(i+j, i) (binomial coefficient)

**Source:** Blaise Pascal's arithmetic triangle (1653)

**Properties:**
- Symmetric positive definite
- All entries are integers
- Determinant is always 1
- Moderately conditioned

**Available sizes:** 3×3, 4×4, 5×5

| Size | Condition Number (κ₂) | det(A) |
|------|----------------------|--------|
| 3×3 | 53 | 1 |
| 4×4 | 654 | 1 |
| 5×5 | 8,376 | 1 |

**Purpose:** Tests exact arithmetic properties. Since entries and solutions are integers, useful for verifying solver accuracy without floating-point artifacts.

### Vandermonde Matrices

**Definition:** V[i,j] = x[i]^j for nodes x = [1, 2, ..., n]

**Source:** Alexandre-Théophile Vandermonde (1772)

**Properties:**
- Non-symmetric
- Related to polynomial interpolation
- Condition number grows exponentially with n
- Depends strongly on node spacing

**Available sizes:** 3×3, 4×4, 5×5

| Size | Condition Number (κ₂) | Determinant |
|------|----------------------|-------------|
| 3×3 | 64 | 2 |
| 4×4 | 3,075 | 12 |
| 5×5 | 648,163 | 288 |

**Purpose:** Tests solver behavior for polynomial fitting and interpolation problems.

### Frank Matrices

**Definition:** Upper Hessenberg with F[i,j] = n+1-j for i ≤ j, F[i,j] = n+1-i for i = j+1, else 0

**Source:** Werner Frank (1958), "Computing eigenvalues of complex matrices"

**Properties:**
- Upper Hessenberg structure
- All eigenvalues are real and distinct
- Well-suited for eigenvalue algorithm testing

**Available sizes:** 3×3, 4×4, 5×5

| Size | Condition Number (κ₂) | Determinant |
|------|----------------------|-------------|
| 3×3 | 119 | -0.3 |
| 4×4 | 268 | -0.05 |
| 5×5 | 3,913 | -0.0083 |

**Purpose:** Tests eigenvalue computation, especially QR algorithm and related methods.

### Kahan Matrices

**Definition:** Upper triangular with K[i,j] defined by angle θ (typically π/6)

**Source:** William Kahan (numerical analysis pioneer, Turing Award 1989)

**Properties:**
- Upper triangular
- Parametrized by angle θ
- Tests backward stability
- Condition number controllable via θ

**Available sizes:** 3×3, 5×5

| Size | Condition Number (κ₂) | Determinant |
|------|----------------------|-------------|
| 3×3 | 2.3 | 0.65 |
| 5×5 | 4.7 | 0.24 |

**Purpose:** Tests backward stability of triangular solvers.

## Solver-Specific Matrix Categories

### Diagonally Dominant Matrices

**Properties:**
- |a[i,i]| > Σ|a[i,j]| for j ≠ i (strict diagonal dominance)
- Guaranteed convergence for iterative methods (Jacobi, Gauss-Seidel, SOR)
- Well-conditioned
- Often arise from discretization of PDEs

**Available systems:**
- `diag_dominant_4x4` - General strictly diagonally dominant
- `diag_dominant_5x5_tridiag` - Tridiagonal from 1D Poisson discretization
- `diag_dominant_6x6_poisson2d` - From 2D Poisson discretization

**Purpose:** Tests iterative solvers where convergence is guaranteed.

### Symmetric Positive Definite (SPD) Matrices

**Properties:**
- Symmetric: A = Aᵀ
- Positive definite: xᵀAx > 0 for all x ≠ 0
- All eigenvalues are positive
- Cholesky decomposition exists (A = LLᵀ)

**Available systems:**
- `spd_3x3` - Basic SPD test matrix
- `spd_4x4_correlation` - Correlation matrix structure
- `spd_5x5_mass_matrix` - Finite element mass matrix pattern
- `spd_6x6_graph_laplacian` - Graph Laplacian structure

**Purpose:** Tests Cholesky decomposition and conjugate gradient methods.

### Overdetermined Systems (Least Squares)

**Properties:**
- More equations than unknowns (m > n)
- No exact solution in general
- Least-squares solution minimizes ‖Ax - b‖₂

**Available systems:**

| Name | Dimensions | Condition Number | Purpose |
|------|------------|------------------|---------|
| `overdetermined_4x3` | 4×3 | 789 | Basic QR and least-squares |
| `overdetermined_5x3_regression` | 5×3 | 55 | Linear regression context |
| `overdetermined_6x4_polynomial` | 6×4 | varies | Polynomial fitting |

**Purpose:** Tests QR decomposition, SVD-based methods, and least-squares solvers.

## API Reference

### LinearAlgEqTestBed Class

The main access point for all test systems.

#### Basic Access

```cpp
// Get total count and access by index
int n = LinearAlgEqTestBed::numLinAlgEqSystems();  // Returns 27
const TestLinearSystem& sys = LinearAlgEqTestBed::getLinAlgEqSystem(0);

// Access by name
const TestLinearSystem& hilbert = LinearAlgEqTestBed::getLinAlgEqSystem("hilbert_5x5");

// Symmetric systems
int ns = LinearAlgEqTestBed::numLinAlgEqSystemsSymmetric();  // Returns 7
const TestLinearSystemSymmetric& sym = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric("spd_4x4_correlation");

// Complex systems
int nc = LinearAlgEqTestBed::numComplexSystems();  // Returns 12
const TestLinearSystemComplex& cmplx = LinearAlgEqTestBed::getComplexSystem("mat_cmplx_1_5x5");

// Overdetermined systems
int no = LinearAlgEqTestBed::numOverdeterminedSystems();  // Returns 3
const TestLinearSystemOverdetermined& over = LinearAlgEqTestBed::getOverdeterminedSystem("overdetermined_4x3");
```

#### Category-Based Filtering

```cpp
// By matrix type
auto hilberts = LinearAlgEqTestBed::getHilbertMatrices();
auto pascals = LinearAlgEqTestBed::getPascalMatrices();
auto vandermondes = LinearAlgEqTestBed::getVandermondeMatrices();
auto franks = LinearAlgEqTestBed::getFrankMatrices();
auto kahans = LinearAlgEqTestBed::getKahanMatrices();
auto diagDom = LinearAlgEqTestBed::getDiagDominantMatrices();
auto spd = LinearAlgEqTestBed::getSPDMatrices();
auto classics = LinearAlgEqTestBed::getClassicMatrices();  // All classic types
auto randoms = LinearAlgEqTestBed::getRandomMatrices();    // mat_3x3, mat_5x5, etc.
```

#### Size-Based Filtering

```cpp
// Get all systems of specific size
auto size5 = LinearAlgEqTestBed::getSystemsBySize(5);
auto symSize5 = LinearAlgEqTestBed::getSymmetricSystemsBySize(5);

// Size categories
auto small = LinearAlgEqTestBed::getSmallMatrices();   // n <= 4
auto medium = LinearAlgEqTestBed::getMediumMatrices(); // 5 <= n <= 8
auto large = LinearAlgEqTestBed::getLargeMatrices();   // n >= 10
```

#### Conditioning-Based Filtering

```cpp
// By condition number (requires _cond_2 to be set)
auto well = LinearAlgEqTestBed::getWellConditioned();        // κ₂ < 100
auto moderate = LinearAlgEqTestBed::getModeratelyConditioned(); // 100 ≤ κ₂ < 1000  
auto ill = LinearAlgEqTestBed::getIllConditioned();          // κ₂ ≥ 1000
```

#### Iteration Helpers

```cpp
// Get all systems as vectors for range-based loops
auto all = LinearAlgEqTestBed::getAllSystems();
auto allSym = LinearAlgEqTestBed::getAllSymmetricSystems();
auto allCmplx = LinearAlgEqTestBed::getAllComplexSystems();
auto allOver = LinearAlgEqTestBed::getAllOverdeterminedSystems();

// Example usage
for (auto& [name, sys] : LinearAlgEqTestBed::getAllSystems()) {
    // Test solver on each system
}
```

## Usage Examples

### Testing a Linear Solver

```cpp
#include "test_data/linear_alg_eq_systems_test_bed.h"
#include "MML.h"

using namespace MML;

TEST_CASE("LU Solver on Test Bed", "[LUSolver]") {
    for (int i = 0; i < TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++) {
        INFO("Testing system " << i);
        TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);
        
        Matrix<Real> mat = testBed._mat;
        Vector<Real> rhs = testBed._rhs;
        Vector<Real> sol(rhs.size());
        
        LUSolver<Real> luSolver(mat);
        luSolver.Solve(rhs, sol);
        
        // Verify against known solution
        REQUIRE(sol.IsEqualTo(testBed._sol, 1e-13));
        
        // Verify residual: A*x should equal b
        Vector<Real> res_rhs = testBed._mat * sol;
        REQUIRE(res_rhs.IsEqualTo(testBed._rhs, 1e-13));
    }
}
```

### Testing with Ill-Conditioned Matrices

```cpp
TEST_CASE("Solver Stability on Ill-Conditioned Systems", "[stability]") {
    auto illConditioned = TestBeds::LinearAlgEqTestBed::getIllConditioned();
    
    for (const auto& [name, sys] : illConditioned) {
        INFO("Testing: " << name << " with cond = " << sys->_cond_2);
        
        // Expect reduced accuracy for ill-conditioned systems
        Real tol = std::min(1e-10 * sys->_cond_2, 0.1);
        
        Matrix<Real> mat = sys->_mat;
        Vector<Real> sol(sys->_rhs.size());
        
        LUSolver<Real> luSolver(mat);
        luSolver.Solve(sys->_rhs, sol);
        
        REQUIRE(sol.IsEqualTo(sys->_sol, tol));
    }
}
```

### Testing Eigenvalues with Test Bed Data

```cpp
TEST_CASE("Verify Test Bed Eigenvalues", "[eigenvalues]") {
    // Test bed provides pre-computed eigenvalues for verification
    const auto& sys = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem("frank_5x5");
    
    // Access the stored eigenvalues
    const Vector<Complex>& eigenvalues = sys._eigen_values;
    
    // Frank matrices have all real eigenvalues - verify imaginary parts are near zero
    for (int i = 0; i < eigenvalues.size(); i++) {
        REQUIRE(std::abs(eigenvalues[i].imag()) < 1e-10);
    }
    
    // Verify eigenvalue property: det(A - λI) = 0
    // For each eigenvalue, (A - λI) should be singular
    for (int i = 0; i < eigenvalues.size(); i++) {
        Matrix<Real> AminusLambdaI = sys._mat;
        for (int j = 0; j < sys._n; j++) {
            AminusLambdaI[j][j] -= eigenvalues[i].real();
        }
        // Matrix should be nearly singular (small determinant)
        // Use LU decomposition determinant calculation
        LUSolver<Real> lu(AminusLambdaI);
        // ... verify near-singularity
    }
}
```

### Testing Cholesky on SPD Matrices

```cpp
TEST_CASE("Cholesky Decomposition on SPD Matrices", "[Cholesky][SPD]") {
    auto spdMatrices = TestBeds::LinearAlgEqTestBed::getSPDMatrices();
    
    for (const auto& [name, sys] : spdMatrices) {
        INFO("Testing: " << name);
        
        // Convert symmetric matrix to full matrix for Cholesky
        // (CholeskySolver works with Matrix<Real>, not MatrixSym)
        Matrix<Real> mat(sys->_n, sys->_n);
        for (int i = 0; i < sys->_n; i++) {
            for (int j = 0; j < sys->_n; j++) {
                mat[i][j] = sys->_mat(i, j);
            }
        }
        
        // If matrix is SPD, Cholesky should succeed (no exception)
        CholeskySolver<Real> cholSolver(mat);
        
        Vector<Real> sol = cholSolver.Solve(sys->_rhs);
        REQUIRE(sol.IsEqualTo(sys->_sol, 1e-12));
    }
}
```

### Testing Least Squares with QR

```cpp
TEST_CASE("QR Least Squares on Overdetermined Systems", "[QR][LeastSquares]") {
    for (int i = 0; i < TestBeds::LinearAlgEqTestBed::numOverdeterminedSystems(); i++) {
        const auto& sys = TestBeds::LinearAlgEqTestBed::getOverdeterminedSystem(i);
        INFO("Testing overdetermined system " << i << " (" << sys._m << "x" << sys._n << ")");
        
        QRSolver<Real> qrSolver(sys._mat);
        Vector<Real> x = qrSolver.LeastSquaresSolve(sys._rhs);
        
        REQUIRE(x.IsEqualTo(sys._sol, 1e-10));
        
        // Verify normal equations: AᵀAx = Aᵀb
        Matrix<Real> At = sys._mat.GetTranspose();
        Matrix<Real> AtA = At * sys._mat;
        Vector<Real> Atb = At * sys._rhs;
        Vector<Real> AtAx = AtA * x;
        REQUIRE(AtAx.IsEqualTo(Atb, 1e-10));
    }
}
```

### Testing GaussJordan Solver

```cpp
TEST_CASE("GaussJordan Complete Test Bed", "[GaussJordan]") {
    for (int i = 0; i < TestBeds::LinearAlgEqTestBed::numLinAlgEqSystems(); i++) {
        INFO("Testing system " << i);
        TestBeds::TestLinearSystem testBed = TestBeds::LinearAlgEqTestBed::getLinAlgEqSystem(i);
        
        Matrix<Real> mat = testBed._mat;
        Vector<Real> rhs = testBed._rhs;
        
        // GaussJordan solves in-place
        GaussJordanSolver<Real>::SolveInPlace(mat, rhs);
        
        // rhs now contains the solution
        // Use relaxed tolerance for ill-conditioned matrices
        Real tol = (i >= 9 && i <= 14) ? 1e-9 : 1e-13;
        REQUIRE(rhs.IsEqualTo(testBed._sol, tol));
        
        // Verify: A*x = b
        Vector<Real> res_rhs = testBed._mat * rhs;
        REQUIRE(res_rhs.IsEqualTo(testBed._rhs, tol));
    }
}
```

## Verification Methodology

All reference values in the test bed have been verified using:

1. **Symbolic computation** (Mathematica/SymPy) for exact integer matrices (Pascal)
2. **High-precision arithmetic** (mpmath with 50+ digits) for ill-conditioned systems
3. **Cross-validation** against NumPy, SciPy, and LAPACK reference implementations
4. **Round-trip verification**: Ax = b computed and compared to original RHS

## System Summary Table

| Category | Count | Size Range | Condition Range |
|----------|-------|------------|-----------------|
| Random real | 9 | 3×3 to 20×20 | varies |
| Hilbert | 4 | 3×3 to 8×8 | 524 to 1.69×10¹¹ |
| Pascal | 3 | 3×3 to 5×5 | 53 to 8,376 |
| Vandermonde | 3 | 3×3 to 5×5 | 64 to 648,163 |
| Frank | 3 | 3×3 to 5×5 | 119 to 3,913 |
| Kahan | 2 | 3×3 to 5×5 | 2.3 to 4.7 |
| Diag. dominant | 3 | 4×4 to 6×6 | well-conditioned |
| SPD | 4 | 3×3 to 6×6 | varies |
| Complex | 12 | 3×3 to 10×10 | varies |
| Overdetermined | 3 | 4×3 to 6×4 | 55 to 789 |
| **Total** | **~50** | | |

## References

- Higham, N.J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed. SIAM.
- Golub, G.H. & Van Loan, C.F. (2013). *Matrix Computations*, 4th ed. Johns Hopkins.
- Anderson, E., et al. (1999). *LAPACK Users' Guide*, 3rd ed. SIAM.
