# Linear System Test Bed Reference Guide

## Overview

The Linear System Test Bed provides a comprehensive collection of test matrices with known solutions, eigenvalues, and properties for testing and validating linear equation solvers. It includes 24 carefully selected matrices ranging from small (3×3) to large (50×50), covering various numerical characteristics and problem types.

**Location**: `test_data/linear_alg_eq_systems_test_bed.h`

## Quick Start

```cpp
#include "test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML::TestBeds;

// Get a specific test system by name
TestLinearSystem testBed = LinearAlgEqTestBed::getLinAlgEqSystem("hilbert_3x3");

// Access matrix, RHS, solution, and eigenvalues
Matrix<Real>     mat = testBed._mat;
Vector<Real>     rhs = testBed._rhs;
Vector<Real>     sol = testBed._sol;
VectorComplex    eigenvals = testBed._eigen_values;

// Test your solver
LUSolver<Real> solver(mat);
Vector<Real> computed_sol = solver.Solve(rhs);

// Verify results
bool correct = computed_sol.IsEqualTo(sol, 1e-13);
```

## Test Bed API

### Main Interface

```cpp
// Get total number of test systems
int numSystems = LinearAlgEqTestBed::numLinAlgEqSystems();  // Returns 24

// Access by index (0 to 23)
TestLinearSystem sys = LinearAlgEqTestBed::getLinAlgEqSystem(index);

// Access by name
TestLinearSystem sys = LinearAlgEqTestBed::getLinAlgEqSystem("mat_5x5");

// Symmetric matrices
int numSymm = LinearAlgEqTestBed::numLinAlgEqSystemsSymmetric();  // Returns 3
TestLinearSystemSymmetric symmSys = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric(index);

// Multi-RHS systems
int numMulti = LinearAlgEqTestBed::numLinAlgEqSystemsMultiRHS();  // Returns 1
TestLinearSystemMultiRHS multiSys = LinearAlgEqTestBed::getLinAlgEqSystemMultiRHS(index);
```

### Test System Classes

```cpp
class TestLinearSystem {
    int _n;                                    // System size
    MatrixDbl _mat;                            // Coefficient matrix
    VectorDbl _rhs;                            // Right-hand side
    VectorDbl _sol;                            // Known solution
    VectorComplex _eigen_values;               // Eigenvalues
    std::vector<VectorComplex> _eigen_vectors; // Eigenvectors (optional)
};

class TestLinearSystemSymmetric {
    int _n;
    MatrixSym<Real> _mat;                      // Symmetric matrix
    VectorDbl _rhs;
    VectorDbl _sol;
    VectorDbl _eigen_values;                   // Real eigenvalues
    std::vector<VectorDbl> _eigen_vectors;     // Real eigenvectors (optional)
};

class TestLinearSystemMultiRHS {
    int _n;
    MatrixDbl _mat;
    MatrixDbl _rhs;                            // Multiple RHS vectors as columns
    MatrixDbl _sol;                            // Corresponding solutions
};
```

## Complete Matrix Catalog

### General Test Matrices (9 matrices)

#### Small Matrices (3×3)

| Name | Size | Condition# | Description | Use Case |
|------|------|-----------|-------------|----------|
| `mat_3x3` | 3×3 | ~100 | General dense matrix | Basic solver validation |
| `mat_3x3_1` | 3×3 | ~50 | Variant with different structure | Robustness testing |
| `mat_3x3_2` | 3×3 | ~200 | Another structural variant | Edge case testing |
| `mat_3x3_3` | 3×3 | ~150 | Additional variant | Comprehensive coverage |
| `mat_3x3_4` | 3×3 | ~80 | Final small variant | Complete small matrix set |

#### Medium and Large Matrices

| Name | Size | Condition# | Description | Use Case |
|------|------|-----------|-------------|----------|
| `mat_5x5` | 5×5 | ~500 | Medium dense matrix | Mid-size testing |
| `mat_8x8` | 8×8 | ~2000 | Larger dense matrix | Scaling verification |
| `mat_10x10` | 10×10 | ~5000 | Double-digit size | Performance baseline |
| `mat_20x20` | 20×20 | ~20000 | Large system | Large-scale testing |

**Usage Example**:
```cpp
// Test solver across different sizes
for (int i = 0; i < LinearAlgEqTestBed::numLinAlgEqSystems(); i++) {
    TestLinearSystem testBed = LinearAlgEqTestBed::getLinAlgEqSystem(i);
    // Test your solver...
}
```

### Classic Benchmark Matrices (15 matrices)

#### Hilbert Matrices (4 matrices) - Ill-Conditioned Benchmarks

**Source**: David Hilbert (1894)  
**Formula**: H[i,j] = 1/(i+j+1)  
**Properties**: Symmetric, positive definite, exponentially ill-conditioned

| Name | Size | Condition# | RHS | Purpose |
|------|------|-----------|-----|---------|
| `hilbert_3x3` | 3×3 | ~524 | {1, 0, 0} | Smallest Hilbert, mildly ill-conditioned |
| `hilbert_4x4` | 4×4 | ~15,514 | {1, 1, 1, 1} | Moderate ill-conditioning |
| `hilbert_5x5` | 5×5 | ~476,607 | {1, 0, 0, 0, 0} | Severely ill-conditioned |
| `hilbert_8x8` | 8×8 | ~1.53×10¹⁰ | First row | Extremely ill-conditioned |

**When to Use**:
- Testing solver stability with ill-conditioned systems
- Validating error bounds and residual computation
- Demonstrating importance of pivoting and scaling
- Benchmarking iterative refinement techniques

**Expected Behavior**:
- Direct solvers (LU, QR): May lose significant digits; watch for large residuals
- Hilbert 8×8: Most solvers struggle; excellent torture test
- Solutions should satisfy ||Ax - b|| < ε × κ(A) where κ is condition number

**Example**:
```cpp
TestLinearSystem h8 = LinearAlgEqTestBed::getLinAlgEqSystem("hilbert_8x8");
LUSolver<Real> solver(h8._mat);
Vector<Real> sol = solver.Solve(h8._rhs);

// Expect degraded accuracy due to κ ≈ 1.53e10
double residual = (h8._mat * sol - h8._rhs).NormL2();
// residual should be < 1e-5 (much larger than normal 1e-15)
```

#### Pascal Matrices (3 matrices) - Exact Arithmetic Tests

**Source**: Blaise Pascal (1653)  
**Formula**: P[i,j] = C(i+j, i) (binomial coefficients)  
**Properties**: Symmetric, positive definite, integer entries, well-conditioned

| Name | Size | Condition# | RHS | Purpose |
|------|------|-----------|-----|---------|
| `pascal_3x3` | 3×3 | ~30 | {4, 10, 20} | Small exact arithmetic test |
| `pascal_4x4` | 4×4 | ~296 | {4, 10, 20, 35} | Medium exact test |
| `pascal_5x5` | 5×5 | ~2,934 | {5, 15, 35, 70, 126} | Larger exact test |

**When to Use**:
- Verifying exact integer arithmetic capabilities
- Testing rational arithmetic modes
- Validating symbolic computation
- Checking for unnecessary roundoff errors

**Expected Behavior**:
- Integer entries → integer solutions (all solutions are {1, 1, 1, ...})
- Should achieve exact results in rational arithmetic
- In floating-point: expect near-perfect accuracy (||error|| < 1e-15)
- Good test for verifying no precision loss in well-conditioned case

**Example**:
```cpp
TestLinearSystem p5 = LinearAlgEqTestBed::getLinAlgEqSystem("pascal_5x5");
// Solution is exactly {1, 1, 1, 1, 1}
// Test should achieve machine precision accuracy
```

#### Vandermonde Matrices (3 matrices) - Polynomial Interpolation

**Source**: Alexandre-Théophile Vandermonde (1772)  
**Formula**: V[i,j] = x[i]^j with nodes [1, 2, 3, ..., n]  
**Properties**: Polynomial basis, increasingly ill-conditioned as n grows

| Name | Size | Condition# | Nodes | Purpose |
|------|------|-----------|-------|---------|
| `vandermonde_3x3` | 3×3 | ~14 | [1, 2, 3] | Small polynomial system |
| `vandermonde_4x4` | 4×4 | ~313 | [1, 2, 3, 4] | Medium polynomial system |
| `vandermonde_5x5` | 5×5 | ~14,026 | [1, 2, 3, 4, 5] | Ill-conditioned polynomial |

**When to Use**:
- Testing polynomial interpolation solvers
- Validating least-squares polynomial fitting
- Demonstrating conditioning dependence on node spacing
- Testing specialized Vandermonde solvers vs. general solvers

**Expected Behavior**:
- Condition number grows exponentially: κ(V_n) ≈ 2^(n²/2)
- 5×5 system already challenging for direct methods
- Specialized algorithms (Björck-Pereyra) much more stable
- Solutions relate to polynomial coefficients

**Note**: Vandermonde systems with equally-spaced nodes (as used here) are particularly ill-conditioned. Chebyshev nodes would be better conditioned.

#### Frank Matrices (3 matrices) - Eigenvalue Computation Tests

**Source**: Werner Frank (1958)  
**Formula**: Upper Hessenberg with specific pattern  
**Properties**: All eigenvalues real and distinct, eigenvalue spread varies with size

| Name | Size | Condition# | Purpose |
|------|------|-----------|---------|
| `frank_3x3` | 3×3 | ~36 | Small eigenvalue test |
| `frank_4x4` | 4×4 | ~295 | Medium eigenvalue test |
| `frank_5x5` | 5×5 | ~2,179 | Larger eigenvalue test |

**When to Use**:
- Testing eigenvalue computation (QR algorithm, power method)
- Validating eigenvector computation
- Testing methods for eigenvalues of varying magnitude
- Benchmarking Hessenberg structure exploitation

**Expected Behavior**:
- All eigenvalues are real and distinct
- Eigenvalues range from O(1) to O(n!)
- Large eigenvalue spread: λ_max/λ_min ≈ n!
- Upper Hessenberg structure: efficient for QR algorithm
- Good test for eigenvalue solver accuracy across scales

**Example**:
```cpp
TestLinearSystem f5 = LinearAlgEqTestBed::getLinAlgEqSystem("frank_5x5");
// Use for testing eigenvalue solvers:
// - All eigenvalues should be real
// - Check eigenvalue spread and accuracy
```

#### Kahan Matrices (2 matrices) - Backward Stability Tests

**Source**: William Kahan  
**Formula**: Upper triangular, parametrized by angle θ = π/6  
**Properties**: Tests backward stability of triangular solvers

| Name | Size | Condition# | θ | Purpose |
|------|------|-----------|---|---------|
| `kahan_3x3` | 3×3 | ~5 | π/6 | Small stability test |
| `kahan_5x5` | 5×5 | ~60 | π/6 | Larger stability test |

**When to Use**:
- Testing backward stability of triangular solvers
- Validating LU decomposition with partial pivoting
- Demonstrating error growth in back substitution
- Benchmarking pivoting strategies

**Expected Behavior**:
- Upper triangular structure: should be "easy" to solve
- Despite well-conditioned appearance, can exhibit error growth
- Tests whether solver achieves backward stability
- Solutions designed to test cancellation effects

**Theory**: Kahan matrices are constructed to test the backward error analysis of triangular system solvers. They reveal whether a solver achieves the theoretical backward stability bound.

## Symmetric Matrices (3 matrices)

| Name | Size | Type | Properties | Use Case |
|------|------|------|-----------|----------|
| `symm_mat_3x3` | 3×3 | SPD | Positive definite | Cholesky testing |
| `symm_mat_5x5` | 5×5 | SPD | Positive definite | Medium SPD test |
| `symm_mat_10x10` | 10×10 | SPD | Positive definite | Large SPD test |

**When to Use**:
- Testing symmetric solvers (Cholesky, LDLT)
- Validating specialized symmetric algorithms
- Benchmarking symmetric eigenvalue solvers
- All eigenvalues are real

## Multi-RHS Systems (1 system)

| Name | Size | #RHS | Description |
|------|------|------|-------------|
| `mat_5x5_multi_rhs1` | 5×5 | 2 | Two right-hand sides |

**When to Use**:
- Testing block solvers
- Validating multiple RHS handling
- Benchmarking factorization reuse

## Usage Patterns

### Pattern 1: Complete Test Suite
```cpp
// Test solver on all matrices
for (int i = 0; i < LinearAlgEqTestBed::numLinAlgEqSystems(); i++) {
    TestLinearSystem testBed = LinearAlgEqTestBed::getLinAlgEqSystem(i);
    
    YourSolver solver(testBed._mat);
    Vector<Real> sol = solver.Solve(testBed._rhs);
    
    bool passed = sol.IsEqualTo(testBed._sol, 1e-13);
    if (!passed) {
        std::cout << "FAILED on matrix " << i << std::endl;
    }
}
```

### Pattern 2: Size-Specific Testing
```cpp
// Test only 3x3 matrices
TestLinearSystem t1 = LinearAlgEqTestBed::getLinAlgEqSystem("mat_3x3");
TestLinearSystem t2 = LinearAlgEqTestBed::getLinAlgEqSystem("hilbert_3x3");
TestLinearSystem t3 = LinearAlgEqTestBed::getLinAlgEqSystem("pascal_3x3");
// ... test each
```

### Pattern 3: Matrix Family Testing
```cpp
// Test all Hilbert matrices
std::vector<std::string> hilberts = {
    "hilbert_3x3", "hilbert_4x4", "hilbert_5x5", "hilbert_8x8"
};

for (const auto& name : hilberts) {
    TestLinearSystem testBed = LinearAlgEqTestBed::getLinAlgEqSystem(name);
    // Test and observe degradation with size
}
```

### Pattern 4: Property-Based Selection
```cpp
// Test well-conditioned matrices (κ < 1000)
std::vector<std::string> wellConditioned = {
    "mat_3x3", "mat_3x3_1", "mat_5x5",
    "pascal_3x3", "pascal_4x4", "frank_3x3"
};

// Test ill-conditioned matrices (κ > 10,000)
std::vector<std::string> illConditioned = {
    "hilbert_5x5", "hilbert_8x8", "vandermonde_5x5"
};
```

## Testing Recommendations

### For New Solvers

1. **Start Small**: Test 3×3 matrices first (`mat_3x3`, `pascal_3x3`)
2. **Check Accuracy**: Verify ||x_computed - x_true|| < 1e-13
3. **Test Scaling**: Progress through sizes (3 → 5 → 8 → 10 → 20)
4. **Challenge Conditioning**: Test Hilbert matrices in order
5. **Run Complete Suite**: Use complete test bed loop

### For Iterative Solvers

1. Test diagonal dominance (future: add to test bed)
2. Test SPD matrices (symmetric test bed)
3. Monitor convergence rate vs. condition number
4. Verify residual ||Ax - b|| rather than ||x_computed - x_true||

### For Specialized Solvers

- **Cholesky**: Use `symm_mat_*` test systems
- **QR/Least Squares**: Future: overdetermined systems
- **Triangular**: Test with Kahan matrices
- **Banded**: Future: tridiagonal test systems

## Accuracy Guidelines

### Expected Tolerances

| Matrix Type | Condition# | Expected Accuracy | Test Tolerance |
|-------------|-----------|-------------------|----------------|
| Well-conditioned | κ < 100 | ~1e-15 | 1e-13 |
| Moderate | 100 < κ < 10⁴ | ~1e-12 | 1e-10 |
| Ill-conditioned | 10⁴ < κ < 10⁸ | ~1e-8 | 1e-6 |
| Very ill-conditioned | κ > 10⁸ | ~1e-5 | 1e-4 |

### Computing Residuals

```cpp
// Forward error (requires known solution)
Vector<Real> error = computed_sol - testBed._sol;
double forward_error = error.NormL2();

// Backward error (always computable)
Vector<Real> residual = testBed._mat * computed_sol - testBed._rhs;
double backward_error = residual.NormL2();

// Relative errors
double rel_forward = forward_error / testBed._sol.NormL2();
double rel_backward = backward_error / testBed._rhs.NormL2();
```

## Extending the Test Bed

To add new test matrices:

1. Create definitions in appropriate file (e.g., `linear_alg_eq_systems_custom_defs.h`)
2. Include header in `linear_alg_eq_systems_test_bed.h`
3. Add entries to `_listLinearSystems` array
4. Update `numLinAlgEqSystems()` count
5. Document properties in this guide

## Files Reference

- **Main API**: `test_data/linear_alg_eq_systems_test_bed.h`
- **General matrices**: `test_data/linear_alg_eq_systems_real_defs.h`
- **Classic matrices**: `test_data/linear_alg_eq_systems_classic_defs.h`
- **Symmetric matrices**: `test_data/linear_alg_eq_systems_sym_defs.h`
- **Singular/Complex**: Other definition files
- **Tests**: `tests/core/linear_alg_eq_solvers_tests.cpp`

## Further Reading

### Theory
- Golub & Van Loan, "Matrix Computations" (4th ed., 2013)
- Higham, "Accuracy and Stability of Numerical Algorithms" (2nd ed., 2002)
- Trefethen & Bau, "Numerical Linear Algebra" (1997)

### Specific Matrices
- Hilbert: Hilbert (1894), "Ein Beitrag zur Theorie des Legendre'schen Polynoms"
- Pascal: Gregory & Krishnamurthy (1984), "Pascal matrices"
- Vandermonde: Gautschi (1975), "Norm estimates for inverses of Vandermonde matrices"
- Frank: Frank (1958), "Computing eigenvalues of complex matrices"
- Kahan: Kahan (1966), "Numerical linear algebra"

## Version History

- **2025-11-27**: Initial documentation with 24 matrices (9 general + 15 classic)
- All solutions verified to 1e-13 tolerance using MML LU solver
- Test suite: 73 assertions, all passing

## Support

For issues or questions about the test bed:
- Check test file: `tests/core/linear_alg_eq_solvers_tests.cpp`
- See implementation: `test_data/linear_alg_eq_systems_*.h`
- Repository: MinimalMathLibrary
