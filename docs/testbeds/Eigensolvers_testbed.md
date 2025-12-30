# Eigenvalue Solvers Test Bed Documentation

## Overview

The Eigenvalue Test Bed provides comprehensive test infrastructure for validating MML's eigenvalue/eigenvector solvers. It includes symmetric matrices, defective matrices, matrices with complex eigenvalues, and challenging pathological cases with known analytical eigenvalues.

**Key Features:**
- **27 test cases** across 8 categories
- Matrices with **known eigenvalues** and eigenvectors
- Metadata: symmetry, defectiveness, condition number, eigenvalue spread
- Support for **real and complex eigenvalues**
- Tests for **repeated, clustered, and pathological** cases

---

## Core Files

| File | Purpose |
|------|---------|
| [eigenvalue_test_bed.h](../../test_data/eigenvalue_test_bed.h) | `TestEigenSystem` struct, retrieval functions, verification utilities |
| [eigenvalue_defs.h](../../test_data/eigenvalue_defs.h) | All test matrix definitions and known eigenvalues |

---

## Data Structure

```cpp
namespace MML::TestBeds {

/// Test case for eigenvalue solvers
struct TestEigenSystem {
    std::string name;                                  ///< Descriptive name
    std::function<Matrix<Real>()> getMatrix;          ///< Factory to get test matrix
    
    // Eigenvalue data
    Vector<Real> realEigenvalues;                     ///< Real eigenvalues (if all real)
    std::vector<std::complex<Real>> complexEigenvalues;  ///< Complex eigenvalues
    std::vector<Vector<Real>> eigenvectors;           ///< Eigenvectors (if known)
    
    // Matrix properties
    int dimension;                                     ///< Matrix dimension n×n
    bool isSymmetric = false;                         ///< Symmetric ⟹ real eigenvalues
    bool isPositiveDefinite = false;                  ///< All eigenvalues positive
    bool hasRealEigenvalues = true;                   ///< All eigenvalues are real
    bool isDefective = false;                         ///< Fewer eigenvectors than eigenvalues
    bool hasRepeatedEigenvalues = false;              ///< Multiple eigenvalues exist
    bool hasClusteredEigenvalues = false;             ///< Eigenvalues very close together
    
    // Numerical properties
    Real conditionNumber = 1.0;                       ///< Condition number (if known)
    Real eigenvalueSpread = 1.0;                      ///< Ratio of |λ_max|/|λ_min|
    int multiplicityMax = 1;                          ///< Maximum algebraic multiplicity
    
    // Metadata
    std::string category;                             ///< standard, repeated, defective, etc.
    std::string description;                          ///< Detailed description
    int difficulty = 1;                               ///< 1=easy, 2=medium, 3=hard, 4=extreme
};

}
```

---

## Test Cases Catalog (27 total)

### Category 1: Standard/Classic Matrices (3)

Well-conditioned reference cases with known eigenvalues.

| # | Name | Dim | Type | Eigenvalues | Condition | Difficulty |
|---|------|-----|------|-------------|-----------|------------|
| 1 | **Pascal 5×5** | 5 | Sym PD | 0.011, 0.18, 1.0, 5.5, 92.3 | ~8500 | Easy |
| 2 | **Hilbert 10×10** | 10 | Sym PD | 1.8e-10 to 1.75 | ~10¹³ | Hard |
| 3 | **Toeplitz Tridiag 20×20** | 20 | Sym PD | 2+2cos(kπ/21) | Low | Easy |

### Category 2: Repeated Eigenvalues (3)

Tests handling of multiplicity.

| # | Name | Dim | Eigenvalues | Max Multiplicity |
|---|------|-----|-------------|------------------|
| 4 | **Repeated Double 4×4** | 4 | 1, 2, 2, 3 | 2 |
| 5 | **Repeated Triple 5×5** | 5 | 1, 2, 2, 2, 3 | 3 |
| 6 | **Repeated Quadruple 4×4** | 4 | 0, 2, 2, 4 | 2 |

### Category 3: Defective Matrices (4)

Non-diagonalizable matrices with Jordan block structure.

| # | Name | Dim | Eigenvalues | Structure | Difficulty |
|---|------|-----|-------------|-----------|------------|
| 7 | **Jordan Block 2×2** | 2 | 2, 2 | J₂(2) | Hard |
| 8 | **Jordan Block 3×3** | 3 | 1, 1, 1 | J₃(1) | Hard |
| 9 | **Mixed Jordan 3×3** | 3 | 2, 2, 3 | J₂(2) ⊕ J₁(3) | Hard |
| 10 | **Two Jordan Blocks 4×4** | 4 | 1, 1, 2, 2 | J₂(1) ⊕ J₂(2) | Hard |

### Category 4: Nearly Defective Matrices (2)

Ill-conditioned eigenvector problems.

| # | Name | Dim | Description | Difficulty |
|---|------|-----|-------------|------------|
| 11 | **Nearly Defective 3×3** | 3 | Perturbed Jordan, ε=10⁻⁸ | Extreme |
| 12 | **Wilkinson W⁻ 21×21** | 21 | Famous close eigenvalue pairs | Extreme |

### Category 5: Clustered Eigenvalues (3)

Distinct but very close eigenvalues.

| # | Name | Dim | Cluster Description | Difficulty |
|---|------|-----|---------------------|------------|
| 13 | **Clustered Tight 4×4** | 4 | 0.999, 1.0, 1.001, 2.0 | Medium |
| 14 | **Clustered Very Tight 5×5** | 5 | 1.0, 1.0001, 1.0002, 1.0003, 2.0 | Hard |
| 15 | **Graded Clusters 6×6** | 6 | Groups at scales 0.1, 1, 10 | Medium |

### Category 6: Complex Eigenvalues (4)

Non-symmetric matrices with complex eigenvalues.

| # | Name | Dim | Eigenvalues | Description |
|---|------|-----|-------------|-------------|
| 16 | **Rotation 2×2** | 2 | e^{±iθ} | Pure rotation, |λ|=1 |
| 17 | **Complex Pair 3×3** | 3 | 1±i, 3 | One real + conjugate pair |
| 18 | **Two Complex Pairs 4×4** | 4 | ±i, ±2i | All purely imaginary |
| 19 | **Companion x⁴-1** | 4 | ±1, ±i | 4th roots of unity |

### Category 7: Special Structures (3)

Matrices with special algebraic properties.

| # | Name | Dim | Property | Eigenvalues |
|---|------|-----|----------|-------------|
| 20 | **Circulant 4×4** | 4 | Circulant | 10, -2, ±2i |
| 21 | **Block Diagonal 6×6** | 6 | Block diagonal | Union of block eigenvalues |
| 22 | **Stochastic 4×4** | 4 | Row stochastic | Dominant λ=1 |

### Category 8: Pathological Cases (5)

Extreme numerical challenges.

| # | Name | Dim | Challenge | Difficulty |
|---|------|-----|-----------|------------|
| 23 | **Nilpotent 3×3** | 3 | All λ=0 but A≠0 | Hard |
| 24 | **Tiny Eigenvalue 3×3** | 3 | λ = 10⁻¹⁵, 1, 2 | Hard |
| 25 | **Huge Spread 3×3** | 3 | λ = 10⁻¹⁵, 1, 10¹⁵ | Extreme |
| 26 | **Almost Singular 4×4** | 4 | Three eigenvalues ~10⁻¹² | Extreme |
| 27 | **Zero Matrix** | n | All λ=0 | Easy |

---

## Matrix Definitions

### Classic Reference Matrices

```cpp
// Pascal matrix (symmetric positive definite)
const Matrix<Real> classic_pascal_5x5{5, 5, {
     1,  1,  1,  1,  1,
     1,  2,  3,  4,  5,
     1,  3,  6, 10, 15,
     1,  4, 10, 20, 35,
     1,  5, 15, 35, 70
}};

// Hilbert matrix H(i,j) = 1/(i+j+1)
Matrix<Real> getHilbert10x10();

// Toeplitz tridiagonal: tridiag(1, 2, 1)
// Eigenvalues: λ_k = 2 + 2·cos(kπ/(n+1))
Matrix<Real> getToeplitz20x20();
```

### Defective (Jordan) Matrices

```cpp
// Classic 2×2 Jordan block
const Matrix<Real> defective_jordan_2x2{2, 2, {
    2, 1,
    0, 2
}};
// Eigenvalue 2 (multiplicity 2), but only ONE eigenvector

// 3×3 Jordan block - maximally defective
const Matrix<Real> defective_jordan_3x3{3, 3, {
    1, 1, 0,
    0, 1, 1,
    0, 0, 1
}};
```

### Complex Eigenvalue Matrices

```cpp
// Rotation matrix: eigenvalues e^{±iθ}
Matrix<Real> getRotation2x2(Real theta = π/4);

// One real + complex pair: λ = 1±i, 3
const Matrix<Real> complex_pair_3x3{3, 3, {
    1.0, -1.0, 0.0,
    1.0,  1.0, 0.0,
    0.0,  0.0, 3.0
}};
```

### Pathological Matrices

```cpp
// Nilpotent: A³=0, all eigenvalues=0, but A≠0
const Matrix<Real> nilpotent_3x3{3, 3, {
    0, 1, 0,
    0, 0, 1,
    0, 0, 0
}};

// Extreme eigenvalue spread
Matrix<Real> getHugeSpread3x3(Real scale = 1e15);
// Eigenvalues: 1/scale, 1, scale
```

---

## API Reference

### Retrieval Functions

```cpp
namespace MML::TestBeds {

// Get all tests
std::vector<TestEigenSystem> getAllEigenTests();           // All 27 tests

// By matrix property
std::vector<TestEigenSystem> getSymmetricEigenTests();     // Symmetric only
std::vector<TestEigenSystem> getDefectiveEigenTests();     // Defective matrices
std::vector<TestEigenSystem> getComplexEigenTests();       // Complex eigenvalues

// By difficulty
std::vector<TestEigenSystem> getEasyEigenTests();          // Difficulty 1
std::vector<TestEigenSystem> getChallengingEigenTests();   // Difficulty 3+

}
```

### Individual Test Getters

```cpp
namespace MML::TestBeds {

// Standard/Classic
TestEigenSystem getPascal5x5Test();
TestEigenSystem getHilbert10x10Test();
TestEigenSystem getToeplitz20x20Test();

// Repeated eigenvalues
TestEigenSystem getRepeatedDouble4x4Test();
TestEigenSystem getRepeatedTriple5x5Test();
TestEigenSystem getRepeatedQuadruple4x4Test();

// Defective
TestEigenSystem getDefectiveJordan2x2Test();
TestEigenSystem getDefectiveJordan3x3Test();
TestEigenSystem getDefectiveMixed3x3Test();
TestEigenSystem getDefectiveTwoBlocks4x4Test();

// Nearly defective
TestEigenSystem getNearlyDefective3x3Test(Real epsilon = 1e-8);
TestEigenSystem getWilkinsonTest(int n = 21);

// Clustered
TestEigenSystem getClusteredTight4x4Test();
TestEigenSystem getClusteredVeryTight5x5Test();
TestEigenSystem getGradedClusters6x6Test();

// Complex eigenvalues
TestEigenSystem getRotation2x2Test(Real theta = π/4);
TestEigenSystem getComplexPair3x3Test();
TestEigenSystem getComplexTwoPairs4x4Test();
TestEigenSystem getCompanionTest();

// Special structure
TestEigenSystem getCirculant4x4Test();
TestEigenSystem getBlockDiagonal6x6Test();
TestEigenSystem getStochastic4x4Test();

// Pathological
TestEigenSystem getNilpotent3x3Test();
TestEigenSystem getTinyEigenvalue3x3Test(Real eps = 1e-15);
TestEigenSystem getHugeSpread3x3Test(Real scale = 1e15);
TestEigenSystem getAlmostSingular4x4Test(Real eps = 1e-12);

}
```

### Verification Utilities

```cpp
namespace MML::TestBeds {

/// Verify eigenvalues match expected (sorts both before comparison)
bool verifyRealEigenvalues(Vector<Real> computed, Vector<Real> expected, 
                            Real tolerance = 1e-10);

/// Compute max absolute eigenvalue error
Real computeEigenvalueError(Vector<Real> computed, Vector<Real> expected);

/// Verify eigenvector orthonormality (symmetric matrices)
bool verifyOrthonormality(const Matrix<Real>& V, Real tolerance = 1e-10);

/// Compute residual ||A*v - λ*v|| for eigenpair
Real computeResidual(const Matrix<Real>& A, const Vector<Real>& v, Real lambda);

}
```

---

## Usage Examples

### Example 1: Test Jacobi Solver on Symmetric Matrices

```cpp
#include "test_data/eigenvalue_test_bed.h"
#include "algorithms/EigenSystemSolvers.h"

using namespace MML;
using namespace MML::TestBeds;

void testJacobiSolver() {
    auto tests = getSymmetricEigenTests();
    
    for (const auto& test : tests) {
        Matrix<Real> A = test.getMatrix();
        MatrixSym<Real> Asym(A);  // Convert to symmetric type
        
        auto result = SymmMatEigenSolverJacobi::Solve(Asym);
        
        if (result.converged) {
            Real error = computeEigenvalueError(result.eigenvalues, test.realEigenvalues);
            bool ortho = verifyOrthonormality(result.eigenvectors);
            
            std::cout << test.name << ": error = " << error
                      << ", orthonormal = " << (ortho ? "YES" : "NO")
                      << ", sweeps = " << result.iterations << std::endl;
        } else {
            std::cout << test.name << ": FAILED to converge" << std::endl;
        }
    }
}
```

### Example 2: Test QR Solver on Tridiagonal Matrix

```cpp
void testQRSolver() {
    auto test = getToeplitz20x20Test();
    Matrix<Real> A = test.getMatrix();
    MatrixSym<Real> Asym(A);
    
    auto result = SymmMatEigenSolverQR::Solve(Asym);
    
    Real error = computeEigenvalueError(result.eigenvalues, test.realEigenvalues);
    std::cout << "Toeplitz 20x20 (QR): max eigenvalue error = " << error << std::endl;
    
    // Verify A*V = V*D
    Matrix<Real>& V = result.eigenvectors;
    for (int i = 0; i < 20; ++i) {
        Vector<Real> v(20);
        for (int j = 0; j < 20; ++j) v[j] = V(j, i);
        
        Real residual = computeResidual(A, v, result.eigenvalues[i]);
        std::cout << "  λ[" << i << "] = " << result.eigenvalues[i]
                  << ", residual = " << residual << std::endl;
    }
}
```

### Example 3: Test General Eigensolver on Non-Symmetric Matrix

```cpp
void testGeneralEigensolver() {
    auto tests = getComplexEigenTests();
    
    for (const auto& test : tests) {
        Matrix<Real> A = test.getMatrix();
        
        auto result = EigenSolver::Solve(A);
        
        std::cout << test.name << " eigenvalues:\n";
        for (int i = 0; i < test.dimension; ++i) {
            auto& computed = result.eigenvalues[i];
            std::cout << "  λ[" << i << "] = " << computed.real() 
                      << " + " << computed.imag() << "i" << std::endl;
        }
    }
}
```

### Example 4: Verify Eigenvector Equation A*v = λ*v

```cpp
void verifyEigenpairs() {
    auto test = getPascal5x5Test();
    Matrix<Real> A = test.getMatrix();
    MatrixSym<Real> Asym(A);
    
    auto result = SymmMatEigenSolverJacobi::Solve(Asym);
    
    std::cout << "Pascal 5x5 eigenpair verification:\n";
    for (int i = 0; i < 5; ++i) {
        Vector<Real> v(5);
        for (int j = 0; j < 5; ++j) v[j] = result.eigenvectors(j, i);
        
        Real lambda = result.eigenvalues[i];
        Real residual = computeResidual(A, v, lambda);
        
        std::cout << "  λ = " << std::setw(12) << lambda 
                  << ", ||A*v - λ*v|| = " << residual << std::endl;
    }
}
```

### Example 5: Handle Repeated Eigenvalues

```cpp
void testRepeatedEigenvalues() {
    auto tests = {
        getRepeatedDouble4x4Test(),
        getRepeatedTriple5x5Test(),
        getRepeatedQuadruple4x4Test()
    };
    
    for (const auto& test : tests) {
        Matrix<Real> A = test.getMatrix();
        MatrixSym<Real> Asym(A);
        
        auto result = SymmMatEigenSolverJacobi::Solve(Asym);
        
        // Verify eigenvalues
        bool eigMatch = verifyRealEigenvalues(result.eigenvalues, test.realEigenvalues);
        
        // Verify orthonormality (crucial for repeated eigenvalues)
        bool ortho = verifyOrthonormality(result.eigenvectors);
        
        std::cout << test.name << ": eig_match = " << (eigMatch ? "YES" : "NO")
                  << ", orthonormal = " << (ortho ? "YES" : "NO")
                  << " (max multiplicity = " << test.multiplicityMax << ")" << std::endl;
    }
}
```

### Example 6: Test Ill-Conditioned Hilbert Matrix

```cpp
void testHilbertMatrix() {
    auto test = getHilbert10x10Test();
    Matrix<Real> A = test.getMatrix();
    MatrixSym<Real> Asym(A);
    
    // Jacobi is more stable for ill-conditioned matrices
    auto result = SymmMatEigenSolverJacobi::Solve(Asym, 1e-12, 200);
    
    std::cout << "Hilbert 10x10 (condition ~10^13):\n";
    std::cout << "  Converged: " << (result.converged ? "YES" : "NO") << std::endl;
    std::cout << "  Iterations: " << result.iterations << std::endl;
    
    // Compare eigenvalues
    for (int i = 0; i < 10; ++i) {
        Real computed = result.eigenvalues[i];
        Real expected = test.realEigenvalues[i];
        Real relError = std::abs(computed - expected) / std::abs(expected);
        
        std::cout << "  λ[" << i << "]: " << std::scientific << computed
                  << ", rel_error = " << relError << std::endl;
    }
}
```

### Example 7: Compare Jacobi vs QR Performance

```cpp
void compareJacobiVsQR() {
    auto tests = getSymmetricEigenTests();
    
    std::cout << std::setw(25) << "Matrix" 
              << std::setw(15) << "Jacobi (ms)"
              << std::setw(15) << "QR (ms)"
              << std::setw(15) << "Max Error" << std::endl;
    
    for (const auto& test : tests) {
        Matrix<Real> A = test.getMatrix();
        MatrixSym<Real> Asym(A);
        
        // Time Jacobi
        auto t1 = std::chrono::high_resolution_clock::now();
        auto jacobiResult = SymmMatEigenSolverJacobi::Solve(Asym);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto jacobiTime = std::chrono::duration<double, std::milli>(t2 - t1).count();
        
        // Time QR
        t1 = std::chrono::high_resolution_clock::now();
        auto qrResult = SymmMatEigenSolverQR::Solve(Asym);
        t2 = std::chrono::high_resolution_clock::now();
        auto qrTime = std::chrono::duration<double, std::milli>(t2 - t1).count();
        
        Real error = computeEigenvalueError(jacobiResult.eigenvalues, qrResult.eigenvalues);
        
        std::cout << std::setw(25) << test.name
                  << std::setw(15) << jacobiTime
                  << std::setw(15) << qrTime
                  << std::setw(15) << error << std::endl;
    }
}
```

---

## MML Eigensolver Reference

### Available Solvers

| Solver | Matrix Type | Algorithm | Eigenvalues | Eigenvectors |
|--------|-------------|-----------|-------------|--------------|
| `SymmMatEigenSolverJacobi` | Symmetric | Jacobi rotations | Real | Yes |
| `SymmMatEigenSolverQR` | Symmetric | Tridiag + QR | Real | Yes |
| `EigenSolver` | General | Hessenberg + QR | Complex | Yes |

### Algorithm Selection Guide

```
Matrix type?
├── Symmetric (A = Aᵀ)
│   ├── Small (n < 50)? → Jacobi (simple, stable)
│   └── Large (n ≥ 50)? → QR (faster asymptotically)
│
└── Non-symmetric
    └── Use EigenSolver (handles complex eigenvalues)
```

### Complexity

| Algorithm | Time | Space | Notes |
|-----------|------|-------|-------|
| **Jacobi** | O(n³) per sweep | O(n²) | 5-10 sweeps typical |
| **QR (symmetric)** | O(n³) total | O(n²) | Tridiag + iteration |
| **QR (general)** | O(n³) | O(n²) | Hessenberg + shifts |

---

## Test Categories by Difficulty

| Difficulty | Description | Examples |
|------------|-------------|----------|
| **Easy (1)** | Well-conditioned, small | Pascal, Toeplitz, Diagonal |
| **Medium (2)** | Repeated/clustered eigenvalues | Clustered, Rotation, Circulant |
| **Hard (3)** | Defective, ill-conditioned | Jordan blocks, Hilbert, Nilpotent |
| **Extreme (4)** | Near-singular, huge spread | Wilkinson, Nearly defective, Huge spread |

---

## Important Concepts

### Symmetric vs Non-Symmetric

| Property | Symmetric | General |
|----------|-----------|---------|
| Eigenvalues | Always real | May be complex |
| Eigenvectors | Orthonormal | May not be orthogonal |
| Diagonalizable | Always | May be defective |
| Recommended solver | Jacobi or QR | General EigenSolver |

### Defective Matrices

A matrix is **defective** if it has fewer linearly independent eigenvectors than its dimension:

```
Algebraic multiplicity (# times λ appears) > Geometric multiplicity (# eigenvectors for λ)
```

Example: Jordan block J₂(2) = [2, 1; 0, 2]
- Eigenvalue 2 with algebraic multiplicity 2
- Only 1 eigenvector: [1, 0]ᵀ
- Cannot be diagonalized!

### Condition Number and Accuracy

For eigenvalue computation, accuracy is limited by:

```
|Δλ| / |λ| ≈ κ(V) × ε_machine
```

where κ(V) is the condition number of the eigenvector matrix.

For symmetric matrices, κ(V) = 1 (eigenvectors are orthonormal), so eigenvalue accuracy is excellent.

---

## See Also

- [LinAlgSystems_testbed.md](LinAlgSystems_testbed.md) - Linear system test cases
- [Optimization_testbed.md](Optimization_testbed.md) - Optimization test cases
- [EigenSystemSolvers.h](../../mml/algorithms/EigenSystemSolvers.h) - MML eigensolver implementations
