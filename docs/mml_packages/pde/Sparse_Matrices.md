# Sparse Matrix Storage Formats

Comprehensive guide to sparse matrix representations in MML's PDE module.

## Table of Contents
1. [Introduction](#introduction)
2. [Why Sparse Matrices?](#why-sparse-matrices)
3. [COO Format](#coo-format-coordinate-list)
4. [CSR Format](#csr-format-compressed-sparse-row)
5. [CSC Format](#csc-format-compressed-sparse-column)
6. [Format Comparison](#format-comparison)
7. [Construction Patterns](#construction-patterns)
8. [Common Operations](#common-operations)
9. [PDE Assembly Example](#pde-assembly-example)
10. [Best Practices](#best-practices)

---

## Introduction

When discretizing PDEs using finite differences, the resulting linear systems produce **sparse matrices** - matrices where most entries are zero. For a 2D Poisson equation on an N×N grid:
- Matrix size: N² × N² (can be millions × millions!)
- Non-zeros per row: ~5 (5-point stencil)
- Sparsity: 99.99%+ for large grids

Storing all entries (including zeros) is prohibitively expensive. Sparse formats store only non-zero entries, achieving massive memory savings and computational speedups.

---

## Why Sparse Matrices?

### Memory Comparison

**Dense storage** (2D Laplacian, 100×100 grid):
- Matrix size: 10,000 × 10,000
- Entries: 100 million
- Memory (double): 800 MB

**Sparse storage** (CSR format):
- Non-zero entries: ~50,000
- Memory: 0.4 MB (2000× less!)

### Performance Comparison

**Matrix-vector product** (Ax = y):

| Format | Operation Count | Memory Access |
|--------|----------------|---------------|
| Dense  | O(n²)          | 100M reads    |
| Sparse CSR | O(nnz)     | 50K reads     |
| **Speedup** | **2000×** | **2000×**     |

For iterative solvers that perform hundreds of SpMV operations, this difference is transformative.

---

## COO Format (Coordinate List)

### Concept

COO stores each non-zero as a **(row, col, value)** triplet. Simple and flexible!

```
Example Matrix:        COO Storage:
┌           ┐          ┌─────┬─────┬──────┐
│ 2  0 -1 │          │ row │ col │value │
│ 0  5  3 │    →     ├─────┼─────┼──────┤
│-1  0  1 │          │  0  │  0  │  2.0 │
└           ┘          │  0  │  2  │ -1.0 │
                       │  1  │  1  │  5.0 │
                       │  1  │  2  │  3.0 │
                       │  2  │  0  │ -1.0 │
                       │  2  │  2  │  1.0 │
                       └─────┴─────┴──────┘
```

### Class: `SparseMatrixCOO<T>`

```cpp
#include "pde/sparse/SparseMatrixCOO.h"
using namespace MML::PDE::Sparse;

// Create empty 100×100 sparse matrix
SparseMatrixCOO<double> A(100, 100);

// Add entries (can be in any order!)
A.addEntry(0, 0, 2.0);
A.addEntry(0, 1, -1.0);
A.addEntry(1, 0, -1.0);
A.addEntry(1, 1, 2.0);
// ... more entries

// Query properties
int n = A.rows();       // 100
int m = A.cols();       // 100
int nnz = A.nnz();      // Number of non-zeros added
```

### When to Use COO

✅ **Ideal for:**
- Matrix construction/assembly
- Adding entries incrementally
- Finite difference stencil application
- FEM assembly loops
- Building from arbitrary order

❌ **Not efficient for:**
- Matrix-vector products (use CSR instead)
- Element access (slow linear search)
- Iterative solvers (convert to CSR first)

### Construction Methods

**Method 1: Individual entries**
```cpp
SparseMatrixCOO<double> A(n, n);

for (int i = 0; i < n; ++i) {
    if (i > 0)  A.addEntry(i, i-1, -1.0);  // Sub-diagonal
    A.addEntry(i, i, 2.0);                 // Diagonal
    if (i < n-1) A.addEntry(i, i+1, -1.0); // Super-diagonal
}
```

**Method 2: Batch construction**
```cpp
std::vector<Triplet<double>> triplets;
triplets.reserve(3 * n);

for (int i = 0; i < n; ++i) {
    if (i > 0)  triplets.emplace_back(i, i-1, -1.0);
    triplets.emplace_back(i, i, 2.0);
    if (i < n-1) triplets.emplace_back(i, i+1, -1.0);
}

A.addEntries(triplets);
```

**Method 3: Stencil helpers**
```cpp
// 1D Laplacian stencil: [-1, 2, -1] / h²
SparseMatrixCOO<double> A(n, n);
double h_sq_inv = 1.0 / (h * h);

for (int i = 0; i < n; ++i) {
    int row = i;
    A.addLaplacian1D(row, i, h_sq_inv);
}
```

### Advanced: Consolidation

COO can store duplicate entries. Consolidate to sum them:

```cpp
A.addEntry(0, 0, 1.0);
A.addEntry(0, 0, 2.0);  // Duplicate!

A.consolidate();  // Sums duplicates: (0,0) = 3.0
```

---

## CSR Format (Compressed Sparse Row)

### Concept

CSR stores non-zeros **row-by-row** in compressed form using three arrays:
- **`values`**: Non-zero values in row-major order
- **`col_indices`**: Column index for each value
- **`row_pointers`**: Where each row starts in values array

```
Example Matrix:        CSR Storage:
┌           ┐          values:       [ 2  -1  |  5   3  | -1   1 ]
│ 2  0 -1 │          col_indices:  [ 0   2  |  1   2  |  0   2 ]
│ 0  5  3 │    →     row_pointers: [ 0        2         4        6 ]
│-1  0  1 │                          ^row 0   ^row 1    ^row 2  ^end
└           ┘
```

**Row i spans** `row_pointers[i]` to `row_pointers[i+1]-1`.

### Class: `SparseMatrixCSR<T>`

```cpp
#include "pde/sparse/SparseMatrixCSR.h"
using namespace MML::PDE::Sparse;

// Method 1: Construct from COO (recommended!)
SparseMatrixCOO<double> coo(n, n);
// ... build COO ...
SparseMatrixCSR<double> A(coo);  // Automatic conversion

// Method 2: Direct construction (advanced)
SparseMatrixCSR<double> A(rows, cols);
// Requires knowledge of CSR internals
```

### When to Use CSR

✅ **Ideal for:**
- **Matrix-vector products** (Ax = y) - most important!
- Iterative solvers (CG, BiCGSTAB, GMRES)
- Row-wise operations
- Preconditioner application

❌ **Less efficient for:**
- Column-wise access (use CSC)
- Adding new entries (immutable after construction)
- Random access (slow search)

### Matrix-Vector Product

CSR makes SpMV extremely efficient:

```cpp
SparseMatrixCSR<double> A(coo);
std::vector<double> x(n, 1.0);  // Input vector
std::vector<double> y(n);       // Output vector

// y = A*x (SpMV)
A.multiply(x, y);

// y = A^T * x (transpose)
A.multiplyTranspose(x, y);

// y = alpha*A*x + beta*y (BLAS-style gemv)
double alpha = 1.0, beta = 0.0;
A.gemv(alpha, x, beta, y);
```

**Why CSR is Fast for SpMV:**
- Sequential memory access (cache-friendly!)
- No branch misprediction
- Vectorization-friendly
- Minimal pointer chasing

### Accessors and Properties

```cpp
// Dimensions
int m = A.rows();
int n = A.cols();
int nnz = A.nnz();

// Element access (slow, for testing only!)
double a_ij = A(i, j);  // Returns 0 if (i,j) not stored

// Diagonal extraction (fast)
double a_ii = A.diagonal(i);

// Row info
int row_nnz = A.rowNnz(i);  // Non-zeros in row i

// Matrix properties
bool square = A.isSquare();
bool symmetric = A.isSymmetric();  // Checks values
bool sym_pattern = A.hasSymmetricPattern();  // Checks structure only
```

### Norms and Analysis

```cpp
// Matrix norms
double frob = A.normFrobenius();  // ||A||_F = sqrt(sum(a_ij²))
double inf_norm = A.normInf();    // ||A||_∞ = max row sum
double one_norm = A.norm1();      // ||A||_1 = max column sum

// Use for condition number estimation
```

---

## CSC Format (Compressed Sparse Column)

### Concept

CSC is the **column-oriented** version of CSR. Stores non-zeros **column-by-column**.

```
Example Matrix:        CSC Storage:
┌           ┐          values:       [ 2  -1  |  5  | -1  3  1 ]
│ 2  0 -1 │          row_indices:  [ 0   2  |  1  |  0  1  2 ]
│ 0  5  3 │    →     col_pointers: [ 0        2      3          6 ]
│-1  0  1 │                          ^col 0   ^col 1 ^col 2    ^end
└           ┘
```

### Class: `SparseMatrixCSC<T>`

```cpp
#include "pde/sparse/SparseMatrixCSC.h"

// Method 1: From COO
SparseMatrixCSC<double> A(coo);

// Method 2: From CSR (transpose structure!)
SparseMatrixCSR<double> A_csr(coo);
SparseMatrixCSC<double> A_csc(A_csr);  // A_csc stores A^T structure
```

### When to Use CSC

✅ **Ideal for:**
- Column-wise operations
- Matrix transpose (CSC of A = CSR of A^T)
- Some direct solvers (LU, Cholesky)
- Column extraction

❌ **Less common in PDE solvers** (most use CSR)

### Operations

```cpp
// Column access (fast)
int col_nnz = A.colNnz(j);

// Matrix-vector product (also supported)
A.multiply(x, y);  // y = A*x
```

---

## Format Comparison

### Memory Layout

| Format | Storage            | Extras       | Total Memory      |
|--------|-------------------|--------------|-------------------|
| COO    | 3 arrays (r,c,v)  | None         | 2·nnz·sizeof(int) + nnz·sizeof(T) |
| CSR    | values, col_idx   | row_ptrs (m+1) | (m+1)·sizeof(int) + nnz·(sizeof(int)+sizeof(T)) |
| CSC    | values, row_idx   | col_ptrs (n+1) | (n+1)·sizeof(int) + nnz·(sizeof(int)+sizeof(T)) |

**Memory winner:** CSR/CSC (slightly more efficient)

### Operation Performance

| Operation            | COO       | CSR       | CSC       |
|---------------------|-----------|-----------|-----------|
| Construction        | ⭐⭐⭐      | ⚠️        | ⚠️        |
| Add entry           | O(1)      | ❌        | ❌        |
| SpMV (Ax)          | O(nnz)    | ⭐⭐⭐ O(nnz) | O(nnz)    |
| Row access          | O(nnz)    | ⭐⭐⭐ O(1) | O(nnz)    |
| Column access       | O(nnz)    | O(nnz)    | ⭐⭐⭐ O(1) |
| Element (i,j)       | O(nnz)    | O(log nnz)| O(log nnz)|
| Conversion to CSR   | -         | -         | O(nnz)    |
| Conversion to CSC   | -         | O(nnz)    | -         |

### When to Use Each Format

**Use COO when:**
- Building matrix from scratch
- Adding entries in arbitrary order
- Assembling finite difference stencils
- Don't know sparsity pattern upfront

**Use CSR when:**
- Solving linear systems with iterative methods
- Need fast SpMV (matrix-vector products)
- Using CG, BiCGSTAB, GMRES solvers
- Row-wise algorithms

**Use CSC when:**
- Need column-wise access
- Implementing certain direct solvers
- Working with matrix transposes

### Recommended Workflow

```cpp
// 1. Build in COO (flexible)
SparseMatrixCOO<double> builder(n, n);
for (/* assembly loop */) {
    builder.addEntry(i, j, value);
}

// 2. Convert to CSR for computation (fast)
SparseMatrixCSR<double> A(builder);

// 3. Use in iterative solver
ConjugateGradient<double> cg;
cg.solve(A, b, x);
```

---

## Construction Patterns

### Pattern 1: 1D Laplacian (Tridiagonal)

Assembles `-u''(x) = f(x)` on [0,1] with Dirichlet BCs.

```cpp
#include "pde/sparse/SparseMatrixCOO.h"

SparseMatrixCOO<double> buildLaplacian1D(int n, double h) {
    SparseMatrixCOO<double> A(n, n);
    double h_sq_inv = 1.0 / (h * h);
    
    for (int i = 0; i < n; ++i) {
        if (i > 0)
            A.addEntry(i, i-1, -h_sq_inv);  // Lower diagonal
        A.addEntry(i, i, 2.0 * h_sq_inv);   // Main diagonal
        if (i < n-1)
            A.addEntry(i, i+1, -h_sq_inv);  // Upper diagonal
    }
    
    return A;
}

// Usage:
int n = 100;
double h = 1.0 / (n + 1);
auto A_coo = buildLaplacian1D(n, h);
SparseMatrixCSR<double> A(A_coo);
```

### Pattern 2: 2D Laplacian (5-Point Stencil)

Assembles `-∇²u = f` on unit square with Dirichlet BCs.

```cpp
SparseMatrixCOO<double> buildLaplacian2D(int nx, int ny, double hx, double hy) {
    int n = nx * ny;  // Total unknowns
    SparseMatrixCOO<double> A(n, n);
    
    double hx_sq_inv = 1.0 / (hx * hx);
    double hy_sq_inv = 1.0 / (hy * hy);
    double center_coeff = 2.0 * (hx_sq_inv + hy_sq_inv);
    
    auto index = [nx](int i, int j) { return i + j * nx; };
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int row = index(i, j);
            
            // 5-point stencil:
            //        (-hy²)
            //          |
            //  (-hx²)--[ c ]--(- hx²)
            //          |
            //        (-hy²)
            //  where c = 2(1/hx² + 1/hy²)
            
            if (i > 0)      A.addEntry(row, index(i-1, j), -hx_sq_inv);  // West
            if (i < nx-1)   A.addEntry(row, index(i+1, j), -hx_sq_inv);  // East
            if (j > 0)      A.addEntry(row, index(i, j-1), -hy_sq_inv);  // South
            if (j < ny-1)   A.addEntry(row, index(i, j+1), -hy_sq_inv);  // North
            
            A.addEntry(row, row, center_coeff);  // Center
        }
    }
    
    return A;
}

// Usage:
int nx = 50, ny = 50;
double hx = 1.0 / (nx + 1);
double hy = 1.0 / (ny + 1);
auto A_coo = buildLaplacian2D(nx, ny, hx, hy);
SparseMatrixCSR<double> A(A_coo);
```

### Pattern 3: Using Stencil Helpers

```cpp
#include "pde/sparse/SparseMatrixCOO.h"

// 2D grid: nx × ny interior points
int nx = 50, ny = 50;
int n = nx * ny;
double h = 1.0 / (nx + 1);
double h_sq_inv = 1.0 / (h * h);

SparseMatrixCOO<double> A(n, n);

for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
        int row = i + j * nx;
        
        // Automatically adds 5-point stencil!
        A.addLaplacian2D_5pt(row, i, j, nx, h_sq_inv);
    }
}

SparseMatrixCSR<double> A_csr(A);
```

---

## Common Operations

### Diagonal Extraction

```cpp
SparseMatrixCSR<double> A(coo);

// Extract full diagonal
std::vector<double> diag(n);
for (int i = 0; i < n; ++i) {
    diag[i] = A.diagonal(i);
}
```

### Row Scaling

```cpp
// Scale row i by factor alpha
void scaleRow(SparseMatrixCSR<double>& A, int row, double alpha) {
    // Warning: This modifies internal structure!
    // Better to do during construction or use immutable operations
}
```

### Matrix-Vector Operations

```cpp
SparseMatrixCSR<double> A(coo);
std::vector<double> x(n, 1.0);
std::vector<double> y(n);

// Standard SpMV: y = A*x
A.multiply(x, y);

// Transpose: y = A^T*x
A.multiplyTranspose(x, y);

// Scaled: y = 2.0*A*x + 0.5*y
A.gemv(2.0, x, 0.5, y);
```

### Symmetry Checking

```cpp
// Check if A == A^T (structure only)
bool sym_structure = A.hasSymmetricPattern();

// Check if A == A^T (values too, with tolerance)
bool sym_values = A.isSymmetric(1e-10);
```

---

## PDE Assembly Example

Complete example: 2D Poisson equation `-∇²u = f` on [0,1]².

```cpp
#include "pde/sparse/SparseMatrixCOO.h"
#include "pde/sparse/SparseMatrixCSR.h"
#include "pde/solvers/ConjugateGradient.h"
#include <vector>
#include <cmath>

using namespace MML::PDE;

int main() {
    // Grid parameters
    int nx = 50;     // Points in x-direction
    int ny = 50;     // Points in y-direction
    int n = nx * ny; // Total unknowns
    
    double hx = 1.0 / (nx + 1);
    double hy = 1.0 / (ny + 1);
    
    // 1. Build matrix in COO format
    Sparse::SparseMatrixCOO<double> A_coo(n, n);
    A_coo.reserve(5 * n);  // Estimate: 5 non-zeros per row
    
    double hx_sq_inv = 1.0 / (hx * hx);
    double hy_sq_inv = 1.0 / (hy * hy);
    double center = 2.0 * (hx_sq_inv + hy_sq_inv);
    
    auto index = [nx](int i, int j) { return i + j * nx; };
    
    // Assemble 5-point Laplacian stencil
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int row = index(i, j);
            
            if (i > 0)    A_coo.addEntry(row, index(i-1, j), -hx_sq_inv);
            if (i < nx-1) A_coo.addEntry(row, index(i+1, j), -hx_sq_inv);
            if (j > 0)    A_coo.addEntry(row, index(i, j-1), -hy_sq_inv);
            if (j < ny-1) A_coo.addEntry(row, index(i, j+1), -hy_sq_inv);
            
            A_coo.addEntry(row, row, center);
        }
    }
    
    // 2. Convert to CSR for efficient computation
    Sparse::SparseMatrixCSR<double> A(A_coo);
    
    std::cout << "Matrix assembled: " << n << "×" << n 
              << ", nnz = " << A.nnz() << "\n";
    std::cout << "Sparsity: " << (100.0 * A.nnz() / (n * n)) << "%\n";
    
    // 3. Setup right-hand side: f(x,y) = sin(πx)sin(πy)
    std::vector<double> b(n);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double x = (i + 1) * hx;
            double y = (j + 1) * hy;
            b[index(i, j)] = std::sin(M_PI * x) * std::sin(M_PI * y);
        }
    }
    
    // 4. Solve with Conjugate Gradient
    Solvers::ConjugateGradient<double> cg;
    cg.setTolerance(1e-8);
    cg.setMaxIterations(1000);
    
    std::vector<double> x(n, 0.0);  // Initial guess
    
    auto result = cg.solve(A, b, x);
    
    std::cout << "CG converged in " << result.iterations << " iterations\n";
    std::cout << "Residual: " << result.residual << "\n";
    
    return 0;
}
```

**Output:**
```
Matrix assembled: 2500×2500, nnz = 12300
Sparsity: 0.197%
CG converged in 127 iterations
Residual: 9.8e-09
```

---

## Best Practices

### ✅ DO:

1. **Build in COO, compute in CSR**
   ```cpp
   SparseMatrixCOO<double> builder(n, n);
   // ... assemble ...
   SparseMatrixCSR<double> A(builder);  // Convert once
   ```

2. **Reserve space when possible**
   ```cpp
   int estimated_nnz = 5 * n;  // 5-point stencil
   A_coo.reserve(estimated_nnz);
   ```

3. **Use stencil helpers for common patterns**
   ```cpp
   A.addLaplacian2D_5pt(row, i, j, nx, h_sq_inv);
   ```

4. **Check matrix properties before solving**
   ```cpp
   if (!A.isSymmetric()) {
       // Use BiCGSTAB instead of CG
   }
   ```

### ❌ DON'T:

1. **Don't use COO for SpMV in tight loops**
   ```cpp
   // BAD: Slow!
   for (int iter = 0; iter < 1000; ++iter) {
       coo.multiply(x, y);  // O(nnz log nnz) per iteration
   }
   
   // GOOD: Convert once
   SparseMatrixCSR<double> csr(coo);
   for (int iter = 0; iter < 1000; ++iter) {
       csr.multiply(x, y);  // O(nnz) per iteration
   }
   ```

2. **Don't access individual elements frequently**
   ```cpp
   // BAD: O(log nnz) per access
   for (int i = 0; i < n; ++i) {
       for (int j = 0; j < n; ++j) {
           double a_ij = A(i, j);  // Slow!
       }
   }
   
   // GOOD: Use SpMV or direct data access
   ```

3. **Don't forget to consolidate COO with duplicates**
   ```cpp
   // If you might have duplicates:
   A_coo.consolidate();  // Sums them
   ```

### Performance Tips

1. **Minimize format conversions**
   - Convert COO → CSR once, not repeatedly
   - Reuse CSR matrix for multiple solves

2. **Choose right format for algorithm**
   - Iterative solvers: CSR
   - Column operations: CSC
   - Assembly: COO

3. **Monitor sparsity**
   ```cpp
   double sparsity = 100.0 * A.nnz() / (n * n);
   if (sparsity > 10.0) {
       // Consider dense methods instead!
   }
   ```

4. **Preallocate vectors**
   ```cpp
   std::vector<double> x(n), y(n);  // Reuse
   A.multiply(x, y);  // No allocation inside
   ```

---

## See Also

- [Grid Infrastructure](Grid_Infrastructure.md) - Structured grids for PDEs
- [Iterative Solvers](Iterative_Solvers.md) - CG, BiCGSTAB, GMRES
- [Elliptic PDEs](Elliptic_PDEs.md) - Poisson/Laplace examples
- [API Reference](API_Reference.md) - Complete class documentation

---

**Next:** [Grid Infrastructure →](Grid_Infrastructure.md)
