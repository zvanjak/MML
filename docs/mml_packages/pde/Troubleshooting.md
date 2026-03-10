# Troubleshooting Guide for PDE Solvers

> **Diagnose and fix common problems in MML PDE solvers**

This guide helps you identify and resolve issues when working with PDE solvers. Each section describes symptoms, root causes, and step-by-step solutions.

---

## Table of Contents

1. [Compilation Issues](#compilation-issues)
2. [Solver Not Converging](#solver-not-converging)
3. [Wrong or Inaccurate Results](#wrong-or-inaccurate-results)
4. [Segmentation Faults and Crashes](#segmentation-faults-and-crashes)
5. [Performance Problems](#performance-problems)
6. [Grid and Indexing Errors](#grid-and-indexing-errors)
7. [Boundary Condition Errors](#boundary-condition-errors)
8. [Matrix Assembly Issues](#matrix-assembly-issues)
9. [Numerical Instability](#numerical-instability)
10. [Quick Reference Checklists](#quick-reference-checklists)

---

## Compilation Issues

### Missing Headers

**Symptom:** Compiler error `'Class' was not declared in this scope` or `No such file or directory`

**Common Causes:**
- Wrong include path
- Missing namespace declaration
- Incomplete header chain

**Solution:**

```cpp
// ❌ Wrong - missing path or wrong header
#include "PoissonSolver.h"
#include "Grid2D"

// ✅ Correct - full path from mml root
#include "mml/pde/elliptic/PoissonSolver.h"
#include "mml/pde/grid/Grid2D.h"
```

**Required Headers by Feature:**

| Feature | Required Header |
|---------|----------------|
| 1D/2D/3D Grids | `mml/pde/grid/Grid1D.h`, `Grid2D.h`, `Grid3D.h` |
| Boundary Conditions | `mml/pde/grid/BoundaryConditions.h` |
| Sparse Matrices | `mml/pde/sparse/SparseMatrixCSR.h` |
| CG Solver | `mml/pde/solvers/ConjugateGradient.h` |
| BiCGSTAB | `mml/pde/solvers/BiCGSTAB.h` |
| GMRES | `mml/pde/solvers/GMRES.h` |
| Poisson Solver | `mml/pde/elliptic/PoissonSolver.h` |
| Heat Solver | `mml/pde/parabolic/HeatSolver.h` |

---

### Template Instantiation Errors

**Symptom:** Long error messages mentioning template arguments, `no matching function`, or `could not deduce template argument`

**Common Causes:**
- Type mismatch between `float` and `double`
- Missing template argument
- Incompatible function signatures

**Solution:**

```cpp
// ❌ Wrong - type mismatch
Grid2D<double> grid(0, 1, 0, 1, 50, 50);
SparseMatrixCSR<float> A;  // float vs double!

// ✅ Correct - consistent types
using Real = double;
Grid2D<Real> grid(0, 1, 0, 1, 50, 50);
SparseMatrixCSR<Real> A;
```

💡 **Tip:** Always define a type alias at the top of your file:
```cpp
using Real = double;  // Change this one line to switch precision
```

---

### Namespace Issues

**Symptom:** `'XYZ' is not a member of 'MML'` or `use of undeclared identifier`

**Solution:**

```cpp
// Option 1: Using directive (convenient but pollutes namespace)
using namespace MML::PDE;

// Option 2: Namespace alias (recommended)
namespace PDE = MML::PDE;
PDE::Grid2D<double> grid(...);

// Option 3: Fully qualified names (verbose but explicit)
MML::PDE::Grid2D<double> grid(...);
```

---

## Solver Not Converging

### Problem: CG Solver Fails to Converge

**Symptom:** 
- Solver reaches `maxIterations` without converging
- Residual oscillates or grows instead of decreasing
- `SolverResult.converged == false`

**Diagnostic Steps:**

```cpp
// 1. Enable verbose output to see convergence history
ConjugateGradient<double> cg;
cg.setVerbose(true);
cg.setMaxIterations(1000);
cg.setTolerance(1e-8);

auto result = cg.solve(A, b, x);

// Check result
if (!result.converged) {
    std::cout << "Failed after " << result.iterations << " iterations\n";
    std::cout << "Final residual: " << result.residual << "\n";
}
```

**Root Cause 1: Matrix is Not Symmetric Positive Definite (SPD)**

CG **requires** SPD matrices. If your matrix is not SPD, CG may not converge.

**Check if matrix is SPD:**
```cpp
// Quick test: Is matrix symmetric?
bool symmetric = true;
for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
        if (std::abs(A(i,j) - A(j,i)) > 1e-10) {
            symmetric = false;
            std::cout << "Asymmetry at (" << i << "," << j << "): "
                      << A(i,j) << " vs " << A(j,i) << "\n";
        }
    }
}
```

**Solutions:**
- ✅ For non-symmetric systems, use **BiCGSTAB** or **GMRES**
- ✅ Check boundary condition implementation (common cause of asymmetry)
- ✅ Verify matrix assembly code

---

**Root Cause 2: Poorly Conditioned Matrix**

**Symptom:** Very slow convergence or stagnation

**Solution: Add a Preconditioner**

```cpp
// ❌ Without preconditioner - slow convergence
ConjugateGradient<double> cg;
auto result = cg.solve(A, b, x);  // May take 1000+ iterations

// ✅ With Jacobi preconditioner - much faster
JacobiPreconditioner<double> jacobi(A);
cg.setPreconditioner(&jacobi);
auto result = cg.solve(A, b, x);  // Often 50-100 iterations

// ✅ For very ill-conditioned systems, try SSOR
SSORPreconditioner<double> ssor(A, 1.5);  // omega = 1.5
cg.setPreconditioner(&ssor);
```

**Preconditioner Selection:**

| Preconditioner | When to Use | Cost |
|----------------|-------------|------|
| None | Well-conditioned, small systems | O(0) |
| Jacobi | Quick improvement, easy to implement | O(n) |
| SSOR | Better than Jacobi, still simple | O(n) |
| ILU(0) | Best for difficult problems | O(n) setup |

---

**Root Cause 3: Tolerance Too Tight**

**Solution:**
```cpp
// ❌ Too tight - may never converge due to floating-point limits
cg.setTolerance(1e-16);

// ✅ Reasonable for double precision
cg.setTolerance(1e-10);

// ✅ More relaxed if solution quality is acceptable
cg.setTolerance(1e-6);
```

---

**Root Cause 4: Boundary Conditions Not Applied Correctly**

Incorrect boundary conditions can create a singular or near-singular matrix.

**Checklist:**
- ☐ All Dirichlet boundaries have rows set to identity
- ☐ RHS vector has correct boundary values
- ☐ Neumann BCs properly modify matrix coefficients
- ☐ At least one Dirichlet BC exists (pure Neumann → singular matrix!)

⚠️ **Warning:** A pure Neumann problem (no Dirichlet BCs) has a singular matrix because the solution is only determined up to a constant!

**Fix for Pure Neumann Problems:**
```cpp
// Pin one node to fix the constant
// e.g., set u[0] = 0
A.setRow(0, {{0, 1.0}});  // Identity row
b[0] = 0.0;               // Fixed value
```

---

### Problem: BiCGSTAB Stagnates or Oscillates

**Symptom:** Residual reaches a plateau or oscillates without decreasing

**Solutions:**

1. **Try GMRES instead** - more robust for difficult problems:
```cpp
GMRES<double> gmres;
gmres.setRestartSize(50);  // Increase if still stagnating
auto result = gmres.solve(A, b, x);
```

2. **Use a better preconditioner:**
```cpp
ILU0Preconditioner<double> ilu(A);
bicgstab.setPreconditioner(&ilu);
```

3. **Check for near-zero diagonal elements:**
```cpp
for (int i = 0; i < n; ++i) {
    if (std::abs(A(i, i)) < 1e-12) {
        std::cout << "Warning: Near-zero diagonal at row " << i << "\n";
    }
}
```

---

## Wrong or Inaccurate Results

### Problem: Solution Has Wrong Values

**Diagnostic:** Compare with known analytical solution

```cpp
// Analytical solution for -u'' = sin(πx), u(0)=u(1)=0
// Exact: u(x) = sin(πx) / π²
auto exact = [](double x) { return std::sin(M_PI * x) / (M_PI * M_PI); };

// Compute error
double maxError = 0.0;
for (int i = 0; i < grid.numNodes(); ++i) {
    double x = grid.x(i);
    double error = std::abs(u[i] - exact(x));
    maxError = std::max(maxError, error);
}
std::cout << "Max error: " << maxError << "\n";
```

---

**Root Cause 1: Sign Convention Error**

The discrete Laplacian can be defined as `-∇²` or `∇²`. Mixing conventions leads to wrong signs.

```cpp
// ❌ Wrong sign - solver expects -∇² but you assembled ∇²
// Symptom: Solution has wrong sign or diverges

// ✅ Check your stencil sign convention
// Standard 1D Laplacian stencil for -u'':
// [-1/h², 2/h², -1/h²]  (positive on diagonal)
```

---

**Root Cause 2: Incorrect Source Term Scaling**

```cpp
// ❌ Wrong - forgot h² factor
b[i] = f(x_i);

// ✅ Correct - multiply by h² for standard stencil
b[i] = h * h * f(x_i);
```

---

**Root Cause 3: Boundary Values Applied Incorrectly**

```cpp
// ❌ Wrong - boundary value not set in RHS
A.setRow(0, {{0, 1.0}});  // Row is correct...
// ...but forgot to set b[0]!

// ✅ Correct - set both matrix row AND RHS
A.setRow(0, {{0, 1.0}});
b[0] = boundaryValue;
```

---

**Root Cause 4: Grid Resolution Too Coarse**

Finite difference methods have O(h²) error. If h is too large, results are inaccurate.

**Solution:** Perform grid convergence study

```cpp
for (int n : {10, 20, 40, 80, 160}) {
    Grid1D<double> grid(0, 1, n);
    // ... solve ...
    std::cout << "n=" << n << ", h=" << grid.dx() 
              << ", error=" << maxError << "\n";
}
```

**Expected:** Error should decrease by factor of 4 when h halves (O(h²) convergence).

---

### Problem: Solution Oscillates or Has Spurious Values

**Root Cause:** Stability violation in time-dependent problems

For explicit schemes (forward Euler), must satisfy CFL condition:

**1D Heat Equation:**
```
α · Δt / Δx² ≤ 0.5  (stability limit)
```

**Solution:**
```cpp
double alpha = 1.0;  // Thermal diffusivity
double dx = grid.dx();
double dt_max = 0.5 * dx * dx / alpha;  // Maximum stable time step

// Use smaller time step
double dt = 0.4 * dt_max;  // Safety factor
```

---

## Segmentation Faults and Crashes

### Problem: Segmentation Fault in Grid Access

**Root Cause 1: Confusing `nx()` with `numNodesX()` (⚠️ Most Common!)**

```cpp
// ❌ WRONG - nx() is CELL count, not node count!
Grid2D<double> grid(0, 1, 0, 1, 50, 50);
for (int i = 0; i <= grid.nx(); ++i) {  // Off by one!
    for (int j = 0; j <= grid.ny(); ++j) {
        // grid(i, j) crashes when i=50 or j=50!
    }
}

// ✅ CORRECT - use numNodesX() for loop bounds
for (int i = 0; i < grid.numNodesX(); ++i) {
    for (int j = 0; j < grid.numNodesY(); ++j) {
        // Safe: i ∈ [0, 50], j ∈ [0, 50]
    }
}
```

**Memory Aid:**
- `nx()`, `ny()`, `nz()` → **C**ell count (think **C** for Count)
- `numNodesX()`, `numNodesY()`, `numNodesZ()` → **N**ode count (think **N** for Nodes)
- **Nodes = Cells + 1** (always!)

---

**Root Cause 2: Linear Index Out of Bounds**

```cpp
// ❌ Wrong formula - easy to mix up
int idx = i * ny + j;  // Wrong if grid is row-major with different convention

// ✅ Use grid's built-in linear index
int idx = grid.linearIndex(i, j);
```

---

**Root Cause 3: Uninitialized Vectors**

```cpp
// ❌ Vector not sized - undefined behavior
std::vector<double> u;
u[0] = 1.0;  // Crash!

// ✅ Properly sized vector
std::vector<double> u(grid.numNodes(), 0.0);
u[0] = 1.0;  // Safe
```

---

### Problem: Crash in Sparse Matrix Operations

**Root Cause:** Accessing non-existent entries in CSR format

```cpp
// ❌ CSR matrix - cannot access arbitrary (i,j) efficiently
double val = A(5, 10);  // May crash or return garbage if entry doesn't exist

// ✅ Use iterators for CSR access
for (auto it = A.rowBegin(5); it != A.rowEnd(5); ++it) {
    int col = it->column;
    double val = it->value;
}
```

---

## Performance Problems

### Problem: Solver is Very Slow

**Diagnostic:** Profile where time is spent

```cpp
auto start = std::chrono::high_resolution_clock::now();
auto result = solver.solve(A, b, x);
auto end = std::chrono::high_resolution_clock::now();

std::cout << "Solve time: " 
          << std::chrono::duration<double>(end - start).count() 
          << " seconds\n";
std::cout << "Iterations: " << result.iterations << "\n";
std::cout << "Time per iteration: " 
          << std::chrono::duration<double>(end - start).count() / result.iterations 
          << " seconds\n";
```

---

**Root Cause 1: Using Wrong Sparse Format**

```cpp
// ❌ COO format - slow for matrix-vector products
SparseMatrixCOO<double> A;  // O(nnz) per SpMV, no cache locality

// ✅ CSR format - optimal for SpMV
SparseMatrixCSR<double> A;  // O(nnz) per SpMV, cache-friendly

// Convert if needed
SparseMatrixCSR<double> A_csr = A_coo.toCSR();
```

**Format Selection:**

| Format | Best For | SpMV Speed |
|--------|----------|------------|
| COO | Construction, modification | Slow |
| CSR | Row-wise access, SpMV | Fast |
| CSC | Column-wise access | Fast (col ops) |

---

**Root Cause 2: No Preconditioner**

Without preconditioning, convergence can be O(n) iterations for ill-conditioned systems.

**Rule of Thumb:**
- Small/well-conditioned: No preconditioner needed
- Medium systems: Jacobi preconditioner
- Large/ill-conditioned: SSOR or ILU preconditioner

---

**Root Cause 3: Grid Too Fine for Problem Requirements**

```cpp
// ❌ Overkill for many problems - millions of unknowns
Grid3D<double> grid(0, 1, 0, 1, 0, 1, 200, 200, 200);  // 8M nodes!

// ✅ Start coarse, refine only where needed
Grid3D<double> grid(0, 1, 0, 1, 0, 1, 50, 50, 50);  // 125K nodes

// Or use grid convergence study to find adequate resolution
```

---

**Root Cause 4: Debug Build**

Debug builds disable optimizations and add bounds checking.

```bash
# ❌ Debug build - 10-100x slower
cmake -DCMAKE_BUILD_TYPE=Debug ..

# ✅ Release build for performance
cmake -DCMAKE_BUILD_TYPE=Release ..
```

---

## Grid and Indexing Errors

### The nx() vs numNodesX() Distinction (⚠️ Critical!)

This is the **most common source of errors** in PDE codes!

```
Grid with nx = 4 cells:

Cell view:     |  0  |  1  |  2  |  3  |      ← 4 cells
               0     1     2     3     4
Node view:     *-----*-----*-----*-----*      ← 5 nodes

nx() = 4          (cell count)
numNodesX() = 5   (node count = nx + 1)
```

**Correct Loop Patterns:**

```cpp
// Iterate over all nodes
for (int i = 0; i < grid.numNodesX(); ++i) {
    // Process node i
}

// Iterate over all cells
for (int i = 0; i < grid.nx(); ++i) {
    // Process cell i (covers nodes i and i+1)
}

// Iterate over interior nodes only (exclude boundaries)
for (int i = 1; i < grid.numNodesX() - 1; ++i) {
    // Process interior node i
}
```

---

### Linear Indexing in 2D/3D

**2D Grid Indexing:**
```cpp
// Row-major (i varies slowest): idx = i * numNodesY + j
// ✅ Always use the grid's method:
int idx = grid.linearIndex(i, j);

// Reverse mapping:
auto [i, j] = grid.nodeIndices(idx);
```

**3D Grid Indexing:**
```cpp
// idx = i * numNodesY * numNodesZ + j * numNodesZ + k
int idx = grid.linearIndex(i, j, k);

// Reverse mapping:
auto [i, j, k] = grid.nodeIndices(idx);
```

---

## Boundary Condition Errors

### Problem: Wrong Boundary Values in Solution

**Root Cause 1: Using Wrong API**

```cpp
// ❌ Wrong - this API doesn't exist
double val = bc.getValue(Side::Left);

// ✅ Correct - use the proper accessor
auto& leftBC = bc.get(Side::Left);
double val = leftBC.value();  // For constant Dirichlet

// Or for spatially-varying BCs:
double val = leftBC.value(x, y);  // 2D
```

---

**Root Cause 2: Boundary Sides Confusion**

```cpp
// 1D boundaries
bc.get(Side::Left);   // x = xmin
bc.get(Side::Right);  // x = xmax

// 2D boundaries
bc.get(Side::Left);   // x = xmin (west)
bc.get(Side::Right);  // x = xmax (east)
bc.get(Side::Bottom); // y = ymin (south)
bc.get(Side::Top);    // y = ymax (north)

// 3D boundaries
bc.get(Side::Left);   // x = xmin
bc.get(Side::Right);  // x = xmax
bc.get(Side::Bottom); // y = ymin
bc.get(Side::Top);    // y = ymax
bc.get(Side::Front);  // z = zmin
bc.get(Side::Back);   // z = zmax
```

---

**Root Cause 3: Dirichlet vs Neumann Confusion**

```cpp
// Dirichlet: u = value (fixed solution value)
auto bc = BoundaryCondition1D<double>::Dirichlet(0.0);

// Neumann: du/dn = value (fixed gradient)
auto bc = BoundaryCondition1D<double>::Neumann(0.0);  // Zero flux

// ⚠️ They have VERY different effects!
// Dirichlet → Adds row to fix value
// Neumann → Modifies stencil at boundary
```

---

**Root Cause 4: Pure Neumann Problem (Singular Matrix)**

If all boundaries have Neumann conditions, the solution is only determined up to a constant.

**Symptom:** CG fails with "matrix appears singular"

**Solution:**
```cpp
// Option 1: Pin one node
bc.setLeft(BoundaryCondition1D<double>::Dirichlet(0.0));

// Option 2: Use compatibility condition
// ∫f dx = ∫(∂u/∂n) ds  (integral of source = integral of flux)
```

---

## Matrix Assembly Issues

### Problem: Matrix is Not Symmetric

**Symptom:** CG fails or produces wrong results

**Diagnostic:**
```cpp
// Check symmetry
double maxAsym = 0.0;
for (int i = 0; i < n; ++i) {
    for (auto& [j, aij] : A.row(i)) {
        double aji = A(j, i);
        maxAsym = std::max(maxAsym, std::abs(aij - aji));
    }
}
std::cout << "Max asymmetry: " << maxAsym << "\n";
```

**Common Causes:**
1. Neumann BC applied only on one side
2. Robin BC coefficients don't match
3. Variable coefficient handling error

---

### Problem: Matrix Has Zero on Diagonal

**Symptom:** Division by zero, NaN in solution

**Diagnostic:**
```cpp
for (int i = 0; i < n; ++i) {
    if (std::abs(A(i, i)) < 1e-14) {
        std::cout << "Zero diagonal at row " << i << "\n";
    }
}
```

**Solution:** Check row corresponds to properly applied boundary condition

---

## Numerical Instability

### Problem: Solution Contains NaN or Inf

**Diagnostic:**
```cpp
for (int i = 0; i < u.size(); ++i) {
    if (!std::isfinite(u[i])) {
        std::cout << "Non-finite value at index " << i << ": " << u[i] << "\n";
    }
}
```

**Common Causes:**
1. Division by zero (check diagonal)
2. Overflow in matrix-vector product
3. CFL violation in explicit time-stepping
4. Uninitialized values

---

### Problem: Explicit Scheme Blows Up

**Symptom:** Solution grows exponentially with each time step

**Root Cause:** CFL condition violated

**Solution for 1D Heat Equation:**
```cpp
double alpha = 1.0;  // Diffusivity
double dx = grid.dx();
double dt_stable = 0.5 * dx * dx / alpha;

// Use smaller time step
double dt = 0.4 * dt_stable;  // 80% of limit

// Or switch to implicit scheme (unconditionally stable)
```

**Solution for 2D Heat Equation:**
```cpp
double dt_stable = 0.25 * dx * dx / alpha;  // More restrictive in 2D
```

---

## Quick Reference Checklists

### ✅ Before Running Solver

- [ ] Grid dimensions correct (cells vs nodes)?
- [ ] Loop bounds use `numNodesX()`, not `nx()`?
- [ ] Boundary conditions applied to matrix AND RHS?
- [ ] At least one Dirichlet BC (or solution pinned)?
- [ ] Source term scaled by h² (if needed)?
- [ ] Matrix format is CSR for iterative solver?
- [ ] Preconditioner selected for large problems?

### ✅ When Solver Doesn't Converge

1. [ ] Is matrix SPD? If not, use BiCGSTAB/GMRES
2. [ ] Try adding Jacobi preconditioner
3. [ ] Increase max iterations
4. [ ] Relax tolerance (try 1e-6 instead of 1e-10)
5. [ ] Check boundary conditions are not creating singular matrix
6. [ ] Enable verbose output to see convergence history

### ✅ When Results Are Wrong

1. [ ] Compare with analytical solution
2. [ ] Check sign convention (-∇² vs ∇²)
3. [ ] Verify source term scaling
4. [ ] Run grid convergence study
5. [ ] Check boundary values in solution
6. [ ] Verify Dirichlet rows are identity with correct RHS

### ✅ When Code Crashes

1. [ ] Check loop bounds (numNodesX, not nx)
2. [ ] Verify vector sizes before access
3. [ ] Use grid.linearIndex() for 2D/3D indexing
4. [ ] Check sparse matrix entry existence before access
5. [ ] Initialize all vectors before use

### ✅ When Code Is Slow

1. [ ] Build in Release mode
2. [ ] Use CSR format for matrices
3. [ ] Add preconditioner
4. [ ] Check if grid is finer than needed
5. [ ] Profile to find bottleneck

---

## Getting Help

If you've tried the solutions above and still have issues:

1. **Create minimal reproducible example** - smallest code that shows the problem
2. **Include diagnostic output** - solver iterations, residual, error values
3. **Describe expected vs actual behavior** - what should happen vs what happens
4. **Check the examples** - see [Examples_Gallery.md](Examples_Gallery.md) for working code

---

## Related Documentation

- [Quick_Start_Guide.md](Quick_Start_Guide.md) - Getting started tutorial
- [Grid_Infrastructure.md](Grid_Infrastructure.md) - Grid indexing details
- [Boundary_Conditions.md](Boundary_Conditions.md) - BC types and API
- [Iterative_Solvers.md](Iterative_Solvers.md) - Solver configuration
- [Elliptic_PDEs.md](Elliptic_PDEs.md) - Poisson/Laplace equations
- [Parabolic_PDEs.md](Parabolic_PDEs.md) - Heat equation and stability
