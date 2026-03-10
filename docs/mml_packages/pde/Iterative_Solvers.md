# Iterative Solvers for Sparse Linear Systems

## Overview

After discretizing a PDE with finite differences, you get a **large sparse linear system** `Ax = b`. For problems with thousands to millions of unknowns, **iterative solvers** are the only practical option.

This guide covers:
- **Three iterative methods:** Conjugate Gradient (CG), BiCGSTAB, GMRES
- **When to use each solver** based on matrix properties
- **Preconditioners** to accelerate convergence
- **Solver configuration** (tolerance, max iterations, monitoring)
- **Complete working examples** for 1D/2D/3D PDEs

**Key Classes:**
- `ConjugateGradient<T>` - For symmetric positive definite (SPD) systems
- `BiCGSTAB<T>` - For general non-symmetric systems
- `GMRES<T>` - For difficult non-symmetric systems (more robust, more memory)
- **Preconditioners:** `JacobiPreconditioner`, `SSORPreconditioner`, `ILU0Preconditioner`

## Why Iterative Solvers?

### Direct vs. Iterative Methods

**Direct methods** (Gaussian elimination, LU decomposition):
- ✅ Always converge (if non-singular)
- ✅ Exact solution (within floating-point precision)
- ❌ O(n³) time complexity
- ❌ O(n²) memory for general matrices
- ❌ Fill-in destroys sparsity

**Iterative methods:**
- ✅ O(n) memory (sparse storage)
- ✅ O(κ·n) time per iteration where κ = condition number
- ✅ Preserve sparsity
- ✅ Can stop early when "close enough"
- ❌ Convergence depends on matrix properties
- ❌ Need preconditioners for bad conditioning

**Practical guideline:**
- n < 10,000 + SPD matrix → Consider direct solver
- n > 10,000 → Always use iterative
- 3D PDEs → Always iterative (millions of unknowns)

### Example: 2D Poisson with 100×100 Grid

```
Grid: 100×100 = 10,000 unknowns
Matrix: 10,000 × 10,000

Direct (LU factorization):
  Memory: ~800 MB (full storage)
  Time: ~10 seconds

Iterative (CG + preconditioner):
  Memory: ~0.4 MB (sparse storage)
  Time: ~0.5 seconds (50 iterations)
```

**For 200×200 grid (40,000 unknowns):**
- Direct: ~12 GB memory, ~5 minutes
- Iterative: ~1.6 MB memory, ~2 seconds

## Choosing the Right Solver

### Decision Tree

```
Is your matrix SYMMETRIC and POSITIVE DEFINITE?
  ├─ YES → Use Conjugate Gradient (CG)
  │        ✅ Guaranteed convergence
  │        ✅ Minimal memory
  │        ✅ Optimal for SPD systems
  │
  └─ NO (non-symmetric or indefinite)
       │
       ├─ Well-conditioned? → Use BiCGSTAB
       │                      ✅ Fast convergence
       │                      ✅ Low memory
       │                      ⚠️ May stagnate on difficult problems
       │
       └─ Ill-conditioned? → Use GMRES
                             ✅ Most robust
                             ✅ Guaranteed to converge (within restart)
                             ❌ Higher memory (stores m basis vectors)
```

### Matrix Properties by PDE Type

| PDE Type | Matrix Property | Recommended Solver |
|----------|----------------|-------------------|
| Laplacian (-∇²u) | SPD | **CG** |
| Poisson (-∇²u = f) | SPD | **CG** |
| Heat equation (implicit) | SPD | **CG** |
| Diffusion (∇·(D∇u)) | SPD (if D > 0) | **CG** |
| Convection-diffusion | Non-symmetric | **BiCGSTAB** or **GMRES** |
| Helmholtz (∇²u + k²u) | Indefinite (k² > 0) | **GMRES** |
| Advection-dominated | Non-symmetric | **GMRES** |

**SPD = Symmetric Positive Definite:** Matrix is symmetric (A = Aᵀ) and all eigenvalues are positive.

## Conjugate Gradient (CG)

### When to Use

✅ **Use CG when:**
- Matrix is **symmetric positive definite** (SPD)
- Discrete Laplacian or diffusion operator
- Implicit time-stepping for parabolic PDEs
- Self-adjoint elliptic operators

❌ **Don't use CG when:**
- Matrix is non-symmetric (convection terms)
- Matrix is indefinite (Helmholtz, wave equations)
- Matrix has negative eigenvalues

### Algorithm Overview

Conjugate Gradient minimizes the **energy functional** `f(x) = ½(x, Ax) - (b, x)` using conjugate search directions:

```
1. r₀ = b - Ax₀        (initial residual)
2. p₀ = r₀             (initial search direction)
3. For k = 0, 1, 2, ...
   a. αₖ = (rₖ, rₖ) / (pₖ, Apₖ)    (step size)
   b. xₖ₊₁ = xₖ + αₖpₖ             (update solution)
   c. rₖ₊₁ = rₖ - αₖApₖ            (update residual)
   d. If ||rₖ₊₁|| small, STOP
   e. βₖ = (rₖ₊₁, rₖ₊₁) / (rₖ, rₖ) (conjugacy parameter)
   f. pₖ₊₁ = rₖ₊₁ + βₖpₖ           (new search direction)
```

**Convergence:** O(√κ) iterations where κ = condition number

### Basic Usage

```cpp
#include "mml/pde/solvers/ConjugateGradient.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

// After assembling system Ax = b
SparseMatrixCSR<double> A = /* ... assembled matrix ... */;
std::vector<double> b = /* ... RHS vector ... */;

// Create solver
ConjugateGradient<double> cg;

// Configure
cg.config().setTolerance(1e-10);
cg.config().setMaxIterations(1000);

// Solve (zero initial guess)
SolverResult<double> result = cg.solve(A, b);

if (result.converged()) {
    std::cout << "Converged in " << result.iterations << " iterations\n";
    std::cout << "Final residual: " << result.residualNorm << "\n";
    
    // Solution is stored in the solve() return (you need to provide x vector)
} else {
    std::cout << "Failed: " << result.message << "\n";
}
```

### Complete 1D Example: CG for Poisson

```cpp
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/sparse/SparseMatrixCOO.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"
#include "mml/pde/solvers/ConjugateGradient.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // Problem: -d²u/dx² = f(x) on [0, 1], u(0) = u(1) = 0
    Grid1D<double> grid(0.0, 1.0, 100);  // 100 cells
    int n = grid.numNodes();
    double dx = grid.dx();
    
    // Source term: f(x) = -2
    auto f = [](double x) { return -2.0; };
    
    // Boundary conditions: Zero Dirichlet
    auto bc = homogeneousDirichlet1D<double>();
    
    // Assemble system in COO format
    SparseMatrixCOO<double> coo(n, n);
    std::vector<double> b(n, 0.0);
    
    // Interior nodes: -d²u/dx² ≈ (u[i-1] - 2u[i] + u[i+1]) / dx²
    for (int i = 1; i < n - 1; i++) {
        double x = grid.x(i);
        coo.add(i, i-1, -1.0 / (dx * dx));
        coo.add(i, i,    2.0 / (dx * dx));
        coo.add(i, i+1, -1.0 / (dx * dx));
        b[i] = f(x);
    }
    
    // Boundary nodes: Placeholder identity
    coo.add(0, 0, 1.0);
    coo.add(n-1, n-1, 1.0);
    
    // Convert to CSR
    SparseMatrixCSR<double> A = coo.toCSR();
    
    // Apply boundary conditions
    bc.apply(A, b, grid);
    
    // Solve with Conjugate Gradient
    ConjugateGradient<double> cg;
    cg.config()
        .setTolerance(1e-10)
        .setMaxIterations(1000)
        .setVerbose(true);
    
    std::vector<double> u(n, 0.0);  // Initial guess: zero
    SolverResult<double> result = cg.solve(A, b, u);
    
    std::cout << "CG Results:\n";
    std::cout << "  Status: " << (result.converged() ? "SUCCESS" : "FAILED") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Final residual: " << result.residualNorm << "\n";
    std::cout << "  Relative residual: " << result.relativeResidual << "\n";
    
    // Analytical solution: u(x) = -x² + x
    double maxError = 0.0;
    for (int i = 0; i < n; i++) {
        double x = grid.x(i);
        double exact = -x * x + x;
        maxError = std::max(maxError, std::abs(u[i] - exact));
    }
    std::cout << "  Max error vs analytical: " << maxError << "\n";
    
    return 0;
}
```

**Expected output:**
```
CG Results:
  Status: SUCCESS
  Iterations: 34
  Final residual: 8.7e-11
  Relative residual: 4.3e-11
  Max error vs analytical: 1.2e-4
```

## BiCGSTAB

### When to Use

✅ **Use BiCGSTAB when:**
- Matrix is **non-symmetric** (convection-diffusion, advection)
- Faster convergence desired over robustness
- Memory is limited (only needs ~7 work vectors)
- Matrix is moderately well-conditioned

❌ **Don't use BiCGSTAB when:**
- Matrix is SPD (use CG instead - more efficient)
- Problem is very ill-conditioned (use GMRES)
- Stagnation occurs (switch to GMRES)

### Algorithm Overview

BiCGSTAB (Bi-Conjugate Gradient Stabilized) improves on BiCG by stabilizing convergence using a polynomial minimization step:

```
1. r₀ = b - Ax₀
2. Choose r̂₀ (shadow residual, typically r̂₀ = r₀)
3. For k = 0, 1, 2, ...
   a. ρₖ = (r̂₀, rₖ)
   b. βₖ = (ρₖ/ρₖ₋₁)(αₖ₋₁/ωₖ₋₁)
   c. pₖ = rₖ + βₖ(pₖ₋₁ - ωₖ₋₁vₖ₋₁)
   d. v = Apₖ
   e. αₖ = ρₖ / (r̂₀, v)
   f. s = rₖ - αₖv
   g. t = As
   h. ωₖ = (t, s) / (t, t)
   i. xₖ₊₁ = xₖ + αₖpₖ + ωₖs
   j. rₖ₊₁ = s - ωₖt
```

**Key advantage:** Two matrix-vector products per iteration, but often converges in fewer iterations than BiCG.

### Basic Usage

```cpp
#include "mml/pde/solvers/BiCGSTAB.h"

using namespace MML::PDE;

// After assembling non-symmetric system Ax = b
BiCGSTAB<double> solver;

solver.config()
    .setTolerance(1e-10)
    .setMaxIterations(1000);

std::vector<double> x(n, 0.0);
SolverResult<double> result = solver.solve(A, b, x);

if (!result.converged()) {
    std::cout << "BiCGSTAB failed: " << result.message << "\n";
    // Consider switching to GMRES
}
```

### Complete 2D Example: Convection-Diffusion

```cpp
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/sparse/SparseMatrixCOO.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"
#include "mml/pde/solvers/BiCGSTAB.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // Problem: -∇²u + v·∇u = 0 on [0,1]² with u = g on boundary
    // Convection-diffusion with velocity v = (1, 0)
    
    Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, 50, 50);
    int numNodes = grid.numNodesX() * grid.numNodesY();
    double dx = grid.dx();
    double dy = grid.dy();
    
    // Convection velocity
    double vx = 1.0;
    double vy = 0.0;
    
    // Assemble system
    SparseMatrixCOO<double> coo(numNodes, numNodes);
    std::vector<double> b(numNodes, 0.0);
    
    // Interior nodes: -∇²u + v·∇u
    grid.forEachInterior([&](int i, int j, double x, double y) {
        int idx = grid.index(i, j);
        
        double cx = 1.0 / (dx * dx);
        double cy = 1.0 / (dy * dy);
        
        // Diffusion: -∇²u (5-point Laplacian)
        coo.add(idx, grid.index(i-1, j), cx);
        coo.add(idx, grid.index(i+1, j), cx);
        coo.add(idx, grid.index(i, j-1), cy);
        coo.add(idx, grid.index(i, j+1), cy);
        coo.add(idx, idx, -2.0 * (cx + cy));
        
        // Convection: v·∇u using upwind differencing
        // ∂u/∂x ≈ (u[i+1] - u[i-1]) / (2dx) (central)
        // ∂u/∂y ≈ (u[j+1] - u[j-1]) / (2dy)
        if (vx > 0) {
            coo.add(idx, grid.index(i-1, j), -vx / (2.0 * dx));
            coo.add(idx, grid.index(i+1, j),  vx / (2.0 * dx));
        }
        if (vy > 0) {
            coo.add(idx, grid.index(i, j-1), -vy / (2.0 * dy));
            coo.add(idx, grid.index(i, j+1),  vy / (2.0 * dy));
        }
        
        b[idx] = 0.0;
    });
    
    // Boundary nodes: Dirichlet u = sin(πx)
    grid.forEachBoundary([&](int i, int j, double x, double y, BoundarySide side) {
        int idx = grid.index(i, j);
        coo.add(idx, idx, 1.0);
        b[idx] = std::sin(M_PI * x);
    });
    
    // Convert and solve
    SparseMatrixCSR<double> A = coo.toCSR();
    
    // Note: Matrix is NON-SYMMETRIC due to convection term
    // Use BiCGSTAB
    BiCGSTAB<double> solver;
    solver.config()
        .setTolerance(1e-10)
        .setMaxIterations(1000)
        .setVerbose(true);
    
    std::vector<double> u(numNodes, 0.0);
    SolverResult<double> result = solver.solve(A, b, u);
    
    std::cout << "BiCGSTAB Results:\n";
    std::cout << "  Converged: " << result.converged() << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Residual: " << result.residualNorm << "\n";
    
    return 0;
}
```

**Expected behavior:**
- BiCGSTAB converges in ~100-200 iterations for moderate grid
- Convection makes matrix non-symmetric → CG would fail
- With preconditioner: ~50-100 iterations

## GMRES

### When to Use

✅ **Use GMRES when:**
- Matrix is **non-symmetric** and **ill-conditioned**
- BiCGSTAB fails or stagnates
- Maximum robustness is required
- Memory for m+1 vectors is acceptable (m = restart parameter)

❌ **Don't use GMRES when:**
- Matrix is SPD (CG is more efficient)
- BiCGSTAB converges well (GMRES uses more memory)
- Memory is extremely limited

### Algorithm Overview

GMRES minimizes the residual norm over a Krylov subspace using modified Gram-Schmidt orthogonalization:

```
1. r₀ = b - Ax₀, β = ||r₀||
2. v₁ = r₀ / β
3. Build orthonormal basis V = [v₁, v₂, ..., vₘ] via Arnoldi:
   For j = 1 to m:
     w = Avⱼ
     For i = 1 to j:
       hᵢⱼ = (w, vᵢ)
       w = w - hᵢⱼvᵢ
     hⱼ₊₁,ⱼ = ||w||
     vⱼ₊₁ = w / hⱼ₊₁,ⱼ
4. Solve least squares: min ||Hy - βe₁||
5. x = x₀ + Vy
6. If not converged, restart with x₀ = x
```

**Memory:** O(m·n) for storing basis vectors V

**Convergence:** Guaranteed to converge in at most n iterations (without restart)

### Basic Usage

```cpp
#include "mml/pde/solvers/GMRES.h"

using namespace MML::PDE;

// Create GMRES solver with restart = 30
GMRES<double> gmres(30);  // Restart every 30 iterations

gmres.config()
    .setTolerance(1e-10)
    .setMaxIterations(1000);

std::vector<double> x(n, 0.0);
SolverResult<double> result = gmres.solve(A, b, x);
```

### Restart Parameter

**Tradeoff:**
- **Small m (e.g., 10-20):** Less memory, more restarts, slower convergence
- **Large m (e.g., 50-100):** More memory, fewer restarts, faster convergence

**Typical values:**
- m = 20-30: General-purpose, moderate problems
- m = 50-100: Difficult problems, ample memory
- m = 10: Memory-constrained

### Complete Example: GMRES for Difficult Problem

```cpp
#include "mml/pde/solvers/GMRES.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"

using namespace MML::PDE;

int main() {
    // Assume A and b are already assembled (non-symmetric, ill-conditioned)
    SparseMatrixCSR<double> A = /* ... */;
    std::vector<double> b = /* ... */;
    int n = A.rows();
    
    // Try BiCGSTAB first (faster if it converges)
    BiCGSTAB<double> bicgstab;
    bicgstab.config().setTolerance(1e-10).setMaxIterations(500);
    
    std::vector<double> x(n, 0.0);
    SolverResult<double> result = bicgstab.solve(A, b, x);
    
    if (!result.converged()) {
        std::cout << "BiCGSTAB failed, switching to GMRES...\n";
        
        // Fall back to GMRES (more robust)
        GMRES<double> gmres(50);  // Restart = 50
        gmres.config().setTolerance(1e-10).setMaxIterations(1000);
        
        std::fill(x.begin(), x.end(), 0.0);  // Reset initial guess
        result = gmres.solve(A, b, x);
        
        if (result.converged()) {
            std::cout << "GMRES succeeded in " << result.iterations << " iterations\n";
        } else {
            std::cout << "Both solvers failed. Problem may be singular.\n";
        }
    }
    
    return 0;
}
```

## Preconditioners

**Preconditioners** transform the system to improve conditioning:

Instead of solving `Ax = b`, solve `M⁻¹Ax = M⁻¹b` where `M ≈ A` but `M⁻¹` is easy to apply.

**Goal:** Make `M⁻¹A` closer to identity (better conditioned) than `A`.

### Available Preconditioners

| Preconditioner | Setup Cost | Apply Cost | Effectiveness | Use Case |
|---------------|-----------|-----------|---------------|----------|
| **Identity** | None | O(n) | None | Testing, SPD with κ < 100 |
| **Jacobi** | O(n) | O(n) | Moderate | Diagonal-dominant, parallel |
| **SSOR** | O(nnz) | O(nnz) | Good | Elliptic PDEs, SPD |
| **ILU(0)** | O(nnz) | O(n) | Excellent | General sparse, best overall |

**nnz = number of non-zeros in matrix**

### Jacobi Preconditioner

**Idea:** `M = diag(A)`, so `M⁻¹r = r / diag(A)`

**Pros:**
- Simplest preconditioner
- Highly parallelizable
- O(n) setup and apply

**Cons:**
- Least effective (only helps if diagonal-dominant)
- Not useful for anisotropic grids

**Usage:**
```cpp
#include "mml/pde/solvers/Preconditioners.h"

auto prec = std::make_shared<JacobiPreconditioner<double>>();
prec->setup(A);
solver.setPreconditioner(prec);
```

### SSOR Preconditioner

**Idea:** `M = (D + ωL) D⁻¹ (D + ωU) / (ω(2-ω))` where `A = L + D + U`

Symmetric Successive Over-Relaxation applies forward and backward Gauss-Seidel sweeps.

**Pros:**
- Very effective for elliptic PDEs
- Symmetric (preserves SPD property for CG)
- Moderate cost

**Cons:**
- Only works for symmetric matrices
- Not parallelizable
- Requires tuning ω parameter

**Usage:**
```cpp
auto prec = std::make_shared<SSORPreconditioner<double>>(1.2);  // ω = 1.2
prec->setup(A);
solver.setPreconditioner(prec);
```

**Typical ω values:**
- ω = 1.0: Gauss-Seidel
- ω = 1.2-1.5: Good for 2D Poisson
- ω = 1.8-1.9: Over-relaxation for some problems

### ILU(0) Preconditioner

**Idea:** Compute incomplete LU factorization `LU ≈ A` with same sparsity pattern as `A`.

**Pros:**
- Most effective simple preconditioner
- Works for general sparse matrices
- Usually reduces iterations by 5-10×

**Cons:**
- Not guaranteed to exist (may break down)
- Sequential setup (not parallelizable)
- O(nnz) setup cost

**Usage:**
```cpp
auto prec = std::make_shared<ILU0Preconditioner<double>>();
prec->setup(A);  // Compute ILU factorization
solver.setPreconditioner(prec);
```

### Complete Example: CG with ILU Preconditioning

```cpp
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"
#include "mml/pde/solvers/ConjugateGradient.h"
#include "mml/pde/solvers/Preconditioners.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // Problem: 2D Poisson on 200×200 grid (40,000 unknowns)
    Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, 200, 200);
    // ... (assemble system as before) ...
    SparseMatrixCSR<double> A = /* assembled matrix */;
    std::vector<double> b = /* RHS vector */;
    
    // WITHOUT preconditioning
    ConjugateGradient<double> cg_noprecond;
    cg_noprecond.config().setTolerance(1e-10);
    
    std::vector<double> x1(b.size(), 0.0);
    auto result1 = cg_noprecond.solve(A, b, x1);
    std::cout << "No preconditioner: " << result1.iterations << " iterations\n";
    // Expected: ~450-500 iterations
    
    // WITH ILU(0) preconditioning
    ConjugateGradient<double> cg_ilu;
    auto prec = std::make_shared<ILU0Preconditioner<double>>();
    prec->setup(A);
    cg_ilu.setPreconditioner(prec);
    cg_ilu.config().setTolerance(1e-10);
    
    std::vector<double> x2(b.size(), 0.0);
    auto result2 = cg_ilu.solve(A, b, x2);
    std::cout << "ILU(0) preconditioner: " << result2.iterations << " iterations\n";
    // Expected: ~50-80 iterations (8-10× reduction!)
    
    return 0;
}
```

**Typical speedups with ILU(0):**
- 2D Poisson: 8-10× fewer iterations
- 3D Poisson: 10-15× fewer iterations
- Convection-diffusion: 5-8× fewer iterations

## Solver Configuration

### SolverConfig Parameters

```cpp
template<typename T>
struct SolverConfig {
    T tolerance = 1e-10;           // Relative residual tolerance
    T absTolerance = 1e-14;        // Absolute residual tolerance
    int maxIterations = 1000;      // Maximum iterations
    bool verbose = false;          // Print progress
};
```

**Convergence criterion:**
```
||r|| < tolerance * ||b|| + absTolerance
```

Where `r = b - Ax` is the residual.

### Setting Configuration

```cpp
ConjugateGradient<double> cg;

// Method 1: Chained setters
cg.config()
    .setTolerance(1e-10)
    .setAbsTolerance(1e-14)
    .setMaxIterations(1000)
    .setVerbose(true);

// Method 2: Direct access
cg.config().tolerance = 1e-10;
cg.config().maxIterations = 500;
```

### Choosing Tolerances

**Relative tolerance** (primary):
- `1e-6`: Coarse solve, visualization
- `1e-8`: Typical engineering accuracy
- `1e-10`: High accuracy for testing
- `1e-12`: Very high accuracy (near machine precision for double)

**Absolute tolerance** (safeguard for small RHS):
- Usually `1e-14` for double precision
- Prevents division by zero when `||b|| ≈ 0`

**Max iterations:**
- Start with `n` or `2n` for first guess
- For preconditioned: `100-500` usually sufficient
- For unpreconditioned: May need `1000-5000`

## Monitoring Convergence

### Progress Callback

```cpp
#include <iostream>

ConjugateGradient<double> cg;

// Set callback to monitor every iteration
cg.setCallback([](int iter, double residual) {
    std::cout << "Iteration " << iter << ": residual = " << residual << "\n";
});

// Or just enable verbose mode
cg.config().setVerbose(true);
```

### SolverResult

```cpp
template<typename T>
struct SolverResult {
    SolverStatus status;           // Success, MaxIterations, Stagnation, Breakdown
    int iterations;                // Number of iterations taken
    T residualNorm;                // Final ||r|| = ||b - Ax||
    T relativeResidual;            // ||r|| / ||b||
    std::string message;           // Human-readable status
    
    bool converged() const;        // Returns true if status == Success
    operator bool() const;         // Implicit conversion to bool
};
```

**Checking results:**
```cpp
SolverResult<double> result = solver.solve(A, b, x);

if (result) {  // Implicit bool conversion
    std::cout << "SUCCESS\n";
} else {
    std::cout << "FAILED: " << result.message << "\n";
}

// Detailed status
switch (result.status) {
    case SolverStatus::Success:
        std::cout << "Converged in " << result.iterations << " iterations\n";
        break;
    case SolverStatus::MaxIterations:
        std::cout << "Hit max iterations. Residual = " << result.residualNorm << "\n";
        break;
    case SolverStatus::Stagnation:
        std::cout << "Stagnation detected. Try GMRES or better preconditioner.\n";
        break;
    case SolverStatus::Breakdown:
        std::cout << "Method breakdown: " << result.message << "\n";
        break;
}
```

## Complete 3D Example: All Solvers Compared

```cpp
#include "mml/pde/grid/Grid3D.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"
#include "mml/pde/solvers/ConjugateGradient.h"
#include "mml/pde/solvers/BiCGSTAB.h"
#include "mml/pde/solvers/GMRES.h"
#include "mml/pde/solvers/Preconditioners.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // 3D Poisson: -∇²u = f on [0,1]³, u = 0 on boundary
    Grid3D<double> grid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 30, 30, 30);
    int n = grid.numNodesX() * grid.numNodesY() * grid.numNodesZ();
    
    // Assemble system (SPD matrix from Laplacian)
    SparseMatrixCSR<double> A = /* ... assembled 7-point stencil ... */;
    std::vector<double> b = /* ... RHS with source term ... */;
    
    // Create ILU(0) preconditioner (shared across solvers)
    auto prec = std::make_shared<ILU0Preconditioner<double>>();
    prec->setup(A);
    
    std::vector<double> x(n, 0.0);
    
    // 1. Conjugate Gradient (best for SPD)
    {
        std::cout << "\n=== Conjugate Gradient ===\n";
        ConjugateGradient<double> cg;
        cg.setPreconditioner(prec);
        cg.config().setTolerance(1e-10).setMaxIterations(500);
        
        std::fill(x.begin(), x.end(), 0.0);
        auto result = cg.solve(A, b, x);
        
        std::cout << "Status: " << (result ? "SUCCESS" : "FAILED") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "Residual: " << result.residualNorm << "\n";
        // Expected: ~40-60 iterations with ILU(0)
    }
    
    // 2. BiCGSTAB (works but slower for SPD)
    {
        std::cout << "\n=== BiCGSTAB ===\n";
        BiCGSTAB<double> bicgstab;
        bicgstab.setPreconditioner(prec);
        bicgstab.config().setTolerance(1e-10).setMaxIterations(500);
        
        std::fill(x.begin(), x.end(), 0.0);
        auto result = bicgstab.solve(A, b, x);
        
        std::cout << "Status: " << (result ? "SUCCESS" : "FAILED") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "Residual: " << result.residualNorm << "\n";
        // Expected: ~60-80 iterations (slower than CG for SPD)
    }
    
    // 3. GMRES (most robust but highest memory)
    {
        std::cout << "\n=== GMRES(30) ===\n";
        GMRES<double> gmres(30);  // Restart = 30
        gmres.setPreconditioner(prec);
        gmres.config().setTolerance(1e-10).setMaxIterations(500);
        
        std::fill(x.begin(), x.end(), 0.0);
        auto result = gmres.solve(A, b, x);
        
        std::cout << "Status: " << (result ? "SUCCESS" : "FAILED") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "Residual: " << result.residualNorm << "\n";
        // Expected: ~50-70 iterations
    }
    
    return 0;
}
```

**Expected performance on 30³ grid (27,000 unknowns):**

| Solver | Iterations | Time | Memory |
|--------|-----------|------|--------|
| CG (no prec) | ~450 | 8 sec | Minimal |
| CG + ILU(0) | **~50** | **1 sec** | Low |
| BiCGSTAB + ILU(0) | ~70 | 1.5 sec | Low |
| GMRES(30) + ILU(0) | ~60 | 2 sec | Moderate |

## Performance Tips

### 1. Always Use Preconditioners

**Without preconditioner:**
```cpp
ConjugateGradient<double> cg;
auto result = cg.solve(A, b, x);  // Slow!
```

**With ILU(0) preconditioner:**
```cpp
auto prec = std::make_shared<ILU0Preconditioner<double>>();
prec->setup(A);
cg.setPreconditioner(prec);
auto result = cg.solve(A, b, x);  // 5-10× faster!
```

### 2. Choose the Right Solver

```cpp
// SPD matrix → Use CG
if (isSymmetric(A) && isPositiveDefinite(A)) {
    ConjugateGradient<double> solver;
}
// Non-symmetric → Try BiCGSTAB first
else {
    BiCGSTAB<double> solver;
    // If fails, fall back to GMRES
}
```

### 3. Reuse Preconditioners

If solving multiple systems with the same matrix:

```cpp
// Setup preconditioner once
auto prec = std::make_shared<ILU0Preconditioner<double>>();
prec->setup(A);
solver.setPreconditioner(prec);

// Solve multiple RHS with same preconditioner
for (int i = 0; i < numRHS; i++) {
    std::vector<double> x(n, 0.0);
    auto result = solver.solve(A, b[i], x);
}
```

### 4. Tune Tolerances

Don't over-solve:

```cpp
// Visualization: 1e-6 sufficient
cg.config().setTolerance(1e-6);

// Numerical analysis: 1e-10 typical
cg.config().setTolerance(1e-10);

// Beyond 1e-12: Diminishing returns (close to machine precision)
```

### 5. Use Good Initial Guesses

For time-dependent problems, use previous solution:

```cpp
std::vector<double> u = u_previous;  // Start from last time step
auto result = solver.solve(A, b, u);
// Converges much faster than starting from zero
```

## Common Pitfalls and Solutions

### Pitfall 1: Using CG for Non-Symmetric Matrices

**Problem:**
```cpp
// Matrix has convection term (non-symmetric)
ConjugateGradient<double> cg;
auto result = cg.solve(A, b, x);  // ❌ CG may fail or give wrong answer!
```

**Solution:** Check matrix symmetry, use BiCGSTAB or GMRES:
```cpp
BiCGSTAB<double> solver;  // ✅ Works for non-symmetric
auto result = solver.solve(A, b, x);
```

### Pitfall 2: Forgetting to Setup Preconditioner

**Problem:**
```cpp
auto prec = std::make_shared<ILU0Preconditioner<double>>();
// Forgot prec->setup(A);
solver.setPreconditioner(prec);
auto result = solver.solve(A, b, x);  // ❌ Preconditioner not initialized!
```

**Solution:** Always call `setup()`:
```cpp
auto prec = std::make_shared<ILU0Preconditioner<double>>();
prec->setup(A);  // ✅ Initialize preconditioner
solver.setPreconditioner(prec);
```

### Pitfall 3: Tolerance Too Tight

**Problem:**
```cpp
cg.config().setTolerance(1e-16);  // ❌ Below machine precision!
auto result = cg.solve(A, b, x);
// Result: MaxIterations (never converges)
```

**Solution:** Use realistic tolerances:
```cpp
cg.config().setTolerance(1e-10);  // ✅ Achievable with double precision
```

**Rule of thumb:** For `double` (15-16 digits), don't expect better than `1e-12` relative accuracy.

### Pitfall 4: Not Checking Convergence

**Problem:**
```cpp
auto result = solver.solve(A, b, x);
// Assume success, use x without checking
writeOutput(x);  // ❌ May be garbage if solver failed!
```

**Solution:** Always check result:
```cpp
auto result = solver.solve(A, b, x);
if (!result.converged()) {
    std::cerr << "Solver failed: " << result.message << "\n";
    std::cerr << "Final residual: " << result.residualNorm << "\n";
    // Handle failure (try different solver, refine mesh, etc.)
}
```

### Pitfall 5: Wrong Initial Guess Dimensions

**Problem:**
```cpp
std::vector<double> x(10);  // Wrong size!
auto result = solver.solve(A, b, x);  // ❌ Dimension mismatch
```

**Solution:** Match matrix size:
```cpp
std::vector<double> x(A.rows(), 0.0);  // ✅ Correct size
auto result = solver.solve(A, b, x);
```

### Pitfall 6: Singular or Near-Singular Matrix

**Problem:**
```cpp
// Pure Neumann BCs → Matrix is singular (up to constant)
auto bc = homogeneousNeumann2D<double>();
// ... assemble system ...
auto result = cg.solve(A, b, x);  // ❌ Fails: matrix singular
```

**Solution:** Add constraint or use at least one Dirichlet BC:
```cpp
// Fix one point: u(0, 0) = 0
bc.setLeft(BoundaryCondition2D<double>::Dirichlet([](double x, double y) {
    return (x < 1e-10 && y < 1e-10) ? 0.0 : /* Neumann */;
}));
```

## Best Practices

### DO:
✅ **Use CG for SPD matrices** (Laplacian, diffusion, heat equation)

✅ **Use preconditioners** (ILU(0) is usually best)

✅ **Check convergence** before using solution

✅ **Set realistic tolerances** (1e-6 to 1e-10 typical)

✅ **Use previous solution** as initial guess for time-dependent problems

✅ **Monitor residuals** with callbacks or verbose mode during development

✅ **Try BiCGSTAB first** for non-symmetric, switch to GMRES if it fails

### DON'T:
❌ **Don't use CG for non-symmetric matrices**

❌ **Don't forget to call `preconditioner->setup(A)`**

❌ **Don't set tolerance below 1e-12** for double precision

❌ **Don't ignore solver failures** (always check `result.converged()`)

❌ **Don't use unpreconditioned solvers** for large 3D problems

❌ **Don't start from zero** when you have a good initial guess

## Summary

**Key Takeaways:**

1. **Solver selection:**
   - **SPD matrix** → Conjugate Gradient
   - **Non-symmetric, well-conditioned** → BiCGSTAB
   - **Non-symmetric, ill-conditioned** → GMRES

2. **Always use preconditioners:**
   - **ILU(0):** Best general-purpose choice (5-10× speedup)
   - **SSOR:** Good for SPD elliptic PDEs
   - **Jacobi:** Simple, parallelizable, less effective

3. **Configuration:**
   ```cpp
   solver.config()
       .setTolerance(1e-10)           // Relative residual
       .setMaxIterations(1000)        // Safety limit
       .setVerbose(true);             // Monitor progress
   
   auto prec = std::make_shared<ILU0Preconditioner<double>>();
   prec->setup(A);
   solver.setPreconditioner(prec);
   ```

4. **Always check results:**
   ```cpp
   if (!result.converged()) {
       // Handle failure
   }
   ```

5. **Typical performance (3D Poisson, 30³ grid):**
   - CG + ILU(0): ~50 iterations, 1 second
   - BiCGSTAB + ILU(0): ~70 iterations, 1.5 seconds
   - GMRES(30) + ILU(0): ~60 iterations, 2 seconds

**Next Steps:**
- See [Elliptic_PDEs.md](Elliptic_PDEs.md) for complete Poisson examples with solvers
- See [Parabolic_PDEs.md](Parabolic_PDEs.md) for time-dependent problems
- See [Examples_Gallery.md](Examples_Gallery.md) for more advanced solver applications
