# API Reference - PDE Solvers Module

> **Complete class and function reference for MML PDE solvers**

This document provides a comprehensive API reference for all classes in the PDE solver module, organized by category.

---

## Table of Contents

1. [Quick Reference Tables](#quick-reference-tables)
2. [Sparse Matrices](#sparse-matrices)
3. [Grid Infrastructure](#grid-infrastructure)
4. [Boundary Conditions](#boundary-conditions)
5. [Iterative Solvers](#iterative-solvers)
6. [Preconditioners](#preconditioners)
7. [PDE Solvers](#pde-solvers)
8. [Utility Types](#utility-types)
9. [Function Index](#function-index)

---

## Quick Reference Tables

### Classes by Category

| Category | Classes | Header |
|----------|---------|--------|
| **Sparse Matrices** | `SparseMatrixCOO<T>`, `SparseMatrixCSR<T>`, `SparseMatrixCSC<T>` | `pde/sparse/*.h` |
| **Grids** | `Grid1D<T>`, `Grid2D<T>`, `Grid3D<T>` | `pde/grid/Grid*.h` |
| **Grid Functions** | `GridFunction1D<T>`, `GridFunction2D<T>`, `GridFunction3D<T>` | `pde/grid/GridFunction.h` |
| **Boundary Conditions** | `BoundaryCondition1D<T>`, `BoundaryConditions1D<T>`, etc. | `pde/grid/BoundaryConditions.h` |
| **Solvers** | `ConjugateGradient<T>`, `BiCGSTAB<T>`, `GMRES<T>` | `pde/solvers/*.h` |
| **Preconditioners** | `JacobiPreconditioner<T>`, `SSORPreconditioner<T>`, `ILU0Preconditioner<T>` | `pde/solvers/Preconditioners.h` |
| **Elliptic PDEs** | `PoissonSolver1D<T>`, `PoissonSolver2D<T>`, `PoissonSolver3D<T>` | `pde/elliptic/PoissonSolver.h` |
| **Parabolic PDEs** | `HeatSolver1D<T>`, `HeatSolver2D<T>`, `HeatSolver3D<T>` | `pde/parabolic/HeatSolver.h` |

### Namespace

All PDE classes are in namespace `MML::PDE`:
```cpp
using namespace MML::PDE;
// or
namespace PDE = MML::PDE;
```

Sparse matrix classes are in `MML::PDE::Sparse` but are also available in `MML::PDE`.

---

## Sparse Matrices

### SparseMatrixCOO<T>

**Header:** `mml/pde/sparse/SparseMatrixCOO.h`

Coordinate (COO) format - stores triplets (row, col, value). Best for matrix construction.

#### Constructors

```cpp
SparseMatrixCOO();                      // Empty 0×0 matrix
SparseMatrixCOO(int rows, int cols);    // Empty matrix with dimensions
```

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `addEntry` | `void addEntry(int row, int col, T value)` | Add entry (accumulates if exists) |
| `setEntry` | `void setEntry(int row, int col, T value)` | Set entry (replaces if exists) |
| `consolidate` | `void consolidate()` | Sort and sum duplicate entries |
| `clear` | `void clear()` | Remove all entries |
| `rows` | `int rows() const` | Number of rows |
| `cols` | `int cols() const` | Number of columns |
| `nnz` | `int nnz() const` | Number of non-zeros |
| `triplets` | `const std::vector<Triplet<T>>& triplets() const` | Access raw triplets |

#### Example

```cpp
SparseMatrixCOO<double> coo(100, 100);
for (int i = 0; i < 100; ++i) {
    coo.addEntry(i, i, 2.0);           // Diagonal
    if (i > 0) coo.addEntry(i, i-1, -1.0);
    if (i < 99) coo.addEntry(i, i+1, -1.0);
}
auto csr = SparseMatrixCSR<double>(coo);  // Convert to CSR
```

---

### SparseMatrixCSR<T>

**Header:** `mml/pde/sparse/SparseMatrixCSR.h`

Compressed Sparse Row (CSR) format - optimal for matrix-vector products and iterative solvers.

#### Constructors

```cpp
SparseMatrixCSR();                           // Empty 0×0 matrix
SparseMatrixCSR(int rows, int cols);         // Empty matrix with dimensions
SparseMatrixCSR(SparseMatrixCOO<T>& coo);    // Convert from COO
SparseMatrixCSR(int rows, int cols,          // From raw CSR data
                std::vector<T> values,
                std::vector<int> col_indices,
                std::vector<int> row_pointers);
```

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `rows` | `int rows() const` | Number of rows |
| `cols` | `int cols() const` | Number of columns |
| `nnz` | `int nnz() const` | Number of non-zeros |
| `operator*` | `std::vector<T> operator*(const std::vector<T>& x) const` | Matrix-vector product |
| `matvec` | `void matvec(const std::vector<T>& x, std::vector<T>& y) const` | y = A*x (no allocation) |
| `operator()` | `T operator()(int i, int j) const` | Access element (slow for CSR) |
| `getDiagonal` | `std::vector<T> getDiagonal() const` | Extract diagonal |
| `setRow` | `void setRow(int row, const std::vector<std::pair<int, T>>& entries)` | Replace row entries |
| `values` | `const std::vector<T>& values() const` | Raw value array |
| `colIndices` | `const std::vector<int>& colIndices() const` | Raw column index array |
| `rowPointers` | `const std::vector<int>& rowPointers() const` | Raw row pointer array |
| `fromCOO` | `void fromCOO(SparseMatrixCOO<T>& coo)` | Build from COO matrix |

#### Example

```cpp
SparseMatrixCSR<double> A(coo);  // From COO
std::vector<double> x(A.cols(), 1.0);
std::vector<double> y = A * x;  // Matrix-vector product

// Efficient iteration over row i
for (int k = A.rowPointers()[i]; k < A.rowPointers()[i+1]; ++k) {
    int col = A.colIndices()[k];
    double val = A.values()[k];
}
```

---

### SparseMatrixCSC<T>

**Header:** `mml/pde/sparse/SparseMatrixCSC.h`

Compressed Sparse Column (CSC) format - optimal for column-wise access.

Similar API to CSR but with column-major storage.

---

## Grid Infrastructure

### Grid1D<T>

**Header:** `mml/pde/grid/Grid1D.h`

One-dimensional uniform grid on interval [xmin, xmax].

#### Constructors

```cpp
Grid1D(int n);                              // [0,1] with n cells
Grid1D(T xmin, T xmax, int n);              // [xmin,xmax] with n cells
Grid1D(const Interval<T>& domain, int n);   // From interval with n cells
```

#### Properties

| Method | Signature | Description |
|--------|-----------|-------------|
| `xmin` | `T xmin() const` | Left boundary |
| `xmax` | `T xmax() const` | Right boundary |
| `n` | `int n() const` | Number of cells |
| `numNodes` | `int numNodes() const` | Number of nodes (= n + 1) |
| `numInteriorNodes` | `int numInteriorNodes() const` | Interior nodes (= n - 1) |
| `dx` | `T dx() const` | Grid spacing |

#### Coordinate Access

| Method | Signature | Description |
|--------|-----------|-------------|
| `x` | `T x(int i) const` | x-coordinate of node i |
| `locate` | `void locate(T xp, int& i, T& alpha) const` | Find cell containing point |

#### Boundary Detection

| Method | Signature | Description |
|--------|-----------|-------------|
| `isLeftBoundary` | `bool isLeftBoundary(int i) const` | True if i == 0 |
| `isRightBoundary` | `bool isRightBoundary(int i) const` | True if i == n |
| `isBoundary` | `bool isBoundary(int i) const` | True if left or right boundary |
| `isInterior` | `bool isInterior(int i) const` | True if not boundary |

#### Example

```cpp
Grid1D<double> grid(0.0, 1.0, 100);  // [0,1] with 100 cells, 101 nodes
std::cout << "Nodes: " << grid.numNodes() << "\n";     // 101
std::cout << "Spacing: " << grid.dx() << "\n";         // 0.01

for (int i = 0; i < grid.numNodes(); ++i) {
    double xi = grid.x(i);
    // Process node at xi
}
```

---

### Grid2D<T>

**Header:** `mml/pde/grid/Grid2D.h`

Two-dimensional uniform grid on rectangle [xmin,xmax] × [ymin,ymax].

#### Constructors

```cpp
Grid2D(int n);                                        // [0,1]² with n×n cells
Grid2D(int nx, int ny);                               // [0,1]² with nx×ny cells
Grid2D(T xmin, T xmax, T ymin, T ymax, int nx, int ny); // Custom domain
Grid2D(const Rectangle<T>& rect, int nx, int ny);     // From rectangle
```

#### Properties

| Method | Signature | Description |
|--------|-----------|-------------|
| `xmin`, `xmax`, `ymin`, `ymax` | `T ...() const` | Domain bounds |
| `nx`, `ny` | `int ...() const` | Cell counts |
| `numNodesX`, `numNodesY` | `int ...() const` | Node counts (= n + 1) |
| `numNodes` | `int numNodes() const` | Total nodes |
| `numCells` | `int numCells() const` | Total cells |
| `numInteriorNodes` | `int numInteriorNodes() const` | Interior nodes |
| `dx`, `dy` | `T ...() const` | Grid spacing |
| `h` | `T h() const` | Minimum spacing |

#### Index Conversion

| Method | Signature | Description |
|--------|-----------|-------------|
| `index` | `int index(int i, int j) const` | Multi-index to linear |
| `index` | `void index(int idx, int& i, int& j) const` | Linear to multi-index |
| `indexI` | `int indexI(int idx) const` | Extract i from linear index |
| `indexJ` | `int indexJ(int idx) const` | Extract j from linear index |

#### Coordinate Access

| Method | Signature | Description |
|--------|-----------|-------------|
| `x` | `T x(int i) const` | x-coordinate of column i |
| `y` | `T y(int j) const` | y-coordinate of row j |
| `coords` | `std::pair<T,T> coords(int i, int j) const` | (x,y) coordinates |

#### Boundary Detection

| Method | Signature | Description |
|--------|-----------|-------------|
| `isLeftBoundary` | `bool isLeftBoundary(int i, int j) const` | True if i == 0 |
| `isRightBoundary` | `bool isRightBoundary(int i, int j) const` | True if i == nx |
| `isBottomBoundary` | `bool isBottomBoundary(int i, int j) const` | True if j == 0 |
| `isTopBoundary` | `bool isTopBoundary(int i, int j) const` | True if j == ny |
| `isBoundary` | `bool isBoundary(int i, int j) const` | True if any boundary |
| `isInterior` | `bool isInterior(int i, int j) const` | True if not boundary |

#### Example

```cpp
Grid2D<double> grid(0, 1, 0, 1, 50, 50);  // Unit square, 50×50 cells

// Iterate over all nodes
for (int j = 0; j < grid.numNodesY(); ++j) {
    for (int i = 0; i < grid.numNodesX(); ++i) {
        int idx = grid.index(i, j);  // Linear index
        double xi = grid.x(i);
        double yj = grid.y(j);
    }
}

// Iterate over interior nodes only
for (int j = 1; j < grid.numNodesY() - 1; ++j) {
    for (int i = 1; i < grid.numNodesX() - 1; ++i) {
        // Process interior node
    }
}
```

---

### Grid3D<T>

**Header:** `mml/pde/grid/Grid3D.h`

Three-dimensional uniform grid. Similar API to Grid2D with additional z-dimension methods.

#### Additional Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `nz` | `int nz() const` | Cell count in z |
| `numNodesZ` | `int numNodesZ() const` | Node count in z |
| `dz` | `T dz() const` | Grid spacing in z |
| `z` | `T z(int k) const` | z-coordinate |
| `index` | `int index(int i, int j, int k) const` | 3D to linear |
| `isFrontBoundary` | `bool isFrontBoundary(int i, int j, int k) const` | k == 0 |
| `isBackBoundary` | `bool isBackBoundary(int i, int j, int k) const` | k == nz |

---

### GridFunction1D<T>, GridFunction2D<T>, GridFunction3D<T>

**Header:** `mml/pde/grid/GridFunction.h`

Store scalar values on grid nodes.

#### Constructors

```cpp
GridFunction1D(const Grid1D<T>& grid);              // Zero-initialized
GridFunction1D(const Grid1D<T>& grid, T value);     // Constant value
GridFunction1D(const Grid1D<T>& grid, std::function<T(T)> f);  // From function
```

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `operator()` | `T& operator()(int i)` | Access value at node i |
| `operator()` | `T& operator()(int i, int j)` | Access value (2D) |
| `operator()` | `T& operator()(int i, int j, int k)` | Access value (3D) |
| `data` | `std::vector<T>& data()` | Raw data access |
| `grid` | `const GridND<T>& grid() const` | Associated grid |
| `max` | `T max() const` | Maximum value |
| `min` | `T min() const` | Minimum value |

---

## Boundary Conditions

### BCType Enum

**Header:** `mml/pde/grid/BoundaryConditions.h`

```cpp
enum class BCType {
    Dirichlet,   // u = g(x)           - Fixed value
    Neumann,     // ∂u/∂n = g(x)       - Fixed flux
    Robin,       // αu + β(∂u/∂n) = g  - Mixed
    Periodic     // u(xmin) = u(xmax)  - Wrapping
};
```

---

### BoundaryCondition1D<T>

**Header:** `mml/pde/grid/BoundaryConditions.h`

Single boundary condition specification.

#### Factory Methods

```cpp
// Dirichlet
static BoundaryCondition1D Dirichlet(T value);
static BoundaryCondition1D Dirichlet(std::function<T(T)> g);

// Neumann
static BoundaryCondition1D Neumann(T flux);
static BoundaryCondition1D Neumann(std::function<T(T)> g);

// Robin: αu + β(∂u/∂n) = g
static BoundaryCondition1D Robin(T alpha, T beta, std::function<T(T)> g);

// Periodic
static BoundaryCondition1D Periodic();
```

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `type` | `BCType` | Type of boundary condition |
| `value` | `std::function<T(T)>` | Value or flux function g(x) |
| `alpha` | `T` | Robin coefficient for u |
| `beta` | `T` | Robin coefficient for ∂u/∂n |

---

### BoundaryConditions1D<T>

Container for left and right boundary conditions.

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `left` | `const BoundaryCondition1D<T>& left() const` | Left BC |
| `right` | `const BoundaryCondition1D<T>& right() const` | Right BC |
| `setLeft` | `void setLeft(const BoundaryCondition1D<T>& bc)` | Set left BC |
| `setRight` | `void setRight(const BoundaryCondition1D<T>& bc)` | Set right BC |
| `apply` | `void apply(SparseMatrixCSR<T>& A, std::vector<T>& b, const Grid1D<T>& grid)` | Apply BCs to system |
| `applyToFunction` | `void applyToFunction(GridFunction1D<T>& u)` | Apply Dirichlet values |

#### Factory Functions

```cpp
// All boundaries zero Dirichlet
BoundaryConditions1D<T> homogeneousDirichlet1D();

// All boundaries zero Neumann (insulated)
BoundaryConditions1D<T> homogeneousNeumann1D();

// Periodic boundaries
BoundaryConditions1D<T> periodic1D();
```

#### Example

```cpp
// Zero Dirichlet on both ends
auto bc = homogeneousDirichlet1D<double>();

// Custom: u(0)=0, u'(1)=0
BoundaryConditions1D<double> bc;
bc.setLeft(BoundaryCondition1D<double>::Dirichlet(0.0));
bc.setRight(BoundaryCondition1D<double>::Neumann(0.0));

// Spatially varying
bc.setLeft(BoundaryCondition1D<double>::Dirichlet([](double x) {
    return std::sin(M_PI * x);
}));
```

---

### BoundaryCondition2D<T>, BoundaryConditions2D<T>

Similar to 1D versions with additional sides: `Left`, `Right`, `Bottom`, `Top`.

Value functions take `(T x, T y)` arguments.

---

### BoundaryCondition3D<T>, BoundaryConditions3D<T>

Similar with six sides: `Left`, `Right`, `Bottom`, `Top`, `Front`, `Back`.

Value functions take `(T x, T y, T z)` arguments.

---

## Iterative Solvers

### SolverConfig<T>

**Header:** `mml/pde/solvers/IterativeSolverBase.h`

Configuration for iterative solvers.

#### Fields and Methods

| Field/Method | Type/Signature | Default | Description |
|--------------|----------------|---------|-------------|
| `tolerance` | `T` | 1e-10 | Relative tolerance |
| `absTolerance` | `T` | 1e-14 | Absolute tolerance |
| `maxIterations` | `int` | 1000 | Maximum iterations |
| `verbose` | `bool` | false | Print progress |
| `setTolerance` | `SolverConfig& setTolerance(T tol)` | | Fluent setter |
| `setMaxIterations` | `SolverConfig& setMaxIterations(int n)` | | Fluent setter |
| `setVerbose` | `SolverConfig& setVerbose(bool v)` | | Fluent setter |

**Convergence criterion:** `‖r‖ < tol × ‖b‖ + absTol`

---

### SolverResult<T>

**Header:** `mml/pde/solvers/IterativeSolverBase.h`

Result of iterative solve.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `status` | `SolverStatus` | Success, MaxIterations, Breakdown, etc. |
| `iterations` | `int` | Number of iterations performed |
| `residualNorm` | `T` | Final residual ‖b - Ax‖ |
| `relativeResidual` | `T` | residualNorm / ‖b‖ |
| `message` | `std::string` | Status message |

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `converged` | `bool converged() const` | True if status == Success |
| `operator bool` | `operator bool() const` | Same as converged() |

---

### SolverStatus Enum

```cpp
enum class SolverStatus {
    Success,         // Converged within tolerance
    MaxIterations,   // Did not converge
    Stagnation,      // Progress stalled
    Breakdown,       // Method breakdown
    InvalidInput     // Invalid dimensions
};
```

---

### ConjugateGradient<T>

**Header:** `mml/pde/solvers/ConjugateGradient.h`

Conjugate Gradient solver for symmetric positive definite (SPD) systems.

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `solve` | `SolverResult<T> solve(const SparseMatrixCSR<T>& A, const std::vector<T>& b, std::vector<T>& x)` | Solve Ax=b |
| `solve` | `SolverResult<T> solve(const SparseMatrixCSR<T>& A, const std::vector<T>& b)` | Solve with zero initial guess |
| `setConfig` | `void setConfig(const SolverConfig<T>& cfg)` | Set configuration |
| `setPreconditioner` | `void setPreconditioner(std::shared_ptr<Preconditioner<T>> prec)` | Set preconditioner |
| `setCallback` | `void setCallback(SolverCallback<T> cb)` | Set progress callback |

#### Convenience Function

```cpp
template<typename T>
SolverResult<T> solveCG(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    const SolverConfig<T>& config = SolverConfig<T>());
```

#### Example

```cpp
ConjugateGradient<double> cg;
cg.setConfig(SolverConfig<double>().setTolerance(1e-8).setMaxIterations(500));

auto prec = std::make_shared<JacobiPreconditioner<double>>();
prec->setup(A);
cg.setPreconditioner(prec);

std::vector<double> x(n, 0.0);
auto result = cg.solve(A, b, x);

if (result.converged()) {
    std::cout << "Converged in " << result.iterations << " iterations\n";
}
```

---

### BiCGSTAB<T>

**Header:** `mml/pde/solvers/BiCGSTAB.h`

Biconjugate Gradient Stabilized for non-symmetric systems.

Same interface as `ConjugateGradient<T>`.

**Use when:** Matrix is non-symmetric (e.g., convection-diffusion).

#### Convenience Function

```cpp
template<typename T>
SolverResult<T> solveBiCGSTAB(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    const SolverConfig<T>& config = SolverConfig<T>());
```

---

### GMRES<T>

**Header:** `mml/pde/solvers/GMRES.h`

Generalized Minimal Residual for difficult non-symmetric systems.

Same base interface plus restart parameter.

#### Additional Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `setRestartSize` | `void setRestartSize(int m)` | Set restart parameter (default 30) |

**Use when:** BiCGSTAB stagnates or problem is very ill-conditioned.

#### Convenience Function

```cpp
template<typename T>
SolverResult<T> solveGMRES(
    const SparseMatrixCSR<T>& A,
    const std::vector<T>& b,
    std::vector<T>& x,
    int restart = 30,
    const SolverConfig<T>& config = SolverConfig<T>());
```

---

## Preconditioners

### Preconditioner<T> (Abstract Base)

**Header:** `mml/pde/solvers/Preconditioners.h`

```cpp
template<typename T>
class Preconditioner {
public:
    virtual void apply(const std::vector<T>& r, std::vector<T>& z) const = 0;
    virtual void setup(const SparseMatrixCSR<T>& A) {}
};
```

---

### JacobiPreconditioner<T>

**Header:** `mml/pde/solvers/Preconditioners.h`

Diagonal preconditioner: M = diag(A).

**Complexity:** O(n) setup, O(n) apply  
**Use:** Quick improvement, highly parallelizable

```cpp
auto prec = std::make_shared<JacobiPreconditioner<double>>();
prec->setup(A);
solver.setPreconditioner(prec);
```

---

### SSORPreconditioner<T>

**Header:** `mml/pde/solvers/Preconditioners.h`

Symmetric Successive Over-Relaxation.

**Complexity:** O(n) setup, O(nnz) apply  
**Use:** Good for elliptic PDEs, better than Jacobi

#### Additional Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| constructor | `SSORPreconditioner(T omega = 1.0)` | Create with relaxation parameter |
| `setOmega` | `void setOmega(T omega)` | Set ω ∈ (0, 2) |
| `omega` | `T omega() const` | Get ω |

```cpp
auto prec = std::make_shared<SSORPreconditioner<double>>(1.5);
prec->setup(A);
```

---

### ILU0Preconditioner<T>

**Header:** `mml/pde/solvers/Preconditioners.h`

Incomplete LU factorization with no fill-in.

**Complexity:** O(nnz) setup, O(nnz) apply  
**Use:** Best for general sparse systems

```cpp
auto prec = std::make_shared<ILU0Preconditioner<double>>();
prec->setup(A);  // Computes ILU(0) factorization
```

---

## PDE Solvers

### PoissonSolver1D<T>

**Header:** `mml/pde/elliptic/PoissonSolver.h`

Solves the 1D Poisson equation: -u'' = f(x)

#### Constructor

```cpp
PoissonSolver1D(const Grid1D<T>& grid, const BoundaryConditions1D<T>& bc);
```

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `setSource` | `void setSource(std::function<T(T)> f)` | Set source term f(x) |
| `solve` | `GridFunction1D<T> solve(T tol = 1e-10, int maxIter = 10000)` | Solve and return solution |
| `assemble` | `std::pair<SparseMatrixCSR<T>, std::vector<T>> assemble() const` | Get discrete system |
| `getMatrix` | `SparseMatrixCSR<T> getMatrix() const` | Get matrix only |
| `grid` | `const Grid1D<T>& grid() const` | Access grid |

#### Example

```cpp
Grid1D<double> grid(0, 1, 100);
auto bc = homogeneousDirichlet1D<double>();

PoissonSolver1D<double> solver(grid, bc);
solver.setSource([](double x) { return std::sin(M_PI * x); });

auto u = solver.solve();
std::cout << "Solution max: " << u.max() << "\n";
```

---

### PoissonSolver2D<T>

**Header:** `mml/pde/elliptic/PoissonSolver.h`

Solves the 2D Poisson equation: -∇²u = f(x,y)

Uses the 5-point stencil.

#### Constructor

```cpp
PoissonSolver2D(const Grid2D<T>& grid, const BoundaryConditions2D<T>& bc);
```

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `setSource` | `void setSource(std::function<T(T, T)> f)` | Set source term f(x,y) |
| `solve` | `GridFunction2D<T> solve(T tol = 1e-10, int maxIter = 10000)` | Solve |
| `assemble` | `std::pair<SparseMatrixCSR<T>, std::vector<T>> assemble() const` | Get discrete system |
| `grid` | `const Grid2D<T>& grid() const` | Access grid |

#### Example

```cpp
Grid2D<double> grid(0, 1, 0, 1, 50, 50);
auto bc = homogeneousDirichlet2D<double>();

PoissonSolver2D<double> solver(grid, bc);
solver.setSource([](double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
});

auto u = solver.solve();
```

---

### PoissonSolver3D<T>

**Header:** `mml/pde/elliptic/PoissonSolver.h`

Solves the 3D Poisson equation using the 7-point stencil.

Similar API with 3D types.

---

### HeatSolver1D<T>

**Header:** `mml/pde/parabolic/HeatSolver.h`

Solves the 1D heat equation: ∂u/∂t = α·∂²u/∂x²

#### Constructor

```cpp
HeatSolver1D(const Grid1D<T>& grid, const BoundaryConditions1D<T>& bc, T alpha);
```

#### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `setInitialCondition` | `void setInitialCondition(std::function<T(T)> u0)` | Set u(x,0) = u0(x) |
| `setSource` | `void setSource(std::function<T(T, T)> f)` | Set source f(x,t) |
| `setObserver` | `void setObserver(std::function<void(T, const std::vector<T>&)> obs)` | Callback per step |
| `setTheta` | `void setTheta(T theta)` | θ=0 (FE), 0.5 (CN), 1 (BE) |
| `getTime` | `T getTime() const` | Current simulation time |
| `getSolution` | `const std::vector<T>& getSolution() const` | Current solution |
| `getGrid` | `const Grid1D<T>& getGrid() const` | Access grid |
| `computeCFLLimit` | `T computeCFLLimit() const` | Max stable dt for explicit |
| `isStable` | `bool isStable(T dt) const` | Check if dt is stable |
| `step` | `bool step(T dt, TimeScheme scheme)` | Single time step |
| `solve` | `HeatSolveResult<T> solve(T dt, T finalTime, TimeScheme scheme)` | Integrate to final time |

#### Example

```cpp
Grid1D<double> grid(0, 1, 100);
auto bc = homogeneousDirichlet1D<double>();

HeatSolver1D<double> solver(grid, bc, 0.01);  // α = 0.01
solver.setInitialCondition([](double x) {
    return std::sin(M_PI * x);
});

// Implicit Crank-Nicolson (unconditionally stable)
auto result = solver.solve(0.001, 1.0, TimeScheme::CrankNicolson);

if (result.stable) {
    std::cout << "Final temperature range: [" 
              << result.minTemperature << ", " 
              << result.maxTemperature << "]\n";
}
```

---

### HeatSolver2D<T>, HeatSolver3D<T>

**Header:** `mml/pde/parabolic/HeatSolver.h`

Similar API for 2D and 3D heat equations.

Note: CFL condition becomes more restrictive in higher dimensions:
- 1D: Δt ≤ h²/(2α)
- 2D: Δt ≤ h²/(4α)
- 3D: Δt ≤ h²/(6α)

---

### TimeScheme Enum

```cpp
enum class TimeScheme {
    ForwardEuler,    // Explicit, O(Δt), conditionally stable
    BackwardEuler,   // Implicit, O(Δt), unconditionally stable
    CrankNicolson,   // Implicit, O(Δt²), unconditionally stable (recommended)
    Theta            // General θ-method
};
```

---

### HeatSolveResult<T>

```cpp
template<typename T>
struct HeatSolveResult {
    bool stable;           // Solution remained stable
    int steps;             // Number of time steps
    T finalTime;           // Final simulation time
    T maxTemperature;      // Max u at final time
    T minTemperature;      // Min u at final time
    T totalEnergy;         // ∫u dx (conservation check)
    std::string message;   // Status message
};
```

---

## Utility Types

### Triplet<T>

**Header:** `mml/pde/sparse/SparseMatrixCOO.h`

```cpp
template<typename T>
struct Triplet {
    int row;
    int col;
    T value;
};
```

---

### Interval<T>

**Header:** `mml/pde/grid/Grid.h`

```cpp
template<typename T>
class Interval {
public:
    Interval();                      // [0, 1]
    Interval(T min, T max);
    T min() const;
    T max() const;
    T length() const;
    bool contains(T x) const;
};
```

---

### Rectangle<T>

**Header:** `mml/pde/grid/Grid.h`

```cpp
template<typename T>
class Rectangle {
public:
    Rectangle();                                  // [0,1]×[0,1]
    Rectangle(T xmin, T xmax, T ymin, T ymax);
    T xMin(), xMax(), yMin(), yMax() const;
    T width(), height() const;
    bool contains(T x, T y) const;
};
```

---

### Box<T>

**Header:** `mml/pde/grid/Grid.h`

```cpp
template<typename T>
class Box {
public:
    Box();  // [0,1]³
    Box(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax);
    // Similar accessors to Rectangle plus z dimension
};
```

---

## Function Index

### Factory Functions

| Function | Header | Description |
|----------|--------|-------------|
| `homogeneousDirichlet1D<T>()` | BoundaryConditions.h | Zero Dirichlet BCs |
| `homogeneousDirichlet2D<T>()` | BoundaryConditions.h | Zero Dirichlet BCs |
| `homogeneousDirichlet3D<T>()` | BoundaryConditions.h | Zero Dirichlet BCs |
| `homogeneousNeumann1D<T>()` | BoundaryConditions.h | Zero flux BCs |
| `periodic1D<T>()` | BoundaryConditions.h | Periodic BCs |

### Convenience Solver Functions

| Function | Header | Description |
|----------|--------|-------------|
| `solveCG<T>(A, b, x, config)` | ConjugateGradient.h | CG solve |
| `solveBiCGSTAB<T>(A, b, x, config)` | BiCGSTAB.h | BiCGSTAB solve |
| `solveGMRES<T>(A, b, x, restart, config)` | GMRES.h | GMRES solve |

### Utility Functions (in Sparse namespace)

| Function | Description |
|----------|-------------|
| `dot(x, y)` | Vector dot product |
| `norm2(x)` | Euclidean norm |
| `axpy(a, x, y)` | y += a*x |
| `residual(A, x, b)` | Return b - A*x |

---

## See Also

- [Quick_Start_Guide.md](Quick_Start_Guide.md) - Getting started tutorial
- [Iterative_Solvers.md](Iterative_Solvers.md) - Solver selection guide
- [Boundary_Conditions.md](Boundary_Conditions.md) - BC types and usage
- [Troubleshooting.md](Troubleshooting.md) - Common issues and solutions
