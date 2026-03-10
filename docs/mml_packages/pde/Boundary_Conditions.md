# Boundary Conditions for PDEs

## Overview

**Boundary conditions** (BCs) specify the behavior of the solution at the domain boundaries and are **essential** for well-posed PDE problems. Without proper BCs, a PDE may have infinitely many solutions (or none at all).

This guide covers:
- **Four BC types:** Dirichlet, Neumann, Robin, Periodic
- **How to specify BCs** in 1D, 2D, and 3D
- **How BCs modify the discrete system** (matrix and RHS)
- **Complete working examples** for all BC types

**Key Classes:**
- `BoundaryCondition1D<T>`, `BoundaryCondition2D<T>`, `BoundaryCondition3D<T>` - Individual BC specification
- `BoundaryConditions1D<T>`, `BoundaryConditions2D<T>`, `BoundaryConditions3D<T>` - BC containers
- **Factory functions:** `homogeneousDirichlet1D<T>()`, `homogeneousNeumann1D<T>()`, `periodic1D<T>()`

## Boundary Condition Types

### 1. Dirichlet Boundary Conditions

**Definition:** Fixed value on the boundary

**Mathematical form:** `u = g(x)` on boundary

**Physical interpretation:**
- **Temperature:** Surface held at fixed temperature
- **Displacement:** Fixed support in structural mechanics
- **Potential:** Known voltage at electrode

**Discrete implementation:**
- Replace matrix row: `u[i] = g(x_i)` → coefficient matrix row becomes identity
- Set RHS: `b[i] = g(x_i)`

**Example:**
```cpp
// Zero Dirichlet (u = 0 on boundary)
auto bc = BoundaryCondition1D<double>::Dirichlet(0.0);

// Non-zero constant (u = 100 on boundary)
auto bc = BoundaryCondition1D<double>::Dirichlet(100.0);

// Spatially varying (u = sin(x) on boundary)
auto bc = BoundaryCondition1D<double>::Dirichlet([](double x) {
    return std::sin(x);
});
```

### 2. Neumann Boundary Conditions

**Definition:** Fixed flux/gradient normal to the boundary

**Mathematical form:** `∂u/∂n = g(x)` on boundary

**Physical interpretation:**
- **Heat flux:** Insulated boundary (g=0) or specified heat flow
- **Fluid flow:** Prescribed velocity gradient
- **Diffusion:** Known mass flux

**Discrete implementation:**
- Use **one-sided finite difference** to approximate derivative
- Left boundary: `(u[1] - u[0])/dx = g` → `-u[0] + u[1] = g·dx`
- Right boundary: `(u[n-1] - u[n-2])/dx = g` → `-u[n-2] + u[n-1] = g·dx`

**Example:**
```cpp
// Zero flux (insulated boundary)
auto bc = BoundaryCondition1D<double>::Neumann(0.0);

// Constant influx
auto bc = BoundaryCondition1D<double>::Neumann(5.0);

// Spatially varying flux
auto bc = BoundaryCondition2D<double>::Neumann([](double x, double y) {
    return x * y;  // Flux varies with position
});
```

### 3. Robin Boundary Conditions

**Definition:** Linear combination of value and gradient

**Mathematical form:** `α·u + β·(∂u/∂n) = g(x)` on boundary

**Physical interpretation:**
- **Convective cooling:** `h·(u - u_ambient) + k·(∂u/∂n) = 0` (Newton's law of cooling)
- **Radiation:** Stefan-Boltzmann boundary conditions
- **Absorbing boundaries:** Prevent wave reflection

**Discrete implementation:**
- Combine value and one-sided difference
- Left: `α·u[0] + β·(u[1] - u[0])/dx = g`
  - Rearranged: `(α - β/dx)·u[0] + (β/dx)·u[1] = g`

**Example:**
```cpp
// Convective BC: h*(u - u_ambient) = 0 at boundary
// Rewritten: h*u + k*(du/dn) = h*u_ambient
double h = 10.0;           // Convection coefficient
double k = 1.0;            // Thermal conductivity
double u_ambient = 20.0;   // Ambient temperature

auto bc = BoundaryCondition1D<double>::Robin(h, k, [=](double x) {
    return h * u_ambient;
});
```

### 4. Periodic Boundary Conditions

**Definition:** Solution wraps around (left boundary = right boundary)

**Mathematical form:** `u(x_left) = u(x_right)` and `∂u/∂x(x_left) = ∂u/∂x(x_right)`

**Physical interpretation:**
- **Circular domains:** Heat flow on a ring
- **Wave propagation:** Toroidal geometry
- **Repeating patterns:** Crystallography, Fourier analysis

**Discrete implementation:**
- Constraint: `u[0] = u[n-1]`
- Stencils at boundaries wrap around to opposite side

**Example:**
```cpp
auto bc = BoundaryCondition1D<double>::Periodic();
```

## 1D Boundary Conditions

### Basic Setup

```cpp
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"

using namespace MML::PDE;

// Create grid
Grid1D<double> grid(0.0, 1.0, 50);  // [0, 1], 50 cells

// Method 1: Individual boundary specification
auto bcLeft = BoundaryCondition1D<double>::Dirichlet(0.0);
auto bcRight = BoundaryCondition1D<double>::Dirichlet(1.0);
BoundaryConditions1D<double> bc(bcLeft, bcRight);

// Method 2: Using setters
BoundaryConditions1D<double> bc;
bc.setLeft(BoundaryCondition1D<double>::Dirichlet(0.0));
bc.setRight(BoundaryCondition1D<double>::Neumann(0.0));

// Method 3: Convenience factory functions
auto bc = homogeneousDirichlet1D<double>();  // u = 0 both ends
auto bc = homogeneousNeumann1D<double>();    // du/dx = 0 both ends
auto bc = periodic1D<double>();              // Periodic wrapping
```

### Complete 1D Example: Mixed Boundary Conditions

Solve: `-d²u/dx² = f(x)` on [0, 1]

With BCs: `u(0) = 0` (Dirichlet left), `du/dx(1) = 0` (Neumann right)

```cpp
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/sparse/SparseMatrixCOO.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // Grid setup
    Grid1D<double> grid(0.0, 1.0, 50);
    int n = grid.numNodes();
    double dx = grid.dx();
    
    // Boundary conditions: Dirichlet left, Neumann right
    BoundaryConditions1D<double> bc(
        BoundaryCondition1D<double>::Dirichlet(0.0),  // u(0) = 0
        BoundaryCondition1D<double>::Neumann(0.0)     // du/dx(1) = 0
    );
    
    // Assemble system: -d²u/dx² = f
    SparseMatrixCOO<double> coo(n, n);
    std::vector<double> b(n, 0.0);
    
    // Source term
    auto f = [](double x) { return std::sin(M_PI * x); };
    
    // Interior nodes: Standard 3-point stencil
    for (int i = 1; i < n - 1; i++) {
        double x = grid.x(i);
        coo.add(i, i-1, -1.0 / (dx * dx));
        coo.add(i, i,    2.0 / (dx * dx));
        coo.add(i, i+1, -1.0 / (dx * dx));
        b[i] = f(x);
    }
    
    // Boundary nodes (initial assembly)
    coo.add(0, 0, 1.0);                    // u[0] (will be replaced by BC)
    coo.add(n-1, n-1, 1.0);                // u[n-1] (will be replaced by BC)
    
    // Convert to CSR for efficient solving
    SparseMatrixCSR<double> A = coo.toCSR();
    
    // Apply boundary conditions (modifies A and b)
    bc.apply(A, b, grid);
    
    // After bc.apply():
    // Row 0: u[0] = 0  (identity row, b[0] = 0)
    // Row n-1: -u[n-2] + u[n-1] = 0  (Neumann condition)
    
    // Solve system (use iterative solver)
    // std::vector<double> u = solver.solve(A, b);
    
    return 0;
}
```

**Key points:**
- Assemble **interior nodes first** with standard stencils
- Add **placeholder rows** for boundary nodes
- Call `bc.apply(A, b, grid)` to **replace** boundary rows
- Matrix A and RHS b are modified **in-place**

## 2D Boundary Conditions

### Four Boundaries Setup

```cpp
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/grid/BoundaryConditions.h"

using namespace MML::PDE;

// Create 2D grid: [0,1] × [0,1]
Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, 50, 50);

// Method 1: Set each boundary individually
BoundaryConditions2D<double> bc;
bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));      // x = 0
bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0));     // x = 1
bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));    // y = 0
bc.setTop(BoundaryCondition2D<double>::Dirichlet(100.0));     // y = 1

// Method 2: Set all at once
BoundaryConditions2D<double> bc;
bc.setAll(BoundaryCondition2D<double>::Dirichlet(0.0));

// Method 3: Factory function
auto bc = homogeneousDirichlet2D<double>();  // u = 0 on all boundaries
```

### Complete 2D Example: Heated Plate

Solve: `-∇²u = 0` (Laplace equation) on [0,1] × [0,1]

With BCs:
- Bottom (y=0): u = 0
- Top (y=1): u = 100
- Left/Right (x=0, x=1): Insulated (du/dx = 0)

```cpp
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/sparse/SparseMatrixCOO.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // Grid: 50×50 cells
    Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, 50, 50);
    int numNodes = grid.numNodesX() * grid.numNodesY();
    double dx = grid.dx();
    double dy = grid.dy();
    
    // Boundary conditions
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));      // Insulated
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));     // Insulated
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));  // Cold
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(100.0));   // Hot
    
    // Assemble system: -∇²u = 0
    SparseMatrixCOO<double> coo(numNodes, numNodes);
    std::vector<double> b(numNodes, 0.0);
    
    // Interior nodes: 5-point Laplacian stencil
    grid.forEachInterior([&](int i, int j, double x, double y) {
        int idx = grid.index(i, j);
        
        double cx = 1.0 / (dx * dx);
        double cy = 1.0 / (dy * dy);
        
        coo.add(idx, grid.index(i-1, j), cx);          // West
        coo.add(idx, grid.index(i+1, j), cx);          // East
        coo.add(idx, grid.index(i, j-1), cy);          // South
        coo.add(idx, grid.index(i, j+1), cy);          // North
        coo.add(idx, idx, -2.0 * (cx + cy));           // Center
        
        b[idx] = 0.0;  // No source term
    });
    
    // Boundary nodes: Add placeholder identity rows
    grid.forEachBoundary([&](int i, int j, double x, double y, BoundarySide side) {
        int idx = grid.index(i, j);
        coo.add(idx, idx, 1.0);
    });
    
    // Convert to CSR
    SparseMatrixCSR<double> A = coo.toCSR();
    
    // Apply boundary conditions (modifies A and b)
    bc.apply(A, b, grid);
    
    // After bc.apply():
    // - Top/bottom boundaries: Identity rows with b = 0 or 100
    // - Left/right boundaries: One-sided difference for du/dx = 0
    
    // Solve system
    // std::vector<double> u = solver.solve(A, b);
    
    return 0;
}
```

### 2D Spatially Varying BCs

```cpp
// Non-uniform temperature on top boundary: u = 100*sin(πx)
bc.setTop(BoundaryCondition2D<double>::Dirichlet([](double x, double y) {
    return 100.0 * std::sin(M_PI * x);
}));

// Heat flux varying with position: q = -k*(∂u/∂n) = 10*x
bc.setRight(BoundaryCondition2D<double>::Neumann([](double x, double y) {
    return 10.0 * x;  // Flux increases with height
}));

// Convective cooling on left boundary
double h = 5.0;
double T_ambient = 20.0;
bc.setLeft(BoundaryCondition2D<double>::Robin(h, 1.0, [=](double x, double y) {
    return h * T_ambient;
}));
```

## 3D Boundary Conditions

### Six Boundaries Setup

```cpp
#include "mml/pde/grid/Grid3D.h"
#include "mml/pde/grid/BoundaryConditions.h"

using namespace MML::PDE;

// 3D grid: [0,1] × [0,1] × [0,1]
Grid3D<double> grid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 30, 30, 30);

// Set all six faces
BoundaryConditions3D<double> bc;
bc.setLeft(BoundaryCondition3D<double>::Dirichlet(0.0));      // x = 0
bc.setRight(BoundaryCondition3D<double>::Dirichlet(0.0));     // x = 1
bc.setBottom(BoundaryCondition3D<double>::Dirichlet(0.0));    // y = 0
bc.setTop(BoundaryCondition3D<double>::Dirichlet(0.0));       // y = 1
bc.setFront(BoundaryCondition3D<double>::Neumann(0.0));       // z = 0
bc.setBack(BoundaryCondition3D<double>::Neumann(0.0));        // z = 1

// Or set all at once
bc.setAll(BoundaryCondition3D<double>::Dirichlet(0.0));

// Factory function
auto bc = homogeneousDirichlet3D<double>();  // u = 0 on all faces
```

### Complete 3D Example: Heat in a Cube

Solve: `-∇²u = f` on [0,1]³

With BCs: All six faces held at u = 0

```cpp
#include "mml/pde/grid/Grid3D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/sparse/SparseMatrixCOO.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"

using namespace MML::PDE;
using namespace MML::PDE::Sparse;

int main() {
    // 3D grid
    Grid3D<double> grid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 30, 30, 30);
    int numNodes = grid.numNodesX() * grid.numNodesY() * grid.numNodesZ();
    double dx = grid.dx(), dy = grid.dy(), dz = grid.dz();
    
    // Boundary conditions: Zero Dirichlet on all faces
    auto bc = homogeneousDirichlet3D<double>();
    
    // Assemble system
    SparseMatrixCOO<double> coo(numNodes, numNodes);
    std::vector<double> b(numNodes, 0.0);
    
    // Source term: Gaussian heat source at center
    auto f = [](double x, double y, double z) {
        double r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
        return 1000.0 * std::exp(-100.0 * r2);
    };
    
    // Interior nodes: 7-point Laplacian stencil
    grid.forEachInterior([&](int i, int j, int k, double x, double y, double z) {
        int idx = grid.index(i, j, k);
        
        double cx = 1.0 / (dx * dx);
        double cy = 1.0 / (dy * dy);
        double cz = 1.0 / (dz * dz);
        
        coo.add(idx, grid.index(i-1, j, k), cx);       // West
        coo.add(idx, grid.index(i+1, j, k), cx);       // East
        coo.add(idx, grid.index(i, j-1, k), cy);       // South
        coo.add(idx, grid.index(i, j+1, k), cy);       // North
        coo.add(idx, grid.index(i, j, k-1), cz);       // Front
        coo.add(idx, grid.index(i, j, k+1), cz);       // Back
        coo.add(idx, idx, -2.0 * (cx + cy + cz));     // Center
        
        b[idx] = f(x, y, z);
    });
    
    // Boundary nodes: Placeholder identity
    grid.forEachBoundary([&](int i, int j, int k, double x, double y, double z, BoundarySide side) {
        int idx = grid.index(i, j, k);
        coo.add(idx, idx, 1.0);
    });
    
    // Convert and apply BCs
    SparseMatrixCSR<double> A = coo.toCSR();
    bc.apply(A, b, grid);
    
    // Solve
    // std::vector<double> u = solver.solve(A, b);
    
    return 0;
}
```

## How Boundary Conditions Modify the System

Understanding what `bc.apply(A, b, grid)` does internally:

### Before BC Application

Assembly creates a system with **placeholder rows** for boundary nodes:

```
For boundary node i:
    A[i, i] = 1.0
    b[i] = 0.0
```

### After Dirichlet BC Application

```
For boundary node i with Dirichlet BC u = g(x_i):
    A: Row i becomes identity (all other entries zeroed)
    A[i, i] = 1.0
    b[i] = g(x_i)

Result: u[i] = g(x_i) is enforced exactly
```

**Example (1D):**
```
Before:        After Dirichlet u[0] = 5:
[1 0 0] [u0]   [1 0 0] [u0]   [5]
[* * *] [u1] = [* * *] [u1] = [*]
[* * *] [u2]   [* * *] [u2]   [*]
         ↓ apply BC ↓
       u[0] = 5
```

### After Neumann BC Application

```
For boundary node i with Neumann BC ∂u/∂n = g:
    A: Row i becomes one-sided difference formula
    
Left boundary (1D):
    du/dx ≈ (u[1] - u[0])/dx = g
    A[0, 0] = -1, A[0, 1] = 1
    b[0] = g * dx
    
Right boundary (1D):
    du/dx ≈ (u[n-1] - u[n-2])/dx = g
    A[n-1, n-2] = -1, A[n-1, n-1] = 1
    b[n-1] = g * dx
```

**Example (1D, right boundary):**
```
Before:              After Neumann du/dx = 2:
[* * 0] [u0]         [* * 0] [u0]   [*]
[* * *] [u1] =  →    [* * *] [u1] = [*]
[0 0 1] [u2]         [0 -1 1] [u2]  [2*dx]
         ↓ apply BC ↓
    -u[1] + u[2] = 2*dx
```

### After Robin BC Application

```
For boundary node i with Robin BC α*u + β*(∂u/∂n) = g:
    Combine identity and one-sided difference
    
Left boundary (1D):
    α*u[0] + β*(u[1] - u[0])/dx = g
    (α - β/dx)*u[0] + (β/dx)*u[1] = g
    
    A[0, 0] = α - β/dx
    A[0, 1] = β/dx
    b[0] = g
```

**Example (1D, convective BC at x=0):**
```
α = 5, β = 1, g = 100, dx = 0.1:
    (5 - 1/0.1)*u[0] + (1/0.1)*u[1] = 100
    -5*u[0] + 10*u[1] = 100

A[0, 0] = -5
A[0, 1] = 10
b[0] = 100
```

## Applying BCs to GridFunctions

Sometimes you need to **apply Dirichlet values directly** to a solution vector (e.g., for initial conditions or visualization):

```cpp
// 1D
GridFunction1D<double> u(grid);
u.fill(0.0);  // Initialize to zero
bc.applyToFunction(u);  // Sets boundary values from Dirichlet BCs

// 2D
GridFunction2D<double> u(grid);
u.fill(0.0);
bc.applyToFunction(u);  // Sets all boundary nodes with Dirichlet values

// 3D
GridFunction3D<double> u(grid);
u.fill(0.0);
bc.applyToFunction(u);
```

**Note:** `applyToFunction()` only applies **Dirichlet** boundary conditions (fixed values). Neumann/Robin BCs are only meaningful in the context of the PDE and are handled by `apply(A, b, grid)`.

## Advanced Examples

### Example 1: Time-Dependent BCs

For time-dependent parabolic PDEs, BCs can change at each time step:

```cpp
// Heat solver with time-varying boundary temperature
BoundaryConditions1D<double> bc;

// Temperature at left boundary oscillates: u(0,t) = 50 + 30*sin(ωt)
double omega = 2.0 * M_PI;  // Frequency

// Update BC value at each time step in solver loop
for (double t = 0; t <= tFinal; t += dt) {
    // Update boundary condition
    bc.setLeft(BoundaryCondition1D<double>::Dirichlet([=](double x) {
        return 50.0 + 30.0 * std::sin(omega * t);
    }));
    
    // Solve time step with updated BCs
    // solver.step(dt, bc);
}
```

**Note:** Some solvers (like `HeatSolver1D`) accept BCs in their constructor and evaluate them at each time step internally. Check the solver documentation.

### Example 2: Mixed BC Types on Different Boundaries

```cpp
// 2D heat equation: Different BC on each side
BoundaryConditions2D<double> bc;

// Left: Fixed temperature u = 100
bc.setLeft(BoundaryCondition2D<double>::Dirichlet(100.0));

// Right: Insulated du/dn = 0
bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));

// Bottom: Convective cooling (Robin)
double h = 10.0;  // Convection coefficient
double T_amb = 20.0;
bc.setBottom(BoundaryCondition2D<double>::Robin(h, 1.0, [=](double x, double y) {
    return h * T_amb;
}));

// Top: Heat flux q = 50 W/m²
bc.setTop(BoundaryCondition2D<double>::Neumann(50.0));
```

### Example 3: Piecewise Boundary Conditions

When you need different BC values on different parts of the same boundary:

```cpp
// Top boundary: u = 100 for x < 0.5, u = 0 for x >= 0.5
bc.setTop(BoundaryCondition2D<double>::Dirichlet([](double x, double y) {
    return (x < 0.5) ? 100.0 : 0.0;
}));

// Or smooth transition
bc.setTop(BoundaryCondition2D<double>::Dirichlet([](double x, double y) {
    return 100.0 * std::exp(-50.0 * (x - 0.5) * (x - 0.5));
}));
```

### Example 4: Symmetric Problems

For problems with symmetry, use **Neumann BCs on symmetry lines**:

```cpp
// Solve only on [0, 1] × [0, 1] (quarter domain)
// Physical domain is [-1, 1] × [-1, 1] with symmetry at x=0 and y=0

BoundaryConditions2D<double> bc;
bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));    // Symmetry line x=0
bc.setBottom(BoundaryCondition2D<double>::Neumann(0.0));  // Symmetry line y=0
bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0)); // Physical boundary x=1
bc.setTop(BoundaryCondition2D<double>::Dirichlet(0.0));   // Physical boundary y=1
```

**Rationale:** At symmetry lines, the solution is **mirror-symmetric**, so the normal derivative is zero: `∂u/∂n = 0`.

## Common Pitfalls and Solutions

### Pitfall 1: Forgetting to Apply BCs

**Problem:**
```cpp
SparseMatrixCSR<double> A = coo.toCSR();
// Forgot to call bc.apply(A, b, grid);
std::vector<double> u = solver.solve(A, b);  // ❌ Wrong! BCs not enforced
```

**Solution:** Always call `bc.apply()` after assembling the matrix:
```cpp
SparseMatrixCSR<double> A = coo.toCSR();
bc.apply(A, b, grid);  // ✅ Correct
std::vector<double> u = solver.solve(A, b);
```

### Pitfall 2: Applying BCs Before Converting to CSR

**Problem:**
```cpp
SparseMatrixCOO<double> coo(n, n);
// ... assembly ...
bc.apply(coo, b, grid);  // ❌ Error! bc.apply() expects CSR, not COO
```

**Solution:** Always convert to CSR first:
```cpp
SparseMatrixCOO<double> coo(n, n);
// ... assembly ...
SparseMatrixCSR<double> A = coo.toCSR();  // Convert first
bc.apply(A, b, grid);  // ✅ Now apply BCs
```

**Reason:** COO format is for construction; BCs need to modify specific rows, which requires CSR's row-based structure.

### Pitfall 3: Wrong Neumann Sign Convention

**Problem:**
```cpp
// Trying to model heat flux OUT of domain (q = -k*du/dn > 0)
// But using positive Neumann value
bc.setRight(BoundaryCondition1D<double>::Neumann(10.0));  // ⚠️ Check sign!
```

**Solution:** Understand the sign convention:
- **Neumann BC:** `∂u/∂n = g` (normal derivative)
- For **outward flux** `q = -k*(∂u/∂n)`:
  - Positive `g` means gradient points **outward** (heat **leaving**)
  - Negative `g` means gradient points **inward** (heat **entering**)

```cpp
// Heat flux OUT: ∂u/∂n > 0
bc.setRight(BoundaryCondition1D<double>::Neumann(10.0));  // ✅

// Heat flux IN: ∂u/∂n < 0
bc.setRight(BoundaryCondition1D<double>::Neumann(-10.0));  // ✅
```

### Pitfall 4: Forgetting Robin Coefficient Signs

**Problem:**
```cpp
// Convective BC: h*(u - u_ambient) = -k*(du/dn)
// Rearranged: h*u + k*(du/dn) = h*u_ambient
// But using wrong sign for alpha or beta
auto bc = BoundaryCondition1D<double>::Robin(-h, k, [=](double x) {  // ❌ Wrong sign!
    return h * u_ambient;
});
```

**Solution:** Standard Robin form is `α*u + β*(∂u/∂n) = g`:
```cpp
// For convective BC: h*u + k*(du/dn) = h*u_ambient
auto bc = BoundaryCondition1D<double>::Robin(h, k, [=](double x) {  // ✅ Both positive
    return h * u_ambient;
});
```

### Pitfall 5: Not Assembling Boundary Nodes First

**Problem:**
```cpp
// Assemble interior nodes only
grid.forEachInterior([&](int i, int j, ...) {
    // Add stencil
});
// Then convert to CSR
SparseMatrixCSR<double> A = coo.toCSR();
bc.apply(A, b, grid);  // ❌ Boundary rows don't exist yet!
```

**Solution:** Add placeholder entries for **all nodes** (including boundaries):
```cpp
// Interior nodes: Standard stencil
grid.forEachInterior([&](int i, int j, ...) {
    // Add stencil
});

// Boundary nodes: Placeholder identity
grid.forEachBoundary([&](int i, int j, ...) {
    int idx = grid.index(i, j);
    coo.add(idx, idx, 1.0);  // ✅ Ensures row exists
});

// Now convert and apply BCs
SparseMatrixCSR<double> A = coo.toCSR();
bc.apply(A, b, grid);  // ✅ Boundary rows exist
```

### Pitfall 6: Modifying Solution After Applying BCs

**Problem:**
```cpp
bc.apply(A, b, grid);
std::vector<double> u = solver.solve(A, b);
// Then modify solution
u[0] = 999.0;  // ⚠️ Violates Dirichlet BC!
```

**Solution:** If you modify the solution, you may need to **reapply BCs**:
```cpp
GridFunction1D<double> u(grid, solverResult);
// Modify solution
u(5) = 123.0;
// Restore boundary values
bc.applyToFunction(u);  // ✅ Re-enforces Dirichlet values
```

## Best Practices

### DO:
✅ **Always call `bc.apply()`** after assembling the matrix and before solving

✅ **Convert COO to CSR first** before applying boundary conditions

✅ **Add placeholder rows** for boundary nodes during assembly

✅ **Use factory functions** for common BC types: `homogeneousDirichlet1D<T>()`

✅ **Check signs** for Neumann and Robin BCs (especially for flux problems)

✅ **Use lambdas** for spatially varying or time-dependent BCs

✅ **Test with simple problems** (e.g., u = 0 everywhere) to verify BC implementation

### DON'T:
❌ **Don't forget** to apply boundary conditions (`bc.apply()`)

❌ **Don't apply BCs to COO matrices** (convert to CSR first)

❌ **Don't mix up** normal derivatives and fluxes (check sign convention)

❌ **Don't modify** solution vector without re-applying Dirichlet BCs

❌ **Don't use Neumann BCs everywhere** in elliptic problems (solution becomes non-unique)

❌ **Don't forget** that `applyToFunction()` only works for Dirichlet BCs

## Summary

**Key Takeaways:**

1. **Four BC types:** Dirichlet (fixed value), Neumann (fixed gradient), Robin (mixed), Periodic (wrapping)

2. **Standard workflow:**
   ```cpp
   // 1. Create BCs
   BoundaryConditions1D<double> bc(left, right);
   
   // 2. Assemble system with placeholders
   SparseMatrixCOO<double> coo(n, n);
   // ... add stencils for interior ...
   // ... add identity for boundaries ...
   
   // 3. Convert to CSR
   SparseMatrixCSR<double> A = coo.toCSR();
   
   // 4. Apply BCs (modifies A and b)
   bc.apply(A, b, grid);
   
   // 5. Solve
   std::vector<double> u = solver.solve(A, b);
   ```

3. **BC modifications:** `bc.apply()` **replaces boundary rows** in the matrix with BC-specific formulas

4. **Dirichlet BCs:** Identity rows, `u[i] = g(x_i)`

5. **Neumann BCs:** One-sided differences, `(u[i+1] - u[i])/dx = g`

6. **Robin BCs:** Combination of value and derivative

7. **applyToFunction():** Sets Dirichlet values directly in solution vector (useful for initialization)

8. **Always** add placeholder boundary rows during assembly, then let `bc.apply()` replace them

**Next Steps:**
- See [Iterative_Solvers.md](Iterative_Solvers.md) for solving the resulting linear systems
- See [Elliptic_PDEs.md](Elliptic_PDEs.md) for complete Poisson/Laplace examples
- See [Parabolic_PDEs.md](Parabolic_PDEs.md) for time-dependent BCs in heat equation

