# Elliptic PDEs: Poisson and Laplace Equations

> **Comprehensive Guide to Steady-State Boundary Value Problems**

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Theory](#mathematical-theory)
3. [Poisson Equation](#poisson-equation)
4. [Laplace Equation](#laplace-equation)
5. [Finite Difference Discretization](#finite-difference-discretization)
6. [1D Implementation](#1d-implementation)
7. [2D Implementation](#2d-implementation)
8. [3D Implementation](#3d-implementation)
9. [Verification and Validation](#verification-and-validation)
10. [Complete Examples](#complete-examples)
11. [Advanced Topics](#advanced-topics)

---

## Overview

Elliptic PDEs describe **steady-state** phenomena where the solution depends only on spatial coordinates, not time. These equations arise naturally in:

- **Electrostatics**: Electric potential from charge distribution
- **Heat conduction**: Steady-state temperature distribution
- **Fluid dynamics**: Potential flow, stream functions
- **Structural mechanics**: Membrane deflection, stress analysis
- **Image processing**: Smoothing, inpainting

### Key Characteristics

- **Boundary value problems**: Solution determined by boundary conditions
- **No time evolution**: Steady-state equilibrium
- **Symmetric matrices**: Resulting linear systems are symmetric positive definite (SPD)
- **Efficient solvers**: CG works beautifully for these problems

---

## Mathematical Theory

### The Poisson Equation

The general form in domain Ω with boundary ∂Ω:

```
-∇²u = f   in Ω
    u = g   on ∂Ω  (Dirichlet BC)
```

Where:
- `u(x)` is the unknown solution
- `f(x)` is the source term (forcing function)
- `g(x)` is the boundary condition
- `∇²` is the Laplacian operator

**In Different Dimensions:**

**1D:**
```
-d²u/dx² = f(x),  x ∈ [a, b]
```

**2D:**
```
-∂²u/∂x² - ∂²u/∂y² = f(x, y),  (x,y) ∈ Ω ⊂ ℝ²
```

**3D:**
```
-∂²u/∂x² - ∂²u/∂y² - ∂²u/∂z² = f(x, y, z),  (x,y,z) ∈ Ω ⊂ ℝ³
```

### The Laplace Equation

Special case of Poisson equation with **f = 0**:

```
∇²u = 0   in Ω
   u = g   on ∂Ω
```

**Physical Interpretation:**
- Solution is **harmonic** (no internal sources)
- Represents equilibrium state
- Satisfies **maximum principle**: max/min on boundary
- Solution is infinitely smooth in the interior

---

## Poisson Equation

### Physical Meaning

The Poisson equation models:

1. **Electrostatics**: `-∇²φ = ρ/ε₀`
   - φ: electric potential
   - ρ: charge density
   - ε₀: permittivity

2. **Heat Conduction**: `-k∇²T = q`
   - T: temperature
   - q: heat source density
   - k: thermal conductivity

3. **Gravitational Potential**: `∇²Φ = 4πGρ`
   - Φ: gravitational potential
   - ρ: mass density
   - G: gravitational constant

### Well-Posedness

For Poisson equation with Dirichlet BC:
- **Existence**: Solution exists if f ∈ L²(Ω) and g is continuous on ∂Ω
- **Uniqueness**: Solution is unique
- **Stability**: Solution depends continuously on data (f, g)

### Regularity

If f is smooth and boundary is smooth:
- Solution u is C^∞ in interior
- Regularity near boundary depends on BC and geometry

---

## Laplace Equation

### Harmonic Functions

Solutions to Laplace equation are called **harmonic functions**. They satisfy:

1. **Mean Value Property**: Value at point = average over any sphere centered at that point
2. **Maximum Principle**: Maximum and minimum occur on boundary
3. **Uniqueness**: Dirichlet problem has unique solution
4. **Smoothness**: Infinitely differentiable in interior

### Applications

**1. Electrostatics (No Charge)**
```
∇²φ = 0   (no charges in region)
```

**2. Steady-State Heat (No Sources)**
```
∇²T = 0   (equilibrium, no heat generation)
```

**3. Potential Flow**
```
∇²ψ = 0   (stream function, incompressible irrotational flow)
```

**4. Membrane Deflection**
```
∇²w = 0   (elastic membrane, no external load)
```

---

## Finite Difference Discretization

### 1D: Three-Point Stencil

**Continuous equation:**
```
-u''(x) = f(x),  x ∈ [a, b]
```

**Taylor series at x_i:**
```
u(x_i+1) = u(x_i) + h·u'(x_i) + (h²/2)·u''(x_i) + (h³/6)·u'''(x_i) + O(h⁴)
u(x_i-1) = u(x_i) - h·u'(x_i) + (h²/2)·u''(x_i) - (h³/6)·u'''(x_i) + O(h⁴)
```

**Add equations:**
```
u(x_i+1) + u(x_i-1) = 2u(x_i) + h²·u''(x_i) + O(h⁴)
```

**Solve for u'':**
```
u''(x_i) = (u(x_i-1) - 2u(x_i) + u(x_i+1))/h² + O(h²)
```

**Discrete equation:**
```
-(u_{i-1} - 2u_i + u_{i+1})/h² = f_i
```

**Stencil notation:**
```
[-1  2  -1]/h² · u = f
```

### 2D: Five-Point Stencil

**Continuous equation:**
```
-∂²u/∂x² - ∂²u/∂y² = f(x, y)
```

**Discrete equation (uniform grid h):**
```
-(u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4u_{i,j})/h² = f_{i,j}
```

**Stencil notation:**
```
       [ 0  -1   0 ]
1/h² · [-1   4  -1] · u = f
       [ 0  -1   0 ]
```

**Matrix structure:**
For nx×ny grid, system is (nx-1)×(ny-1) dimensional with:
- **Diagonal**: 4/h²
- **Off-diagonal**: -1/h² (4 non-zeros per row)
- **Bandwidth**: O(nx)
- **Sparsity**: ~5 non-zeros per row

### 3D: Seven-Point Stencil

**Continuous equation:**
```
-∂²u/∂x² - ∂²u/∂y² - ∂²u/∂z² = f(x, y, z)
```

**Discrete equation:**
```
-(u_{i-1,j,k} + u_{i+1,j,k} + u_{i,j-1,k} + u_{i,j+1,k} + 
  u_{i,j,k-1} + u_{i,j,k+1} - 6u_{i,j,k})/h² = f_{i,j,k}
```

**Stencil notation:**
```
Layer k-1:     Layer k:           Layer k+1:
[ 0  0  0 ]    [ 0  -1   0 ]      [ 0  0  0 ]
[ 0 -1  0 ]    [-1   6  -1]       [ 0 -1  0 ]
[ 0  0  0 ]    [ 0  -1   0 ]      [ 0  0  0 ]
```
All scaled by 1/h².

### Truncation Error

All three schemes are **second-order accurate**:
```
Local truncation error: τ = O(h²)
Global error: ||u - u_h|| = O(h²)
```

**Convergence study:** Halving h should reduce error by factor of 4.

---

## 1D Implementation

### Basic Usage

```cpp
#include <MML.h>
using namespace MML;

// Define domain: [0, 1]
Interval<double> domain(0.0, 1.0);

// Create grid with 100 cells
Grid1D<double> grid(domain, 100);

// Homogeneous Dirichlet boundary conditions: u(0) = u(1) = 0
auto bc = homogeneousDirichlet1D<double>();

// Create solver
PoissonSolver1D<double> solver(grid, bc);

// Define source term: f(x) = sin(πx)
auto f = [](double x) { return std::sin(M_PI * x); };

// Set source and solve
solver.setSource(f);
auto solution = solver.solve();

// Access solution at grid points
for (int i = 0; i < grid.numNodesX(); ++i) {
    double x = grid.nodeX(i);
    double u = solution[i];
    std::cout << x << " " << u << "\n";
}
```

### Non-Homogeneous Boundary Conditions

```cpp
// u(0) = 1.0, u(1) = 2.0
BoundaryConditions1D<double> bc;
bc.setLeft(BoundaryCondition1D<double>::Dirichlet(1.0));
bc.setRight(BoundaryCondition1D<double>::Dirichlet(2.0));

PoissonSolver1D<double> solver(grid, bc);
```

### Analytical Solution Verification

**Problem:** `-u'' = sin(πx)`, `u(0) = u(1) = 0`

**Analytical solution:** `u(x) = sin(πx)/π²`

```cpp
#include <MML.h>
#include <iostream>
#include <cmath>

int main() {
    using namespace MML;
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    auto bc = homogeneousDirichlet1D<double>();
    
    PoissonSolver1D<double> solver(grid, bc);
    
    auto f = [](double x) { return std::sin(M_PI * x); };
    solver.setSource(f);
    auto solution = solver.solve();
    
    // Compare with analytical solution
    double maxError = 0.0;
    for (int i = 0; i < grid.numNodesX(); ++i) {
        double x = grid.nodeX(i);
        double numerical = solution[i];
        double analytical = std::sin(M_PI * x) / (M_PI * M_PI);
        double error = std::abs(numerical - analytical);
        maxError = std::max(maxError, error);
    }
    
    std::cout << "Maximum error: " << maxError << "\n";
    // Expected: ~8e-6 for 100 cells (O(h²) accuracy)
    
    return 0;
}
```

---

## 2D Implementation

### Basic Laplace Equation

**Problem:** `∇²u = 0` in unit square, `u = g` on boundary

```cpp
#include <MML.h>
using namespace MML;

// Define domain: [0,1] × [0,1]
Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);

// Create 50×50 grid
Grid2D<double> grid(domain, 50, 50);

// Boundary conditions: u = sin(πx) on top, u = 0 elsewhere
BoundaryConditions2D<double> bc;

bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));
bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0));
bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));
bc.setTop(BoundaryCondition2D<double>::Dirichlet(
    [](double x, double y) { return std::sin(M_PI * x); }
));

// Create solver (Laplace: f = 0)
PoissonSolver2D<double> solver(grid, bc);

// Solve with zero source
auto f = [](double x, double y) { return 0.0; };
solver.setSource(f);
auto solution = solver.solve();

// Access solution
for (int j = 0; j < grid.numNodesY(); ++j) {
    for (int i = 0; i < grid.numNodesX(); ++i) {
        int idx = grid.index(i, j);
        double u = solution[idx];
        // Use u...
    }
}
```

### Poisson with Source Term

**Problem:** `-∇²u = exp(-((x-0.5)² + (y-0.5)²)/0.1²)` (Gaussian source)

```cpp
// Same grid and BC setup...

// Gaussian source centered at (0.5, 0.5)
auto f = [](double x, double y) {
    double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
    return std::exp(-r2 / (0.1 * 0.1));
};

PoissonSolver2D<double> solver(grid, bc);
solver.setSource(f);
auto solution = solver.solve();
```

### Mixed Boundary Conditions

```cpp
// u = 0 on left, right, bottom
// ∂u/∂n = 0 on top (Neumann - insulated)

BoundaryConditions2D<double> bc;
bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));
bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0));
bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));
bc.setTop(BoundaryCondition2D<double>::Neumann(0.0));  // Insulated

PoissonSolver2D<double> solver(grid, bc);
```

### CSV Output for Visualization

```cpp
#include <fstream>

void saveSolution2D(const Grid2D<double>& grid, 
                    const std::vector<double>& solution,
                    const std::string& filename) {
    std::ofstream file(filename);
    file << "x,y,u\n";
    
    for (int j = 0; j < grid.numNodesY(); ++j) {
        for (int i = 0; i < grid.numNodesX(); ++i) {
            double x = grid.nodeX(i);
            double y = grid.nodeY(j);
            double u = solution[grid.index(i, j)];
            file << x << "," << y << "," << u << "\n";
        }
    }
}

// Usage
solver.setSource(f);
auto solution = solver.solve();
saveSolution2D(grid, solution, "laplace_solution.csv");
```

**Python plotting:**
```python
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = pd.read_csv('laplace_solution.csv')
x = data['x'].unique()
y = data['y'].unique()
X, Y = np.meshgrid(x, y)
Z = data['u'].values.reshape(len(y), len(x))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
plt.show()
```

---

## 3D Implementation

### Basic 3D Poisson

```cpp
#include <MML.h>
using namespace MML;

// Define 3D domain: [0,1]³
Box<double> domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

// Create 30×30×30 grid (27,000 unknowns!)
Grid3D<double> grid(domain, 30, 30, 30);

// Homogeneous Dirichlet on all faces
auto bc = homogeneousDirichlet3D<double>();

PoissonSolver3D<double> solver(grid, bc);

// Point source at center
auto f = [](double x, double y, double z) {
    double r2 = (x - 0.5) * (x - 0.5) + 
                (y - 0.5) * (y - 0.5) + 
                (z - 0.5) * (z - 0.5);
    return (r2 < 0.01) ? 100.0 : 0.0;  // Localized source
};

solver.setSource(f);
auto solution = solver.solve();

// Access solution
for (int k = 0; k < grid.numNodesZ(); ++k) {
    for (int j = 0; j < grid.numNodesY(); ++j) {
        for (int i = 0; i < grid.numNodesX(); ++i) {
            int idx = grid.index(i, j, k);
            double u = solution[idx];
            // Use u...
        }
    }
}
```

### 3D Grid Memory Requirements

**Storage:**
- Grid: `nx × ny × nz` nodes
- Solution vector: `(nx+1) × (ny+1) × (nz+1)` doubles
- Sparse matrix: ~7 non-zeros per row

**Examples:**
- 50³ = 125,000 unknowns ≈ 1 MB solution + ~7 MB matrix
- 100³ = 1,000,000 unknowns ≈ 8 MB solution + ~56 MB matrix
- 200³ = 8,000,000 unknowns ≈ 64 MB solution + ~448 MB matrix

**Practical limits on desktop:**
- **Direct solver**: ~100³ (1M unknowns)
- **Iterative solver (CG)**: ~300³ (27M unknowns)

---

## Verification and Validation

### Method of Manufactured Solutions (MMS)

**Strategy:** Choose analytical solution, derive source term, verify implementation.

**1D Example:**

**Choose:** `u(x) = sin(2πx)`

**Compute:** `f(x) = -u''(x) = 4π² sin(2πx)`

**Boundary:** `u(0) = 0`, `u(1) = 0` ✓

```cpp
// Manufactured solution test
auto u_exact = [](double x) { return std::sin(2 * M_PI * x); };
auto f_manufactured = [](double x) { 
    return 4 * M_PI * M_PI * std::sin(2 * M_PI * x); 
};

Grid1D<double> grid(Interval<double>(0.0, 1.0), 100);
auto bc = homogeneousDirichlet1D<double>();
PoissonSolver1D<double> solver(grid, bc);

solver.setSource(f_manufactured);
auto solution = solver.solve();

// Compute L∞ error
double maxError = 0.0;
for (int i = 0; i < grid.numNodesX(); ++i) {
    double x = grid.nodeX(i);
    double error = std::abs(solution[i] - u_exact(x));
    maxError = std::max(maxError, error);
}

std::cout << "Max error: " << maxError << "\n";
```

**2D Example:**

**Choose:** `u(x,y) = sin(πx) sin(πy)`

**Compute:** `f(x,y) = -∇²u = 2π² sin(πx) sin(πy)`

**Boundary:** `u = 0` on all sides ✓

```cpp
auto u_exact = [](double x, double y) { 
    return std::sin(M_PI * x) * std::sin(M_PI * y); 
};

auto f = [](double x, double y) {
    return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
};

Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
Grid2D<double> grid(domain, 50, 50);
auto bc = homogeneousDirichlet2D<double>();
PoissonSolver2D<double> solver(grid, bc);

solver.setSource(f);
auto solution = solver.solve();

// Verify
double maxError = 0.0;
for (int j = 0; j < grid.numNodesY(); ++j) {
    for (int i = 0; i < grid.numNodesX(); ++i) {
        double x = grid.nodeX(i);
        double y = grid.nodeY(j);
        double error = std::abs(solution[grid.index(i, j)] - u_exact(x, y));
        maxError = std::max(maxError, error);
    }
}

std::cout << "Max error: " << maxError << "\n";
```

### Grid Convergence Study

**Test O(h²) convergence:**

```cpp
#include <vector>

std::vector<int> gridSizes = {10, 20, 40, 80, 160};
std::vector<double> errors;
std::vector<double> h_values;

for (int n : gridSizes) {
    Grid2D<double> grid(Rectangle<double>(0,1,0,1), n, n);
    auto bc = homogeneousDirichlet2D<double>();
    PoissonSolver2D<double> solver(grid, bc);
    
    solver.setSource(f_manufactured);
    auto solution = solver.solve();
    
    double maxError = 0.0;
    // ... compute maxError ...
    
    errors.push_back(maxError);
    h_values.push_back(1.0 / n);
}

// Check convergence rate
for (size_t i = 1; i < errors.size(); ++i) {
    double rate = std::log(errors[i-1] / errors[i]) / 
                  std::log(h_values[i-1] / h_values[i]);
    std::cout << "Rate: " << rate << " (expect ~2.0)\n";
}
```

**Expected output:**
```
Grid 10:  error = 1.2e-3
Grid 20:  error = 3.0e-4  (rate = 2.0)
Grid 40:  error = 7.5e-5  (rate = 2.0)
Grid 80:  error = 1.9e-5  (rate = 2.0)
Grid 160: error = 4.7e-6  (rate = 2.0)
```

---

## Complete Examples

### Example 1: Electrostatic Potential

**Problem:** Charged conducting plates in 2D

```cpp
#include <MML.h>
#include <iostream>

int main() {
    using namespace MML;
    
    // Domain: [-1, 1] × [-1, 1]
    Rectangle<double> domain(-1.0, 1.0, -1.0, 1.0);
    Grid2D<double> grid(domain, 100, 100);
    
    // Boundary conditions:
    // - Left plate (x = -1):  φ = +100 V
    // - Right plate (x = +1): φ = -100 V
    // - Top/bottom: φ = 0 (grounded)
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Dirichlet(100.0));
    bc.setRight(BoundaryCondition2D<double>::Dirichlet(-100.0));
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(0.0));
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));
    
    PoissonSolver2D<double> solver(grid, bc);
    
    // No charges in region (Laplace equation)
    auto f = [](double x, double y) { return 0.0; };
    solver.setSource(f);
    
    auto phi = solver.solve();
    
    // Compute electric field: E = -∇φ
    std::cout << "Electric field at center:\n";
    int i_mid = grid.numNodesX() / 2;
    int j_mid = grid.numNodesY() / 2;
    double h = grid.dx();
    
    int idx = grid.index(i_mid, j_mid);
    int idx_left = grid.index(i_mid - 1, j_mid);
    int idx_right = grid.index(i_mid + 1, j_mid);
    
    double Ex = -(phi[idx_right] - phi[idx_left]) / (2 * h);
    std::cout << "Ex = " << Ex << " V/m\n";
    
    return 0;
}
```

### Example 2: Steady-State Heat in Plate

**Problem:** Square plate with different edge temperatures

```cpp
#include <MML.h>
#include <iostream>
#include <fstream>

int main() {
    using namespace MML;
    
    Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 80, 80);
    
    // Temperature boundary conditions:
    // - Left: 0°C (cold)
    // - Right: 0°C (cold)
    // - Bottom: 0°C (cold)
    // - Top: 100°C (hot)
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Dirichlet(0.0));
    bc.setRight(BoundaryCondition2D<double>::Dirichlet(0.0));
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(100.0));
    
    PoissonSolver2D<double> solver(grid, bc);
    
    // No internal heat sources
    auto f = [](double x, double y) { return 0.0; };
    solver.setSource(f);
    
    auto T = solver.solve();
    
    // Save temperature distribution
    std::ofstream file("plate_temperature.csv");
    file << "x,y,T\n";
    for (int j = 0; j < grid.numNodesY(); ++j) {
        for (int i = 0; i < grid.numNodesX(); ++i) {
            file << grid.nodeX(i) << "," 
                 << grid.nodeY(j) << ","
                 << T[grid.index(i, j)] << "\n";
        }
    }
    
    std::cout << "Temperature at center: " 
              << T[grid.index(40, 40)] << "°C\n";
    
    return 0;
}
```

### Example 3: Point Source (Green's Function)

**Problem:** Point charge in 2D domain

```cpp
#include <MML.h>
#include <iostream>

int main() {
    using namespace MML;
    
    Rectangle<double> domain(-1.0, 1.0, -1.0, 1.0);
    Grid2D<double> grid(domain, 100, 100);
    
    // Grounded boundary
    auto bc = homogeneousDirichlet2D<double>();
    
    PoissonSolver2D<double> solver(grid, bc);
    
    // Point source at origin with Gaussian approximation
    double sigma = 0.05;  // Width of Gaussian
    auto f = [sigma](double x, double y) {
        double r2 = x * x + y * y;
        return 100.0 * std::exp(-r2 / (2 * sigma * sigma)) / 
               (2 * M_PI * sigma * sigma);
    };
    
    solver.setSource(f);
    auto phi = solver.solve();
    
    // Check radial symmetry
    std::cout << "Radial profile:\n";
    for (int i = 50; i < grid.numNodesX(); ++i) {
        double x = grid.nodeX(i);
        double r = x;  // Along x-axis from center
        double phi_r = phi[grid.index(i, 50)];
        std::cout << "r = " << r << ", φ = " << phi_r << "\n";
    }
    
    return 0;
}
```

---

## Advanced Topics

### Matrix Properties

The discrete Poisson system `Au = f` has special structure:

**1D:**
```
A = 1/h² · tridiag(-1, 2, -1)
```
- **Dimension**: (n-1) × (n-1) for n+1 grid points
- **Bandwidth**: 3
- **Condition number**: κ(A) ≈ 4/h² = O(n²)

**2D:**
```
A = 1/h² · (I ⊗ T + T ⊗ I)
```
where T = tridiag(-1, 2, -1)
- **Dimension**: (nx-1)(ny-1) × (nx-1)(ny-1)
- **Bandwidth**: O(nx)
- **Sparsity**: 5 non-zeros per row
- **Condition number**: κ(A) ≈ 8/h² = O(n²)

**Properties:**
- **Symmetric**: Aᵀ = A
- **Positive Definite**: xᵀAx > 0 for all x ≠ 0
- **Diagonally Dominant**: |aᵢᵢ| > Σⱼ≠ᵢ |aᵢⱼ|
- **M-matrix**: Inverse has all positive entries

### Solver Choice

**Direct Solvers:**
- LU decomposition: O(n³) time, O(n²) space
- Cholesky (for SPD): O(n³) time, more stable
- Banded solver: O(bn²) for bandwidth b
- **Good for**: Small problems (n < 10,000)

**Iterative Solvers:**
- Conjugate Gradient: O(n^1.5) with preconditioning
- Multigrid: O(n) (optimal!)
- **Good for**: Large problems (n > 10,000)

**MML Implementation:**
- Uses Conjugate Gradient with Jacobi preconditioner
- Typical iterations: 10-50 for 2D, 20-100 for 3D
- Tolerance: 1e-10 (adjustable)

### Boundary Condition Implementation

**Dirichlet BC (u = g):**
- **Strong enforcement**: Remove boundary nodes from system
- Modify RHS for interior nodes adjacent to boundary
- Our implementation uses this approach

**Neumann BC (∂u/∂n = g):**
- **One-sided difference** at boundary
- Example (left boundary): `(-3u₀ + 4u₁ - u₂)/(2h) = g`
- Modifies stencil at boundary

**Robin BC (αu + β∂u/∂n = g):**
- Combination of Dirichlet and Neumann
- Common in heat transfer (convective cooling)

### Performance Optimization

**Grid Size Selection:**
- Balance accuracy vs. cost
- Rule of thumb: 10-20 points per characteristic length
- For smooth solutions: 50×50 often sufficient
- For sharp gradients: 200×200 or adaptive refinement

**Sparse Matrix Format:**
- Use CSR for matrix-vector products
- Build with COO, convert to CSR before solving
- Memory: ~40 bytes per non-zero in CSR

**Preconditioner Impact:**
- No preconditioner: 100-500 iterations
- Jacobi: 50-200 iterations
- SSOR: 30-100 iterations
- ILU: 10-50 iterations

### Extensions

**Non-Uniform Grids:**
- Refine near boundaries or singularities
- Use coordinate transformations
- Requires non-uniform stencils

**Higher-Order Schemes:**
- 9-point stencil in 2D: O(h⁴) accuracy
- Compact schemes: better for smooth solutions
- Trade-off: accuracy vs. complexity

**Adaptive Mesh Refinement:**
- Refine grid where error is large
- Requires error estimators
- Complex but powerful for multi-scale problems

**Nonlinear Poisson:**
- `-∇²u = f(u)` (semilinear)
- Use Newton iteration or fixed-point
- Each iteration solves linear Poisson

---

## Summary

### Key Takeaways

✅ **Elliptic PDEs** describe steady-state equilibrium (no time dependence)

✅ **Poisson equation**: `-∇²u = f` with source term

✅ **Laplace equation**: `∇²u = 0` (special case, no sources)

✅ **Finite differences** convert PDE to sparse linear system

✅ **Matrix is SPD**: Conjugate Gradient is perfect solver

✅ **Second-order accurate**: O(h²) convergence

✅ **Verification**: Use manufactured solutions and grid refinement

### Quick Reference

| Dimension | Stencil | Matrix Size | Non-zeros/row | Bandwidth |
|-----------|---------|-------------|---------------|-----------|
| 1D        | 3-point | n-1         | 3             | 3         |
| 2D        | 5-point | (nx-1)(ny-1)| 5             | O(nx)     |
| 3D        | 7-point | (nx-1)(ny-1)(nz-1)| 7       | O(nx·ny)  |

### Related Documentation

- [Quick Start Guide](Quick_Start_Guide.md) - Get started quickly
- [Sparse Matrices](Sparse_Matrices.md) - Understanding CSR/COO formats
- [Grid Infrastructure](Grid_Infrastructure.md) - Grid1D/2D/3D details
- [Boundary Conditions](Boundary_Conditions.md) - BC types and API
- [Iterative Solvers](Iterative_Solvers.md) - CG, BiCGSTAB, GMRES
- [Examples Gallery](Examples_Gallery.md) - More complete examples
- [Troubleshooting](Troubleshooting.md) - Common issues

---

*Elliptic PDEs documentation - Part of MinimalMathLibrary PDE Solver Module*
