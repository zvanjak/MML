# Quick Start Guide - PDE Solvers in 5 Minutes

Get up and running with MML PDE solvers fast! This guide walks you through solving your first elliptic and parabolic PDEs.

---

## 📋 Prerequisites

```cpp
// Required C++17 headers (all standard library)
#include <iostream>
#include <cmath>
#include <functional>

// MML PDE headers
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/sparse/SparseMatrixCSR.h"
#include "mml/pde/solvers/ConjugateGradient.h"
#include "mml/pde/elliptic/PoissonSolver.h"
#include "mml/pde/parabolic/HeatSolver.h"

using namespace MML::PDE;
using Real = double;  // or float for single precision
```

💡 **Tip:** All MML headers are in the `mml/` directory. The PDE module lives in `mml/pde/`.

---

## ⚡ Example 1: Solve 1D Poisson Equation

**Problem:** Solve $-u''(x) = \sin(\pi x)$ on $[0,1]$ with $u(0) = u(1) = 0$

### Step 1: Create Grid

```cpp
// Create 1D grid on [0,1] with 100 cells (101 nodes)
Interval<Real> domain(0.0, 1.0);
Grid1D<Real> grid(domain, 100);

std::cout << "Grid: " << grid.numNodes() << " nodes, h = " << grid.dx() << "\n";
// Output: Grid: 101 nodes, h = 0.01
```

### Step 2: Set Boundary Conditions

```cpp
// Homogeneous Dirichlet: u(0) = 0, u(1) = 0
auto bc = homogeneousDirichlet1D<Real>();

// Or create custom boundary conditions:
// BoundaryConditions1D<Real> bc;
// bc.setLeft(BoundaryCondition1D<Real>::Dirichlet(0.0));
// bc.setRight(BoundaryCondition1D<Real>::Dirichlet(0.0));
```

### Step 3: Create Solver and Set Source

```cpp
PoissonSolver1D<Real> solver(grid, bc);

// Set source term f(x) = sin(πx)
solver.setSource([](Real x) {
    return std::sin(M_PI * x);
});
```

### Step 4: Solve

```cpp
// Solve the system (uses built-in Conjugate Gradient)
auto u = solver.solve();

std::cout << "Solution computed!\n";

// Find max value manually
Real maxVal = 0.0;
for (int i = 0; i < grid.numNodes(); ++i) {
    maxVal = std::max(maxVal, u(i));
}
std::cout << "Max value: " << maxVal << "\n";
// Output: Max value: ~0.101 (exact: 1/π² ≈ 0.10132)
```

### Complete Code

```cpp
#include <iostream>
#include <cmath>
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/elliptic/PoissonSolver.h"
#include "mml/pde/solvers/ConjugateGradient.h"

int main() {
    using namespace MML::PDE;
    using Real = double;
    
    // 1. Create grid
    Interval<Real> domain(0.0, 1.0);
    Grid1D<Real> grid(domain, 100);
    
    // 2. Set boundary conditions
    auto bc = homogeneousDirichlet1D<Real>();
    
    // 3. Create solver and set source
    PoissonSolver1D<Real> solver(grid, bc);
    solver.setSource([](Real x) { return std::sin(M_PI * x); });
    
    // 4. Solve
    auto u = solver.solve();
    
    // Find max value
    Real maxVal = 0.0;
    for (int i = 0; i < grid.numNodes(); ++i) {
        maxVal = std::max(maxVal, u(i));
    }
    std::cout << "Max value: " << maxVal << " (exact: " << 1.0/(M_PI*M_PI) << ")\n";
    
    return 0;
}
```

---

## 🔥 Example 2: Solve 1D Heat Equation

**Problem:** Solve $\frac{\partial u}{\partial t} = 0.01 \frac{\partial^2 u}{\partial x^2}$ on $[0,1]$ with $u(x,0) = \sin(\pi x)$

### Step 1: Create Grid and BCs

```cpp
Interval<Real> domain(0.0, 1.0);
Grid1D<Real> grid(domain, 100);
auto bc = homogeneousDirichlet1D<Real>();
```

### Step 2: Create Heat Solver

```cpp
Real alpha = 0.01;  // Thermal diffusivity
HeatSolver1D<Real> solver(grid, bc, alpha);
```

### Step 3: Set Initial Condition

```cpp
solver.setInitialCondition([](Real x) {
    return std::sin(M_PI * x);
});
```

### Step 4: Solve for Time Interval

```cpp
Real dt = 0.001;     // Time step
Real tFinal = 1.0;   // Solve until t = 1.0

// Use Crank-Nicolson (2nd order, unconditionally stable)
auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);

std::cout << "Heat equation solved!\n";
std::cout << "Steps taken: " << result.steps << "\n";
std::cout << "Final time: " << result.finalTime << "\n";
std::cout << "Final max temp: " << result.maxTemperature << "\n";
// Output: Final max temp: ~0.0004 (exponential decay: exp(-π²αt))
```

### Complete Code

```cpp
#include <iostream>
#include <cmath>
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/parabolic/HeatSolver.h"

int main() {
    using namespace MML::PDE;
    using Real = double;
    
    // 1. Create grid and boundary conditions
    Interval<Real> domain(0.0, 1.0);
    Grid1D<Real> grid(domain, 100);
    auto bc = homogeneousDirichlet1D<Real>();
    
    // 2. Create heat solver
    Real alpha = 0.01;
    HeatSolver1D<Real> solver(grid, bc, alpha);
    
    // 3. Set initial condition
    solver.setInitialCondition([](Real x) { return std::sin(M_PI * x); });
    
    // 4. Solve
    Real dt = 0.001, tFinal = 1.0;
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    
    std::cout << "Steps: " << result.steps << ", Final max: " << result.maxTemperature << "\n";
    
    return 0;
}
```

---

## 📐 Example 3: Solve 2D Laplace Equation

**Problem:** Solve $\nabla^2 u = 0$ on unit square with $u = g$ on boundary

### Code

```cpp
#include <iostream>
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/elliptic/PoissonSolver.h"
#include "mml/pde/solvers/ConjugateGradient.h"

int main() {
    using namespace MML::PDE;
    using Real = double;
    
    // 1. Create 50×50 grid on [0,1]×[0,1]
    Rectangle<Real> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<Real> grid(domain, 50, 50);
    
    // 2. Set non-homogeneous Dirichlet boundary conditions
    BoundaryConditions2D<Real> bc;
    bc.setBottom(BoundaryCondition2D<Real>::Dirichlet(0.0));
    bc.setTop(BoundaryCondition2D<Real>::Dirichlet(
        std::function<Real(Real, Real)>([](Real x, Real) { 
            return std::sin(M_PI * x); 
        })
    ));
    bc.setLeft(BoundaryCondition2D<Real>::Dirichlet(0.0));
    bc.setRight(BoundaryCondition2D<Real>::Dirichlet(0.0));
    
    // 3. Create Poisson solver (f = 0 → Laplace)
    PoissonSolver2D<Real> solver(grid, bc);
    solver.setSource([](Real x, Real y) { return 0.0; });  // Laplace: ∇²u = 0
    
    // 4. Solve (uses built-in CG)
    auto u = solver.solve();
    
    std::cout << "2D Laplace equation solved on " << grid.numNodes() << " nodes\n";
    
    // Find max value
    Real maxVal = 0.0;
    for (int i = 0; i < grid.numNodesX(); ++i) {
        for (int j = 0; j < grid.numNodesY(); ++j) {
            maxVal = std::max(maxVal, u(i, j));
        }
    }
    std::cout << "Max value: " << maxVal << "\n";
    
    return 0;
}
```

**Result:** Smooth interpolation from bottom (u=0) to top (u=sin(πx)) satisfying Laplace equation!

---

## 🏠 Example 4: 2D Heat Equation - Room Cooling

**Problem:** Simulate room temperature when you open a window (cold air enters from top)

### Code

```cpp
#include <iostream>
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include "mml/pde/parabolic/HeatSolver.h"

int main() {
    using namespace MML::PDE;
    using Real = double;
    
    // 1. Create grid: 5m × 3m room, 50×30 cells
    Rectangle<Real> domain(0.0, 5.0, 0.0, 3.0);
    Grid2D<Real> room(domain, 50, 30);
    
    // 2. Boundary conditions
    BoundaryConditions2D<Real> bc;
    // Top: open window at 0°C
    bc.setTop(BoundaryCondition2D<Real>::Dirichlet(0.0));
    // Other sides: insulated (Neumann with zero flux)
    bc.setBottom(BoundaryCondition2D<Real>::Neumann(0.0));
    bc.setLeft(BoundaryCondition2D<Real>::Neumann(0.0));
    bc.setRight(BoundaryCondition2D<Real>::Neumann(0.0));
    
    // 3. Create solver (air thermal diffusivity)
    Real alpha = 0.01;
    HeatSolver2D<Real> solver(room, bc, alpha);
    
    // 4. Initial condition: room at 20°C
    solver.setInitialCondition([](Real, Real) { return 20.0; });
    
    // 5. Solve for 100 seconds
    Real dt = 0.1;  // 0.1 second time steps
    Real tFinal = 100.0;
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    
    std::cout << "Room cooling simulation complete!\n";
    std::cout << "Time simulated: " << result.finalTime << " seconds\n";
    std::cout << "Temperature now ranges from " << result.minTemperature 
              << "°C to " << result.maxTemperature << "°C\n";
    std::cout << "Total energy: " << result.totalEnergy << "\n";
    
    return 0;
}
```

**Output Example:**
```
Room cooling simulation complete!
Time simulated: 100 seconds
Temperature now ranges from 5.2°C to 18.7°C
Total energy: 1234.5
```

💡 **Visualization Tip:** Export solution at different times using observer callback (see Advanced section).

---

## 🎯 Time Stepping Schemes Comparison

```cpp
// Forward Euler: Explicit, fast, conditionally stable
auto result1 = solver.solve(tFinal, dt, TimeScheme::ForwardEuler);

// Backward Euler: Implicit, unconditionally stable, O(Δt) accurate
auto result2 = solver.solve(tFinal, dt, TimeScheme::BackwardEuler);

// Crank-Nicolson: Implicit, O(Δt²) accurate, RECOMMENDED
auto result3 = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);

// Theta method: Generalized (θ=0→FE, θ=0.5→CN, θ=1→BE)
solver.setTheta(0.5);  // Same as Crank-Nicolson
auto result4 = solver.solve(tFinal, dt, TimeScheme::Theta);
```

| Scheme | Order | Stability | Cost/Step | Recommended? |
|--------|-------|-----------|-----------|--------------|
| Forward Euler | O(Δt) | CFL: Δt ≤ h²/(2αd) | Cheap | Education only |
| Backward Euler | O(Δt) | Unconditional | Moderate | Large time steps |
| Crank-Nicolson | O(Δt²) | Unconditional | Moderate | ⭐ Production |

---

## 🔧 Advanced: Monitoring Solution Evolution

Use observer callback to capture solution at intermediate times:

```cpp
#include <vector>

std::vector<Real> times;
std::vector<std::vector<Real>> snapshots;

// Observer function called at each time step
auto observer = [&](Real t, const std::vector<Real>& u) {
    times.push_back(t);
    snapshots.push_back(u);  // Store solution
    
    // Find max value
    Real maxVal = *std::max_element(u.begin(), u.end());
    std::cout << "t = " << t << ", max = " << maxVal << "\n";
};

// Set observer before solving
solver.setObserver(observer);

// Solve
auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);

std::cout << "Captured " << snapshots.size() << " snapshots\n";

// Now you can analyze or export snapshots for visualization
```

---

## 📊 Exporting Solutions for Visualization

```cpp
#include <fstream>

// Export 2D solution to CSV
void exportToCSV(const GridFunction2D<Real>& u, const Grid2D<Real>& grid, 
                 const std::string& filename) {
    std::ofstream file(filename);
    file << "x,y,u\n";
    
    for (int j = 0; j <= grid.ny(); ++j) {
        for (int i = 0; i <= grid.nx(); ++i) {
            Real x = grid.x(i);
            Real y = grid.y(j);
            int idx = i + j * (grid.nx() + 1);
            file << x << "," << y << "," << u[idx] << "\n";
        }
    }
    file.close();
}

// Usage
exportToCSV(u, grid, "solution.csv");
```

Then visualize with Python:
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('solution.csv')
plt.tricontourf(df['x'], df['y'], df['u'], levels=20)
plt.colorbar(label='u')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution')
plt.show()
```

---

## 🚀 Performance Tips

### 1. Choose the Right Solver

```cpp
// PoissonSolver uses built-in Conjugate Gradient (optimal for SPD matrices)
auto u = poissonSolver.solve();

// For custom tolerance:
auto u = poissonSolver.solve(1e-10, 5000);  // tolerance, max iterations

// The solver automatically chooses the best method based on problem type
```

### 2. Add a Preconditioner

```cpp
// Solvers use built-in optimizations and preconditioners
// For custom settings, you can adjust solver tolerance:
auto u = poissonSolver.solve(1e-8, 2000);  // tighter tolerance, more iterations

// For heat equations, choose time scheme carefully:
// Crank-Nicolson is usually optimal (2nd order + unconditionally stable)
auto result = heatSolver.solve(tFinal, dt, TimeScheme::CrankNicolson);
```

### 3. Configure Solver Tolerance

```cpp
// Configure solver tolerance and max iterations
auto u = poissonSolver.solve(1e-8,   // Relative residual tolerance
                             1000);   // Maximum iterations

// Default is solve() which uses tolerance=1e-10, maxIter=10000
auto u = poissonSolver.solve();  // Good defaults for most problems
```

### 4. Grid Resolution Guidelines

| Problem Size | 2D Grid | 3D Grid | Memory | Solve Time* |
|--------------|---------|---------|--------|-------------|
| Small | 50×50 | 25×25×25 | ~2 MB | 0.01 s |
| Medium | 200×200 | 50×50×50 | ~30 MB | 1 s |
| Large | 500×500 | 100×100×100 | ~200 MB | 30 s |
| Very Large | 1000×1000 | 200×200×200 | ~1 GB | 5 min |

*Approximate, without preconditioner, single-threaded

---

## ⚠️ Common Mistakes

### Mistake 1: Using `nx()` Instead of `numNodesX()`

```cpp
// ❌ WRONG: nx() returns CELL count
for (int i = 0; i < grid.nx(); ++i) {  // Misses last node!
    // ...
}

// ✅ CORRECT: numNodesX() returns NODE count (nx + 1)
for (int i = 0; i <= grid.nx(); ++i) {  // Correct
    // ...
}

// Or better yet:
for (int i = 0; i < grid.numNodesX(); ++i) {  // Also correct
    // ...
}
```

### Mistake 2: CFL Violation with Forward Euler

```cpp
Rectangle<Real> domain(0.0, 1.0, 0.0, 1.0);
Grid2D<Real> grid(domain, 100, 100);
Real h = grid.dx();  // 0.01
Real alpha = 0.01;
Real dt = 0.01;  // ❌ TOO LARGE!

// CFL condition: dt ≤ h²/(2αd) = 0.01²/(2×0.01×2) = 0.0025
// Using dt=0.01 will cause instability!

// ✅ CORRECT:
Real dtMax = h * h / (2 * alpha * 2);  // Factor of 2 for 2D
Real dt = 0.8 * dtMax;  // Use 80% of limit for safety

// Or just use implicit scheme (no CFL restriction):
auto result = solver.solve(tFinal, 0.01, TimeScheme::CrankNicolson);  // ✅
```

### Mistake 3: Wrong Boundary Condition API

```cpp
// ❌ WRONG: No such method
auto value = bc.getValue(x, y, t);

// ✅ CORRECT: Use get(side).value(x, y)
auto bcTop = bc.get(BoundarySide2D::Top);
auto value = bcTop.value(x, y);
```

---

## 🔍 Troubleshooting

| Problem | Possible Cause | Solution |
|---------|---------------|----------|
| **CG not converging** | Matrix not SPD | Check BCs, add preconditioner, increase max iterations |
| **Wrong results** | BC mismatch | Verify BC setup, check sign convention |
| **Segmentation fault** | Grid bounds error | Use `numNodesX()` not `nx()`, check array access |
| **Too slow** | No preconditioner | Add SSOR or ILU preconditioner |
| **Instability (heat)** | CFL violation | Reduce dt or use implicit scheme |

👉 **[Full Troubleshooting Guide](Troubleshooting.md)**

---

## 📚 Next Steps

You now know the basics! Here's where to go next:

### Learn More About Components
- **[Sparse Matrices](Sparse_Matrices.md)** - CSR/CSC/COO formats in depth
- **[Grid Infrastructure](Grid_Infrastructure.md)** - Grid internals and stencils
- **[Boundary Conditions](Boundary_Conditions.md)** - All BC types explained
- **[Iterative Solvers](Iterative_Solvers.md)** - CG/BiCGSTAB/GMRES theory and practice

### Dive into PDE Types
- **[Elliptic PDEs](Elliptic_PDEs.md)** - Complete Poisson/Laplace guide
- **[Parabolic PDEs](Parabolic_PDEs.md)** - Complete heat equation guide

### Practice with Examples
- **[Examples Gallery](Examples_Gallery.md)** - 10+ complete working examples
- **[Performance Guide](Performance_Guide.md)** - Optimization techniques

---

## 🎓 Summary

You learned how to:
- ✅ Create 1D and 2D grids
- ✅ Set up Dirichlet and Neumann boundary conditions
- ✅ Solve Poisson equations (elliptic PDEs)
- ✅ Solve heat equations (parabolic PDEs)
- ✅ Choose appropriate time-stepping schemes
- ✅ Monitor solution evolution with observers
- ✅ Export solutions for visualization

**You're ready to solve PDEs!** 🎉

---

## 💡 Quick Reference Card

```cpp
// Elliptic (Poisson): -∇²u = f
Rectangle<Real> domain(xmin, xmax, ymin, ymax);
Grid2D<Real> grid(domain, nx, ny);
auto bc = homogeneousDirichlet2D<Real>();
PoissonSolver2D<Real> solver(grid, bc);
solver.setSource([](Real x, Real y) { return f(x,y); });
auto u = solver.solve();

// Parabolic (Heat): ∂u/∂t = α∇²u
HeatSolver2D<Real> solver(grid, bc, alpha);
solver.setInitialCondition([](Real x, Real y) { return initial(x,y); });
auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
```

---

*Ready for more? Check out the **[Examples Gallery](Examples_Gallery.md)** for complete working code!*

---

*Last updated: January 1, 2026*
