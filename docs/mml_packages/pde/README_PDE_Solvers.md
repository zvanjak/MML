# PDE Solvers - Finite Difference Methods for Partial Differential Equations

**MinimalMathLibrary (MML) - PDE Module**

Comprehensive toolkit for solving partial differential equations using finite difference methods on structured grids with iterative linear solvers.

---

## 🎯 Overview

The MML PDE module provides production-ready solvers for the three fundamental types of PDEs:

| PDE Type | Equation | Applications | Status |
|----------|----------|--------------|--------|
| **Elliptic** | $-\nabla^2 u = f$ | Steady-state heat, electrostatics, stress | ✅ Complete |
| **Parabolic** | $\frac{\partial u}{\partial t} = \alpha\nabla^2 u$ | Transient diffusion, heat flow | ✅ Complete |
| **Hyperbolic** | $\frac{\partial^2 u}{\partial t^2} = c^2\nabla^2 u$ | Wave propagation, acoustics | 🚧 Planned |

---

## ✨ Key Features

- ✅ **Sparse Matrix Infrastructure** - CSR, CSC, COO formats optimized for PDE discretization
- ✅ **Structured Grid System** - Uniform 1D, 2D, and 3D grids with flexible boundary conditions
- ✅ **Iterative Solvers** - CG, BiCGSTAB, GMRES with preconditioners (Jacobi, SSOR, ILU)
- ✅ **Elliptic PDE Solvers** - Poisson and Laplace equations in 1D/2D/3D
- ✅ **Parabolic PDE Solvers** - Heat equation with multiple time-stepping schemes
- ✅ **Boundary Conditions** - Dirichlet, Neumann, Robin, and periodic BCs
- ✅ **Header-Only** - Easy integration, no linking required
- ✅ **Thoroughly Tested** - 2500+ unit tests including PDE-specific validation

---

## 🚀 Quick Start

Solve the 2D Poisson equation $-\nabla^2 u = \sin(\pi x)\sin(\pi y)$ with zero Dirichlet boundaries:

```cpp
#include "mml/pde/elliptic/PoissonSolver.h"

using namespace MML::PDE;
using Real = double;

int main() {
    // Create 50×50 grid on [0,1]×[0,1]
    Grid2D<Real> grid(0, 1, 0, 1, 50, 50);
    
    // Zero Dirichlet boundary conditions
    auto bc = homogeneousDirichlet2D<Real>();
    
    // Create solver and set source term
    PoissonSolver2D<Real> solver(grid, bc);
    solver.setSource([](Real x, Real y) {
        return std::sin(M_PI * x) * std::sin(M_PI * y);
    });
    
    // Solve using Conjugate Gradient
    ConjugateGradient<Real> cg;
    auto u = solver.solve(cg);
    
    std::cout << "Solution computed on " << grid.numNodes() << " nodes\n";
    std::cout << "Max value: " << u.max() << "\n";
    
    return 0;
}
```

**Result:** Solution computed in milliseconds with O(h²) spatial accuracy!

👉 **[See Quick Start Guide](Quick_Start_Guide.md)** for more examples and step-by-step tutorial.

---

## 📐 Module Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        PDE SOLVER MODULE                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐         │
│  │   Sparse     │  │     Grid     │  │  Boundary    │         │
│  │   Matrices   │  │ Infrastructure│  │  Conditions  │         │
│  │              │  │              │  │              │         │
│  │ • CSR/CSC/COO│  │ • Grid1D/2D/3D│  │ • Dirichlet  │         │
│  │ • SpMV       │  │ • Stencils   │  │ • Neumann    │         │
│  │ • Assembly   │  │ • Iteration  │  │ • Robin      │         │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘         │
│         │                 │                 │                 │
│         └─────────┬───────┴─────────┬───────┘                 │
│                   │                 │                         │
│         ┌─────────▼─────────────────▼─────────┐               │
│         │    Iterative Linear Solvers          │               │
│         │  • CG (SPD systems)                  │               │
│         │  • BiCGSTAB (non-symmetric)          │               │
│         │  • GMRES (general)                   │               │
│         │  • Preconditioners (Jacobi/SSOR/ILU) │               │
│         └─────────┬─────────────────┬──────────┘               │
│                   │                 │                         │
│         ┌─────────▼─────────┐  ┌────▼──────────┐              │
│         │   Elliptic PDEs   │  │ Parabolic PDEs│              │
│         │                   │  │               │              │
│         │ • PoissonSolver1D │  │ • HeatSolver1D│              │
│         │ • PoissonSolver2D │  │ • HeatSolver2D│              │
│         │ • PoissonSolver3D │  │ • HeatSolver3D│              │
│         │                   │  │ • Time schemes│              │
│         └───────────────────┘  └───────────────┘              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 📚 Documentation Index

### Getting Started
- **[Quick Start Guide](Quick_Start_Guide.md)** - 5-minute tutorial with examples ⭐
- **[Examples Gallery](Examples_Gallery.md)** - Complete working examples with visualization

### Core Components
- **[Sparse Matrices](Sparse_Matrices.md)** - CSR/CSC/COO formats, when to use each
- **[Grid Infrastructure](Grid_Infrastructure.md)** - Structured grids, indexing, stencils
- **[Boundary Conditions](Boundary_Conditions.md)** - Dirichlet, Neumann, Robin, periodic
- **[Iterative Solvers](Iterative_Solvers.md)** - CG, BiCGSTAB, GMRES, preconditioners

### PDE Solvers
- **[Elliptic PDEs](Elliptic_PDEs.md)** - Poisson and Laplace equations ⚡
- **[Parabolic PDEs](Parabolic_PDEs.md)** - Heat and diffusion equations 🔥
- **[Hyperbolic PDEs](Hyperbolic_PDEs.md)** - Wave and advection (planned) 🌊

### Advanced Topics
- **[Performance Guide](Performance_Guide.md)** - Optimization, benchmarks, scaling 📊
- **[Mathematical Background](Mathematical_Background.md)** - Theory and fundamentals 📐
- **[API Reference](API_Reference.md)** - Complete class and function reference 📚

### Support
- **[Troubleshooting](Troubleshooting.md)** - Common issues and solutions 🔧 ⭐

---

## 🎓 What Can You Solve?

### Elliptic Problems (Steady-State)
- **Electrostatics:** Compute electric potential from charge distributions
- **Steady-State Heat:** Temperature distribution in solids
- **Structural Analysis:** Stress and displacement in materials
- **Fluid Flow:** Stream functions and potential flow
- **Image Processing:** Poisson image editing

### Parabolic Problems (Time-Dependent Diffusion)
- **Transient Heat Transfer:** Temperature evolution over time
- **Mass Diffusion:** Concentration gradients in materials
- **Option Pricing:** Black-Scholes PDE for financial derivatives
- **Image Processing:** Anisotropic diffusion filtering
- **Biological Processes:** Reaction-diffusion systems

### Example: Room Cooling Simulation 🏠

```cpp
// Simulate temperature in a room when you open a window
Grid2D<Real> room(0, 5, 0, 3, 50, 30);  // 5m × 3m room

// Initial condition: warm room at 20°C
auto u0 = [](Real x, Real y) { return Real(20.0); };

// Boundary: cold window at top, insulated walls elsewhere
BoundaryConditions2D<Real> bc;
bc.setTop(BoundaryType::Dirichlet, [](Real x, Real) { return Real(0.0); });  // 0°C window
bc.setBottom(BoundaryType::Neumann, [](Real, Real) { return Real(0.0); });   // Insulated floor
bc.setLeft(BoundaryType::Neumann, [](Real, Real) { return Real(0.0); });     // Insulated wall
bc.setRight(BoundaryType::Neumann, [](Real, Real) { return Real(0.0); });    // Insulated wall

// Solve for 100 seconds
HeatSolver2D<Real> solver(room, Real(0.01), bc);  // α = 0.01 (air thermal diffusivity)
auto result = solver.solve(u0, Real(0.1), Real(100.0), TimeScheme::CrankNicolson);

std::cout << "Room cooling simulation complete!\n";
std::cout << "Temperature range: [" << result.minTemperature 
          << ", " << result.maxTemperature << "] °C\n";
```

---

## ⚙️ Performance Characteristics

| Grid Size | DOF | Matrix NNZ | CG Iterations | Solve Time* | Memory |
|-----------|-----|------------|---------------|-------------|--------|
| 2D: 50×50 | 2,500 | ~12K | ~50 | 5 ms | 1 MB |
| 2D: 100×100 | 10,000 | ~50K | ~100 | 20 ms | 4 MB |
| 2D: 500×500 | 250,000 | ~1.2M | ~500 | 2 sec | 80 MB |
| 3D: 50×50×50 | 125,000 | ~800K | ~200 | 1 sec | 50 MB |

*Single-threaded on Intel i7, double precision, without preconditioner

💡 **Tip:** Add a preconditioner (SSOR or ILU) to reduce iterations by 2-5×!

---

## 🔬 Numerical Accuracy

| Method | Spatial Order | Temporal Order | Stability |
|--------|---------------|----------------|-----------|
| Poisson (5-point) | O(h²) | N/A | Unconditional |
| Heat (Forward Euler) | O(h²) | O(Δt) | CFL: Δt ≤ h²/(2α·d) |
| Heat (Backward Euler) | O(h²) | O(Δt) | Unconditional |
| Heat (Crank-Nicolson) | O(h²) | O(Δt²) | Unconditional ⭐ |

---

## 🛠️ System Requirements

### Required
- **C++17** or later
- Standard library (`<vector>`, `<functional>`, `<cmath>`, etc.)

### Optional
- **CMake 3.14+** (for building tests and examples)
- **Catch2** (for running unit tests, auto-downloaded by CMake)

### Tested Compilers
- ✅ MSVC 2022 (Windows)
- ✅ GCC 13+ (Linux)
- ✅ Clang 14+ (macOS)

---

## 📦 Installation

### Header-Only Integration

```cpp
// Include the headers you need:
#include "mml/pde/sparse/SparseMatrixCSR.h"
#include "mml/pde/grid/Grid2D.h"
#include "mml/pde/solvers/ConjugateGradient.h"
#include "mml/pde/elliptic/PoissonSolver.h"
#include "mml/pde/parabolic/HeatSolver.h"
```

### CMake Integration

```cmake
# Add MML to your project
add_subdirectory(path/to/MinimalMathLibrary)

# Link your target
target_link_libraries(your_target PRIVATE MML)
```

---

## 🧪 Testing

The PDE module includes comprehensive tests:

```bash
# Build and run all tests
cd build
cmake --build . --target MML_Tests
ctest -R pde -V

# Run specific PDE test suite
./tests/Release/MML_Tests.exe "[poisson]"
./tests/Release/MML_Tests.exe "[heat]"
```

**Test Coverage:**
- 📊 Sparse matrices: ~150 tests
- 📐 Grid infrastructure: ~80 tests
- 🔁 Iterative solvers: ~120 tests
- ⚡ Elliptic PDEs: ~50 tests
- 🔥 Parabolic PDEs: ~40 tests

---

## 🎯 Design Philosophy

### 1. **Correctness First**
- Verified against manufactured solutions
- Extensive unit tests for every component
- O(h²) spatial accuracy guaranteed

### 2. **Performance Matters**
- Optimized sparse matrix formats (CSR for SpMV)
- Minimal memory allocations
- Cache-friendly data structures

### 3. **Easy to Use**
- Clear, intuitive API
- Comprehensive documentation
- Working examples for every use case

### 4. **Flexible and Extensible**
- Template-based for any floating-point type
- Mix-and-match components (solvers, preconditioners, BCs)
- Easy to add custom boundary conditions

---

## 🤝 Common Workflows

### Workflow 1: Solve Elliptic PDE
1. Create `Grid2D` or `Grid3D` with desired resolution
2. Set up `BoundaryConditions2D/3D` (Dirichlet/Neumann/Robin)
3. Create `PoissonSolver2D/3D` with grid and BCs
4. Set source term `f(x,y)` or `f(x,y,z)`
5. Choose iterative solver (`ConjugateGradient`, `BiCGSTAB`, `GMRES`)
6. Solve and get `GridFunction2D/3D` result

### Workflow 2: Solve Parabolic PDE
1. Create `Grid1D/2D/3D` with spatial resolution
2. Set up `BoundaryConditions` 
3. Create `HeatSolver1D/2D/3D` with grid, diffusivity α, and BCs
4. Optionally set source term `f(x,t)` or `f(x,y,t)`
5. Choose time integration scheme (`CrankNicolson` recommended)
6. Set initial condition `u₀(x)` or `u₀(x,y)`
7. Solve for time interval `[0, T]` with time step Δt
8. Get `HeatSolveResult` with final solution and diagnostics

---

## ⚠️ Common Pitfalls

| Issue | Problem | Solution |
|-------|---------|----------|
| **Grid bounds** | Using `nx()` instead of `numNodesX()` | Always use `numNodesX()` for loop bounds! |
| **CG divergence** | Matrix not SPD | Check BC assembly, use preconditioner |
| **CFL violation** | Explicit scheme unstable | Use Δt ≤ h²/(2αd) or switch to implicit |
| **BC API** | Using wrong BC accessor | Use `bc.get(side).value(x, y)` pattern |
| **Wrong sign** | Sign error in PDE | We use $-\nabla^2 u = f$ (negative for SPD) |

👉 **[Full Troubleshooting Guide](Troubleshooting.md)**

---

## 🗺️ Roadmap

### ✅ Completed (v1.0)
- Sparse matrix infrastructure
- Structured grid system
- Iterative solvers (CG, BiCGSTAB, GMRES)
- Preconditioners (Jacobi, SSOR, ILU)
- Elliptic PDE solvers (Poisson/Laplace)
- Parabolic PDE solvers (Heat equation)
- Comprehensive documentation

### 🚧 In Progress
- Hyperbolic PDE solvers (Wave, Advection)
- ADI (Alternating Direction Implicit) methods
- Multigrid solvers

### 🔮 Future
- Non-uniform grids
- Adaptive mesh refinement
- Parallel solvers (OpenMP, MPI)
- GPU acceleration
- Nonlinear PDE solvers (Newton-Krylov)
- Finite element methods (FEM)

---

## 📖 References

### Textbooks
- **Numerical Recipes** (Press et al.) - Classic reference for scientific computing
- **Finite Difference Methods for PDEs** (LeVeque) - Excellent FD introduction
- **Iterative Methods for Sparse Linear Systems** (Saad) - Comprehensive solver reference

### Related MML Modules
- **Linear Algebra** - Dense matrix operations, decompositions
- **Numerical Methods** - ODE solvers, integration, root finding
- **Optimization** - Minimization, least squares

---

## 👥 Contributing

Found a bug? Have a feature request? Want to improve documentation?

- 📝 **Documentation:** Help us improve these guides!
- 🐛 **Bug Reports:** Include minimal reproducible example
- ✨ **Feature Requests:** Describe use case and expected API
- 🔧 **Pull Requests:** Follow existing code style

---

## 📄 License

MinimalMathLibrary is released under the MIT License — free for personal, academic, and commercial use.

See [LICENSE.md](../../LICENSE.md) for details.

---

## 🙏 Acknowledgments

The PDE module was inspired by:
- Numerical Recipes implementations
- PETSc design philosophy
- SciPy sparse matrix conventions
- Deal.II finite element library

Special thanks to the numerical analysis community for decades of research in PDE solvers!

---

## 🚀 Next Steps

1. **[Quick Start Guide](Quick_Start_Guide.md)** - Solve your first PDE in 5 minutes
2. **[Examples Gallery](Examples_Gallery.md)** - Browse complete working examples
3. **[Elliptic PDEs](Elliptic_PDEs.md)** - Deep dive into Poisson solvers
4. **[Parabolic PDEs](Parabolic_PDEs.md)** - Master heat equation solvers

**Happy solving!** ⚡🔥🌊

---

*Last updated: January 1, 2026*  
*MML PDE Module v1.0*
