<div align="center">

# 🔢 MML — Minimal Math Library

### **The Complete C++ Numerical Computing Toolkit**

*Single-header • Cross-platform • Visualization included*

[![Ubuntu](https://github.com/zvanjak/MML/workflows/Ubuntu/badge.svg)](https://github.com/zvanjak/MML/actions?query=workflow%3AUbuntu)
[![Windows](https://github.com/zvanjak/MML/workflows/Windows/badge.svg)](https://github.com/zvanjak/MML/actions?query=workflow%3AWindows)
[![macOS](https://github.com/zvanjak/MML/workflows/macOS/badge.svg)](https://github.com/zvanjak/MML/actions?query=workflow%3AmacOS)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)
[![Tests](https://img.shields.io/badge/tests-4540%20passing-brightgreen.svg)](tests/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE.md)

---

##  Mission 

**🚀 Just `#include "MML.h"` and compute** — vectors, matrices, tensors, ODE solvers, eigenvalues, integration, and more!

[Quick Start](#-quick-start) • [Features](#-key-features) • [Documentation](#-documentation) • [Examples](#-real-examples) • [Visualization](#-visualization-suite)

</div>

---

## 🎯 What is MML?

MML is a **comprehensive, single-header C++ mathematical library** for numerical computing. With just one `#include "MML.h"` directive, you get access to a complete toolkit of mathematical objects and algorithms — from basic vectors and matrices to ODE solvers and field operations.
And as a bonus, you also get cross-platform visualization tools for plotting functions, curves, surfaces, and vector fields!

<table>
<tr>
<td width="50%">

### The Problem

Most C++ math libraries require:
- Complex build systems and dependencies
- Linking against multiple libraries
- Platform-specific configuration
- Steep learning curves

**Result:** Possibly hours of setup before writing actual code.

</td>
<td width="50%">

### The MML Solution

**Zero-friction integration:**

```cpp
#include "MML.h"
// That's it. Start computing.
```

- ✅ Single header file
- ✅ Pure C++17, no dependencies
- ✅ Works on Windows, Linux, Mac
- ✅ 4,540 tests ensure correctness

</td>
</tr>
</table>

If needed, you can also use it piece-wise by including only selected headers from `mml/` directory.

---

## 🏛️ Design Philosophy

MML is built on three core principles:

<table>
<tr>
<td align="center" width="33%">

### 🎯 Completeness & Simplicity

A **comprehensive toolkit** covering vectors, matrices, tensors, ODEs, integration, field operations, and more. **Intuitive syntax** makes mathematical objects first-class citizens in C++.

</td>
<td align="center" width="33%">

### 🔬 Correctness & Precision

**4,540 unit tests** validate every algorithm against known analytical solutions. Numerical methods **report achieved precision**, and dedicated testbeds demonstrate accuracy.

</td>
<td align="center" width="33%">

### ⚡ Trivial Integration & Visualization

**Single header, zero dependencies** — just `#include "MML.h"` and start computing immediately. **Cross-platform visualizers** for functions, surfaces, vector fields, and particle systems on Windows, Linux and Mac.

</td>
</tr>
</table>

### Flagship Example: Verify Gauss's Divergence Theorem

📄 *[View full source](src/readme_examples/readme00_fundamental_theorems.cpp)*

```cpp
// Verify ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS over a unit cube
// Define vector field F(x,y,z) = (x², y², z²)
VectorFunction<3> F([](const VectorN<Real, 3>& p) {
    return VectorN<Real, 3>{ p[0]*p[0], p[1]*p[1], p[2]*p[2] };
});

// Compute divergence NUMERICALLY - MML calculates ∇·F automatically!
ScalarFunctionFromStdFunc<3> divF([&F](const VectorN<Real, 3>& p) {
    return VectorFieldOperations::DivCart<3>(F, p);
});

// Volume integral: ∫∫∫(∇·F)dV
auto y_lo = [](Real) { return 0.0; };  auto y_hi = [](Real) { return 1.0; };
auto z_lo = [](Real,Real) { return 0.0; };  auto z_hi = [](Real,Real) { return 1.0; };
Real volIntegral = Integrate3D(divF, GAUSS10, 0, 1, y_lo, y_hi, z_lo, z_hi).value;

// Surface integral: ∮∮(F·n̂)dS through all 6 faces
Cube3D unitCube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
Real surfIntegral = SurfaceIntegration::SurfaceIntegral(F, unitCube, 1e-8);

std::cout << "Volume integral:  " << volIntegral << "\n";   // 3.0000000000
std::cout << "Surface integral: " << surfIntegral << "\n";  // 3.0000000000
std::cout << "Error: " << std::abs(volIntegral - surfIntegral) << "\n";  // 9.77e-15 ✓
```

---

## 🚀 Quick Start

### Installation

**Option 1: Simply get single header (Recommended)**
```bash
# Download MML.h directly from the repository
curl -O https://raw.githubusercontent.com/zvanjak/MML/master/mml/single_header/MML.h

# Include in your project
#include "MML.h"
```

**Option 2: Get whole repository (if you need individual modules)**
```bash
git clone https://github.com/zvanjak/MML.git
cd MML
cmake -B build
cmake --build build
```


**Option 3: Visual Studio Code**

1. Open VS Code and press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
2. Type **"Git: Clone"** and press Enter
3. Paste the repository URL:
   ```
   https://github.com/zvanjak/MML.git
   ```
4. Select a folder and open the cloned repository
5. When prompted, install recommended extensions (C/C++, CMake Tools)
6. Press `Ctrl+Shift+P` → **"CMake: Configure"** to set up the build
7. Press `F7` to build or use the CMake status bar

> 💡 **Tip:** The [CMake Tools extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools) provides IntelliSense, build buttons, and test integration out of the box.

### First Program

```cpp
#include "MML.h"
using namespace MML;

int main() {
    // Create a matrix and solve a linear system
    Matrix<Real> A{3, 3, {4, 1, 2, 
                          1, 5, 1, 
                          2, 1, 6}};
    Vector<Real> b{1, 2, 3};
    
    LUSolver<Real> solver(A);
    Vector<Real> x = solver.Solve(b);
    
    std::cout << "Solution: " << x << std::endl;
    std::cout << "Residual: " << (A * x - b).NormL2() << std::endl;
    
    return 0;
}
```

**Compile:**
```bash
g++ -std=c++17 -O3 myprogram.cpp -o myprogram
```

---

## ✨ Key Features

### 🏗️ Library Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                          MML.h (single header)                      │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│   mml/                                                              │
│   ├── base/        Vectors, Matrices, Tensors, Functions, Geometry  │
│   ├── core/        Derivation, Integration, Linear Solvers, Fields  │
│   ├── algorithms/  ODE, Root Finding, Interpolation, Eigen Solvers  │
│   ├── systems/     Dynamical Systems, Linear Systems, Attractors    │
│   ├── interfaces/  Abstract interfaces                              │
│   └── tools/       Visualization, Serialization, Console Printing   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

MML is organized into four main layers, each building on the previous one:

- **[Base](docs/base/README_Base.md)** — The mathematical foundation. Vectors, matrices, tensors, functions, polynomials, quaternions, and geometry primitives. These are the objects you compute with — designed with intuitive syntax so mathematical expressions in code read naturally.
- **[Core](docs/core/README_Core.md)** — Numerical operations on base objects. Numerical differentiation (up to 8th order accuracy), integration (1D/2D/3D with adaptive quadrature), linear system solvers (LU, QR, SVD, Cholesky), vector field operations (gradient, divergence, curl, Laplacian), coordinate transformations, and metric tensors.
- **[Algorithms](docs/algorithms/README_Algorithms.md)** — Higher-level problem-solving methods. ODE solvers (fixed and adaptive step), root finding (Bisection, Newton, Brent), eigenvalue decomposition, interpolation, curve fitting, path and surface integration, differential geometry, and function analysis.
- **[Tools](docs/tools/)** — Bridging computation and presentation. Console printing with publication-quality formatting (6 export formats), file serialization for all mathematical objects, and cross-platform visualizers for functions, surfaces, vector fields, and particle systems.


### 📐 Mathematical Objects

| Category | Types | Description |
|----------|-------|-------------|
| [**Vectors**](docs/base/Vectors.md) | `Vector<T>`, `VectorN<T,N>`, 2D/3D Cartesian, Polar, Spherical | Dynamic and fixed-size vectors in multiple coordinate systems |
| [**Matrices**](docs/base/Matrices.md) | `Matrix<T>`, `MatrixNM<N,M>`, Symmetric, Tridiagonal, Band | Full matrix algebra with specialized storage |
| [**Tensors**](docs/base/Tensors.md) | `Tensor2<N>` through `Tensor5<N>` | Rank 2-5 tensors for advanced computations |
| [**Functions**](docs/base/Functions.md) | `IRealFunction`, `IScalarFunction<N>`, `IVectorFunction<N>` | First-class function objects with calculus support |
| [**Curves & Surfaces**](docs/core/Curves_and_surfaces.md) | `ParametricCurve<N>`, `ParametricSurface<N>`, predefined 2D/3D | Parametric curves, surfaces, arc length, Frenet frames |
| [**Geometry**](docs/base/Geometry.md) | Points, Lines, Planes, Triangles, Bodies | 2D and 3D geometric primitives |
| [**Polynomials**](docs/base/Polynoms.md) | `Polynom<T>` | Generic polynomials over any field |
| [**Quaternions**](docs/base/Quaternions.md) | Full quaternion algebra | 3D rotations, SLERP interpolation |

### 🔢 Numerical Algorithms

| Category | Algorithms | Description |
|----------|------------|-------------|
| [**Linear Algebra**](docs/core/Linear_equations_solvers.md) | LUSolver, QRSolver, SVD, Cholesky | Matrix decompositions and solvers |
| [**Eigensolvers**](docs/algorithms/Eigen_solvers.md) | EigenSolver (symmetric & general) | Real matrix eigenvalue computation |
| [**Derivation**](docs/core/Derivation.md) | NDer1-8, NSecDer, NThirdDer, Gradient, Jacobian | 1st/2nd/3rd derivatives, O(h) to O(h⁸) accuracy |
| [**Integration 1D**](docs/core/Integration.md) | Trap, Simpson, Romberg, Gauss-Kronrod (G7K15, G10K21) | Adaptive quadrature with error estimates |
| [**Integration 2D/3D**](docs/core/Multidim_integration.md) | Integrate2D, Integrate3D, Monte Carlo | Multidimensional integration |
| [**Improper Integrals**](docs/core/Integration.md) | IntegrateUpperInf, IntegrateLowerInf, IntegrateInfInf | Semi-infinite and infinite bounds |
| [**Path & Surface**](docs/algorithms/Path_integration.md) | PathIntegration, SurfaceIntegration | Curve and surface integrals, flux |
| [**ODE Solvers**](docs/algorithms/Differential_equations_solvers.md) | ODESystemFixedStepSolver, Euler, RK4, Adaptive | Ordinary differential equations |
| [**Dynamical Systems**](docs/systems/) | DynamicalSystem, fixed points, Lyapunov exponents | Phase space analysis and stability |
| [**Root Finding**](docs/algorithms/Root_finding.md) | Bisection, Newton, Secant, Brent | Equation solving |
| [**Interpolation**](docs/base/Interpolated_functions.md) | LinearInterpRealFunc, SplineInterpRealFunc, PolynomInterpRealFunc | Function approximation |

### 🌀 Advanced Capabilities

| Feature | Description |
|---------|-------------|
| [**Field Operations**](docs/core/Field_operations.md) | `GradientCart`, `DivCart`, `CurlCart`, `LaplacianCart` in Cartesian, cylindrical, spherical |
| [**Coordinate Transforms**](docs/core/Coordinate_transformations.md) | `CoordTransfSphericalToCartesian`, `CoordTransfLorentzXAxis`, rotations |
| [**Metric Tensors**](docs/core/Metric_tensor.md) | General coordinate system support |
| [**Differential Geometry**](docs/algorithms/Differential_geometry.md) | Curvature, torsion, Frenet frames |
| [**Dynamical Systems**](docs/systems/) | Fixed point classification, Lyapunov exponents, bifurcation analysis |
| [**Function Analysis**](docs/algorithms/Function_analyzer.md) | Find roots, extrema, inflection points |

## 📚 Tools & Utilities

**Practical tools** for working with mathematical objects — the fourth layer of MML.

| Tool | Description |
|------|-------------|
| [**ConsolePrinter**](docs/tools/ConsolePrinter.md) | Beautiful table formatting, 6 export formats (TXT, CSV, JSON, HTML, LaTeX, Markdown), 5 border styles |
| [**Serializer**](docs/tools/Serializer.md) | Save functions, ODE solutions, particle simulations, vector fields to files |
| [**Visualizers**](docs/tools/Visualizers.md) | Cross-platform plotting for functions, fields, curves, surfaces, particles |
| **Random** | Random number generators for simulations and Monte Carlo |

**Why Tools Matter:**

The Tools layer bridges the gap between computation and presentation. Whether you're debugging numerical algorithms, preparing results for publication, or creating animations for teaching — these utilities make it effortless:

- **ConsolePrinter**: Format any vector, matrix, or table into publication-quality output with a single call
- **Serializer**: Persist simulation results for later analysis or visualization in external tools
- **Visualizers**: Launch real-time viewers directly from code — no manual data export needed



See the full [Visualization Gallery](#-visualization-gallery) below for all visualizer types.

---

## 🧪 Real Examples — Physics Simulations You Can Run

> **Self-contained, production-ready examples** demonstrating MML's capabilities in real-world physics simulations. Each example includes all physics code in-directory — no external dependencies.

### 🌌 [Example 00: N-Body Gravity](docs/examples/Example_00_N_body_gravity.md) — *Flagship Example*

**Solar System & Star Cluster Simulations** — Newton's law of universal gravitation with **7 integrators** (Euler, RK4, Verlet, Leapfrog, RK5, DP5, DP8). Symplectic integrators for long-term orbital stability, adaptive methods for highest accuracy. Self-contained physics engine (~870 lines).

<table>
<tr>
<td align="center" width="33%">

![Solar System](docs/images/readme/examples/00_N_body_gravity/_Example00_solar_system.png)

*Solar system orbital mechanics*

</td>
<td align="center" width="33%">

![Particle Sim](docs/images/readme/examples/00_N_body_gravity/_Example00_solar_system_particle_sim.png)

*Real-time particle visualization*

</td>
<td align="center" width="33%">

![Cluster Overview](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_overview.png)

*Star cluster collision overview*

</td>
</tr>
<tr>
<td align="center">

![Cluster Step 1](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_step%201.png)

*Cluster approach*

</td>
<td align="center">

![Cluster Step 3](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_step%203.png)

*Gravitational interaction*

</td>
<td align="center">

![Cluster Trajectories](docs/images/readme/examples/00_N_body_gravity/_Example00_star_clusters_trajectories_visualization.png)

*Full trajectory visualization*

</td>
</tr>
</table>

```cpp
NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config1_Solar_system();
NBodyGravitySimulator solver(config);

auto results_verlet = solver.SolveVerlet(0.01, 365000);       // 10 years, symplectic
auto results_dp8 = solver.SolveDP8(10.0, 1e-12, 0.01, 0.001); // 10 years, adaptive DP8
```

---

### 🎾 [Example 02: Double Pendulum](docs/examples/Example_02_double_pendulum.md) — *Chaos Theory*

Deterministic chaos in action — infinitely small differences in initial conditions lead to completely different outcomes.

<table>
<tr>
<td align="center" width="33%">

![Trajectory](docs/images/readme/examples/02_double_pendulum/01_trajectory_for%20both%20angles.png)

*Angle trajectories*

</td>
<td align="center" width="33%">

![Phase Space](docs/images/readme/examples/02_double_pendulum/02_phase_space.png)

*Phase space portrait*

</td>
<td align="center" width="33%">

![Butterfly Effect](docs/images/readme/examples/02_double_pendulum/03_butterfly_effect.png)

*Butterfly effect divergence*

</td>
</tr>
</table>

---

### 🏎️ [Example 03: Formula 1 G-Force Analysis](docs/examples/Example_03_F1_GForce_analysis.md) — *Parametric Curves*

Real F1 telemetry data (Silverstone, Monza) analyzed using MML's parametric curve and curvature calculations. Lateral G = v²κ/g, longitudinal G = (1/g)·dv/dt.

<table>
<tr>
<td align="center" width="33%">

![Track Path](docs/images/readme/examples/03_formula_1_sim/01_track_path.png)

*Track layout from telemetry*

</td>
<td align="center" width="33%">

![G-Forces](docs/images/readme/examples/03_formula_1_sim/02_g_forces.png)

*G-force profile around lap*

</td>
<td align="center" width="33%">

![Speed Profile](docs/images/readme/examples/03_formula_1_sim/03_speed_profile.png)

*Speed profile analysis*

</td>
</tr>
</table>

---

### 💥 [Example 04: 2D Collision Simulator](docs/examples/Example_04_collision_simulator_2d.md) — *Kinetic Theory*

**30,000+ particles** with exact elastic collision physics, spatial partitioning for O(N) performance, and multi-threaded execution. Watch shock waves propagate!

<table>
<tr>
<td align="center" width="33%">

![Shock 1](docs/images/readme/examples/04_collision_simulator_2d/01_shock_wave.png)

*Initial shock front*

</td>
<td align="center" width="33%">

![Shock 2](docs/images/readme/examples/04_collision_simulator_2d/02_shock_wave.png)

*Wave propagation*

</td>
<td align="center" width="33%">

![Shock 3](docs/images/readme/examples/04_collision_simulator_2d/03_shock_wave.png)

*Shock wave dispersion*

</td>
</tr>
</table>

---

### 📦 [Example 05: Rigid Body Collisions](docs/examples/Example_05_rigid_body.md) — *3D Dynamics*

Two parallelepipeds and a sphere in a cubic container with elastic collisions, full rotational dynamics using quaternions and inertia tensors.

<table>
<tr>
<td align="center" width="50%">

![Start](docs/images/readme/examples/05_rigid_body/01_rigid_body_start.png)

*Initial configuration*

</td>
<td align="center" width="50%">

![Collision](docs/images/readme/examples/05_rigid_body/02_rigid_body.png)

*Mid-collision dynamics*

</td>
</tr>
</table>

---

### More Examples

| # | Example | Description |
|---|---------|-------------|
| 01 | [**Projectile Launch**](docs/examples/Example_01_projectile_launch.md) | Ballistic trajectory with air resistance — drag models, vacuum vs air comparison |
| 06 | [**Lorentz Transformations**](docs/examples/Example_06_Lorentz_transformations.md) | Special relativity — time dilation, length contraction, Twin Paradox |

### 🚀 Try It Now

```bash
cmake -B build && cmake --build build

# Run the flagship N-Body simulation
./build/src/examples/Release/Example00_NBodyGravity      # Windows
./build/src/examples/Example00_NBodyGravity              # Linux
```

All examples produce visualization output viewable with the included Qt-based viewers.

---

## 📊 Visualization Suite


MML includes a powerful **cross-platform visualization suite** for functions, fields, curves, surfaces, and particle systems. All visualizers work on Windows (WPF), Linux (Qt), and macOS (Qt).

### Visualization Gallery

#### Windows — WPF Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![WPF Real](docs/images/readme/visualization_suite/win/wpf_real_func_multi_damped_oscillations.png)

</td>
<td align="center" width="25%">

**Real Functions (Lorentz)**

![WPF Lorentz](docs/images/readme/visualization_suite/win/wpf_real_func_multi_Lorentz.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![WPF Scalar 2D](docs/images/readme/visualization_suite/win/wpf_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![WPF Scalar 3D](docs/images/readme/visualization_suite/win/wpf_scalar_func_3d.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![WPF Curve 2D](docs/images/readme/visualization_suite/win/wpf_param_curve_2d_butterfly.png)

</td>
<td align="center">

**Parametric Curve 3D**

![WPF Curve 3D](docs/images/readme/visualization_suite/win/wpf_param_curve_3d.png)

</td>
<td align="center">

**Parametric Surface**

![WPF Surface](docs/images/readme/visualization_suite/win/wpf_param_surface.png)

</td>
<td align="center">

**Vector Field 3D**

![WPF VecField](docs/images/readme/visualization_suite/win/wpf_vector_field_3d.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![WPF Particle 2D](docs/images/readme/visualization_suite/win/wpf_particle_vis_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![WPF Particle 3D](docs/images/readme/visualization_suite/win/wpf_particle_vis_3d.png)

</td>
<td align="center">

**Rigid Body Simulation**

![WPF Rigid](docs/images/readme/visualization_suite/win/wpf_rigid_body.png)

</td>
<td align="center">

</td>
</tr>
</table>

#### Windows — Qt Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![Qt Real](docs/images/readme/visualization_suite/win/win_qt_real_func_multi.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![Qt Scalar 2D](docs/images/readme/visualization_suite/win/win_qt_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![Qt Scalar 3D](docs/images/readme/visualization_suite/win/win_qt_scalar_func_3d.png)

</td>
<td align="center" width="25%">

**Parametric Surface**

![Qt Surface](docs/images/readme/visualization_suite/win/win_qt_param_surface.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![Qt Curve 2D](docs/images/readme/visualization_suite/win/win_qt_param_curve_2d.png)

</td>
<td align="center">

**Parametric Curve 3D**

![Qt Curve 3D](docs/images/readme/visualization_suite/win/win_qt_param_curve_3d.png)

</td>
<td align="center">

**Vector Field 3D**

![Qt VecField](docs/images/readme/visualization_suite/win/win_qt_vector_field_3d.png)

</td>
<td align="center">

**Rigid Body Simulation**

![Qt Rigid](docs/images/readme/visualization_suite/win/win_qt_rigid_body.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![Qt Particle 2D](docs/images/readme/visualization_suite/win/win_qt_particle_vis_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![Qt Particle 3D](docs/images/readme/visualization_suite/win/win_qt_particle_vis_3d.png)

</td>
<td align="center">

</td>
<td align="center">

</td>
</tr>
</table>

#### Linux — Qt Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![Linux Real](docs/images/readme/visualization_suite/linux/linux_qt_real_func_multri.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![Linux Scalar 2D](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_2d.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![Linux Scalar 3D](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_3d.png)

</td>
<td align="center" width="25%">

**Parametric Surface**

![Linux Surface](docs/images/readme/visualization_suite/linux/linux_qt_param_surface.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![Linux Curve 2D](docs/images/readme/visualization_suite/linux/linux_qt_param_curve_2d.png)

</td>
<td align="center">

**Parametric Curve 3D**

![Linux Curve 3D](docs/images/readme/visualization_suite/linux/linux_qt_param_curve_3d.png)

</td>
<td align="center">

**Vector Field 3D**

![Linux VecField](docs/images/readme/visualization_suite/linux/linux_qt_vector_field_3d.png)

</td>
<td align="center">

**Rigid Body Simulation**

![Linux Rigid](docs/images/readme/visualization_suite/linux/linux_qt_rigid_body_vis.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![Linux Particle 2D](docs/images/readme/visualization_suite/linux/linux_qt_particle_vis_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![Linux Particle 3D](docs/images/readme/visualization_suite/linux/linux_qt_particle_vis_3d.png)

</td>
<td align="center">

**Scalar 2D (Dark Theme)**

![Linux Scalar Dark](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_2d_dark.png)

</td>
<td align="center">

**Parametric Surface (Dark)**

![Linux Surface Dark](docs/images/readme/visualization_suite/linux/linux_qt_param_surface_dark.png)

</td>
</tr>
</table>

#### macOS — Qt Visualizers

<table>
<tr>
<td align="center" width="25%">

**Real Functions**

![Mac Real](docs/images/readme/visualization_suite/mac/mac_qt_real_func_multi.png)

</td>
<td align="center" width="25%">

**Real Functions (Lorentz)**

![Mac Lorentz](docs/images/readme/visualization_suite/mac/mac_qt_real_func_multi_Lorentz_system.png)

</td>
<td align="center" width="25%">

**Scalar Function 2D**

![Mac Scalar 2D](docs/images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_monkey_saddle.png)

</td>
<td align="center" width="25%">

**Scalar Function 3D**

![Mac Scalar 3D](docs/images/readme/visualization_suite/mac/mac_qt_scalar_function_3d_gyroid.png)

</td>
</tr>
<tr>
<td align="center">

**Parametric Curve 2D**

![Mac Curve 2D](docs/images/readme/visualization_suite/mac/mac_qt_param_curve_2de_butterfly.png)

</td>
<td align="center">

**Parametric Curve 3D**

![Mac Curve 3D](docs/images/readme/visualization_suite/mac/mac_qt_param_curve_3d_trefoil.png)

</td>
<td align="center">

**Parametric Surface**

![Mac Surface](docs/images/readme/visualization_suite/mac/mac_qt_param_surface_klein.png)

</td>
<td align="center">

**Vector Field 3D**

![Mac VecField](docs/images/readme/visualization_suite/mac/mac_qt_vector_field_3d_gravity.png)

</td>
</tr>
<tr>
<td align="center">

**Particle Visualizer 2D**

![Mac Particle 2D](docs/images/readme/visualization_suite/mac/mac_qt_particle_visualizer_2d.png)

</td>
<td align="center">

**Particle Visualizer 3D**

![Mac Particle 3D](docs/images/readme/visualization_suite/mac/mac_qt_partice_visualizer_3d.png)

</td>
<td align="center">

**Scalar Function 2D (Ripple)**

![Mac Ripple](docs/images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_ripple.png)

</td>
<td align="center">

**Rigid Body Simulation**

![Mac Rigid](docs/images/readme/visualization_suite/mac/mac_qt_rigid_body.png)

</td>
</tr>
</table>

### Available Visualizers

All visualizers have **complete cross-platform support** with multiple backend options:
- **Windows**: WPF-based visualizers (primary), Qt and FLTK also available
- **Linux**: Qt-based visualizers (primary), FLTK for lightweight 2D
- **macOS**: Qt-based visualizers (primary), FLTK for lightweight 2D

| Visualizer | Purpose | Output | Platform |
|------------|---------|--------|----------|
| **RealFunctionVisualizer** | Plot 1D functions | 2D graphs | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **MultiRealFunctionVisualizer** | Compare multiple 1D functions | 2D overlaid graphs | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **ScalarFunction2DVisualizer** | 3D surface plots | 3D surfaces | Windows (WPF), Linux/macOS (Qt) |
| **ParametricCurve2DVisualizer** | 2D parametric curves | 2D curves | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **ParametricCurve3DVisualizer** | 3D parametric curves | 3D curves | Windows (WPF), Linux/macOS (Qt) |
| **VectorField2DVisualizer** | 2D vector field arrows | 2D arrow plots | Windows (WPF), Linux/macOS (Qt/FLTK) |
| **VectorField3DVisualizer** | 3D vector field arrows | 3D arrow plots | Windows (WPF), Linux/macOS (Qt) |
| **ParticleVisualizer2D** | Animated 2D particle systems | 2D animations with playback | Windows (WPF), Linux/macOS (Qt) |
| **ParticleVisualizer3D** | Animated 3D particle systems | 3D animations with playback | Windows (WPF), Linux/macOS (Qt) |


### Visualization Examples
- `src/visualization_examples/` — Ready-to-run demos for every visualizer type
- `src/book/chapters/Chapter_02_vizualization/` — Comprehensive examples from the companion book
- Examples of custom styling, multi-plot layouts, and animation controls

**Real Function Plotting — Lorenz System Time Series:**

```cpp
// Lorenz attractor — chaotic time evolution of x(t), y(t), z(t)
ODESystem lorenz_system(3, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
    const Real sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
});

Vector<Real> initial_state({ 1.0, 1.0, 1.0 });
CashKarpIntegrator solver(lorenz_system);
ODESystemSolution sol = solver.integrate(initial_state, 0.0, 50.0, 0.001, 1e-10, 0.001);

// Extract solution components as real functions via spline interpolation
Vector<Real> t_vals = sol.getTValues();
SplineInterpRealFunc x_t(t_vals, sol.getXValues(0));
SplineInterpRealFunc y_t(t_vals, sol.getXValues(1));
SplineInterpRealFunc z_t(t_vals, sol.getXValues(2));

std::vector<IRealFunction*> time_series = { &x_t, &y_t, &z_t };
Visualizer::VisualizeMultiRealFunction(
    time_series, "Lorenz System: Chaotic Time Evolution",
    { "x(t)", "y(t)", "z(t)" }, 0.0, 50.0, 1000, "lorenz_time_series.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Lorenz Win](docs/images/readme/visualization_suite/win/wpf_real_func_multi_Lorentz.png) | ![Lorenz Linux](docs/images/readme/visualization_suite/linux/linux_qt_real_func_multri.png) | ![Lorenz Mac](docs/images/readme/visualization_suite/mac/mac_qt_real_func_multi_Lorentz_system.png) |

**Scalar Function 2D — Ripple (2D Sinc / Sombrero):**

```cpp
// 2D Sinc function — concentric ripples with central peak
ScalarFunction<2> sinc_2d{ [](const VectorN<Real, 2>& v) {
    Real r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
    if (r < 1e-10) return Real(50.0);
    return 50.0 * std::sin(r) / r;
}};

Visualizer::VisualizeScalarFunc2DCartesian(
    sinc_2d, "2D Sinc (Sombrero): z = 50·sin(r)/r",
    -15.0, 15.0, 80, -15.0, 15.0, 80, "sinc_2d.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Ripple Win](docs/images/readme/visualization_suite/win/wpf_scalar_func_2d.png) | ![Ripple Linux](docs/images/readme/visualization_suite/linux/linux_qt_scalar_func_2d.png) | ![Ripple Mac](docs/images/readme/visualization_suite/mac/mac_qt_scalar_func_2d_ripple.png) |

**Parametric Curves — Trefoil Knot:**

```cpp
// Trefoil knot: x = sin(t) + 2sin(2t), y = cos(t) - 2cos(2t), z = -sin(3t)
auto trefoil = [](Real t) {
    Real x = 50.0 * (std::sin(t) + 2.0*std::sin(2.0*t));
    Real y = 50.0 * (std::cos(t) - 2.0*std::cos(2.0*t));
    Real z = 50.0 * (-std::sin(3.0*t));
    return VectorN<Real, 3>{x, y, z};
};

ParametricCurve<3> knot(trefoil);
Visualizer::VisualizeParamCurve3D(
    knot, "Trefoil Knot", 0.0, 2.0*Constants::PI, 500, "trefoil.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Trefoil Win](docs/images/readme/visualization_suite/win/wpf_param_curve_3d.png) | ![Trefoil Linux](docs/images/readme/visualization_suite/linux/linux_qt_param_curve_3d.png) | ![Trefoil Mac](docs/images/readme/visualization_suite/mac/mac_qt_param_curve_3d_trefoil.png) |

**Vector Fields — Two-Body Gravity:**

```cpp
// Gravity field of two masses
VectorFunction<3> gravity{[](const VectorN<Real, 3> &x) {
    const VectorN<Real, 3> x1{100, 0, 0}, x2{-100, 0, 0};
    const Real m1 = 1000, m2 = 1000, G = 10;
    return -G * m1 * (x - x1) / pow((x - x1).NormL2(), 3)
           -G * m2 * (x - x2) / pow((x - x2).NormL2(), 3);
}};

Visualizer::VisualizeVectorField3DCartesian(
    gravity, "Gravity Field",
    -200, 200, 15, -200, 200, 15, -200, 200, 15, "gravity.mml");
```

| Windows (WPF) | Linux (Qt) | macOS (Qt) |
|:-------------:|:----------:|:----------:|
| ![Field Win](docs/images/readme/visualization_suite/win/wpf_vector_field_3d.png) | ![Field Linux](docs/images/readme/visualization_suite/linux/linux_qt_vector_field_3d.png) | ![Field Mac](docs/images/readme/visualization_suite/mac/mac_qt_vector_field_3d_gravity.png) |

### Data Export

All visualizers use serialized data files that can also be loaded by external tools (Python, MATLAB, etc.):

```cpp
// Serialize function data
Serializer::SaveRealFunc(f, "Function", 0, 10, 100, "data.txt");

// Serialize ODE solution
Serializer::SaveODESolutionAsMultiFunc(
    solution, "ODE Solution", {"x", "v"}, "ode_data.txt");

// Serialize 2D vector field
Serializer::SaveVectorFunc2DCartesian(
    field, "Vector Field", -10, 10, 20, -10, 10, 20, "field_data.txt");
```


---

## 📝 Code Examples

### Vectors & Matrices

📄 *[View full source](src/readme_examples/readme02_vectors_matrices.cpp)*

```cpp
#include "MML.h"
using namespace MML;

// Real vectors and matrices
Vector<Real> vec1{1.5, -2.0, 0.5}, vec2{1.0, 1.0, -3.0};
Matrix<Real> mat_3x3{3, 3, {1.0,  2.0, -1.0,
                           -1.0,  5.0,  6.0,
                            3.0, -2.0,  1.0}};

// Vector and matrix arithmetic
Vector<Real> result = Real{2.0} * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
std::cout << "Result: " << result << std::endl;

// Complex vectors and matrices
VectorComplex vec_cmplx{{Complex(1,1), Complex(-1,2)}};
MatrixComplex mat_cmplx{2, 2, {Complex(0.5,1),  Complex(-1,2),
                               Complex(-1,-2), Complex(-2,2)}};

// Matrix properties
std::cout << "IsOrthogonal:  " << Utils::IsOrthogonal(mat_3x3) << std::endl;
std::cout << "IsHermitian:   " << Utils::IsHermitian(mat_cmplx) << std::endl;

/* Expected OUTPUT:
    Result: [   -3.137858162,     3.922322703,    -8.629109946]
    IsOrthogonal:  0
    IsHermitian:   1
*/
```

---

### Linear Systems & Eigenvalues

📄 *[View full source](src/readme_examples/readme03_linear_systems.cpp)*

```cpp
// Define a linear system Ax = b
Matrix<Real> A{5, 5, {0.2,  4.1, -2.1, -7.4,  1.6,
                      1.6,  1.5, -1.1,  0.7,  5.0,
                     -3.8, -8.0,  9.6, -5.4, -7.8,
                      4.6, -8.2,  8.4,  0.4,  8.0,
                     -2.6,  2.9,  0.1, -9.6, -2.7}};
Vector<Real> b{1.1, 4.7, 0.1, 9.3, 0.4};

// Solve using LU decomposition
LUSolver<Real> luSolver(A);
Vector<Real> x = luSolver.Solve(b);

std::cout << "Solution: " << x << std::endl;
std::cout << "Verification A*x: " << (A * x) << std::endl;
std::cout << "Residual norm: " << (A*x - b).NormL2() << std::endl;

// Compute eigenvalues and eigenvectors
auto eigenResult = EigenSolver::Solve(A);

std::cout << "Eigenvalues:" << std::endl;
for (const auto& ev : eigenResult.eigenvalues)
    std::cout << "  λ = " << ev.real << " + " << ev.imag << "i" << std::endl;

// QR decomposition
QRSolver<Real> qr(A);
Matrix<Real> Q = qr.GetQ();
Matrix<Real> R = qr.GetR();
std::cout << "QR valid: " << Matrix<Real>::AreEqual(A, Q * R, 1e-10) << std::endl;
```

---

### Polynomials & Algebra

📄 *[View full source](src/readme_examples/readme09_polynomials.cpp)*

```cpp
// Create polynomials: p(x) = 2x³ - 3x² + x - 5
PolynomReal p{-5, 1, -3, 2};  // Coefficients: [constant, x, x², x³]

// Evaluate at a point
Real val = p(2.0);  // p(2) = 2(8) - 3(4) + 2 - 5 = 16 - 12 + 2 - 5 = 1
std::cout << "p(2) = " << val << std::endl;

// Polynomial arithmetic
PolynomReal q{1, 2};           // q(x) = 2x + 1
PolynomReal sum = p + q;       // Addition
PolynomReal prod = p * q;      // Multiplication
std::cout << "Degree of p*q: " << prod.degree() << std::endl;

// Calculus on polynomials
PolynomReal dp = p.derivative();        // p'(x) = 6x² - 6x + 1
PolynomReal ip = p.integral();          // ∫p(x)dx (constant = 0)
std::cout << "p'(x) at x=1: " << dp(1.0) << std::endl;

// Solve quadratic: x² - 5x + 6 = 0 → x = 2, 3
Complex r1, r2;
int numReal = SolveQuadratic(1.0, -5.0, 6.0, r1, r2);
std::cout << "Roots: " << r1.real() << ", " << r2.real() << std::endl;

// Solve cubic: x³ - 6x² + 11x - 6 = 0 → x = 1, 2, 3
Complex c1, c2, c3;
SolveCubic(1.0, -6.0, 11.0, -6.0, c1, c2, c3);
std::cout << "Cubic roots: " << c1.real() << ", " << c2.real() 
          << ", " << c3.real() << std::endl;
```

---

### Defining Functions

📄 *[View full source](src/readme_examples/readme04_functions_interpolation.cpp)*

MML provides multiple ways to create function objects that can be derived, integrated, and analyzed:

```cpp
// =====================================================================
// CASE 1: From standalone function
// =====================================================================
double MyFunction(double x) { return sin(x) * (1.0 + 0.5*x*x); }
RealFunction f1(MyFunction);

// =====================================================================
// CASE 2: Direct lambda creation (most common and recommended)
// =====================================================================
RealFunction f2{[](Real x) { return sin(x) * (1.0 + 0.5*x*x); }};

// Different function types — all created with lambdas
ScalarFunction<3> scalar([](const VectorN<Real,3>& x) { 
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];   // scalar field: r²
});

VectorFunction<3> vector([](const VectorN<Real,3>& x) -> VectorN<Real,3> { 
    return {x[1]*x[2], x[0]*x[2], x[0]*x[1]};   // vector field
});

ParametricCurve<3> helix([](Real t) -> VectorN<Real,3> { 
    return {cos(t), sin(t), 0.2*t};             // 3D parametric curve
});

ParametricSurface<3> torus([](Real u, Real v) -> VectorN<Real,3> { 
    Real R = 3.0, r = 1.0;                      // parametric surface
    return {(R + r*cos(v))*cos(u), (R + r*cos(v))*sin(u), r*sin(v)}; 
});

// =====================================================================
// CASE 3a: From class with operator() — useful for stateful functions
// =====================================================================
class FunctionFromClassOperator {
    double _amplitude;
public:
    FunctionFromClassOperator(double amp) : _amplitude(amp) {}
    double operator()(double x) const { return _amplitude * std::sin(x); }
};

FunctionFromClassOperator obj(2.5);  // amplitude = 2.5
RealFunctionFromStdFunc f3(std::function<double(double)>{obj});

// =====================================================================
// CASE 3b: Class inheriting IRealFunction directly
// =====================================================================
class MyDerivedFunction : public IRealFunction {
    double _frequency;
public:
    MyDerivedFunction(double freq) : _frequency(freq) {}
    double operator()(double x) const override { return std::cos(_frequency * x); }
};
MyDerivedFunction f4(2.0);  // frequency = 2.0 → cos(2x)

// =====================================================================
// CASE 4: Wrapper for external/legacy class you can't modify
// =====================================================================
class ExternalComplexClass {
public:
    double ComputeValue(double x) const { return std::exp(-x*x); }
};

class ExternalClassWrapper : public IRealFunction {
    const ExternalComplexClass& _ref;
public:
    ExternalClassWrapper(const ExternalComplexClass& obj) : _ref(obj) {}
    double operator()(double x) const override { return _ref.ComputeValue(x); }
};

ExternalComplexClass externalObj;
ExternalClassWrapper f5(externalObj);  // Now usable with all MML algorithms!

// All these functions can now be derived, integrated, analyzed...
std::cout << "f2(1) = " << f2(1.0) << std::endl;           // 1.26147
std::cout << "f2'(1) = " << Derivation::NDer4(f2, 1.0) << std::endl;  // derivative
std::cout << "∫f2 = " << IntegrateSimpson(f2, 0, 1).value << std::endl;  // integral
```

---

### Interpolation

📄 *[View full source](src/readme_examples/readme04_functions_interpolation.cpp)*

Create smooth functions from discrete data points:

```cpp
// Data points (could be from experiment, simulation, file...)
Vector<Real> x_data{0, 2.5, 5.0, 7.5, 10.0};
Vector<Real> y_data{0, 1.5, 4.0, 9.5, 16.0};

// Three interpolation methods — each creates a callable function
LinearInterpRealFunc   linear_interp(x_data, y_data);      // Piecewise linear
SplineInterpRealFunc   spline_interp(x_data, y_data);      // Cubic spline (smooth)
PolynomInterpRealFunc  poly_interp(x_data, y_data, 3);     // Polynomial degree 3

// Evaluate anywhere in the interpolation range
Real x = 3.7;
std::cout << "Linear:     " << linear_interp(x) << std::endl;  // 2.70
std::cout << "Spline:     " << spline_interp(x) << std::endl;  // 2.41 (smoother)
std::cout << "Polynomial: " << poly_interp(x) << std::endl;    // 2.43

// Comparison at multiple points:
//     x     Linear   Spline   Polynom
//   -----  -------  -------  -------
//    1.0    0.600    0.576    0.480
//    3.0    2.000    1.842    1.760
//    5.0    4.000    4.000    4.000
//    7.0    7.400    6.954    7.680
//    9.0   12.200   12.654   14.240
```

---

### Numerical Derivatives

📄 *[View full source](src/readme_examples/readme05_numerical_calculus.cpp)*

```cpp
// Compare derivative orders on f(x) = sin(x) at x = 1.0
// Analytical derivative: f'(1) = cos(1) ≈ 0.5403023058681398

RealFunction f{[](Real x) { return std::sin(x); }};
Real x = 1.0;
Real analytical = std::cos(1.0);

Real der1 = Derivation::NDer1(f, x);   // O(h) - forward difference
Real der2 = Derivation::NDer2(f, x);   // O(h²) - central difference
Real der4 = Derivation::NDer4(f, x);   // O(h⁴) - 5-point stencil
Real der6 = Derivation::NDer6(f, x);   // O(h⁶) - 7-point stencil
Real der8 = Derivation::NDer8(f, x);   // O(h⁸) - 9-point stencil

std::cout << std::scientific << std::setprecision(10);
std::cout << "NDer1 error: " << std::abs(der1 - analytical) << std::endl;  // ~1e-8
std::cout << "NDer2 error: " << std::abs(der2 - analytical) << std::endl;  // ~1e-11
std::cout << "NDer8 error: " << std::abs(der8 - analytical) << std::endl;  // ~1e-14
```

---

### Numerical Integration

📄 *[View full source](src/readme_examples/readme05_numerical_calculus.cpp)*

```cpp
// Integrate f(x) = sin(x) from 0 to π — Analytical result: 2.0
RealFunction f{[](Real x) { return std::sin(x); }};
Real a = 0.0, b = Constants::PI;

auto trap = IntegrateTrap(f, a, b, 1e-8);
auto simp = IntegrateSimpson(f, a, b, 1e-8);
auto romb = IntegrateRomberg(f, a, b, 1e-10);

std::cout << "Trapezoid: " << trap.value << " (error: " << trap.error_estimate << ")\n";
std::cout << "Simpson:   " << simp.value << " (error: " << simp.error_estimate << ")\n";
std::cout << "Romberg:   " << romb.value << " (error: " << romb.error_estimate << ")\n";

// Monte Carlo for high-dimensional integrals
ScalarFunction<3> volume_func([](const VectorN<Real, 3>& v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
});
MonteCarloIntegrator<3> mc_integrator;
VectorN<Real, 3> lower{0, 0, 0}, upper{1, 1, 1};
auto mc_result = mc_integrator.integrate(volume_func, lower, upper, 
                                          MonteCarloConfig().samples(100000));
std::cout << "MC estimate: " << mc_result.value 
          << " +/- " << mc_result.error_estimate << std::endl;
```

---

### Root Finding Algorithms

📄 *[View full source](src/readme_examples/readme10_root_finding.cpp)*

```cpp
// Find root of f(x) = x³ - 2x - 5 (has root near x ≈ 2.0945)
RealFunction f{[](Real x) { return x*x*x - 2*x - 5; }};

// Compare different methods
Real root_bisect = RootFinding::FindRootBisection(f, 2.0, 3.0, 1e-12);
Real root_brent  = RootFinding::FindRootBrent(f, 2.0, 3.0, 1e-12);
Real root_newton = RootFinding::FindRootNewton(f, 2.0, 3.0, 1e-12);
Real root_ridder = RootFinding::FindRootRidders(f, 2.0, 3.0, 1e-12);

std::cout << std::setprecision(15);
std::cout << "Bisection: " << root_bisect << std::endl;
std::cout << "Brent:     " << root_brent << std::endl;
std::cout << "Newton:    " << root_newton << std::endl;
std::cout << "Ridders:   " << root_ridder << std::endl;

// Get detailed convergence info using config
RootFinding::RootFindingConfig config;
config.tolerance = 1e-14;
config.max_iterations = 100;

auto result = RootFinding::FindRootBrent(f, 2.0, 3.0, config);
std::cout << "Iterations: " << result.iterations_used << std::endl;
std::cout << "f(root) =   " << result.function_value << std::endl;

// Find multiple roots by bracketing
Vector<Real> brackets_lo, brackets_hi;
int numRoots = RootFinding::FindRootBrackets(f, -10.0, 10.0, 100, 
                                              brackets_lo, brackets_hi);
std::cout << "Found " << numRoots << " root bracket(s)" << std::endl;
```

---

### Field Operations

📄 *[View full source](src/readme_examples/readme06_field_operations.cpp)*

```cpp
// Scalar field: gravitational potential φ(x,y,z) = -1/r
ScalarFunction<3> potential([](const VectorN<Real, 3>& x) {
    return -1.0 / x.NormL2();
});
VectorN<Real, 3> pos{1.0, 2.0, 2.0};

// Gradient ∇φ (gives force direction)
auto grad = ScalarFieldOperations::GradientCart<3>(potential, pos);
std::cout << "Gradient at (1,2,2): " << grad << std::endl;

// Laplacian ∇²φ (zero outside mass for gravity!)
Real laplacian = ScalarFieldOperations::LaplacianCart<3>(potential, pos);
std::cout << "Laplacian at (1,2,2): " << laplacian << std::endl;

// Vector field: rotating velocity field v = (y, -x, z)
VectorFunction<3> velocity([](const VectorN<Real, 3>& x) -> VectorN<Real, 3> {
    return {x[1], -x[0], x[2]};
});

// Divergence ∇·v (compression/expansion rate)
Real div = VectorFieldOperations::DivCart<3>(velocity, pos);
std::cout << "Divergence at (1,2,2): " << div << std::endl;

// Curl ∇×v (rotation/vorticity)
auto curl = VectorFieldOperations::CurlCart(velocity, pos);
std::cout << "Curl at (1,2,2): " << curl << std::endl;
```

---

### Coordinate Transformations

📄 *[View full source](src/readme_examples/readme11_coord_transforms.cpp)*

```cpp
// Cartesian point
Vector3Cartesian cart_pos{1.0, 1.0, 1.0};

// Convert to Spherical (r, θ, φ) - Math/ISO convention
// θ = polar angle from z-axis, φ = azimuthal angle in xy-plane
Vector3Spherical sph_pos = CoordTransfCartToSpher.transf(cart_pos);
std::cout << "Spherical: r=" << sph_pos[0] << ", θ=" << sph_pos[1] 
          << ", φ=" << sph_pos[2] << std::endl;

// Convert to Cylindrical (r, φ, z)
Vector3Cylindrical cyl_pos = CoordTransfCartToCyl.transf(cart_pos);
std::cout << "Cylindrical: r=" << cyl_pos[0] << ", φ=" << cyl_pos[1] 
          << ", z=" << cyl_pos[2] << std::endl;

// Convert back to Cartesian
Vector3Cartesian back = CoordTransfSpherToCart.transf(sph_pos);
std::cout << "Back to Cartesian: " << back << std::endl;

// Field operations in spherical coordinates
// Inverse-square potential: φ = -1/r
ScalarFunction<3> pot_spher([](const VectorN<Real, 3>& x) { return -1.0/x[0]; });

// Use a point away from origin for gradient
Vector3Spherical test_sph{2.0, Constants::PI/4, Constants::PI/4};
auto grad_spher = ScalarFieldOperations::GradientSpher(pot_spher, test_sph);
std::cout << "Gradient of -1/r in spherical at r=2:" << std::endl;
std::cout << "  ∂φ/∂r = " << grad_spher[0] << " (analytical: 1/r² = 0.25)" << std::endl;

// Covariant vector transformation (transform gradient to Cartesian coords)
Vector3Cartesian test_cart = CoordTransfSpherToCart.transf(test_sph);
auto force_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, test_cart);
std::cout << "Force vector in Cartesian: " << force_cart << std::endl;
```

---

### Parametric Curves

📄 *[View full source](src/readme_examples/readme12_parametric_curves.cpp)*

```cpp
// Define a 3D helix: r(t) = (cos(t), sin(t), 0.2t)
Curves::CurveCartesian3D helix([](Real t) -> VectorN<Real, 3> { 
    return {cos(t), sin(t), 0.2*t}; 
});

Real t = Constants::PI / 4;

// Curve properties at parameter t
auto pos = helix(t);                    // Position on curve
auto tangent = helix.getTangent(t);     // Tangent vector dr/dt
auto unit_tan = helix.getTangentUnit(t);// Unit tangent T
auto normal = helix.getNormal(t);       // Normal vector (acceleration)
auto binormal = helix.getBinormal(t);   // Binormal B = T × N

std::cout << std::setprecision(6);
std::cout << "Position:       " << pos << std::endl;
std::cout << "Tangent dr/dt:  " << tangent << std::endl;
std::cout << "Unit tangent:   " << unit_tan << std::endl;
std::cout << "Normal (accel): " << normal << std::endl;
std::cout << "Binormal:       " << binormal << std::endl;

// Curvature κ (Frenet-Serret apparatus)
Real curvature = helix.getCurvature(t);
std::cout << "Curvature κ:    " << curvature << std::endl;

// Predefined curves
Curves::LemniscateCurve lemniscate;       // Figure-eight curve
Curves::ToroidalSpiralCurve torus(5, 2);  // Spiral on torus
Curves::Circle3DXZCurve circle(3.0);      // Circle in XZ plane
```

---

### Path & Line Integrals

📄 *[View full source](src/readme_examples/readme13_path_integrals.cpp)*

```cpp
// Line integral of vector field along a curve
// ∫ F·dr where F = (y, -x, 0) along unit circle

VectorFunction<3> F([](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
    return {p[1], -p[0], 0.0};
});

// Unit circle in XY plane: r(t) = (cos(t), sin(t), 0), t ∈ [0, 2π]
ParametricCurve<3> circle([](Real t) -> VectorN<Real, 3> {
    return {cos(t), sin(t), 0.0};
});

// Work integral: W = ∫₀^{2π} F(r(t)) · r'(t) dt
Real work = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, 1e-8);
std::cout << "Work integral ∮ F·dr:" << std::endl;
std::cout << "  Numerical:  " << work << std::endl;
std::cout << "  Analytical: -2π = " << -2*Constants::PI << std::endl;

// Arc length computation: ∫ ds
Real arc_length = PathIntegration::ParametricCurveLength<3>(circle, 0.0, 2*Constants::PI);
std::cout << "Arc length of unit circle:" << std::endl;
std::cout << "  Numerical:  " << arc_length << std::endl;
std::cout << "  Analytical: 2π = " << 2*Constants::PI << std::endl;
```

---

### Function Analysis

📄 *[View full source](src/readme_examples/readme14_function_analysis.cpp)*

```cpp
// Analyze function behavior over an interval
RealFunction f{[](Real x) { return x*x*x - 3*x + 1; }};

RealFunctionAnalyzer analyzer(f, "x³ - 3x + 1");
analyzer.PrintIntervalAnalysis(-3.0, 3.0, 100, 1e-6);
/*  Output:
    f(x) = x³ - 3x + 1 - Analysis in interval [-3.00, 3.00]:
      Defined    : yes
      Continuous : yes
      Monotonic  : no
      Min        : -0.999...
      Max        : 2.999...
*/

// Find roots using root bracket search + bisection
Vector<Real> xb1, xb2;
int numBrackets = RootFinding::FindRootBrackets(f, -3.0, 3.0, 100, xb1, xb2);
std::cout << "Found " << numBrackets << " root brackets:" << std::endl;
std::cout << std::setprecision(10);
for (int i = 0; i < numBrackets; i++) {
    Real root = RootFinding::FindRootBisection(f, xb1[i], xb2[i], 1e-10);
    std::cout << "  x = " << root << ", f(x) = " << f(root) << std::endl;
}

// Analyze a function with discontinuity
RealFunctionFromStdFunc step([](Real x) -> Real { 
    if (x < 0) return 0.0;
    else if (x > 0) return 1.0;
    else return 0.5;
});
RealFunctionAnalyzer step_analyzer(step, "step(x)");
step_analyzer.PrintIntervalAnalysis(-2.0, 2.0, 100, 1e-6);
// Detects discontinuity at x = 0
```

---

### Differential Equations

📄 *[View full source](src/readme_examples/readme07_ode_solvers.cpp)*

```cpp
// Simple harmonic oscillator: d²x/dt² = -ω²x
// As system: dx/dt = v, dv/dt = -ω²x
const Real omega = 2.0;

ODESystem system(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    const Real omega = 2.0;
    dydt[0] = y[1];                    // dx/dt = v
    dydt[1] = -omega*omega * y[0];     // dv/dt = -ω²x
});

Vector<Real> initial_cond{1.0, 0.0};  // x(0) = 1, v(0) = 0

// Solve with RK4 stepper
RungeKutta4_StepCalculator stepper;
ODESystemFixedStepSolver solver(system, stepper);
ODESystemSolution solution = solver.integrate(initial_cond, 0.0, 10.0, 1000);

// Access final values (index 999 for 1000 points)
int last = solution.size() - 1;
std::cout << "Final position: " << solution.getXValue(last, 0) << std::endl;
std::cout << "Final velocity: " << solution.getXValue(last, 1) << std::endl;

// Analytical solution at t=10: x(t) = cos(ωt), v(t) = -ω sin(ωt)
Real t_final = 10.0;
std::cout << "Analytical x:   " << std::cos(omega * t_final) << std::endl;
std::cout << "Analytical v:   " << -omega * std::sin(omega * t_final) << std::endl;
```

---

### 🦋 Dynamical Systems Analysis — *The Crown Jewel*

📄 *[View full source](src/readme_examples/readme15_dynamical_systems.cpp)*

```cpp
// THE LORENZ SYSTEM - The butterfly effect in action
// dx/dt = σ(y - x), dy/dt = x(ρ - z) - y, dz/dt = xy - βz
using namespace MML::Systems;

LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // Classic chaotic parameters

std::cout << "Dissipative: " << (lorenz.isDissipative() ? "yes" : "no") << std::endl;
std::cout << "Flow divergence: " << lorenz.getDivergence() << " (always negative -> attractor exists)" << std::endl;

// FIXED POINT ANALYSIS - Find equilibria and classify stability
std::vector<Vector<Real>> guesses = {
    Vector<Real>{0.0, 0.0, 0.0},      // Origin
    Vector<Real>{8.0, 8.0, 27.0},     // Near C+
    Vector<Real>{-8.0, -8.0, 27.0}    // Near C-
};

auto fixedPoints = FixedPointFinder::FindMultiple(lorenz, guesses);

std::cout << "Found " << fixedPoints.size() << " fixed points:" << std::endl;
for (const auto& fp : fixedPoints) {
    std::cout << "  Fixed Point: (" << fp.location[0] << ", " 
              << fp.location[1] << ", " << fp.location[2] << ")" << std::endl;
    std::cout << "    Type: " << ToString(fp.type) << std::endl;
    std::cout << "    Stable: " << (fp.isStable ? "yes" : "NO (unstable)") << std::endl;
    std::cout << "    Eigenvalues: ";
    for (const auto& ev : fp.eigenvalues)
        std::cout << "(" << ev.real() << "+" << ev.imag() << "i) ";
    std::cout << std::endl;
}

// LYAPUNOV EXPONENTS - Quantify chaos via sensitivity to initial conditions
Vector<Real> x0 = lorenz.getDefaultInitialCondition();
auto lyapResult = LyapunovAnalyzer::Compute(lorenz, x0, 
                                             500.0,   // Total integration time
                                             1.0,     // Orthonormalization interval
                                             0.01);   // Step size

std::cout << "Lyapunov spectrum: [" << lyapResult.exponents[0] << ", "
          << lyapResult.exponents[1] << ", " << lyapResult.exponents[2] << "]" << std::endl;
std::cout << "Maximum exponent: " << lyapResult.maxExponent;
if (lyapResult.maxExponent > 0) std::cout << " (POSITIVE -> CHAOS!)";
std::cout << std::endl;
std::cout << "Kaplan-Yorke dimension: " << lyapResult.kaplanYorkeDimension 
          << " (fractal attractor!)" << std::endl;
std::cout << "System is " << (lyapResult.isChaotic ? "CHAOTIC" : "regular") << std::endl;

// COMPARE CLASSIC CHAOTIC SYSTEMS
RosslerSystem rossler(0.2, 0.2, 5.7);       // Spiral chaos
VanDerPolSystem vanderpol(1.0);             // Limit cycle (NOT chaotic)
DoublePendulumSystem pendulum(1.0, 1.0, 9.81);  // Mechanical chaos

// Compute their Lyapunov spectra...
// Lorenz:       λ₁ ≈ +0.91  (CHAOTIC)
// Rössler:      λ₁ ≈ +0.07  (CHAOTIC)
// Van der Pol:  λ₁ ≈  0.00  (periodic)
// Double Pend:  λ₁ ≈ +1.5   (CHAOTIC)

// BIFURCATION ANALYSIS - Route to chaos via parameter sweeps
LorenzSystem lorenzSweep;
Vector<Real> sweepIC{1.0, 1.0, 1.0};

auto bifurcation = BifurcationAnalyzer::Sweep(
    lorenzSweep,
    1,              // Parameter index (rho)
    20.0, 30.0,     // Parameter range
    6,              // Number of parameter values
    sweepIC,
    2,              // Record z-component maxima
    50.0,           // Transient time
    20.0,           // Recording time
    0.01            // Step size
);

// Reveals: rho < 24: periodic → period-doubling cascade → rho ≈ 24.74: chaos onset!
```

---

## 🔬 Precision & Testing

MML takes numerical accuracy seriously. Every algorithm is rigorously validated against analytical solutions, and precision characteristics are documented.

### Test Suite Overview

**4,540 unit tests** across **111 test files**, organized by domain:

| Domain | Files | Key Validations |
|--------|-------|-----------------|
| **Linear Algebra** | 18 | LU/QR/SVD decomposition, eigensolvers, condition numbers up to 10¹⁵ |
| **Calculus** | 14 | Derivation (orders 1-8), integration (1D/2D/3D), Gauss-Kronrod |
| **ODE Solvers** | 6 | All steppers, event detection, stiff systems (λ = -10⁶) |
| **Geometry** | 22 | 2D/3D primitives, convex hull, Voronoi, KD-tree, triangulation |
| **Functions** | 12 | Real, scalar, vector, parametric curves/surfaces, interpolation |
| **Field Operations** | 8 | Gradient, divergence, curl in Cartesian/spherical/cylindrical |
| **Root Finding** | 4 | Bisection, Newton, Brent, Ridders with convergence verification |

### Predefined Test Beds

The `TestBeds::` namespace provides **battle-tested inputs** for algorithm validation:

| Test Bed | Contents | Purpose |
|----------|----------|---------|
| **Linear Systems** | Hilbert, Vandermonde, Kahan matrices | Stress-test solvers (κ = 10³ to 10¹⁵) |
| **ODE Systems** | Lorenz, Van der Pol, Robertson stiff problem | Validate adaptive stepping |
| **Integration** | Oscillatory, singular (1/√x), discontinuous | Test quadrature robustness |
| **Eigenvalues** | Symmetric, non-symmetric, repeated roots | Verify decomposition accuracy |
| **Functions** | sin, exp, polynomials with known derivatives | Precision benchmarks |
| **Parametric Curves** | Helix, Lemniscate, Torus Spiral | Differential geometry (analytical curvature) |

### Precision Benchmarks

**Numerical Derivation** — Error vs analytical derivative of sin(x) at x = 1.0:

| Method | Stencil | Error |
|--------|---------|-------|
| `NDer1` | 2-point forward | ~10⁻⁸ |
| `NDer2` | 3-point central | ~10⁻¹¹ |
| `NDer4` | 5-point | ~10⁻¹³ |
| `NDer6` | 7-point | ~10⁻¹⁴ |
| `NDer8` | 9-point | ~10⁻¹⁵ (machine ε!) |

**ODE Solver Energy Conservation** — Kepler orbit after 1000 periods:

| Solver | Relative Energy Drift |
|--------|----------------------|
| RK4 (fixed step) | ~10⁻⁶ |
| RKF45 (adaptive) | ~10⁻¹² |
| Dormand-Prince 8(5,3) | ~10⁻¹⁴ |

**Ill-Conditioned Problems** — MML handles edge cases:

```cpp
// Hilbert matrix — notoriously ill-conditioned (κ ≈ 10¹⁵ for n=10)
Matrix<Real> H = TestBeds::hilbert_10x10();   // condition number ~10¹³
LUSolver<Real> solver(H);
Vector<Real> x = solver.Solve(b);
// Achieves residual ||Ax - b|| / ||b|| < 10⁻³ despite extreme ill-conditioning

// Stiff ODE — Robertson chemical kinetics (λ ratio = 10⁸)
TestBeds::RobertsonStiffODE stiff;  // λ = {-0.04, -3×10⁴, -10⁸}
// Implicit solvers handle this; explicit solvers would need dt < 10⁻⁸
```

### Running the Test Suite

```bash
# Build and run all tests
cd build && ctest -j8 --output-on-failure

# Run specific test category
ctest -R "integration" -V

# Output: 4540 tests pass in ~30 seconds
```

📚 **Precision Analysis Reports:** [Overview](docs/testing_precision/README.md) • [Derivation](docs/testing_precision/DERIVATION_ANALYSIS.md) • [Integration](docs/testing_precision/INTEGRATION_ANALYSIS.md) • [ODE Solvers](docs/testing_precision/ODE_SOLVER_ANALYSIS.md)

📁 **Test Bed Documentation:** [Functions](docs/testbeds/Functions_testbed.md) • [ODE Systems](docs/testbeds/ODESystems_testbed.md) • [Linear Systems](docs/testbeds/LinAlgSystems_testbed.md) • [Curves & Surfaces](docs/testbeds/ParametricCurvesSurfaces_testbed.md)

---

## 📚 Documentation

| Resource | Description |
|----------|-------------|
| [Base Types](docs/base/README_Base.md) | Vectors, Matrices, Tensors, Functions |
| [Core Operations](docs/core/README_Core.md) | Derivation, Integration, Field Operations |
| [Algorithms](docs/algorithms/README_Algorithms.md) | ODE Solvers, Root Finding, Interpolation |
| [Systems](docs/systems/) | Dynamical Systems, Phase Space, Stability Analysis |
| [Tools](docs/tools/) | Visualization, Serialization, Console Output |
| [Examples](examples/) | Complete working examples |
| [Test Suites](docs/testbeds/) | Physics simulations and validation |
| [📖 Book References](references/book_references.md) | Foundational textbooks that informed MML |
| [📄 Paper References](references/paperes_references.md) | Academic papers behind the algorithms |

---

## 🛠️ Building & Testing

```bash
# Configure and build
cmake -B build
cmake --build build

# Run tests
cd build && ctest -j8 --output-on-failure

# Build examples
cmake --build build --target examples
```

---

## 📜 License

MML is released under the **MIT License** — free for personal, academic, and commercial use.

See [LICENSE.md](LICENSE.md) for details.

---

## ☕ Support MML

If MML has been useful to you, consider supporting its continued development:

<p align="center">
  <a href="https://github.com/sponsors/zvanjak">
    <img src="https://img.shields.io/badge/Sponsor-❤️-pink?logo=github&style=for-the-badge" alt="GitHub Sponsors">
  </a>

</p>

Your support helps maintain and improve MML! 🚀

---

<div align="center">

**Made with ❤️ for the C++ scientific computing community**

[⬆ Back to Top](#-mml--minimal-math-library)

</div>
