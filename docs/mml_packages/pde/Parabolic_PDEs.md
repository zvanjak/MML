# Parabolic PDEs: Heat and Diffusion Equations

> **Comprehensive Guide to Time-Dependent Diffusion Processes**

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Theory](#mathematical-theory)
3. [The Heat Equation](#the-heat-equation)
4. [Time Discretization Schemes](#time-discretization-schemes)
5. [Stability Analysis](#stability-analysis)
6. [1D Implementation](#1d-implementation)
7. [2D Implementation](#2d-implementation)
8. [3D Implementation](#3d-implementation)
9. [Observer Pattern](#observer-pattern)
10. [Verification and Validation](#verification-and-validation)
11. [Complete Examples](#complete-examples)
12. [Advanced Topics](#advanced-topics)

---

## Overview

Parabolic PDEs describe **time-dependent diffusion** processes where the solution evolves from an initial state according to diffusion laws. These equations model:

- **Heat conduction**: Temperature evolution in materials
- **Mass diffusion**: Concentration of chemicals, pollutants
- **Momentum diffusion**: Viscous fluid flow (Navier-Stokes)
- **Option pricing**: Black-Scholes equation in finance
- **Image processing**: Smoothing, denoising filters

### Key Characteristics

- **Initial-boundary value problems**: Need initial condition u(x, 0) and boundary conditions
- **Time evolution**: Solution progresses forward in time
- **Smoothing effect**: Initial discontinuities smooth out instantly
- **Irreversible**: Cannot "un-diffuse" (entropy increases)
- **Maximum principle**: No new extrema created in interior

---

## Mathematical Theory

### The Heat Equation

The general form in domain Ω over time interval [0, T]:

```
∂u/∂t = α∇²u + f(x, t)   in Ω × (0, T)
   u(x, 0) = u₀(x)        in Ω         (initial condition)
   u(x, t) = g(x, t)      on ∂Ω × (0, T) (boundary condition)
```

Where:
- `u(x, t)` is the unknown solution (e.g., temperature)
- `α` is the diffusion coefficient (thermal diffusivity, diffusion constant)
- `f(x, t)` is the source term (heat generation, chemical reaction)
- `u₀(x)` is the initial condition
- `g(x, t)` is the boundary condition

**In Different Dimensions:**

**1D:**
```
∂u/∂t = α ∂²u/∂x² + f(x, t),  x ∈ [a, b], t > 0
```

**2D:**
```
∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²) + f(x, y, t),  (x,y) ∈ Ω ⊂ ℝ², t > 0
```

**3D:**
```
∂u/∂t = α(∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²) + f(x, y, z, t),  (x,y,z) ∈ Ω ⊂ ℝ³, t > 0
```

### Physical Interpretation

**1. Heat Conduction (Fourier's Law)**

Temperature T satisfies:
```
ρc ∂T/∂t = k∇²T + q
```

Divide by ρc to get standard form:
```
∂T/∂t = α∇²T + f
```

where α = k/(ρc) is **thermal diffusivity**:
- k: thermal conductivity [W/(m·K)]
- ρ: density [kg/m³]
- c: specific heat capacity [J/(kg·K)]
- q: heat source [W/m³]

**Typical values:**
- Copper: α ≈ 1.1 × 10⁻⁴ m²/s
- Aluminum: α ≈ 9.7 × 10⁻⁵ m²/s
- Steel: α ≈ 1.4 × 10⁻⁵ m²/s
- Water: α ≈ 1.4 × 10⁻⁷ m²/s
- Air: α ≈ 2.2 × 10⁻⁵ m²/s

**2. Fick's Law (Mass Diffusion)**

Concentration c satisfies:
```
∂c/∂t = D∇²c + R
```

where D is diffusion coefficient and R is reaction rate.

**3. Black-Scholes Equation**

Option price V satisfies (transformed):
```
∂V/∂τ = (σ²/2)∂²V/∂S² + rS∂V/∂S - rV
```

Can be transformed to heat equation form.

### Well-Posedness

For heat equation with Dirichlet BC:
- **Existence**: Solution exists for u₀ ∈ L²(Ω), f ∈ L²(Ω × [0,T])
- **Uniqueness**: Solution is unique
- **Stability**: Solution depends continuously on initial/boundary data
- **Regularity**: Smooth for t > 0 even if u₀ is discontinuous

### Maximum Principle

**Parabolic maximum principle:**
```
min(u₀, min(g on ∂Ω × [0,T])) ≤ u(x, t) ≤ max(u₀, max(g on ∂Ω × [0,T]))
```

**Consequence:** No new maxima or minima created in interior!

---

## The Heat Equation

### Fundamental Solution (1D)

For infinite domain with point source at origin:
```
∂u/∂t = α ∂²u/∂x²,  x ∈ ℝ, t > 0
u(x, 0) = δ(x)
```

**Solution:**
```
u(x, t) = 1/√(4παt) exp(-x²/(4αt))
```

This is a **Gaussian** that spreads with time:
- Width: σ(t) = √(2αt)
- Height: 1/√(4παt) → 0 as t → ∞
- Conservation: ∫u dx = 1 for all t

### Separation of Variables (1D)

For homogeneous BC and no source:
```
∂u/∂t = α ∂²u/∂x²,  x ∈ [0, L], t > 0
u(0, t) = u(L, t) = 0
u(x, 0) = u₀(x)
```

**Solution:**
```
u(x, t) = Σₙ₌₁^∞ bₙ sin(nπx/L) exp(-α(nπ/L)²t)
```

where bₙ are Fourier sine coefficients of u₀.

**Key observation:** Higher modes decay faster! Mode n decays as exp(-n²t).

### Energy Decay

Define **energy**:
```
E(t) = ∫_Ω u²(x, t) dx
```

For homogeneous Dirichlet BC and no source:
```
dE/dt = -2α ∫_Ω |∇u|² dx ≤ 0
```

Energy **monotonically decreases** → equilibrium at u = 0.

---

## Time Discretization Schemes

### Method of Lines

**Strategy:**
1. Discretize space with finite differences → system of ODEs
2. Discretize time with time-stepping scheme

**Semi-discrete form (after spatial discretization):**
```
du/dt = Au + f(t)
```

where A is the spatial discretization matrix (e.g., Laplacian).

### Forward Euler (Explicit)

**Formula:**
```
(uⁿ⁺¹ - uⁿ)/Δt = Auⁿ + fⁿ
```

**Rearrange:**
```
uⁿ⁺¹ = uⁿ + Δt(Auⁿ + fⁿ)
```

**Properties:**
- ✅ **Explicit**: No linear system to solve (cheap per step)
- ✅ **Easy to implement**: One matrix-vector multiply
- ❌ **Conditionally stable**: Requires Δt ≤ Δt_CFL
- ❌ **First-order**: Local error O(Δt²), global error O(Δt)

**CFL Condition (1D):**
```
Δt ≤ h²/(2α)
```

**CFL Condition (2D):**
```
Δt ≤ h²/(4α)
```

**CFL Condition (3D):**
```
Δt ≤ h²/(6α)
```

**When to use:** Very small time steps acceptable, α very small, research/prototyping.

### Backward Euler (Implicit)

**Formula:**
```
(uⁿ⁺¹ - uⁿ)/Δt = Auⁿ⁺¹ + fⁿ⁺¹
```

**Rearrange:**
```
(I - ΔtA)uⁿ⁺¹ = uⁿ + Δt fⁿ⁺¹
```

**Properties:**
- ✅ **Unconditionally stable**: Any Δt allowed
- ✅ **First-order**: O(Δt) accuracy
- ✅ **Robust**: Good for stiff problems
- ❌ **Implicit**: Must solve linear system each step
- ❌ **Dissipative**: Overdamps oscillations

**System matrix:** `M = I - ΔtA`
- Still **SPD** (symmetric positive definite)
- Can use Conjugate Gradient
- Sparse structure preserved

**When to use:** Large time steps needed, long simulations, production code.

### Crank-Nicolson (Implicit)

**Formula (θ = 0.5):**
```
(uⁿ⁺¹ - uⁿ)/Δt = (A uⁿ⁺¹ + A uⁿ)/2 + (fⁿ⁺¹ + fⁿ)/2
```

**Rearrange:**
```
(I - Δt/2 A)uⁿ⁺¹ = (I + Δt/2 A)uⁿ + Δt/2(fⁿ⁺¹ + fⁿ)
```

**Properties:**
- ✅ **Unconditionally stable**: Any Δt allowed
- ✅ **Second-order**: O(Δt²) accuracy (best!)
- ✅ **No artificial dissipation**: Preserves oscillations
- ❌ **Implicit**: Linear system each step
- ❌ **Can oscillate**: For very large Δt

**When to use:** Production simulations, accuracy critical, moderate time steps.

### θ-Method (Generalized)

**Formula:**
```
(uⁿ⁺¹ - uⁿ)/Δt = θAuⁿ⁺¹ + (1-θ)Auⁿ + θfⁿ⁺¹ + (1-θ)fⁿ
```

**Special cases:**
- θ = 0: Forward Euler (explicit)
- θ = 0.5: Crank-Nicolson (second-order)
- θ = 1: Backward Euler (most dissipative)

**Stability:**
- θ < 0.5: Conditionally stable
- θ ≥ 0.5: Unconditionally stable

### Comparison Table

| Scheme | Order | Stability | Cost/step | Dissipation | Best Use |
|--------|-------|-----------|-----------|-------------|----------|
| Forward Euler | O(Δt) | CFL: Δt ≤ h²/(2αd) | Low | None | Research, small α |
| Backward Euler | O(Δt) | Unconditional | High | High | Robust, stiff |
| Crank-Nicolson | O(Δt²) | Unconditional | High | None | Production, accuracy |
| θ = 0.6 | O(Δt) | Unconditional | High | Some | Robust, less oscillation |

---

## Stability Analysis

### Von Neumann Analysis (1D)

**Discrete equation (Forward Euler):**
```
uⱼⁿ⁺¹ = uⱼⁿ + r(uⱼ₋₁ⁿ - 2uⱼⁿ + uⱼ₊₁ⁿ)
```

where r = αΔt/h² is the **mesh Fourier number**.

**Ansatz:** uⱼⁿ = ξⁿ exp(ikjh)

**Amplification factor:**
```
ξ = 1 + r(e⁻ⁱᵏʰ - 2 + eⁱᵏʰ)
  = 1 + r(-2 + 2cos(kh))
  = 1 - 4r sin²(kh/2)
```

**Stability condition:** |ξ| ≤ 1 for all k
```
-1 ≤ 1 - 4r sin²(kh/2) ≤ 1
```

Worst case: sin²(kh/2) = 1 (highest mode):
```
1 - 4r ≥ -1  ⟹  r ≤ 1/2  ⟹  Δt ≤ h²/(2α)
```

**This is the CFL condition for 1D heat equation!**

### CFL Derivation (General)

For d-dimensional heat equation with grid spacing h:
```
Δt ≤ h²/(2αd)
```

where d is spatial dimension:
- 1D: Δt ≤ h²/(2α)
- 2D: Δt ≤ h²/(4α)
- 3D: Δt ≤ h²/(6α)

**Physical interpretation:**
- Information travels distance ~√(αΔt) per time step
- Must be less than grid spacing: √(αΔt) < h
- Gives: Δt < h²/α (order of magnitude)

**Example (2D):**
- Grid: h = 0.01 m
- Material: steel, α = 1.4 × 10⁻⁵ m²/s
- CFL: Δt ≤ (0.01)² / (4 × 1.4 × 10⁻⁵) = 1.79 seconds
- To simulate 1 hour: 3600/1.79 ≈ 2000 steps!

This is why **implicit methods** are essential for parabolic PDEs.

### Implicit Scheme Stability

**Backward Euler amplification factor:**
```
ξ = 1 / (1 + 4r sin²(kh/2))
```

For r > 0: always |ξ| < 1 → **unconditionally stable**!

**Crank-Nicolson amplification factor:**
```
ξ = (1 - 2r sin²(kh/2)) / (1 + 2r sin²(kh/2))
```

Always |ξ| ≤ 1 → **unconditionally stable**, and ξ → -1 for large r (can oscillate).

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

// Homogeneous Dirichlet: u(0,t) = u(1,t) = 0
auto bc = homogeneousDirichlet1D<double>();

// Thermal diffusivity (e.g., steel)
double alpha = 1.4e-5;  // m²/s

// Create solver
HeatSolver1D<double> solver(grid, bc, alpha);

// Initial condition: u(x, 0) = sin(πx)
auto u0 = [](double x) { return std::sin(M_PI * x); };
solver.setInitialCondition(u0);

// Time parameters
double tFinal = 10.0;  // seconds
double dt = 0.1;       // seconds

// Solve with Crank-Nicolson
auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);

// Access final solution
auto u_final = solver.getSolution();
```

### Exponential Decay Example

**Problem:** `-∂u/∂t = α ∂²u/∂x²`, `u(x,0) = sin(πx)`, `u(0,t) = u(1,t) = 0`

**Analytical solution:** `u(x, t) = sin(πx) exp(-απ²t)`

```cpp
#include <MML.h>
#include <iostream>
#include <cmath>

int main() {
    using namespace MML;
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    auto bc = homogeneousDirichlet1D<double>();
    
    double alpha = 0.01;
    HeatSolver1D<double> solver(grid, bc, alpha);
    
    // Initial condition
    auto u0 = [](double x) { return std::sin(M_PI * x); };
    solver.setInitialCondition(u0);
    
    // Solve to t = 1.0
    double tFinal = 1.0;
    double dt = 0.01;
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    auto u_numerical = solver.getSolution();
    
    // Compare with analytical solution
    double decayFactor = std::exp(-alpha * M_PI * M_PI * tFinal);
    
    double maxError = 0.0;
    for (int i = 0; i < grid.numNodesX(); ++i) {
        double x = grid.nodeX(i);
        double u_exact = std::sin(M_PI * x) * decayFactor;
        double error = std::abs(u_numerical[i] - u_exact);
        maxError = std::max(maxError, error);
    }
    
    std::cout << "Time: " << tFinal << " s\n";
    std::cout << "Decay factor: " << decayFactor << "\n";
    std::cout << "Maximum error: " << maxError << "\n";
    std::cout << "Steps taken: " << result.numSteps << "\n";
    
    return 0;
}
```

### CFL Stability Demonstration

```cpp
#include <MML.h>
#include <iostream>

void testStability(double dt, const std::string& scheme_name) {
    using namespace MML;
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 50);  // h = 0.02
    auto bc = homogeneousDirichlet1D<double>();
    
    double alpha = 0.01;
    double h = 1.0 / 50;
    double dt_CFL = h * h / (2 * alpha);  // = 0.02
    
    std::cout << "\n" << scheme_name << " with dt = " << dt << "\n";
    std::cout << "CFL limit: " << dt_CFL << "\n";
    std::cout << "dt/dt_CFL = " << dt / dt_CFL << "\n";
    
    HeatSolver1D<double> solver(grid, bc, alpha);
    
    auto u0 = [](double x) { return std::sin(M_PI * x); };
    solver.setInitialCondition(u0);
    
    try {
        auto result = solver.solve(1.0, dt, TimeScheme::ForwardEuler);
        auto u = solver.getSolution();
        
        // Check for instability (unbounded growth)
        double maxVal = 0.0;
        for (double val : u) {
            maxVal = std::max(maxVal, std::abs(val));
        }
        
        if (maxVal > 10.0) {
            std::cout << "❌ UNSTABLE! Max value: " << maxVal << "\n";
        } else {
            std::cout << "✓ Stable. Max value: " << maxVal << "\n";
        }
    } catch (const std::exception& e) {
        std::cout << "❌ Failed: " << e.what() << "\n";
    }
}

int main() {
    // Below CFL: stable
    testStability(0.01, "Forward Euler");
    
    // Exactly at CFL: marginally stable
    testStability(0.02, "Forward Euler");
    
    // Above CFL: unstable!
    testStability(0.03, "Forward Euler");
    
    // Implicit: always stable
    testStability(0.1, "Backward Euler");
    
    return 0;
}
```

---

## 2D Implementation

### Room Cooling Simulation

**Problem:** Room with one open window, temperature equilibrating.

```cpp
#include <MML.h>
#include <iostream>
#include <fstream>

int main() {
    using namespace MML;
    
    // Room: 5m × 4m
    Rectangle<double> room(0.0, 5.0, 0.0, 4.0);
    Grid2D<double> grid(room, 50, 40);
    
    // Air thermal diffusivity
    double alpha = 2.2e-5;  // m²/s
    
    // Boundary conditions:
    // - Left wall (x=0): insulated, ∂T/∂n = 0
    // - Right wall (x=5): insulated, ∂T/∂n = 0
    // - Bottom (y=0): insulated, ∂T/∂n = 0
    // - Top (y=4): open window, T = 0°C
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));
    bc.setBottom(BoundaryCondition2D<double>::Neumann(0.0));
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(0.0));  // Open window!
    
    HeatSolver2D<double> solver(grid, bc, alpha);
    
    // Initial condition: room at 25°C
    auto T0 = [](double x, double y) { return 25.0; };
    solver.setInitialCondition(T0);
    
    // Simulate cooling for 100 seconds
    double tFinal = 100.0;
    double dt = 1.0;  // 1 second steps (well above CFL, but implicit!)
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    
    // Get final temperature
    auto T_final = solver.getSolution();
    
    // Save to CSV
    std::ofstream file("room_cooling.csv");
    file << "x,y,T\n";
    for (int j = 0; j < grid.numNodesY(); ++j) {
        for (int i = 0; i < grid.numNodesX(); ++i) {
            file << grid.nodeX(i) << ","
                 << grid.nodeY(j) << ","
                 << T_final[grid.index(i, j)] << "\n";
        }
    }
    
    std::cout << "Room cooling simulation complete!\n";
    std::cout << "Steps: " << result.numSteps << "\n";
    std::cout << "Time: " << result.tFinal << " s\n";
    
    // Temperature at center after 100s
    int i_mid = grid.numNodesX() / 2;
    int j_mid = grid.numNodesY() / 2;
    std::cout << "Center temperature: " 
              << T_final[grid.index(i_mid, j_mid)] << "°C\n";
    
    return 0;
}
```

### 2D with Source Term

**Problem:** Plate with internal heat generation.

```cpp
// Same setup...

// Heat source: q = 100 W/m³ in center region
auto q = [](double x, double y, double t) {
    double dx = x - 2.5;
    double dy = y - 2.0;
    double r2 = dx * dx + dy * dy;
    return (r2 < 0.25) ? 100.0 : 0.0;  // Source in circle of radius 0.5m
};

// Note: Source term support (setSourceTerm) not yet implemented in HeatSolver
// This example shows the intended API for future implementation
HeatSolver2D<double> solver(grid, bc, alpha);
solver.setInitialCondition(T0);
// solver.setSourceTerm(q);  // TODO: Not yet implemented

auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
```

### Time Scheme Comparison

```cpp
void compareSchemes() {
    using namespace MML;
    
    Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 50, 50);
    auto bc = homogeneousDirichlet2D<double>();
    double alpha = 0.01;
    
    auto u0 = [](double x, double y) {
        return std::sin(M_PI * x) * std::sin(M_PI * y);
    };
    
    // Analytical decay
    double t = 1.0;
    double decay = std::exp(-2 * alpha * M_PI * M_PI * t);
    
    std::cout << "Analytical decay: " << decay << "\n\n";
    
    // Test each scheme
    std::vector<std::pair<TimeScheme, std::string>> schemes = {
        {TimeScheme::ForwardEuler, "Forward Euler"},
        {TimeScheme::BackwardEuler, "Backward Euler"},
        {TimeScheme::CrankNicolson, "Crank-Nicolson"}
    };
    
    for (auto [scheme, name] : schemes) {
        HeatSolver2D<double> solver(grid, bc, alpha);
        solver.setInitialCondition(u0);
        
        double dt = (scheme == TimeScheme::ForwardEuler) ? 0.001 : 0.01;
        auto result = solver.solve(t, dt, scheme);
        auto u = solver.getSolution();
        
        double numerical_decay = u[grid.index(25, 25)];  // Center value
        double error = std::abs(numerical_decay - decay);
        
        std::cout << name << " (dt=" << dt << "):\n";
        std::cout << "  Steps: " << result.numSteps << "\n";
        std::cout << "  Decay: " << numerical_decay << "\n";
        std::cout << "  Error: " << error << "\n\n";
    }
}
```

---

## 3D Implementation

### 3D Heat Diffusion

```cpp
#include <MML.h>
using namespace MML;

int main() {
    // Domain: 1m × 1m × 1m cube
    Box<double> cube(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    
    // 30×30×30 grid (27,000 unknowns)
    Grid3D<double> grid(cube, 30, 30, 30);
    
    // All faces at 0°C
    auto bc = homogeneousDirichlet3D<double>();
    
    // Thermal diffusivity
    double alpha = 1.0e-5;
    
    HeatSolver3D<double> solver(grid, bc, alpha);
    
    // Initial: hot center
    auto T0 = [](double x, double y, double z) {
        double r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
        return 100.0 * std::exp(-r2 / 0.01);
    };
    solver.setInitialCondition(T0);
    
    // Solve
    double tFinal = 10.0;
    double dt = 0.1;
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    
    std::cout << "3D diffusion complete\n";
    std::cout << "Grid: 30×30×30 = 27,000 nodes\n";
    std::cout << "Steps: " << result.numSteps << "\n";
    
    return 0;
}
```

---

## Observer Pattern

### Monitoring Temperature Evolution

The observer pattern allows you to monitor the solution at each time step without storing all intermediate states.

```cpp
#include <MML.h>
#include <iostream>
#include <fstream>

class TemperatureMonitor {
private:
    std::ofstream file;
    const Grid1D<double>& grid;
    int sampleIndex;  // Which point to monitor

public:
    TemperatureMonitor(const std::string& filename, 
                      const Grid1D<double>& g, 
                      int idx)
        : file(filename), grid(g), sampleIndex(idx) {
        file << "time,temperature\n";
    }
    
    void operator()(double t, const std::vector<double>& u) {
        file << t << "," << u[sampleIndex] << "\n";
    }
};

int main() {
    using namespace MML;
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    auto bc = homogeneousDirichlet1D<double>();
    
    double alpha = 0.01;
    HeatSolver1D<double> solver(grid, bc, alpha);
    
    auto u0 = [](double x) { return std::sin(M_PI * x); };
    solver.setInitialCondition(u0);
    
    // Create observer to monitor center point
    int centerIdx = grid.numNodesX() / 2;
    TemperatureMonitor monitor("temperature_history.csv", grid, centerIdx);
    
    // Set observer
    solver.setObserver([&monitor](double t, const std::vector<double>& u) {
        monitor(t, u);
    });
    
    // Solve
    solver.solve(10.0, 0.1, TimeScheme::CrankNicolson);
    
    std::cout << "Temperature history saved to temperature_history.csv\n";
    
    return 0;
}
```

### Energy Tracking Observer

```cpp
class EnergyTracker {
private:
    std::ofstream file;
    const Grid2D<double>& grid;

public:
    EnergyTracker(const std::string& filename, const Grid2D<double>& g)
        : file(filename), grid(g) {
        file << "time,energy,max_temp,min_temp\n";
    }
    
    void operator()(double t, const std::vector<double>& u) {
        double energy = 0.0;
        double maxT = u[0];
        double minT = u[0];
        
        double dx = grid.dx();
        double dy = grid.dy();
        double dA = dx * dy;
        
        for (double T : u) {
            energy += T * T * dA;  // ∫T² dA
            maxT = std::max(maxT, T);
            minT = std::min(minT, T);
        }
        
        file << t << "," << energy << "," << maxT << "," << minT << "\n";
    }
};

// Usage
EnergyTracker tracker("energy_history.csv", grid);
solver.setObserver([&tracker](double t, const std::vector<double>& u) {
    tracker(t, u);
});
```

---

## Verification and Validation

### Method of Manufactured Solutions (1D)

**Choose:** `u(x, t) = exp(-t) sin(πx)`

**Compute source:** 
```
f = ∂u/∂t - α∂²u/∂x²
  = -exp(-t) sin(πx) - α(-π²) exp(-t) sin(πx)
  = exp(-t) sin(πx)(-1 + απ²)
```

```cpp
auto u_exact = [](double x, double t) {
    return std::exp(-t) * std::sin(M_PI * x);
};

auto f = [alpha = 0.01](double x, double t) {
    return std::exp(-t) * std::sin(M_PI * x) * (-1.0 + alpha * M_PI * M_PI);
};

// Test
HeatSolver1D<double> solver(grid, bc, alpha);
solver.setInitialCondition([](double x) { return std::sin(M_PI * x); });
// Note: Source term support not yet implemented
// solver.setSourceTerm(f);  // TODO: Not yet implemented

auto result = solver.solve(1.0, 0.01, TimeScheme::CrankNicolson);
auto u_numerical = solver.getSolution();

// Verify at t = 1.0
double maxError = 0.0;
for (int i = 0; i < grid.numNodesX(); ++i) {
    double x = grid.nodeX(i);
    double error = std::abs(u_numerical[i] - u_exact(x, 1.0));
    maxError = std::max(maxError, error);
}

std::cout << "MMS verification error: " << maxError << "\n";
```

### Convergence in Time

```cpp
void timeConvergenceStudy() {
    using namespace MML;
    
    // Fix spatial resolution
    Grid1D<double> grid(Interval<double>(0, 1), 100);
    auto bc = homogeneousDirichlet1D<double>();
    double alpha = 0.01;
    
    auto u0 = [](double x) { return std::sin(M_PI * x); };
    double tFinal = 1.0;
    
    // Reference solution (very fine dt)
    HeatSolver1D<double> ref_solver(grid, bc, alpha);
    ref_solver.setInitialCondition(u0);
    auto ref_result = ref_solver.solve(tFinal, 0.0001, TimeScheme::CrankNicolson);
    auto u_ref = ref_solver.getSolution();
    
    // Test various dt
    std::vector<double> dt_values = {0.1, 0.05, 0.025, 0.0125};
    std::vector<double> errors;
    
    for (double dt : dt_values) {
        HeatSolver1D<double> solver(grid, bc, alpha);
        solver.setInitialCondition(u0);
        auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
        auto u = solver.getSolution();
        
        double maxError = 0.0;
        for (size_t i = 0; i < u.size(); ++i) {
            maxError = std::max(maxError, std::abs(u[i] - u_ref[i]));
        }
        errors.push_back(maxError);
        
        std::cout << "dt = " << dt << ", error = " << maxError << "\n";
    }
    
    // Check convergence rate (should be ~2 for CN)
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        std::cout << "Convergence rate: " << rate << " (expect ~2.0)\n";
    }
}
```

---

## Complete Examples

### Example 1: 1D Rod Cooling

**Problem:** Hot rod cooling in ambient temperature.

```cpp
#include <MML.h>
#include <iostream>
#include <fstream>

int main() {
    using namespace MML;
    
    // Steel rod: 1 meter long
    Interval<double> rod(0.0, 1.0);
    Grid1D<double> grid(rod, 100);
    
    // Fixed temperatures at ends: 0°C
    auto bc = homogeneousDirichlet1D<double>();
    
    // Steel thermal diffusivity
    double alpha = 1.4e-5;  // m²/s
    
    HeatSolver1D<double> solver(grid, bc, alpha);
    
    // Initial: hot in center
    auto T0 = [](double x) {
        return 100.0 * std::exp(-100.0 * (x - 0.5) * (x - 0.5));
    };
    solver.setInitialCondition(T0);
    
    // Track center temperature
    std::ofstream file("rod_cooling.csv");
    file << "time,T_center\n";
    
    int centerIdx = grid.numNodesX() / 2;
    solver.setObserver([&file, centerIdx](double t, const std::vector<double>& u) {
        file << t << "," << u[centerIdx] << "\n";
    });
    
    // Simulate 1000 seconds
    solver.solve(1000.0, 10.0, TimeScheme::CrankNicolson);
    
    std::cout << "Rod cooling simulation complete!\n";
    
    return 0;
}
```

### Example 2: 2D Diffusion with Reaction

**Problem:** Chemical reaction-diffusion.

```cpp
#include <MML.h>
#include <iostream>

int main() {
    using namespace MML;
    
    Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 80, 80);
    
    // Zero concentration on boundary
    auto bc = homogeneousDirichlet2D<double>();
    
    // Diffusion coefficient
    double D = 0.001;  // m²/s
    
    HeatSolver2D<double> solver(grid, bc, D);
    
    // Initial concentration: spot at center
    auto c0 = [](double x, double y) {
        double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
        return (r2 < 0.01) ? 1.0 : 0.0;
    };
    solver.setInitialCondition(c0);
    
    // Reaction term: -kc (first-order decay)
    double k = 0.1;  // 1/s
    auto reaction = [k](double x, double y, double t) {
        // This would need access to current concentration
        // In practice, use nonlinear solver or operator splitting
        return 0.0;  // Simplified: no reaction for linear solver
    };
    
    // Solve
    solver.solve(10.0, 0.1, TimeScheme::CrankNicolson);
    
    auto c_final = solver.getSolution();
    
    // Compute total mass
    double mass = 0.0;
    for (double c : c_final) {
        mass += c * grid.dx() * grid.dy();
    }
    
    std::cout << "Final total mass: " << mass << "\n";
    
    return 0;
}
```

### Example 3: Comparison of All Schemes

```cpp
void compareAllSchemes() {
    using namespace MML;
    
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    auto bc = homogeneousDirichlet1D<double>();
    double alpha = 0.01;
    
    auto u0 = [](double x) { return std::sin(M_PI * x); };
    
    double tFinal = 1.0;
    double decay_exact = std::exp(-alpha * M_PI * M_PI * tFinal);
    
    std::cout << "Exact decay factor: " << decay_exact << "\n\n";
    
    // Forward Euler (small dt required)
    {
        HeatSolver1D<double> solver(grid, bc, alpha);
        solver.setInitialCondition(u0);
        
        double h = 1.0 / 100;
        double dt_CFL = 0.9 * h * h / (2 * alpha);  // Safety factor 0.9
        
        auto result = solver.solve(tFinal, dt_CFL, TimeScheme::ForwardEuler);
        auto u = solver.getSolution();
        
        std::cout << "Forward Euler:\n";
        std::cout << "  dt = " << dt_CFL << " (CFL-limited)\n";
        std::cout << "  Steps: " << result.numSteps << "\n";
        std::cout << "  Decay: " << u[50] << "\n";
        std::cout << "  Error: " << std::abs(u[50] - decay_exact) << "\n\n";
    }
    
    // Backward Euler (large dt allowed)
    {
        HeatSolver1D<double> solver(grid, bc, alpha);
        solver.setInitialCondition(u0);
        
        double dt = 0.1;  // 100x larger than FE!
        auto result = solver.solve(tFinal, dt, TimeScheme::BackwardEuler);
        auto u = solver.getSolution();
        
        std::cout << "Backward Euler:\n";
        std::cout << "  dt = " << dt << "\n";
        std::cout << "  Steps: " << result.numSteps << "\n";
        std::cout << "  Decay: " << u[50] << "\n";
        std::cout << "  Error: " << std::abs(u[50] - decay_exact) << "\n\n";
    }
    
    // Crank-Nicolson (best accuracy)
    {
        HeatSolver1D<double> solver(grid, bc, alpha);
        solver.setInitialCondition(u0);
        
        double dt = 0.1;
        auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
        auto u = solver.getSolution();
        
        std::cout << "Crank-Nicolson:\n";
        std::cout << "  dt = " << dt << "\n";
        std::cout << "  Steps: " << result.numSteps << "\n";
        std::cout << "  Decay: " << u[50] << "\n";
        std::cout << "  Error: " << std::abs(u[50] - decay_exact) << "\n\n";
    }
}
```

---

## Advanced Topics

### Adaptive Time Stepping

**Strategy:** Adjust dt based on local truncation error estimate.

**Error estimate (CN vs BE):**
```
error ≈ ||u_CN - u_BE|| / 2
```

**Algorithm:**
```
1. Take step with CN and BE
2. Estimate error
3. If error > tolerance: reject, reduce dt
4. If error < tolerance/10: accept, increase dt
5. Otherwise: accept, keep dt
```

**Implementation:**
```cpp
double adaptive_dt = dt_initial;
double t = 0.0;

while (t < tFinal) {
    // Take step with both schemes
    auto u_CN = solver.step(adaptive_dt, TimeScheme::CrankNicolson);
    auto u_BE = solver.step(adaptive_dt, TimeScheme::BackwardEuler);
    
    // Estimate error
    double error = computeError(u_CN, u_BE);
    
    if (error > tolerance) {
        // Reject step, reduce dt
        adaptive_dt *= 0.5;
    } else {
        // Accept step
        t += adaptive_dt;
        
        if (error < tolerance / 10) {
            // Increase dt
            adaptive_dt = std::min(adaptive_dt * 1.5, dt_max);
        }
    }
}
```

### Operator Splitting

For reaction-diffusion: `∂u/∂t = D∇²u + R(u)`

**Strang splitting (2nd order):**
```
1. React for Δt/2:     u* = u^n + (Δt/2) R(u^n)
2. Diffuse for Δt:     u** solves ∂u/∂t = D∇²u from u* to u**
3. React for Δt/2:     u^{n+1} = u** + (Δt/2) R(u**)
```

This allows using specialized solvers for each term.

### Nonlinear Heat Equation

**Stefan problem (phase change):**
```
∂u/∂t = ∇·(k(u)∇u) + f
```

where k(u) depends on temperature (e.g., accounts for latent heat).

**Solution:** Newton iteration or lagged nonlinearity
```
u^{n+1,m+1} = u^{n+1,m} + δu
where δu solves:
(I/Δt - ∇·(k(u^{n+1,m})∇)) δu = residual
```

### Preconditioning for Implicit Methods

Matrix to invert: `M = I - ΔtA`

**Preconditioners:**
- **Jacobi**: P = diag(M)
- **SSOR**: Symmetric Successive Over-Relaxation
- **ILU**: Incomplete LU factorization

**Impact:** 100-500 iterations → 10-50 iterations

### Multigrid Methods

For very fine grids (100³ and larger), multigrid achieves O(n) complexity:

**Idea:**
1. Smooth error on fine grid
2. Restrict residual to coarse grid
3. Solve coarse grid problem
4. Interpolate correction to fine grid
5. Smooth again

Can solve 1M unknowns in seconds!

---

## Summary

### Key Takeaways

✅ **Parabolic PDEs** describe time-dependent diffusion (heat, mass, momentum)

✅ **Heat equation**: `∂u/∂t = α∇²u + f` with initial and boundary conditions

✅ **Explicit schemes** (Forward Euler): CFL-limited, cheap per step

✅ **Implicit schemes** (Backward Euler, Crank-Nicolson): Unconditionally stable

✅ **Crank-Nicolson**: O(Δt²) accuracy, recommended for production

✅ **Observer pattern**: Monitor solution evolution without storing everything

✅ **Verification**: Manufactured solutions, analytical comparisons, convergence studies

### Quick Reference

| Scheme | Order | Stability | dt Restriction | Best Use |
|--------|-------|-----------|----------------|----------|
| Forward Euler | O(Δt) | Conditional | Δt ≤ h²/(2αd) | Research, prototypes |
| Backward Euler | O(Δt) | Unconditional | None | Robust, stiff problems |
| Crank-Nicolson | O(Δt²) | Unconditional | None | **Production** |
| θ-method (θ>0.5) | O(Δt) | Unconditional | None | Damped oscillations |

### CFL Conditions

- **1D**: Δt ≤ h²/(2α)
- **2D**: Δt ≤ h²/(4α)
- **3D**: Δt ≤ h²/(6α)

**Only applies to explicit schemes!**

### Related Documentation

- [Quick Start Guide](Quick_Start_Guide.md) - Get started quickly
- [Elliptic PDEs](Elliptic_PDEs.md) - Poisson and Laplace equations
- [Sparse Matrices](Sparse_Matrices.md) - Understanding CSR/COO formats
- [Grid Infrastructure](Grid_Infrastructure.md) - Grid1D/2D/3D details
- [Boundary Conditions](Boundary_Conditions.md) - BC types and API
- [Iterative Solvers](Iterative_Solvers.md) - CG for implicit systems
- [Examples Gallery](Examples_Gallery.md) - More complete examples
- [Troubleshooting](Troubleshooting.md) - Common issues

---

*Parabolic PDEs documentation - Part of MinimalMathLibrary PDE Solver Module*
