# Hyperbolic PDEs: Wave and Advection Equations

> **Comprehensive Guide to Wave Propagation and Transport Problems**

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Theory](#mathematical-theory)
3. [The Wave Equation](#the-wave-equation)
4. [The Advection Equation](#the-advection-equation)
5. [Time Discretization Schemes](#time-discretization-schemes)
6. [CFL Stability Condition](#cfl-stability-condition)
7. [1D Wave Solver](#1d-wave-solver)
8. [2D Wave Solver](#2d-wave-solver)
9. [1D Advection Solver](#1d-advection-solver)
10. [2D Advection Solver](#2d-advection-solver)
11. [Numerical Issues](#numerical-issues)
12. [Complete Examples](#complete-examples)
13. [Advanced Topics](#advanced-topics)

---

## Overview

Hyperbolic PDEs describe **wave propagation** and **transport** phenomena where information travels at finite speeds along characteristic curves. These equations model:

- **Mechanical waves**: Vibrating strings, membranes, sound waves
- **Electromagnetic waves**: Light, radio waves, radar
- **Fluid dynamics**: Shallow water waves, acoustics
- **Traffic flow**: Vehicle density waves
- **Chemical transport**: Pollutant advection, tracer transport

### Key Characteristics

- **Information propagates at finite speed**: The wave speed c determines how fast disturbances travel
- **Initial value problems**: Solution depends on initial conditions, not just boundaries
- **No smoothing effect**: Unlike parabolic PDEs, sharp features remain sharp (in exact solution)
- **Reversible**: Solutions can propagate both forward and backward in time
- **Conservation properties**: Energy (wave equation) or mass (advection) should be conserved

### Comparison with Other PDE Types

| Property | Elliptic | Parabolic | **Hyperbolic** |
|----------|----------|-----------|----------------|
| Time dependence | None (steady) | First-order | Second-order (wave) or First-order (advection) |
| Physical domain | Boundary-driven | Diffusion | Wave propagation |
| Information flow | Global | Forward in time | Along characteristics |
| Smoothing | Harmonic | Strong smoothing | **None** |
| Stability challenge | SPD system | Time step limit | CFL condition |

---

## Mathematical Theory

### The Wave Equation

The general form of the wave equation in domain Ω over time interval [0, T]:

```
∂²u/∂t² = c²∇²u + f(x, t)   in Ω × (0, T)
    u(x, 0) = u₀(x)          in Ω         (initial displacement)
∂u/∂t(x, 0) = v₀(x)          in Ω         (initial velocity)
```

Where:
- `u(x, t)` is the displacement (e.g., string deflection, pressure)
- `c` is the wave speed (m/s)
- `f(x, t)` is the source/forcing term
- `u₀(x)` is the initial displacement
- `v₀(x)` is the initial velocity

**In Different Dimensions:**

**1D:**
```
∂²u/∂t² = c² ∂²u/∂x² + f(x, t),   x ∈ [a, b], t > 0
```

**2D:**
```
∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y²) + f(x, y, t),   (x,y) ∈ Ω ⊂ ℝ², t > 0
```

**3D:**
```
∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²) + f(x, y, z, t)
```

### The Advection Equation

The linear advection (transport) equation:

```
∂u/∂t + c·∂u/∂x = 0,   x ∈ [a, b], t > 0
   u(x, 0) = u₀(x)     (initial condition)
```

**Multi-dimensional form:**
```
∂u/∂t + a·∇u = 0
```

where `a = (a_x, a_y, ...)` is the velocity field.

**Exact Solution:** The advection equation has a simple exact solution:
```
u(x, t) = u₀(x - ct)
```
The initial profile simply translates at speed c without changing shape.

### Physical Interpretations

**Wave Equation Applications:**

1. **Vibrating String** (1D): Guitar, violin, piano
   - c² = T/ρ where T is tension, ρ is linear density
   
2. **Vibrating Membrane** (2D): Drum head, speaker cone
   - c² = T/σ where T is surface tension, σ is area density
   
3. **Acoustic Waves**: Sound propagation
   - c = √(K/ρ) where K is bulk modulus, ρ is density
   - Air at 20°C: c ≈ 343 m/s

4. **Electromagnetic Waves**:
   - c = 1/√(εμ), in vacuum c ≈ 3×10⁸ m/s

**Advection Equation Applications:**

1. **Pollutant Transport**: Tracer in flowing river
2. **Traffic Flow**: Vehicle density waves
3. **Weather Patterns**: Temperature fronts moving with wind
4. **Numerical Test Problem**: Standard test for scheme accuracy

---

## The Wave Equation

### Energy Conservation

The wave equation conserves total energy:

```
E(t) = ∫_Ω [½(∂u/∂t)² + ½c²|∇u|²] dx
```

Where:
- `½(∂u/∂t)²` is the kinetic energy density
- `½c²|∇u|²` is the potential (strain) energy density

For homogeneous boundary conditions and f = 0, dE/dt = 0.

### Standing Waves

Solutions of the form u(x,t) = X(x)T(t) lead to standing wave modes:

**1D with Dirichlet BCs** (u(0,t) = u(L,t) = 0):
```
uₙ(x, t) = sin(nπx/L)·[Aₙcos(ωₙt) + Bₙsin(ωₙt)]
```

where ωₙ = nπc/L are the natural frequencies.

**2D rectangular membrane** (Dirichlet on all sides):
```
u_{mn}(x, y, t) = sin(mπx/Lₓ)·sin(nπy/L_y)·cos(ω_{mn}t)
```

where ω_{mn} = πc√((m/Lₓ)² + (n/L_y)²).

### d'Alembert's Solution (1D)

The general solution to the 1D wave equation (infinite domain):
```
u(x, t) = ½[u₀(x - ct) + u₀(x + ct)] + 1/(2c)∫_{x-ct}^{x+ct} v₀(s) ds
```

This shows:
- Initial displacement splits into left and right traveling waves
- Initial velocity creates a wave packet

---

## The Advection Equation

### Characteristics

The advection equation is a first-order hyperbolic PDE. Its characteristics are lines:
```
dx/dt = c
```

Along characteristics, u is constant: Du/Dt = 0.

This means information travels at speed c from left to right (if c > 0).

### Numerical Challenges

The advection equation is the canonical test problem for numerical schemes because:

1. **Centered differences are UNSTABLE**: FTCS (Forward Time, Central Space) is unconditionally unstable
2. **Upwinding is essential**: Information must come from upwind direction
3. **Numerical diffusion**: First-order schemes smear sharp features
4. **Numerical dispersion**: Second-order schemes create spurious oscillations

---

## Time Discretization Schemes

### For Wave Equation

#### Leapfrog Scheme (Default)

Second-order in time and space, non-dissipative:

```
(u^{n+1} - 2u^n + u^{n-1})/Δt² = c²·(u_{i-1} - 2u_i + u_{i+1})/h²
```

Rearranged:
```
u_i^{n+1} = 2u_i^n - u_i^{n-1} + r²·(u_{i-1}^n - 2u_i^n + u_{i+1}^n)
```

where r = cΔt/h is the Courant number.

**Properties:**
- ✅ Second-order accurate O(Δt², h²)
- ✅ Non-dissipative (preserves energy)
- ✅ Time-reversible
- ❌ Requires two time levels (startup needed)
- ❌ Conditionally stable (CFL: r ≤ 1)

**Startup:** Uses Taylor expansion for u^{-1}:
```
u^{-1} ≈ u^0 - Δt·v^0 + (Δt²/2)·c²·∇²u^0
```

### For Advection Equation

#### 1. Upwind Scheme (First-Order)

Uses information from the upwind direction:

**For c > 0 (flow to right):**
```
u_i^{n+1} = u_i^n - r·(u_i^n - u_{i-1}^n)
```

**For c < 0 (flow to left):**
```
u_i^{n+1} = u_i^n - r·(u_{i+1}^n - u_i^n)
```

**Properties:**
- ✅ Unconditionally stable for |r| ≤ 1
- ✅ Simple and robust
- ❌ First-order accurate O(Δt, h)
- ❌ High numerical diffusion

#### 2. Lax-Friedrichs Scheme

Centered differences with averaging:
```
u_i^{n+1} = ½(u_{i-1}^n + u_{i+1}^n) - (r/2)·(u_{i+1}^n - u_{i-1}^n)
```

**Properties:**
- ✅ Simple, symmetric
- ✅ Stable for |r| ≤ 1
- ❌ Very diffusive

#### 3. Lax-Wendroff Scheme (Second-Order)

Second-order accurate in time and space:
```
u_i^{n+1} = u_i^n - (r/2)·(u_{i+1}^n - u_{i-1}^n) + (r²/2)·(u_{i-1}^n - 2u_i^n + u_{i+1}^n)
```

**Properties:**
- ✅ Second-order accurate O(Δt², h²)
- ✅ Less diffusive than upwind
- ❌ Dispersive (creates oscillations near discontinuities)
- ❌ Stable for |r| ≤ 1

#### 4. Beam-Warming Scheme (Second-Order Upwind)

Second-order using upwind-biased stencil:

**For c > 0:**
```
u_i^{n+1} = u_i^n - (r/2)·(3u_i^n - 4u_{i-1}^n + u_{i-2}^n) 
          + (r²/2)·(u_i^n - 2u_{i-1}^n + u_{i-2}^n)
```

**Properties:**
- ✅ Second-order accurate
- ✅ Less dispersive than Lax-Wendroff for smooth solutions
- ❌ More complex boundary treatment

#### 5. FTCS Scheme (FOR REFERENCE ONLY - UNSTABLE!)

Forward Time, Central Space:
```
u_i^{n+1} = u_i^n - (r/2)·(u_{i+1}^n - u_{i-1}^n)
```

**⚠️ WARNING:** This scheme is **unconditionally unstable** for advection! 
It is included only for educational purposes to demonstrate instability.

---

## CFL Stability Condition

The Courant-Friedrichs-Lewy (CFL) condition is essential for explicit schemes.

### Wave Equation

**1D:** Δt ≤ h/c   (Courant number r = cΔt/h ≤ 1)

**2D:** Δt ≤ h/(c√2)   (diagonal wave propagation)

**3D:** Δt ≤ h/(c√3)

### Advection Equation

**1D:** Δt ≤ h/|c|   (Courant number r = |c|Δt/h ≤ 1)

**2D:** Δt ≤ 1/(|a_x|/hₓ + |a_y|/h_y)

### Physical Interpretation

The CFL condition ensures that **numerical information travels at least as fast as physical information**. If the time step is too large, information "outruns" the numerical stencil, causing instability.

```cpp
// Check CFL before solving
WaveSolver1D<double> solver(grid, bc, waveSpeed);
double dt = 0.5 * solver.computeCFLLimit();  // 50% of limit for safety
if (!solver.isStable(dt)) {
    throw std::runtime_error("Time step violates CFL condition!");
}
```

---

## 1D Wave Solver

### Class: `WaveSolver1D<T>`

```cpp
#include "mml/pde/hyperbolic/WaveSolver.h"

namespace MML::PDE {

template<typename T>
class WaveSolver1D {
public:
    // Construction
    WaveSolver1D(const Grid1D<T>& grid, 
                 const BoundaryConditions1D<T>& bc, 
                 T waveSpeed);
    
    // Initial conditions
    void setInitialDisplacement(std::function<T(T)> u0);
    void setInitialVelocity(std::function<T(T)> v0);
    
    // Optional source term
    void setSource(std::function<T(T, T)> f);
    
    // Observer for snapshots
    void setObserver(std::function<void(T, const std::vector<T>&, 
                                        const std::vector<T>&)> obs);
    
    // Time stepping
    bool step(T dt, WaveScheme scheme = WaveScheme::Leapfrog);
    WaveSolveResult<T> solve(T finalTime, T dt, 
                              WaveScheme scheme = WaveScheme::Leapfrog);
    
    // Stability
    T computeCFLLimit() const;
    T computeCourantNumber(T dt) const;
    bool isStable(T dt) const;
    
    // Energy
    T computeEnergy(T dt) const;
    
    // Accessors
    T getTime() const;
    const std::vector<T>& getDisplacement() const;
    const std::vector<T>& getVelocity() const;
    const Grid1D<T>& getGrid() const;
};

}  // namespace MML::PDE
```

### Basic Usage

```cpp
#include "mml/pde/hyperbolic/WaveSolver.h"
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"

using namespace MML::PDE;

// Create grid: string of length L=1 with 100 nodes
int n = 100;
Grid1D<double> grid(0.0, 1.0, n);

// Fixed ends (Dirichlet)
auto bc = homogeneousDirichlet1D<double>();

// Wave speed (determines frequency)
double c = 1.0;

// Create solver
WaveSolver1D<double> solver(grid, bc, c);

// Initial pluck: triangular shape
solver.setInitialDisplacement([](double x) {
    if (x < 0.5) return 2.0 * x;
    else return 2.0 * (1.0 - x);
});

// String initially at rest
solver.setInitialVelocity([](double x) { return 0.0; });

// Time step (80% of CFL limit)
double dt = 0.8 * solver.computeCFLLimit();

// Solve for one period
double period = 2.0 / c;  // T = 2L/c for fundamental mode
auto result = solver.solve(period, dt);

// Check results
std::cout << "Steps: " << result.steps << "\n";
std::cout << "Energy error: " << result.energyError * 100 << "%\n";
```

### Standing Wave Example

```cpp
// Fundamental mode: sin(πx)·cos(πct)
solver.setInitialDisplacement([](double x) {
    return std::sin(M_PI * x);
});
solver.setInitialVelocity([](double x) { return 0.0; });

// After time T/4 = 1/(2c), string should be flat (all energy kinetic)
auto result = solver.solve(0.5 / c, dt);

// After time T/2 = 1/c, string should be inverted
result = solver.solve(1.0 / c, dt);
```

### Traveling Wave Example

```cpp
// Gaussian pulse moving to the right
double x0 = 0.3;  // Initial center
double sigma = 0.05;  // Width

solver.setInitialDisplacement([x0, sigma](double x) {
    return std::exp(-(x - x0) * (x - x0) / (sigma * sigma));
});

// Initial velocity for right-traveling wave: v₀ = c·∂u₀/∂x
solver.setInitialVelocity([x0, sigma, c](double x) {
    double dudx = -2.0 * (x - x0) / (sigma * sigma) * 
                  std::exp(-(x - x0) * (x - x0) / (sigma * sigma));
    return c * dudx;
});
```

---

## 2D Wave Solver

### Class: `WaveSolver2D<T>`

```cpp
template<typename T>
class WaveSolver2D {
public:
    WaveSolver2D(const Grid2D<T>& grid, 
                 const BoundaryConditions2D<T>& bc, 
                 T waveSpeed);
    
    void setInitialDisplacement(std::function<T(T, T)> u0);
    void setInitialVelocity(std::function<T(T, T)> v0);
    void setSource(std::function<T(T, T, T)> f);
    
    bool step(T dt, WaveScheme scheme = WaveScheme::Leapfrog);
    WaveSolveResult<T> solve(T finalTime, T dt, WaveScheme scheme = WaveScheme::Leapfrog);
    
    T computeCFLLimit() const;  // h / (c * sqrt(2))
    bool isStable(T dt) const;
    T computeEnergy(T dt) const;
    
    const std::vector<T>& getDisplacement() const;
    const Grid2D<T>& getGrid() const;
};
```

### Vibrating Membrane Example

```cpp
// Square membrane [0,1] × [0,1] with fixed edges
Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, 50, 50);
auto bc = homogeneousDirichlet2D<double>();
double c = 1.0;

WaveSolver2D<double> solver(grid, bc, c);

// Fundamental mode: sin(πx)·sin(πy)·cos(√2·πt)
solver.setInitialDisplacement([](double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
});
solver.setInitialVelocity([](double x, double y) { return 0.0; });

// Time step (80% of 2D CFL limit)
double dt = 0.8 * solver.computeCFLLimit();

// Solve for quarter period
double omega = M_PI * c * std::sqrt(2.0);
double quarterPeriod = M_PI / (2.0 * omega);
auto result = solver.solve(quarterPeriod, dt);
```

---

## 1D Advection Solver

### Class: `AdvectionSolver1D<T>`

```cpp
template<typename T>
class AdvectionSolver1D {
public:
    AdvectionSolver1D(const Grid1D<T>& grid, T velocity, bool periodic = true);
    
    void setInitialCondition(std::function<T(T)> u0);
    void setPeriodic(bool periodic);
    void setObserver(std::function<void(T, const std::vector<T>&)> obs);
    
    bool step(T dt, AdvectionScheme scheme = AdvectionScheme::Upwind);
    AdvectionSolveResult<T> solve(T finalTime, T dt, 
                                   AdvectionScheme scheme = AdvectionScheme::Upwind);
    
    T computeCFLLimit() const;
    T computeCourantNumber(T dt) const;
    bool isStable(T dt) const;
    T computeMass() const;
    
    T getTime() const;
    const std::vector<T>& getSolution() const;
    T getVelocity() const;
};
```

### Translation Test Example

```cpp
#include "mml/pde/hyperbolic/AdvectionSolver.h"

// Create grid with periodic BCs
int n = 100;
Grid1D<double> grid(0.0, 10.0, n);
double velocity = 1.0;  // Move right at speed 1

AdvectionSolver1D<double> solver(grid, velocity, true);  // periodic

// Gaussian pulse
double x0 = 2.0;
double sigma = 0.3;
solver.setInitialCondition([x0, sigma](double x) {
    return std::exp(-(x - x0) * (x - x0) / (sigma * sigma));
});

// Time step
double dt = 0.8 * solver.computeCFLLimit();

// Advect for time t = 3.0, pulse should move to x = 5.0
auto result = solver.solve(3.0, dt, AdvectionScheme::Upwind);

// Check mass conservation
std::cout << "Mass error: " << result.massError * 100 << "%\n";
```

### Scheme Comparison Example

```cpp
// Compare different schemes on same problem
std::vector<AdvectionScheme> schemes = {
    AdvectionScheme::Upwind,
    AdvectionScheme::LaxFriedrichs,
    AdvectionScheme::LaxWendroff,
    AdvectionScheme::BeamWarming
};

for (auto scheme : schemes) {
    AdvectionSolver1D<double> solver(grid, velocity, true);
    solver.setInitialCondition(gaussian);
    
    auto result = solver.solve(3.0, dt, scheme);
    
    // Measure numerical diffusion (peak height reduction)
    double peakHeight = *std::max_element(
        solver.getSolution().begin(), 
        solver.getSolution().end()
    );
    std::cout << "Scheme peak: " << peakHeight << "\n";
}
```

---

## 2D Advection Solver

### Class: `AdvectionSolver2D<T>`

```cpp
template<typename T>
class AdvectionSolver2D {
public:
    AdvectionSolver2D(const Grid2D<T>& grid, 
                      T vx, T vy,
                      bool periodicX = true, 
                      bool periodicY = true);
    
    void setInitialCondition(std::function<T(T, T)> u0);
    AdvectionSolveResult<T> solve(T finalTime, T dt, 
                                   AdvectionScheme scheme = AdvectionScheme::Upwind);
    
    T computeCFLLimit() const;  // 1 / (|vx|/hx + |vy|/hy)
    bool isStable(T dt) const;
    
    const std::vector<T>& getSolution() const;
};
```

### Diagonal Transport Example

```cpp
// Transport in diagonal direction
Grid2D<double> grid(0.0, 10.0, 0.0, 10.0, 50, 50);
double vx = 1.0, vy = 1.0;  // Diagonal velocity

AdvectionSolver2D<double> solver(grid, vx, vy, true, true);

// Gaussian blob
double x0 = 2.0, y0 = 2.0;
double sigma = 0.5;
solver.setInitialCondition([x0, y0, sigma](double x, double y) {
    double dx = x - x0, dy = y - y0;
    return std::exp(-(dx*dx + dy*dy) / (sigma * sigma));
});

// Advect - blob should move to (x0+vx*t, y0+vy*t)
double dt = 0.5 * solver.computeCFLLimit();
auto result = solver.solve(2.0, dt);
```

---

## Numerical Issues

### Numerical Diffusion

First-order schemes (Upwind, Lax-Friedrichs) introduce **artificial viscosity** that smears sharp features:

```
Original: ___/‾‾‾‾\___
After diffusion: ___/‾‾‾‾\___  →  ___/~~‾‾~~\___
```

**Mitigation:**
- Use second-order schemes (Lax-Wendroff, Beam-Warming)
- Increase grid resolution
- Use high-resolution schemes (flux limiters) for shocks

### Numerical Dispersion

Second-order schemes introduce **phase errors** causing oscillations near discontinuities:

```
Original: ___/‾‾‾‾\___
With dispersion: ~~_/‾‾‾‾\~~~ (wiggles)
```

**Mitigation:**
- Avoid sharp discontinuities in initial conditions
- Use flux limiters or TVD schemes
- Smooth initial data if possible

### Instability Detection

The solvers automatically detect instability:

```cpp
auto result = solver.solve(finalTime, dt);
if (!result.stable) {
    std::cerr << "Solver failed: " << result.message << "\n";
    // Check CFL, reduce dt, or choose different scheme
}
```

---

## Complete Examples

### Example 1: Vibrating Guitar String

```cpp
#include "mml/pde/hyperbolic/WaveSolver.h"
#include "mml/pde/grid/Grid1D.h"
#include "mml/pde/grid/BoundaryConditions.h"
#include <fstream>

using namespace MML::PDE;

int main() {
    // Guitar string: L = 0.65m, c ≈ 400 m/s (A4 = 440 Hz)
    double L = 0.65;
    double c = 440.0 * 2.0 * L;  // c = f * 2L for fundamental
    int n = 200;
    
    Grid1D<double> grid(0.0, L, n);
    auto bc = homogeneousDirichlet1D<double>();
    
    WaveSolver1D<double> solver(grid, bc, c);
    
    // Pluck at 1/5 of length (typical guitar pluck position)
    double pluckPos = L / 5.0;
    solver.setInitialDisplacement([L, pluckPos](double x) {
        if (x < pluckPos) 
            return x / pluckPos;
        else 
            return (L - x) / (L - pluckPos);
    });
    solver.setInitialVelocity([](double x) { return 0.0; });
    
    // Record waveform at pickup position (near bridge)
    double pickupPos = 0.9 * L;
    int pickupIdx = static_cast<int>(pickupPos / grid.dx());
    
    std::ofstream waveform("guitar_waveform.dat");
    solver.setObserver([&waveform, pickupIdx](double t, 
                        const std::vector<double>& u,
                        const std::vector<double>& v) {
        waveform << t << " " << u[pickupIdx] << "\n";
    });
    
    double dt = 0.8 * solver.computeCFLLimit();
    double duration = 0.01;  // 10ms of sound
    
    auto result = solver.solve(duration, dt);
    
    std::cout << "Simulated " << result.steps << " time steps\n";
    std::cout << "Energy conservation error: " << result.energyError * 100 << "%\n";
    
    return 0;
}
```

### Example 2: Pollutant Transport in River

```cpp
#include "mml/pde/hyperbolic/AdvectionSolver.h"

using namespace MML::PDE;

int main() {
    // River section: 1 km long, flow velocity 0.5 m/s
    double L = 1000.0;  // meters
    double v = 0.5;     // m/s
    int n = 500;
    
    Grid1D<double> grid(0.0, L, n);
    
    // Non-periodic (outflow at right boundary)
    AdvectionSolver1D<double> solver(grid, v, false);
    
    // Initial pollutant concentration: localized spill at x = 100m
    double spillCenter = 100.0;
    double spillWidth = 20.0;
    solver.setInitialCondition([spillCenter, spillWidth](double x) {
        double d = x - spillCenter;
        return std::exp(-d * d / (spillWidth * spillWidth));
    });
    
    // Track plume location over time
    std::vector<double> times = {0, 200, 400, 600, 800, 1000};  // seconds
    double dt = 0.8 * solver.computeCFLLimit();
    
    for (double t : times) {
        auto result = solver.solve(t, dt, AdvectionScheme::LaxWendroff);
        
        // Find peak location
        const auto& u = solver.getSolution();
        auto maxIt = std::max_element(u.begin(), u.end());
        int maxIdx = std::distance(u.begin(), maxIt);
        double peakX = grid.x(maxIdx);
        
        std::cout << "t=" << t << "s: peak at x=" << peakX 
                  << "m (expected: " << spillCenter + v * t << "m)\n";
    }
    
    return 0;
}
```

### Example 3: 2D Drum Membrane

```cpp
#include "mml/pde/hyperbolic/WaveSolver.h"
#include <fstream>

using namespace MML::PDE;

int main() {
    // Square drum: 30cm × 30cm
    double L = 0.3;
    int n = 60;
    double c = 100.0;  // Wave speed in membrane
    
    Grid2D<double> grid(0.0, L, 0.0, L, n, n);
    auto bc = homogeneousDirichlet2D<double>();
    
    WaveSolver2D<double> solver(grid, bc, c);
    
    // Strike near center
    double x0 = 0.4 * L, y0 = 0.5 * L;
    double sigma = 0.02;  // 2cm radius impact
    
    solver.setInitialDisplacement([x0, y0, sigma](double x, double y) {
        double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        return std::exp(-r2 / (sigma * sigma));
    });
    solver.setInitialVelocity([](double x, double y) { return 0.0; });
    
    double dt = 0.8 * solver.computeCFLLimit();
    double duration = 0.01;  // 10ms
    
    auto result = solver.solve(duration, dt);
    
    // Output final displacement for visualization
    std::ofstream out("drum_displacement.dat");
    const auto& u = solver.getDisplacement();
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            out << grid.x(i) << " " << grid.y(j) << " " 
                << u[j * n + i] << "\n";
        }
        out << "\n";
    }
    
    std::cout << "Energy error: " << result.energyError * 100 << "%\n";
    return 0;
}
```

---

## Advanced Topics

### Observer Pattern

Both solvers support observers for recording solution history:

```cpp
// Record snapshots every 0.1 time units
std::vector<std::vector<double>> history;
std::vector<double> times;

solver.setObserver([&](double t, const std::vector<double>& u, 
                        const std::vector<double>& v) {
    if (times.empty() || t - times.back() >= 0.1) {
        history.push_back(u);
        times.push_back(t);
    }
});
```

### Custom Boundary Conditions

```cpp
// Time-dependent forcing at left boundary (driven oscillation)
double omega = 2.0 * M_PI;  // Driving frequency
auto bc = BoundaryConditions1D<double>(
    BoundaryCondition<double>(BCType::Dirichlet, [omega](double x, double t) {
        return std::sin(omega * t);
    }),
    BoundaryCondition<double>(BCType::Dirichlet, 0.0)  // Fixed right end
);
```

### Convergence Analysis

Verify second-order accuracy with grid refinement:

```cpp
std::vector<int> resolutions = {50, 100, 200, 400};
std::vector<double> errors;

for (int n : resolutions) {
    Grid1D<double> grid(0.0, 1.0, n);
    WaveSolver1D<double> solver(grid, bc, c);
    // ... set up problem with known solution ...
    
    double dt = 0.8 * solver.computeCFLLimit();
    solver.solve(finalTime, dt);
    
    // Compute error against exact solution
    double error = computeL2Error(solver.getDisplacement(), exactSolution);
    errors.push_back(error);
}

// Check convergence rate
for (size_t i = 1; i < errors.size(); ++i) {
    double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
    std::cout << "Refinement " << i << ": rate = " << rate << "\n";
    // Should be approximately 2.0 for second-order scheme
}
```

---

## API Reference

### Enumerations

```cpp
enum class WaveScheme {
    Leapfrog,      // Second-order, non-dissipative (default)
    LaxWendroff    // Second-order, some dissipation
};

enum class AdvectionScheme {
    Upwind,        // First-order, stable, diffusive (default)
    LaxFriedrichs, // First-order, simple, very diffusive
    LaxWendroff,   // Second-order, dispersive
    BeamWarming,   // Second-order upwind
    FTCS           // Unstable! For demonstration only
};
```

### Result Structures

```cpp
template<typename T>
struct WaveSolveResult {
    bool stable;           // Did solution remain stable?
    int steps;             // Number of time steps taken
    T finalTime;           // Final simulation time
    T maxDisplacement;     // Maximum value at final time
    T minDisplacement;     // Minimum value at final time
    T totalEnergy;         // Total energy at final time
    T initialEnergy;       // Initial total energy
    T energyError;         // Relative energy error |E-E₀|/E₀
    std::string message;   // Status message
};

template<typename T>
struct AdvectionSolveResult {
    bool stable;
    int steps;
    T finalTime;
    T maxValue;
    T minValue;
    T totalMass;          // Integral of u (conservation check)
    T initialMass;
    T massError;          // Relative mass error
    T maxAbsValue;        // For oscillation detection
    std::string message;
};
```

---

## See Also

- [Grid Infrastructure](Grid_Infrastructure.md) - Grid1D, Grid2D classes
- [Boundary Conditions](Boundary_Conditions.md) - BC types and setup
- [Elliptic PDEs](Elliptic_PDEs.md) - Steady-state problems
- [Parabolic PDEs](Parabolic_PDEs.md) - Diffusion equations
- [Examples Gallery](Examples_Gallery.md) - More complete examples

---

*Last updated: January 2026*
