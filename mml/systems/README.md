# Systems Package

**Dynamical Systems Analysis**

The Systems package provides comprehensive tools for analyzing dynamical systems, including fixed point analysis, Lyapunov exponent computation, bifurcation diagrams, phase space tools, and a unified linear algebra interface.

## Features

### Dynamical Systems Analysis
- **Fixed Point Detection** - Find equilibria using Newton-Raphson
- **Stability Classification** - Nodes, foci, saddles, centers
- **Lyapunov Exponents** - Full spectrum computation
- **Bifurcation Diagrams** - Parameter sweep analysis
- **Poincaré Sections** - Return maps and phase space slicing

### Phase Space Tools
- **Trajectory Integration** - Adaptive ODE solvers
- **Phase Portraits** - Vector field visualization data
- **Limit Cycles** - Periodic orbit detection
- **Chaos Detection** - Lyapunov-based chaos indicators

### Linear Systems Interface
- **Unified Facade** - One class for all linear algebra
- **Smart Solver Selection** - Automatic algorithm choice
- **Matrix Decompositions** - LU, QR, SVD, Cholesky
- **Condition Analysis** - Stability assessment and diagnostics

## Quick Start

### Fixed Point Analysis

```cpp
#include "systems/DynamicalSystem.h"

using namespace MML::Systems;

// Define a 2D system: dx/dt = f(x, y)
class VanDerPol : public IDynamicalSystem {
public:
    VanDerPol(double mu) : mu_(mu) {}
    
    int getDim() const override { return 2; }
    
    void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
        dxdt[0] = x[1];
        dxdt[1] = mu_ * (1 - x[0]*x[0]) * x[1] - x[0];
    }
    
private:
    double mu_;
};

// Create system
VanDerPol system(1.0);  // μ = 1.0

// Find fixed points
DynamicalSystemAnalyzer analyzer;
auto fixedPoints = analyzer.findFixedPoints(system, searchRegion);

for (const auto& fp : fixedPoints) {
    std::cout << "Fixed point at: " << fp.location << "\n";
    std::cout << "Type: " << ToString(fp.type) << "\n";
    std::cout << "Eigenvalues: ";
    for (const auto& ev : fp.eigenvalues)
        std::cout << ev << " ";
    std::cout << "\n";
}
```

### Lyapunov Exponents

```cpp
#include "systems/DynamicalSystem.h"

using namespace MML::Systems;

// Lorenz system (classic chaos)
class Lorenz : public IDynamicalSystem {
public:
    Lorenz(double sigma = 10.0, double rho = 28.0, double beta = 8.0/3.0)
        : sigma_(sigma), rho_(rho), beta_(beta) {}
    
    int getDim() const override { return 3; }
    
    void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
        dxdt[0] = sigma_ * (x[1] - x[0]);
        dxdt[1] = x[0] * (rho_ - x[2]) - x[1];
        dxdt[2] = x[0] * x[1] - beta_ * x[2];
    }
    
private:
    double sigma_, rho_, beta_;
};

Lorenz lorenz;
Vector<Real> x0 = {1.0, 1.0, 1.0};  // Initial condition
double totalTime = 1000.0;

DynamicalSystemAnalyzer analyzer;
auto result = analyzer.computeLyapunovSpectrum(lorenz, x0, totalTime);

std::cout << "Lyapunov exponents: " << result.exponents << "\n";
std::cout << "Max exponent (λ₁): " << result.maxExponent << "\n";
std::cout << "Sum (should be negative for dissipative): " << result.sum << "\n";
std::cout << "Kaplan-Yorke dimension: " << result.kaplanYorkeDimension << "\n";
std::cout << "Is chaotic: " << (result.isChaotic ? "Yes" : "No") << "\n";
```

### Bifurcation Diagram

```cpp
#include "systems/DynamicalSystem.h"

using namespace MML::Systems;

// Logistic map: x_{n+1} = r * x_n * (1 - x_n)
class LogisticMap : public IDiscreteSystem {
public:
    LogisticMap(double r) : r_(r) {}
    
    void step(const Vector<Real>& x, Vector<Real>& next) const override {
        next[0] = r_ * x[0] * (1 - x[0]);
    }
    
    void setParameter(int i, Real val) override { r_ = val; }
    
private:
    double r_;
};

DynamicalSystemAnalyzer analyzer;

// Sweep r from 2.5 to 4.0
auto diagram = analyzer.bifurcationDiagram(
    LogisticMap(3.0),   // System template
    0,                  // Parameter index
    2.5, 4.0,           // Parameter range
    1000,               // Number of parameter values
    {0.5},              // Initial condition
    500,                // Transient steps to discard
    100                 // Points to record per parameter
);

// diagram.parameterValues[] and diagram.attractorValues[][] 
// can be plotted to visualize period-doubling cascade to chaos
```

### Poincaré Section

```cpp
#include "systems/DynamicalSystem.h"

using namespace MML::Systems;

Lorenz lorenz;
Vector<Real> x0 = {1.0, 1.0, 20.0};

// Define section: z = 27 (crossing from below)
PoincareSection<Real> section(2, 27.0, +1);  // variable 2 (z), value 27, positive direction

DynamicalSystemAnalyzer analyzer;
auto crossings = analyzer.poincareMap(lorenz, x0, section, 10000.0, 1000);

// crossings contains (x, y) values at each z=27 crossing
// Plot to see the characteristic Lorenz attractor cross-section
```

### Linear System Interface

```cpp
#include "systems/LinearSystem.h"

using namespace MML::Systems;

// Create system from matrix
Matrix<Real> A = {{4, 1, 2}, {1, 5, 1}, {2, 1, 6}};
LinearSystem<Real> system(A);

// Solve Ax = b
Vector<Real> b = {1, 2, 3};
auto solution = system.solve(b);

std::cout << "Solution: " << solution.x << "\n";
std::cout << "Residual: " << solution.verification.absoluteResidual << "\n";

// Get matrix properties
auto props = system.analyze();
std::cout << "Condition number: " << props.conditionNumber << "\n";
std::cout << "Is symmetric: " << props.isSymmetric << "\n";
std::cout << "Is positive definite: " << props.isPositiveDefinite << "\n";
std::cout << "Stability: " << ToString(props.stability) << "\n";

// Decompositions (cached)
auto lu = system.getLU();
auto qr = system.getQR();
auto svd = system.getSVD();

// Eigenvalue analysis
auto eigen = system.getEigenvalues();
std::cout << "Eigenvalues: " << eigen << "\n";
```

## API Reference

### Fixed Point Types

| Type | Eigenvalue Pattern | Behavior |
|------|-------------------|----------|
| `StableNode` | All λ < 0 (real) | Trajectories converge |
| `UnstableNode` | All λ > 0 (real) | Trajectories diverge |
| `Saddle` | Mixed signs | Stable/unstable manifolds |
| `StableFocus` | Re(λ) < 0 (complex) | Spiral inward |
| `UnstableFocus` | Re(λ) > 0 (complex) | Spiral outward |
| `Center` | Re(λ) = 0 (imaginary) | Periodic orbits |

### Result Structures

```cpp
// Fixed point result
struct FixedPoint<Type> {
    Vector<Type> location;              // Position in state space
    std::vector<std::complex<Type>> eigenvalues;  // Jacobian eigenvalues
    Matrix<Type> jacobian;              // Jacobian at fixed point
    FixedPointType type;                // Stability classification
    bool isStable;                      // Overall stability
    Type convergenceResidual;           // ||f(x*)|| at convergence
    int iterations;                     // Newton iterations used
};

// Lyapunov result
struct LyapunovResult<Type> {
    Vector<Type> exponents;             // λ₁ ≥ λ₂ ≥ ... ≥ λₙ
    Type maxExponent;                   // λ₁
    Type sum;                           // Σλᵢ
    Type kaplanYorkeDimension;          // Fractal dimension estimate
    bool isChaotic;                     // True if λ₁ > 0
    int numOrthonormalizations;         // QR steps performed
    Type totalTime;                     // Integration time
};

// Bifurcation diagram
struct BifurcationDiagram<Type> {
    std::string parameterName;
    std::vector<Type> parameterValues;
    std::vector<std::vector<Type>> attractorValues;
    int numTransientSteps;
    int numRecordedPoints;
};
```

### IDynamicalSystem Interface

| Method | Description |
|--------|-------------|
| `getDim()` | Get state space dimension |
| `derivs(t, x, dxdt)` | Compute dx/dt = f(x, t) |
| `jacobian(t, x, J)` | Compute Jacobian matrix |
| `getStateName(i)` | Get name of state variable i |
| `getParamName(i)` | Get name of parameter i |
| `isAutonomous()` | True if no explicit time dependence |
| `isHamiltonian()` | True if energy-conserving |

### LinearSystem Interface

| Method | Description |
|--------|-------------|
| `solve(b)` | Solve Ax = b |
| `analyze()` | Get matrix properties |
| `getLU()` | LU decomposition A = PLU |
| `getQR()` | QR decomposition A = QR |
| `getSVD()` | SVD A = UΣVᵀ |
| `getCholesky()` | Cholesky A = LLᵀ (SPD only) |
| `getEigenvalues()` | Compute eigenvalues |
| `getConditionNumber()` | Condition number κ(A) |
| `inverse()` | Compute A⁻¹ |

## File Structure

```
systems/
├── README.md              # This file
├── CMakeLists.txt         # Build configuration
│
├── include/
│   ├── DynamicalSystem.h  # Main dynamical systems analysis
│   │   - IDynamicalSystem interface
│   │   - Fixed point analysis
│   │   - Lyapunov exponents
│   │   - Bifurcation diagrams
│   │   - Poincaré sections
│   │
│   └── LinearSystem.h     # Unified linear algebra facade
│       - Smart solver selection
│       - Matrix decompositions
│       - Condition analysis
│       - Solution verification
│
└── tests/
    └── systems_tests.cpp
```

## Mathematical Background

### Lyapunov Exponents

The Lyapunov exponent measures the rate of separation of infinitesimally close trajectories:

$$\lambda = \lim_{t \to \infty} \frac{1}{t} \ln \frac{\|\delta x(t)\|}{\|\delta x(0)\|}$$

Computed via QR decomposition of the variational equation solution.

### Fixed Point Stability

At equilibrium $x^*$ where $f(x^*) = 0$, linearize:

$$\dot{y} = Df(x^*) \cdot y$$

Stability determined by eigenvalues of Jacobian $Df(x^*)$.

### Kaplan-Yorke Dimension

$$D_{KY} = j + \frac{\sum_{i=1}^{j} \lambda_i}{|\lambda_{j+1}|}$$

where $j$ is largest integer such that $\sum_{i=1}^{j} \lambda_i \geq 0$.

### Condition Number

$$\kappa(A) = \|A\| \cdot \|A^{-1}\| = \frac{\sigma_{\max}}{\sigma_{\min}}$$

Measures sensitivity of linear system solution to perturbations.

## Example Systems

### Classic Chaotic Systems

| System | Dimension | Behavior |
|--------|-----------|----------|
| Lorenz | 3D | Strange attractor, chaos |
| Rössler | 3D | Simpler chaotic attractor |
| Chua | 3D | Double-scroll attractor |
| Hénon | 2D (map) | Fractal structure |
| Logistic | 1D (map) | Period-doubling to chaos |

### Classic Oscillators

| System | Dimension | Behavior |
|--------|-----------|----------|
| Van der Pol | 2D | Limit cycle |
| Duffing | 2D | Nonlinear oscillator |
| Pendulum | 2D | Periodic/chaotic |
| Lotka-Volterra | 2D | Predator-prey cycles |

## References

- Strogatz, S.H. (2015). "Nonlinear Dynamics and Chaos" (2nd ed.)
- Ott, E. (2002). "Chaos in Dynamical Systems" (2nd ed.)
- Parker, T.S., and Chua, L.O. (1989). "Practical Numerical Algorithms for Chaotic Systems"
- Numerical Recipes in C++, 3rd Edition, Chapter 17: Integration of ODEs

## See Also

- [mml_packages/README.md](../README.md) - Package overview
- [Optimization Package](../optimization/README.md) - Parameter optimization
- [Symbolic Package](../symbolic/README.md) - Jacobian computation via AD
- [MML algorithms/ODESolverAdaptive.h](../../../mml/algorithms/) - ODE solvers
