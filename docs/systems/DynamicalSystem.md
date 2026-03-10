# DynamicalSystem - Comprehensive Dynamical Systems Analysis

`DynamicalSystem.h` provides MML's complete framework for analyzing continuous and discrete dynamical systems, featuring chaos detection, stability analysis, bifurcation diagrams, and Lyapunov exponent computation.

## Overview

```cpp
#include "systems/DynamicalSystem.h"
using namespace MML;

// Create a Lorenz system
LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // σ, ρ, β parameters

// Compute Lyapunov exponents (chaos detector)
auto lyap = LyapunovExponents::Compute(lorenz, 
    Vector<Real>({1.0, 1.0, 1.0}),  // Initial condition
    1000.0,   // Total time
    0.01,     // Time step
    100);     // Orthonormalization interval

std::cout << "Max Lyapunov exponent: " << lyap.maxExponent << "\n";
std::cout << "Chaotic? " << (lyap.isChaotic ? "Yes" : "No") << "\n";
```

## Design Philosophy

1. **Unified Interface** - `IDynamicalSystem` extends `IODESystemParametrized` with analysis capabilities
2. **Built-in Systems** - Classic examples ready to use (Lorenz, Rössler, Hénon, Logistic Map)
3. **Comprehensive Analysis** - Fixed points, stability, bifurcations, Lyapunov exponents
4. **Discrete & Continuous** - Both ODE systems and iterative maps
5. **Research-Ready** - Publication-quality tools for chaos and nonlinear dynamics

## Features

### IDynamicalSystem Interface

Extended interface for systems that support analysis:

```cpp
class MySystem : public DynamicalSystemBase<3, 2>  // 3D state, 2 parameters
{
public:
    void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
    {
        // Define your system dynamics
        dydt[0] = /* ... */;
        dydt[1] = /* ... */;
        dydt[2] = /* ... */;
    }
    
    // Optional: Analytical Jacobian (faster than numerical)
    bool hasAnalyticalJacobian() const override { return true; }
    
    void jacobian(Real t, const Vector<Real>& y, Matrix<Real>& J) const override
    {
        // Compute ∂f/∂x
        J(0,0) = /* ... */;
        // ... etc
    }
    
    // System properties
    bool isAutonomous() const override { return true; }
    bool isDissipative() const override { return true; }
    bool isHamiltonian() const override { return false; }
};
```

**Key Methods:**

- `getStateName(i)` - Name of state variable i
- `getParamName(i)` - Name of parameter i
- `getParamRange(i)` - Valid range for parameter i
- `getDefaultInitialCondition()` - Suggested starting state
- `hasAnalyticalJacobian()` - Does system provide analytical Jacobian?
- `jacobian(t, x, J)` - Compute Jacobian matrix ∂f/∂x
- `isAutonomous()` - No explicit time dependence?
- `isHamiltonian()` - Energy conserving?
- `isDissipative()` - Contracting phase space?

## Fixed Point Analysis

### Finding Fixed Points

Fixed points are equilibrium states where dx/dt = 0:

```cpp
VanDerPolSystem vdp(1.0);  // μ = 1.0

// Find fixed point near origin
Vector<Real> guess({0.1, 0.1});
auto fp = FixedPointFinder::Find(vdp, guess);

if (fp.convergenceResidual < 1e-10) {
    std::cout << "Found fixed point at: " << fp.location << "\n";
    std::cout << "Type: " << ToString(fp.type) << "\n";
    std::cout << "Stable? " << (fp.isStable ? "Yes" : "No") << "\n";
    std::cout << "Eigenvalues: ";
    for (const auto& ev : fp.eigenvalues)
        std::cout << ev.real() << "+" << ev.imag() << "i ";
    std::cout << "\n";
}
```

**FixedPoint Result Structure:**

```cpp
struct FixedPoint<Real> {
    Vector<Real> location;                      // Position in state space
    std::vector<std::complex<Real>> eigenvalues;// Jacobian eigenvalues
    Matrix<Real> jacobian;                      // Jacobian at fixed point
    FixedPointType type;                        // Stability classification
    bool isStable;                              // Overall stability
    Real convergenceResidual;                   // ||f(x*)|| at convergence
    int iterations;                             // Newton iterations used
};
```

### Fixed Point Types

```cpp
enum class FixedPointType {
    StableNode,        // All eigenvalues negative real (sink)
    UnstableNode,      // All eigenvalues positive real (source)
    Saddle,            // Mixed signs (saddle)
    StableFocus,       // Complex with negative real (spiral in)
    UnstableFocus,     // Complex with positive real (spiral out)
    Center,            // Purely imaginary (periodic orbits)
    StableStarNode,    // Degenerate stable node
    UnstableStarNode,  // Degenerate unstable node
    Unknown            // Could not classify
};
```

### Multiple Fixed Points

Find all fixed points from a grid of initial guesses:

```cpp
std::vector<Vector<Real>> guesses = {
    Vector<Real>({0.0, 0.0}),
    Vector<Real>({1.0, 1.0}),
    Vector<Real>({-1.0, -1.0}),
    Vector<Real>({1.0, -1.0})
};

auto fixedPoints = FixedPointFinder::FindMultiple(vdp, guesses);

std::cout << "Found " << fixedPoints.size() << " unique fixed points:\n";
for (const auto& fp : fixedPoints) {
    std::cout << "  Location: " << fp.location 
              << " Type: " << ToString(fp.type) << "\n";
}
```

## Lyapunov Exponents

Lyapunov exponents quantify chaos by measuring exponential divergence rates:

```cpp
LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);

auto result = LyapunovExponents::Compute(
    lorenz,
    Vector<Real>({1.0, 1.0, 1.0}),  // Initial condition
    1000.0,    // Total integration time
    0.01,      // Time step
    100        // Steps between QR decompositions
);

std::cout << "Lyapunov exponents: " << result.exponents << "\n";
std::cout << "Max exponent: " << result.maxExponent << "\n";
std::cout << "Sum: " << result.sum << "\n";
std::cout << "Kaplan-Yorke dimension: " << result.kaplanYorkeDimension << "\n";
std::cout << "Chaotic? " << (result.isChaotic ? "Yes" : "No") << "\n";
```

**LyapunovResult Structure:**

```cpp
struct LyapunovResult<Real> {
    Vector<Real> exponents;    // Lyapunov exponents (λ₁ ≥ λ₂ ≥ ... ≥ λₙ)
    Real maxExponent;          // Largest exponent (λ₁)
    Real sum;                  // Sum of all exponents
    Real kaplanYorkeDimension; // Kaplan-Yorke dimension
    bool isChaotic;            // True if λ₁ > 0
    int numOrthonormalizations;// Number of QR steps performed
    Real totalTime;            // Total integration time
};
```

**Interpretation:**

- **λ₁ > 0**: Chaotic dynamics (exponential divergence)
- **λ₁ ≈ 0**: Quasi-periodic or marginally stable
- **λ₁ < 0**: Stable fixed point or limit cycle
- **Sum of λᵢ < 0**: Dissipative system (volume contraction)
- **Kaplan-Yorke dimension**: Fractal dimension of attractor

### Algorithm Details

Uses **Continuous QR Decomposition** method:
1. Evolve tangent space along with trajectory
2. Periodically orthonormalize using QR decomposition
3. Accumulate logarithms of stretching factors
4. Average over total time

```cpp
// Control computation parameters
auto lyap = LyapunovExponents::Compute(
    system,
    x0,
    10000.0,  // Longer time → better convergence
    0.01,     // Smaller dt → more accurate
    50        // More frequent QR → better numerical stability
);
```

## Bifurcation Analysis

Sweep parameters and detect transitions in system behavior:

```cpp
LogisticMap logistic(3.0);

auto diagram = BifurcationAnalyzer::Sweep(
    logistic,
    0,              // Parameter index (r)
    2.8,            // Min parameter value
    4.0,            // Max parameter value
    500,            // Number of parameter steps
    Vector<Real>({0.5}),  // Initial condition
    0,              // Component to record (x)
    100.0,          // Transient time
    50.0,           // Recording time
    0.01            // Time step
);

// Plot bifurcation diagram
for (size_t i = 0; i < diagram.parameterValues.size(); ++i) {
    Real r = diagram.parameterValues[i];
    for (Real x : diagram.attractorValues[i]) {
        std::cout << r << " " << x << "\n";
    }
}
```

**Typical Bifurcations in Logistic Map:**
- r < 3.0: Fixed point
- r ≈ 3.0: Period-doubling (period-2)
- r ≈ 3.449: Period-4
- r ≈ 3.544: Period-8
- r ≈ 3.5699: Onset of chaos (Feigenbaum point)
- r = 4.0: Full chaos

### 2D System Bifurcations

```cpp
VanDerPolSystem vdp(0.1);

auto diagram = BifurcationAnalyzer::Sweep(
    vdp,
    0,              // Parameter index (μ)
    0.1,            // Min μ
    5.0,            // Max μ
    100,            // Steps
    Vector<Real>({2.0, 0.0}),  // Initial condition
    0,              // Record x component
    500.0,          // Longer transient
    200.0,          // Recording time
    0.01
);
```

## Discrete Maps

MML includes classic discrete maps with analytical properties:

### Logistic Map

```cpp
LogisticMap logistic(3.9);  // r = 3.9 (chaotic)

Vector<Real> x({0.5});
for (int i = 0; i < 100; ++i) {
    x = logistic.iterate(x);
    std::cout << x[0] << "\n";
}

// Compute Lyapunov exponent
auto lyap = DiscreteMapLyapunov<1>::compute(logistic, Vector<Real>({0.5}), 10000, 1000);
std::cout << "Lyapunov exponent: " << lyap.maxExponent << "\n";

// For r=4, there's an analytical formula
logistic.setR(4.0);
std::cout << "Analytical λ: " << logistic.analyticalLyapunov() << "\n";  // ≈ 0.693
```

### Hénon Map

Classic 2D chaotic map:

```cpp
HenonMap henon(1.4, 0.3);  // Classic parameters

Vector<Real> x({0.0, 0.0});
auto traj = henon.trajectory(x, 10000);

// Analyze attractor
for (const auto& point : traj) {
    std::cout << point[0] << " " << point[1] << "\n";
}

// Compute Lyapunov spectrum
auto lyap = DiscreteMapLyapunov<2>::compute(henon, x, 50000, 5000);
std::cout << "Lyapunov exponents: ";
for (Real exp : lyap.exponents)
    std::cout << exp << " ";
std::cout << "\n";
```

### Standard Map (Chirikov-Taylor)

Area-preserving kicked rotor:

```cpp
StandardMap standard(1.0);  // K = 1.0

// Generate orbit
Vector<Real> x({0.5, 1.0});  // (p, θ)
for (int i = 0; i < 1000; ++i) {
    x = standard.iterate(x);
    std::cout << x[0] << " " << x[1] << "\n";
}

// Check area preservation
Matrix<Real> J(2, 2);
standard.jacobian(x, J);
Real det = J(0,0)*J(1,1) - J(0,1)*J(1,0);
std::cout << "Determinant of Jacobian: " << det << "\n";  // Should be ≈ 1
```

### Tent Map

Exactly solvable chaotic map:

```cpp
TentMap tent(2.0);  // μ = 2.0

Vector<Real> x({0.3});
auto traj = tent.trajectory(x, 100);

// Exact Lyapunov exponent for μ ∈ (1, 2]
std::cout << "Analytical λ: " << tent.analyticalLyapunov() << "\n";  // log(2) ≈ 0.693

// Compare with numerical
auto lyap = DiscreteMapLyapunov<1>::compute(tent, x, 10000, 1000);
std::cout << "Numerical λ: " << lyap.maxExponent << "\n";
```

## Built-in Continuous Systems

### Lorenz System

The quintessential chaotic system:

```cpp
LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // σ, ρ, β

// System properties
std::cout << "Dimension: " << lorenz.getDim() << "\n";  // 3
std::cout << "Parameters: " << lorenz.getNumParam() << "\n";  // 3
std::cout << "Dissipative? " << lorenz.isDissipative() << "\n";  // true
std::cout << "Divergence: " << lorenz.getDivergence() << "\n";  // -(σ + 1 + β)

// Classic chaotic parameters
auto ic = lorenz.getDefaultInitialCondition();  // (1, 1, 1)

// Integrate trajectory
ODESystemSolver<Real> solver(lorenz);
auto [tVals, xVals] = solver.Solve_RungeKutta4(ic, 0.0, 50.0, 0.01);
```

### Rössler System

Simpler chaotic attractor:

```cpp
RosslerSystem rossler(0.2, 0.2, 5.7);  // a, b, c

// Classic parameters produce chaos with λ₁ ≈ 0.07
auto lyap = LyapunovExponents::Compute(rossler, 
    rossler.getDefaultInitialCondition(),
    5000.0, 0.01, 100);
```

### Van der Pol Oscillator

Self-sustained oscillations:

```cpp
VanDerPolSystem vdp(1.0);  // μ = 1.0

// Find the stable limit cycle
auto ic = vdp.getDefaultInitialCondition();  // (2, 0)

// For μ << 1: nearly sinusoidal
// For μ >> 1: relaxation oscillations
vdp.setParam(0, 5.0);  // Increase μ for relaxation behavior
```

### Duffing Oscillator

Driven nonlinear oscillator (non-autonomous):

```cpp
DuffingSystem duffing(1.0, 1.0, 0.3, 0.5, 1.2);  // α, β, δ, γ, ω

std::cout << "Autonomous? " << duffing.isAutonomous() << "\n";  // false

// State: [x, v, θ] where θ = ωt (phase of driving force)
// d²x/dt² + δ dx/dt + αx + βx³ = γ cos(θ)
```

### Chua's Circuit

Electronic chaos generator with double scroll:

```cpp
ChuaCircuit chua(15.6, 28.0, -1.143, -0.714);  // α, β, a, b

// Piecewise-linear system
// Famous for its double-scroll attractor
auto ic = Vector<Real>({0.1, 0.1, 0.1});
```

## Parameter Manipulation

All systems inherit parameter management:

```cpp
LorenzSystem lorenz;

// Set individual parameters
lorenz.setParam(0, 10.0);   // σ
lorenz.setParam(1, 28.0);   // ρ
lorenz.setParam(2, 8.0/3.0);// β

// Get individual parameters
Real sigma = lorenz.getParam(0);

// Set all parameters at once
lorenz.setParams(Vector<Real>({10.0, 28.0, 8.0/3.0}));

// Get all parameters
Vector<Real> params = lorenz.getParams();

// Query parameter names and ranges
std::cout << "Param 0: " << lorenz.getParamName(0) << "\n";  // "sigma"
auto [min, max] = lorenz.getParamRange(1);
std::cout << "ρ range: [" << min << ", " << max << "]\n";  // [0, 200]
```

## Advanced Features

### Poincaré Sections

Define sections through phase space:

```cpp
// Poincaré section at z = 27 with positive crossing
PoincareSection<Real> section(2, 27.0, +1);

// Use with BifurcationAnalyzer for return maps
// (Advanced feature - implementation details in header)
```

### Invariant Quantities

For Hamiltonian systems:

```cpp
class MyHamiltonianSystem : public DynamicalSystemBase<4, 1>
{
public:
    bool isHamiltonian() const override { return true; }
    
    int getNumInvariants() const override { return 1; }
    
    std::string getInvariantName(int i) const override { 
        return "Energy"; 
    }
    
    Real computeInvariant(int i, const Vector<Real>& x) const override {
        // H(q,p) = kinetic + potential energy
        return 0.5 * (x[2]*x[2] + x[3]*x[3]) + potential(x[0], x[1]);
    }
};
```

### Custom System Example

Complete example of defining a new system:

```cpp
// Predator-Prey (Lotka-Volterra)
class LotkaVolterraSystem : public DynamicalSystemBase<2, 4>
{
public:
    LotkaVolterraSystem(Real alpha = 1.0, Real beta = 0.1, 
                        Real gamma = 1.5, Real delta = 0.075)
    {
        _params[0] = alpha; _params[1] = beta;
        _params[2] = gamma; _params[3] = delta;
        
        _stateNames[0] = "prey";
        _stateNames[1] = "predator";
        
        _paramNames[0] = "alpha";  // Prey growth
        _paramNames[1] = "beta";   // Predation rate
        _paramNames[2] = "gamma";  // Predator death
        _paramNames[3] = "delta";  // Predator efficiency
        
        _paramRanges[0] = {0.0, 5.0};
        _paramRanges[1] = {0.0, 1.0};
        _paramRanges[2] = {0.0, 5.0};
        _paramRanges[3] = {0.0, 1.0};
    }
    
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
    {
        Real x = y[0], z = y[1];  // prey, predator
        Real alpha = _params[0], beta = _params[1];
        Real gamma = _params[2], delta = _params[3];
        
        dydt[0] = alpha * x - beta * x * z;         // dx/dt
        dydt[1] = delta * x * z - gamma * z;        // dz/dt
    }
    
    bool hasAnalyticalJacobian() const override { return true; }
    
    void jacobian(Real t, const Vector<Real>& y, Matrix<Real>& J) const override
    {
        Real x = y[0], z = y[1];
        Real alpha = _params[0], beta = _params[1];
        Real gamma = _params[2], delta = _params[3];
        
        J.Resize(2, 2);
        J(0,0) = alpha - beta * z;   J(0,1) = -beta * x;
        J(1,0) = delta * z;          J(1,1) = delta * x - gamma;
    }
    
    Vector<Real> getDefaultInitialCondition() const override
    {
        return Vector<Real>({10.0, 5.0});  // Initial populations
    }
};

// Usage
LotkaVolterraSystem lv;
auto fp = FixedPointFinder::Find(lv, Vector<Real>({5.0, 5.0}));
std::cout << "Equilibrium populations: " << fp.location << "\n";
std::cout << "Type: " << ToString(fp.type) << "\n";  // Likely Center
```

## Integration with ODE Solvers

All continuous systems work with MML's ODE solvers:

```cpp
LorenzSystem lorenz(10, 28, 8.0/3.0);
ODESystemSolver<Real> solver(lorenz);

auto ic = lorenz.getDefaultInitialCondition();

// Use any ODE method
auto sol1 = solver.Solve_RungeKutta4(ic, 0.0, 50.0, 0.01);
auto sol2 = solver.Solve_DormandPrince5_Adaptive(ic, 0.0, 50.0, 1e-6, 1e-8);
auto sol3 = solver.Solve_BackwardEuler(ic, 0.0, 50.0, 0.01);
```

## Performance Tips

1. **Analytical Jacobian** - Override `jacobian()` for 10x speedup in Lyapunov/fixed point computations
2. **Adaptive Stepsize** - Use Dormand-Prince for trajectory integration, fixed step for Lyapunov
3. **Transient Removal** - Always discard initial transient before analysis
4. **Long Time** - Lyapunov exponents need ~1000 orbital periods for convergence
5. **QR Frequency** - Balance numerical stability (small interval) vs speed (large interval)

## Common Patterns

### Complete Analysis Workflow

```cpp
// 1. Create system
LorenzSystem lorenz(10, 28, 8.0/3.0);

// 2. Find fixed points
auto fp1 = FixedPointFinder::Find(lorenz, Vector<Real>({0, 0, 0}));
auto fp2 = FixedPointFinder::Find(lorenz, Vector<Real>({8, 8, 27}));

// 3. Compute Lyapunov exponents
auto lyap = LyapunovExponents::Compute(lorenz,
    Vector<Real>({1, 1, 1}), 5000.0, 0.01, 100);

// 4. Generate bifurcation diagram
auto bif = BifurcationAnalyzer::Sweep(lorenz,
    1,      // Sweep ρ parameter
    0.0, 30.0, 300,
    Vector<Real>({1, 1, 1}),
    2,      // Record z component
    500.0, 200.0, 0.01);

// 5. Integrate trajectory for visualization
ODESystemSolver<Real> solver(lorenz);
auto [t, x] = solver.Solve_RungeKutta4(
    lorenz.getDefaultInitialCondition(), 0.0, 100.0, 0.01);
```

### Chaos Identification Checklist

A system is **chaotic** if:
1. ✓ Largest Lyapunov exponent > 0
2. ✓ Sensitive dependence on initial conditions
3. ✓ Trajectory is bounded (not escaping to infinity)
4. ✓ Aperiodic (never exactly repeats)

```cpp
// Test for chaos
auto lyap = LyapunovExponents::Compute(system, x0, 5000.0, 0.01, 100);
if (lyap.isChaotic && lyap.sum < 0) {
    std::cout << "System exhibits dissipative chaos\n";
    std::cout << "Attractor dimension ≈ " << lyap.kaplanYorkeDimension << "\n";
}
```

## Implementation Details

| Feature | Algorithm | Reference |
|---------|-----------|-----------|
| Fixed Points | Newton-Raphson with Jacobian | Numerical Recipes |
| Stability | Eigenvalues of Jacobian | Linear Algebra |
| Lyapunov Exponents | Continuous QR method | Benettin et al. 1980 |
| Bifurcation | Parameter continuation | Strogatz, Nonlinear Dynamics |
| Discrete Maps | Direct iteration | Ott, Chaos in Dynamical Systems |

## See Also

- [IODESystem.h](../interfaces/IODESystem.md) - Base ODE system interface
- [ODESystemSolver.h](../algorithms/ODESystemSolver.md) - ODE integration methods
- [EigenSystemSolvers.h](../algorithms/EigenSystemSolvers.md) - Eigenvalue computation
- [RootFinding.h](../algorithms/RootFinding.md) - Newton's method for fixed points

## References

1. Strogatz, S. H. (2018). *Nonlinear Dynamics and Chaos*
2. Ott, E. (2002). *Chaos in Dynamical Systems*
3. Sprott, J. C. (2003). *Chaos and Time-Series Analysis*
4. Benettin, G. et al. (1980). "Lyapunov Characteristic Exponents for smooth dynamical systems"
