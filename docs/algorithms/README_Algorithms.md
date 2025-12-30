# MML Algorithms Layer Documentation

## Overview

The Algorithms layer implements numerical methods for solving mathematical problems. These are production-ready implementations of classical and modern algorithms for root finding, solving differential equations, optimization, and statistical analysis.

**Files**: `mml/algorithms/*.h`

### Algorithms Overview and Quick Navigation
- **Monte Carlo integration**: See [MonteCarlo_integration.md](MonteCarlo_integration.md) for plain/stratified/hit-or-miss methods and error behavior.
- **Root finding**: See [Root_finding.md](Root_finding.md) for bracketing, bisection, Newton, and secant (1D and systems).
- **ODE solvers**: See [Differential_equations_solvers.md](Differential_equations_solvers.md) for fixed-step vs. adaptive steppers and Jacobian usage.
- **Eigen solvers**: See [Eigen_solvers.md](Eigen_solvers.md) for symmetric vs. nonsymmetric approaches and verification tips.
- **Optimization**: See [Function_optimization.md](Function_optimization.md) for line search, gradient, conjugate gradient, quasi-Newton.
- **Function analysis**: See [Function_analyzer.md](Function_analyzer.md) for critical points and properties.
- **Statistics**: See [Statistics.md](Statistics.md) for descriptive stats, distributions, and hypothesis testing.
- **Fourier**: See [Fourier_transformation.md](Fourier_transformation.md) for spectral analysis workflows.
- **Path integration**: See [Path_integration.md](Path_integration.md) for line/surface/volume integrals.

## Contents

1. [Root Finding](#root-finding)
   - [Bracketing Methods](#bracketing-methods)
   - [Newton-Raphson](#newton-raphson)
   - [Secant Method](#secant-method)
   - [Multidimensional Roots](#multidimensional-roots)
2. [ODE Solvers](#ode-solvers)
   - [Fixed-Step Methods](#fixed-step-methods)
   - [Adaptive Methods](#adaptive-methods)
   - [Specialized Solvers](#specialized-solvers)
3. [Optimization](#optimization)
   - [Line Search](#line-search)
   - [Gradient Methods](#gradient-methods)
   - [Conjugate Gradient](#conjugate-gradient)
   - [Quasi-Newton Methods](#quasi-newton-methods)
4. [Heuristic Optimization](#heuristic-optimization)
   - [Simulated Annealing](#simulated-annealing)
5. [Function Analysis](#function-analysis)
   - [Critical Points](#critical-points)
   - [Extrema Finding](#extrema-finding)
   - [Function Properties](#function-properties)
6. [Statistics](#statistics)
   - [Descriptive Statistics](#descriptive-statistics)
   - [Order Statistics](#order-statistics)
   - [Robust Statistics](#robust-statistics)
   - [Correlation and Covariance](#correlation-and-covariance)
   - [Distributions](#distributions)
   - [Hypothesis Testing](#hypothesis-testing)
   - [Confidence Intervals](#confidence-intervals)

---

## Root Finding

**File**: `mml/algorithms/RootFinding.h`

**Purpose**: Find values x where f(x) = 0.

### Bracketing Methods

#### Bracket Root

**Algorithm**: Geometrically expand interval until root is bracketed.

```cpp
#include "algorithms/RootFinding.h"

RealFunction f([](Real x) { return x*x - 2.0; });  // Find √2

Real x1 = 0.0, x2 = 1.0;
bool found = RootFinding::BracketRoot(f, x1, x2);
if (found) {
    // x1 and x2 now bracket a root: f(x1)*f(x2) < 0
}
```

**Properties**:
- **Complexity**: O(log n) expansions
- **Robustness**: Geomet ric scaling (factor 1.6)
- **Failure Mode**: Returns false if range becomes too large

#### Find Root Brackets

**Algorithm**: Subdivide interval and find all sign changes.

```cpp
RealFunction f([](Real x) { 
    return std::sin(x);  // Multiple roots
});

Vector<Real> x1, x2;
int numRoots = RootFinding::FindRootBrackets(f, 0.0, 10.0, 100, x1, x2);

// numRoots = number of bracketed intervals found
for (int i = 0; i < numRoots; ++i) {
    // Root exists in [x1[i], x2[i]]
}
```

**Properties**:
- **Complexity**: O(n) function evaluations
- **Use Case**: Find all roots in interval
- **Note**: May miss roots if sampling is too coarse

---

### Bisection Method

**Algorithm**: Repeatedly halve interval until convergence.

```cpp
RealFunction f([](Real x) { return x*x - 2.0; });

Real x1 = 1.0, x2 = 2.0;  // Must bracket root
Real tol = 1e-10;
Real root = RootFinding::FindRootBisection(f, x1, x2, tol);
// root ≈ 1.414213562 (√2)
```

**Properties**:
- **Convergence**: Linear, guaranteed
- **Rate**: Halves error each iteration
- **Complexity**: O(log(ε⁻¹)) iterations
- **Robustness**: Most reliable method
- **Speed**: Slowest of all methods

**Pros**:
✅ Always converges if root is bracketed  
✅ No derivative needed  
✅ Predictable iterations  

**Cons**:
❌ Slow convergence  
❌ Requires bracket  
❌ Only one root at a time  

---

### Newton-Raphson Method

**Algorithm**: Use tangent line approximation.

```
x_{n+1} = x_n - f(x_n)/f'(x_n)
```

```cpp
RealFunction f([](Real x) { return x*x - 2.0; });
// Derivative computed numerically internally

Real x1 = 1.0, x2 = 2.0;  // Initial bracket
Real tol = 1e-10;
Real root = RootFinding::FindRootNewton(f, x1, x2, tol);
// Converges in ~5 iterations
```

**Properties**:
- **Convergence**: Quadratic (near root)
- **Rate**: Error ~ error²
- **Complexity**: O(log log ε⁻¹) iterations
- **Speed**: Very fast near root
- **Requirement**: f'(x) ≠ 0

**Pros**:
✅ Quadratic convergence  
✅ Very fast near root  
✅ Few iterations needed  

**Cons**:
❌ Can diverge if initial guess poor  
❌ Fails if f'(x) = 0  
❌ Requires derivative  

**Example: Finding nth Root**

```cpp
// Find ⁵√32 (fifth root of 32)
RealFunction f([](Real x) { 
    return std::pow(x, 5) - 32.0; 
});

Real root = RootFinding::FindRootNewton(f, 1.0, 3.0, 1e-10);
// root = 2.0
```

---

### Secant Method

**Algorithm**: Approximate derivative using secant line.

```
x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
```

```cpp
RealFunction f([](Real x) { return std::cos(x) - x; });

Real x1 = 0.0, x2 = 1.0;
Real tol = 1e-10;
Real root = RootFinding::FindRootSecant(f, x1, x2, tol);
// root ≈ 0.739085 (solution to cos(x) = x)
```

**Properties**:
- **Convergence**: Superlinear (rate ≈ 1.618)
- **Speed**: Between bisection and Newton-Raphson
- **Advantage**: No derivative needed
- **Complexity**: O(log ε⁻¹) iterations

**Pros**:
✅ No derivative needed  
✅ Faster than bisection  
✅ More robust than Newton  

**Cons**:
❌ Slower than Newton-Raphson  
❌ Can still diverge  
❌ Needs two initial points  

---

### Multidimensional Root Finding

**File**: `mml/algorithms/RootFinding.h` (Newton-Raphson for systems)

**Problem**: Find **x** where **F**(**x**) = **0**, **F**: ℝⁿ → ℝⁿ

```cpp
// System: f₁(x,y) = x² + y² - 1 = 0
//         f₂(x,y) = x - y = 0
// Solution: intersection of circle and line

VectorFunction<2> F([](const VectorN<Real, 2>& v) {
    Real x = v[0], y = v[1];
    return VectorN<Real, 2>{
        x*x + y*y - 1.0,  // f₁
        x - y             // f₂
    };
});

VectorN<Real, 2> initial_guess{0.5, 0.5};
VectorN<Real, 2> root = RootFinding::NewtonRaphsonSystem(F, initial_guess, 1e-10);
// root ≈ [0.707, 0.707] = [1/√2, 1/√2]
```

**Algorithm**:
1. Compute Jacobian **J** = ∂**F**/∂**x**
2. Solve **J** Δ**x** = -**F**(**x**)
3. Update **x** ← **x** + Δ**x**
4. Repeat until ||**F**(**x**)|| < tol

**Properties**:
- **Convergence**: Quadratic (if Jacobian nonsingular)
- **Complexity**: O(n³) per iteration (Jacobian solve)
- **Robustness**: Damped Newton or trust region for global convergence

---

## ODE Solvers

**Files**: `mml/algorithms/ODESystemSolver.h`, `ODESystemStepCalculators.h`, `ODEAdaptiveIntegrator.h`

**Problem**: Solve dy/dt = f(t, y), given y(t₀) = y₀

### Architecture

MML uses a modular ODE solver design:

- **ODESystem**: Defines the system dy/dt = f(t, y)
- **StepCalculator**: Computes single step (Euler, RK4, etc.) for fixed-step integration
- **ODESystemFixedStepSolver**: Fixed-step integration using step calculators
- **ODEAdaptiveIntegrator**: Adaptive step-size integration with error control
- **Stepper**: Adaptive steppers (DormandPrince5, CashKarp, DormandPrince8)

### Fixed-Step Methods

#### Euler Method

**Algorithm**: Simplest first-order method.

```
y_{n+1} = y_n + h * f(t_n, y_n)
```

```cpp
#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"

// Define system: dy/dt = -y (exponential decay)
ODESystem system(1, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    dydt[0] = -y[0];
});

// Create Euler step calculator
EulerStep_Calculator euler;

// Create fixed-step solver
ODESystemFixedStepSolver solver(system, euler);

// Solve
Vector<Real> y0({1.0});  // Initial condition
ODESystemSolution sol = solver.integrate(y0, 0.0, 5.0, 100);

// sol contains solution at 101 points (including initial)
```

**Properties**:
- **Order**: 1st order (error ~ h)
- **Stability**: Conditionally stable
- **Speed**: Fastest per step
- **Accuracy**: Poorest

**Use When**: Quick approximation, educational purposes

#### RK4 (Runge-Kutta 4th Order)

**Algorithm**: Classic 4th-order method.

```
k₁ = h*f(t, y)
k₂ = h*f(t + h/2, y + k₁/2)
k₃ = h*f(t + h/2, y + k₂/2)
k₄ = h*f(t + h, y + k₃)
y_{n+1} = y_n + (k₁ + 2k₂ + 2k₃ + k₄)/6
```

```cpp
RungeKutta4_StepCalculator rk4;
ODESystemFixedStepSolver solver(system, rk4);

Vector<Real> y0({1.0});
ODESystemSolution sol = solver.integrate(y0, 0.0, 5.0, 50);
// Much more accurate than Euler with same number of steps
```

**Properties**:
- **Order**: 4th order (error ~ h⁴)
- **Stability**: Good
- **Speed**: 4 function evaluations per step
- **Accuracy**: Excellent

**Pros**:
✅ High accuracy  
✅ Stable  
✅ Industry standard  

**Cons**:
❌ Fixed step size (may waste computation)  
❌ 4× slower than Euler per step  

#### Midpoint Method

**Algorithm**: 2nd-order method.

```cpp
Midpoint_StepCalculator midpoint;
ODESystemFixedStepSolver solver(system, midpoint);
```

**Properties**:
- **Order**: 2nd order (error ~ h²)
- **Evaluations**: 2 per step
- **Accuracy**: Better than Euler, worse than RK4

---

### Adaptive Methods

#### Cash-Karp 5(4) (Recommended)

**Algorithm**: 5th-order method with embedded 4th-order error estimate.

```cpp
#include "algorithms/ODEAdaptiveIntegrator.h"

// Create adaptive integrator with Cash-Karp stepper
ODEAdaptiveIntegrator<CashKarp_Stepper> integrator(system);

// Solve with error tolerance
Real eps = 1e-6;        // Error tolerance
Real saveInterval = 0.1; // Save points at this interval
Vector<Real> y0({1.0});

ODESystemSolution sol = integrator.integrate(y0, 0.0, 5.0, saveInterval, eps);
// Automatically adjusts step size to maintain accuracy

// Access solution statistics
SolutionStatistics stats = integrator.getStatistics();
std::cout << "Accepted steps: " << stats.acceptedSteps << std::endl;
std::cout << "Rejected steps: " << stats.rejectedSteps << std::endl;
```

**Properties**:
- **Order**: 5th order with 4th order error estimate
- **Adaptivity**: Embedded error control with PI controller
- **Efficiency**: Fewer steps than fixed-step for same accuracy
- **Evaluations**: 6 per step

**Error Control**:
```
scale_i = eps * (|y_i| + |h * dxdt_i| + tiny)
error_i = |y_5th_i - y_4th_i| / scale_i
acceptable if max(error_i) ≤ 1.0
```

**Step Size Adjustment**:
```
h_new = h_old * safety_factor * (1.0/error)^0.2
```

#### Dormand-Prince 5(4) (FSAL)

**Algorithm**: Industry-standard method used by MATLAB ode45 and SciPy RK45.

```cpp
ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(system);
ODESystemSolution sol = integrator.integrate(y0, 0.0, 5.0, 0.1, 1e-6);
```

**Properties**:
- **FSAL**: First Same As Last - reuses final derivative (only 6 evals per step)
- **Dense Output**: 4th-order interpolation within steps
- **PI Control**: Smooth step-size adaptation

#### Dormand-Prince 8(7) (High-Order)

**Algorithm**: 8th-order method for high-precision requirements.

```cpp
ODEAdaptiveIntegrator<DormandPrince8_Stepper> integrator(system);
ODESystemSolution sol = integrator.integrate(y0, 0.0, 5.0, 0.1, 1e-10);
```

**Properties**:
- **Order**: 8th order with 7th order error estimate
- **Use Case**: High-precision long-time integrations
- **Evaluations**: 13 per step

---

### Specialized Solvers

#### Leapfrog (Verlet Integration)

**Purpose**: Energy-conserving integrator for Hamiltonian systems.

**Algorithm**: Symplectic integrator.

```
v_{n+1/2} = v_n + (h/2) * a_n
x_{n+1} = x_n + h * v_{n+1/2}
v_{n+1} = v_{n+1/2} + (h/2) * a_{n+1}
```

```cpp
// Use Leapfrog step calculator for fixed-step symplectic integration
Leapfrog_StepCalculator leapfrog;
ODESystemFixedStepSolver solver(system, leapfrog);

// System must be: y = [x, v]ᵀ  (position, velocity)
Vector<Real> y0({0.0, 1.0});  // x=0, v=1
ODESystemSolution sol = solver.integrate(y0, 0.0, 10.0, 1000);
```

**Properties**:
- **Order**: 2nd order
- **Conservation**: Preserves energy (symplectic)
- **Use Case**: Molecular dynamics, celestial mechanics
- **Advantage**: Long-time stability

**Example: Simple Harmonic Oscillator**

```cpp
// ẍ + ω²x = 0  →  ẋ = v, v̇ = -ω²x
ODESystem sho(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    Real omega = 1.0;
    dydt[0] = y[1];                    // ẋ = v
    dydt[1] = -omega*omega * y[0];     // v̇ = -ω²x
});

Leapfrog_StepCalculator leapfrog;
ODESystemFixedStepSolver solver(sho, leapfrog);
Vector<Real> y0({1.0, 0.0});  // x=1, v=0
ODESystemSolution sol = solver.integrate(y0, 0.0, 100.0, 1000);

// Energy E = ½v² + ½ω²x² remains constant!
```

---

### ODE Solver Selection Guide

| Method | Order | When to Use |
|--------|-------|-------------|
| **Euler** | 1 | Quick demos, educational |
| **Midpoint** | 2 | Better than Euler, simple |
| **RK4** | 4 | General purpose, fixed step |
| **Cash-Karp** | 5(4) | Adaptive, general purpose |
| **Dormand-Prince 5** | 5(4) | Adaptive with FSAL, dense output |
| **Dormand-Prince 8** | 8(7) | High precision requirements |
| **Leapfrog** | 2 | Hamiltonian, energy conservation |

**Recommendations**:
- **Default**: Dormand-Prince 5 for most problems (FSAL, dense output)
- **Stiff systems**: Consider implicit methods (future)
- **Long-time**: Leapfrog for Hamiltonian
- **High precision**: Dormand-Prince 8 or RK4 with small step

---

## Optimization

**Files**: `mml/algorithms/Optimization.h`, `mml/algorithms/Optimization/OptimizationMultidim.h`

**Problem**: Find **x*** that minimizes f(**x**)

### 1D Minimization

**Purpose**: Find minimum of a function in one dimension.

```cpp
#include "algorithms/Optimization.h"

// Minimize f(x) = (x-2)² + 1
RealFunctionFromStdFunc f([](Real x) { return (x-2)*(x-2) + 1; });

// First bracket the minimum
MinimumBracket bracket = Minimization::BracketMinimum(f, 0.0, 1.0);

// Then refine with golden section or Brent
MinimizationResult result = Minimization::GoldenSectionSearch(f, bracket);
std::cout << "x_min = " << result.xmin << ", f_min = " << result.fmin << std::endl;
// x_min ≈ 2.0, f_min ≈ 1.0
```

#### Golden Section Search

**Algorithm**: Fibonacci-based interval reduction.

```cpp
// Convenience function: bracket + minimize in one call
MinimizationResult result = Minimization::GoldenSectionSearch(f, 0.0, 10.0, 1e-6);
Real x_min = result.xmin;
Real f_min = result.fmin;
```

**Properties**:
- **Convergence**: Linear with ratio 0.618
- **Ratio**: Golden ratio φ = 1.618...
- **No derivatives needed**
- **Robust**: Always converges in bracket

#### Brent's Method

**Algorithm**: Combines golden section with parabolic interpolation.

```cpp
MinimizationResult result = Minimization::BrentMinimize(f, bracket, 1e-6);
// Or with auto-bracketing:
MinimizationResult result = Minimization::BrentMinimize(f, 0.0, 10.0, 1e-6);
```

**Properties**:
- **Convergence**: Superlinear to quadratic
- **Speed**: Faster than golden section for smooth functions
- **Robustness**: Falls back to golden section if needed

---

### Multidimensional Optimization

**File**: `mml/algorithms/Optimization/OptimizationMultidim.h`

#### Powell's Method (Derivative-Free)

**Algorithm**: Direction set method that doesn't require derivatives.

```cpp
#include "algorithms/Optimization/OptimizationMultidim.h"

ScalarFunction<2> f([](const VectorN<Real, 2>& v) {
    Real x = v[0], y = v[1];
    return x*x + 4*y*y;  // Elliptical paraboloid
});

VectorN<Real, 2> x0{1.0, 1.0};
MultidimMinimizationResult result = PowellMinimize(f, x0, 1e-6);
// result.xmin ≈ [0, 0]
```

**Properties**:
- **Convergence**: Quadratic for quadratic functions
- **Speed**: Good for smooth functions
- **Pro**: No derivatives required
- **Con**: Can be slow for ill-conditioned problems

#### Conjugate Gradient

**Algorithm**: Use conjugate directions to avoid zigzagging.

```cpp
ConjugateGradient cg(1e-8, 1e-8, 200, ConjugateGradient::Method::PolakRibiere);
MultidimMinimizationResult result = cg.Minimize(f, x0);
// Or use convenience function:
MultidimMinimizationResult result = ConjugateGradientMinimize(f, x0, 1e-6);
// Converges in ≤ n iterations for quadratic f (n = dimension)
```

**Properties**:
- **Convergence**: Superlinear
- **Speed**: Much faster than steepest descent
- **Ideal for**: Quadratic functions, large dimensions
- **Memory**: O(n) storage

**Comparison**:
- Steepest descent: 50+ iterations
- Conjugate gradient: 2 iterations (for 2D quadratic)

---

### Quasi-Newton Methods

#### BFGS (Broyden-Fletcher-Goldfarb-Shanno)

**Algorithm**: Approximate Hessian using gradient history.

```cpp
BFGS bfgs(1e-8, 1e-8, 200);
MultidimMinimizationResult result = bfgs.Minimize(f, x0);
// Or use convenience function:
MultidimMinimizationResult result = BFGSMinimize(f, x0, 1e-6);
```

**Properties**:
- **Convergence**: Superlinear to quadratic
- **Speed**: Fastest gradient-based method
- **Memory**: O(n²) for Hessian approximation
- **Pro**: Quadratic convergence without computing Hessian

**Update Formula**:
```
B_{k+1} = B_k + (y_k y_k^T)/(y_k^T s_k) - (B_k s_k s_k^T B_k)/(s_k^T B_k s_k)
where s_k = x_{k+1} - x_k, y_k = ∇f_{k+1} - ∇f_k
```

---

### Constrained Optimization

*Note: Full constrained optimization methods (penalty method, Lagrange multipliers) are planned for future versions. Currently, constraints can be handled by transforming the problem or using box constraints with simulated annealing.*

---

## Heuristic Optimization

**File**: `mml/algorithms/Optimization/SimulatedAnnealing.h`

**Purpose**: Global optimization for multimodal and non-convex problems.

### Simulated Annealing

**Algorithm**: Probabilistic method inspired by metallurgy.

```cpp
#include "algorithms/Optimization/SimulatedAnnealing.h"

// Find global minimum of Rastrigin function (many local minima)
auto rastrigin = [](const Vector<Real>& v) {
    Real x = v[0], y = v[1];
    Real A = 10.0;
    return 2*A + (x*x - A*std::cos(2*Constants::PI*x)) 
               + (y*y - A*std::cos(2*Constants::PI*y));
};

Vector<Real> x0({4.0, 4.0});

// Configure simulated annealing
SimulatedAnnealing sa(
    100.0,   // Initial temperature T0
    0.95,    // Cooling rate alpha
    10000,   // Max iterations
    0.5      // Step size
);

// Run optimization
HeuristicOptimizationResult result = sa.Minimize(rastrigin, x0);
std::cout << "Best x: " << result.xbest << std::endl;
std::cout << "Best f: " << result.fbest << std::endl;
std::cout << "Iterations: " << result.iterations << std::endl;
std::cout << "Accepted moves: " << result.acceptedMoves << std::endl;
// x_best ≈ [0, 0] (global minimum)

// Or use convenience function:
HeuristicOptimizationResult result = SimulatedAnnealingMinimize(
    rastrigin, x0, 100.0, 0.95, 10000, 0.5
);
```

**Algorithm**:
1. Generate random neighbor within step size
2. If better, accept
3. If worse, accept with probability exp(-ΔE/T)
4. Decrease temperature according to cooling schedule

**Cooling Schedules**:
- **Exponential**: T(k) = T₀ × αᵏ (most common)
- **Linear**: T(k) = T₀ - k × (T₀ - T_final) / max_iter
- **Logarithmic**: T(k) = T₀ / (1 + log(1 + k))

```cpp
// Custom cooling schedule
ExponentialCooling cooling(0.95);
SimulatedAnnealing sa(100.0, cooling, 10000, 0.5);
```

**Properties**:
- **Global**: Can escape local minima
- **Stochastic**: Different runs give different results
- **No derivatives**: Works for non-smooth functions
- **Flexible**: Works with bounds constraints

**Tuning**:
- **T₀**: Start high enough to accept ~80% of uphill moves
- **α**: Typical 0.9-0.99 (slower cooling = better results)
- **Step size**: Problem-dependent, affects neighbor generation

**With Box Constraints**:
```cpp
Vector<Real> lowerBounds({-5.12, -5.12});
Vector<Real> upperBounds({5.12, 5.12});

HeuristicOptimizationResult result = SimulatedAnnealingMinimize(
    rastrigin, x0, 100.0, 0.95, 10000, 0.5, lowerBounds, upperBounds
);
```

---

### Future Heuristic Methods

Planned for future versions:
- **Genetic Algorithms**: Population-based evolutionary optimization
- **Particle Swarm Optimization**: Swarm intelligence method
- **Differential Evolution**: Robust global optimizer

---

## Function Analysis

**File**: `mml/algorithms/FunctionsAnalyzer.h`

**Purpose**: Analyze properties of functions.

### RealFunctionAnalyzer Class

The `RealFunctionAnalyzer` class provides comprehensive analysis of real-valued functions.

```cpp
#include "algorithms/FunctionsAnalyzer.h"

RealFunctionFromStdFunc f([](Real x) { return x*x*x - 3*x; });
RealFunctionAnalyzer analyzer(f);
```

### Finding Roots

```cpp
// Find all roots in an interval
std::vector<Real> roots = analyzer.GetRoots(-5.0, 5.0, 1e-6);
// roots = {-√3, 0, √3}  (where x³ - 3x = 0)
```

### Finding Local Optima

```cpp
// Find local optimums (where f'(x) = 0 and it's a true min/max)
Vector<Real> optimums = analyzer.GetLocalOptimums(-5.0, 5.0, 1e-6);
// optimums contains x values where there are local min/max

// Get classified critical points with type information
std::vector<CriticalPoint> classified = analyzer.GetLocalOptimumsClassified(-5.0, 5.0);
for (const auto& cp : classified) {
    std::cout << "x = " << cp.x << ", f(x) = " << cp.value;
    switch (cp.type) {
        case CriticalPointType::LOCAL_MINIMUM:
            std::cout << " (local minimum)";
            break;
        case CriticalPointType::LOCAL_MAXIMUM:
            std::cout << " (local maximum)";
            break;
        case CriticalPointType::SADDLE_POINT:
            std::cout << " (inflection point)";
            break;
    }
    std::cout << std::endl;
}
```

### Finding Inflection Points

```cpp
// Find inflection points (where f''(x) changes sign)
Vector<Real> inflections = analyzer.GetInflectionPoints(-5.0, 5.0, 1e-6);
```

### Point Analysis

```cpp
// Analyze function properties at a specific point
Real x = 1.0;
bool defined = analyzer.isDefinedAtPoint(x);
bool continuous = analyzer.isContinuousAtPoint(x, 1e-6);
bool derivDefined = analyzer.isDerivativeDefinedAtPoint(x, 1e-6);
bool isLocalOpt = analyzer.isLocalOptimum(x, 1e-6);
bool isInflection = analyzer.isInflectionPoint(x, 1e-6);

// Print comprehensive analysis
analyzer.PrintPointAnalysis(x);
```

### Interval Analysis

```cpp
// Analyze function over an interval
analyzer.PrintIntervalAnalysis(-5.0, 5.0, 100);
// Prints: defined regions, continuity, critical points, inflections
```

---

## Statistics

**Files**: `mml/algorithms/Statistics.h`, `mml/algorithms/Statistics/*.h`

**Purpose**: Statistical analysis of data.

### Descriptive Statistics

```cpp
#include "algorithms/Statistics.h"

Vector<Real> data({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});

// Central tendency
Real mean = Statistics::Mean(data);           // 5.5  (alias for Avg)
Real avg = Statistics::Avg(data);             // 5.5
Real median = Statistics::Median(data);       // 5.5

// Dispersion
Real variance = Statistics::Variance(data);   // 9.17
Real stddev = Statistics::StdDev(data);       // 3.03

// Combined computation (more efficient)
Real outAvg, outVar;
Statistics::AvgVar(data, outAvg, outVar);     // Computes both in one pass

Real outAvg2, outStdDev;
Statistics::AvgStdDev(data, outAvg2, outStdDev);

// Full moments analysis
Real ave, adev, sdev, var, skew, curt;
Statistics::Moments(data, ave, adev, sdev, var, skew, curt);
// ave = mean, adev = mean absolute deviation
// sdev = standard deviation, var = variance
// skew = skewness (0 = symmetric)
// curt = kurtosis (-1.2 for uniform-ish, 0 for normal)
```

### Order Statistics

```cpp
// Median (middle value)
Real median = Statistics::Median(data);           // 5.5

// Percentile (0-100 scale)
Real p25 = Statistics::Percentile(data, 25.0);   // 25th percentile (Q1)
Real p75 = Statistics::Percentile(data, 75.0);   // 75th percentile (Q3)
Real iqr = p75 - p25;                             // Interquartile range

// Or use dedicated IQR function
Real iqr2 = Statistics::IQR(data);

// Full quartiles at once
Real q1, q2, q3;
Statistics::Quartiles(data, q1, q2, q3);          // Q1=25%, Q2=50%(median), Q3=75%

// Range and MinMax
Real range = Statistics::Range(data);             // max - min
Real minVal, maxVal;
Statistics::MinMax(data, minVal, maxVal);
```

### Robust Statistics

```cpp
// Mode (most frequent value)
Real mode = Statistics::Mode(data);

// Trimmed mean (exclude extreme values)
Real trimmed = Statistics::TrimmedMean(data, 10.0);  // Trim 10% from each end

// Median Absolute Deviation (robust alternative to stddev)
Real mad = Statistics::MAD(data);                    // mad * 1.4826 ≈ stddev for normal

// Alternative means
Real gmean = Statistics::GeometricMean(data);        // For multiplicative data
Real hmean = Statistics::HarmonicMean(data);         // For rates/ratios

// Weighted statistics
Vector<Real> weights({1, 1, 2, 2, 3, 3, 2, 2, 1, 1});
Real wmean = Statistics::WeightedMean(data, weights);
Real wvar = Statistics::WeightedVariance(data, weights);
```

### Correlation and Covariance

```cpp
Vector<Real> x({1, 2, 3, 4, 5});
Vector<Real> y({2, 4, 5, 4, 5});

// Covariance
Real cov = Statistics::Covariance(x, y);

// Pearson correlation coefficient
Real r = Statistics::PearsonCorrelation(x, y);      // r ∈ [-1, 1]
// 1 = perfect positive, -1 = perfect negative, 0 = no linear correlation

// Correlation with significance testing
Statistics::CorrelationResult result = Statistics::PearsonCorrelationWithTest(x, y);
// result.r          - correlation coefficient
// result.pValue     - p-value for testing r=0
// result.tStatistic - t-statistic for significance

// Coefficient of determination
Real r_squared = Statistics::RSquared(x, y);        // R² = r²

// For multivariate data (n observations × p variables)
Matrix<Real> multiData(100, 5);  // 100 observations, 5 variables
Matrix<Real> covMatrix = Statistics::CovarianceMatrix(multiData);
Matrix<Real> corrMatrix = Statistics::CorrelationMatrix(multiData);
```

### Rank Correlation

```cpp
#include "algorithms/Statistics/RankCorrelation.h"

Vector<Real> x({1, 2, 3, 4, 5});
Vector<Real> y({5, 6, 7, 8, 7});

// Spearman rank correlation (non-parametric)
Real rho = Statistics::SpearmanCorrelation(x, y);

// Kendall's tau (ordinal association)
Real tau = Statistics::KendallTau(x, y);

// With significance testing
Statistics::CorrelationResult spearman = Statistics::SpearmanCorrelationWithTest(x, y);
Statistics::CorrelationResult kendall = Statistics::KendallTauWithTest(x, y);
```

### Distributions

**Files**: `mml/algorithms/Statistics/CoreDistributions.h`, `mml/algorithms/Statistics/Distributions.h`

#### Normal Distribution

```cpp
#include "algorithms/Statistics/CoreDistributions.h"

// Create standard normal (μ=0, σ=1)
Statistics::NormalDistribution stdNormal;             // Default: μ=0, σ=1
Statistics::NormalDistribution normal(5.0, 2.0);      // μ=5, σ=2

// PDF: Probability density function
Real pdf = normal.pdf(6.0);                           // f(6) = 0.176

// CDF: Cumulative distribution function P(X ≤ x)
Real cdf = normal.cdf(5.0);                           // P(X ≤ 5) = 0.5

// Inverse CDF (quantile function) - find x for given probability
Real x = normal.inverseCdf(0.975);                    // 95% upper bound

// Z-score (standardize a value)
Real z = normal.zScore(7.0);                          // (7-5)/2 = 1.0
```

#### Student's t-Distribution

```cpp
// For small samples and unknown population variance
Statistics::TDistribution t_dist(10);                 // 10 degrees of freedom

Real pdf = t_dist.pdf(1.5);                           // t-density at 1.5
Real cdf = t_dist.cdf(2.0);                           // P(T ≤ 2.0)
Real critical = t_dist.inverseCdf(0.975);             // Critical value for 95% CI
```

#### Chi-Square Distribution

```cpp
Statistics::ChiSquareDistribution chi2(5);            // 5 degrees of freedom

Real pdf = chi2.pdf(4.0);                             // χ² density at 4.0
Real cdf = chi2.cdf(11.07);                           // P(χ² ≤ 11.07) ≈ 0.95
Real critical = chi2.inverseCdf(0.95);                // 95th percentile
```

#### F-Distribution

```cpp
Statistics::FDistribution F(5, 10);                   // df1=5, df2=10

Real pdf = F.pdf(2.0);                                // F-density at 2.0
Real cdf = F.cdf(3.33);                               // P(F ≤ 3.33)
Real critical = F.inverseCdf(0.95);                   // Critical value for ANOVA
```

#### Other Distributions

```cpp
#include "algorithms/Statistics/Distributions.h"

// Cauchy (heavy-tailed, no mean/variance)
Statistics::CauchyDistribution cauchy(0.0, 1.0);      // location=0, scale=1

// Exponential (waiting times, lifetimes)
Statistics::ExponentialDistribution expo(2.0);        // rate λ=2 (mean=0.5)

// Logistic (S-shaped curve, logistic regression)
Statistics::LogisticDistribution logistic(0.0, 1.0);  // location=0, scale=1
```

### Hypothesis Testing

**File**: `mml/algorithms/Statistics/HypothesisTesting.h`

#### T-Tests

```cpp
#include "algorithms/Statistics/HypothesisTesting.h"

// One-sample t-test: H₀: μ = μ₀
Real mu0 = 5.0;
Real alpha = 0.05;  // 95% confidence

Statistics::HypothesisTestResult result = Statistics::OneSampleTTest(data, mu0, alpha);
// result.testStatistic  - t-statistic
// result.pValue         - p-value
// result.criticalValue  - critical value at alpha
// result.rejectNull     - true if should reject H₀
// result.confidenceLevel - 1 - alpha (0.95)
// result.degreesOfFreedom
// result.testName       - "One-Sample t-Test"

if (result.rejectNull) {
    // Reject null hypothesis at 95% confidence
}

// Two-sample t-test: H₀: μ₁ = μ₂
Vector<Real> sample1({...});
Vector<Real> sample2({...});
Statistics::HypothesisTestResult twoSample = 
    Statistics::TwoSampleTTest(sample1, sample2, alpha);

// Paired t-test (before/after measurements)
Vector<Real> before({...});
Vector<Real> after({...});
Statistics::HypothesisTestResult paired = 
    Statistics::PairedTTest(before, after, alpha);

// Welch's t-test (unequal variances)
Statistics::HypothesisTestResult welch = 
    Statistics::WelchTTest(sample1, sample2, alpha);
```

#### One-Way ANOVA

```cpp
// One-way ANOVA: Compare means of multiple groups
std::vector<Vector<Real>> groups = {group1, group2, group3};
Statistics::HypothesisTestResult anova = Statistics::OneWayANOVA(groups, alpha);

if (anova.rejectNull) {
    // At least one group mean differs significantly
}
```

---

### Confidence Intervals

**File**: `mml/algorithms/Statistics/ConfidenceIntervals.h`

```cpp
#include "algorithms/Statistics/ConfidenceIntervals.h"

// Confidence interval for mean (known variance - Z-interval)
Real knownSigma = 2.0;
Real confidence = 0.95;
Statistics::ConfidenceInterval zInterval = 
    Statistics::MeanConfidenceInterval_Z(data, knownSigma, confidence);
// zInterval.lower, zInterval.upper, zInterval.center

// Confidence interval for mean (unknown variance - t-interval)
Statistics::ConfidenceInterval tInterval = 
    Statistics::MeanConfidenceInterval_T(data, confidence);

// Confidence interval for proportion
int successes = 75, trials = 100;
Statistics::ConfidenceInterval propInterval = 
    Statistics::ProportionConfidenceInterval(successes, trials, confidence);

// Confidence interval for variance (chi-square based)
Statistics::ConfidenceInterval varInterval = 
    Statistics::VarianceConfidenceInterval(data, confidence);
```

---

### Time Series Analysis

**File**: `mml/algorithms/Statistics/TimeSeries.h`

```cpp
#include "algorithms/Statistics/TimeSeries.h"

Vector<Real> series({...});  // Time series data

// Autocorrelation at lag k
Real acf_1 = Statistics::Autocorrelation(series, 1);   // Lag 1
Real acf_5 = Statistics::Autocorrelation(series, 5);   // Lag 5

// Partial autocorrelation
Real pacf_1 = Statistics::PartialAutocorrelation(series, 1);

// Full ACF up to maxLag
int maxLag = 20;
Vector<Real> acf = Statistics::AutocorrelationFunction(series, maxLag);
```

---

## Algorithm Selection Guide

### Root Finding

| Function Type | Recommended Method | Reason |
|---------------|-------------------|---------|
| Smooth, bracketed | Newton-Raphson | Fast convergence |
| Non-smooth | Bisection | Robustness |
| No derivative | Secant | Good compromise |
| Multiple roots | Bracket + Bisection | Find all roots |
| Multidimensional | Newton system | Quadratic convergence |

### ODE Solving

| Problem Type | Recommended Solver | Reason |
|--------------|-------------------|---------|
| General, smooth | RKCK adaptive | Efficiency + accuracy |
| Stiff system | Implicit method* | Stability |
| Hamiltonian | Leapfrog | Energy conservation |
| High precision | RK4 small step | Accuracy |
| Quick demo | Euler | Simplicity |

*Note: Stiff solvers (implicit methods) not yet implemented

### Optimization

| Problem Type | Recommended Method | Reason |
|--------------|-------------------|--------|
| Smooth, gradient available | BFGS | Quasi-Newton, fast convergence |
| Quadratic-like | ConjugateGradient | Optimal for quadratic forms |
| Derivative-free | PowellMinimize | No gradient required |
| Non-smooth/many local minima | SimulatedAnnealing | Global search |
| 1D optimization | BrentMinimize | Robust, efficient |
| Initial bracketing | BracketMinimum | Find containing interval |

---

## Numerical Considerations

### Accuracy vs. Speed Trade-offs

```cpp
// High accuracy (slow)
ODEAdaptiveIntegrator<DormandPrince8_Stepper> precise(system);
auto sol_accurate = precise.Integrate(y0, t0, tf, 1e-12, 1e-12);

// Moderate accuracy (balanced)
ODEAdaptiveIntegrator<DormandPrince5_Stepper> balanced(system);
auto sol_moderate = balanced.Integrate(y0, t0, tf, 1e-6, 1e-6);

// Low accuracy (fast)
ODEAdaptiveIntegrator<CashKarp_Stepper> fast(system);
auto sol_fast = fast.Integrate(y0, t0, tf, 1e-3, 1e-3);
```

### Ill-Conditioned Problems

```cpp
// Problem: Nearly singular Hessian in optimization
// Solution: Add regularization
auto regularized_f = [&](const VectorN<Real, n>& x) {
    Real lambda = 1e-6;  // Regularization parameter
    return objective(x) + lambda * x.NormL2();
};
```

### Stiff ODE Systems

Systems where explicit methods require very small steps:
```cpp
// Characteristic: Rapidly decaying transients + slow dynamics
// Example: Chemical kinetics, electrical circuits

// Future: Implicit methods (BDF, Radau)
// Current workaround: Use small fixed step or very strict tolerance
```

---

## Error Handling

```cpp
try {
    Real root = RootFinding::FindRootNewton(f, x1, x2, tol);
} catch (const RootFindingError& e) {
    // Derivative zero, divergence, or max iterations
    // Fall back to bisection
    root = RootFinding::FindRootBisection(f, x1, x2, tol);
}

try {
    Minimization::BFGS optimizer(f, grad);
    auto result = optimizer.Minimize(x0);
} catch (const OptimizationError& e) {
    // Non-descent direction or line search failure
    // Try different initial guess or use global method
    SimulatedAnnealing sa(...);
    auto global_result = sa.Optimize(bounds);
}
```

---

## Performance Tips

1. **Choose appropriate method** for problem structure
2. **Use adaptive methods** when accuracy requirements vary
3. **Vectorize** when possible (batch operations)
4. **Warm start** iterative methods with good initial guess
5. **Cache** expensive function evaluations
6. **Use cooling schedules wisely** in simulated annealing

---

## Testing

Run algorithm tests:
```bash
cd build
.\tests\Debug\MML_Tests.exe "[algorithms]"
.\tests\Debug\MML_Tests.exe "[root-finding]"
.\tests\Debug\MML_Tests.exe "[ode-solvers]"
.\tests\Debug\MML_Tests.exe "[optimization]"
```

---

## Summary

The Algorithms layer provides:
- **Root Finding**: 5 methods (bisection, Newton-Raphson, secant, Ridder, Brent)
- **ODE Solvers**: Fixed-step (Euler, Midpoint, RK4, Leapfrog) + Adaptive (Cash-Karp, DormandPrince5/8)
- **Optimization 1D**: Bracketing + Golden section + Brent's method
- **Optimization Multidim**: Conjugate gradient, BFGS, Powell
- **Heuristic Optimization**: Simulated Annealing with cooling schedules
- **Function Analysis**: Root finding, critical points, inflection points
- **Statistics**: Descriptive stats, order stats, robust stats, correlation, distributions, hypothesis testing

All algorithms:
✅ Production-ready implementations  
✅ Well-tested and validated  
✅ Configurable parameters  
✅ Error handling and robustness  
✅ Documented complexity and convergence  

**Previous Layer**: [Core Layer Documentation](README_Core.md)  
**Next Layer**: [Testbeds Layer Documentation](README_Testbeds.md)

---

*Last Updated: December 27, 2025*  
*MML Version: 1.0*
