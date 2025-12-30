# Function Optimization

Comprehensive toolkit for **unconstrained optimization** - from univariate golden section to sophisticated multidimensional quasi-Newton methods.

## Overview

**Optimization** finds minima (or maxima) of functions: given f(x), find x* such that f(x*) ≤ f(x) for all nearby x.

**Applications**:
- **Machine Learning**: Training neural networks, logistic regression, SVM
- **Engineering**: Design optimization, parameter fitting, control systems
- **Physics**: Energy minimization, least action principle, equilibrium states
- **Finance**: Portfolio optimization, risk minimization, model calibration
- **Statistics**: Maximum likelihood estimation, curve fitting

**Library Architecture**:
- **1D Methods**: Bracketing + refinement (Golden Section, Brent)
- **Derivative-Free Multidim**: Nelder-Mead (simplex), Powell (direction set)
- **Gradient-Based Multidim**: Conjugate Gradient, BFGS quasi-Newton
- **Configuration**: Fluent builder API for termination criteria and observers

## Quick Reference

### 1D Minimization (Optimization.h)

| Method | Convergence | Derivatives | Best For |
|--------|-------------|-------------|----------|
| **BracketMinimum** | N/A | No | Finding initial bracket |
| **GoldenSectionSearch** | Linear (0.618) | No | Robust, simple |
| **BrentMinimize** | Superlinear | No | **General purpose** |
| **BrentMinimizeWithDeriv** | Superlinear | Yes | Fast when derivatives cheap |

### Multidimensional Minimization (OptimizationMultidim.h)

| Method | Derivatives | Convergence | Storage | Best For |
|--------|-------------|-------------|---------|----------|
| **NelderMead** | No | Linear | O(N²) | **Black-box functions** |
| **Powell** | No | Superlinear | O(N²) | Smooth functions, moderate N |
| **ConjugateGradient** | Yes | Superlinear | O(N) | **Large-scale smooth** |
| **BFGS** | Yes | Superlinear | O(N²) | **General smooth, N < 1000** |

### Convergence Rates

| Method | Rate | Iterations to ε |
|--------|------|-----------------|
| Golden Section | 0.618 per iter | ~40 for 10⁻⁸ |
| Bisection | 0.5 per iter | ~27 for 10⁻⁸ |
| Brent (parabolic) | ~1.324 | ~15-20 for 10⁻⁸ |
| Newton (1D) | Quadratic | ~4-5 for 10⁻⁸ |
| CG/BFGS | Superlinear | ~10-50 (problem dependent) |

---

## 1D Optimization

### BracketMinimum - Find Initial Bracket

**Purpose**: Locate an interval [a, c] guaranteed to contain a minimum.

**Algorithm**: Geometric expansion using golden ratio (φ ≈ 1.618) with parabolic acceleration.

```cpp
#include "algorithms/Optimization.h"
using namespace MML;

// Function to minimize: f(x) = (x - 2)² + 1
RealFunctionFromStdFunc func([](Real x) { 
    return (x - 2.0) * (x - 2.0) + 1.0; 
});

// Find bracket starting from [0, 1]
MinimumBracket bracket = Minimization::BracketMinimum(func, 0.0, 1.0);

if (bracket.valid) {
    std::cout << "Bracket: [" << bracket.ax << ", " << bracket.cx << "]\n";
    std::cout << "Middle point: " << bracket.bx << ", f(bx) = " << bracket.fb << "\n";
}
// Output: Bracket contains x=2 with fb < fa and fb < fc
```

**Postcondition**: `bx` is between `ax` and `cx`, and `fb ≤ min(fa, fc)`.

### GoldenSectionSearch - Robust 1D Minimization

**Algorithm**: Divide interval by golden ratio, reusing one evaluation per iteration.

**Convergence**: Linear with ratio 0.618 - guaranteed for any unimodal function.

```cpp
// Golden section search (robust, derivative-free)
auto result = Minimization::GoldenSectionSearch(func, bracket);

std::cout << "Minimum at x = " << result.xmin << "\n";      // ≈ 2.0
std::cout << "f(x*) = " << result.fmin << "\n";             // ≈ 1.0
std::cout << "Iterations: " << result.iterations << "\n";   // ~40

// Or use convenience function (brackets automatically)
auto result2 = Minimization::GoldenSectionSearch(func, 0.0, 5.0);
```

**When to Use**:
- Function evaluation is expensive
- Derivatives not available
- Robustness more important than speed
- Function may not be smooth

### BrentMinimize - Fast 1D Minimization

**Algorithm**: Combines golden section (safe) with parabolic interpolation (fast).

**Convergence**: Superlinear (~1.324) for smooth functions, degrades gracefully.

```cpp
// Brent's method (preferred for general use)
auto result = Minimization::BrentMinimize(func, bracket);
// OR with automatic bracketing:
auto result = Minimization::BrentMinimize(func, 0.0, 5.0);

std::cout << "Minimum: x* = " << result.xmin << "\n";
std::cout << "f(x*) = " << result.fmin << "\n";
std::cout << "Converged: " << result.converged << "\n";
```

**Why Brent is Preferred**:
- Near-optimal convergence for smooth functions
- Never slower than golden section
- Handles non-smooth functions gracefully
- Industry-standard algorithm

### BrentMinimizeWithDeriv - When Derivatives Available

**Algorithm**: Brent's method using derivative sign to narrow bracket + secant extrapolation.

```cpp
// Function and its derivative
RealFunctionFromStdFunc func([](Real x) { return x*x - 4*x + 5; });
RealFunctionFromStdFunc dfunc([](Real x) { return 2*x - 4; });

// Brent with derivatives (faster convergence)
auto result = Minimization::BrentMinimizeWithDeriv(func, dfunc, 0.0, 5.0);
// Minimum at x* = 2.0
```

**When to Use**: Derivative is cheap (closed-form, automatic differentiation).

### Maximization Wrappers

```cpp
// Find maximum (internally negates function)
auto maxResult = Minimization::BrentMaximize(func, a, b);
std::cout << "Maximum at x = " << maxResult.xmin << "\n";
std::cout << "f(x) = " << maxResult.fmin << "\n";  // Actual max value (not negated)
```

---

## Multidimensional Optimization

### Result Structure

All multidimensional methods return `MultidimMinimizationResult`:

```cpp
struct MultidimMinimizationResult {
    Vector<Real> xmin;      // Location of minimum
    Real         fmin;      // Function value at minimum
    int          iterations; // Iterations or function evaluations
    bool         converged;  // True if converged within tolerance
};
```

### NelderMead - Derivative-Free Simplex Method

**Algorithm**: The "amoeba" method - a simplex (N+1 points in N dimensions) that adapts through reflection, expansion, contraction, and shrinking.

**Convergence**: Linear but very robust. Works on non-smooth, noisy functions.

```cpp
#include "algorithms/Optimization/OptimizationMultidim.h"
using namespace MML;

// Rosenbrock function: f(x,y) = (1-x)² + 100(y-x²)²
class Rosenbrock : public IScalarFunction<2> {
public:
    Real operator()(const VectorN<Real, 2>& x) const override {
        Real a = 1.0 - x[0];
        Real b = x[1] - x[0] * x[0];
        return a * a + 100.0 * b * b;
    }
};

Rosenbrock func;
VectorN<Real, 2> start = {-1.0, 1.0};

// Method 1: Simple call with uniform delta
auto result = NelderMeadMinimize(func, start, /*delta=*/1.0);

// Method 2: Using optimizer object for more control
NelderMead optimizer(/*ftol=*/1e-8, /*maxIter=*/5000);
auto result2 = optimizer.Minimize(func, start, /*delta=*/1.0);

std::cout << "Minimum at: (" << result.xmin[0] << ", " << result.xmin[1] << ")\n";
std::cout << "f(x*) = " << result.fmin << "\n";
std::cout << "Function evals: " << result.iterations << "\n";
// Output: (1.0, 1.0), f = 0, ~500-1000 evals
```

**Advanced: Custom Initial Simplex**

```cpp
// Per-dimension deltas (scale each axis differently)
Vector<Real> deltas(2);
deltas[0] = 0.5;   // x-axis perturbation
deltas[1] = 2.0;   // y-axis perturbation

auto result = optimizer.Minimize(func, start, deltas);

// Or provide explicit simplex (3 points for 2D)
Matrix<Real> simplex(3, 2);
simplex(0, 0) = 0.0;  simplex(0, 1) = 0.0;   // Vertex 0
simplex(1, 0) = 1.0;  simplex(1, 1) = 0.0;   // Vertex 1
simplex(2, 0) = 0.5;  simplex(2, 1) = 1.0;   // Vertex 2

auto result2 = optimizer.Minimize(func, simplex);
```

**When to Use NelderMead**:
- ✅ Derivatives unavailable or expensive
- ✅ Non-smooth or noisy functions
- ✅ Black-box optimization
- ✅ Moderate dimensions (N < 10-20)
- ❌ Large-scale problems (slow convergence)
- ❌ High-precision requirements

### Powell - Direction Set Method

**Algorithm**: Successive 1D line minimizations along directions that become conjugate.

**Convergence**: Superlinear for quadratic-like functions.

```cpp
// Powell's method - no derivatives needed
Powell optimizer(/*ftol=*/3e-8, /*maxIter=*/200);
auto result = optimizer.Minimize(func, start);

// Or convenience wrapper
auto result2 = PowellMinimize(func, start, /*ftol=*/1e-8);

std::cout << "Powell converged: " << result.converged << "\n";
std::cout << "Iterations: " << result.iterations << "\n";
```

**Advanced: Custom Initial Directions**

```cpp
// Initialize with custom direction matrix (columns are directions)
Matrix<Real> directions(2, 2, 0.0);
directions(0, 0) = 1.0;  // First direction: (1, 0)
directions(1, 1) = 1.0;  // Second direction: (0, 1)

auto result = optimizer.Minimize(func, start, directions);
```

**When to Use Powell**:
- ✅ Smooth functions, derivatives unavailable
- ✅ Moderate dimensions (N < 20)
- ✅ Better than NelderMead for smooth functions
- ❌ Non-smooth functions (use NelderMead)

---

## Gradient-Based Methods

### Differentiable Function Interface

Gradient-based methods require implementing `IDifferentiableScalarFunction<N>`:

```cpp
template<int N>
class MyFunction : public IDifferentiableScalarFunction<N> {
public:
    Real operator()(const VectorN<Real, N>& x) const override {
        // Compute function value
    }
    
    void Gradient(const VectorN<Real, N>& x, VectorN<Real, N>& grad) const override {
        // Fill grad with ∂f/∂xᵢ
    }
};
```

**Example: Rosenbrock with Gradient**

```cpp
class RosenbrockWithGrad : public IDifferentiableScalarFunction<2> {
public:
    Real operator()(const VectorN<Real, 2>& x) const override {
        Real a = 1.0 - x[0];
        Real b = x[1] - x[0] * x[0];
        return a * a + 100.0 * b * b;
    }
    
    void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& g) const override {
        // ∂f/∂x = -2(1-x) - 400x(y-x²)
        // ∂f/∂y = 200(y-x²)
        g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
        g[1] = 200.0 * (x[1] - x[0] * x[0]);
    }
};
```

### ConjugateGradient - Efficient Large-Scale Optimization

**Algorithm**: Generate search directions that are conjugate (H-orthogonal) to previous directions.

**Variants**:
- **Fletcher-Reeves**: γ = (g_new · g_new) / (g_old · g_old)
- **Polak-Ribière**: γ = (g_new · (g_new - g_old)) / (g_old · g_old) — **recommended**

**Convergence**: Superlinear, theoretically N iterations for quadratic.

```cpp
RosenbrockWithGrad func;
VectorN<Real, 2> start = {-1.0, 1.0};

// Conjugate gradient with Polak-Ribière (default)
ConjugateGradient cg(/*ftol=*/3e-8, /*gtol=*/1e-8, /*maxIter=*/200);
auto result = cg.Minimize(func, start);

// Use Fletcher-Reeves variant
cg.setMethod(ConjugateGradient::Method::FletcherReeves);
auto result2 = cg.Minimize(func, start);

// Or convenience wrapper
auto result3 = ConjugateGradientMinimize(func, start, 1e-8, 
                                          ConjugateGradient::Method::PolakRibiere);

std::cout << "CG converged in " << result.iterations << " iterations\n";
```

**When to Use CG**:
- ✅ Large-scale problems (N > 100)
- ✅ Memory-limited (O(N) storage)
- ✅ Gradients available
- ✅ Moderately ill-conditioned problems
- ❌ Very ill-conditioned (use BFGS with preconditioning)

### BFGS - Quasi-Newton Method

**Algorithm**: Build inverse Hessian approximation from gradient differences.

**Convergence**: Superlinear, approaching quadratic near solution.

```cpp
RosenbrockWithGrad func;
VectorN<Real, 2> start = {-1.0, 1.0};

// BFGS quasi-Newton
BFGS bfgs(/*ftol=*/3e-8, /*gtol=*/1e-8, /*maxIter=*/200);
auto result = bfgs.Minimize(func, start);

// Or convenience wrapper
auto result2 = BFGSMinimize(func, start);

std::cout << "BFGS: " << result.iterations << " iterations\n";
// Typically 20-30 iterations for Rosenbrock
```

**BFGS Update Formula**:
```
H_{k+1} = H_k + ρ_k s_k s_k^T - ρ_k H_k y_k y_k^T H_k / (y_k^T H_k y_k) + ρ_k ...
where:
  s_k = x_{k+1} - x_k
  y_k = g_{k+1} - g_k
  ρ_k = 1 / (y_k^T s_k)
```

**When to Use BFGS**:
- ✅ Smooth, well-behaved functions
- ✅ Moderate dimensions (N < 1000)
- ✅ Need fast convergence
- ✅ Gradients available (or cheap to approximate)
- ❌ Very large N (use L-BFGS or CG)
- ❌ Limited memory (O(N²) storage)

---

## Line Search

All multidimensional methods use **line minimization** along search directions.

### Line Function Wrapper

```cpp
// Converts N-dim function to 1D: g(t) = f(p + t*xi)
LineFunction<N> lineFunc(func, point, direction);

// For differentiable functions:
DLineFunction<N> dlineFunc(diffFunc, point, direction);
Real dirDerivative = dlineFunc.derivative(t);  // ∇f · direction
```

### LineMinimizer

```cpp
// Minimize along direction, updating point and direction in-place
VectorN<Real, 2> p = {0.0, 0.0};
VectorN<Real, 2> xi = {1.0, 0.5};  // Search direction

Real fmin = LineMinimizer::Minimize(func, p, xi);
// p is now updated to minimum along line
// xi is scaled to actual displacement
```

---

## Configuration API

### OptimizationConfig - Fluent Builder

```cpp
#include "algorithms/Optimization/OptimizationConfig.h"

// Configure optimization with fluent API
OptimizationConfig<double> config;
config.WithMaxIterations(5000)
      .WithTolerance(1e-8)
      .WithConsoleOutput(50)      // Print every 50 iterations
      .WithTrajectory(10);        // Record every 10 iterations

// Use with optimizer
// SimulatedAnnealing<double> sa(config);
// auto result = sa.Minimize(func, x0);
```

### Termination Criteria

| Criterion | Description | Usage |
|-----------|-------------|-------|
| `MaxIterationsCriterion` | Stop after N iterations | Safety net |
| `FunctionToleranceCriterion` | |f_{k+1} - f_k| < ftol | Value convergence |
| `GradientNormCriterion` | ‖∇f‖ < gtol | Gradient convergence |
| `StagnationCriterion` | No improvement for N iters | Plateau detection |
| `TimeLimitCriterion` | Wall-clock time limit | Time-constrained |
| `TargetValueCriterion` | f(x) ≤ target | Known optimum |

**Composite Criteria**:

```cpp
// Stop on ANY: maxIter OR stagnation
config.WithAnyCriterion()
      .OrCriterion<MaxIterationsCriterion<double>>(1000)
      .OrCriterion<StagnationCriterion<double>>(100);

// Stop on ALL: ftol AND gtol (then maxIter as safety)
config.WithAllCriteria()
      .AndCriterion<FunctionToleranceCriterion<double>>(1e-8)
      .AndCriterion<GradientNormCriterion<double>>(1e-6);
```

**Convenience Presets**:

```cpp
// Standard: maxIter OR stagnation
config.WithStandardTermination(/*maxIter=*/1000, /*stagnation=*/100);

// Precise: (ftol AND gtol) OR maxIter
config.WithPreciseTermination(/*ftol=*/1e-8, /*gtol=*/1e-6, /*maxIter=*/10000);

// Timed: time limit with safety maxIter
config.WithTimedTermination(/*seconds=*/60.0, /*maxIter=*/100000);
```

### Observers

**Console Output**:
```cpp
config.WithConsoleOutput(/*printEvery=*/50, /*verbose=*/true);
```

**Trajectory Recording**:
```cpp
config.WithTrajectory(/*saveEvery=*/10, /*maxSize=*/1000);
// Access via TrajectoryObserver::GetTrajectory()
```

**Custom Callbacks**:
```cpp
config.WithCallback([](const OptimizationState<double>& state) {
    std::cout << "Iter " << state.iteration << ": f = " << state.fval << "\n";
    return true;  // false to stop early
});
```

---

## Algorithm Selection Guide

### Decision Tree

```
Is gradient available?
├─ NO → Is function smooth?
│       ├─ YES → Powell (moderate N) or NelderMead (any N)
│       └─ NO  → NelderMead (robust to noise/discontinuities)
│
└─ YES → How large is N?
         ├─ N < 1000 → BFGS (best general choice)
         └─ N ≥ 1000 → Conjugate Gradient (memory efficient)
```

### Method Comparison

| Scenario | Best Method | Alternative |
|----------|-------------|-------------|
| Black-box function | NelderMead | Powell |
| Smooth, small N (<20) | BFGS | Powell |
| Smooth, large N (>100) | CG | Limited-memory BFGS |
| Noisy function | NelderMead | - |
| High precision needed | BFGS | CG with restart |
| Memory constrained | CG | NelderMead |
| Function evals expensive | BFGS | CG |
| 1D optimization | Brent | Golden Section |

### Performance Guidelines

1. **Start with BFGS** if gradients available — it's the best general-purpose method
2. **Use NelderMead** for truly black-box functions or noisy evaluations
3. **Switch to CG** for large N (>1000) or memory constraints
4. **Always bracket first** in 1D before refinement
5. **Scale variables** to similar magnitudes for better conditioning
6. **Check gradient norm** for convergence, not just function value
7. **Use restarts** if CG stalls (reset to steepest descent)

---

## Mathematical Background

### Optimality Conditions

**First-Order (Necessary)**: At minimum x*, ∇f(x*) = 0

**Second-Order (Sufficient)**: At x* with ∇f(x*) = 0:
- If H(x*) positive definite → x* is local minimum
- If H(x*) negative definite → x* is local maximum
- If H(x*) indefinite → x* is saddle point

### Convergence Theory

**Linear Convergence**: ‖x_{k+1} - x*‖ ≤ c‖x_k - x*‖, c < 1
- Golden section: c = 0.618

**Superlinear Convergence**: ‖x_{k+1} - x*‖ / ‖x_k - x*‖ → 0
- CG, BFGS: Eventually superlinear

**Quadratic Convergence**: ‖x_{k+1} - x*‖ ≤ c‖x_k - x*‖²
- Newton's method (with exact Hessian)

### Condition Number

**Definition**: κ = λ_max / λ_min (ratio of Hessian eigenvalues)

**Effect on Convergence**:
- κ ≈ 1: Fast convergence (nearly spherical contours)
- κ >> 1: Slow convergence (elongated valleys)

**Remedies**:
- Preconditioning
- Variable scaling
- Quasi-Newton methods (adapt to curvature)

---

## Error Handling

```cpp
try {
    auto result = Minimization::BrentMinimize(func, a, b);
    if (!result.converged) {
        std::cerr << "Warning: Did not converge in " 
                  << result.iterations << " iterations\n";
    }
} catch (const OptimizationError& e) {
    std::cerr << "Optimization failed: " << e.what() << "\n";
    std::cerr << "Iterations: " << e.iterations() << "\n";
}

// Multidimensional
try {
    auto result = BFGSMinimize(func, start);
} catch (const MultidimOptimizationError& e) {
    std::cerr << "Multidim optimization failed: " << e.what() << "\n";
}
```

---

## Complete Examples

### Example 1: Simple 1D Minimization

```cpp
#include "algorithms/Optimization.h"
using namespace MML;

int main() {
    // Minimize f(x) = x⁴ - 3x³ + 2
    RealFunctionFromStdFunc f([](Real x) {
        return x*x*x*x - 3*x*x*x + 2;
    });
    
    // Find and refine minimum
    auto result = Minimization::BrentMinimize(f, 0.0, 3.0);
    
    std::cout << "Minimum at x = " << result.xmin << "\n";
    std::cout << "f(x*) = " << result.fmin << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    
    return 0;
}
```

### Example 2: Multidimensional Without Derivatives

```cpp
#include "algorithms/Optimization/OptimizationMultidim.h"
using namespace MML;

// Himmelblau's function: has 4 local minima
class Himmelblau : public IScalarFunction<2> {
public:
    Real operator()(const VectorN<Real, 2>& x) const override {
        Real a = x[0]*x[0] + x[1] - 11;
        Real b = x[0] + x[1]*x[1] - 7;
        return a*a + b*b;
    }
};

int main() {
    Himmelblau func;
    
    // Find different minima from different starting points
    std::vector<VectorN<Real, 2>> starts = {
        {0.0, 0.0}, {-4.0, 0.0}, {4.0, 0.0}, {0.0, 4.0}
    };
    
    for (const auto& start : starts) {
        auto result = NelderMeadMinimize(func, start);
        std::cout << "From (" << start[0] << ", " << start[1] << "): ";
        std::cout << "min at (" << result.xmin[0] << ", " << result.xmin[1] << ")";
        std::cout << ", f = " << result.fmin << "\n";
    }
    
    return 0;
}
```

### Example 3: BFGS with Gradients

```cpp
#include "algorithms/Optimization/OptimizationMultidim.h"
using namespace MML;

// Quadratic bowl: f(x) = x^T A x + b^T x
class QuadraticBowl : public IDifferentiableScalarFunction<3> {
    Matrix<Real> A;
    Vector<Real> b;
    
public:
    QuadraticBowl() : A(3, 3, 0.0), b(3) {
        // Positive definite matrix
        A(0,0) = 4; A(0,1) = 1; A(0,2) = 0;
        A(1,0) = 1; A(1,1) = 3; A(1,2) = 0.5;
        A(2,0) = 0; A(2,1) = 0.5; A(2,2) = 2;
        b[0] = 1; b[1] = 2; b[2] = 1;
    }
    
    Real operator()(const VectorN<Real, 3>& x) const override {
        Real result = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
                result += x[i] * A(i,j) * x[j];
            result += b[i] * x[i];
        }
        return result;
    }
    
    void Gradient(const VectorN<Real, 3>& x, VectorN<Real, 3>& g) const override {
        // ∇f = 2Ax + b
        for (int i = 0; i < 3; ++i) {
            g[i] = b[i];
            for (int j = 0; j < 3; ++j)
                g[i] += 2 * A(i,j) * x[j];
        }
    }
};

int main() {
    QuadraticBowl func;
    VectorN<Real, 3> start = {10.0, 10.0, 10.0};
    
    // BFGS should find minimum quickly for quadratic
    BFGS optimizer(1e-10, 1e-10, 100);
    auto result = optimizer.Minimize(func, start);
    
    std::cout << "BFGS found minimum in " << result.iterations << " iterations\n";
    std::cout << "x* = (" << result.xmin[0] << ", " 
              << result.xmin[1] << ", " << result.xmin[2] << ")\n";
    std::cout << "f(x*) = " << result.fmin << "\n";
    
    return 0;
}
```

---

## References

1. **Numerical Recipes** (Press et al.) - Chapters 10.1-10.8
2. **Nocedal & Wright** - "Numerical Optimization" - Comprehensive reference
3. **Fletcher** - "Practical Methods of Optimization"
4. **Gill, Murray, Wright** - "Practical Optimization"

---

## See Also

- [Root_finding.md](Root_finding.md) - Find zeros of functions
- [Differential_equations_solvers.md](Differential_equations_solvers.md) - Solve ODEs  
- [Statistics.md](Statistics.md) - Least squares fitting (optimization application)
- [MonteCarlo_integration.md](MonteCarlo_integration.md) - Stochastic optimization connection
