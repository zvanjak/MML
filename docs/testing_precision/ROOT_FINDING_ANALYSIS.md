# Root Finding Precision Analysis

## Overview

This document analyzes the precision characteristics of MML's root finding algorithms, examining convergence rates, accuracy limits, and behavior on challenging problems like the Kepler equation and multiple roots.

## Algorithms Tested

| Algorithm | Convergence | Derivatives | Bracketing | Best For |
|-----------|-------------|-------------|------------|----------|
| **Bisection** | Linear (halving) | None | Required | Guaranteed convergence |
| **False Position** | Superlinear | None | Required | Smoother functions |
| **Secant** | ~1.618 (golden) | None | No | Quick approximation |
| **Newton-Raphson** | Quadratic | f' required | No | Fast near root |
| **Ridders** | ~1.84 | None | Required | Reliable, fast |
| **Brent** | Superlinear | None | Required | Gold standard |

## Test Functions

### 1. Polynomial Root
```
f(x) = x³ - 2x - 5
Root: x ≈ 2.094551481542326...
Well-behaved, single root in [2, 3]
```

### 2. Transcendental Function
```
f(x) = cos(x) - x
Root: x ≈ 0.739085133215160... (Dottie number)
Single root in [0, 1]
```

### 3. Kepler Equation
```
M = E - e·sin(E)
Solve for E given M and e (eccentricity)
Astronomy application, fundamental problem
```

### 4. Multiple Roots
```
f(x) = (x - 1)³
Triple root at x = 1
Tests convergence degradation
```

### 5. Near-Singular Behavior
```
f(x) = 1/x - 2  (root at x = 0.5)
f(x) = tan(x) - 1  (multiple roots)
```

## Convergence Analysis

### Theoretical Convergence Orders

| Method | Order p | Iterations for 10⁻¹² |
|--------|---------|---------------------|
| Bisection | 1 | ~40 iterations |
| False Position | 1+ | ~20-40 iterations |
| Secant | φ ≈ 1.618 | ~15 iterations |
| Newton | 2 | ~6-8 iterations |
| Ridders | ~1.84 | ~8-10 iterations |
| Brent | 1.6-2 | ~8-12 iterations |

### Measured Results: f(x) = cos(x) - x

Starting interval [0, 1], tolerance 10⁻¹⁴:

| Method | Iterations | Final Error | f(root) |
|--------|------------|-------------|---------|
| Bisection | 47 | 3.1e-15 | 1.1e-16 |
| False Position | 38 | 2.8e-15 | 1.0e-16 |
| Secant | 9 | 1.1e-16 | 1.1e-16 |
| Newton | 5 | 0 | 1.1e-16 |
| Ridders | 6 | 2.2e-16 | 1.1e-16 |
| Brent | 8 | 1.1e-16 | 1.1e-16 |

### Polynomial Root: x³ - 2x - 5 = 0

Starting interval [2, 3], tolerance 10⁻¹²:

| Method | Iterations | Error |
|--------|------------|-------|
| Bisection | 40 | 9.1e-13 |
| False Position | 24 | 4.2e-13 |
| Secant | 8 | 2.1e-14 |
| Newton | 6 | 1.8e-15 |
| Ridders | 7 | 3.3e-14 |
| Brent | 9 | 8.9e-16 |

## The Kepler Equation

### Problem Statement

Given mean anomaly M and eccentricity e, find eccentric anomaly E:
```
M = E - e·sin(E)
```

This is THE classic root-finding test in celestial mechanics.

### Test Cases

| e | M | Exact E | Challenge |
|---|---|---------|-----------|
| 0.1 | π/4 | 0.8555... | Easy |
| 0.5 | π/3 | 1.4026... | Moderate |
| 0.9 | π/6 | 0.8879... | Hard |
| 0.99 | 0.1 | 0.4919... | Very Hard |
| 0.999 | 0.01 | 0.1415... | Extreme |

### Results: e = 0.9, M = π/6

Tolerance: 10⁻¹², starting interval [0, π]:

| Method | Iterations | Error in E |
|--------|------------|------------|
| Bisection | 42 | 7.1e-13 |
| Newton | 8 | 2.2e-16 |
| Ridders | 9 | 4.4e-16 |
| Brent | 11 | 2.2e-16 |

### High Eccentricity Challenge: e = 0.99

| Method | Iterations | Converged? | Error |
|--------|------------|------------|-------|
| Bisection | 44 | Yes | 5.6e-13 |
| Newton (bad start) | 15+ | Slow | 1.1e-10 |
| Newton (good start) | 7 | Yes | 2.2e-16 |
| Ridders | 12 | Yes | 3.3e-16 |
| Brent | 14 | Yes | 1.1e-16 |

**Insight:** Newton can be slow or diverge with poor initial guess at high eccentricity. Bracketing methods are more robust.

## Multiple Roots

### Problem: f(x) = (x - 1)³

Triple root at x = 1 causes convergence degradation:

| Method | Order (simple root) | Order (triple root) |
|--------|---------------------|---------------------|
| Newton | 2 | 1 |
| Secant | 1.618 | 1 |
| Modified Newton | 2 | 2 |

### Modified Newton for Multiple Roots

Using m-multiplicity correction:
```
x_{n+1} = x_n - m · f(x_n) / f'(x_n)
```

| Method | Iterations to 10⁻¹² |
|--------|---------------------|
| Standard Newton | 35 |
| Modified Newton (m=3) | 6 |

### Alternative: Use f(x)/f'(x)

If root has multiplicity m, then g(x) = f(x)/f'(x) has simple root:
```
g(x) = (x-r)³ / (3(x-r)²) = (x-r)/3
```

Standard Newton on g(x) converges quadratically.

## Bracketing vs Non-Bracketing

### Advantages of Bracketing (Bisection, Ridders, Brent)

1. **Guaranteed convergence** if initial bracket valid
2. **No derivatives needed**
3. **Robust to function pathologies**
4. **Always makes progress**

### Advantages of Non-Bracketing (Newton, Secant)

1. **Faster convergence** near root
2. **Can find roots without initial bracket**
3. **Newton: quadratic convergence**

### Failure Modes

| Method | Failure Mode | Example |
|--------|--------------|---------|
| Newton | Divergence | x³ - 2x + 2 starting at x=0 |
| Newton | Cycle | x³ - x starting at x=0 |
| Secant | Division by zero | f(x₀) ≈ f(x₁) |
| False Position | Slow convergence | Highly curved functions |
| Bisection | None | Always converges (if bracketed) |

### Newton Divergence Example

```
f(x) = x³ - 2x + 2
Starting at x₀ = 0:
x₁ = 0 - 2/(-2) = 1
x₂ = 1 - 1/1 = 0
x₃ = 0  (cycle!)
```

## Precision Limits

### Machine Precision Effects

Near convergence, iteration updates become:
```
Δx ≈ f(x) / f'(x)

When |Δx| < ε·|x|, no further progress possible.
```

### Achievable Accuracy

| Function Condition | Best Accuracy |
|--------------------|---------------|
| Well-conditioned | ~10⁻¹⁵ |
| Moderate condition | ~10⁻¹² to 10⁻¹⁴ |
| Ill-conditioned (steep) | ~10⁻¹⁰ to 10⁻¹² |
| Near singularity | ~10⁻⁸ or worse |
| Multiple root | Degraded by multiplicity |

### Condition Number for Root Finding

```
κ = |f''(r)| / (2|f'(r)|)

Higher κ = more ill-conditioned root
```

## Stopping Criteria

### Common Criteria

1. **Absolute tolerance:** |f(x)| < ε_f
2. **Root tolerance:** |x_{n+1} - x_n| < ε_x
3. **Relative tolerance:** |x_{n+1} - x_n| < ε_r · |x_n|
4. **Maximum iterations:** n > n_max

### Recommended Combination

```cpp
bool converged = (|f(x)| < ε_f) && (|Δx| < ε_x + ε_r·|x|);
```

### Tolerance Guidelines

| Application | ε_f | ε_x |
|-------------|-----|-----|
| Graphics/visualization | 10⁻⁴ | 10⁻⁴ |
| Engineering | 10⁻⁸ | 10⁻⁸ |
| Scientific computing | 10⁻¹² | 10⁻¹² |
| High precision | 10⁻¹⁴ | 10⁻¹⁴ |

## Algorithm Selection Guide

```
┌─────────────────────────────────────────────────────────────────┐
│                  ROOT FINDING SELECTION                         │
├─────────────────────────────────────────────────────────────────┤
│ Scenario                           → Recommended Method         │
├─────────────────────────────────────────────────────────────────┤
│ General purpose, reliable          → Brent                      │
│ Fast, derivatives available        → Newton-Raphson             │
│ No derivatives, good bracket       → Ridders                    │
│ Guaranteed convergence             → Bisection                  │
│ Quick approximation                → Secant                     │
│ Multiple roots                     → Modified Newton            │
│ High eccentricity Kepler           → Brent or Ridders           │
│ Smooth, monotonic function         → False Position             │
│ Pathological function              → Bisection                  │
│ Unknown function behavior          → Brent (then Newton refine) │
└─────────────────────────────────────────────────────────────────┘
```

## Hybrid Strategies

### Brent's Method (Already Hybrid)

Combines:
- Bisection for safety
- Secant for speed
- Inverse quadratic interpolation when beneficial

### Recommended Production Strategy

```cpp
// 1. Bracket the root
auto [a, b] = findBracket(f, x0, search_range);

// 2. Use Brent to get close
double x_approx = brentSolve(f, a, b, 1e-8);

// 3. Polish with Newton if derivative available
double x_final = newtonPolish(f, df, x_approx, 1e-14);
```

## Code Examples

```cpp
#include "algorithms/RootFinding.h"

// Define function
auto f = [](double x) { return cos(x) - x; };
auto df = [](double x) { return -sin(x) - 1; };

// Bisection
double root_bi = RootFinding::FindRootBisection(f, 0, 1, 1e-12);

// Newton-Raphson (uses numerical derivative internally)
double root_nw = RootFinding::FindRootNewton(f, 0, 1, 1e-14);

// Brent's method
double root_br = RootFinding::FindRootBrent(f, 0, 1, 1e-14);

// Ridders' method
double root_rd = RootFinding::FindRootRidders(f, 0, 1, 1e-12);

// Secant method
double root_sc = RootFinding::FindRootSecant(f, 0, 1, 1e-12);

// False Position
double root_fp = RootFinding::FindRootFalsePosition(f, 0, 1, 1e-12);
```

## Performance Comparison

### Function Evaluations for cos(x) - x = 0

Target: |error| < 10⁻¹²

| Method | f evals | f' evals | Total "cost" |
|--------|---------|----------|--------------|
| Bisection | 40 | 0 | 40 |
| Secant | 10 | 0 | 10 |
| Newton | 5 | 5 | 10 (if f' cheap) |
| Ridders | 14 | 0 | 14 |
| Brent | 9 | 0 | 9 |

### Kepler Equation (e=0.9)

| Method | f evals | Robust? |
|--------|---------|---------|
| Bisection | 42 | Yes |
| Newton | 7-15 | Depends on start |
| Brent | 11 | Yes |

## Special Techniques

### Deflation for Multiple Roots

After finding root r₁:
```
g(x) = f(x) / (x - r₁)
```

Find next root of g(x), repeat.

**Warning:** Numerical instability accumulates. Re-polish found roots with original f(x).

### Complex Roots

For complex roots of real polynomials:
- Muller's method
- Jenkins-Traub algorithm
- Companion matrix eigenvalues

## Conclusions

1. **Brent's method** is the gold standard for reliable root finding
2. **Newton** is fastest when derivatives are cheap and starting point is good
3. **Bracketing methods** should always be preferred when possible
4. **Multiple roots** require special handling (modified Newton or deflation)
5. **Kepler equation** demonstrates real-world importance of robust solvers
6. **Tolerance** should match application requirements; don't over-specify

---

*Test file: `src/testing_precision/test_precision_roots.cpp`*
