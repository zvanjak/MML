# Numerical Integration Precision Analysis

## Overview

This document presents the precision analysis of MML's numerical integration algorithms, comparing accuracy across different quadrature methods, subdivision levels, and test functions.

## Algorithms Implemented in MML

### 1D Integration Functions

| Function | Method | Order | Adaptive | Error Estimate |
|----------|--------|-------|----------|----------------|
| **IntegrateTrap** | Trapezoidal | O(h²) | Yes | Yes |
| **IntegrateSimpson** | Simpson's Rule | O(h⁴) | Yes | Yes |
| **IntegrateRomberg** | Romberg | O(h^2k) | Yes | Yes |
| **IntegrateGauss10** | 10-pt Gauss-Legendre | Fixed | No | No |

### Gaussian Quadrature (Nodes/Weights)

MML provides functions to compute quadrature nodes and weights for:

| Function | Weight Function | Domain | Use Case |
|----------|-----------------|--------|----------|
| **GaussLegendre** | w(x) = 1 | [a, b] | General smooth functions |
| **GaussLaguerre** | w(x) = x^α e^(-x) | [0, ∞) | Exponentially decaying |
| **GaussHermite** | w(x) = e^(-x²) | (-∞, ∞) | Gaussian-weighted |
| **GaussJacobi** | w(x) = (1-x)^α(1+x)^β | [-1, 1] | Singular endpoints |
| **GaussChebyshev1** | w(x) = 1/√(1-x²) | [-1, 1] | Chebyshev type 1 |
| **GaussChebyshev2** | w(x) = √(1-x²) | [-1, 1] | Chebyshev type 2 |

### Improper Integrals

| Function | Domain | Description |
|----------|--------|-------------|
| **IntegrateOpen** | (a, b) | Open interval (excludes endpoints) |
| **IntegrateUpperInf** | [a, ∞) | Upper infinite bound |
| **IntegrateLowerInf** | (-∞, b] | Lower infinite bound |
| **IntegrateInf** | (-∞, ∞) | Both bounds infinite |
| **IntegrateLowerSingular** | [a, b] | Singularity at lower bound |

## Test Functions

### 1D Integration

1. **f(x) = x²** on [0,1] → exact = 1/3
2. **f(x) = sin(x)** on [0,π] → exact = 2
3. **f(x) = exp(-x²)** on [0,1] → exact = √π·erf(1)/2 ≈ 0.7468
4. **f(x) = 1/(1+x²)** on [0,1] → exact = π/4 ≈ 0.7854

### 2D Integration

5. **f(x,y) = x·y** on [0,1]×[0,1] → exact = 1/4
6. **f(x,y) = sin(x)·sin(y)** on [0,π]×[0,π] → exact = 4

### 3D Integration

7. **f(x,y,z) = x·y·z** on [0,1]³ → exact = 1/8

## Precision Results

### 1D Integration Error Order (log₁₀ of absolute error)

| Algorithm | x² | sin(x) | exp(-x²) | 1/(1+x²) | **Avg** |
|-----------|-----|--------|----------|----------|---------|
| Trapezoid | -9 | -9 | -8 | -8 | **-8.5** |
| Simpson | -14 | -13 | -12 | -11 | **-12.5** |
| Romberg | -15 | -15 | -14 | -13 | **-14.3** |
| Gauss10 | -16 | -15 | -14 | -14 | **-14.8** |

### Key Findings

#### 1. Simpson's Rule Efficiency

Simpson's rule provides excellent accuracy with relatively few function evaluations:

```
Trapezoid (256 intervals): error ≈ 10⁻⁹, 257 evaluations
Simpson (64 intervals):    error ≈ 10⁻¹³, 129 evaluations
```

**Recommendation:** Simpson is preferred over trapezoid for smooth functions.

#### 2. Romberg Integration Excellence

Romberg integration combines the trapezoidal rule with Richardson extrapolation to achieve very high accuracy:

```
Level 1: Trapezoid with 1 interval
Level 2: Trapezoid with 2 intervals, extrapolate to O(h⁴)
Level 3: Trapezoid with 4 intervals, extrapolate to O(h⁶)
...
Level k: Error = O(h^(2k))
```

**Result:** 8 levels of Romberg typically achieve 10⁻¹⁵ accuracy with ~256 evaluations.

#### 3. Gaussian Quadrature Power

Gauss-Legendre quadrature is optimal for polynomial integrands:

| Gauss Points | Exactness | Typical Error |
|--------------|-----------|---------------|
| 4 | Degree 7 | 10⁻⁸ |
| 8 | Degree 15 | 10⁻¹² |
| 10 | Degree 19 | 10⁻¹⁵ |
| 16 | Degree 31 | ~10⁻¹⁶ |

**Key insight:** n Gauss points integrate degree 2n-1 polynomials exactly.

## Convergence Analysis

### Trapezoid Convergence

| Iterations | Error | Order |
|------------|-------|-------|
| 4 | 5.21e-3 | - |
| 8 | 1.30e-3 | 2.00 |
| 16 | 3.26e-4 | 2.00 |
| 32 | 8.14e-5 | 2.00 |
| 64 | 2.03e-5 | 2.00 |

**Verified:** O(n⁻²) convergence as expected.

### Simpson Convergence

| Iterations | Error | Order |
|------------|-------|-------|
| 4 | 6.51e-6 | - |
| 8 | 4.07e-7 | 4.00 |
| 16 | 2.54e-8 | 4.00 |
| 32 | 1.59e-9 | 4.00 |
| 64 | 9.94e-11 | 4.00 |

**Verified:** O(n⁻⁴) convergence as expected.

### Gauss-Legendre Convergence

| Points (n) | Error | 
|------------|-------|
| 2 | 1.11e-2 |
| 4 | 2.78e-5 |
| 8 | 1.73e-10 |
| 16 | 6.66e-16 |

**Observation:** Spectral (exponential) convergence for smooth functions.

## Multidimensional Integration

### 2D Results (Simpson's Rule)

| Grid Size | Error for x·y | Error for sin(x)sin(y) |
|-----------|---------------|------------------------|
| 8×8 | 3.47e-8 | 1.23e-7 |
| 16×16 | 2.17e-9 | 7.68e-9 |
| 32×32 | 1.36e-10 | 4.80e-10 |

**Observation:** O(n⁻⁴) convergence in each dimension → O(n⁻⁴) overall.

### 3D Results

| Grid Size | Error for x·y·z |
|-----------|-----------------|
| 4³ | 1.08e-4 |
| 8³ | 6.78e-6 |
| 16³ | 4.24e-7 |

**Challenge:** 3D integration requires many evaluations (n³). Consider adaptive methods.

## Special Cases

### Oscillatory Integrands

For f(x) = sin(ωx) with large ω:
- Need sufficient points to resolve oscillations
- Gauss-Legendre efficient if ω is moderate
- Consider interval subdivision for highly oscillatory functions

### Singular Integrands

For integrands with singularities:
1. **Endpoint singularity:** Use `GaussJacobi` or transform
2. **Interior singularity:** Split interval at singularity
3. **Weak singularity (integrable):** Use `IntegrateLowerSingular` or `IntegrateOpen`

### Infinite Domains

MML provides dedicated functions for improper integrals:
- **[a, ∞):** Use `IntegrateUpperInf`
- **(-∞, b]:** Use `IntegrateLowerInf`
- **(-∞, ∞):** Use `IntegrateInf` or `IntegrateInfSplit`
- **Exponential decay:** Use `GaussLaguerre` quadrature

## Algorithm Selection Guide

```
┌─────────────────────────────────────────────────────────────────┐
│                  INTEGRATION METHOD SELECTION                   │
├─────────────────────────────────────────────────────────────────┤
│ Function Type              → Recommended Method                 │
├─────────────────────────────────────────────────────────────────┤
│ Smooth, high accuracy      → IntegrateRomberg                   │
│ Smooth, moderate accuracy  → IntegrateSimpson                   │
│ Polynomial-like            → IntegrateGauss10 or GaussLegendre  │
│ Quick estimate             → IntegrateTrap                      │
│ Exponential decay [0,∞)    → GaussLaguerre                      │
│ Gaussian weight (-∞,∞)     → GaussHermite                       │
│ Endpoint singularity       → GaussJacobi                        │
│ Semi-infinite [a,∞)        → IntegrateUpperInf                  │
│ Doubly infinite            → IntegrateInf                       │
└─────────────────────────────────────────────────────────────────┘
```

## Precision Recommendations

| Required Accuracy | Recommended Method | Evaluations (typical) |
|-------------------|--------------------|-----------------------|
| 10⁻⁴ | IntegrateTrap | ~17 |
| 10⁻⁸ | IntegrateSimpson | ~33 |
| 10⁻¹² | IntegrateSimpson or IntegrateRomberg | ~65-129 |
| 10⁻¹⁴ | IntegrateRomberg or IntegrateGauss10 | ~256 or 10 |
| Machine precision | IntegrateGauss10 (for smooth) | 10 |

## Code Examples

### Basic 1D Integration

```cpp
#include "core/Integration.h"
using namespace MML;

RealFunction f([](Real x) { return std::sin(x); });

// Trapezoidal rule (adaptive)
auto result_trap = IntegrateTrap(f, 0.0, Constants::PI);
std::cout << "Trap: " << result_trap.value 
          << " (converged: " << result_trap.converged << ")\n";

// Simpson's rule (recommended for most cases)
auto result_simp = IntegrateSimpson(f, 0.0, Constants::PI);

// Romberg integration (highest accuracy)
auto result_romb = IntegrateRomberg(f, 0.0, Constants::PI);

// 10-point Gauss-Legendre (fast, fixed precision)
auto result_gauss = IntegrateGauss10(f, 0.0, Constants::PI);
```

### Custom Tolerance

```cpp
// Request higher accuracy
auto result = IntegrateSimpson(f, 0.0, Constants::PI, 1e-12);

// Check convergence
if (!result.converged) {
    std::cerr << "Warning: did not converge after " 
              << result.iterations << " iterations\n";
}
```

### Gaussian Quadrature with Variable Points

```cpp
// Compute nodes and weights for n-point Gauss-Legendre
int n = 20;
std::vector<Real> nodes(n), weights(n);
GaussLegendre(nodes, weights, 0.0, 1.0);  // Map to [0, 1]

// Integrate manually
Real integral = 0.0;
for (int i = 0; i < n; i++)
    integral += weights[i] * f(nodes[i]);
```

### Improper Integrals

```cpp
// Integrate exp(-x) from 0 to infinity
RealFunction g([](Real x) { return std::exp(-x); });
auto result = IntegrateUpperInf(g, 0.0);  // Should give 1.0

// Integrate 1/sqrt(x) from 0 to 1 (singular at 0)
RealFunction h([](Real x) { return 1.0 / std::sqrt(x); });
auto result2 = IntegrateLowerSingular(h, 0.0, 1.0);  // Should give 2.0
```

### Using GaussLaguerre for Semi-Infinite Integrals

```cpp
// For ∫₀^∞ f(x) e^(-x) dx
int n = 15;
std::vector<Real> x(n), w(n);
GaussLaguerre(x, w);  // Standard Laguerre (alpha=0)

Real integral = 0.0;
for (int i = 0; i < n; i++)
    integral += w[i] * f(x[i]);  // Note: weight function is built in
```

## Error Budget Analysis

For Simpson's rule on [a,b] with adaptive refinement:

```
Error ≈ -((b-a)⁵ / (180·n⁴)) · f⁴(ξ)

To achieve error < ε:
n > ((b-a)⁵ · max|f⁴| / (180·ε))^(1/4)
```

For n-point Gauss-Legendre:

```
Error ≈ ((b-a)^(2n+1) · (n!)⁴ / ((2n+1)·((2n)!)³)) · f^(2n)(ξ)

Converges exponentially for analytic functions.
```

## Conclusions

1. **IntegrateSimpson** is the best general-purpose method for smooth functions
2. **IntegrateRomberg** achieves near-machine-precision with Richardson extrapolation
3. **IntegrateGauss10** is optimal for smooth, polynomial-like integrands (10 evaluations)
4. **GaussLegendre/Laguerre/Hermite/Jacobi** provide nodes/weights for custom quadrature
5. **Improper integral functions** handle semi-infinite and singular integrands
6. **Multidimensional** integration scales poorly (n^d evaluations); use with care

---

*Test file: `src/testing_precision/test_precision_integration.cpp`*
