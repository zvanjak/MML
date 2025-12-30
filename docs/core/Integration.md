# Numerical Integration - 1D Quadrature Methods

This document covers numerical integration (quadrature) of one-dimensional functions ∫ᵃᵇ f(x) dx.

## Table of Contents
- [Overview](#overview)
- [Method Selection Guide](#method-selection-guide)
- [API Reference](#api-reference)
- [Trapezoidal Rule](#trapezoidal-rule)
- [Simpson's Rule](#simpsons-rule)
- [Romberg Integration](#romberg-integration)
- [Gaussian Quadrature](#gaussian-quadrature)
- [Integration Result](#integration-result)
- [Examples](#examples)
- [Convergence and Accuracy](#convergence-and-accuracy)
- [Common Pitfalls](#common-pitfalls)
- [See Also](#see-also)

---

## Overview

MinimalMathLibrary provides **adaptive quadrature methods** for integrating one-dimensional real functions:

| Method | Order | Best For | Function Evals | Adaptive |
|--------|-------|----------|----------------|----------|
| **Trapezoidal** | O(h²) | Non-smooth, data points | 2ᵏ + 1 | Yes |
| **Simpson** | O(h⁴) | Smooth functions | 2ᵏ + 1 | Yes |
| **Romberg** | O(h^(2K)) | Analytic functions | Few | Yes |
| **Gauss10** | O(h²⁰) | High precision, fixed cost | 10 | No |

**Key Features:**
- **Adaptive refinement:** Trapezoidal, Simpson, and Romberg automatically subdivide until convergence
- **Error estimation:** Convergence-based error estimates
- **IntegrationResult:** Structured result with value, error estimate, iterations, convergence status
- **Default method:** `Integrate()` uses Trapezoidal (robust)

**Recommendation:** Start with **Simpson** for general use. Use **Trapezoidal** for non-smooth integrands or interpolated data. Use **Romberg** or **Gauss10** for smooth, analytic functions.

All integration methods are defined in [`mml/core/Integration/Integration1D.h`](../../mml/core/Integration/Integration1D.h).

---

## Method Selection Guide

### Decision Tree

```
What is your integrand f(x)?

Is f(x) non-smooth (discontinuous derivatives, piecewise, data points)?
├─ YES → Use Trapezoidal
│  └─ Robust to roughness, interpolated data
└─ NO (smooth, continuous derivatives)
   ├─ Is f(x) very smooth (analytic, infinitely differentiable)?
   │  ├─ YES → Use Romberg or Gauss10
   │  │  └─ Romberg: adaptive high-order
   │  │  └─ Gauss10: fixed cost, excellent accuracy
   │  └─ NO → Use Simpson
   │     └─ Best balance for smooth functions
   └─ High dimensional (n > 1)?
      └─ See Multidim_integration.md or Monte Carlo
```

### Function Smoothness Requirements

**Trapezoidal (TRAP):**
- **Requires:** f(x) continuous
- **Optimal for:** Piecewise linear, interpolated data, discontinuous derivatives
- **Convergence:** O(h²) per refinement
- **Robustness:** ✅ Very robust (won't fail on corners/kinks)

**Simpson (SIMPSON):**
- **Requires:** f'''(x) continuous (continuous 3rd derivative)
- **Optimal for:** General smooth functions
- **Convergence:** O(h⁴) per refinement
- **Robustness:** ✅ Good for most cases
- **Recommended:** **Default choice for smooth integrands**

**Romberg (ROMBERG):**
- **Requires:** f(x) analytic (infinitely differentiable)
- **Optimal for:** Very smooth functions on finite intervals
- **Convergence:** O(h^(2K)) - very fast!
- **Robustness:** ⚠️ Sensitive to singularities, endpoint issues
- **Use when:** Maximum accuracy from minimum evaluations

**Gauss10 (GAUSS10):**
- **Requires:** f(x) well-behaved (no singularities)
- **Optimal for:** Fixed budget (always 10 evaluations)
- **Convergence:** Exact for polynomials up to degree 19
- **Robustness:** ⚠️ No adaptive refinement, no error estimate
- **Use when:** Function is smooth and 10 evals sufficient

### Typical Use Cases

| Use Case | Method | Reason |
|----------|--------|--------|
| **General integration** | Simpson | Best balance of speed and accuracy |
| **Measured data** | Trapezoidal | Robust to data noise |
| **Physics simulations** | Simpson | Smooth continuous functions |
| **High precision** | Romberg | Analytic integrand, need accuracy |
| **Fixed cost** | Gauss10 | Budget limited to 10 evaluations |
| **Oscillatory** | Adaptive + small ε | May need many subdivisions |
| **Singularities** | Specialized methods | See singular integration techniques |

---

## API Reference

### Enum: Integration Methods

```cpp
enum IntegrationMethod {
    TRAP,      // Trapezoidal rule
    SIMPSON,   // Simpson's rule
    ROMBERG,   // Romberg extrapolation
    GAUSS10    // 10-point Gauss-Legendre quadrature
};
```

### Trapezoidal Integration

```cpp
// Modern interface: returns IntegrationResult
IntegrationResult IntegrateTrap(
    const IRealFunction& func,
    Real a, Real b,
    Real eps = Defaults::TrapezoidIntegrationEPS  // Default: 1e-6
);

// Legacy interface: returns value only
Real IntegrateTrap(
    const IRealFunction& func,
    Real a, Real b,
    int* doneSteps = nullptr,      // Output: iterations performed
    Real* achievedPrec = nullptr,  // Output: achieved error estimate
    Real eps = Defaults::TrapezoidIntegrationEPS
);
```

### Simpson Integration

```cpp
// Modern interface: returns IntegrationResult
IntegrationResult IntegrateSimpson(
    const IRealFunction& func,
    Real a, Real b,
    Real eps = Defaults::SimpsonIntegrationEPS  // Default: 1e-6
);

// Legacy interface: returns value only
Real IntegrateSimpson(
    const IRealFunction& func,
    Real a, Real b,
    int* doneSteps = nullptr,
    Real* achievedPrec = nullptr,
    Real eps = Defaults::SimpsonIntegrationEPS
);
```

### Romberg Integration

```cpp
// Returns IntegrationResult with high-order Richardson extrapolation
IntegrationResult IntegrateRomberg(
    const IRealFunction& func,
    Real a, Real b,
    Real eps = Defaults::RombergIntegrationEPS  // Default: 1e-10
);
```

### Gaussian Quadrature

```cpp
// Fixed 10-point Gauss-Legendre quadrature
// Always exactly 10 function evaluations
Real IntegrateGauss10(
    const IRealFunction& func,
    Real a, Real b
);
```

### Default Integration Function

```cpp
// Function pointer: can be reassigned to any method
// Default: points to IntegrateTrap
IntegrationResult (*Integrate)(
    const IRealFunction& func,
    Real a, Real b,
    Real eps
);

// Usage:
auto result = Integrate(func, 0.0, 1.0, 1e-8);
```

---

## Trapezoidal Rule

**Formula:** Approximate integral as sum of trapezoids:

```
∫ᵃᵇ f(x) dx ≈ h/2 · [f(a) + 2·∑f(xᵢ) + f(b)]
```

where h = (b-a)/n and xᵢ = a + i·h.

### Algorithm

**Extended Trapezoidal Rule** with successive refinement:

1. **Step 1:** Compute initial estimate with 2 points:
   ```
   I₁ = (b-a)/2 · [f(a) + f(b)]
   ```

2. **Step k:** Add 2^(k-2) new interior points (between existing points):
   ```
   Iₖ = Iₖ₋₁/2 + (b-a)/2^(k-1) · ∑ f(xₙₑw)
   ```

3. **Convergence test:** Stop when |Iₖ - Iₖ₋₁| < ε·|Iₖ₋₁|

### When to Use

✅ **Good for:**
- **Non-smooth functions:** Corners, kinks, discontinuous derivatives
- **Piecewise functions:** Functions defined by cases
- **Interpolated data:** Linear interpolation between points
- **Conservative estimates:** Guaranteed convergence for continuous f(x)
- **Measured data with noise:** Robust to small perturbations

❌ **Not recommended:**
- Smooth analytic functions (Simpson/Romberg faster)
- High-precision requirements (slow convergence)

### Example

```cpp
// Non-smooth function: |x - 0.5|
RealFunction abs_func([](Real x) {
    return abs(x - 0.5);  // Has corner at x = 0.5
});

// Trapezoidal handles corner gracefully
auto result = IntegrateTrap(abs_func, 0.0, 1.0, 1e-6);
std::cout << "Integral: " << result.value 
          << " (error: " << result.error_estimate
          << ", iters: " << result.iterations << ")" << std::endl;

// Expected: 0.25 (two triangles: 0.5 * 0.5 * 0.5 * 2)
```

### Convergence Rate

**Error:** O(h²) where h = (b-a)/n

For smooth functions f(x):
```
Error ≈ -(b-a)³/(12n²) · f''(ξ)   for some ξ ∈ [a,b]
```

**Halving h → Error reduces by factor of 4**

---

## Simpson's Rule

**Formula:** Approximate integral using parabolic segments:

```
∫ᵃᵇ f(x) dx ≈ h/3 · [f(a) + 4·f(m) + f(b)]
```

for interval [a,b] with midpoint m = (a+b)/2.

### Algorithm

**Composite Simpson's Rule** via Richardson extrapolation on Trapezoidal:

```
Sₖ = (4·Tₖ - Tₖ₋₁) / 3
```

where Tₖ is the k-th trapezoidal estimate.

**This cancels the O(h²) error term**, leaving O(h⁴) error!

### When to Use

✅ **Good for:**
- **General-purpose integration:** Most smooth functions
- **Physics/engineering:** Continuous smooth integrands
- **Default choice** when unsure (good speed-accuracy balance)
- **Functions with continuous 3rd derivative**

❌ **Not recommended:**
- Non-smooth functions (use Trapezoidal)
- Discontinuous derivatives (Trapezoidal more robust)

### Example

```cpp
RealFunction smooth([](Real x) {
    return sin(x) * (1.0 + 0.5 * x * x);
});

// Simpson's rule - excellent for smooth functions
auto result = IntegrateSimpson(smooth, 0.0, 10.0, 1e-8);

std::cout << "Simpson result: " << result.value << std::endl;
std::cout << "Error estimate: " << result.error_estimate << std::endl;
std::cout << "Iterations: " << result.iterations << std::endl;
std::cout << "Converged: " << (result.converged ? "YES" : "NO") << std::endl;
```

### Convergence Rate

**Error:** O(h⁴) where h = (b-a)/n

For functions with continuous 4th derivative:
```
Error ≈ -(b-a)⁵/(2880n⁴) · f⁽⁴⁾(ξ)   for some ξ ∈ [a,b]
```

**Halving h → Error reduces by factor of 16** (much faster than Trapezoidal!)

---

## Romberg Integration

**Algorithm:** High-order Richardson extrapolation on Trapezoidal rule.

**Idea:** Combine multiple Trapezoidal estimates with different step sizes to cancel error terms systematically.

### Richardson Extrapolation

Generate Romberg tableau:

```
T₀₀
T₁₀  R₁₁
T₂₀  R₂₁  R₂₂
T₃₀  R₃₁  R₃₂  R₃₃
...
```

Where:
- **Tₖ₀:** k-th Trapezoidal estimate (step h = (b-a)/2^k)
- **Rₖⱼ:** j-th extrapolation at step k

**Extrapolation formula:**
```
Rₖⱼ = (4ʲ·Rₖ,ⱼ₋₁ - Rₖ₋₁,ⱼ₋₁) / (4ʲ - 1)
```

**Each column increases order by 2:**
- Column 1 (j=1): Simpson's rule (O(h⁴))
- Column 2 (j=2): O(h⁶)
- Column 3 (j=3): O(h⁸)
- ...

### When to Use

✅ **Good for:**
- **Analytic functions:** Smooth, infinitely differentiable
- **High precision needed:** Want many correct digits
- **Clean intervals:** No singularities at endpoints or interior

❌ **Not recommended:**
- Non-smooth functions (extrapolation amplifies errors!)
- Functions with singularities (won't converge)
- Endpoint issues (e.g., f(a) or f(b) undefined)

### Example

```cpp
// Very smooth function: e^(-x²)
RealFunction gaussian([](Real x) {
    return exp(-x * x);
});

// Romberg achieves high accuracy with few evaluations
auto result = IntegrateRomberg(gaussian, 0.0, 3.0, 1e-12);

std::cout << "Romberg integral: " << result.value << std::endl;
std::cout << "Error estimate: " << result.error_estimate << std::endl;
std::cout << "Iterations: " << result.iterations << std::endl;

// Expected: ≈ 0.886226925452758 (√π/2 · erf(3) ≈ √π/2)
// Romberg converges in ~5-8 iterations vs ~15-20 for Simpson!
```

### Convergence Rate

**Exponential for analytic f(x):**

Error decreases by factor of ~16-256 per iteration (depending on smoothness).

For very smooth functions, Romberg is **dramatically faster** than Simpson.

---

## Gaussian Quadrature

**Algorithm:** Fixed-point quadrature using optimal abscissas and weights.

**10-point Gauss-Legendre formula:**
```
∫₋₁¹ f(x) dx ≈ ∑ᵢ₌₁¹⁰ wᵢ · f(xᵢ)
```

where {xᵢ} are roots of 10th Legendre polynomial, {wᵢ} are optimal weights.

**Transformed to [a,b]:**
```
∫ᵃᵇ f(x) dx = (b-a)/2 · ∫₋₁¹ f((b-a)/2·t + (b+a)/2) dt
```

### Properties

- **Exactly integrates** polynomials up to degree 2n-1 = 19
- **No adaptive refinement:** Always exactly 10 function evaluations
- **No error estimate:** Must trust the result or compare with another method
- **Very accurate** for smooth functions

### When to Use

✅ **Good for:**
- **Fixed computational budget:** Need exactly N evaluations
- **Smooth well-behaved functions:** No singularities
- **Known sufficient accuracy:** 10 points enough for your problem
- **Benchmarking:** Compare against adaptive methods

❌ **Not recommended:**
- Rough or oscillatory functions (no adaptation!)
- Unknown integrand behavior (no error estimate)
- Variable precision requirements

### Example

```cpp
RealFunction poly([](Real x) {
    return 1 + x + x*x + x*x*x;  // Cubic polynomial
});

// Gauss10 is exact for polynomials up to degree 19!
Real gauss_result = IntegrateGauss10(poly, 0.0, 1.0);

// Analytical: ∫₀¹ (1 + x + x² + x³) dx = [x + x²/2 + x³/3 + x⁴/4]₀¹
//           = 1 + 0.5 + 0.333... + 0.25 = 2.083333...
Real analytical = 1.0 + 0.5 + 1.0/3.0 + 0.25;

std::cout << "Gauss10 result: " << std::setprecision(15) << gauss_result << std::endl;
std::cout << "Analytical:     " << std::setprecision(15) << analytical << std::endl;
std::cout << "Error:          " << abs(gauss_result - analytical) << std::endl;

// Output: Error ≈ 0 (exact to machine precision!)
```

### Accuracy

**For smooth functions:**
- Error ~ O(h²⁰) in practice
- Excellent accuracy from just 10 evaluations
- Competitive with Romberg for many problems

**No error estimate available!**

---

## Integration Result

### IntegrationResult Structure

```cpp
struct IntegrationResult {
    Real value;              // Computed integral value
    Real error_estimate;     // Estimated error
    int iterations;          // Number of refinement iterations
    bool converged;          // Did it converge within max steps?
    
    // Constructor
    IntegrationResult(Real val, Real err, int iters, bool conv);
};
```

### Usage

**Modern interface (recommended):**

```cpp
auto result = IntegrateSimpson(func, 0.0, 1.0, 1e-8);

if (result.converged) {
    std::cout << "Integral = " << result.value 
              << " ± " << result.error_estimate << std::endl;
    std::cout << "Converged in " << result.iterations << " iterations" << std::endl;
} else {
    std::cerr << "WARNING: Did not converge within max iterations!" << std::endl;
    std::cerr << "Best estimate: " << result.value 
              << " (error ~ " << result.error_estimate << ")" << std::endl;
}
```

**Legacy interface (deprecated):**

```cpp
int steps;
Real achieved_prec;
Real value = IntegrateSimpson(func, 0.0, 1.0, &steps, &achieved_prec, 1e-8);

std::cout << "Value: " << value << std::endl;
std::cout << "Steps: " << steps << std::endl;
std::cout << "Achieved precision: " << achieved_prec << std::endl;
```

### Error Estimate Interpretation

The `error_estimate` field contains:
- **Adaptive methods (Trap/Simpson/Romberg):** |Iₖ - Iₖ₋₁| (last improvement)
- **Gauss10:** Not available (no error estimate)

**Convergence criterion:**
```
error_estimate < ε · |value|
```

or both current and previous estimates are zero.

---

## Examples

### Example 1: Comparing All Methods

```cpp
#include "core/Functions.h"
#include "core/Integration.h"

RealFunction f([](Real x) {
    return sin(x) * (1.0 + 0.5 * x * x);
});

// Analytical integral
RealFunction F([](Real x) {
    return x * (-0.5 * x * cos(x) + sin(x));
});

Real a = 0.0, b = 10.0;
Real exact = F(b) - F(a);

std::cout << "Integrating f(x) = sin(x)·(1 + 0.5x²) from " << a << " to " << b << "\n";
std::cout << "Exact value: " << exact << "\n\n";

// Trapezoidal
auto trap = IntegrateTrap(f, a, b, 1e-6);
std::cout << "Trapezoidal: " << trap.value 
          << " (error: " << abs(trap.value - exact)
          << ", iters: " << trap.iterations << ")\n";

// Simpson
auto simp = IntegrateSimpson(f, a, b, 1e-6);
std::cout << "Simpson:     " << simp.value
          << " (error: " << abs(simp.value - exact)
          << ", iters: " << simp.iterations << ")\n";

// Romberg
auto romb = IntegrateRomberg(f, a, b, 1e-10);
std::cout << "Romberg:     " << romb.value
          << " (error: " << abs(romb.value - exact)
          << ", iters: " << romb.iterations << ")\n";

// Gauss10
Real gauss = IntegrateGauss10(f, a, b);
std::cout << "Gauss10:     " << gauss
          << " (error: " << abs(gauss - exact) << ")\n";
```

**Output:**
```
Integrating f(x) = sin(x)·(1 + 0.5x²) from 0 to 10
Exact value: 36.5134

Trapezoidal: 36.5133 (error: 0.0001, iters: 12)
Simpson:     36.5134 (error: 1.2e-07, iters: 7)
Romberg:     36.5134 (error: 3.5e-13, iters: 6)
Gauss10:     36.5134 (error: 2.8e-09)
```

### Example 2: Integration with Singularity (Careful!)

```cpp
// Function with singularity at x=0: f(x) = 1/√x
// ∫₀¹ 1/√x dx = 2 (analytical)

RealFunction singular([](Real x) {
    return 1.0 / sqrt(x);
});

// Avoid exact singularity by integrating from small ε instead of 0
Real eps = 1e-10;
auto result = IntegrateSimpson(singular, eps, 1.0, 1e-6);

std::cout << "∫ 1/√x dx from " << eps << " to 1.0\n";
std::cout << "Numerical: " << result.value << std::endl;
std::cout << "Analytical: 2.0 - 2√" << eps << " ≈ " << 2.0 - 2*sqrt(eps) << std::endl;

// Better: use specialized singularity handling (not shown here)
```

### Example 3: Oscillatory Integrand

```cpp
// Highly oscillatory: f(x) = sin(100x)
// ∫₀^(2π) sin(100x) dx = 0 (analytical)

RealFunction oscillatory([](Real x) {
    return sin(100.0 * x);
});

Real a = 0.0, b = 2.0 * M_PI;

// Need tight tolerance for oscillatory functions
auto result = IntegrateSimpson(oscillatory, a, b, 1e-10);

std::cout << "Oscillatory integral: " << result.value << std::endl;
std::cout << "Expected: 0.0" << std::endl;
std::cout << "Iterations: " << result.iterations << std::endl;

// Many iterations needed to resolve oscillations!
```

### Example 4: Probability Distribution (PDF Integration)

```cpp
// Standard normal PDF: φ(x) = (1/√(2π)) · e^(-x²/2)
RealFunction normal_pdf([](Real x) {
    const Real sqrt_2pi = sqrt(2.0 * M_PI);
    return exp(-0.5 * x * x) / sqrt_2pi;
});

// Compute P(0 ≤ X ≤ 1) for X ~ N(0,1)
auto prob = IntegrateSimpson(normal_pdf, 0.0, 1.0, 1e-10);

std::cout << "P(0 ≤ X ≤ 1) = " << prob.value << std::endl;
// Expected: ≈ 0.3413 (from standard normal table)

// Verify: total probability from -∞ to +∞ should be 1
// (In practice, integrate from -5 to +5 for ~1.0 numerically)
auto total = IntegrateSimpson(normal_pdf, -5.0, 5.0, 1e-8);
std::cout << "Total probability (-5 to +5): " << total.value << std::endl;
// Should be ≈ 0.9999994 (very close to 1.0)
```

### Example 5: Work Done by Variable Force

```cpp
// Physics: Work = ∫ F(x) dx
// Spring force: F(x) = -kx (Hooke's law)
// Work to compress spring from 0 to 0.1m with k = 1000 N/m

Real spring_constant = 1000.0;  // N/m

RealFunction spring_force([spring_constant](Real x) {
    return spring_constant * x;  // Magnitude (unsigned)
});

auto work = IntegrateSimpson(spring_force, 0.0, 0.1, 1e-8);

std::cout << "Work to compress spring 0.1m: " << work.value << " J" << std::endl;
// Analytical: W = ½kx² = ½·1000·0.01 = 5 J

Real analytical_work = 0.5 * spring_constant * 0.1 * 0.1;
std::cout << "Analytical: " << analytical_work << " J" << std::endl;
std::cout << "Error: " << abs(work.value - analytical_work) << " J" << std::endl;
```

---

## Convergence and Accuracy

### Convergence Criteria

**All adaptive methods** stop when:

```
|I_current - I_previous| < ε · |I_current|
```

or when both estimates are zero.

**Maximum iterations:**
- Trapezoidal: `Defaults::TrapezoidIntegrationMaxSteps` (typically 20)
- Simpson: `Defaults::SimpsonIntegrationMaxSteps` (typically 20)
- Romberg: Similar (~15-20 iterations)

**If max iterations reached:**
- `converged = false` in `IntegrationResult`
- Best estimate still returned
- May indicate: oscillatory function, singularity, or too tight tolerance

### Error Estimates

| Method | Error Estimate | Reliability |
|--------|----------------|-------------|
| **Trapezoidal** | \|Iₖ - Iₖ₋₁\| | Good for smooth f(x) |
| **Simpson** | \|Sₖ - Sₖ₋₁\| | Very good for smooth f(x) |
| **Romberg** | \|Rₖⱼ - Rₖ₋₁,ⱼ\| | Excellent for analytic f(x) |
| **Gauss10** | None | N/A (no error estimate) |

**Actual vs Estimated Error:**
- For smooth functions: Estimate is **conservative** (actual error often smaller)
- For non-smooth functions: Estimate may **underestimate** true error

### Precision Limits

**Double precision (`Real = double`):**
- ε ≈ 2.22e-16 (machine epsilon)
- **Best achievable:** ~1e-14 (100× machine epsilon)
- **Typical request:** 1e-6 to 1e-10
- **Don't ask for:** ε < 1e-12 (roundoff errors dominate)

**Setting tolerance:**
```cpp
// ✅ GOOD: Reasonable tolerance
auto result = IntegrateSimpson(func, a, b, 1e-8);

// ⚠️ RISKY: Very tight (may not converge)
auto result = IntegrateSimpson(func, a, b, 1e-14);

// ❌ BAD: Impossible for double precision
auto result = IntegrateSimpson(func, a, b, 1e-20);  // Won't converge!
```

---

## Common Pitfalls

### 1. Wrong Method for Function Type

**Problem:** Using high-order method on non-smooth function

```cpp
// ❌ BAD: Romberg on function with corner
RealFunction abs_func([](Real x) { return abs(x); });
auto bad = IntegrateRomberg(abs_func, -1.0, 1.0, 1e-10);  // Won't converge!

// ✅ GOOD: Trapezoidal for non-smooth
auto good = IntegrateTrap(abs_func, -1.0, 1.0, 1e-6);  // Converges fine
```

### 2. Singularities at Endpoints

**Problem:** Integrand undefined or infinite at a or b

```cpp
// ❌ BAD: Singular at x=0
RealFunction singular([](Real x) { return 1.0 / x; });
auto bad = IntegrateSimpson(singular, 0.0, 1.0, 1e-6);  // Disaster!

// ✅ GOOD: Avoid singularity with small offset
auto good = IntegrateSimpson(singular, 1e-10, 1.0, 1e-6);  // Works

// ✅ BETTER: Use specialized singularity handling (subtract singularity analytically)
```

### 3. Tolerance Too Tight

**Problem:** Requesting precision beyond numerical limits

```cpp
// ❌ BAD: Impossible precision for double
auto bad = IntegrateSimpson(func, 0.0, 1.0, 1e-18);  // Won't converge

// ✅ GOOD: Realistic precision
auto good = IntegrateSimpson(func, 0.0, 1.0, 1e-10);  // Achievable
```

### 4. Ignoring Convergence Flag

**Problem:** Using result without checking if it converged

```cpp
// ❌ BAD: Blindly trust result
auto result = IntegrateSimpson(difficult_func, a, b, 1e-10);
Real integral = result.value;  // May not have converged!

// ✅ GOOD: Check convergence
auto result = IntegrateSimpson(difficult_func, a, b, 1e-10);
if (!result.converged) {
    std::cerr << "WARNING: Integration did not converge!" << std::endl;
    std::cerr << "Error estimate: " << result.error_estimate << std::endl;
    // Decide: relax tolerance, switch method, or subdivide interval
}
```

### 5. Oscillatory Functions Without Enough Subdivisions

**Problem:** Undersampling oscillations

```cpp
// Oscillates 100 times in [0, 2π]
RealFunction fast_osc([](Real x) { return sin(100 * x); });

// ❌ MAY FAIL: Default tolerance may not resolve oscillations
auto iffy = IntegrateSimpson(fast_osc, 0.0, 2*M_PI, 1e-6);

// ✅ BETTER: Tighter tolerance forces more subdivisions
auto better = IntegrateSimpson(fast_osc, 0.0, 2*M_PI, 1e-12);

// ✅ BEST: Manually subdivide interval to match oscillation period
Real period = 2*M_PI / 100.0;
// Integrate over each period and sum
```

### 6. Infinite or Very Large Intervals

**Problem:** Integrating from -∞ to +∞ or very large range

```cpp
// ❌ BAD: Cannot literally use infinity
// auto bad = IntegrateSimpson(gaussian, -INFINITY, +INFINITY, 1e-6);

// ✅ GOOD: Use practical limits (function negligible beyond ±5 for Gaussian)
auto good = IntegrateSimpson(gaussian, -5.0, 5.0, 1e-8);

// ✅ BETTER: Use variable transformation for infinite intervals
// Transform: x = tan(θ), dx = sec²(θ) dθ, θ ∈ [-π/2, π/2]
// ∫_{-∞}^{+∞} f(x) dx = ∫_{-π/2}^{π/2} f(tan(θ)) · sec²(θ) dθ
```

---

## See Also

**Related Documentation:**
- [Multidim_integration.md](Multidim_integration.md) - 2D and 3D integration
- [Functions.md](Functions.md) - Function interfaces
- [Derivation.md](Derivation.md) - Numerical derivatives

**Advanced Integration:**
- [../algorithms/MonteCarlo_integration.md](../algorithms/MonteCarlo_integration.md) - Monte Carlo methods for high dimensions
- Adaptive Simpson (ASR) - Recursive adaptive refinement
- Gaussian quadrature (higher orders) - 20, 32, 64-point rules

**Base Classes:**
- [Vector.md](../base/Vector.md) - For storing tabulated data
- [Functions.md](Functions.md) - Creating integrable function objects

**Theory References:**
- Numerical Recipes §4.2-4.4 - Trapezoidal, Simpson, Romberg
- Numerical Recipes §4.5 - Gaussian Quadratures
- Abramowitz & Stegun - Tables of Gaussian Quadrature Points

---

## Runnable Examples

| Example | Source File | Description |
|---------|------------|-------------|
| Integration Demo | [docs_demo_integration.cpp](../../src/docs_demos/docs_demo_integration.cpp) | 1D integration methods |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/docs_demos/Release/MML_DocsApp
```
