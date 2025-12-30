# Root Finding Algorithms

Comprehensive guide to finding zeros of real-valued functions using bracketing and derivative-based methods.

## Overview

Root finding algorithms solve the fundamental problem: **given function f(x), find x* such that f(x*) = 0**. This appears throughout:
- **Physics**: Equilibrium points, resonant frequencies, orbital mechanics
- **Engineering**: Circuit design, structural analysis, control systems
- **Economics**: Break-even points, Nash equilibria, pricing models
- **Numerical analysis**: Implicit equations, nonlinear systems preparation

The library provides **4 robust methods** spanning bracketing (guaranteed convergence) to Newton-Raphson (fast but needs care).

## Quick Reference

| Method | Convergence | Requires Bracket | Needs Derivative | Speed | Robustness | Use Case |
|--------|-------------|------------------|------------------|-------|------------|----------|
| **BracketRoot** | N/A | No (creates bracket) | No | Fast search | ⭐⭐⭐ | **Find initial bracket** |
| **FindRootBrackets** | N/A | No (finds all) | No | Systematic | ⭐⭐⭐⭐ | **Locate all roots** |
| **FindRootBisection** | Linear | Yes | No | Slow | ⭐⭐⭐⭐⭐ | **Guaranteed convergence** |
| **FindRootFalsePosition** | Superlinear (~1.6) | Yes | No | Medium | ⭐⭐⭐⭐ | **Faster than bisection** |
| **FindRootSecant** | Superlinear (~1.618) | No | No | Fast | ⭐⭐⭐ | **No derivative Newton** |
| **FindRootNewton** | Quadratic | Optional | Yes (numerical) | **Fastest** | ⭐⭐ | **High-precision refinement** |
| **FindRootRidders** | Quadratic | Yes | No | Fast | ⭐⭐⭐⭐ | **Robust fast bracketing** |
| **FindRootBrent** | Superlinear-Quadratic | Yes | No | Fast | ⭐⭐⭐⭐⭐ | **Best all-around** |

**Key Insight**: Always bracket first (BracketRoot/FindRootBrackets), then refine with your method of choice. Brent is the recommended default.

## Mathematical Background

### Root Finding Problem

Given continuous function f: ℝ → ℝ, find x* where:
```
f(x*) = 0
```

**Challenges**:
- Multiple roots or no roots may exist
- Numerical precision limits (~10⁻¹⁵ for double precision)
- Ill-conditioned functions (near-horizontal near root)
- Discontinuities or infinite derivatives

### Intermediate Value Theorem

**Foundation for bracketing**: If f is continuous on [a,b] and f(a)·f(b) < 0, then ∃x* ∈ (a,b) such that f(x*) = 0.

**Critical**: Sign change guarantees root existence in interval.

### Convergence Rates

Let eₙ = |xₙ - x*| be error at iteration n.

**Linear convergence**: eₙ₊₁ ≈ C·eₙ where 0 < C < 1
- Bisection: C = 0.5 (halves interval each step)
- Predictable: n steps → 2⁻ⁿ smaller interval

**Quadratic convergence**: eₙ₊₁ ≈ C·eₙ²
- Newton-Raphson: Near root, error squares each iteration
- Example: 10⁻² → 10⁻⁴ → 10⁻⁸ → 10⁻¹⁶ (4 steps!)
- **Requirement**: Good initial guess, f'(x*) ≠ 0

## Methods

### BracketRoot - Find Initial Bracket

**Purpose**: Expand initial guess range until sign change detected.

**Algorithm**:
```cpp
static bool BracketRoot(const IRealFunction& func, Real& x1, Real& x2, int MaxTry = 50)
```

**Strategy**:
1. Start with [x₁, x₂]
2. Evaluate f(x₁), f(x₂)
3. If f(x₁)·f(x₂) < 0: **Done!** Root bracketed
4. Else: Expand interval geometrically (scaleFactor = 1.6)
   - Expand in direction of smaller |f|
   - x₁ ← x₁ + 1.6·(x₁ - x₂) if |f(x₁)| < |f(x₂)|
   - x₂ ← x₂ + 1.6·(x₂ - x₁) otherwise
5. Repeat up to MaxTry iterations

**Returns**: 
- `true` if bracket found (x₁, x₂ updated with bracket)
- `false` if search failed (no sign change found)

**Use When**:
- Have rough estimate of root location
- Need single root bracket for refinement
- Exploratory search in unknown region

**Limitations**:
- May miss roots if expanding in wrong direction
- Not systematic (depends on initial guess quality)
- Can bracket wrong root if multiple exist

---

### FindRootBrackets - Systematic Root Search

**Purpose**: Systematically search interval, find ALL root brackets.

**Algorithm**:
```cpp
static int FindRootBrackets(const IRealFunction& func, const Real x1, const Real x2, 
                           const int numPoints, Vector<Real>& xb1, Vector<Real>& xb2)
```

**Strategy**:
1. Subdivide [x₁, x₂] into numPoints equal segments
2. Evaluate function at each segment endpoint
3. Check for sign changes: f(xᵢ)·f(xᵢ₊₁) ≤ 0
4. Store each bracket pair in (xb1, xb2) vectors
5. Dynamic resize if more roots found than expected

**Returns**: Number of brackets found (numRoots)

**Output**: 
- `xb1[0..numRoots-1]`: Left bracket endpoints
- `xb2[0..numRoots-1]`: Right bracket endpoints

**Resolution**: dx = (x₂ - x₁) / numPoints

**Use When**:
- Need to find ALL roots in interval
- Function behavior unknown
- Preparing for batch processing of roots
- Verification (counting expected roots)

**Sampling Strategy**:
- **numPoints = 100**: Quick survey, may miss narrow features
- **numPoints = 1000**: Detailed search, captures most roots
- **numPoints ~ 10/δx**: Where δx = expected root spacing

**Limitations**:
- Can miss roots if spacing too coarse
- Sign changes only (misses tangent roots f(x*) = 0, f'(x*) = 0)
- Brackets may contain multiple roots

---

### FindRootBisection - Guaranteed Convergence

**Purpose**: Robust, guaranteed root finding via interval halving.

**Algorithm**:
```cpp
static Real FindRootBisection(const IRealFunction& func, Real x1, Real x2, Real xacc)
```

**Classical Bisection**:
```
Given: [a, b] with f(a)·f(b) < 0
Repeat:
  c = (a + b)/2
  Evaluate f(c)
  If f(a)·f(c) < 0: b ← c  (root in [a,c])
  Else:             a ← c  (root in [c,b])
Until: |b - a| < xacc
```

**Implementation Details**:
- Maintains rtb (root bound) and dx (interval width)
- Reduces dx by factor 2 each iteration
- **Safety checks**: Validates f(x₁)·f(x₂) < 0 initially
- **Robustness**: Detects NaN/Inf function values
- **Max iterations**: Defaults::BisectionMaxSteps (typically 100)

**Convergence**:
- **Rate**: Linear, exact reduction factor 0.5
- **Iterations needed**: n ≈ log₂((x₂-x₁)/xacc)
- **Example**: [0, 1] to 10⁻¹² accuracy → ~40 iterations

**Error Bound**: After n iterations:
```
|xₙ - x*| ≤ (x₂ - x₁) / 2ⁿ
```

**Guarantees**:
- ✅ **Always converges** if root bracketed
- ✅ **Predictable** iteration count
- ✅ **No derivative** needed
- ✅ **Works with discontinuous f** (if sign change exists)

**Use When**:
- **Robustness critical** (production code, safety systems)
- Function evaluation cheap
- Moderate accuracy needed (10⁻⁶ to 10⁻¹⁰)
- Derivative unavailable or unreliable

**Limitations**:
- ❌ **Slow** compared to Newton (40× more iterations typical)
- ❌ Requires valid bracket
- ❌ Cannot exploit function smoothness

---

### FindRootNewton - Fast Quadratic Convergence

**Purpose**: Rapid root refinement using function derivative information.

**Algorithm**:
```cpp
static Real FindRootNewton(const IRealFunction& func, Real x1, Real x2, Real xacc)
```

**Newton-Raphson Method**:
```
Given: Initial guess x₀ ∈ [x₁, x₂]
Repeat:
  xₙ₊₁ = xₙ - f(xₙ)/f'(xₙ)
Until: |xₙ₊₁ - xₙ| < xacc
```

**Geometric Interpretation**: 
- Linearize f at current point: f(x) ≈ f(xₙ) + f'(xₙ)·(x - xₙ)
- Find where linearization crosses x-axis
- Use that as next guess

**Derivative Computation**: Uses `Derivation::NDer4()` (4-point numerical)
```cpp
f'(x) ≈ [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)
```
Accuracy: O(h⁴), automatically chooses h

**Implementation Safeguards**:

1. **Bracket checking**: Ensures xₙ stays in [x₁, x₂]
   - Throws if Newton step jumps outside
   - Prevents wild divergence

2. **Derivative validation**:
   ```cpp
   if (|f'(x)| < √ε) throw "Near-zero derivative"
   ```
   - Prevents division by tiny numbers
   - Detects horizontal tangent situations

3. **NaN/Inf detection**: Checks f'(x) and step dx for non-finite values

4. **Max iterations**: Defaults::NewtonRaphsonMaxSteps (typically 100)

**Convergence Analysis**:

**Near root** (x close to x*):
```
eₙ₊₁ ≈ |f''(x*)|/(2|f'(x*)|) · eₙ²
```

**Example convergence** (good conditions):
```
e₀ = 10⁻²  →  e₁ ≈ 10⁻⁴  →  e₂ ≈ 10⁻⁸  →  e₃ ≈ 10⁻¹⁶
```
**3 iterations** achieve machine precision!

**When Convergence Fails**:
- **f'(x*) = 0**: Local extremum, no convergence (division by small f')
- **Poor initial guess**: May diverge or find wrong root
- **Multiple roots**: f'(x) = 0 between roots causes issues
- **Inflection near root**: f''(x) large → slow convergence

**Use When**:
- ✅ **Speed critical** (real-time systems, optimization loops)
- ✅ Good initial guess available (from bracketing)
- ✅ High precision needed (10⁻¹² to 10⁻¹⁵)
- ✅ Smooth function with f'(x*) ≠ 0
- ✅ Can afford derivative calculation

**Limitations**:
- ❌ **No convergence guarantee**
- ❌ Sensitive to initial guess
- ❌ Fails at extrema (f'(x*) ≈ 0)
- ❌ Can diverge (oscillate or escape bracket)

---

### FindRootFalsePosition - Faster Bracketing

**Purpose**: Root finding via linear interpolation while maintaining bracket.

**Algorithm**:
```cpp
static Real FindRootFalsePosition(const IRealFunction& func, Real x1, Real x2, Real xacc)
```

**Method** (Regula Falsi):
1. Start with bracket [x₁, x₂] where f(x₁)·f(x₂) < 0
2. Linear interpolation: x_new = (x₁·f₂ - x₂·f₁) / (f₂ - f₁)
3. Evaluate f(x_new)
4. Update bracket: replace x₁ or x₂ with x_new to maintain sign change
5. Continue until |f(x_new)| < xacc or bracket width < xacc

**Convergence**:
- **Superlinear**: Order ~1.618 (golden ratio φ)
- Faster than bisection, slower than Newton
- Guaranteed convergence (maintains bracket)

**Use When**:
- ✅ Bisection is too slow
- ✅ Newton is unreliable (poor initial guess)
- ✅ Function has discontinuous derivatives
- ✅ Need guaranteed convergence with better speed

**Limitations**:
- ❌ Can be slow with one-sided convergence (one endpoint stays fixed)
- ❌ Slower than Newton or secant for smooth functions

---

### FindRootSecant - Derivative-Free Newton

**Purpose**: Newton-like method without explicit derivative calculation.

**Algorithm**:
```cpp
static Real FindRootSecant(const IRealFunction& func, Real x1, Real x2, Real xacc)
```

**Method**:
```
Given: Two initial guesses x₀, x₁
Iteration: xₙ₊₁ = xₙ - f(xₙ)·(xₙ - xₙ₋₁) / (f(xₙ) - f(xₙ₋₁))
```

**Derivative Approximation**: 
- Approximates f'(xₙ) ≈ (f(xₙ) - f(xₙ₋₁)) / (xₙ - xₙ₋₁)
- Only **one function evaluation per iteration** (vs 2 for numerical Newton)

**Convergence**:
- **Superlinear**: Order φ = (1 + √5)/2 ≈ 1.618
- Faster than bisection, slightly slower than Newton
- **Example**: 5-10 iterations for machine precision

**Use When**:
- ✅ Derivatives expensive or unavailable
- ✅ Have two good initial guesses
- ✅ Smooth function behavior
- ✅ Newton is overkill for the problem

**Limitations**:
- ❌ May diverge with poor initial guesses
- ❌ No bracket maintenance (can jump away from root)
- ❌ Fails if f(xₙ) = f(xₙ₋₁)

---

### FindRootRidders - Quadratic Without Derivatives

**Purpose**: Quadratic convergence without requiring derivatives.

**Algorithm**:
```cpp
static Real FindRootRidders(const IRealFunction& func, Real x1, Real x2, Real xacc)
```

**Method** (Exponential Interpolation):
1. Start with bracket [x₁, x₂]
2. Compute midpoint: xₘ = (x₁ + x₂) / 2
3. Evaluate: f₁ = f(x₁), fₘ = f(xₘ), f₂ = f(x₂)
4. Compute: s = √(fₘ² - f₁·f₂)
5. New estimate: x₄ = xₘ + (xₘ - x₁)·sign(f₁ - f₂)·fₘ / s
6. Update bracket with x₄
7. Continue until convergence

**Mathematical Basis**: Fits exponential function through three points:
- f(x) ≈ A·e^(λx) + B·e^(-λx)
- More sophisticated than linear interpolation

**Convergence**:
- **Quadratic**: Same as Newton!
- Error: eₙ₊₁ ≈ C·eₙ²
- 2-3× faster than false position
- **Guaranteed convergence** (maintains bracket)

**Use When**:
- ✅ **Best general-purpose choice**
- ✅ Need quadratic convergence without derivatives
- ✅ Want robustness AND speed
- ✅ Function may have difficult curvature

**Advantages over Newton**:
- ✅ No derivative needed
- ✅ More robust (always converges)
- ✅ Often outperforms Newton on difficult functions

---

### FindRootBrent - The Gold Standard

**Purpose**: Industry-standard hybrid method combining inverse quadratic interpolation, secant, and bisection.

**Algorithm**:
```cpp
static Real FindRootBrent(const IRealFunction& func, Real x1, Real x2, Real xacc)
```

**Method** (Adaptive Hybrid):
- **Primary**: Inverse quadratic interpolation (superquadratic near root)
- **Secondary**: Secant method (superlinear)
- **Fallback**: Bisection (guaranteed progress)
- Automatically switches between methods based on progress

**Key Features**:
- Falls back to bisection if:
  - Interpolated point outside bracket
  - Step too small (stagnation)
  - Insufficient bracket reduction
- Maintains |f(b)| ≤ |f(a)| invariant
- Tracks previous steps to detect poor progress

**Convergence**:
- **Superlinear to quadratic** depending on function behavior
- Typical: 6-12 iterations for double precision
- **Guaranteed convergence** (bracket always maintained)
- Root guaranteed to be in original bracket

**Use When**:
- ✅ **Default choice for production code**
- ✅ Both speed and reliability matter
- ✅ General-purpose root finding
- ✅ Unsure which method to use

**Why Brent?**:
- 📚 Widely trusted: scipy, GNU GSL, MATLAB, Boost
- 🏆 Best all-around performance
- 🔒 Extremely robust (never fails if root exists)
- 🚀 Fast: comparable to Newton without derivatives

**Comparison**:
- vs Bisection: Much faster, same reliability
- vs Newton: No derivatives, more robust, similar speed
- vs Ridders: Similar performance, more widely used and tested

---

## Integration with Other Modules

### IRealFunction Interface
```cpp
namespace MML {
  class IRealFunction {
    virtual Real operator()(Real x) const = 0;
  };
}
```

**All root finding functions accept IRealFunction**:
- Function objects
- Lambda expressions (via wrapper)
- Functors with state

**Example with lambda**:
```cpp
auto f = [](Real x) { return x*x - 2.0; };
// Wrap in IRealFunction implementation
Real root = RootFinding::FindRootBisection(funcWrapper, 0.0, 2.0, 1e-10);
```

### Derivation Module

**FindRootNewton** uses `Derivation::NDer4()`:
```cpp
Real df = Derivation::NDer4(func, x);
```

**Automatic step size**: Adapts to function scale, typically h ~ x·∛ε

**Accuracy**: O(h⁴) error, sufficient for Newton convergence

### Statistics Module

**Root quality assessment**:
```cpp
Vector<Real> residuals(n);
for (int i = 0; i < n; i++)
  residuals[i] = func(roots[i]);  // Should be near zero

Real maxError = Statistics::MaxAbs(residuals);
Real avgError = Statistics::Mean(residuals.Apply([](Real r) { return std::abs(r); }));
```

### Function Analysis

**Combine with FunctionsAnalyzer**:
```cpp
// Find all roots systematically
Vector<Real> xb1, xb2;
int n = RootFinding::FindRootBrackets(f, -10.0, 10.0, 1000, xb1, xb2);

Vector<Real> roots(n);
for (int i = 0; i < n; i++)
  roots[i] = RootFinding::FindRootNewton(f, xb1[i], xb2[i], 1e-12);

// Classify roots using FunctionsAnalyzer
RealFunctionAnalyzer analyzer(f);
for (auto root : roots) {
  Real fpp = analyzer.GetSecondDerivative(root);
  // fpp > 0: local minimum, fpp < 0: local maximum
}
```

## Examples

### Example 1: Basic Root Finding - Quadratic

Find root of f(x) = x² - 2 (√2 ≈ 1.414213562373095):

```cpp
#include "algorithms/RootFinding.h"

class QuadraticFunc : public IRealFunction {
public:
  Real operator()(Real x) const override { return x*x - 2.0; }
};

void Example1() {
  QuadraticFunc f;
  
  // Method 1: Direct bisection if bracket known
  Real root1 = RootFinding::FindRootBisection(f, 1.0, 2.0, 1e-10);
  std::cout << "Bisection: " << std::setprecision(15) << root1 << "\n";
  // Output: 1.41421356237310
  
  // Method 2: Newton-Raphson for speed
  Real root2 = RootFinding::FindRootNewton(f, 1.0, 2.0, 1e-14);
  std::cout << "Newton:    " << std::setprecision(15) << root2 << "\n";
  // Output: 1.41421356237310 (typically 4-5 iterations)
  
  // Verify
  std::cout << "f(root1) = " << f(root1) << "\n";  // ~1e-11
  std::cout << "f(root2) = " << f(root2) << "\n";  // ~1e-15
}
```

### Example 2: Bracket Expansion

Find root when only rough estimate available:

```cpp
class TranscendentalFunc : public IRealFunction {
public:
  Real operator()(Real x) const override { 
    return std::cos(x) - x;  // Root at x ≈ 0.739085
  }
};

void Example2() {
  TranscendentalFunc f;
  
  // Guess: root somewhere near 1.0
  Real x1 = 0.5, x2 = 1.0;
  
  if (RootFinding::BracketRoot(f, x1, x2)) {
    std::cout << "Bracket found: [" << x1 << ", " << x2 << "]\n";
    // Output: [0.3125, 1.0] (or similar)
    
    // Now refine
    Real root = RootFinding::FindRootNewton(f, x1, x2, 1e-12);
    std::cout << "Root: " << std::setprecision(12) << root << "\n";
    // Output: 0.739085133215
    
    // Verify: cos(root) should equal root
    std::cout << "cos(root) = " << std::cos(root) << "\n";
    std::cout << "root      = " << root << "\n";
  } else {
    std::cout << "Bracket expansion failed\n";
  }
}
```

### Example 3: Finding All Roots - Polynomial

Find all real roots of P(x) = (x-1)(x-3)(x-5) = x³ - 9x² + 23x - 15:

```cpp
class CubicPoly : public IRealFunction {
public:
  Real operator()(Real x) const override {
    return x*x*x - 9*x*x + 23*x - 15;
  }
};

void Example3() {
  CubicPoly f;
  
  // Search interval [0, 6] with 600 sample points
  Vector<Real> xb1, xb2;
  int numBrackets = RootFinding::FindRootBrackets(f, 0.0, 6.0, 600, xb1, xb2);
  
  std::cout << "Found " << numBrackets << " brackets:\n";
  
  Vector<Real> roots(numBrackets);
  for (int i = 0; i < numBrackets; i++) {
    // Refine each bracket with Newton
    roots[i] = RootFinding::FindRootNewton(f, xb1[i], xb2[i], 1e-12);
    std::cout << "Root " << i+1 << ": " << std::setprecision(12) << roots[i];
    std::cout << "  f(x) = " << f(roots[i]) << "\n";
  }
  
  // Expected output:
  // Root 1: 1.000000000000  f(x) = 0e+00
  // Root 2: 3.000000000000  f(x) = 0e+00
  // Root 3: 5.000000000000  f(x) = 0e+00
}
```

### Example 4: Resonant Frequency (Physics Application)

Find resonant frequency of damped harmonic oscillator:

Impedance: Z(ω) = √[R² + (ωL - 1/(ωC))²]

Resonance at dZ/dω = 0, leads to: ω²LC - 1 = 0

```cpp
class ResonanceEquation : public IRealFunction {
  Real L, C;
public:
  ResonanceEquation(Real inductance, Real capacitance) 
    : L(inductance), C(capacitance) {}
  
  Real operator()(Real omega) const override {
    return omega*omega * L*C - 1.0;
  }
};

void Example4() {
  Real L = 10e-3;   // 10 mH
  Real C = 100e-9;  // 100 nF
  
  ResonanceEquation f(L, C);
  
  // Theoretical: ω₀ = 1/√(LC)
  Real omega_theory = 1.0 / std::sqrt(L * C);
  std::cout << "Theoretical ω₀ = " << omega_theory << " rad/s\n";
  
  // Numerical: bracket near expected value
  Real w1 = omega_theory * 0.5;
  Real w2 = omega_theory * 1.5;
  
  Real omega_numerical = RootFinding::FindRootBisection(f, w1, w2, 1e-6);
  std::cout << "Numerical ω₀   = " << omega_numerical << " rad/s\n";
  
  // Convert to frequency: f = ω/(2π)
  Real freq_Hz = omega_numerical / (2 * Constants::PI);
  std::cout << "Resonant frequency = " << freq_Hz << " Hz\n";
  // Output: ~15915.5 Hz (15.9 kHz)
}
```

### Example 5: Convergence Comparison

Compare convergence rates of Bisection vs Newton:

```cpp
class TestFunc : public IRealFunction {
public:
  Real operator()(Real x) const override {
    return std::exp(x) - 3*x*x;  // Root near x ≈ 3.733
  }
};

void Example5() {
  TestFunc f;
  Real x1 = 3.0, x2 = 4.0;
  
  // Custom bisection with iteration tracking
  std::cout << "BISECTION:\n";
  Real a = x1, b = x2;
  for (int iter = 1; iter <= 20; iter++) {
    Real c = 0.5*(a + b);
    Real fc = f(c);
    
    if (f(a) * fc < 0) b = c;
    else               a = c;
    
    if (iter % 5 == 0)
      std::cout << "Iter " << iter << ": x = " << std::setprecision(12) << c 
                << "  f(x) = " << fc << "\n";
  }
  
  // Custom Newton with tracking
  std::cout << "\nNEWTON:\n";
  Real x = 0.5*(x1 + x2);
  for (int iter = 1; iter <= 6; iter++) {
    Real fx = f(x);
    Real dfx = Derivation::NDer4(f, x);
    x -= fx/dfx;
    
    std::cout << "Iter " << iter << ": x = " << std::setprecision(12) << x 
              << "  f(x) = " << f(x) << "\n";
  }
  
  /* Output:
     BISECTION:
     Iter 5:  x = 3.71875      f(x) = 0.7xxx
     Iter 10: x = 3.73291      f(x) = 0.02xxx
     Iter 15: x = 3.73322      f(x) = 0.0007xxx
     Iter 20: x = 3.73324      f(x) = 2e-5
     
     NEWTON:
     Iter 1:  x = 3.72541      f(x) = -0.15
     Iter 2:  x = 3.73323      f(x) = 0.003
     Iter 3:  x = 3.73324      f(x) = 1.2e-6
     Iter 4:  x = 3.73324      f(x) = 2e-13
     
     Newton achieves in 4 iterations what Bisection needs 20+ for!
  */
}
```

### Example 6: Multiple Roots with Visualization

Systematically find and classify all roots of a complex function:

```cpp
class ComplexFunc : public IRealFunction {
public:
  Real operator()(Real x) const override {
    return std::sin(x) * std::exp(-x*x/10) - 0.3;
  }
};

void Example6() {
  ComplexFunc f;
  
  // Survey wide range
  Vector<Real> xb1, xb2;
  int n = RootFinding::FindRootBrackets(f, -5.0, 5.0, 1000, xb1, xb2);
  
  std::cout << "Found " << n << " sign changes\n\n";
  
  // Refine each and report
  for (int i = 0; i < n; i++) {
    // Use Newton for refinement (faster)
    Real root = RootFinding::FindRootNewton(f, xb1[i], xb2[i], 1e-12);
    
    // Compute derivative at root for classification
    Real df = Derivation::NDer4(f, root);
    std::string crossing = (df > 0) ? "ascending" : "descending";
    
    std::cout << "Root " << i+1 << ": x = " << std::setw(10) << root 
              << "  (" << crossing << " crossing)\n";
    std::cout << "  Residual: f(x) = " << f(root) << "\n";
    std::cout << "  Initial bracket: [" << xb1[i] << ", " << xb2[i] << "]\n\n";
  }
}
```

### Example 7: Orbital Mechanics - Kepler's Equation

Solve Kepler's equation for eccentric anomaly E:
```
M = E - e·sin(E)
```
where M = mean anomaly, e = eccentricity

```cpp
class KeplerEquation : public IRealFunction {
  Real M, e;
public:
  KeplerEquation(Real meanAnomaly, Real eccentricity)
    : M(meanAnomaly), e(eccentricity) {}
  
  Real operator()(Real E) const override {
    return E - e * std::sin(E) - M;
  }
};

void Example7() {
  Real M = Constants::PI / 3;  // Mean anomaly = 60°
  Real e = 0.3;                 // Moderate eccentricity
  
  KeplerEquation kepler(M, e);
  
  // E must be in range [0, 2π], bracket around M
  Real E = RootFinding::FindRootNewton(kepler, M - 1.0, M + 1.0, 1e-14);
  
  std::cout << "Mean anomaly M = " << M << " rad (" << M*180/Constants::PI << "°)\n";
  std::cout << "Eccentricity e = " << e << "\n";
  std::cout << "Eccentric anomaly E = " << E << " rad (" << E*180/Constants::PI << "°)\n";
  
  // Verify solution
  std::cout << "\nVerification: M = E - e·sin(E)\n";
  std::cout << "  LHS (M) = " << M << "\n";
  std::cout << "  RHS     = " << (E - e*std::sin(E)) << "\n";
  std::cout << "  Error   = " << std::abs(M - (E - e*std::sin(E))) << "\n";
}
```

### Example 8: Implicit Equation - Circle-Line Intersection

Find intersection of circle x² + y² = r² with line y = mx + b:

Substituting: x² + (mx + b)² = r²  →  (1+m²)x² + 2mbx + (b²-r²) = 0

```cpp
class CircleLineIntersection : public IRealFunction {
  Real m, b, r;
public:
  CircleLineIntersection(Real slope, Real intercept, Real radius)
    : m(slope), b(intercept), r(radius) {}
  
  Real operator()(Real x) const override {
    Real y = m*x + b;
    return x*x + y*y - r*r;  // = 0 at intersection
  }
};

void Example8() {
  Real m = 0.5;   // Line slope
  Real b = 1.0;   // y-intercept
  Real r = 2.0;   // Circle radius
  
  CircleLineIntersection f(m, b, r);
  
  // Find both intersection points
  Vector<Real> xb1, xb2;
  int n = RootFinding::FindRootBrackets(f, -3.0, 3.0, 100, xb1, xb2);
  
  std::cout << "Line: y = " << m << "x + " << b << "\n";
  std::cout << "Circle: x² + y² = " << r*r << "\n";
  std::cout << "Intersections: " << n << "\n\n";
  
  for (int i = 0; i < n; i++) {
    Real x = RootFinding::FindRootBisection(f, xb1[i], xb2[i], 1e-10);
    Real y = m*x + b;
    
    std::cout << "Point " << i+1 << ": (" << x << ", " << y << ")\n";
    
    // Verify both equations
    Real circleCheck = x*x + y*y;
    Real lineCheck = y - (m*x + b);
    std::cout << "  Circle check: x²+y² = " << circleCheck 
              << " (should be " << r*r << ")\n";
    std::cout << "  Line check: y-mx-b = " << lineCheck << "\n\n";
  }
}
```

## Best Practices

### Solver Selection Strategy

**Decision Tree**:

```
Do you have a bracket [a,b] with f(a)·f(b) < 0?
├─ YES: Root guaranteed in [a,b]
│   ├─ Need guaranteed convergence? → FindRootBisection
│   ├─ Want speed & high precision? → FindRootNewton (with bracket check)
│   └─ Don't know which? → Start Newton, fall back to Bisection if issues
│
└─ NO: Need to find bracket first
    ├─ Have good guess? → BracketRoot (expand from guess)
    ├─ Unknown region? → FindRootBrackets (systematic search)
    └─ Multiple roots expected? → FindRootBrackets (finds all)
```

### Tolerance Selection

| Application | Accuracy | Method | Tolerance | Iterations (typical) |
|-------------|----------|--------|-----------|---------------------|
| Engineering rough | 10⁻⁴ | Bisection | 1e-4 | ~13 (on [0,1]) |
| Scientific standard | 10⁻⁸ | Newton | 1e-8 | ~5 |
| High precision | 10⁻¹² | Newton | 1e-12 | ~6 |
| Machine limit | 10⁻¹⁵ | Newton | 1e-14 | ~7 |

**Rule of thumb**: 
- Bisection: n_iter ≈ log₂((b-a)/tol)
- Newton: n_iter ≈ log₂(log(1/tol)) + constant

### Numerical Stability

**Function Scaling**:
```cpp
// Bad: Different scales
Real f_bad(Real x) { return 1e10 * x - 1e-5; }  // Root at x ≈ 1e-15

// Good: Normalized
Real f_good(Real x) { return x - 1e-15; }  // Direct representation
```

**Ill-Conditioned Roots** (f'(x*) ≈ 0):
```cpp
// Problematic: Multiple root f(x) = (x-1)²
// f'(1) = 0 → Newton fails

// Solution: Reformulate
// Instead of f(x) = (x-1)², use g(x) = x - 1
// Or deflate: f(x)/(x-1) if first root known
```

**Derivative Issues**:
```cpp
// If f'(x) has discontinuities, Newton can jump erratically
// Solution: Use bisection near discontinuities, Newton in smooth regions
```

### Common Pitfalls

❌ **Pitfall 1**: Using Newton without bracket
```cpp
// DON'T:
Real root = FindRootNewton(f, x1, x2, 1e-10);  // May diverge!

// DO:
Real root;
if (RootFinding::BracketRoot(f, x1, x2))
  root = RootFinding::FindRootNewton(f, x1, x2, 1e-10);
else
  std::cerr << "Could not bracket root\n";
```

❌ **Pitfall 2**: Coarse sampling missing roots
```cpp
// DON'T (may miss narrow peaks):
Vector<Real> xb1, xb2;
int n = FindRootBrackets(f, 0, 10, 10, xb1, xb2);  // Only 10 points!

// DO (adequate resolution):
int n = FindRootBrackets(f, 0, 10, 1000, xb1, xb2);  // 1000 points
```

❌ **Pitfall 3**: Ignoring function evaluation cost
```cpp
// If f(x) is expensive (numerical integration, simulation):
// Prefer Newton (few evaluations) over Bisection (many evaluations)
// Even if Newton needs derivative (2× cost per iteration, but 10× fewer iterations)
```

❌ **Pitfall 4**: Not verifying convergence
```cpp
// ALWAYS check residual:
Real root = FindRootNewton(f, a, b, 1e-10);
Real residual = std::abs(f(root));
if (residual > 1e-8)
  std::cerr << "Warning: Poor convergence, residual = " << residual << "\n";
```

✅ **Best Practice Pattern**:
```cpp
// 1. Bracket
Real x1 = initial_guess - delta;
Real x2 = initial_guess + delta;
if (!RootFinding::BracketRoot(f, x1, x2)) {
  // Fallback: Systematic search
  Vector<Real> xb1, xb2;
  int n = RootFinding::FindRootBrackets(f, search_min, search_max, 1000, xb1, xb2);
  if (n == 0) throw std::runtime_error("No roots found");
  x1 = xb1[0]; x2 = xb2[0];  // Use first bracket
}

// 2. Refine
Real root;
try {
  root = RootFinding::FindRootNewton(f, x1, x2, tolerance);
} catch (const RootFindingError&) {
  // Newton failed, fall back to bisection
  root = RootFinding::FindRootBisection(f, x1, x2, tolerance * 10);
}

// 3. Verify
Real residual = std::abs(f(root));
assert(residual < 1e-6);  // Application-specific threshold
```

## Performance Considerations

### Complexity Analysis

| Method | Per-Iteration | Total to 10⁻¹⁰ | Function Evals | Derivative Evals |
|--------|---------------|----------------|----------------|------------------|
| **BracketRoot** | O(1) | ~30 expansions | ~60 | 0 |
| **FindRootBrackets** | O(n) | O(n) | n | 0 |
| **Bisection** | O(1) | ~33 iterations | 33 | 0 |
| **Newton** | O(1) | ~5 iterations | 5 | 5 |

**Key Insight**: Newton is 6-7× fewer iterations, but each iteration needs derivative (~5 extra function evaluations for NDer4). **Total work similar** if f(x) cheap, **Newton wins** if f(x) expensive.

### When Function Evaluation is Expensive

Examples:
- f(x) = ∫₀ˣ g(t) dt (numerical integration each call)
- f(x) = solution to ODE at time x
- f(x) = result of finite element simulation

**Strategy**: Minimize total function evaluations
- Prefer Newton (even with numerical derivative overhead)
- Cache function values if possible
- Use adaptive tolerance (start loose, refine iteratively)

### Bracketing Strategy for Multiple Roots

```cpp
// Find all roots in [a, b] efficiently:

// Step 1: Coarse survey
Vector<Real> xb1_coarse, xb2_coarse;
int n_coarse = FindRootBrackets(f, a, b, 100, xb1_coarse, xb2_coarse);

// Step 2: Refine each bracket region with finer sampling
for (int i = 0; i < n_coarse; i++) {
  Vector<Real> xb1_fine, xb2_fine;
  int n_fine = FindRootBrackets(f, xb1_coarse[i], xb2_coarse[i], 100, xb1_fine, xb2_fine);
  
  // Step 3: Solve each fine bracket
  for (int j = 0; j < n_fine; j++) {
    Real root = FindRootNewton(f, xb1_fine[j], xb2_fine[j], 1e-12);
    // Store root...
  }
}

// Adaptive resolution based on function complexity!
```

### Parallel Root Finding

```cpp
// If multiple independent brackets, refine in parallel:
std::vector<Real> roots(numBrackets);

#pragma omp parallel for
for (int i = 0; i < numBrackets; i++) {
  roots[i] = RootFinding::FindRootNewton(f, xb1[i], xb2[i], 1e-10);
}
```

## Advanced Topics

### Deflation for Multiple Roots

After finding root r₁, find next root by deflation:

```cpp
class DeflatedFunction : public IRealFunction {
  const IRealFunction& f_original;
  std::vector<Real> known_roots;
public:
  DeflatedFunction(const IRealFunction& f) : f_original(f) {}
  
  void AddRoot(Real root) { known_roots.push_back(root); }
  
  Real operator()(Real x) const override {
    Real result = f_original(x);
    for (Real r : known_roots)
      result /= (x - r);  // Divide out (x - r) factor
    return result;
  }
};

// Usage:
DeflatedFunction f_deflated(original_f);
Real r1 = FindRootNewton(original_f, x1, x2, tol);
f_deflated.AddRoot(r1);

Real r2 = FindRootNewton(f_deflated, y1, y2, tol);
f_deflated.AddRoot(r2);

Real r3 = FindRootNewton(f_deflated, z1, z2, tol);
```

**Caution**: Numerical errors accumulate! Deflation works best for well-separated roots.

### Complex Roots (Future Extension)

Current implementation: **Real roots only**

For complex roots, need:
1. Complex arithmetic support
2. Modified Newton for complex domain
3. Different convergence criteria

Placeholder for future `FindRootComplex()`.

### Hybrid Methods

**Newton with Bisection Backup** (idea for FindRootSafe):
```cpp
Real FindRootHybrid(const IRealFunction& f, Real x1, Real x2, Real tol) {
  Real root = 0.5 * (x1 + x2);
  Real dx_old = std::abs(x2 - x1);
  Real dx = dx_old;
  
  for (int iter = 0; iter < maxIter; iter++) {
    Real fx = f(root);
    Real dfx = Derivation::NDer4(f, root);
    
    // Check if Newton step stays in bounds and converges faster than bisection
    Real dx_newton = fx / dfx;
    if (std::abs(dx_newton) < std::abs(dx) && 
        root - dx_newton > x1 && root - dx_newton < x2) {
      // Take Newton step
      dx = dx_newton;
      root -= dx;
    } else {
      // Fall back to bisection
      dx = 0.5 * (x2 - x1);
      root = x1 + dx;
      
      if (f(x1) * fx < 0) x2 = root;
      else                x1 = root;
    }
    
    if (std::abs(dx) < tol) return root;
  }
  throw RootFindingError("Hybrid method failed to converge");
}
```

Combines **robustness of bisection** with **speed of Newton**.

---

## Runnable Examples

| Topic | Demo File | Function |
|-------|-----------|----------|
| Find Initial Bracket | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_BracketRoot()` |
| Find All Brackets | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootBrackets()` |
| Bisection Method | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootBisection()` |
| False Position | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootFalsePosition()` |
| Secant Method | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootSecant()` |
| Newton-Raphson | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootNewton()` |
| Ridders' Method | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootRidders()` |
| Brent's Method | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootBrent()` |
| Method Comparison | `src/docs_demos/docs_demo_root_finding.cpp` | `Docs_Demo_RootFinding_Comparison()` |

*Comprehensive demos covering all 8 root finding methods with comparison.*

---

## Summary

### Method Comparison

| Feature | Bisection | FalsePos | Secant | Newton | Ridders | Brent |
|---------|-----------|----------|--------|--------|---------|-------|
| **Convergence** | Linear | ~φ | ~φ | Quadratic | Quadratic | Adaptive |
| **Derivative** | No | No | No | Yes | No | No |
| **Bracket** | Yes | Yes | No | Optional | Yes | Yes |
| **Robustness** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Speed** | Slow | Medium | Fast | Fastest | Fast | Fast |
| **Best For** | Robustness | Faster bracket | Simple | Precision | No-deriv fast | Default |

### When to Use What

**Exploration Phase** (unknown function):
1. `FindRootBrackets()` with moderate sampling (100-1000 points)
2. Visualize function if possible
3. Identify approximate root locations

**Refinement Phase** (bracket known):
1. Try `FindRootBrent()` first (recommended default!)
2. Use `FindRootNewton()` for maximum speed with good guesses
3. Fall back to `FindRootBisection()` for guaranteed convergence
4. Verify with residual check: |f(x*)| < tolerance

**Production Code**:
- **Brent** for best all-around performance (default recommendation)
- **Bisection** for guaranteed robustness when reliability critical
- **Newton** for performance-critical paths with good initial guesses
- **Ridders** when quadratic convergence needed without derivatives

### Key Takeaways

1. ✅ **Always bracket first** - guarantees root existence
2. ✅ **Use Brent as default** - best all-around performance
3. ✅ **Verify convergence** - check |f(x*)| after solving
4. ✅ **Scale functions** - avoid extreme magnitude differences
5. ✅ **Newton for speed** - when f smooth and guess good
6. ✅ **Bisection for robustness** - when convergence critical
7. ✅ **Sample adequately** - 10× expected root spacing minimum
8. ✅ **Watch for f'(x*) ≈ 0** - Newton fails, use bracketing methods
9. ✅ **Catch exceptions** - RootFindingError signals algorithm failure

### References

- **Numerical Recipes** (Press et al.): Chapter 9 - Root Finding and Nonlinear Sets of Equations
- **Numerical Analysis** (Burden & Faires): Chapter 2 - Solutions of Equations in One Variable
- **Numerical Methods** (Chapra & Canale): Chapter 5 & 6 - Roots: Bracketing and Open Methods
- **Brent (1973)**: Algorithms for Minimization without Derivatives

---

**Part of MinimalMathLibrary** - `mml/algorithms/RootFinding.h`
