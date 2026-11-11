# Numerical Derivation - Finite Difference Methods

This document covers numerical computation of derivatives using finite difference approximations.

## Table of Contents
- [Overview](#overview)
- [Order Selection Guide](#order-selection-guide)
- [API Reference](#api-reference)
- [Automatic Step Size](#automatic-step-size)
- [Legacy One-Sided Sampling](#legacy-one-sided-sampling)
- [Error Estimation](#error-estimation)
- [Function Types](#function-types)
- [Jacobians, Hessians and Gradients — one vocabulary](#jacobians-hessians-and-gradients--one-vocabulary)
- [Higher-Order Derivatives](#higher-order-derivatives)
- [Examples](#examples)
- [Performance and Accuracy](#performance-and-accuracy)
- [Common Pitfalls](#common-pitfalls)
- [See Also](#see-also)

---

## Overview

MinimalMathLibrary provides **finite-difference numerical differentiation** for functions where analytical derivatives are unavailable or impractical.

**Available Methods:** Orders 1, 2, 4, 6, and 8

| Method | Formula | Error | Function Evaluations | Accuracy |
|--------|---------|-------|---------------------|----------|
| **Order 1** (forward) | f'(x) ≈ (f(x+h) - f(x))/h | O(h) | 2 | Low |
| **Order 2** (central) | f'(x) ≈ (f(x+h) - f(x-h))/(2h) | O(h²) | 2 | Good |
| **Order 4** (central) | Complex 5-point stencil | O(h⁴) | 4-5 | **Recommended** |
| **Order 6** (central) | Complex 7-point stencil | O(h⁶) | 6-7 | High precision |
| **Order 8** (central) | Complex 9-point stencil | O(h⁸) | 8-9 | Very high precision |

**Key Features:**
- **Automatic step sizing:** Optimal `h` based on order and `Real` precision
- **Error estimation:** Optional error output parameter
- **Detailed result honesty:** default results are `NotComputed`; disabled or failed
  error estimates are `NoEstimate`, while a genuinely computed zero remains available
  through `HasEstimate(result.error)` (see [RESULT_POLICY](../RESULT_POLICY.md)).
- **Multiple function types:** Real, Scalar, Vector, Parametric curves, Tensor fields
- **Higher derivatives:** Second and third derivatives available
- **Legacy left/right wrappers:** Shifted stencils; see their actual evaluation centers and sampling reach below

**Default:** `Derive()` uses **Order 4** (best balance of accuracy and robustness)

All derivation methods are defined in:
- [`mml/core/Derivation.h`](../../mml/core/Derivation.h) - Main header
- [`mml/core/Derivation/DerivationRealFunction.h`](../../mml/core/Derivation/DerivationRealFunction.h)
- [`mml/core/Derivation/DerivationScalarFunction.h`](../../mml/core/Derivation/DerivationScalarFunction.h)
- [`mml/core/Derivation/DerivationVectorFunction.h`](../../mml/core/Derivation/DerivationVectorFunction.h)

---

## Order Selection Guide

### Decision Tree

```
Is your function smooth (infinitely differentiable)?
├─ YES (e.g., sin, exp, polynomials)
│  ├─ Need high precision? → Order 6 or 8
│  └─ General use? → Order 4 (default)
└─ NO (noisy, discrete, or limited precision data)
   ├─ Noisy data? → Order 2 (more robust)
   └─ Quick estimate? → Order 1 (fastest)
```

### Recommendations by Use Case

**General Purpose (Default):**
- **Order 4** (`Derive()` or `NDer4()`)
- Best balance of accuracy and stability
- Recommended for most applications

**High Accuracy (Smooth Functions):**
- **Order 6** (`NDer6()`) or **Order 8** (`NDer8()`)
- Use when function is very smooth
- Requires low numerical noise
- Example: analytic test functions, well-conditioned physics

**Noisy Data:**
- **Order 2** (`NDer2()`)
- More robust to small perturbations
- Better for experimental/measured data

**Quick Estimation:**
- **Order 1** (`NDer1()`)
- Cheapest (fewest function evaluations)
- Use only for rough estimates

### Accuracy vs Stability Trade-off

**Higher order = Higher accuracy BUT:**
- More sensitive to numerical noise
- Requires smaller step size `h`
- More function evaluations

**Optimal Order Selection:**
```cpp
// Smooth analytical function → Order 6-8
Real exact_der = Derivation::NDer8(sin_func, x);

// General case → Order 4 (default)
Real balanced_der = Derivation::Derive(some_func, x, &err);

// Noisy data → Order 2
Real robust_der = Derivation::NDer2(noisy_data_func, x);
```

---

## API Reference

### Richardson refinement contracts

`NDerRichardson`, `NSecDerRichardson` and `NThirdDerRichardson` require a finite
evaluation point, a positive finite initial step, a positive tableau depth
`max_iter`, and a finite shrinking factor `con > 1` whose square is finite.
Invalid parameters throw `NumericInputError` before allocating the tableau or
evaluating the function. Zero step is not an automatic-step request in these APIs.

Each initial or refined stencil must have finite, strictly ordered, distinct
coordinates around the requested point; the third derivative also requires
distinct outer samples at twice the step. Collapsed or overflowing coordinates,
step underflow, nonfinite callback values and nonfinite extrapolation arithmetic
throw `NumericalMethodError`. A failed refinement does not silently return an
earlier estimate or write a success-shaped error estimate. Callback exceptions
propagate unchanged. The adaptive aliases share these contracts.

The error output is a heuristic difference between neighboring Richardson
tableau estimates, not a certified error bound. A one-level tableau is allowed
and returns the raw central-difference estimate, with
`std::numeric_limits<Real>::max()` meaning that no extrapolation error estimate
is available. On failure, the caller's error output is left unchanged.

Rank-3 and rank-4 tensor `NDer1Partial`, `NDer2Partial` and `NDer4Partial`
default-step overloads use the same formulas and scaled steps as their
explicit-step counterparts, including optional error estimation.

### Real Functions (f: ℝ → ℝ)

```cpp
namespace Derivation {

// ========== First Derivative ==========

// Order 1 (forward difference)
Real NDer1(const IRealFunction& f, Real x);
Real NDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr);

// Order 2 (central difference)
Real NDer2(const IRealFunction& f, Real x);
Real NDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr);

// Order 4 (5-point stencil) - DEFAULT
Real NDer4(const IRealFunction& f, Real x);
Real NDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr);
Real Derive(const IRealFunction& f, Real x, Real* error = nullptr); // Calls NDer4

// Order 6 (7-point stencil)
Real NDer6(const IRealFunction& f, Real x);
Real NDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr);

// Order 8 (9-point stencil)
Real NDer8(const IRealFunction& f, Real x);
Real NDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr);

// ========== One-Sided Derivatives ==========

// Legacy left/right wrappers: higher orders evaluate at shifted centers, not at x
Real NDer1Left(const IRealFunction& f, Real x, Real* error = nullptr);
Real NDer1Right(const IRealFunction& f, Real x, Real* error = nullptr);
Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr);
Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr);
// ... (available for all orders)

// ========== Second Derivative ==========

Real NSecDer1(const IRealFunction& f, Real x);
Real NSecDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr);
Real NSecDer2(const IRealFunction& f, Real x);
Real NSecDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr);
// ... (NSecDer4, NSecDer6, NSecDer8)

// ========== Third Derivative ==========

Real NThirdDer1(const IRealFunction& f, Real x);
Real NThirdDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr);
Real NThirdDer2(const IRealFunction& f, Real x);
// ... (NThirdDer4, NThirdDer6, NThirdDer8)

}
```

### Scalar Functions (f: ℝⁿ → ℝ)

```cpp
// Partial derivative ∂f/∂xᵢ at point
template<int N>
Real NDer4Partial(
    const IScalarFunction<N>& f,
    int deriv_index,              // Index of variable (0 to N-1)
    const VectorN<Real, N>& point,
    Real* error = nullptr
);

// Gradient: all partial derivatives ∇f = [∂f/∂x₁, ..., ∂f/∂xₙ]
template<int N>
VectorN<Real, N> NDer4PartialByAll(
    const IScalarFunction<N>& f,
    const VectorN<Real, N>& point,
    VectorN<Real, N>* error = nullptr
);

// Second partial derivative ∂²f/(∂xᵢ∂xⱼ)
template<int N>
Real NSecDer4Partial(
    const IScalarFunction<N>& f,
    int der_ind1, int der_ind2,   // Indices of variables
    const VectorN<Real, N>& point,
    Real* error = nullptr
);

// Available for all orders: NDer1Partial, NDer2Partial, NDer6Partial, NDer8Partial
```

### Vector Functions (f: ℝⁿ → ℝⁿ)

```cpp
// Single component partial: ∂fᵢ/∂xⱼ
template<int N>
Real NDer4Partial(
    const IVectorFunction<N>& f,
    int func_index,               // Component of output (0 to N-1)
    int deriv_index,              // Variable to differentiate (0 to N-1)
    const VectorN<Real, N>& point,
    Real* error = nullptr
);

// One row of Jacobian: [∂fᵢ/∂x₁, ..., ∂fᵢ/∂xₙ]
template<int N>
VectorN<Real, N> NDer4PartialByAll(
    const IVectorFunction<N>& f,
    int func_index,
    const VectorN<Real, N>& point,
    VectorN<Real, N>* error = nullptr
);

// Full Jacobian matrix: J[i,j] = ∂fᵢ/∂xⱼ
template<int N>
MatrixNM<Real, N, N> NDer4PartialAllByAll(
    const IVectorFunction<N>& f,
    const VectorN<Real, N>& point,
    MatrixNM<Real, N, N>* error = nullptr
);
```

### Parametric Curves (f: ℝ → ℝⁿ)

```cpp
// First derivative: tangent vector
template<int N>
VectorN<Real, N> NDer4(
    const IParametricCurve<N>& curve,
    Real t,
    Real* error = nullptr
);

// Second derivative: acceleration/curvature
template<int N>
VectorN<Real, N> NSecDer4(
    const IParametricCurve<N>& curve,
    Real t,
    Real* error = nullptr
);

// Third derivative
template<int N>
VectorN<Real, N> NThirdDer4(
    const IParametricCurve<N>& curve,
    Real t,
    Real* error = nullptr
);
```

### Tensor Fields (T: ℝⁿ → ℝⁿˣⁿ or higher)

```cpp
// Partial derivative of tensor component
template<int N>
Real NDer4Partial(
    const ITensorField2<N>& field,
    int i, int j,                 // Tensor indices
    int deriv_index,              // Spatial derivative index
    const VectorN<Real, N>& point,
    Real* error = nullptr
);

// Available for ITensorField3<N>, ITensorField4<N>, ITensorField5<N>
```

---

## Automatic Step Size

An explicitly supplied derivative step must be **finite and strictly positive**. Zero,
negative, NaN and infinite steps throw `NumericInputError` before any callback, with or
without an error pointer, across real/complex, scalar/vector partials, tensor fields,
curves and surfaces. Aggregate gradient/Jacobian-element wrappers inherit this rule;
error outputs are not modified on invalid-step rejection. Detailed real/complex APIs
convert it to `InvalidInput` only when `ConvertToStatus` is selected.

Zero means automatic **only for APIs that explicitly promise it**: `DerivativeConfig::step`,
the numerical `CalcJacobian*`/`CalcHessian` entry points, and `NDer2PartialDyn`/
`NDer4PartialDyn`. It is not an automatic-step sentinel for an explicit `NDer*` step.

When **h is omitted**, the step is selected based on:
1. **Derivative order and stencil accuracy order** (both affect the step scale)
2. **Machine epsilon** (`std::numeric_limits<Real>::epsilon()`)
3. **Balance between truncation and roundoff errors**

### Automatic Step Formulas

The first-derivative base steps in `DerivationBase.h`, with epsilon supplied by `Constants::Eps`, are:

| Accuracy order | Base step |
|----------------|-----------|
| 1 | 2√ε |
| 2 | (3ε)^(1/3) |
| 4 | (11.25ε)^(1/5) |
| 6 | (ε/168)^(1/7) |
| 8 | (551.25ε)^(1/9) |

Ordinary default overloads multiply the relevant base step by `max(1, |x|)`, or the
corresponding input coordinate magnitude for partials (`|z|` for complex arguments).
An explicit positive step is used unchanged. The Jacobian-specific conventions below
are intentionally distinct.

Surface second partials `NDer2_uu`, `NDer2_uw` and `NDer2_ww` instead use
`Precision::FourthRootEps`, approximately epsilon^(1/4). Their O(h²) truncation
balances O(epsilon/h²) roundoff, not the O(epsilon/h) term of a first derivative.
Pure partials retain coordinate-relative scaling; the mixed default retains its
common step scaled by `max(1, |u|, |w|)`. The one-step and independent-u/w-step
explicit overloads and first-partial defaults are unchanged.
This is a smooth, well-scaled-function heuristic, not a universal accuracy
guarantee: translation cancellation, parameter conditioning and higher derivatives
can still require explicit steps.

These centered surface routines do not clip or switch to one-sided stencils at
chart boundaries. All requested samples must be valid: uu/ww error estimation
also evaluates at ±2h, while mixed error estimation stays inside its corner
stencil. A larger automatic step may require more room inside a bounded chart.

Real-function and curve `NSecDer2`/`NSecDer4`, including the detailed real-function
APIs, and scalar `NSecDer2Partial`/`NSecDer4Partial` use derivative-order-aware bases:

| Second-derivative accuracy order | Base step | Pure samples without / with error |
|---|---|---|
| 2 | `Precision::FourthRootEps`, approximately epsilon^(1/4) | 0, ±h / also ±2h |
| 4 | `Precision::SixthRootEps`, epsilon^(1/6) | 0, ±h, ±2h / also ±3h |

The sixth-root constant is computed with typed square/cube roots of the actual
type's epsilon, including platforms where `long double` differs from `double`.
Real and curve defaults retain `max(1, |x|)` scaling. A scalar pure partial scales
by its coordinate, while a mixed partial retains a common h scaled by
`max(1, |x_i|, |x_j|)`. The fourth-order mixed stencil reaches the ±2h corners;
error estimation also uses half-h corners.

Default second-derivative aliases and numerical `CalcHessian` inherit these bases.
Positive explicit steps and nonzero `DerivativeConfig::step` remain absolute;
an explicit step argument takes priority over the config field. First derivatives,
third derivatives, Richardson/adaptive starting steps, analytic/AD dispatch and
error/result policies are unchanged. No centered stencil is clipped or changed to
one-sided at a boundary. Larger automatic steps may expose a callback-domain error
that a smaller step avoided; choose an explicit step when interior room is limited.
These defaults and error estimates are smooth-function heuristics, not certified
accuracy guarantees.

Runnable example: `MML_DocsApp --second-derivative-steps`.

```cpp
RealFunction exponential([](Real x) { return std::exp(x); });
const auto result = Derivation::NSecDer4Detailed(exponential, Real(0.3));
std::cout << "Second derivative: " << result.value << '\n';
std::cout << "Step used: " << result.step_used << '\n';
```

`ComplexStep` and `ComplexStep2` also require finite positive steps and finite real
evaluation points, and reject nonfinite callback values or final quotients with
`NumericalMethodError`. The first method retains support for extremely small positive
steps, including representable subnormals. The second method requires a positive finite
representable `h²`, checked before callbacks; it involves real subtraction and remains
subject to cancellation. Both propagate callback exceptions unchanged.

### When to Use Custom h

**Use automatic h (omit parameter) when:**
- Function is well-scaled (values ≈ 1)
- No prior knowledge of function behavior
- General-purpose computation

**Provide explicit h when:**
- Function has extreme scaling (very large/small values)
- You know local curvature is very high/low
- Specific accuracy requirements
- Testing/debugging derivatives

**Guidelines for custom h:**
```cpp
// Too large h → truncation error dominates
// Too small h → roundoff error dominates
// Optimal: balance both errors

// Rough estimate for order n:
h_optimal ≈ ε^(1/(n+1)) * characteristic_scale
```

**Example: Function with large values**
```cpp
RealFunction huge([](Real x) { return 1e10 * sin(x); });

// Automatic h might be too small (roundoff errors)
Real der_auto = Derivation::NDer4(huge, 1.0);

// Manual h scaled to function magnitude
Real scale = 1e10;
Real h_custom = NDer4_h * pow(scale, 1.0/5.0);
Real der_manual = Derivation::NDer4(huge, 1.0, h_custom);
```

---

## Legacy One-Sided Sampling

The `NDer*Left`/`NDer*Right` wrappers retain their historical formulas. They are
**not a family of endpoint-anchored high-order one-sided formulas**. Orders 2, 4, 6
and 8 call a central stencil at a shifted center. Automatic overloads choose the
step using the original `x`, then shift the center using that same step.

Offsets below are multiples of positive `h` relative to the original `x`. A range
describes the sampling extent, not a promise to evaluate every interior grid point.

| Accuracy order | Underlying center, Left / Right | Value-sample extent, Left / Right | Extent with error, Left / Right |
|----------------|--------------------------------|----------------------------------|--------------------------------|
| 1 | x-h / x (forward stencil) | [-1, 0] / [0, 1] | [-2, 0] / [-1, 1] |
| 2 | x-3h / x+3h | [-4, -2] / [2, 4] | [-5, -1] / [1, 5] |
| 4 | x-4h / x+4h | [-6, -2] / [2, 6] | [-7, -1] / [1, 7] |
| 6 | x-5h / x+5h | [-8, -2] / [2, 8] | [-9, -1] / [1, 9] |
| 8 | x-6h / x+6h | [-10, -2] / [2, 10] | [-11, -1] / [1, 11] |

For example, `NDer4Left` applied to x² at x=1 with h=0.01 returns approximately
**1.92**, the derivative at x=0.96, not the endpoint derivative 2. Its optional error
estimates the underlying stencil error; it does **not** include the bias from moving
the requested evaluation point. The first-order wrappers are ordinary forward/backward
quotients and retain their first-order endpoint bias.

Requesting error estimates adds samples. In particular, `NDer1Right` at x=0 samples
0 and h without error, but also **-h** with error, so a function restricted to
nonnegative inputs can fail only when an error pointer is supplied. Higher-order
shifted wrappers stay on their named side but reach farther with error. Check the
entire required interval, not just the named point. Callback exceptions propagate
unchanged; direct error outputs are written only after a successful calculation.

## Error Estimation

### Error Output Parameter

Most derivation functions accept an **optional error pointer**:

```cpp
Real error;
Real derivative = Derivation::NDer4(func, x, &error);

std::cout << "f'(" << x << ") = " << derivative 
          << " ± " << error << std::endl;
```

### Error Composition

The reported heuristic is **estimated truncation error + a roundoff allowance**,
not a certified error bound. The value returned at the requested step is not
replaced by a finer or extrapolated value just because an error output is requested.

**Truncation error** (from finite h):
- Depends on higher derivatives and h
- Order n method: truncation ∝ h^n

**Roundoff error** (from finite precision):
- From function evaluation errors
- For an m-th derivative, scales roughly as ε·|f|/h^m (grows as h decreases!)

**Optimal h minimizes total error:**
For an m-th derivative with accuracy order p, the model is
C₁·h^p + C₂·ε/h^m, with an optimum proportional to ε^(1/(m+p)).

### Leading Truncation Terms

These smooth-function Taylor terms describe magnitude, not rigorous bounds.
Higher derivatives are approximated by centered finite-difference residuals
with the corresponding powers of h; vector-valued results use Euclidean norms.

| Returned derivative | Accuracy order | Leading truncation magnitude |
|---|---:|---|
| First derivative | 2 | \|f'''(x)\| h² / 6 |
| First derivative | 4 | \|f⁽⁵⁾(x)\| h⁴ / 30 |
| Pure second derivative | 2 | \|f⁽⁴⁾(x)\| h² / 12 |
| Pure second derivative | 4 | \|f⁽⁶⁾(x)\| h⁴ / 90 |
| Third derivative | 2 | \|f⁽⁵⁾(x)\| h² / 4 |
| Third derivative | 4 | 7 \|f⁽⁷⁾(x)\| h⁴ / 120 |

The pure second-derivative formulas are shared by real functions, scalar
partials, curves and surface uu/ww derivatives; third-derivative formulas are
shared by real functions and curves.

### Mixed Partials and Evaluation Cost

Scalar mixed `NSecDer2Partial` / `NSecDer4Partial` and surface `NDer2_uw`
include truncation by comparing the requested-step result Dₕ with Dₕ/₂.
The truncation estimate is 4/3 times their difference for second-order methods
and 16/15 times their difference for fourth-order methods. Roundoff allowances
from both calculations are included; fourth-order scalar estimates include the
outer 2h samples as well. Surface u/w steps can differ; both are halved together.
The derivative returned remains Dₕ.

| Routine | Without error output | With error output |
|---|---:|---:|
| Real/curve pure second, order 2 | 3 | 5 |
| Real/curve pure second, order 4 | 5 | 7 |
| Real/curve third, order 2 | 4 | 6 |
| Real/curve third, order 4 | 6 | 8 |
| Scalar pure second, order 2 / 4 | 3 / 5 | 5 / 7 |
| Scalar mixed second, order 2 / 4 | 4 / 8 | 8 / 12 |
| Surface uu/ww | 3 | 5 |
| Surface uw | 4 | 8 |

Mixed error estimation adds only samples inside the existing stencil region;
fourth-order scalar partials reuse the h-grid in the refined calculation.
Surface uu/ww error estimation extends the stencil from ±h to ±2h.
As with other centered stencils, callers must provide a function valid at all
sampled coordinates; these routines do not clamp to a declared parameter domain.
The affected higher-derivative entry points require positive finite explicit
steps. If mixed error refinement cannot form distinct finite half-step
coordinates, it throws `NumericalMethodError` instead of reporting zero
truncation. Callback exceptions propagate. Disable the error output to avoid
the additional samples, not to obtain a domain-safe one-sided stencil.

### Interpreting Error Estimates

```cpp
Real err;
Real der = Derivation::Derive(f, x, &err);

if (err < 1e-10) {
    // Excellent accuracy
} else if (err < 1e-6) {
    // Good accuracy for most applications
} else if (err < 1e-3) {
    // Moderate accuracy, may need refinement
} else {
    // Poor accuracy:
    // - Try different order
    // - Adjust h manually
    // - Function may be poorly behaved
}
```

---

## Function Types

All derivation methods work with **function interfaces** (see [Functions.md](Functions.md)):

### 1. IRealFunction (f: ℝ → ℝ)

```cpp
#include <mml/core/Functions.h>
#include <mml/core/Derivation.h>

RealFunction f([](Real x) { return sin(x) * exp(-x); });

Real x = 1.0;
Real df = Derivation::Derive(f, x);          // f'(1.0)
Real d2f = Derivation::NSecDer4(f, x);       // f''(1.0)
Real d3f = Derivation::NThirdDer4(f, x);     // f'''(1.0)
```

### 2. IScalarFunction<N> (f: ℝⁿ → ℝ)

**Partial derivatives:**

```cpp
ScalarFunction<3> f([](const VectorN<Real, 3>& x) {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];  // ||x||²
});

VectorN<Real, 3> point{1.0, 2.0, 3.0};

// Single partial: ∂f/∂x₁ at point
Real df_dx1 = Derivation::NDer4Partial(f, 0, point);  // = 2*x[0] = 2.0

// Gradient: ∇f = [∂f/∂x₁, ∂f/∂x₂, ∂f/∂x₃]
VectorN<Real, 3> grad = Derivation::NDer4PartialByAll(f, point);
// grad = [2.0, 4.0, 6.0]

// Second partial: ∂²f/(∂x₁∂x₂)
Real d2f = Derivation::NSecDer4Partial(f, 0, 1, point);  // = 0 (mixed partials)
```

### 3. IVectorFunction<N> (f: ℝⁿ → ℝⁿ)

**Jacobian matrix:**

```cpp
VectorFunction<3> f([](const VectorN<Real, 3>& x) {
    return VectorN<Real, 3>{
        x[0] * x[1],           // f₁ = x·y
        x[1] * x[2],           // f₂ = y·z
        x[0] * x[2]            // f₃ = x·z
    };
});

VectorN<Real, 3> point{1.0, 2.0, 3.0};

// Single Jacobian element: J[1,2] = ∂f₂/∂x₃
Real J_12 = Derivation::NDer4Partial(f, 1, 2, point);  // = y = 2.0

// One row: [∂f₂/∂x₁, ∂f₂/∂x₂, ∂f₂/∂x₃]
VectorN<Real, 3> row2 = Derivation::NDer4PartialByAll(f, 1, point);
// row2 = [0, z, y] = [0, 3, 2]

// Full Jacobian: J[i,j] = ∂fᵢ/∂xⱼ
MatrixNM<Real, 3, 3> J = Derivation::NDer4PartialAllByAll(f, point);
/*
J = [ y  x  0 ]   [ 2  1  0 ]
    [ 0  z  y ] = [ 0  3  2 ]
    [ z  0  x ]   [ 3  0  1 ]
*/
```

For full fixed-size or dynamic vector-function Jacobians, prefer
`CalcJacobian` (fixed size via `IVectorFunction`, dynamic size via `std::function`). They perturb one input coordinate at a time
and fill the complete output column from each vector evaluation, requiring four
vector-function calls per input dimension:

```cpp
int evaluations = 0;
Matrix<Real> J = Derivation::CalcJacobian(function, point, 0.0, &evaluations);
// evaluations == 4 * point.size()
```

With `h == 0`, each column uses
$h_j=\epsilon^{1/5}\max(1,|x_j|)$. A nonzero `h` remains an explicit absolute
step for every coordinate. Square and non-square outputs are supported, output
dimensions are validated, non-finite values are rejected, and the optional
evaluation pointer accumulates attempted vector-function calls.

### 4. IParametricCurve<N> (r: ℝ → ℝⁿ)

**Curve derivatives:**

```cpp
ParametricCurve<3> helix([](Real t) {
    return VectorN<Real, 3>{
        cos(t),        // x = cos(t)
        sin(t),        // y = sin(t)
        t              // z = t
    };
});

Real t = M_PI / 4;

// Tangent vector: dr/dt
VectorN<Real, 3> tangent = Derivation::NDer4(helix, t);
// tangent ≈ [-sin(t), cos(t), 1] = [-0.707, 0.707, 1]

// Acceleration: d²r/dt²
VectorN<Real, 3> accel = Derivation::NSecDer4(helix, t);
// accel ≈ [-cos(t), -sin(t), 0] = [-0.707, -0.707, 0]

// Speed: ||dr/dt||
Real speed = tangent.NormL2();  // ≈ 1.414
```

### 5. ITensorField2<N> (T: ℝⁿ → ℝⁿˣⁿ)

**Tensor field derivatives:**

```cpp
// Stress tensor field σ(x)
TensorField2<3> stress_field([](const VectorN<Real, 3>& x) {
    // Return 3×3 stress tensor at position x
    MatrixNM<Real, 3, 3> sigma;
    sigma[0][0] = x[0] * x[0];  // σₓₓ
    sigma[1][1] = x[1] * x[1];  // σᵧᵧ
    // ...
    return sigma;
});

VectorN<Real, 3> pos{1.0, 1.0, 1.0};

// Derivative of tensor component: ∂σₓₓ/∂x
Real dSxx_dx = Derivation::NDer4Partial(stress_field, 0, 0, 0, pos);
// = 2*x = 2.0
```

---

## Jacobians, Hessians and Gradients — one vocabulary

MML computes Jacobians three ways — finite differences (`MML::Derivation`), forward-mode
automatic differentiation with `Dual<T>` numbers and reverse-mode AD with a `Tape` (both in
`MML::AD`). Since 2.1 they share one naming scheme ([CONTRACT_POLICY §11](../CONTRACT_POLICY.md#11-naming-convention)):
`Calc<Object>[<Method>]`, where the method suffix names the technology and the *argument types*
select fixed vs. dynamic size (no `Dyn` suffix). The 2.0 spellings (`calcJacobian`,
`calcJacobianDyn`, `jacobianForwardAD`, `AD::gradient`, …) remain as `[[deprecated]]` aliases.

| Entry point | Input | Result | Method | Default step | Order / cost |
|---|---|---|---|---|---|
| `Derivation::CalcJacobian(f, pos [, h, &evals])` | `IVectorFunction<N>` / `IVectorFunctionNM<N,M>` + `VectorN` | `MatrixNM` | central differences | `ε^{1/5}·max(1,|x_j|)` per column | 4th order, 4·N vector evaluations |
| `Derivation::CalcJacobian(f, pos [, outputDim] [, h, &evals])` | `std::function<Vector(const Vector&)>` + `Vector` | `Matrix` | central differences | same | 4th order, 4·n; validates output dimension and finiteness |
| `Derivation::CalcJacobianInPlace(f, pos, J [, h, &evals])`, `CalcJacobianColumns(f, pos, J, m, h, &evals)` | as above, pre-allocated `J` | in place | same | same | same — for Newton loops |
| `Derivation::CalcJacobianOrder2(f, pos [, h])` | `std::function` + `Vector` | `Matrix` | central differences | `ε^{1/3}` | 2nd order, 2·n² vector evaluations (per-entry implementation) |
| `Derivation::CalcJacobianODE(derivs, t, x [, J] [, h])` | `derivs(t, x, dxdt)` | `Matrix` | central differences | `√ε·max(1,|x_j|)` | 2nd order, 2·n — what stiff steppers use |
| `Derivation::CalcJacobianMap(map, x, J [, h])` | `x_{n+1} = f(x_n)` | in place | central differences | `√ε·max(1,|x_j|)` | 2nd order, 2·n |
| `Derivation::CalcHessian(f, pos [, H] [, h])` | `IScalarFunction<N>` + `VectorN` | `MatrixNM` | `NSecDer4Partial` per entry | per-entry automatic step | 4th order, N(N+1)/2 mixed partials |
| `AD::CalcJacobianForwardAD(f, point)` | callable on `std::vector<Dual<T>>` | `Matrix<T>` | forward AD | — (exact to rounding) | n+1 evaluations |
| `AD::CalcJacobianReverseAD(f, point)` | callable on `std::vector<ADVar>` | `Matrix<Real>` | reverse AD | — | m+1 evaluations + m backward sweeps |
| `AD::CalcHessianForwardAD(f, point)` | callable on `std::vector<Dual<Dual<T>>>` | `Matrix<T>` | nested forward AD | — | n(n+1)/2 evaluations |
| `AD::CalcGradientForwardAD(f, x)`, `AD::CalcGradientReverseAD(f, x)` | scalar callable | `std::vector` | forward / reverse AD | — | n evaluations / 1 evaluation + 1 sweep |
| `AD::CalcDerivativeForwardAD(f, x)` | scalar callable on `Dual<T>` | `{value, derivative}` | forward AD | — | 1 evaluation |

The finite-difference entries report **no error estimate** — they are raw matrices, and the step
column above is the only accuracy knob. Zero selects the default. An explicit positive step is
absolute for `CalcJacobian`, `CalcJacobianOrder2`, the dynamic component helpers and `CalcHessian`.
For legacy `CalcJacobianODE`/`CalcJacobianMap`, it is instead a **relative factor**:
the actual step is `h·max(1,|x_j|)`. The dynamic component helpers select unscaled
`ε^{1/3}`/`ε^{1/5}` when zero is passed; large coordinates can therefore require an explicit step.

Numerical Jacobians cancel common sample offsets before applying binary-scaled difference
quotients, avoiding unnecessary overflow of weighted samples or denominators. All public
Jacobian paths reject nonfinite input coordinates, nonfinite/negative steps, collapsed or
overflowing stencil coordinates, invalid callback output shapes/nonfinite components and
nonfinite final entries. Invalid input arguments throw `ArgumentError`; unrepresentable stencils
or numerical values throw `NumericalMethodError`; shape mismatches throw `VectorDimensionError`.
Each column's coordinates are checked before its callbacks, and a rejected stencil makes no
callbacks for that column. Callback exceptions and accumulated attempted-call counters are
preserved. In-place outputs need not be complete after an exception. These checks are safety
contracts, not certified accuracy bounds.

Extreme-scale Jacobian tests must use finite, nonzero magnitudes representable in
the active `Real` type; double-range constants can instead become infinity or zero
in float. The regressions separately verify overflowing raw differences and
denominators, representable subnormal steps, cancellation of equal huge outer
samples, and rejection of nonfinite final entries after finite callback samples.
These preconditions distinguish numerical range handling from earlier input or
callback validation failures; they do not relax the derivative accuracy bounds.

Per-element error estimates are available through
`NDer4Partial(..., &err)` and friends. The AD entries are exact up to rounding and do not take a
step at all. For a user-supplied Jacobian on an ODE system see `IODESystemWithJacobian`.

```cpp
// Numerical: fixed size through the interface, dynamic size through std::function
VectorFunction<2> polar([](const VectorN<Real, 2>& q) {
	return VectorN<Real, 2>{ q[0] * std::cos(q[1]), q[0] * std::sin(q[1]) };
});
VectorN<Real, 2> q0{ 2.0, 0.5 };
MatrixNM<Real, 2, 2> Jfix = Derivation::CalcJacobian(polar, q0);

int evals = 0;
Matrix<Real> Jdyn = Derivation::CalcJacobian(
	[](const Vector<Real>& q) { return Vector<Real>{ q[0] * std::cos(q[1]), q[0] * std::sin(q[1]) }; },
	Vector<Real>{ 2.0, 0.5 }, 0.0, &evals);              // evals == 8

// Forward and reverse AD on the same map
auto polarAD = [](const auto& q) {
	using std::cos; using std::sin;
	return std::vector<std::decay_t<decltype(q[0])>>{ q[0] * cos(q[1]), q[0] * sin(q[1]) };
};
Matrix<Real> Jfwd = AD::CalcJacobianForwardAD(polarAD, std::vector<Real>{ 2.0, 0.5 });
Matrix<Real> Jrev = AD::CalcJacobianReverseAD(polarAD, std::vector<Real>{ 2.0, 0.5 });
// Jfix, Jdyn, Jfwd, Jrev agree to ~1e-10 (finite differences) / rounding (AD)
```

---

## Higher-Order Derivatives

For scalar-generic expressions used by existing geometry and field algorithms,
see [Automatic differentiation for geometry and fields](Automatic_differentiation.md).
The opt-in adapters retain the ordinary Real-valued interfaces while supplying
first/second partials, third curve derivatives, and metric jets. Nested elementary
functions support higher-order forward propagation; AD does not supply a
roundoff bound or remove nonsmooth points and chart singularities.

### Second Derivatives

**Real functions (f: ℝ → ℝ):**
```cpp
RealFunction f([](Real x) { return x*x*x; });  // f(x) = x³

Real x = 2.0;
Real f_x = Derivation::Derive(f, x);          // f'(2) = 12
Real f_xx = Derivation::NSecDer4(f, x);       // f''(2) = 12
```

**Scalar functions (Hessian matrix):**
```cpp
ScalarFunction<2> f([](const VectorN<Real, 2>& x) {
    return x[0]*x[0] + x[0]*x[1] + x[1]*x[1];  // x² + xy + y²
});

VectorN<Real, 2> point{1.0, 1.0};

// Second partials (Hessian elements)
Real fxx = Derivation::NSecDer4Partial(f, 0, 0, point);  // ∂²f/∂x² = 2
Real fxy = Derivation::NSecDer4Partial(f, 0, 1, point);  // ∂²f/∂x∂y = 1
Real fyy = Derivation::NSecDer4Partial(f, 1, 1, point);  // ∂²f/∂y² = 2

// Hessian matrix:
// H = [ fxx  fxy ] = [ 2  1 ]
//     [ fxy  fyy ]   [ 1  2 ]
```

### Third Derivatives

```cpp
RealFunction f([](Real x) { return x*x*x*x; });  // f(x) = x⁴

Real x = 2.0;
Real f_x = Derivation::Derive(f, x);           // f'(2) = 32
Real f_xx = Derivation::NSecDer4(f, x);        // f''(2) = 96
Real f_xxx = Derivation::NThirdDer4(f, x);     // f'''(2) = 96
```

### Curvature from Derivatives

**Curvature of plane curve y = f(x):**

```cpp
RealFunction f([](Real x) { return sin(x); });

Real x = M_PI / 4;
Real f_x = Derivation::Derive(f, x);
Real f_xx = Derivation::NSecDer4(f, x);

// Curvature κ = |f''| / (1 + f'²)^(3/2)
Real kappa = abs(f_xx) / pow(1 + f_x*f_x, 1.5);
```

**Curvature of parametric curve:**

```cpp
ParametricCurve<3> curve([](Real t) { /* ... */ });

Real t = 1.0;
VectorN<Real, 3> r_t = Derivation::NDer4(curve, t);      // dr/dt
VectorN<Real, 3> r_tt = Derivation::NSecDer4(curve, t);  // d²r/dt²

// κ = ||r' × r''|| / ||r'||³
VectorN<Real, 3> cross = r_t.Cross(r_tt);
Real kappa = cross.NormL2() / pow(r_t.NormL2(), 3);
```

---

## Examples

### Example 1: Simple Function Derivative

```cpp
#include <mml/core/Functions.h>
#include <mml/core/Derivation.h>

RealFunction f([](Real x) {
    return sin(x) * (1.0 + 0.5 * x * x);
});

Real x = 0.5;

// Order 1 (forward difference) - least accurate
Real der1 = Derivation::NDer1(f, x);

// Order 2 (central difference) - good balance
Real der2 = Derivation::NDer2(f, x);

// Order 4 (default) - recommended
Real der4 = Derivation::Derive(f, x);

// Order 6 - high precision
Real err6;
Real der6 = Derivation::NDer6(f, x, &err6);
std::cout << "f'(0.5) = " << der6 << " ± " << err6 << std::endl;

// Order 8 - very high precision
Real der8 = Derivation::NDer8(f, x);

std::cout << "Comparison of orders:\n"
          << "Order 1: " << der1 << "\n"
          << "Order 2: " << der2 << "\n"
          << "Order 4: " << der4 << "\n"
          << "Order 6: " << der6 << "\n"
          << "Order 8: " << der8 << std::endl;
```

### Example 2: Gradient of Scalar Field

```cpp
// Gravitational potential: φ = -GM/r
ScalarFunction<3> potential([](const VectorN<Real, 3>& pos) {
    Real r = pos.NormL2();
    const Real GM = 1.0;  // G * Mass
    return -GM / r;
});

VectorN<Real, 3> position{3.0, 4.0, 0.0};  // r = 5

// Gradient: gravitational field g = -∇φ
VectorN<Real, 3> grad_phi = Derivation::NDer4PartialByAll(potential, position);

VectorN<Real, 3> g_field = -1.0 * grad_phi;

std::cout << "Position: " << position << std::endl;
std::cout << "∇φ: " << grad_phi << std::endl;
std::cout << "Gravitational field g = -∇φ: " << g_field << std::endl;

// Expected: g points toward origin with magnitude GM/r² = 1/25 = 0.04
// Direction: g ∝ -r̂ = [-0.6, -0.8, 0]
// Magnitude: ||g|| = 0.04
```

### Example 3: Jacobian Matrix for Coordinate Transformation

```cpp
// Polar to Cartesian transformation: (r,θ) → (x,y)
VectorFunction<2> polar_to_cart([](const VectorN<Real, 2>& polar) {
    Real r = polar[0];
    Real theta = polar[1];
    return VectorN<Real, 2>{
        r * cos(theta),  // x
        r * sin(theta)   // y
    };
});

VectorN<Real, 2> polar{2.0, M_PI/4};  // r=2, θ=45°

// Jacobian: J[i,j] = ∂xᵢ/∂polarⱼ
MatrixNM<Real, 2, 2> J = Derivation::NDer4PartialAllByAll(polar_to_cart, polar);

std::cout << "Jacobian at (r=" << polar[0] << ", θ=" << polar[1] << "):\n"
          << J << std::endl;

/*
Expected Jacobian (analytical):
J = [ ∂x/∂r  ∂x/∂θ ]   [ cos(θ)  -r*sin(θ) ]   [ 0.707  -1.414 ]
    [ ∂y/∂r  ∂y/∂θ ] = [ sin(θ)   r*cos(θ) ] = [ 0.707   1.414 ]
*/

// Determinant = r (area scaling factor)
Real det_J = J.Det();
std::cout << "det(J) = " << det_J << " (should equal r = 2)" << std::endl;
```

### Example 4: Curvature and Torsion of 3D Curve

```cpp
// Helical curve: r(t) = (a·cos(t), a·sin(t), b·t)
Real a = 2.0, b = 1.0;
ParametricCurve<3> helix([a, b](Real t) {
    return VectorN<Real, 3>{
        a * cos(t),
        a * sin(t),
        b * t
    };
});

Real t = M_PI / 2;

// First derivative: tangent
VectorN<Real, 3> r_t = Derivation::NDer4(helix, t);

// Second derivative
VectorN<Real, 3> r_tt = Derivation::NSecDer4(helix, t);

// Third derivative (for torsion)
VectorN<Real, 3> r_ttt = Derivation::NThirdDer4(helix, t);

// Curvature: κ = ||r' × r''|| / ||r'||³
VectorN<Real, 3> cross1 = r_t.Cross(r_tt);
Real speed_cubed = pow(r_t.NormL2(), 3);
Real curvature = cross1.NormL2() / speed_cubed;

// Torsion: τ = (r' × r'') · r''' / ||r' × r''||²
Real torsion = cross1.DotProduct(r_ttt) / cross1.NormSquared();

std::cout << "Helix parameters: a=" << a << ", b=" << b << std::endl;
std::cout << "At t=" << t << ":\n";
std::cout << "Curvature κ = " << curvature << std::endl;
std::cout << "Torsion τ = " << torsion << std::endl;

// Analytical values for helix:
Real kappa_exact = a / (a*a + b*b);
Real tau_exact = b / (a*a + b*b);
std::cout << "Expected: κ = " << kappa_exact << ", τ = " << tau_exact << std::endl;
```

### Example 5: Laplacian of Scalar Field

```cpp
// Scalar field: f(x,y,z) = x² + y² - 2z²
ScalarFunction<3> field([](const VectorN<Real, 3>& x) {
    return x[0]*x[0] + x[1]*x[1] - 2*x[2]*x[2];
});

VectorN<Real, 3> point{1.0, 1.0, 1.0};

// Laplacian: ∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²
Real fxx = Derivation::NSecDer4Partial(field, 0, 0, point);  // = 2
Real fyy = Derivation::NSecDer4Partial(field, 1, 1, point);  // = 2
Real fzz = Derivation::NSecDer4Partial(field, 2, 2, point);  // = -4

Real laplacian = fxx + fyy + fzz;
std::cout << "∇²f at " << point << " = " << laplacian << std::endl;
// Expected: 2 + 2 - 4 = 0 (harmonic in this case)
```

### Example 6: Custom Step Size for Scaled Function

```cpp
// Function with large magnitude (ill-scaled)
Real scale = 1e8;
RealFunction huge([scale](Real x) {
    return scale * sin(x);
});

Real x = 1.0;

// Automatic h (may be suboptimal)
Real der_auto = Derivation::NDer4(huge, x);

// Manual h adjustment for better accuracy
Real h_base = pow(Constants::Eps, 1.0/5.0);  // Order 4 base
Real h_adjusted = h_base * pow(scale, -1.0/5.0);  // Scale correction

Real err_manual;
Real der_manual = Derivation::NDer4(huge, x, h_adjusted, &err_manual);

// Analytical derivative
Real der_exact = scale * cos(x);

std::cout << "Exact derivative: " << der_exact << std::endl;
std::cout << "Auto h derivative: " << der_auto 
          << " (error: " << abs(der_auto - der_exact) << ")" << std::endl;
std::cout << "Manual h derivative: " << der_manual 
          << " (error: " << abs(der_manual - der_exact) 
          << ", est: " << err_manual << ")" << std::endl;
```

---

## Performance and Accuracy

### Function Evaluation Counts

| Method | First Derivative | Second Derivative | Gradient (n dims) |
|--------|------------------|-------------------|-------------------|
| Order 1 | 2 | 3 | 2n |
| Order 2 | 2 | 3 | 2n |
| Order 4 | 4-5 | 5-7 | 4n-5n |
| Order 6 | 6-7 | 8-10 | 6n-7n |
| Order 8 | 8-9 | 11-13 | 8n-9n |

**Full vector-function Jacobian:** the `CalcJacobian` overloads use
exactly $4n$ vector-function evaluations for an $m\times n$ fourth-order
Jacobian. Computing individual partial derivatives separately can cost more.

### Accuracy Comparison (Example)

**Test function:** f(x) = sin(x)·exp(-x) at x = 1.0  
**Exact:** f'(1) ≈ 0.003080 (using automatic h)

| Order | Computed | Absolute Error | Relative Error |
|-------|----------|----------------|----------------|
| 1 | 0.003081 | 1.0e-6 | 3.2e-4 |
| 2 | 0.003080 | 2.5e-9 | 8.1e-7 |
| 4 | 0.003080 | 1.2e-12 | 3.9e-10 |
| 6 | 0.003080 | 8.5e-15 | 2.8e-12 |
| 8 | 0.003080 | 4.4e-16 | 1.4e-13 |

**Observation:** Each order gains ~2-3 decimal places of accuracy

### Computational Cost vs Accuracy

**For smooth functions:**
- Order 4 is **sweet spot** (4-5 evaluations, excellent accuracy)
- Order 6-8 costly for marginal gain (unless extreme precision needed)

**For noisy/rough functions:**
- Lower order (2) more robust
- Higher order amplifies noise

### Timing (Relative)

Assuming function evaluation cost = 1 unit:

| Operation | Order 2 | Order 4 | Order 6 | Order 8 |
|-----------|---------|---------|---------|---------|
| f'(x) | 2 | 5 | 7 | 9 |
| ∇f (3D) | 6 | 15 | 21 | 27 |
| Full Jacobian (3×3, optimized order 4) | — | 12 | — | — |

---

## Common Pitfalls

### 1. Wrong Step Size

**Problem:** Manual h too large or too small

```cpp
// ❌ BAD: h too large → truncation error dominates
Real der = Derivation::NDer4(f, x, 0.1);  // Inaccurate

// ❌ BAD: h too small → roundoff error dominates
Real der = Derivation::NDer4(f, x, 1e-15);  // Unstable

// ✅ GOOD: Use automatic h
Real der = Derivation::NDer4(f, x);  // Optimal balance
```

### 2. Ignoring Error Estimates

**Problem:** Blindly trusting result without checking error

```cpp
// ❌ BAD: No error checking
Real der = Derivation::Derive(f, x);

// ✅ GOOD: Monitor error
Real err;
Real der = Derivation::Derive(f, x, &err);
if (err > tolerance) {
    // Try different method or warn user
}
```

### 3. Using High Order for Noisy Functions

**Problem:** Order 8 on discrete/noisy data amplifies noise

```cpp
// Measured data (noisy)
RealFunction noisy_data([](Real x) {
    return measured_value(x);  // Has ±0.01 noise
});

// ❌ BAD: High order amplifies noise
Real der8 = Derivation::NDer8(noisy_data, x);  // Garbage

// ✅ GOOD: Low order more robust
Real der2 = Derivation::NDer2(noisy_data, x);  // Better
```

### 4. Derivatives at Boundaries

**Problem:** A side-named wrapper is not necessarily an endpoint derivative

```cpp
RealFunction bounded([](Real x) {
    if (x < 0 || x > 1) throw DomainError("outside [0,1]");
    return x * x;
});
Real h = 0.01;
Real shifted = Derivation::NDer4Right(bounded, 0.0, h);
std::cout << shifted << '\n';  // About 0.08: derivative at 4h, NOT at zero
```

The samples are inside the domain here, but the requested endpoint derivative is zero.
Choose a justified interior point or an endpoint-specific method if the derivative at
the boundary itself is required. See [Legacy One-Sided Sampling](#legacy-one-sided-sampling)
before supplying an error pointer.

### 5. Not Scaling Inputs

**Problem:** Variables with vastly different scales

```cpp
ScalarFunction<2> f([](const VectorN<Real, 2>& x) {
    // x[0] ~ 1e-3, x[1] ~ 1e6  (poor scaling!)
    return x[0] * x[1];
});

// ❌ BAD: Automatic h suboptimal for both dimensions
VectorN<Real, 2> grad = Derivation::NDer4PartialByAll(f, point);

// ✅ GOOD: Rescale problem first
// Let u = x[0]/1e-3, v = x[1]/1e6
// Then solve in (u,v) space with similar magnitudes
```

### 6. Confusing Partial Derivative Indices

**Problem:** Wrong index order in multivariable functions

```cpp
ScalarFunction<3> f([](const VectorN<Real, 3>& x) {
    return x[0] * x[1] * x[2];  // xyz
});

VectorN<Real, 3> point{1, 2, 3};

// Indices: 0=x, 1=y, 2=z

// ∂f/∂y at point
Real df_dy = Derivation::NDer4Partial(f, 1, point);  // Correct: index 1 = y

// ❌ Common mistake: thinking indices start at 1
Real wrong = Derivation::NDer4Partial(f, 2, point);  // This is ∂f/∂z!
```

---

## See Also

**Related Documentation:**
- [Functions.md](Functions.md) - Function interfaces and creation
- [Integration.md](Integration.md) - Numerical integration methods
- [Vector_field_operations.md](Vector_field_operations.md) - Gradient, divergence, curl, Laplacian
- [ODE_system.md](ODE_system.md) - Systems requiring derivatives

**Base Classes:**
- [Vector.md](../base/Vector.md) - Dynamic vectors
- [VectorN.md](../base/VectorN.md) - Fixed-size vectors
- [Matrix.md](../base/Matrix.md) - Matrix operations

**Theory References:**
- Numerical Recipes §5.7 - Numerical Derivatives
- Fornberg, B. (1988) - "Generation of Finite Difference Formulas"
- Press et al. - "Numerical Differentiation"

---

## Runnable Examples

| Example | Source File | Description |
|---------|------------|-------------|
| Derivation Demo | [docs_demo_derivation.cpp](../../src/docs_demos/docs_demo_derivation.cpp) | Real, scalar, vector function derivatives |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/docs_demos/Release/MML_DocsApp
```
