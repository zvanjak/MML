# Root Finding Algorithms

Comprehensive guide to finding zeros of real-valued functions using bracketing and derivative-based methods.

## Overview

Root finding algorithms solve the fundamental problem: **given function f(x), find x* such that f(x*) = 0**. This appears throughout:
- **Physics**: Equilibrium points, resonant frequencies, orbital mechanics
- **Engineering**: Circuit design, structural analysis, control systems
- **Economics**: Break-even points, Nash equilibria, pricing models
- **Numerical analysis**: Implicit equations, nonlinear systems preparation

The library provides **six scalar refinement methods** plus bracketing,
polynomial, and complex-root utilities. Scalar callables can be passed directly
or through `IRealFunction`.

## Quick Reference

| Method | Convergence | Requires Bracket | Needs Derivative | Speed | Robustness | Use Case |
|--------|-------------|------------------|------------------|-------|------------|----------|
| **BracketRoot** | N/A | No (creates bracket) | No | Fast search | ⭐⭐⭐ | **Find initial bracket** |
| **FindRootBrackets** | N/A | No (finds all) | No | Systematic | ⭐⭐⭐⭐ | **Locate all roots** |
| **FindRootBisection** | Linear | Yes | No | Slow | ⭐⭐⭐⭐⭐ | **Guaranteed convergence** |
| **FindRootFalsePosition** | Linear (may stagnate) | Yes | No | Medium | ⭐⭐⭐⭐ | **Shape-aware bracketing** |
| **FindRootSecant** | Superlinear (~1.618) | Yes in MML | No | Fast | ⭐⭐⭐ | **No derivative Newton** |
| **FindRootNewton** | Quadratic | Yes in MML | Yes (numerical) | **Fastest** | ⭐⭐ | **High-precision refinement** |
| **FindRootRidders** | Quadratic | Yes | No | Fast | ⭐⭐⭐⭐ | **Robust fast bracketing** |
| **FindRootBrent** | Superlinear-Quadratic | Yes | No | Fast | ⭐⭐⭐⭐⭐ | **Best all-around** |

**Key Insight**: Always bracket first (BracketRoot/FindRootBrackets), then refine with your method of choice. Brent is the recommended default.

### Capability Map

| Problem | Preferred API |
|---------|---------------|
| One bracketed scalar root | `FindRootBrent` |
| Guaranteed scalar interval reduction | `FindRootBisection` |
| Candidate intervals, including tangent roots | `FindRootIntervals` |
| All detectable real roots on an interval | `FindAllRealRootsInInterval` |
| Polynomial real/complex roots with diagnostics | `FindPolynomialRootsDetailed` |
| One complex-function root | `FindRootNewtonComplex` or `FindRootMuller` |
| Square nonlinear system $F(x)=0$ | `SolveNonlinearSystemNewton` |

## Configuration and Diagnostics

Every scalar refinement method follows the library-wide plain / Detailed convention
([CONTRACT_POLICY §11](../CONTRACT_POLICY.md#11-naming-convention), [RESULT_POLICY](../RESULT_POLICY.md)):
`FindRootX(f, x1, x2, xacc)` returns the `Real` root and throws `RootFindingError` when it
cannot be found; `FindRootXDetailed(f, x1, x2, xacc)` and `FindRootXDetailed(f, x1, x2, config)`
return a `RootFindingResult` and never throw for numerical failure. (Until 2.1 the result-returning
form was spelled `FindRootX(f, x1, x2, config)`; that overload remains as a deprecated alias.)

```cpp
RootFinding::RootFindingConfig config;
config.x_tolerance = 1e-12;
config.f_tolerance = 1e-12;
config.relative_tolerance = 1e-12;
config.max_iterations = 100;

auto result = RootFinding::FindRootBrentDetailed(
  [](Real x) { return x*x - 2.0; }, 1.0, 2.0, config);

if (result) {
  std::cout << result.root << " residual=" << result.function_value;
}
```

For source compatibility, `tolerance` remains the fallback absolute x and
function tolerance when the explicit fields are zero. Results consistently
report `root`, `function_value`, `x_error`, iteration and function-evaluation
counts, convergence, `AlgorithmStatus`, algorithm name, message, and elapsed
time. Invalid configuration/endpoints, non-finite values, stalls, and iteration
exhaustion are returned as structured failures by detailed overloads; legacy
overloads continue to throw `RootFindingError` on failure.

Scalar refinement accepts either the configured x criterion **or** the function-value
criterion; it does not require both. Bisection returns the evaluated midpoint whose
residual was checked, including on iteration exhaustion. Safeguarded Newton checks
endpoint and initial midpoint residuals before requesting a derivative, and validates
each new residual before accepting a small step. An accepted endpoint or initial
midpoint uses zero iterations and has no x-error estimate (`x_error` is `NoEstimate`,
which is NaN); its achieved tolerance is the absolute residual. Numerical Newton's
function-evaluation count includes its finite-difference samples.

`function_evaluations` counts attempted objective invocations, not iterations:
endpoint, initial-midpoint, numerical-derivative samples and final residual checks
all count. An analytic derivative callback is not an objective invocation.
Invalid configuration or endpoints rejected before evaluation report zero calls.
Standalone callback exceptions still propagate; they do not manufacture a result.
The runnable `Docs_Demo_CallbackAccounting` (`MML_DocsApp --callback-accounting`)
checks, for example, eight objective calls for one numerical-Newton iteration:

```cpp
int rootCalls = 0;
auto function = [&](Real x) { ++rootCalls; return x * x - 2; };
RootFinding::RootFindingConfig config;
config.max_iterations = 1;
const auto root = RootFinding::FindRootNewtonDetailed(function, 1, 2, config);
if (root.function_evaluations != rootCalls || rootCalls != 8 || root.iterations_used != 1)
    throw std::runtime_error("Newton callback accounting failed");
```

`initial_step` and `max_stagnant_iterations` are reserved source-compatibility fields:
scalar refinement solvers currently ignore them. In particular, the latter does not
enable a configurable stagnation detector, although negative values remain invalid
configuration. Method-specific degeneracy checks can still report `Stalled`.
For scalar refinement, `max_iterations == 0` retains its algorithm-specific default.

### Bracket signs and function units

After validating finite endpoint values and accepting endpoints within the function
tolerance (including either signed zero), scalar refiners require strictly opposite
signs. They compare signs directly, never by multiplying endpoint values: a tiny
opposite-sign bracket remains valid, and tiny same-sign values remain invalid.
`BracketRoot` also compares signs directly, but accepts only exact endpoint zeros;
it has no function-tolerance setting.

Secant and false-position interpolation normalize function values before dividing;
Ridders uses an exponent-safe geometric mean and a scaled hypotenuse instead of
squaring or multiplying raw residuals. Safeguarded Newton tests the proposed point
against the current bracket and compares step lengths directly. Brent's interpolation
already uses dimensionless ratios. These protections avoid sign-product and
interpolation overflow/underflow; they do not repair a function that itself returns
nonfinite values, underflows to zero, or loses accuracy through subnormal rounding.

When replacing `f` by `unit * f`, scale the absolute `f_tolerance` by `abs(unit)`;
leave x tolerances unchanged. Tolerances are not automatically rescaled, and a zero
`f_tolerance` still inherits the legacy tolerance. The plain overload has only one
fallback tolerance, so use a Detailed configuration to control x and residual units
independently.

This runnable example is covered by `Docs_Demo_RootSignScaling` in
[`docs_demo_root_finding.cpp`](../../src/docs_demos/algorithms/docs_demo_root_finding.cpp);
run `MML_DocsApp --root-sign-scaling`:

```cpp
const int exponent = std::numeric_limits<Real>::max_exponent * 3 / 4;
for (Real unit : {std::ldexp(Real(1), -exponent), Real(1), std::ldexp(Real(1), exponent)}) {
    RealFunctionFromStdFunc function([=](Real x) { return unit * (x - Real(0.25)); });
    RootFinding::RootFindingConfig config;
    config.x_tolerance = Precision::EqualityTolerance;
    config.f_tolerance = unit * Precision::EqualityTolerance;
    config.max_iterations = 16;
    const auto result = RootFinding::FindRootBrentDetailed(function, 0, 1, config);
    if (!result) throw std::runtime_error(result.error_message);
    if (std::abs(result.root - Real(0.25)) > Precision::DefaultToleranceStrict
        || std::abs(result.function_value) > config.f_tolerance)
        throw std::runtime_error("Scaled root failed its root/residual checks");
    std::cout << "unit=" << unit << " root=" << result.root
              << " residual=" << result.function_value << '\n';
}
```

Precision matters when testing these contracts: a small scalar step must still
change the represented point to exercise a new residual. Complex Newton uses a
numerical derivative, so even a linear function need not reach its exact root in
one floating-point update. Choose a residual tolerance appropriate to `Real`;
the final permitted update is checked against that same tolerance. A nonfinite
residual fixture must also leave the derivative stencil in the finite domain if
it is intended to test rejection of the updated point rather than derivative
failure.

## Nonlinear Systems

`NonlinearSystemSolvers.h` solves square systems

$$
F(x)=0,\qquad F:\mathbb{R}^n\rightarrow\mathbb{R}^n
$$

with dense Newton iteration. This differs from multidimensional optimization:
the solver directly computes $J(x_k)\Delta x=-F(x_k)$ instead of minimizing the
scalar objective $\frac12\lVert F(x)\rVert^2$.

Dynamic `Vector`/`Matrix` and fixed-size `VectorN`/`MatrixNM` entry points are
available. Supply an analytical Jacobian callback or omit it to use a
column-wise fourth-order finite-difference Jacobian:

```cpp
RootFinding::NonlinearSystemConfig config;
config.residual_tolerance = 1e-12;
config.step_tolerance = 1e-12;
config.store_trace = true;

auto function = [](const Vector<Real>& x) {
  return Vector<Real>({x[0]*x[0] + x[1]*x[1] - 1.0, x[0] - x[1]});
};
auto jacobian = [](const Vector<Real>& x) {
  return Matrix<Real>(2, 2, {2*x[0], 2*x[1], 1.0, -1.0});
};

auto result = RootFinding::SolveNonlinearSystemNewton(
  function, jacobian, Vector<Real>({0.8, 0.4}), config);
```

`numerical_jacobian_step == 0` uses the shared derivation layer's
coordinate-scaled automatic steps. A nonzero value is an explicit absolute
finite-difference step. Numerical Jacobian function evaluations are included in
the result's aggregate evaluation count.

Backtracking is enabled by default and only accepts residual-reducing steps.
`NonlinearSystemResult` reports the solution and residual vector, residual and
step norms, function/Jacobian evaluations, linear solves, iteration count,
status/message/timing, and optional iteration traces containing point, norms,
and damping. Singular Jacobians, dimension mismatches, non-finite values,
stagnation, and iteration exhaustion are returned as structured failures.

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

**Purpose**: Expand initial guess range until a finite endpoint root or sign change is detected.

**Algorithm**:
```cpp
static bool BracketRoot(const IRealFunction& func, Real& x1, Real& x2, int MaxTry = 50)
```

**Strategy**:
1. Start with [x₁, x₂]
2. Evaluate f(x₁), f(x₂)
3. If either endpoint is exactly zero, or the finite endpoint values have opposite signs: **Done!**
4. Else: Expand interval geometrically (scaleFactor = 1.6)
   - Expand in direction of smaller |f|
   - x₁ ← x₁ + 1.6·(x₁ - x₂) if |f(x₁)| < |f(x₂)|
   - x₂ ← x₂ + 1.6·(x₂ - x₁) otherwise
5. Repeat up to MaxTry expansions, checking the endpoints after the final expansion too

**Returns**: 
- `true` if a bracket or exact endpoint root is found (x₁, x₂ updated)
- `false` if the search is exhausted or an endpoint function value is non-finite

An expansion budget of zero still checks the initial endpoints. Sign comparison
does not multiply function values, so representable tiny values cannot hide a
sign change through product underflow.

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
- Tangent candidates depend on local sampling and `tangent_tolerance`
- Brackets may contain multiple roots

### FindRootIntervals - Structured Isolation

`FindRootIntervals` is the preferred sampling API. It returns
`std::vector<RootCandidate>` instead of parallel endpoint vectors:

```cpp
RootFinding::RootIsolationConfig config;
config.num_intervals = 500;
config.zero_tolerance = 1e-12;
config.merge_tolerance = 1e-10;

auto candidates = RootFinding::FindRootIntervals(
  [](Real x) { return std::sin(x); }, 0.0, 10.0, config);

for (const auto& candidate : candidates) {
  std::cout << candidate.interval.lower << " "
        << candidate.interval.upper << " "
        << candidate.residual << "\n";
}
```

`RootCandidateType::SignChange` records a finite interval and both endpoint
values. `RootCandidateType::ExactSample` records a zero-width interval at a
sample satisfying `zero_tolerance`. `RootCandidateType::TangentCandidate`
represents a local minimum of $|f|$ that refines below `tangent_tolerance`
without changing sign. Adjacent detections are merged using `merge_tolerance`,
retaining the lowest-residual representative. Distinct roots in neighboring
intervals remain separate.

Tangent detection is enabled by default. Suspicious three-point local minima
are adaptively refined with bounded golden-section minimization of $|f|$.
Neighborhoods already containing a sign change are excluded, preventing an
ordinary root from being reported twice. Set `detect_tangent_roots` to `false`
for sign-change-only scanning, tune `tangent_tolerance` to the expected function
scale, and use `tangent_refinement_iterations` to control local refinement.

An overload accepting both `f` and `df` additionally requires a sampled
derivative sign reversal or a small derivative before refinement:

```cpp
auto candidates = RootFinding::FindRootIntervals(
  [](Real x) { return (x - 0.3)*(x - 0.3); },
  [](Real x) { return 2.0*(x - 0.3); },
  0.0, 1.0, config);
```

`FindRootBrackets` remains available as a compatibility wrapper. Exact interior
sample roots are translated back to small neighboring brackets so older
refinement workflows continue to work. Because a small residual local minimum
can be near-zero without being a true root, tangent candidates are hypotheses;
downstream refinement must still verify the final residual.

### FindAllRealRootsInInterval - Isolate and Refine

`FindAllRealRootsInInterval` combines structured isolation with final
refinement and residual verification:

```cpp
RootFinding::FindAllRealRootsConfig config;
config.isolation.num_intervals = 500;
config.isolation.tangent_tolerance = 1e-10;
config.refinement.x_tolerance = 1e-12;
config.refinement.f_tolerance = 1e-12;
config.merge_tolerance = 1e-10;

auto result = RootFinding::FindAllRealRootsInInterval(
  [](Real x) { return (x - 1.0)*(x - 2.0)*(x - 2.0); },
  0.0, 3.0, config);
```

Sign-change candidates are refined with Brent, exact samples are re-evaluated
against the final function tolerance, and tangent candidates undergo another
bounded residual minimization before acceptance. The result contains sorted,
de-duplicated `roots`, accepted and optionally rejected candidate diagnostics,
candidate-found/accepted/rejected and unique-root counts, aggregate function evaluations, elapsed time, status,
and a `complete` flag. A search with no candidates is a successful complete
result with an empty root vector. If any candidate fails final verification,
the status is `ToleranceUnachievable` and successful roots remain available.
Its aggregate evaluation count includes isolation, refinement and merge-residual
checks. Each accepted/rejected candidate's refinement count covers only that
candidate's refinement, including a throwing invocation or a non-finite result.
Thus a caught refinement failure retains its actual partial count rather than zero.
Tangent-refinement counts are measured from callbacks, independently of the
configured golden-section iteration count. Disabling rejected-candidate storage
does not remove their work from the aggregate count.

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

**Implementation Safeguards** (Numerical Recipes `rtsafe`, §9.4):

1. **Bracket maintenance**: a bracket `[xl, xh]` with `f(xl) < 0 < f(xh)` is kept
   and shrunk at every iteration, so the method can never leave it.

2. **Bisection fallback**: a Newton step is taken only if it stays inside the
   bracket and at least halves the previous step (`|2f| ≤ |dx_old·f'|`);
   otherwise a bisection step `x = (xl + xh)/2` is taken instead. A zero or
   non-finite derivative (stationary point, inflection) therefore costs one
   bisection step rather than a failure.

3. **Guaranteed convergence**: for any continuous `f` with a sign change the
   method converges; the only failure modes are an unbracketed root, a
   non-finite function value, or exhausting `max_iterations`.

4. **Analytic derivative**: `FindRootNewtonDetailed(f, df, x1, x2, config)` uses a
   user-supplied `f'` instead of `NDer4` (4 fewer evaluations per step).

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
- **Linear in general**; can stagnate when one endpoint remains fixed
- Often faster than bisection on balanced functions, but no Illinois scaling is applied
- Guaranteed convergence (maintains bracket)

**Use When**:
- ✅ Bisection is too slow
- ✅ Newton is unreliable (poor initial guess)
- ✅ Function has discontinuous derivatives
- ✅ Need a simple shape-aware bracketed method

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
#include <mml/algorithms/RootFinding.h>

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
| **BracketRoot** | O(1) | ~30 expansions | ~32 | 0 |
| **FindRootBrackets** | O(n) | n sampling cells | n + 1, plus optional tangent refinement | 0 |
| **Bisection** | O(1) | ~33 iterations | ~35 | 0 |
| **Newton (analytic derivative)** | O(1) | ~5 iterations | ~8 | ~5 |
| **Newton (numerical derivative)** | O(1) | ~5 iterations | ~28 | 0 analytic callbacks |

These are illustrative work counts, not guaranteed convergence rates. Numerical
Newton's `NDer4` requests four objective samples per derivative, in addition to
two endpoints, an initial midpoint and each new-iterate residual. Fewer iterations
therefore need not mean fewer objective calls; use `function_evaluations` to compare work.

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

### Complex Roots

`RootFindingComplex.h` provides detailed complex-domain solvers:

- `FindRootNewtonComplex` uses numerical complex differentiation from an
  initial complex guess.
- `FindRootMuller` uses three complex starting points and requires no
  derivative.

Both return `ComplexRootFindingResult` with the root, residual, iteration count,
status, method name, message, and elapsed time. They check the initial current
estimate and every permitted update, including the last one, against
`abs(function_value) < tolerance`. An initial root takes zero iterations.
Unlike the scalar configuration, their explicit `max_iterations == 0` permits
no updates: the initial current estimate is checked, then a non-root reports
`MaxIterationsExceeded`. For Muller the current estimate is the third starting
point, `z2`. Negative budgets, nonpositive/non-finite tolerances, and non-finite
starting points report `InvalidInput`; non-finite sampled values or steps report
`NumericalInstability`. A rejected non-finite step retains the last evaluated
point and its matching residual rather than pairing a new point with an old value.
General polynomial roots are
also available through the Laguerre, Durand-Kerner, companion-matrix, and
low-degree polynomial APIs in `RootFindingPolynoms.h`.

### Polynomial Root Diagnostics and Polishing

Legacy polynomial functions such as `LaguerreRoots`, `EigenvalueRoots`,
`BairstowRoots`, and the quadratic/cubic/quartic convenience solvers remain
available. `FindPolynomialRootsDetailed` adds method selection and structured
verification:

```cpp
PolynomReal polynomial({4.0, 0.0, -5.0, 0.0, 1.0});
RootFinding::PolynomialRootConfig config;
config.method = RootFinding::PolynomialRootMethod::Laguerre;
config.tolerance = 1e-12;
config.polish = true;
config.polishing_iterations = 8;

auto result = RootFinding::FindPolynomialRootsDetailed(polynomial, config);
```

Each `PolynomialRootResult` reports the complex root, absolute and
coefficient-scale-normalized residuals, estimated multiplicity, polishing
iterations/status, convergence status, and message. `PolynomialRootSetResult`
adds the method, polynomial degree, deflation/polishing flags, maximum
residuals, aggregate status, timing, and evaluation count.

Polishing always evaluates Newton corrections against the original polynomial,
not a progressively deflated copy. Nearby estimates are clustered to estimate
multiplicity, and modified Newton steps use that estimate for repeated roots.
`PolishPolynomialRoots` is also public for callers that already have root
estimates from another source.

### Coefficient Scaling and Degree

The quadratic (real and complex), cubic, and quartic solvers reduce degree only when
the leading coefficient is **exactly zero**, not when its absolute magnitude is below
machine epsilon. General root methods require the stored leading coefficient to be
nonzero; trailing zero storage is invalid rather than silently changing the degree.
A nonzero constant, however small, has no roots. The detailed API rejects the zero
polynomial because its root set is not finite.

All polynomial root methods use power-of-two coefficient normalization. Common finite,
nonzero coefficient scaling therefore does not change degree decisions or impose an
absolute coefficient-size convergence floor. The stable quadratic and Cardano branches
are retained; the complex quadratic and Bairstow factor extraction also use the stable
quadratic construction. Common scaling can still round subnormal input coefficients
before the solver receives them; normalization cannot undo that information loss.
Nonfinite inputs are rejected. Normalization that cannot preserve a coefficient
component, unrepresentable arithmetic, and iterative nonconvergence are explicit
failures, not permission to discard a nonzero leading term.

Polishing and final verification use a normalized copy of the **original, undeflated**
polynomial. The relative residual is
`abs(P(z)) / sum(abs(a[k]) * abs(z)^k)` and acceptance uses that scale-independent
quantity. Its numerator and denominator must be finite; an exact zero numerator and
zero denominator at `z == 0` represent an exact root. Absolute residuals are converted
back to the caller's coefficient units and can underflow or overflow there without
changing the normalized verification criterion.

`PolynomialRootConfig.max_iterations` is forwarded to the selected method, including
the companion matrix's total QR iteration budget. `EigenvalueRoots(polynomial, maxIter)`
exposes the same limit directly; the existing one-argument overload retains its
1000-iteration default. Polishing has its separate `polishing_iterations` budget.

The following witness and scaled companion examples run with
`MML_DocsApp --polynomial-scaling`:

```cpp
Complex first, second;
const int count = SolveQuadratic(Real(1e-20), 0, Real(-1e-20), first, second);
std::cout << count << " real roots: " << first << ", " << second << '\n'; // -1 and +1.

for (Real scale : {Real(1e-200), Real(1), Real(1e200)}) {
    const PolynomReal polynomial{4 * scale, 0, -5 * scale, 0, scale};
    RootFinding::PolynomialRootConfig config;
    config.method = RootFinding::PolynomialRootMethod::CompanionMatrix;
    config.max_iterations = 500;
    const auto result = RootFinding::FindPolynomialRootsDetailed(polynomial, config);
    if (!result.IsSuccess())
        throw RootFindingError(result.error_message);
    std::cout << "Relative residual: " << result.max_relative_residual << '\n';
}
```

## Runnable Examples

| Topic | Demo File | Function |
|-------|-----------|----------|
| Scalar methods and comparison | `src/docs_demos/algorithms/docs_demo_root_finding.cpp` | `Docs_Demo_RootFinding_Comparison()` |
| Structured scalar diagnostics | `src/docs_demos/algorithms/docs_demo_root_finding.cpp` | `Docs_Demo_RootFindingDiagnostics()` |
| Bracketing and all-real-roots search | `src/docs_demos/algorithms/docs_demo_root_finding.cpp` | `Docs_Demo_FindRootBrackets()`, `Docs_Demo_FindAllRealRoots()` |
| Polynomial diagnostics and polishing | `src/docs_demos/algorithms/docs_demo_root_finding.cpp` | `Docs_Demo_PolynomialRootDiagnostics()` |
| Nonlinear-system Newton | `src/docs_demos/algorithms/docs_demo_root_finding.cpp` | `Docs_Demo_NonlinearSystemNewton()` |
| Code examples workflow | `src/code_examples/readme10_root_finding.cpp` | `Readme_RootFinding()` |

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
