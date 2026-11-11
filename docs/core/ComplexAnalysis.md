# Complex Analysis

MML's complex-analysis layer provides small, composable tools for numerical work with complex-valued functions: callable wrappers, complex derivatives, contour parameterizations, contour integration, winding numbers, Cauchy integral formulas, residues, and argument-principle counts.

The implementation is header-only and lives mainly in:

- `mml/base/ComplexFunction.h`
- `mml/interfaces/IComplexFunction.h`
- `mml/core/Derivation/DerivationComplex.h`
- `mml/core/ComplexAnalysis.h`

For complex root finding, use `mml/algorithms/RootFinding.h`; those algorithms are documented with the root-finding layer.

## Function Wrappers

Complex-analysis routines consume abstract function interfaces:

| Interface / wrapper | Purpose |
|---------------------|---------|
| `IComplexFunction` | Abstract `Complex -> Complex` callable |
| `IRealToComplexFunction` | Abstract `Real -> Complex` callable, used for contour parameterizations and signals |
| `ComplexFunction` | Wraps a `Complex (*)(Complex)` function pointer |
| `ComplexFunctionFromStdFunc` | Wraps a `std::function<Complex(Complex)>`, including lambdas and functors |
| `RealToComplexFunction` | Wraps a `Complex (*)(Real)` function pointer |
| `RealToComplexFunctionFromStdFunc` | Wraps a `std::function<Complex(Real)>` |

`ComplexFunctionFromStdFunc` and `RealToComplexFunctionFromStdFunc` store a copied `std::function`. If the source lambda captures local values by reference, the usual C++ lifetime rules still apply.

```cpp
#include <mml/base/ComplexFunction.h>
using namespace MML;

ComplexFunctionFromStdFunc f([](Complex z) {
	return z * z + Complex(1.0, 0.0);
});

Complex value = f(Complex(0.0, 1.0));
```

## Complex Derivatives

Complex numerical derivatives are provided through the derivation layer:

| Function | Description |
|----------|-------------|
| `Derivation::NDer1Complex(f, z)` | First-order complex derivative approximation |
| `Derivation::NDer2Complex(f, z)` | Second-order centered complex derivative |
| `Derivation::NDer4Complex(f, z)` | Fourth-order stencil |
| `Derivation::NDer6Complex(f, z)` | Sixth-order stencil |
| `Derivation::NDer4ComplexDetailed(f, z, config)` | Fourth-order derivative with structured status/error reporting |

These functions assume the target function is complex differentiable near the evaluation point.

```cpp
#include <mml/core/Derivation.h>
#include <mml/base/ComplexFunction.h>
using namespace MML;

ComplexFunctionFromStdFunc f([](Complex z) { return std::exp(z); });
Complex dz = Derivation::NDer4Complex(f, Complex(1.0, 0.25));
```

## Contours

`mml/core/ComplexAnalysis.h` defines three built-in contour parameterizations in `MML::ComplexAnalysis`:

| Contour | Parameter range | Notes |
|---------|-----------------|-------|
| `CircleContour(center, radius)` | `0` to `2*pi` | Counter-clockwise circle; exposes `center()`, `radius()`, `t_start()`, `t_end()`, and `derivative(t)` |
| `LineSegmentContour(z1, z2)` | `0` to `1` | Straight segment from `z1` to `z2`; derivative is constant |
| `ArcContour(center, radius, t1, t2)` | `t1` to `t2` | Circular arc with the supplied angle interval |

Centers, endpoints, and angles must be finite. Circle/arc radii must be positive, finite,
and representable at the center; constructors reject collapsed or overflowing geometry.
Reversing an arc's parameter limits reverses orientation.

You can also use the generic contour-integral overload with your own `gamma(t)` and `gamma_deriv(t)` callables.

## Contour Integration

Contour integration returns a `ContourIntegrationResult`:

| Field | Meaning |
|-------|---------|
| `value` | Computed complex integral |
| `error_estimate`, `error_estimated` | Absolute estimated uncertainty, or `NoEstimate` when unavailable/disabled |
| `function_evaluations` | Actual attempted function calls, including numerical derivative probes in count helpers |
| `iterations` | Maximum completed subdivision depth, zero when the initial panels suffice |
| `converged` | Whether the aggregate estimated error met the requested absolute tolerance |
| `status`, `error_message` | Success or an explicit input, arithmetic, depth, budget, or resolution failure |
| `algorithm_name`, `elapsed_time_ms` | Operation and elapsed time, also retained on captured failures |

Conversion to `Complex` is **explicit**, not implicit. `ValueOrThrow(context)` rejects unusable
results; `.value` deliberately bypasses that protection. For compatibility, `ContourIntegral`
still returns a result. Its canonical diagnostic spelling is `ContourIntegralDetailed`, with
the same contour overloads and an optional `ContourIntegrationConfig`.

```cpp
#include <mml/core/ComplexAnalysis.h>
using namespace MML;

ComplexFunctionFromStdFunc f([](Complex z) {
	return MML::Real(1.0) / z;
});

ComplexAnalysis::CircleContour unit(Complex(0.0, 0.0), 1.0);
ComplexAnalysis::ContourIntegrationConfig config;
config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
auto integral = ComplexAnalysis::ContourIntegralDetailed(f, unit, config);
Complex value = integral.ValueOrThrow("contour example"); // Approximately 2*pi*i.
```

Available overloads:

- `ContourIntegral(f, CircleContour, tol)`
- `ContourIntegral(f, LineSegmentContour, tol)`
- `ContourIntegral(f, ArcContour, tol)`
- `ContourIntegral(f, gamma, gamma_deriv, t_start, t_end, tol)`

The implementation integrates `f(gamma(t)) * gamma'(t)` with embedded nonuniform G7/K15
Gauss-Kronrod rules and largest-error-first refinement of the aggregate frontier. It checks
endpoints separately, accumulates complex components with compensation, includes a roundoff
floor, and stops on budget/depth exhaustion, stagnation, or unrepresentable subdivisions.
This replaces the five-node dyadic Simpson convergence test that could alias harmonic four.

`ContourIntegrationConfig` inherits `estimate_error`, `check_finite`, and `exception_policy`.
Its additional fields are:

| Field | Default | Meaning |
|-------|---------|---------|
| `tolerance` | `Precision::ComplexAnalysisTolerance` | Positive finite absolute tolerance in **returned-value units** |
| `max_depth` | `Precision::ComplexAnalysisMaxDepth` (30) | Nonnegative subdivision limit |
| `max_evaluations` | `Precision::ComplexAnalysisMaxEvaluations` (100000) | Positive hard limit on attempted function calls |
| `initial_panels` | `Precision::ComplexAnalysisInitialPanels` (2) | Positive initial panel count; Cauchy derivatives use at least `n+1` panels |

Zero/nonfinite tolerances and invalid configuration are rejected before function evaluation.
The float default absolute budget is chosen above the unchanged quadrature
roundoff floor for an ordinary unit-circle integral. Explicit tighter requests
still fail when their requested accuracy is unachievable; defaults do not
disable roundoff, uncertainty, or resource checks.
Default `Propagate` preserves callback/input exceptions; `ConvertToStatus` captures them with
diagnostics. Ordinary non-convergence is returned as data under either policy. Disabling error
reporting hides the estimate, not the adaptive checks. A zero parameter interval returns exact
zero without evaluating the function.

## Winding, Cauchy, Residues, and Zeros

The higher-level helpers are intentionally focused on circular contours:

| Function | Description |
|----------|-------------|
| `WindingNumber(contour, z0, tol)` | Computes and numerically certifies the integer winding estimate for a circle |
| `CauchyIntegralFormula(f, contour, z0, tol)` | Evaluates `f(z0)` from the Cauchy integral formula |
| `CauchyDerivative(f, contour, z0, n, tol)` | Computes the `n`-th derivative from the generalized Cauchy integral formula |
| `Residue(f, z0, radius, tol)` | Computes a residue by integrating around a small circle centered at `z0` |
| `ResidueSimplePole(f, z0, h)` | Fixed-step simple-pole estimate from opposite real/imaginary directions |
| `ArgumentPrinciple(f, contour, tol)` | Computes zeros minus poles inside a circle by integrating `f'(z) / f(z)` |
| `CountZeros(f, contour, tol)` | Convenience alias for pole-free functions |

Each helper has a `...Detailed` counterpart. Plain complex/int helpers reject unusable results
instead of discarding their status. Config-bearing overloads permit explicit work limits.
`Residue` and argument-principle scalar overloads retain their `DefaultToleranceStrict` defaults.

Cauchy formulas require a finite point **strictly inside** the circle and a nonnegative derivative
order whose factorial is representable. Values and errors are both scaled by `n!/(2*pi*i)`;
unachievable output accuracy is a failure, not apparent convergence in unscaled integral units.

`IntegerContourResult` retains the unrounded complex `value`, `integer_value`,
`integer_certified`, and (for numerical differentiation) `derivative_error_estimate`.
`IntegerOrThrow(context)` requires both usable integration and integer certification.
Certification checks finite/range-safe conversion and requires the complex distance to the
nearest integer **plus estimated uncertainty** to fit the requested tolerance and a unique-integer
window. Argument-principle uncertainty includes the largest sampled normalized derivative-error
estimate; derivative stencils are refined within a bounded work budget. Merely rounding the real
part is insufficient.

```cpp
#include <mml/core/ComplexAnalysis.h>
using namespace MML;

ComplexFunctionFromStdFunc polynomial([](Complex z) {
	return z * z - Complex(1.0, 0.0);
});

ComplexAnalysis::CircleContour circle(Complex(0.0, 0.0), 2.0);
int zeros = ComplexAnalysis::CountZeros(polynomial, circle);
```

For diagnostics on the previously aliased nonvanishing example:

```cpp
ComplexFunctionFromStdFunc nonvanishing([](Complex z) {
	return std::exp(z * z * z * z);
});
ComplexAnalysis::CircleContour unit(Complex(0, 0), 1);
auto count = ComplexAnalysis::CountZerosDetailed(nonvanishing, unit);
int zeros = count.IntegerOrThrow("nonvanishing example"); // 0.
```

`ResidueSimplePoleDetailed` uses the actual representable displacements, with opposite-pair
weights that cancel a linear analytic background even when rounded step sizes differ. It
rejects nonpositive/nonfinite/nonprogressing steps. This is a **fixed finite-step estimate**:
`Completed`, `converged=false`, `error_estimated=false`, and `NoEstimate`. It does not certify
the limiting residue or check that a pole is simple; use `ResidueDetailed` for adaptive
contour-based accuracy.

## Notes

- Numerical integer certification is conditional on the theorem's hypotheses and adequate sampling,
  **not a mathematical proof**. No finite black-box sampling rule can exclude every unresolved
  oscillation, zero, pole, or deliberately aliased function. Quadrature and derivative uncertainty
  are estimates, not rigorous analytic bounds.
- The contour theorem helpers do not verify analyticity or singularity placement. The caller is responsible for choosing a contour that matches the theorem being used.
- `Residue` requires a radius small enough to isolate the target singularity from any others.
- `ArgumentPrinciple` computes `f'(z)` with error-estimated `NDer4ComplexDetailed` stencils;
  these probe a neighborhood off the contour. The function must be analytic there and nonzero
  on the contour. Detected zeros/nonfinite samples or unresolved derivatives fail explicitly.
- Routine tests cover wrapper evaluation, complex derivatives, Newton/Muller complex root finding, contour classes, contour integrals, winding numbers, Cauchy formulas, residues, and argument-principle zero counts.

The examples and original regression witnesses run in `MML_DocsApp --complex-contours`
(`src/docs_demos/docs_demos_cookbook.cpp`).

## See Also

- [Numerical derivation](Derivation.md)
- [Root finding](../algorithms/Root_finding.md)
- [Functions](../base/Functions.md)