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

You can also use the generic contour-integral overload with your own `gamma(t)` and `gamma_deriv(t)` callables.

## Contour Integration

Contour integration returns a `ContourIntegrationResult`:

| Field | Meaning |
|-------|---------|
| `value` | Computed complex integral |
| `error_estimate` | Absolute error estimate used by the adaptive rule |
| `function_evaluations` | Number of integrand evaluations |
| `converged` | Whether the adaptive integration reached its tolerance before maximum depth |

The result converts implicitly to `Complex` for convenience, but using the named fields keeps diagnostics visible.

```cpp
#include <mml/core/ComplexAnalysis.h>
using namespace MML;

ComplexFunctionFromStdFunc f([](Complex z) {
	return REAL(1.0) / z;
});

ComplexAnalysis::CircleContour unit(Complex(0.0, 0.0), 1.0);
auto integral = ComplexAnalysis::ContourIntegral(f, unit);
// integral.value is approximately 2*pi*i.
```

Available overloads:

- `ContourIntegral(f, CircleContour, tol)`
- `ContourIntegral(f, LineSegmentContour, tol)`
- `ContourIntegral(f, ArcContour, tol)`
- `ContourIntegral(f, gamma, gamma_deriv, t_start, t_end, tol)`

The implementation uses adaptive Simpson integration on the real parameter interval after rewriting the integral as `integral f(gamma(t)) * gamma'(t) dt`.

## Winding, Cauchy, Residues, and Zeros

The higher-level helpers are intentionally focused on circular contours:

| Function | Description |
|----------|-------------|
| `WindingNumber(contour, z0, tol)` | Computes `(1 / 2*pi*i) integral dz / (z - z0)` for a `CircleContour`, rounded to an integer |
| `CauchyIntegralFormula(f, contour, z0, tol)` | Evaluates `f(z0)` from the Cauchy integral formula |
| `CauchyDerivative(f, contour, z0, n, tol)` | Computes the `n`-th derivative from the generalized Cauchy integral formula |
| `Residue(f, z0, radius, tol)` | Computes a residue by integrating around a small circle centered at `z0` |
| `ResidueSimplePole(f, z0, h)` | Approximates the simple-pole residue using the limit `(z - z0) f(z)` from several directions |
| `ArgumentPrinciple(f, contour, tol)` | Computes zeros minus poles inside a circle by integrating `f'(z) / f(z)` |
| `CountZeros(f, contour, tol)` | Convenience alias for pole-free functions |

```cpp
#include <mml/core/ComplexAnalysis.h>
using namespace MML;

ComplexFunctionFromStdFunc polynomial([](Complex z) {
	return z * z - Complex(1.0, 0.0);
});

ComplexAnalysis::CircleContour circle(Complex(0.0, 0.0), 2.0);
int zeros = ComplexAnalysis::CountZeros(polynomial, circle);
```

## Notes

- The contour theorem helpers do not verify analyticity or singularity placement. The caller is responsible for choosing a contour that matches the theorem being used.
- `Residue` requires a radius small enough to isolate the target singularity from any others.
- `ArgumentPrinciple` computes `f'(z)` numerically with `Derivation::NDer4Complex`, so it is sensitive to zeros on or very near the contour.
- Routine tests cover wrapper evaluation, complex derivatives, Newton/Muller complex root finding, contour classes, contour integrals, winding numbers, Cauchy formulas, residues, and argument-principle zero counts.

## See Also

- [Numerical derivation](Derivation.md)
- [Root finding](../algorithms/Root_finding.md)
- [Functions](../base/Functions.md)