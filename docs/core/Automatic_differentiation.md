# Automatic differentiation for geometry and fields

**Header:** `mml/core/Derivation/ADAdapters.h`

**Runnable examples:** [docs_demo_automatic_differentiation.cpp](../../src/docs_demos/core/docs_demo_automatic_differentiation.cpp)

MML's opt-in AD adapters retain a scalar-generic expression and expose the normal
`Real`-valued MML interfaces. Define a model once, then use its derivatives in
existing curve, surface, metric, and field operations without writing derivative
formulas or choosing finite-difference steps.

This integration uses **forward-mode AD**, including nested dual numbers for
higher derivatives. The existing low-level forward/reverse Jacobian and gradient
utilities remain available; see [Derivation](Derivation.md). There are no new
nonlinear-system, ODE-system, or optimization adapters.

AD removes finite-difference truncation and subtraction of nearby samples; it
does not remove floating-point roundoff. Nested elementary-function identities
can leave small residuals even when their analytical higher derivatives vanish.
Curvature and connection calculations also accumulate roundoff in products,
normalization, inversion, and contractions. Comparisons must respect the active
`Real` precision and the scale of the result.

A numerical reference has a separate accuracy budget. Differentiating numerical
connections compounds stencil roundoff (proportional to epsilon divided by the
step squared) as well as truncation. Near a cancelling Ricci scalar, an absolute
budget is needed even when the individual terms are well resolved. This does
not justify relaxing AD's analytical checks, exact polynomial identities,
provider dispatch, or failure propagation. The adapter regressions exercise
these contracts with both modular and generated float/double headers.

## Build and run the examples

Build `MML_DocsApp` in Release and run it with `--automatic-differentiation`.
The selected demo runs every C++ example below and checks its reference results;
a failed check throws rather than printing a success-shaped result.

The modular include set used by the demo is:

```cpp
#include <mml/core/Derivation/ADAdapters.h>
#include <mml/core/DifferentialGeometry/InducedMetric.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Integration/PathIntegration.h>
#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>
using namespace MML;
```

Alternatively, the generated `MML.h` contains the adapters and provider contracts.

Field `Detailed` results count calls to the Real-valued field `operator()`, not
derivative-provider operations or the adapter's internal dual-number evaluations.
Provider-only Cartesian operations can therefore have `function_evaluations=0`;
determinant divergence with metric/vector providers counts its one center field
call. See the [Detailed field contract](Field_operations.md#detailed-evaluation-contract)
for order selection, finite checking, error availability, and failure accounting.

## Authoring contract: opt-in, not automatic conversion

AD must evaluate the expression with dual-number inputs. An arbitrary
`Real`-only virtual function, function pointer, or erased `std::function` cannot
be converted into such an expression by an adapter.

- Use generic parameters (`auto` or templates), retaining their scalar type.
- For elementary functions use ADL: `using std::sin; sin(x)`, not `std::sin(x)`
  when `x` may be dual-valued. The same rule applies to `cos`, `exp`, `sqrt`, etc.
- Return scalar-generic containers. The examples use `std::array<S,N>`.
  A `Vec3Cart` contains `Real`, so it is not a suitable generic expression result.
- Give constant components the generic type too: `S(0)`, not a separately
  returned `Real` component.
- Expressions are called through `const` adapters. Their captures must remain
  valid, and any captured mutable external state determines thread safety.
- Curves take one scalar; surfaces take two scalars of the same type. Fields
  and metrics receive a fixed-size array of generic coordinates.
- Outputs must have the declared dimensions; metric outputs must be square,
  symmetric **covariant** component matrices.

All factories own their expression by value. Bounds describe the parameter
domain; as with ordinary MML curves/surfaces, evaluation does not enforce those
bounds. Queries must have finite coordinates. Invalid bounds, wrong output
shapes, nonfinite values/derivatives, and invalid indices are explicit errors.

## 1. Curves: derivatives, Frenet geometry, and arc length

`MakeADCurve<3>` is the default; `MakeADCurve<2>` creates a Cartesian planar
curve. The adapter supplies derivatives through third order. Existing curve
geometry and path length/scalar/vector line integrals acquire these derivatives
through `ICurveDerivativeProvider<N>`.

Runnable counterpart: `CurveExample`.

```cpp
auto helix = AD::MakeADCurve([](auto t) {
    using std::sin;
    using std::cos;
    using S = decltype(t);
    return std::array<S, 3>{2 * cos(t), 2 * sin(t), t};
}, 0, 2 * Constants::PI);

const Curves::ICurveCartesian3D& curve = helix;
Real curvature = curve.getCurvature(0.7);
Real torsion = curve.getTorsion(0.7);
Real length = PathIntegration::ParametricCurveLength(curve, 0, 1);
std::cout << "Helix: curvature=" << curvature
          << ", torsion=" << torsion << ", length=" << length << '\n';
```

For a helix `(a cos(t), a sin(t), b t)`, curvature is `a/(a*a+b*b)` and
torsion is `b/(a*a+b*b)`. Here the results are **0.4**, **0.2**, and
**sqrt(5)** for the length. The base-reference calls use the AD provider, not
the separately implemented analytical `HelixCurve` class.

AD supplies derivatives; it does not replace numerical quadrature or change
existing curve conventions for stationary points/undefined Frenet frames.
The runnable `CurveExample` also prints the absolute curvature error against
the ordinary Real-only curve evaluated at the same parameter:

```cpp
Curves::CurveCartesian3D numerical([](Real t) {
    return VectorN<Real, 3>{2 * std::cos(t), 2 * std::sin(t), t};
});
std::cout << std::scientific << "Helix curvature absolute error: AD=" << std::abs(curvature - Real(0.4))
          << ", finite differences=" << std::abs(numerical.getCurvature(0.7) - Real(0.4)) << std::defaultfloat << '\n';
```

## 2. Surfaces: extrinsic and intrinsic curvature

`MakeADSurface` implements the existing virtual tangent, normal, and
fundamental-form methods. Curvature and principal-direction consumers use
those methods. A degenerate AD surface normal throws `DomainError`.

The adapter also supplies an **induced metric jet**: metric values and their
first/second partials. `InducedMetric2D` consumes this capability for connections
and curvature, including through ordinary surface and metric base references.
This is a separate route from the second fundamental form; intrinsic curvature
is not silently obtained by calling the extrinsic curvature formula.

Runnable counterpart: `SurfaceExample`.

```cpp
auto sphere = AD::MakeADSurface([](auto u, auto w) {
    using std::sin;
    using std::cos;
    using S = decltype(u);
    return std::array<S, 3>{2 * sin(u) * cos(w), 2 * sin(u) * sin(w), 2 * cos(u)};
}, 0, Constants::PI, 0, 2 * Constants::PI);

Surfaces::ISurfaceCartesian& surface = sphere;
Real extrinsicK = surface.GaussianCurvature(0.8, 0.4);
Real intrinsicK = DifferentialGeometry::GaussianCurvatureIntrinsic(surface, 0.8, 0.4);
DifferentialGeometry::InducedMetric2D induced(surface);
Real scalarCurvature = induced.GetRicciScalar({0.8, 0.4});
std::cout << "Sphere: extrinsic K=" << extrinsicK << ", intrinsic K=" << intrinsicK
          << ", scalar curvature=" << scalarCurvature << '\n';
```

A radius-2 sphere has Gaussian curvature **1/4**, scalar curvature **1/2**,
and mean-curvature magnitude **1/2**. The signed mean curvature depends on the
parameterization's normal orientation.
The runnable `SurfaceExample` likewise compares the intrinsic curvature error
against the ordinary numerical surface-to-metric route:

```cpp
Surfaces::Sphere numerical(2);
Real numericalIntrinsicK = DifferentialGeometry::GaussianCurvatureIntrinsic(numerical, 0.8, 0.4);
std::cout << std::scientific << "Sphere intrinsic K absolute error: AD=" << std::abs(intrinsicK - Real(0.25))
          << ", finite differences=" << std::abs(numericalIntrinsicK - Real(0.25)) << std::defaultfloat << '\n';
```

Printed errors depend on the platform and derivative configuration. These are
comparisons at regular points, not a claim that AD always outperforms every
analytical formula or well-tuned stencil.

Extrinsic curvature needs second surface derivatives. Differentiating the
induced metric twice needs **third surface derivatives**; the nested AD route
handles those tangent derivatives too. Chart poles still require a different
regular parameterization. AD does not remove chart singularities.

## 3. Metrics: connections and Riemann without derivative stencils

`MakeADMetric<N>` takes an expression returning covariant components `g_ij`.
It supplies `g`, `dg`, and `ddg`; the connection and its derivative are computed
from those ordinary `Real` values. Matrix inversion remains the existing
checked operation, with its singularity policy:

`d(g^-1) = -g^-1 (dg) g^-1`.

No dual-number matrix inversion is required. Riemann, Ricci, scalar curvature,
Einstein tensors, and existing geodesic consumers use the same connection
contract. Geodesic integration itself is not differentiated.

Runnable counterpart: `MetricExample`.

```cpp
auto metric = AD::MakeADMetric<2>([](const auto& q) {
    using std::sin;
    using S = std::decay_t<decltype(q[0])>;
    return std::array<std::array<S, 2>, 2>{{
        {S(1), S(0)}, {S(0), sin(q[0]) * sin(q[0])}
    }};
});

const MetricTensorField<2>& geometry = metric;
VectorN<Real, 2> point{0.8, 0.4};
Real gamma = geometry.GetChristoffelSymbolSecondKind(0, 1, 1, point);
Real riemann = geometry.GetRiemannCurvatureTensor(0, 1, 0, 1, point);
Real scalar = geometry.GetRicciScalar(point);
std::cout << "Unit sphere metric: Gamma=" << gamma
          << ", R^0_101=" << riemann << ", scalar curvature=" << scalar << '\n';
```

Expected: `Gamma^0_11 = -sin(theta)*cos(theta)`,
`R^0_101 = sin(theta)^2`, and scalar curvature **2**.

The optional factory signature argument is
`MetricTensorField<N>::SignatureType::Lorentzian`; the default is Riemannian.
It declares the intended signature, not a positive-definiteness certification.
Connections require nonsingular covariant metrics. Symmetric off-diagonal
components are supported.

## 4. Fields: gradient, Hessian, Jacobian, divergence, curl, and Laplacian

`MakeADScalarField<N>` supplies first/second partials;
`MakeADVectorField<N>` supplies component first partials and a `Jacobian`.
Ordinary Cartesian, spherical, and cylindrical field operations consume the
providers. General-coordinate gradients, connection-based divergence, and
metric covariant derivatives do too.

Runnable counterpart: `FieldExample`.

```cpp
auto potential = AD::MakeADScalarField<3>([](const auto& p) {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
});
auto flow = AD::MakeADVectorField<3>([](const auto& p) {
    using S = std::decay_t<decltype(p[0])>;
    return std::array<S, 3>{-p[1], p[0], p[2] * p[2]};
});

VectorN<Real, 3> point{1, 2, 3};
auto gradient = ScalarFieldOperations::GradientCart(potential, point);
Real laplacian = ScalarFieldOperations::LaplacianCart(potential, point);
Real divergence = VectorFieldOperations::DivCart(flow, point);
auto curl = VectorFieldOperations::CurlCart(flow, point);
auto hessian = potential.Hessian(point);
auto jacobian = flow.Jacobian(point);
auto detailed = ScalarFieldOperations::GradientCartDetailed(potential, point);
std::cout << "Fields: gradient=" << gradient << ", Laplacian=" << laplacian
          << ", divergence=" << divergence << ", curl=" << curl << '\n';
```

Expected: gradient **(2,4,6)**, Laplacian **6**, divergence **6**, curl
**(0,0,2)**, Hessian **2I**, and `Jacobian(2,2) = 6`.

Coordinate conventions do not change: general-coordinate gradient/divergence
use coordinate-basis components; spherical/cylindrical convenience operations
use physical orthonormal components. See [Field operations](Field_operations.md).

`DivergenceViaDetDetailed` uses the identity
`partial_i(log(sqrt(abs(det(g))))) = trace(g^-1 partial_i(g))/2` when **both**
the vector field and metric supply derivatives and `derivative_order == 0`.
Otherwise it retains the checked numerical volume-ratio method.

Provider-backed `Detailed` results report **NoEstimate (NaN)** for error, not
zero: AD does not supply a roundoff/error bound. Their evaluation count counts
derivative-provider calls plus ordinary field evaluations (not metric
evaluations). The supplied field adapters evaluate their expression once per
partial. Timing and `ConvertToStatus` retain the usual result contract.

## Numerical references and explicit failures

No AD failure triggers an automatic numerical retry. A provider is selected
because it exists; an exception from it propagates, or is reported by the
existing `Detailed` exception policy.

Field derivative axes, vector components, and metric component indices must
be in `[0, N)`. Invalid indices throw `IndexError` before the adapter invokes
the user expression, including either axis of a second partial derivative.

Runnable counterpart: `NumericalReferenceAndFailureExample` (the polynomial
potential is the same definition as in `FieldExample`).

```cpp
VectorN<Real, 3> point{1, 2, 3};
auto numericalGradient = ScalarFieldOperations::GradientCart(potential, point, 4);

auto squareRoot = AD::MakeADScalarField<2>([](const auto& p) {
    using std::sqrt;
    return sqrt(p[0]) + p[1];
});
FieldOperationConfig config;
config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
auto failed = ScalarFieldOperations::GradientCartDetailed(squareRoot, {0, 1}, config);
```

The explicit gradient order uses finite differences even on an AD adapter.
The square-root derivative at zero is not finite: `failed.status` is
`InvalidInput`, with a diagnostic message, not a successful zero gradient.

### Exact routing boundary

| Consumer | AD/provider route | Numerical route retained |
|---|---|---|
| Cartesian curve geometry and `PathIntegration` | `ICurveDerivativeProvider<N>` | No provider: existing configurable derivative functions |
| Surface normals/tangents/forms and curvature | Existing virtual overrides | Ordinary surfaces keep their existing implementations |
| Metric connection and Riemann | `IMetricJetProvider<N>`; induced metrics forward surface capability | Ordinary metrics; explicit `GetChristoffelSymbolSecondKindNumerical` |
| Raw field operations and metric covariant derivatives | Scalar/vector provider | No provider: existing derivative functions |
| Detailed gradients, Cartesian divergence/curl/Laplacian | Provider with `derivative_order == 0` | Nonzero order retains the existing numerical path |
| Spherical/cylindrical detailed divergence/curl/Laplacian | Same provider-aware raw operation as before | These wrappers retain their existing default-stencil behavior; their config does not select a stencil |
| Determinant divergence | Both vector and metric providers, default order | Nonzero order or missing capability |

Explicit `Derivation::NDer*`, `CalcJacobian`, and `CalcHessian` remain numerical.
Typed differential-geometry field/frame/form APIs and surface-flux quadrature
have their own explicit numerical derivative paths; this addition does not
silently change them. Parallel-transport curve tangents retain their configured
numerical step, even if the metric connection is AD-backed.

## Quotient arithmetic and constant reverse outputs

Forward `Dual` division uses direct scalar division for the primal value. Free,
compound and mixed-scalar division share the same tangent implementation; compound
self-division reads both operands before updating either field. Floating-point leaf
tangents use binary-scaled, compensated product differences instead of forming the
denominator square or an overflowing reciprocal. Reverse division uses the same
leaf arithmetic for the denominator contribution, including the incoming adjoint.
This retains representable derivatives such as **-1e-200** at numerator and denominator
values **1e200**. Nested duals remain differentiable, including rational Hessians;
this is not a guarantee that all intermediate arithmetic in every nested expression
is range-safe.

`CalcGradientReverseAD` accepts a finite, untracked `ADVar` constant result and returns
exactly one zero per input coordinate. `CalcJacobianReverseAD` likewise fills constant
output rows with exact zeros, preserving the full output-by-input shape and any
interleaved variable-dependent rows. These zero results do not depend on adjoints left
by callback-side work. A constant scalar callback is evaluated once; the reverse
Jacobian retains its shape-probe evaluation plus one evaluation per output row.

Untracked NaN/infinite constant outputs throw `NumericalMethodError` rather than being
accepted as successful zero derivatives. Callback exceptions propagate unchanged.
Direct `ADVar(7).backward()` still throws `ArgumentError`: only the high-level
mathematical helpers interpret an untracked constant output as having zero derivatives.
Recorded constants remain supported through the ordinary tape sweep. Reverse-mode
division by zero still throws `DivisionByZeroError`. Domain and tape-lifetime
contracts are described below.

## Floating leaves, singular endpoints and real powers

`Dual<T>` requires a floating-point leaf: `float`, `double`, `long double`, or
nested duals over one of these. Integer/bool leaves, including integer-only
class-template deduction, are compile-time errors. Forward derivative and gradient
helpers, and both `std::vector`/MML `Vector` Jacobian and Hessian overloads, reject
integer-deduced element types. Their deprecated aliases follow the same rule. Use
floating inputs such as `Real(2)`; integer constants mixed with floating duals remain
supported.

Raw forward and reverse `asin`/`acos` do not clamp singular derivatives to finite
numbers. At either endpoint **-1** or **1**, a nonzero directional seed produces a
signed infinity; a zero seed can produce NaN. Inputs outside the real domain produce
NaN. Nested derivatives likewise retain nonfinite results, not fabricated finite
higher derivatives. The checked AD adapters reject nonfinite derivatives with their
existing `DomainError` contract.

AD power overloads require finite primal base and exponent values. A fixed zero
exponent retains **x^0 = 1**, including **0^0**, with exactly zero derivatives.
Negative bases require a fixed integer exponent; zero bases reject negative
exponents. Variable exponents require strictly positive bases, including when their
current value is zero or an integer. Invalid cases throw `DomainError` before
recording a reverse power operation.

Forward fixed exponents have no seeded derivative coefficients at any nesting
level. Reverse fixed exponents are untracked constants or recorded expressions
whose ancestry contains no input-variable nodes; this is structural classification,
not symbolic simplification of expressions such as `x - x`. Positive bases do not
require an ancestry walk. Ordinary fixed integer powers, including negative-base
powers and zero-base polynomial Hessians, retain their real derivatives. Singular
or overflowing results from otherwise valid finite inputs remain explicit IEEE
nonfinite results; checked adapters still reject them.

The following runs in `src/docs_demos/core/docs_demo_automatic_differentiation.cpp`:

```cpp
auto squared = AD::CalcDerivativeForwardAD(
    [](auto x) { using std::pow; return pow(x, 2); }, Real(0));
auto constantPower = AD::CalcDerivativeForwardAD(
    [](auto x) { using std::pow; return pow(x, 0); }, Real(0));
auto endpoint = AD::asin(AD::Dual<Real>(1, 1));
// squared is (0, 0); constantPower is (1, 0); endpoint.deriv is +infinity.
```

## Reverse tape lifetime and handle invalidation

The modular tape and handle APIs are declared in
[ReverseAD.h](../../mml/core/Derivation/ReverseAD.h).

`ADVar` is a non-owning value handle, not an owner of its `Tape`. The tape object
must be alive at the same address whenever a recorded handle is used. Destruction
is a caller lifetime precondition: generation checking does not make a dangling
tape pointer safe, detect destroyed tapes, or extend a tape's lifetime.

For a live tape, recorded handles capture its invalidation generation:

| Operation | Effect on existing recorded handles |
|---|---|
| Append nodes, including storage reallocation | Handles remain valid |
| `clear()` | All handles into that tape become invalid, even if indices are reused |
| Copy construction | Independent graph storage; source handles remain valid |
| Copy assignment | Destination handles become invalid; source handles remain valid |
| Move construction | Source handles become invalid; they do not follow the graph |
| Move assignment | Source and replaced destination handles become invalid |
| Direct self-copy or self-move assignment | No change to handles, graph or completed-sweep state |

Moved-from tapes are empty, have no completed sweep, and can record a fresh graph.
The destination retains the transferred entries and completed-sweep state. Copies
retain independent entries and their copied sweep state. Untracked constants are
unaffected by these transitions.

Invalid recorded handles throw `ArgumentError` from value, adjoint, index and tape
accessors, backward sweeps, arithmetic and recording factories. Both operands are
validated before a recording factory inserts constants or nodes. Division validates
handles before its zero-denominator check. Copying a handle preserves its generation;
it does not rebind an invalid handle to a replacement graph.

Raw integer indices and `TapeEntry` references are not generation-aware handles:
callers must use them only with the current graph and respect storage invalidation.
This is not a thread-safety or concurrent-mutation guarantee. Invalidation-counter
exhaustion throws `ResourceLimitError` before modifying a graph rather than wrapping
and resurrecting old handles; copy/move replacement can therefore throw.

The bounded docs demo checks rejection even after a new node reuses the old index:

```cpp
AD::Tape tape;
AD::ADVar old(&tape, 2);
tape.clear();
AD::ADVar replacement(&tape, 5);
bool rejected = false;
try {
    (void)old.value();
} catch (const ArgumentError& error) {
    rejected = true;
    std::cout << "Stale AD handle rejected: " << error.what() << '\n';
}
if (!rejected)
    throw std::runtime_error("Stale AD handle was accepted");
```

## Extension contracts and limitations

The interfaces in [DerivativeProviders.h](../../mml/core/Derivation/DerivativeProviders.h)
are AD-independent and return ordinary MML vectors/matrices. Analytical
implementations can supply them too:

- `ICurveDerivativeProvider<N>`: curve derivatives, orders 1 through 3.
- `IScalarDerivativeProvider<N>`: first and second scalar partials.
- `IVectorDerivativeProvider<N>`: vector component partials.
- `IMetricJetProvider<N>`: `MetricJet<N>` with `value`, `first[a]`,
  and `second[a][b]`. All matrices must be finite and consistent with the
  covariant metric expression.

`MetricTensorField::DeriveChristoffelSymbolSecondKind` is now a protected
virtual extension point. A custom analytical connection can override its
derivative without accepting dual inputs or supplying an entire metric jet.
Provider dispatch requires RTTI. The surface-jet shortcut in `InducedMetric2D`
is exact-type guarded: subclasses that alter metric components retain the
numerical fallback unless they explicitly supply their own derivative route.
The AD adapter classes are final, preventing inherited derivatives from silently
disagreeing with a subclass's changed expression.

Nested elementary functions support higher-order propagation; regression tests
cover trig/inverse trig, exponential/logarithms, roots, powers, hyperbolic
functions, and mixed exponent seeds. Nevertheless:

- AD differentiates the executed expression/branch, not an unavailable
  mathematical model. It does not certify differentiability at `abs(0)`,
  branch boundaries, or discontinuities.
- Rounding, overflow, conditioning, and geometric degeneracy still matter.
- Quadrature and ODE integration remain numerical. Removing derivative stencils
  is not a promise of exact final geometry or universally faster execution.
- Higher-order forward mode can be expensive. There is no shared mutable cache
  or global derivative-backend switch.

Related: [Curves and surfaces](Curves_and_surfaces.md),
[Metric tensor](Metric_tensor.md), [Field operations](Field_operations.md),
[Derivation](Derivation.md).
