# DynamicalSystem - Comprehensive Dynamical Systems Analysis

`DynamicalSystem.h` provides continuous-flow and discrete-map models, local equilibrium
analysis, finite-time Lyapunov estimates, and parameter-sweep data. Numerical indicators
are not certificates of chaos, boundedness, or global stability.

The examples use `MML::Real`, `MML`, and `MML::Systems`. Unless a snippet shows includes,
include the umbrella header and use both namespaces as below. Each example is a function
body; class definitions shown inside a body are local example types.

## Overview

```cpp
#include <mml/systems/DynamicalSystem.h>
using namespace MML;
using namespace MML::Systems;

LorenzSystem lorenz(10, 28, 8.0 / 3.0);
auto lyap = LyapunovAnalyzer::Compute(lorenz, Vector<Real>{1, 1, 1}, 20, 0.25, 0.01);
std::cout << "Lorenz finite-time maximum: " << lyap.maxExponent << "\n";
std::cout << "Positive-growth heuristic: " << lyap.isChaotic << "\n";
```

Arguments after the state are **total physical time**, **maximum physical QR interval**,
and **IVP step setting**, respectively. The interval is not a step count.
The continuous Lyapunov engine uses fixed-step RK4 for `Real=float`, with this
setting as a maximum step; wider precisions use it as the initial adaptive step.
The facade spelling is `DynamicalSystemAnalyzer<>::ComputeLyapunov`; neither engine
uses `LyapunovExponents::Compute`.

Run the bounded examples with `MML_DocsApp --systems-api`. Their counterparts are in
[`docs_demo_dynamical_system.cpp`](../../src/docs_demos/systems/docs_demo_dynamical_system.cpp).
Previously repaired contracts also have the `--fixed-points`, `--lyapunov-contracts`,
`--systems-contracts`, `--report-contracts`, and `--model-contracts` routes.

## Design Philosophy

1. **Unified Interface** - `IDynamicalSystem` extends `IODESystemParametrized` with analysis capabilities
2. **Built-in Systems** - Classic examples ready to use (Lorenz, Rössler, Hénon, Logistic Map)
3. **Comprehensive Analysis** - Fixed points, stability, bifurcations, Lyapunov exponents
4. **Discrete & Continuous** - Both ODE systems and iterative maps
5. **Explicit Contracts** - Status-bearing failures, bounded searches, and qualified numerical estimates

## Features

### IDynamicalSystem Interface

Extended interface for systems that support analysis:

```cpp
class DampedOscillator : public DynamicalSystemBase<2, 2> {
public:
    DampedOscillator(Real frequency = 1, Real damping = 0.2) {
        setParams(Vector<Real>{frequency, damping});
        _stateNames = {"position", "velocity"};
        _paramNames = {"frequency", "damping"};
    }
    void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
        derivative[0] = x[1];
        derivative[1] = -_params[0] * _params[0] * x[0] - _params[1] * x[1];
    }
    bool hasAnalyticalJacobian() const override { return true; }
    void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
        matrix = Matrix<Real>(2, 2, {0, 1, -_params[0] * _params[0], -_params[1]});
    }
    bool isDissipative() const override { return _params[1] > 0; }
};
DampedOscillator oscillator;
Vector<Real> derivative(2);
oscillator.derivs(0, Vector<Real>{1, 0}, derivative);
Matrix<Real> matrix;
oscillator.jacobian(0, Vector<Real>{1, 0}, matrix);
std::cout << "Custom oscillator derivative: " << derivative << "\n";
```

**Key Methods:**

- `getStateName(i)` - Name of state variable i
- `getParamName(i)` - Name of parameter i
- `getParamRange(i)` - Suggested exploration range for parameter i, not an admissibility constraint
- `getDefaultInitialCondition()` - Suggested starting state
- `hasAnalyticalJacobian()` - Does system provide analytical Jacobian?
- `jacobian(t, x, J)` - Compute Jacobian matrix ∂f/∂x
- `isAutonomous()` - No explicit time dependence?
- `isHamiltonian()` - Model's Hamiltonian-structure declaration
- `isDissipative()` - Model's dissipation declaration, distinct from an analyzer's finite-time sum indicator

The direct `IDynamicalSystem` numerical-Jacobian default shares its central
difference kernel with `Derivation::CalcJacobianODE`: the automatic relative
step is `Precision::CbrtEps * max(1, abs(x[j]))`. It requires a positive dimension,
an exactly sized finite state and finite time, and preflights every central stencil
before invoking the RHS. Short states therefore never reach unchecked indexing.
An unrepresentable stencil fails explicitly; it does not switch to a one-sided rule.
Each RHS evaluation must write exactly the system dimension's finite entries.
Wrong output sizes throw `VectorDimensionError`; missing/nonfinite entries or
unrepresentable derivatives throw `NumericalMethodError`. Callback exceptions
propagate. Successful evaluation uses two vector RHS calls per state coordinate;
the caller's matrix is unchanged on failure. Analytic overrides remain responsible
for their own contracts.

## Fixed Point Analysis

### Finding Fixed Points

Fixed points are equilibrium states of **autonomous** continuous systems where
dx/dt = 0. Nonautonomous systems are rejected: a root of f(0,x) alone does not
establish a time-independent equilibrium.

```cpp
VanDerPolSystem vdp(1.0);  // μ = 1.0

// Find fixed point near origin
Vector<Real> guess({0.1, 0.1});
auto fp = FixedPointFinder::Find(vdp, guess);

if (fp.status == AlgorithmStatus::Success) {
    std::cout << "Found fixed point at: " << fp.location << "\n";
    std::cout << "Type: " << to_string(fp.type) << "\n";
    std::cout << "Asymptotic stability certified by linearization? "
              << (fp.isStable ? "Yes" : "No") << "\n";
    std::cout << "Eigenvalues: ";
    for (const auto& ev : fp.eigenvalues)
        std::cout << ev.real() << "+" << ev.imag() << "i ";
    std::cout << "\n";
}
```

**FixedPoint Result Structure:**

| Field | Meaning |
|-------|---------|
| `location`, `convergenceResidual` | Current candidate and its residual norm; residual is NaN until evaluated |
| `jacobian`, `eigenvalues`, `type` | Local Jacobian, spectrum and linearized classification |
| `isStable` | Strictly decaying linear spectrum; false does not prove instability |
| `iterations` | Completed Newton updates |
| `status`, `error_message` | Finding/classification outcome and failure diagnostic |
| `function_evaluations` | Direct Newton RHS calls, excluding calls internal to Jacobians |

Always check `status`, not only a residual comparison. Failures retain the current
candidate location and any computed residual, with a specific status and diagnostic:
`InvalidInput`, `SingularMatrix`, `NumericalInstability`, or
`MaxIterationsExceeded` as appropriate. A failed/unusable Jacobian eigensolve also
makes the overall result unsuccessful, even if the candidate residual met tolerance.
No failed spectrum is published as a stability certificate.

The guess must have the system's positive dimension and finite components;
tolerance must be finite and positive, and the update budget nonnegative.
The convergence test remains the absolute residual norm `< tolerance`.
A zero budget still checks the initial state. A root reached on the final
permitted update receives its Jacobian and spectrum just like an earlier root.
`iterations` counts completed Newton updates; `function_evaluations` counts direct
Newton residual callbacks, including the initial/final checks, but not callbacks
made internally by a user-supplied or numerical Jacobian.
Uncomputed residuals/Jacobian entries are NaN rather than fabricated zeroes.

`FindFixedPointDetailed` copies the payload's status, diagnostic and evaluation
count. Unexpected callback exceptions still propagate by default;
`DynSysConfig::exception_policy = EvaluationExceptionPolicy::ConvertToStatus`
retains the existing exception-to-failure conversion.

### Fixed Point Types

| Label | Linearized spectrum |
|-------|---------------------|
| `StableNode`, `UnstableNode` | Strictly negative or positive real eigenvalues |
| `StableFocus`, `UnstableFocus` | Strictly negative or positive real parts with complex pairs |
| `Saddle`, `SaddleFocus` | Both signs of real parts; the latter includes complex pairs |
| `Center` | Purely imaginary nonzero spectrum; not a nonlinear stability certificate |
| `Unknown` | Unclassified, including mixed neutral and nonneutral modes |

`StableStarNode` and `UnstableStarNode` are retained enum labels; do not infer a
geometric-multiplicity or diagonalizability test from their presence.

Signs use the existing `Precision::DefaultToleranceStrict` band. Strictly negative
real parts certify local asymptotic stability (subject to the usual smoothness
assumptions); `isStable == false` does **not** prove instability.
Mixed neutral and decaying/expanding modes remain `Unknown` rather than being
forced into a hyperbolic focus/node label. `Center` is only a spectral label:
linearization alone does not prove nonlinear stability, bounded motion, or
periodic orbits, and no semisimplicity test is claimed.
Calling `ClassifyFixedPoint` directly validates/classifies only the supplied
Jacobian; its status then describes classification, not verification of a root
at the stored location.
Analyzer reports count linear centers and unclassified points separately from
unstable points, and include saddle-foci in the saddle count.

Runnable counterparts and regression checks are in
[`docs_demo_dynamical_system.cpp`](../../src/docs_demos/systems/docs_demo_dynamical_system.cpp).
Run `MML_DocsApp --fixed-points` for the bounded examples.

### Multiple Fixed Points

Find distinct successful fixed points from several initial guesses (this does
not guarantee finding every equilibrium):

```cpp
VanDerPolSystem vdp(1);
std::vector<Vector<Real>> guesses = {
    Vector<Real>({0.0, 0.0}),
    Vector<Real>({1.0, 1.0}),
    Vector<Real>({-1.0, -1.0}),
    Vector<Real>({1.0, -1.0})
};

auto fixedPoints = FixedPointFinder::FindMultiple(vdp, guesses);

std::cout << "Found " << fixedPoints.size() << " unique fixed points:\n";
for (const auto& fp : fixedPoints) {
    std::cout << "  Location: " << fp.location 
              << " Type: " << to_string(fp.type) << "\n";
}
```

Only successful finite results enter the collection; failed roots and failed
classifications are excluded. Duplicate locations are merged using the finite
positive `uniqueTolerance`. Invalid shared system/tolerance parameters throw
`ArgumentError` because the collection API has no result-status envelope;
individual invalid guesses are skipped. An empty valid guess collection remains
an empty result.

## Lyapunov Exponents

Lyapunov analysis estimates finite-time tangent growth rates along a trajectory.
Both continuous and discrete APIs return exponents in descending order before
computing the maximum or applying their positive-growth heuristics.

```cpp
LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);

auto result = LyapunovAnalyzer::Compute(
    lorenz,
    Vector<Real>({1.0, 1.0, 1.0}),  // Initial condition
    20.0,      // Total integration time
    0.25,      // Maximum interval between QR decompositions
    0.01       // Initial adaptive IVP step
);

std::cout << "Lyapunov exponents: " << result.exponents << "\n";
std::cout << "Max exponent: " << result.maxExponent << "\n";
std::cout << "Sum: " << result.sum << "\n";
std::cout << "Kaplan-Yorke dimension: " << result.kaplanYorkeDimension << "\n";
std::cout << "Positive-growth heuristic: " << result.isChaotic << "\n";
```

**LyapunovResult Structure:**

| Field | Meaning |
|-------|---------|
| `exponents`, `maxExponent` | Sorted finite-time spectrum and its largest entry |
| `sum` | Sum of the estimated exponents |
| `kaplanYorkeDimension` | Kaplan–Yorke estimate derived from that spectrum |
| `isChaotic` | Largest exponent exceeds `Precision::LyapunovChaosThreshold`; not a proof of chaos |
| `numOrthonormalizations` | Actual accepted QR blocks, including internal subdivisions |
| `totalTime` | Physical time covered |
| `function_evaluations` | Direct trajectory RHS calls, including rejected QR trials; excludes calls internal to supplied/default Jacobians |
| `qr_refinements` | Internal block halvings used to recover tangent-frame resolution |

**Interpretation:**

- Positive rates indicate measured tangent expansion, not necessarily a chaotic attractor.
- Near-zero rates alone do not distinguish periodic, quasi-periodic or neutral dynamics.
- Negative rates indicate measured contraction, not a global stability proof.
- A negative sum indicates average local volume contraction along the measured trajectory.
- Kaplan–Yorke is a spectrum-based dimension estimate, not a certified fractal dimension.
- Check finite-time convergence, transients, IVP accuracy and QR-interval sensitivity separately.

### Analyzer Report Contracts

`DynamicalSystemAnalyzer::Analyze` combines non-exhaustive fixed-point searches
with finite-time Lyapunov estimates. The fixed-point part requires an autonomous
system. Counts refer only to returned local numerical linearizations; neutral
and unclassified points are not relabeled as unstable. Finding no point does
not prove that no equilibrium exists.

Legacy fields and quick-query numerical behavior are unchanged:

| Indicator | Existing numerical meaning | Not established |
|-----------|----------------------------|-----------------|
| `isChaotic` / `IsChaotic` | Largest estimated exponent exceeds `Precision::LyapunovChaosThreshold` | Chaos or a bounded chaotic attractor |
| `isDissipative` / `IsDissipative` | Estimated sum is below `-Precision::LooseTolerance` | Global volume contraction or boundedness |
| `isConservative` | Absolute estimated sum is below `Precision::CoarseTolerance` | Hamiltonian structure or exact volume preservation |
| `kaplanYorkeDimension` / `GetFractalDimension` | Spectrum-based Kaplan-Yorke estimate | A measured fractal dimension; invariant-set and asymptotic-spectrum assumptions are not verified |
| `lyapunovTime` / `GetLyapunovTime` | Reciprocal largest exponent above the growth threshold; otherwise infinity | A certified predictability horizon or infinite predictability |

The two sum indicators use independent tolerances and can both be true.
Below-threshold positive growth is not relabeled as a nonpositive exponent or
regular dynamics. No boundedness or convergence conclusion follows from these
indicators. For example, `x' = 1` has a zero exponent but unbounded motion, and
`x' = x, y' = -2*y` has a negative sum yet unbounded expanding trajectories.

Default reports have zero dimensions/counts, false indicator flags, empty
collections and summary, and unavailable (`NoEstimate`-style NaN) report-level
dimension/time estimates. An empty Lyapunov spectrum means the analysis has not
run; false default flags do not classify a system. A populated report contains
computed estimates, but `Analyze` leaves its trajectory and Poincare collections
empty. Request those analyses separately and store their results explicitly if
needed. Underlying analysis exceptions still propagate.

Run `MML_DocsApp --report-contracts` for this drift example, the bounded Lorenz
report example in the systems overview, and the qualified indicator workflow:

```cpp
class ConstantDrift : public DynamicalSystemBase<1, 0> {
public:
    void derivs(Real, const Vector<Real>&, Vector<Real>& derivative) const override {
        derivative[0] = 1;
    }
    void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
        matrix = Matrix<Real>(1, 1, {0});
    }
};
ConstantDrift drift;
auto driftReport = DynamicalSystemAnalyzer<>(drift).Analyze(Vector<Real>{0}, {}, 1);
std::cout << driftReport.summary;
```

### Algorithm Details

Uses **periodic QR normalization of the integrated variational system**:
1. Integrate the complete trajectory and tangent matrix at the same physical clock.
2. Orthonormalize with scaled, twice-applied modified Gram–Schmidt.
3. Accumulate logarithmic stretches without forming an overflowing unscaled norm.
4. Average accepted stretches over the requested total time.

For continuous flows the requested QR interval is a **maximum**. A trial block
whose stretches leave the tolerance-derived resolution window, lose rank, or
overflow the variational derivative is retried from its original state with half
the interval. The shorter interval is retained for later blocks. Refinement is
bounded by `Precision::LyapunovMaxQRRefinements` per accepted block; exhaustion or
inability to advance floating-point time raises `NumericalMethodError`.
This safeguard prevents silently discarding strong contraction; it is **not** a
certified error bound for the returned exponents.

QR refinement and IVP accuracy are separate controls. In the float RK4 engine,
a resolved tangent frame can still have appreciable integration error when the
step times the growth or contraction rate is large. Refine the step setting
and compare spectra against independent references or successive step sizes;
shortening the QR interval alone is not an exponent-accuracy test. Precision-aware
test fixtures retain strict spectral bounds by choosing a step whose RK4
stability-polynomial rate bias is below the requested accuracy.

```cpp
class StrongDecay : public DynamicalSystemBase<1, 0> {
public:
    void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
        derivative[0] = -100 * x[0];
    }
    void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
        matrix.resize(1, 1);
        matrix(0, 0) = -100;
    }
};
StrongDecay decay;
auto contraction = LyapunovAnalyzer::Compute(decay, Vector<Real>{0}, 1.0, 1.0, 0.01);
std::cout << "Expected rate -100, estimated: " << contraction.maxExponent << "\n";
std::cout << "Accepted QR blocks: " << contraction.numOrthonormalizations << "\n";
std::cout << "QR refinements: " << contraction.qr_refinements << "\n";
```

Continuous requests require positive finite time, QR interval and initial step,
and a finite initial state of the system dimension. Callback outputs must have
the required shape and finite entries. `ComputeLyapunovDetailed` propagates errors
by default; with `ConvertToStatus`, it reports specific input, numerical, singular
or resource failures, diagnostics, and direct RHS work even on failure.

For discrete maps, every nonzero representable independent stretch contributes
its finite logarithm, however small. Exactly zero or numerically dependent
tangent columns raise `SingularMatrixError`: this API does **not** continue an
extended spectrum containing negative infinity. Non-finite callback/frame values
also fail explicitly; there is no identity reset or neutral-growth replacement.
These checks apply to the stored Jacobian entries: a tiny constant rounded to
zero is a zero stretch, and a constant converted to infinity is a nonfinite
callback value. Extreme-scale rank tests must also keep their perpendicular
components representable before column-relative rank checks are meaningful.
Measurement iterations must be positive and transient iterations nonnegative.
The discrete `isChaotic` heuristic retains its separate per-iteration threshold,
`Precision::DiscreteLyapunovChaosThreshold`.

The examples in this section and the tent-map example below run with
`MML_DocsApp --lyapunov-contracts`.

## Bifurcation Analysis

The continuous engine accepts an `IDynamicalSystem`, **not** an `IDiscreteMap`.
It records sampled local maxima of one state component after each transient.
Durations are physical times and the last argument is the initial adaptive step
and output sampling cadence, not a requested number of steps. Short runs illustrate
the API, not convergence to an attractor or detection of every bifurcation.

### Discrete parameter samples

For maps, iterate directly. This example restarts every sample from `x = 0.3`,
discards 256 iterations, and records the next eight. It does not classify the
number of periodic points or certify chaos:

```cpp
for (Real parameter : {Real(2.8), Real(3.2), Real(3.5), Real(3.9)}) {
    LogisticMap sample(parameter);
    Vector<Real> x{0.3};
    for (int iteration = 0; iteration < 256; ++iteration)
        x = sample.iterate(x);
    for (int iteration = 0; iteration < 8; ++iteration) {
        x = sample.iterate(x);
        std::cout << parameter << " " << x[0] << "\n";
    }
}
```

The familiar logistic period-doubling route contains periodic windows as well as
chaotic regimes; neither the parameter alone nor finitely many distinct samples
establishes a trajectory's asymptotic behavior. Initial conditions and finite
precision can also select exceptional orbits.

### Continuous-flow sweep

```cpp
VanDerPolSystem vdp(0.5);
auto diagram = BifurcationAnalyzer::Sweep(vdp, 0, 0.5, 1.5, 3, Vector<Real>{2, 0}, 0, 5, 10, 0.02);
for (std::size_t sample = 0; sample < diagram.parameterValues.size(); ++sample)
    for (Real maximum : diagram.attractorValues[sample])
        std::cout << diagram.parameterValues[sample] << " " << maximum << "\n";
```

## Detailed Evaluation Contracts

`DynSysConfig::exception_policy` defaults to `Propagate`. `ConvertToStatus`
preserves a diagnostic and the work already performed when an exception occurs:
invalid requests map to `InvalidInput`, rank failures to `SingularMatrix`,
numerical/ODE failures to `NumericalInstability`, resource bounds to
`ResourceLimitExceeded`, and exhausted searches to `MaxIterationsExceeded`.
Other callback exceptions retain `AlgorithmSpecificFailure` and their message.

Across fixed-point, Lyapunov, trajectory, Poincaré and bifurcation Detailed APIs,
`function_evaluations` counts actual **direct RHS callbacks**, including rejected
integration attempts, transient phases and calls that throw. It does not count
output samples, intersections or parameter values, and cannot observe nested
derivative calls made internally by supplied/default Jacobians.

Input dimensions and finite initial states are always checked before work.
Trajectory time and sweep phase durations must be finite and nonnegative;
initial steps and output intervals must be finite and positive. Poincaré counts
must be nonnegative, the search budget positive, the transient nonnegative, the
section value finite, and direction one of -1, 0 or +1. The section index must
belong to the state dimension and the search must advance a finite physical clock.

Valid zero-time trajectories return exactly the initial state with zero RHS work.
Zero-crossing requests return an empty successful result with zero RHS work;
negative counts are invalid. A bounded Poincaré search that finds too few
crossings retains its partial intersections with `MaxIterationsExceeded` in the
Detailed API; the simple API throws `ConvergenceError`.

`check_finite = false` disables optional integration-stage-state and section-value
checks, **not** mandatory input, callback-shape, finite RHS-output or final-result
validation. Fixed-point and Lyapunov arithmetic also keep the finite checks
needed for trustworthy convergence and spectrum calculations.

The sweep fields `numTransientSteps` and `numRecordedPoints` retain their legacy
rounded `phaseDuration / initialStep` estimates, with checked integer bounds.
Despite their names, they are not adaptive step counts, actual stored-point
counts or RHS work. Use Detailed `function_evaluations` for direct RHS work.
Exceeding the estimates' integer range is an explicit resource failure rather
than an overflowing narrowing conversion.
Sweep parameter bounds must be finite. Overflow-safe interpolation supports
finite endpoints even when their difference overflows, in either direction.
One sample evaluates the starting endpoint; multiple bifurcation samples still
require different endpoints.

This bounded example runs with `MML_DocsApp --systems-contracts`:

```cpp
class CountedDecay : public DynamicalSystemBase<1, 1> {
public:
    mutable int calls = 0;
    CountedDecay() { setParam(0, 1); }
    void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
        ++calls;
        derivative[0] = -_params[0] * x[0];
    }
};
CountedDecay decay;
auto idle = IntegrateTrajectoryDetailed(decay, Vector<Real>{1}, 0, 0.25);
auto trajectory = IntegrateTrajectoryDetailed(decay, Vector<Real>{1}, 1, 0.25);
std::cout << "Zero-time points and RHS calls: "
          << idle.points.size() << ", " << idle.function_evaluations << "\n";
std::cout << "Trajectory RHS calls: " << trajectory.function_evaluations
          << ", observed: " << decay.calls << "\n";

DynSysConfig config;
config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
auto invalid = IntegrateTrajectoryDetailed(decay, Vector<Real>{}, 0, 0.25, 0.01, config);
std::cout << "Invalid input: " << invalid.error_message << "\n";
```

### Parameter scans and box searches

`ScanChaosVsParameter` validates the parameter index, positive sample count,
finite endpoints, initial-state shape and finiteness, representable variational
dimension, and finite positive Lyapunov duration before changing a parameter.
Each sample **restarts from the supplied initial state**. A one-sample scan
evaluates the starting endpoint; constant and descending ranges are supported.
In contrast, a bifurcation sweep uses **continuation**: the final state of both
phases becomes the next sample's starting state. These state-starting policies
remain unchanged.

Bifurcation continuation also carries one **physical clock** through transient,
recording, and subsequent parameter samples. Each separate sweep starts at
`t = 0`; a zero-duration phase preserves the current clock without RHS work.
This repairs the former resets to zero at each phase, which changed the forcing
in explicitly time-dependent systems and could fabricate maxima. Such sweeps
sample successive forcing phases, not independent fixed-phase experiments.
The complete phase schedule is checked before any parameter write or RHS call:
overflowing endpoints or positive durations unable to advance floating-point
time are invalid inputs. This preflight uses a linear pass over the requested
samples without trajectory storage. Nominal step-count metadata remains
per phase, and RHS diagnostics still report actual work.

Both scans restore the original finite parameter on normal completion and on
operation failure. When restoration succeeds, the original exception is rethrown
unchanged. If the restoration setter itself throws, no restoration guarantee
is possible: `ParameterRestorationError` explicitly reports unknown model state
and retains `OperationException()` (null after a successful computation) and
`RestorationException()` as `std::exception_ptr` values. No throwing destructor
or silent success is used. Detailed status conversion reports
`AlgorithmSpecificFailure` with both diagnostics and the actual RHS work;
exception-pointer inspection requires the propagating API/policy.

`FindFixedPointsInBox` requires an autonomous positive-dimensional system,
matching finite bounds ordered on every axis, and positive grid size and budget.
Degenerate axes remain valid. A one-point-per-axis grid retains the midpoint
behavior. The full grid cardinality is checked with overflow-safe integer
arithmetic **before allocating guesses or calling the RHS**. The default limit
is `Precision::FixedPointMaxGridStarts` (100,000); the trailing `maxGridStarts`
argument permits a deliberate caller override. This bounds starting guesses,
not total solver time or memory in bytes. A finite grid is an exploratory
search, not a proof that every fixed point has been found.

These additional bounded examples also run with `MML_DocsApp --systems-contracts`,
using the `CountedDecay` class above:

```cpp
const Real bound = std::numeric_limits<Real>::max();
auto wide = BifurcationAnalyzer::Sweep(
    decay, 0, -bound, bound, 3, Vector<Real>{0}, 0, 0, 0);
DynamicalSystemAnalyzer<> analyzer(decay);
auto scan = analyzer.ScanChaosVsParameter(0, 1, 2, 1, Vector<Real>{0}, 0.02);
auto roots = analyzer.FindFixedPointsInBox(Vector<Real>{-2}, Vector<Real>{2}, 3, 3);
std::cout << "Extreme-span midpoint: " << wide.parameterValues[1] << "\n";
std::cout << "Single scan parameter: " << scan[0].first << "\n";
std::cout << "Grid roots: " << roots.size() << "\n";
```

The same bounded demo includes an analytical clock check. The dummy parameter
deliberately leaves the RHS unchanged, isolating the solution `x(t) = sin(t)`.
There is no interior maximum on `[pi/2, 3*pi/2]`. Three successive recording
windows of length `pi`, starting at zero, contain respectively one, zero, and
one interior maximum:

```cpp
class DrivenClock : public DynamicalSystemBase<1, 1> {
public:
    DrivenClock() { setParam(0, 7); }
    bool isAutonomous() const override { return false; }
    void derivs(Real time, const Vector<Real>&, Vector<Real>& derivative) const override {
        derivative[0] = std::cos(time);
    }
};
DrivenClock driven;
auto afterTransient = BifurcationAnalyzer::Sweep(
    driven, 0, 1, 2, 1, Vector<Real>{0}, 0, Constants::PI / 2, Constants::PI, 0.02);
auto continued = BifurcationAnalyzer::Sweep(
    driven, 0, 1, 3, 3, Vector<Real>{0}, 0, 0, Constants::PI, 0.02);
std::cout << "Maxima after transient: " << afterTransient.attractorValues[0].size() << "\n";
std::cout << "Continued maxima: " << continued.attractorValues[0].size()
          << ", " << continued.attractorValues[1].size()
          << ", " << continued.attractorValues[2].size() << "\n";
```

## Discrete Maps

MML includes classic discrete maps with analytical properties. They implement `IDiscreteMap<N>`
(namespace `MML`), which is an `IFunction<Vector<Real>, const Vector<Real>&>` — so `map(x)` and
`map.iterate(x)` are the same call — and carries the shared `IParametrized` mixin: `getNumParam()`,
`getParam(i)` / `setParam(i, v)` (index-checked) and `getParams()` / `setParams(v)` (size-checked),
exactly as for the continuous `IDynamicalSystem`s above. `iterate()` rejects a state of the wrong
dimension with `VectorDimensionError`.

Built-in map constructors and parameter setters require finite coefficients. Invalid values
raise `ArgumentError`; a rejected built-in bulk update leaves every parameter unchanged.
Finite coefficients outside the usual chaotic regimes remain supported, including a singular
Henon map with `b = 0`. State and parameter name accessors reject invalid indices with `IndexError`.

```cpp
HenonMap henon(1.4, 0.3);
henon.getNumParam();                       // 2
henon.setParams(Vector<Real>({1.2, 0.3})); // VectorDimensionError if the size is not 2
Vector<Real> next = henon(Vector<Real>({0.1, 0.1}));   // == henon.iterate(...)
```

### Logistic Map

```cpp
LogisticMap logistic(3.9);
auto orbit = logistic.orbit(Vector<Real>{0.3}, 32);
auto lyap = DiscreteMapLyapunov<1>::Compute(logistic, Vector<Real>{0.3}, 10000, 1000);
std::cout << "Logistic samples: " << orbit.size() << "\n";
std::cout << "Finite-iteration growth estimate: " << lyap.maxExponent << "\n";
logistic.setParam(0, 4);
std::cout << "Analytical rate at r=4: " << logistic.analyticalLyapunov() << "\n";
```

`Compute` takes measurement and discarded-transient **iteration counts**, unlike
the continuous engine's physical times. The analytic `log(2)` at `r = 4` is a
typical invariant-measure result, not the exponent of every initial condition.
For example, starting exactly at `x = 0.5` immediately reaches an endpoint orbit;
the QR API explicitly rejects its zero tangent stretch.

### Hénon Map

Classic 2D chaotic map:

```cpp
HenonMap henon(1.4, 0.3);
auto henonOrbit = henon.orbit(Vector<Real>{0, 0}, 32);
auto henonLyap = DiscreteMapLyapunov<2>::Compute(henon, Vector<Real>{0, 0}, 10000, 1000);
std::cout << "Henon samples: " << henonOrbit.size() << "\n";
std::cout << "Henon finite-iteration spectrum: " << henonLyap.exponents << "\n";
```

`orbit(x0, n)` returns the next `n` iterates, not the initial state plus `n` iterates.
It is a finite trajectory sample, not by itself an attractor analysis.

### Standard Map (Chirikov-Taylor)

Area-preserving kicked rotor:

```cpp
StandardMap standard(1);
Vector<Real> state{0.5, 1};
for (int iteration = 0; iteration < 32; ++iteration)
    state = standard.iterate(state);
Matrix<Real> jacobian;
standard.jacobian(state, jacobian);
Real determinant = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
std::cout << "Standard-map Jacobian determinant: " << determinant << "\n";
```

### Tent Map

Exactly solvable chaotic map:

```cpp
TentMap tent(2.0);  // μ = 2.0

Vector<Real> x({0.3});
for (int iteration = 0; iteration < 100; ++iteration)
    x = tent.iterate(x);

// Constant tangent-magnitude rate at mu = 2, with the implementation's kink convention.
std::cout << "Analytical λ: " << tent.analyticalLyapunov() << "\n";  // log(2) ≈ 0.693

// Compare with numerical
auto lyap = DiscreteMapLyapunov<1>::Compute(tent, x, 10000, 1000);
std::cout << "Numerical λ: " << lyap.maxExponent << "\n";
```

## Built-in Continuous Systems

Use explicit initial states when the experiment requires them; several models inherit
the zero default. Model property flags are declarations, not recomputed certificates
for every accepted parameter choice. Rossler, Van der Pol and Chua currently supply
Jacobian overrides without advertising `hasAnalyticalJacobian() == true`; this page
does not change those existing capability flags.

### Lorenz System

The quintessential chaotic system:

```cpp
LorenzSystem lorenz(10, 28, 8.0 / 3.0);
DynamicalSystemAnalyzer<> analyzer(lorenz);
auto trajectory = analyzer.IntegrateTrajectory(Vector<Real>{1, 1, 1}, 2, 0.1, 0.01);
std::cout << "Trajectory samples: " << trajectory.size() << "\n";
```

This convenience API returns state vectors, including both endpoints, at the
requested output cadence. Use an ODE solver directly when time coordinates and
solver statistics are needed.

### Rössler System

Simpler chaotic attractor:

```cpp
RosslerSystem rossler(0.2, 0.2, 5.7);
auto rosslerLyap = LyapunovAnalyzer::Compute(rossler, Vector<Real>{1, 1, 0}, 20, 0.25, 0.01);
std::cout << "Rossler finite-time maximum: " << rosslerLyap.maxExponent << "\n";
```

The short interval is for API demonstration, not a converged reference exponent.

### Van der Pol Oscillator

Self-sustained oscillations:

```cpp
VanDerPolSystem vdp(1);
Vector<Real> initial{2, 0};
DynamicalSystemAnalyzer<> vdpAnalyzer(vdp);
auto points = vdpAnalyzer.IntegrateTrajectory(initial, 10, 0.1, 0.01);
vdp.setParam(0, 5);
std::cout << "Van der Pol trajectory points: " << points.size() << "\n";
```

Supply the nonzero initial state explicitly. The current default is the inherited
zero vector, an equilibrium; it does not start at `(2, 0)` or automatically find a cycle.

### Duffing Oscillator

The usual driven two-state equation is represented with a third phase state:

```cpp
DuffingSystem duffing(0.3, -1, 1, 0.5, 1.2); // delta, alpha, beta, gamma, omega
Vector<Real> duffingDerivative(3);
duffing.derivs(0, Vector<Real>{1, 0, 0}, duffingDerivative);
std::cout << "Duffing phase rate: " << duffingDerivative[2] << "\n";
std::cout << "Current autonomy metadata: " << duffing.isAutonomous() << "\n";
```

State order is `[x, v, theta]`; the equations use `cos(theta)` and `theta' = omega`.
There is no explicit time dependence in this lifted RHS. Its current `isAutonomous()`
metadata nevertheless returns false, so the equilibrium-search API rejects it.
This documentation reflects that existing mismatch without changing the model or flag.

### Chua's Circuit

Electronic chaos generator with double scroll:

```cpp
ChuaCircuit chua(15.6, 28, -1.143, -0.714);
Vector<Real> chuaDerivative(3);
chua.derivs(0, Vector<Real>{0.1, 0.1, 0.1}, chuaDerivative);
std::cout << "Chua derivative: " << chuaDerivative << "\n";
```

## Parameter Manipulation

The seven built-in continuous models inherit parameter management from `DynamicalSystemBase`.
Their constructors and public scalar/bulk setters reject NaN and either infinity with
`ArgumentError`. A rejected built-in update leaves every parameter unchanged; wrong-sized
bulk vectors raise `VectorDimensionError` before any parameter is changed.

`getParamRange()` describes a **suggested exploration range**, not a hard model domain.
The Lorenz, Rossler, Van der Pol, Duffing, Chua and Henon-Heiles models accept all finite
coefficients, including those outside their suggested boxes. These checks do not certify
physical realism, bounded trajectories, or representable evaluations for arbitrary states
and extreme coefficients; analysis callback checks still report numerical failures.

The double pendulum additionally requires **strictly positive masses and lengths**.
Gravity can be any finite signed value, including zero; negative gravity reverses the
chosen gravitational direction. There is no arbitrary positive-value cutoff or enforcement
of the suggested mass/length boxes. Floating-point overflow/underflow can still occur for
extreme scales. Its denominator uses the equivalent positive form
`m1 + m2 * sin(theta1 - theta2)^2`, avoiding cancellation for aligned rods at large mass ratios.

State names, parameter names/ranges, and invariant indices are checked: negative indices
and indices at or beyond their respective counts raise `IndexError`. A model with no
parameters or invariants has no valid index in that category.

`DynamicalSystemBase<N, P>` requires `N > 0` and `P >= 0` at compile time;
`IDiscreteMap<N>` likewise requires `N > 0`. Zero-parameter continuous models remain supported.
Custom derived models should initialize through `setParams()` and implement additional
per-parameter constraints in the protected, non-mutating `ValidateParameter()` hook.
Direct writes to protected storage or overrides of public setters remain the author's
responsibility. The shared `IParametrized` interface does not promise transactional behavior
for arbitrary user implementations.

Runnable domain and rejection examples are in `Docs_Demo_ModelDomains()` and can be selected
with `MML_DocsApp --model-contracts`:

```cpp
LorenzSystem lorenz;
lorenz.setParam(1, 150); // Suggested rho range is [0, 100], not an admissibility limit.
std::cout << "Accepted rho outside suggested range: " << lorenz.getParam(1) << "\n";

DoublePendulumSystem pendulum(1, 1, 1, 1, -9.80665);
pendulum.setParam(4, 0); // Zero or reversed gravity is supported.
const auto original = pendulum.getParams();
try {
    pendulum.setParams(Vector<Real>{2, 1, 0, 1, 0});
}
catch (const ArgumentError& error) {
    std::cout << "Rejected zero length: " << error.what() << "\n";
}
std::cout << "Parameters unchanged: " << (pendulum.getParams() == original) << "\n";
```

The usual valid parameter access remains unchanged:

```cpp
LorenzSystem lorenz;

// Set individual parameters
lorenz.setParam(0, 10.0);   // σ
lorenz.setParam(1, 28.0);   // ρ
lorenz.setParam(2, 8.0/3.0);// β

// Get individual parameters
Real sigma = lorenz.getParam(0);

// Set all parameters at once
lorenz.setParams(Vector<Real>({10.0, 28.0, 8.0/3.0}));

// Get all parameters
Vector<Real> params = lorenz.getParams();

// Query parameter names and ranges
std::cout << "Param 0: " << lorenz.getParamName(0) << "\n";  // "sigma"
auto [min, max] = lorenz.getParamRange(1);
std::cout << "Suggested rho range: [" << min << ", " << max << "]\n";  // [0, 100]
```

## Advanced Features

### Poincaré Sections

Define sections through phase space:

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<> analyzer(lorenz);
PoincareSection<Real> section(2, 27, +1);
auto crossings = analyzer.ComputePoincareSection(Vector<Real>{1, 1, 20}, section, 4, 0.01, 50, 1);
std::cout << "Bounded Poincare crossings: " << crossings.size() << "\n";
```

Arguments after the section are the requested crossing count, initial adaptive
step, search-time budget, and discarded transient time. Returned vectors are full
states at interpolated crossings, not just the two coordinates one might plot.
The simple API throws if the budget is exhausted; the Detailed API retains partial
crossings and a failure status, as described above.

### Invariant Quantities

Use the supplied invariant API rather than declaring an incomplete abstract example:

```cpp
HenonHeilesSystem hamiltonian;
auto state = hamiltonian.getDefaultInitialCondition();
std::cout << hamiltonian.getInvariantName(0) << ": " << hamiltonian.computeInvariant(0, state) << "\n";
std::cout << "Independent E parameter: " << hamiltonian.getParam(0) << "\n";
DoublePendulumSystem pendulum;
std::cout << "Pendulum energy: " << pendulum.computeInvariant(0, pendulum.getDefaultInitialCondition()) << "\n";
```

The current Henon-Heiles `E` parameter is independent metadata: it does not choose
the initial condition or alter the RHS. Its default value is `0.1`, while the
default state has computed energy `0.065`. Use `computeInvariant()` or `getEnergy()`
to evaluate the state actually integrated. This example evaluates energies at
individual states; it does not establish conservation by an arbitrary numerical solver.

### Custom System Example

Complete example of defining a new system:

```cpp
class LotkaVolterraSystem : public DynamicalSystemBase<2, 4> {
public:
    LotkaVolterraSystem(Real alpha = 1, Real beta = 0.1, Real gamma = 1.5, Real delta = 0.075) {
        setParams(Vector<Real>{alpha, beta, gamma, delta});
        _stateNames = {"prey", "predator"};
        _paramNames = {"alpha", "beta", "gamma", "delta"};
    }
    void derivs(Real, const Vector<Real>& x, Vector<Real>& derivative) const override {
        derivative[0] = _params[0] * x[0] - _params[1] * x[0] * x[1];
        derivative[1] = _params[3] * x[0] * x[1] - _params[2] * x[1];
    }
    bool hasAnalyticalJacobian() const override { return true; }
    void jacobian(Real, const Vector<Real>& x, Matrix<Real>& matrix) const override {
        matrix = Matrix<Real>(2, 2, {
            _params[0] - _params[1] * x[1], -_params[1] * x[0],
            _params[3] * x[1], _params[3] * x[0] - _params[2]});
    }
};
LotkaVolterraSystem populations;
auto equilibrium = FixedPointFinder::Find(populations, Vector<Real>{18, 12});
if (equilibrium.status != AlgorithmStatus::Success)
    throw std::runtime_error(equilibrium.error_message);
std::cout << "Equilibrium populations: " << equilibrium.location << "\n";
std::cout << "Linearized type: " << to_string(equilibrium.type) << "\n";
```

For these parameters the coexistence equilibrium is `(20, 10)`. A `Center` label
describes its linearization, not a general proof of nonlinear stability or periodicity.
The example is mathematical population dynamics; any desired positivity domain
for parameters or state variables needs its own model-specific contract.

## Integration with ODE Solvers

All continuous systems work with MML's ODE solvers:

```cpp
#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODESolverAdaptive.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>

LorenzSystem lorenz(10, 28, 8.0 / 3.0);
RungeKutta4_StepCalculator rk4;
ODESystemFixedStepSolver fixedSolver(lorenz, rk4);
auto fixedSolution = fixedSolver.integrate(Vector<Real>{1, 1, 1}, 0, 1, 100);

ODEAdaptiveIntegrator<DormandPrince5_Stepper> adaptive(lorenz);
ODEIntegratorConfig config;
config.initial_step_size = 0.01;
config.output_interval = 0.1;
auto adaptiveSolution = adaptive.integrate(Vector<Real>{1, 1, 1}, 0, 1, config);
std::cout << "Fixed-step stored points: " << fixedSolution.size() << "\n";
std::cout << "Adaptive stored points: " << adaptiveSolution.size() << "\n";
```

The fixed-step argument `100` is an **interval count**, producing 101 stored
points. Adaptive `output_interval` is an output cadence, not the adaptive step size;
the example stores 11 points. `size()` counts saved points and `getNumSteps()`
counts saved intervals; neither is the adaptive solver's attempted-step count.
The system and step calculator must outlive the solver that references them.

## Performance Tips

1. **Analytical Jacobian** - A correct override can avoid finite-difference RHS calls; verify it independently and advertise `hasAnalyticalJacobian()` accurately. No fixed speedup is promised.
2. **Adaptive Stepsize** - Trajectory, Lyapunov, bifurcation, and Poincaré analysis use adaptive Dormand-Prince integration; `h` supplies the initial step and bifurcation sampling cadence
3. **Transient Policy** - Discrete QR accepts an explicit transient count; continuous `Compute` starts from its supplied state without an automatic burn-in. Sweep and section APIs have their own explicit phase durations.
4. **Duration** - Check estimates against increasing durations and relevant model time scales; there is no universal number of orbital periods that guarantees convergence.
5. **QR Frequency** - Check interval sensitivity and numerical resolution; the continuous engine can shorten unresolved QR blocks, but does not certify exponent error.

## Common Patterns

### Complete Analysis Workflow

```cpp
LorenzSystem lorenz(10, 28, 8.0 / 3.0);
DynamicalSystemAnalyzer<> analyzer(lorenz);
auto fixedPoint = analyzer.FindFixedPoint(Vector<Real>{0, 0, 0});
if (fixedPoint.status != AlgorithmStatus::Success)
    throw std::runtime_error(fixedPoint.error_message);
std::cout << "Origin linearized type: " << to_string(fixedPoint.type) << "\n";
auto lyap = analyzer.ComputeLyapunov(Vector<Real>{1, 1, 1}, 20, 0.25, 0.01);
auto bifurcation = analyzer.ComputeBifurcation(1, 24, 28, 3, Vector<Real>{1, 1, 1}, 2, 2, 4);
auto trajectory = analyzer.IntegrateTrajectory(Vector<Real>{1, 1, 1}, 2, 0.1, 0.01);
std::cout << "Spectrum: " << lyap.exponents << "\n";
std::cout << "Sweep samples: " << bifurcation.parameterValues.size() << "\n";
std::cout << "Trajectory samples: " << trajectory.size() << "\n";
```

Each call is a separate experiment from its supplied state. The sweep itself
continues state/time between samples and restores `rho` afterward. These bounded
durations demonstrate composition; they are not a complete scientific convergence study.

### Practical Chaos Assessment

Finite-time indicators are not a chaos certificate:
1. Check positive-growth estimates for transient, duration, IVP and QR sensitivity.
2. Establish boundedness or an appropriate invariant set independently.
3. Investigate aperiodicity and sensitive dependence with model-specific evidence;
   finite observations cannot prove that a trajectory never repeats.
4. Do not substitute a negative exponent sum for a boundedness check.

```cpp
LorenzSystem lorenz(10, 28, 8.0 / 3.0);
auto lyap = LyapunovAnalyzer::Compute(lorenz,
    Vector<Real>({1.0, 1.0, 1.0}), 20.0, 0.25, 0.01);
if (lyap.isChaotic && lyap.sum < -Precision::LooseTolerance) {
    std::cout << "Positive tangent growth with a negative finite-time sum\n";
    std::cout << "Kaplan-Yorke estimate: " << lyap.kaplanYorkeDimension << "\n";
    std::cout << "Boundedness and chaotic-attractor assumptions require independent evidence\n";
}
```

## Implementation Details

| Feature | Algorithm | Reference |
|---------|-----------|-----------|
| Fixed Points | Newton-Raphson with Jacobian | Numerical Recipes |
| Stability | Eigenvalues of Jacobian | Linear Algebra |
| Lyapunov Exponents | Continuous QR method | Benettin et al. 1980 |
| Bifurcation | Parameter continuation | Strogatz, Nonlinear Dynamics |
| Discrete Maps | Direct iteration | Ott, Chaos in Dynamical Systems |

## See Also

- [IODESystem.h](../../mml/interfaces/IODESystem.h) - Base ODE system interface
- [Adaptive ODE integration](../../mml/algorithms/ODESolvers/ODESolverAdaptive.h) - Integration methods and controls
- [Eigenvalue solvers](../../mml/algorithms/EigenSystemSolvers.h) - Eigenvalue computation
- [Fixed-point engine](../../mml/systems/DynamicalSystem/FixedPointAnalysis.h) - Newton search and classification

## References

1. Strogatz, S. H. (2018). *Nonlinear Dynamics and Chaos*
2. Ott, E. (2002). *Chaos in Dynamical Systems*
3. Sprott, J. C. (2003). *Chaos and Time-Series Analysis*
4. Benettin, G. et al. (1980). "Lyapunov Characteristic Exponents for smooth dynamical systems"
