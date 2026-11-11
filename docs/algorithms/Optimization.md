# Optimization

## Callback diagnostics: Nelder-Mead and ALM

`MultidimOptimizationResult.function_evaluations` records objective callback
attempts. `iterations` and `iterations_used` record algorithm iterations, not
surrogate evaluation totals. Nelder-Mead counts initial simplex vertices,
reflection, expansion, contraction and every shrink sample at the callback
boundary. An initially converged two-dimensional simplex can therefore take
zero iterations but three objective calls.

The following example is compiled and checked by
`Docs_Demo_CallbackAccounting` in
[`src/docs_demos/algorithms/docs_demo_root_finding.cpp`](../../src/docs_demos/algorithms/docs_demo_root_finding.cpp).
Run `MML_DocsApp --callback-accounting`.

```cpp
int objectiveCalls = 0;
ScalarFunctionFromStdFunc<2> objective([&](const VectorN<Real, 2>&) {
    ++objectiveCalls;
    return Real(1);
});
Optimization::NelderMead optimizer;
const auto minimum = optimizer.Minimize(objective, VectorN<Real, 2>{0, 0});
if (minimum.function_evaluations != objectiveCalls || objectiveCalls != 3
    || minimum.iterations != 0)
    throw std::runtime_error("Nelder-Mead callback accounting failed");
```

Callback exceptions propagate unchanged. After a failure,
`NelderMead::getNumFuncEvals()` includes the throwing invocation. Each public
`Minimize` attempt resets that counter, including input-validation failures
that perform zero objective calls. The counter is also current while the
objective callback executes.

The existing Nelder-Mead budget behavior is unchanged: `maxIter` is an
evaluation threshold checked **between simplex iterations**, and also caps
iterations. Initialization always evaluates the whole simplex; an iteration
already started completes its reflection/expansion/contraction/shrink batch.
The reported callback total can consequently exceed that threshold.
Box-constrained Nelder-Mead additionally counts its final projected-objective
check.

ALM's `function_evaluations` counts the original objective, not constraint
callbacks or simplex iterations. It includes objective work in an exhausted
inner solve, bounded-wrapper final checks and ALM's final objective check;
`inner_iterations` remains a separate aggregate. The same runnable demo verifies
four objective calls despite zero inner iterations in an exhausted unbounded
two-dimensional solve.

## Callback diagnostics: 1-D minimization

`Minimization1DResult.function_evaluations` counts objective invocations made by
Golden Section, Brent and derivative-based Brent. A supplied-bracket overload
counts only fresh samples requested during that call; cached bracket values and
earlier caller work do not count. An interval overload includes its own
`BracketMinimum` work and refinement, without double counting. Failed bracketing
retains its actual callback total even though refinement takes zero iterations.
Configuration, callable and maximization wrappers preserve these counts.

Analytic derivative callbacks remain separate from objective callbacks.
For example, one derivative-Brent iteration after interval bracketing requests
five objective calls and two derivative calls in the demo below. A convergence
check that uses cached values does not invent another objective call.
Input-validation and callback exceptions retain their existing behavior:
exceptions propagate, and no result is manufactured on a throwing path.

The snippets below are compiled and checked by
`Docs_Demo_Optimization_CallbackAccounting` in
[`src/docs_demos/algorithms/docs_demo_optimization.cpp`](../../src/docs_demos/algorithms/docs_demo_optimization.cpp).
Run `MML_DocsApp --optimization-accounting`.

```cpp
int scalarCalls = 0;
auto scalar = [&](Real x) {
    ++scalarCalls;
    const Real offset = x - Real(1) / 3;
    return 1 + offset * offset;
};
Minimization1DConfig config;
config.max_iterations = 1;
const auto minimum = Minimization::BrentMinimize(scalar, -1, 0, config);
if (minimum.function_evaluations != scalarCalls || scalarCalls != 5 || minimum.iterations != 1)
    throw std::runtime_error("Brent bracketing/refinement accounting failed");
```

## Scalar Brent result status

Ordinary and derivative-based Brent report `status == AlgorithmStatus::Success`
whenever a convergence check is met, including initial bracket-width completion
and derivative Brent's accepted short-step minimum. `converged`, `IsSuccess()`
and explicit boolean conversion then agree. Interval, callable, configuration
and maximization wrappers preserve this outcome and the actual callback count.
Default results remain `NotComputed`; exhausted solves remain
`MaxIterationsExceeded`, and failed bracketing remains `AlgorithmSpecificFailure`.
Those states do not become usable merely because a partial value is present.
See [RESULT_POLICY.md](../RESULT_POLICY.md).

Using the `scalar` objective from the earlier snippet:

```cpp
scalarCalls = 0;
Minimization1DConfig successConfig;
const auto successful = Minimization::BrentMinimize(scalar, -1, 0, successConfig);
if (!successful.converged || successful.status != AlgorithmStatus::Success
    || !successful.IsSuccess() || !static_cast<bool>(successful)
    || successful.function_evaluations != scalarCalls)
    throw std::runtime_error("Brent success-status contract failed");
```

## Callback diagnostics: Powell

Powell counts initial and extrapolated objective samples plus every callback
in nested line-search bracketing, retries and refinement. Direction-replacement
line searches are included. `iterations` and `iterations_used` remain direction-set
iterations; the iteration cap is not a callback budget.

```cpp
int objectiveCalls = 0;
ScalarFunctionFromStdFunc<2> objective([&](const VectorN<Real, 2>& x) {
    ++objectiveCalls;
    return 1 + x[0] * x[0] + 2 * x[1] * x[1];
});
Optimization::Powell optimizer(Precision::OptimizationTolerance, 1);
const auto powell = optimizer.Minimize(objective, VectorN<Real, 2>{-5, 5});
if (powell.function_evaluations != objectiveCalls || optimizer.getNumFuncEvals() != objectiveCalls
    || powell.iterations != 1 || powell.converged)
    throw std::runtime_error("Powell nested-callback accounting failed");
```

`Powell::getNumFuncEvals()` is current inside the objective callback and includes
a throwing invocation. It resets before each public `Minimize` attempt, including
validation failures that perform no callbacks. Exceptions propagate unchanged.
Box-constrained Powell additionally counts the final projected-objective check.
The demo also checks partial counts after a propagated failure and the bounded
wrapper's final check.

## Powell custom direction matrices

For an `IScalarFunction<N>`, the custom-direction overload requires an exactly
`N`-by-`N` matrix with finite entries; its columns are the initial search directions.
Empty, undersized, oversized and rectangular matrices are rejected with
`MultidimOptimizationInputError` before any entry is read or objective callback is
issued. The diagnostic includes expected and actual dimensions. Rejected matrices
are unchanged; there is no implicit resizing or truncation. Correctly sized
matrices retain the existing finite-entry checks and direction-update behavior.

The same runnable demo checks rejection and counter reset after earlier solves:

```cpp
objectiveCalls = 0;
Matrix<Real> wrongDirections(3, 3, Real(0));
for (int i = 0; i < 3; ++i) wrongDirections(i, i) = 1;
const auto before = wrongDirections;
bool rejected = false;
try {
    const auto unexpected = optimizer.Minimize(objective, VectorN<Real, 2>{-5, 5}, wrongDirections);
    (void)unexpected;
} catch (const Optimization::MultidimOptimizationInputError& error) {
    if (std::string_view(error.what()).find("expected 2 x 2, got 3 x 3") == std::string_view::npos) throw;
    rejected = true;
}
if (!rejected || objectiveCalls != 0 || optimizer.getNumFuncEvals() != 0 || wrongDirections != before)
    throw std::runtime_error("Powell direction-shape preflight failed");
```

These accounting repairs do not change numerical results, convergence tests,
exception policies or iteration/batch-budget behavior.

## Bound-feasibility input contracts

`BoundConstraints::IsFeasible` and `MaxViolation` require finite coordinates in
both dynamic and fixed-size vectors. NaN and either infinity throw
`BoundConstraintsError`, including in an unbounded component; they never imply
feasibility or zero violation. The coordinate diagnostic identifies its index.
`IsFeasible` additionally requires a finite, nonnegative tolerance; zero means
exact bound comparisons. Dimension mismatches still throw `VectorDimensionError`,
and querying a default empty bound domain is an error.

Bounds may use negative infinity as a lower endpoint and positive infinity as an
upper endpoint (`Unbounded` uses both). Those are bound markers, not permitted
query coordinates. NaN bounds, reversed intervals, positive-infinite lower
endpoints and negative-infinite upper endpoints are rejected. Fixed **finite**
coordinates and one-sided bounds remain supported. ALM uses the same bound-data
validation rather than duplicating the finite-coordinate admission check.

The following snippet is compiled and checked by
`Docs_Demo_Optimization_BoundsContracts` in
[`src/docs_demos/algorithms/docs_demo_optimization.cpp`](../../src/docs_demos/algorithms/docs_demo_optimization.cpp).
Run `MML_DocsApp --bounds-contracts`.

```cpp
Optimization::BoundConstraints unit(1, 0, 1);
const Real nan = std::numeric_limits<Real>::quiet_NaN();
bool invalidRejected = false;
try {
    const auto unexpected = unit.IsFeasible(Vector<Real>{nan});
    (void)unexpected;
} catch (const Optimization::BoundConstraintsError&) {
    invalidRejected = true;
}
if (!invalidRejected)
    throw std::runtime_error("NaN feasibility query was not rejected");

const auto unbounded = Optimization::BoundConstraints::Unbounded(1);
if (!unbounded.IsFeasible(VectorN<Real, 1>{2}) || unbounded.MaxViolation(Vector<Real>{2}) != 0)
    throw std::runtime_error("Unbounded finite-coordinate feasibility changed");
```

Finite projection/clamping and box-wrapper final objective recomputation are
unchanged. Projected-gradient and box-wrapper entry points already reject
nonfinite starts before objective callbacks; this repair closes the independent
direct-query gap, not a claim that every consumer previously admitted NaN.