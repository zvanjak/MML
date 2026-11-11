# MML Result-Reporting Policy

**Status:** Adopted (epic MinimalMathLibrary-wewi.3) · **Scope:** every public API in `mml/` that
returns a result object · **Companion:** [CONTRACT_POLICY.md](CONTRACT_POLICY.md) §6, §8, §10 ·
**Runnable counterpart:** `src/docs_demos/core/docs_demo_result_policy.cpp` (`MML_DocsApp`)

The Detailed API (`XDetailed()` returning a result struct built on
[`AlgorithmTypes.h`](../mml/base/AlgorithmTypes.h)) is only as honest as its weakest populator.
Analysis v4 found `converged = true` on fixed quadrature rules, literal `0.0` residuals, invented
error estimates and implicit `operator Real()` that silently discarded a failed flag. This page
states what every field means, when it may be filled, and which result type each layer returns.
It is normative: a result that violates it is a bug, not a style issue.

---

## 1. The one rule

> **A result claims only what it checked.**

- `status == Success` and `converged == true` **only after a convergence or accuracy criterion
  was evaluated and met** (residual below tolerance, bracket width below tolerance, successive
  estimates agreeing, ...).
- A method that has **no criterion** (fixed Gauss rule, direct LU, closed-form formula, a
  geometric fallback such as fan triangulation, a plain midpoint sum) reports
  `status == Completed`, `converged == false`, and leaves every estimate it did not compute as
  `NoEstimate`.
- A numeric estimate field (`error_estimate`, `achieved_tolerance`, `residual`, `max_residual`,
  `x_error`) holds either a **computed** number or `NoEstimate` (`NaN`). Never a literal `0.0`
  standing in for "unknown". Test with `HasEstimate(x)`.

## 2. Status vocabulary (`AlgorithmStatus`)

| Status | Meaning | `converged` | payload usable |
|---|---|---|---|
| `Success` | criterion checked **and met** | `true` | yes |
| `Completed` | no criterion exists; ran to completion | `false` | yes |
| `MaxIterationsExceeded` | iteration budget exhausted before the criterion was met | `false` | partial — best iterate so far |
| `ResourceLimitExceeded` | storage or allocation limit prevented further progress | `false` | partial, if the algorithm retained one |
| `Stalled` | no progress between iterations | `false` | partial |
| `ToleranceUnachievable` | requested tolerance below what the method/precision can deliver | `false` | partial |
| `NumericalInstability` | NaN/Inf/ill-conditioning detected | `false` | no |
| `SingularMatrix` | singular or numerically singular operator | `false` | no |
| `InvalidInput` | preconditions violated (only when the API chose to report instead of throw) | `false` | no |
| `AlgorithmSpecificFailure` | anything else; `error_message` explains | `false` | no |

`IsUsable(status)` / `result.IsSuccess()` / `if (result)` are all **`Success || Completed`**.
Fixed Gauss-Kronrod rules return `Completed` while retaining their computed embedded
error estimates; having an estimate is not evidence that an accuracy target was
checked. Such results remain usable with `converged=false`. Adaptive GK separately
tests its tolerance and reports `Success` or a nonconvergence failure.
Iterative algorithms that can never be `Completed` may override `IsSuccess()` as
`converged && status == Success` (RootFinding does); the two definitions agree for them.

**Partial payloads** (`MaxIterationsExceeded`, `ResourceLimitExceeded`, `Stalled`, `ToleranceUnachievable`) must still be
the best estimate the algorithm has — never zero-initialised or stale — and the estimate fields
must describe *that* payload (e.g. the residual of the returned `x`, not of some earlier iterate).

## 3. Field contract

| Field | Populate when | Otherwise |
|---|---|---|
| `converged` | criterion met | `false` |
| `status` | always | — |
| `iterations_used` | always (0 for one-shot methods; refinement levels for Romberg-type methods) | — |
| `achieved_tolerance` / `error_estimate` / `residual` | computed from the returned payload | `NoEstimate` |
| `error_message` | status is not `Success`/`Completed`; specific and actionable (`"Failed to converge after 100 iterations (residual: 1.5e-4)"`) | empty |
| `algorithm_name` | always | — |
| `elapsed_time_ms` | always (use `AlgorithmTimer`) | — |
| `function_evaluations` | whenever the algorithm can count them | `0` |
| `error_estimated` (where it exists) | derived convenience, must equal `HasEstimate(error_estimate)` | — |

Helpers: `MakeSuccessResult`, `MakeCompletedResult`, `MakeFailureResult`,
`MakeEvaluationSuccessResult`, `MakeEvaluationFailureResult` in `AlgorithmTypes.h` set these
consistently; prefer them to hand-filling.

## 4. Plain / Detailed pairs

Every algorithm exposes the pair

```cpp
FamilyResult XDetailed(input, const FamilyConfig& config = {});   // never throws on non-convergence
PayloadType  X(input, ...);                                       // thin shim over XDetailed
```

The Detailed form is spelled either `XDetailed(...)` (integration, derivation, DAE, ...) or as
the config overload `X(input, const FamilyConfig&)` (RootFinding, CONTRACT_POLICY §10); which
spelling wins library-wide is decided by the naming epic (wewi.5). The semantics below apply to
both.

- `XDetailed` reports non-convergence as **data** (status + partial payload). It throws only
  for invalid input (CONTRACT_POLICY §6), or for runtime failures when
  `EvaluationExceptionPolicy::Propagate` is selected.
- The plain `X` **returns the payload** and **throws `ConvergenceError`** (or the family's
  specific exception) when the Detailed result is not usable. So `Real I = IntegrateTrap(f, a, b);`
  is the normal call, and it can never hand back a silently non-converged number — that was
  the Poincaré-section and Newton defects fixed in the P0 sprint.
- A result type may offer `explicit operator PayloadType()` or a named `value` member; it is
  **never implicit**. `Real I = IntegrateTrapDetailed(f, a, b);` must not compile;
  `Real I = IntegrateTrapDetailed(f, a, b).value;` is the explicit opt-out that visibly discards
  diagnostics. (This supersedes the former CONTRACT_POLICY §10 allowance for implicit
  conversion.)
- Result types and Detailed functions are `[[nodiscard]]`.
- Integration is the reference implementation: one `BasicIntegrationResult<T>` for 1D, improper,
  Gauss-Kronrod, 2D/3D, adaptive, path, surface, Monte Carlo (`IntegrationResult`) and contour
  integrals (`ContourIntegrationResult = BasicIntegrationResult<Complex>`); plain `IntegrateX(...)`
  returns `Real`, `IntegrateXDetailed(...)` (same arguments, plus a config-bearing overload)
  returns the result; the result-returning cores are `IntegrateXCore` (internal).

## 5. Which result type per layer

| Layer | Returns | Notes |
|---|---|---|
| `mml/base` | **no result objects** — values or exceptions | containers, polynomials, geometry primitives are total functions or throw (`InterpolationStatus` is the one tolerated exception: it is a *query* on a constructed interpolant, not an algorithm outcome) |
| `mml/core`, `mml/algorithms`, `mml/systems` — iterative | `{Family}Result : IterativeResultBase` | eigen, root finding, optimization, ODE/DAE, Krylov, nonlinear systems |
| same — one-shot evaluation | `{Family}Result : EvaluationResultBase` | integration, differentiation, triangulation, analyzers |
| `mml/tools` (serializer, persistence, data_loader, visualizer launcher) | **one** `IOResult` | `success`, `IOError error`, `message`, `filename`, `offset` and `offset_kind`; explicit boolean conversion |

**Domain status enums** (`SolverStatus` in SparseSolvers, `LPStatus`, `SolutionStatus`,
`GraphAlgorithmStatus`) are permitted only when they carry information `AlgorithmStatus` cannot
(e.g. `Infinite` solutions, `Unbounded` LP) **and** the owning header provides
`AlgorithmStatus ToAlgorithmStatus(DomainStatus)` so callers can treat all results uniformly
(task wewi.3.4). New code uses `AlgorithmStatus` directly.

**Public containers in signatures:** algorithm layers take and return `Vector<T>` / `Matrix<T>`,
not `std::vector<T>` (task wewi.3.4 migrates the Krylov solvers).

### Optimization and Krylov migration

The canonical optimization results are `Minimization1DResult`, `MultidimOptimizationResult`
and `LinearProgrammingResult`; the old type names remain deprecated aliases. All inherit
`IterativeResultBase`. Legacy optimization `iterations` mirrors `iterations_used`; a function
evaluation count is never inferred from an iteration count.
Nelder–Mead now reports actual simplex iterations in `iterations`/`iterations_used` and
its separately counted objective calls in `function_evaluations` (the former `iterations`
field incorrectly held objective calls).

LP `.status` is now `AlgorithmStatus`. Use `.lp_status` for `Optimal`, `Infeasible`,
`Unbounded`, etc.; `ToAlgorithmStatus(LPStatus)` supplies the common classification.
Phase-one budget exhaustion is `MaxIterationsExceeded`, not a proof of infeasibility.
Every simplex outcome retains the current original-variable iterate and its measured primal
constraint violation; an uncomputed basis condition number is `NoEstimate`.

Krylov APIs take `Vector<T>`, return `KrylovResult<T>`, and expose `.converged` (a boolean),
`.iterations_used` and `.error_message`. `SolverResult<T>` and `SolverStatus` remain deprecated
type aliases, but callers must migrate the old status enumerators and result member spellings.
`residualNorm` and `achieved_tolerance` are the absolute norm of **the true residual of the
returned iterate**, including iteration-limit outcomes and small nonzero right-hand sides.
Success requires that residual to satisfy the configured criterion, not just a recursive
estimate; BiCGSTAB refreshes and restarts when rounding separates those estimates.

### I/O migration

Serializer, persistence and visualizer functions return `[[nodiscard]] IOResult`.
`SerializeResult`, `VisualizerResult` and `SerializeError` remain deprecated type aliases;
use `IOResult`, `IOError`, `.message` and `.filename` instead of the visualizer-specific
`.errorMessage` and `.dataFilePath`. Launch diagnostics retain `exitCode`, `timedOut` and
`dataOnly`; data-only success is not a claim that a process launched.

Dataset loading offers `IOResult` overloads with a `Dataset&` output parameter: `LoadCSV`,
`LoadCSVSafe`, `LoadTSVSafe`, `LoadJSON`, `LoadJSONSafe`, and auto-format `LoadData`.
Output datasets are replaced only on success. Deprecated no-output-parameter safe wrappers
retain the old payload-bearing `LoadResult` interface for migration; their diagnostics come
from the same `IOResult`, not a separate error model.

Byte offsets are zero-based; line offsets are one-based physical lines, including blank
lines. `offset == -1` and `offset_kind == Unknown` explicitly mean no position is available
(e.g. open failure, process failure, or a nonseekable stream). Do not invent an offset.
JSON syntax errors report the parser byte, CSV schema errors the actual source line,
and failed file reads retain the best available stream position. File-save success includes
successful flush and close, not just buffered insertion.

```cpp
std::stringstream storage;
IOResult saved = Persistence::SaveJson(storage, Vector<Real>{MML::Real(1.0), MML::Real(2.0)});
Vector<Real> restored;
IOResult loaded = Persistence::LoadJson(storage, restored);
if (!saved || !loaded) { /* inspect .error, .message, .filename, .offset */ }
```

This example runs in `MML_DocsApp --result-policy`.

### Exception classification

Checked `Vector::at` and `VectorN::at` access throws `VectorAccessBoundsError`, a subtype of
`IndexError` and `std::out_of_range`, with the attempted index and container size.
Size mismatches still throw `VectorDimensionError`.
Cholesky rejects non-Hermitian values with `DomainError`; nonsquare shapes remain
`MatrixDimensionError`. Exception string constructors are explicit.

Shared core exception definitions are centralized in `mml/MMLExceptions.h`:

| Exception | Owning header | Standard base |
|-----------|---------------|---------------|
| `MML::FourierError` | `mml/MMLExceptions.h` | `std::invalid_argument` |
| `MML::StatisticsError` | `mml/MMLExceptions.h` | `std::runtime_error` |
| `MML::CurveFittingError` | `mml/MMLExceptions.h` | `std::invalid_argument` |

All three retain explicit string constructors and `MMLException` inheritance/message handling.
Include the central header directly when catching or constructing these types without the
corresponding algorithms. The former separate domain-exception headers have been deleted
without forwarding files. The full single-header distribution includes each definition once.

Core still implements FFT, basic statistics and curve fitting, so these shared types cannot
be defined exclusively in packages without an invalid core-to-package dependency. Native
MML-Packages code consumes these exact canonical types; it does not redefine, wrap or alias
them. Package extensions are free to add genuinely package-only exceptions in their own
headers. Centralization does not alter existing thrown classes.

```cpp
#include <mml/MMLExceptions.h>

try {
    throw MML::StatisticsError("same type in core and packages");
} catch (const MML::MMLException& error) {
    std::cout << error.message() << '\n';
}
```

Runnable counterpart: `Docs_Demo_ExceptionOwnership()` in
[docs_demo_result_policy.cpp](../src/docs_demos/core/docs_demo_result_policy.cpp),
selected with `MML_DocsApp --exception-ownership`.

## 6. Worked examples

```cpp
// Iterative: Success only after the check (RootFinding's Detailed form is the config overload)
auto r = RootFinding::FindRootBrentDetailed(f, 0.0, 2.0, cfg);
if (!r) { /* r.status, r.error_message, r.root is the best bracket midpoint */ }

// One-shot fixed rule: Completed, no estimate
auto g = IntegrateGauss10Detailed(f, a, b);
assert(g.status == AlgorithmStatus::Completed);
assert(!g.converged && !HasEstimate(g.error_estimate));

// Adaptive rule: Success with a real estimate
auto t = IntegrateTrapDetailed(f, a, b);
if (t.converged) use(t.value, t.error_estimate);

// Plain call: the value, or ConvergenceError - never a silently non-converged number
Real I = IntegrateSimpson(f, a, b);
// Explicit opt-out when you hold a result and want just the number
Real J = IntegrateSimpsonDetailed(f, a, b).value;     // or static_cast<Real>(...)
```

## 7. Checklist for reviewers

1. Is `Success` set only after a comparison against a tolerance? If the method has no
   tolerance, is it `Completed`?
2. Is every estimate field either computed from the **returned** payload or `NoEstimate`?
3. Does the plain overload throw when the Detailed result is not usable?
4. Is every conversion operator `explicit`; is the result `[[nodiscard]]`?
5. Does the layer return the right type (§5)? Does a domain enum have `ToAlgorithmStatus`?
6. Is there a test asserting the status **and** the estimate honesty for the degenerate path
   (zero integrand, already-converged initial guess, singular input)?

---

## Enforcement

- **Audit task wewi.3.2** sweeps every assignment to `converged`, `residual`, `error`, `status`
  in `mml/core` and `mml/algorithms` against §1–§3.
- **Tasks wewi.3.3–3.6** apply §4–§5 (explicit conversions and unified integration results;
  optimization vocabularies and Krylov signatures; `IOResult`; exception classes).
- Contract tests: `tests/base/result_policy_tests.cpp` (vocabulary) and the per-family
  adversarial tests under `tests/` (populators).
- When a case is not covered here, extend this page in the same change.
