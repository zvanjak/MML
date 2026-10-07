# MML Contract Policy

**Status:** Adopted (epic MinimalMathLibrary-i5u3.2) · **Scope:** every public API in `mml/`
**Audience:** contributors. The user-facing distillation of these guarantees lives in
[Fundamentals.md](Fundamentals.md) §5 — keep the two in sync.

This is the library-wide behavioral contract. Same-named operations must behave the same
way on every type. New code must conform; existing code is brought into conformance by
the layer-by-layer sweeps (i5u3.2.2, i5u3.2.3). Where current code deviates, the policy
wins — the deviation is a bug.

---

## 1. Equality semantics

`operator==` / `operator!=` **never throw**. Comparing objects of different dimensions
returns `false` — a 3-vector simply *is not equal* to a 4-vector.

- Exact element-wise comparison; no tolerance.
- Tolerance-based comparison is a named method: `IsEqualTo(other, eps)` (existing name, keep it).
- *Known deviation:* `Vector::operator==` throws `VectorDimensionError` on size mismatch
  while `Matrix::operator==` returns `false` — Vector must be changed to return `false`.

## 2. Zero queries: exact vs. tolerance

Naming states the semantics:

| Name | Meaning |
|------|---------|
| `IsZero()` | exact: every component `== 0` |
| `IsNearZero(eps = Precision::DefaultTolerance)` | tolerance: norm (or max-abs) `<= eps` |

A method named `isZero`/`IsZero` that internally applies a tolerance is a policy
violation — rename it to `IsNearZero` or make it exact.

## 3. Normalization

`Normalized()` / `GetAsUnitVector()` **throw** (`VectorDimensionError` or
`GeometryError` as appropriate) when the norm is near zero
(`< Precision::DefaultTolerance`). Silent return of a zero or garbage vector is
never acceptable.

- Callers that legitimately want a fallback use the explicit variant
  `GetAsUnitVectorOrZero()` — the name carries the contract.
- Constructors taking a direction (lines, planes, rays) validate via the throwing path.

## 4. Checked vs. unchecked access

Three tiers, consistently:

| Access | Contract |
|--------|----------|
| `operator[]`, `operator()` | **unchecked** — fast path; `assert` in Debug only |
| `at(...)` | **checked** — throws `VectorAccessBoundsError` / `MatrixAccessBoundsError` / `TensorIndexError` |
| high-level operations (solvers, transforms, algorithms) | **validate their inputs** at entry and throw MML exceptions |

No third pattern (silent clamp, silent resize, return-default) is permitted in element access.

## 5. Null callbacks

A callable slot (`std::function`, function pointer, virtual hook) that is required but
not provided **throws `NotImplementedError`** at the point of use. Never a silent no-op,
never a default result. Optional callbacks must be documented as optional and tested
for emptiness with a clearly named guard.

## 6. Exceptions at public boundaries

Every public entry point reports failure through the **MML exception taxonomy**
(`MMLExceptions.h`; all types derive from `MMLException` and a matching `std::` base).

- No raw `std::runtime_error`/`std::logic_error` thrown from `mml/` code.
- No `catch (...)` that swallows errors; catch specific types, add context, rethrow or
  convert to a result object.
- Bool-return failure reporting is allowed only in result objects (§8) and simple
  predicates — never as the sole error channel of an algorithm.

## 7. No `std::cout` in library code

Solvers and algorithms never write to `std::cout`/`std::cerr`. Diagnostics go through:

- an `std::ostream&` parameter (default-able to a null stream), or
- a verbosity callback on the config object.

`verbose` flags on config objects must route through one of those two, not hardcode
`std::cout`. Printing helpers (`Print(...)`, visualizers, serializers) take their
stream explicitly.

## 8. Result objects for iterative algorithms

Every iterative algorithm (root finding, optimization, ODE/DAE, eigen, iterative
linear solvers) returns a result object that reports, at minimum:

```
converged      (bool)
iterations     (int)
residual / achieved_tolerance
status         (AlgorithmStatus enum)
error_message  (string, empty on success)
```

Throwing is reserved for *invalid input*; *failure to converge* is data, reported in
the result. Convenience wrappers that throw on non-convergence must be thin shims over
the result-returning core (pattern already used by `RootFinding::FindRootNewton`).

## 9. Ownership vocabulary

Type and method names state ownership:

| Suffix / prefix | Meaning |
|-----------------|---------|
| `...View` / `...Ref` | non-owning; caller guarantees lifetime |
| plain name / `...Owned` | owning copy |
| `Get...()` returning by value | independent copy, safe to keep |
| `Get...()` returning reference | tied to the source object's lifetime — document it |

A method must not return a reference into internal storage without the name or doc
making the lifetime coupling explicit.

## 10. Algorithm API shape (Config + Result)

Every iterative algorithm follows the Config + Result pattern. Reference implementation:
`mml/algorithms/RootFinding.h`; shared base types: `mml/base/AlgorithmTypes.h`.

**Config structs** — named `{Family}Config`, `snake_case` fields, every field defaulted:
`tolerance` (Real), `max_iterations` (int), `verbose` (bool, positive phrasing) + any
algorithm-specific fields. `verbose` routes through §7, never `std::cout`.

**Result structs** — named `{Family}Result`, extend `IterativeResultBase` (which supplies
`converged`, `iterations_used`, `achieved_tolerance`, `status` (`AlgorithmStatus`),
`error_message`, `algorithm_name`, `elapsed_time_ms`, `function_evaluations`) and add the
payload named for what it is (`root`, `eigenvalues`, `minimizer`, ...). Error messages are
specific and actionable ("Failed to converge after 100 iterations (residual: 1.5e-4)");
empty on success (§6, §8).

**Signatures** — config-based core + thin simple overload that delegates:

```cpp
FamilyResult Algorithm(const InputType& input, const FamilyConfig& config);
FamilyResult Algorithm(const InputType& input, Real tol = 1e-10, int max_iter = 100); // delegates
```

Inputs by `const&` (containers/functions/configs) or value (scalars); result returned by
value. Algorithms taking callables also provide a constrained `Function&&` overload
(`RealFunctionCallable`) so raw lambdas work.

**Convenience conversion** — `operator PrimaryType()` on a Result is allowed only when one
value is mathematically obvious, and must carry a `@warning` that diagnostics are discarded.

**Checklist for a new algorithm:** Config struct · Result struct on `IterativeResultBase` ·
config-based core · simple overload · Doxygen on all fields · tests for config parameters and
convergence reporting · doc snippet backed by a `src/docs_demos/` demo (AGENTS.md contract).

---

## Enforcement

1. **New code:** review against this page; regression tests assert contract behavior
   (see `tests/base/bugfix_epic_i5u3_tests.cpp` for the pattern).
2. **Existing code:** conformance sweeps per layer — base (i5u3.2.2), core + algorithms
   (i5u3.2.3) — fix deviations and add contract tests.
3. **Ambiguity:** when a case is not covered here, extend this document in the same
   change that introduces the behavior, and keep it one page.
