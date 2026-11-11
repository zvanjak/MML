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

*Which* tolerance a comparison may use, where the number lives (`mml/MMLPrecision.h` only) and the
scalar helpers (`IsNearlyEqual`, `IsNearlyZero`, `IsNegligibleRelativeTo`) are specified in
**[TOLERANCE_POLICY.md](TOLERANCE_POLICY.md)**.

## 3. Normalization

`Normalized()` **throws** (`VectorDimensionError` or
`GeometryError` as appropriate) when the norm is near zero
(`< Precision::RoundoffZeroThreshold`, 100 ulp of unit scale — TOLERANCE_POLICY §2). Silent return of a zero or garbage vector is
never acceptable.

- Callers that legitimately want a fallback use the explicit variant
  `NormalizedOrZero()` — the name carries the contract.
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

**What each field may claim** — `Success` only after a checked criterion, `Completed` for
methods without one, `NoEstimate` instead of a literal `0.0` for estimates not computed,
which result type each layer returns — is specified in **[RESULT_POLICY.md](RESULT_POLICY.md)**.
That page is normative for every populator of these fields.

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

**Convenience conversion** — a Result may expose the payload through a named member (`value`,
`root`, ...) or an **`explicit`** `operator PrimaryType()`. Implicit conversion is not allowed:
it silently discards `converged`/`status` (analysis v4, CORE C6). Result types are
`[[nodiscard]]`. See [RESULT_POLICY.md](RESULT_POLICY.md) §4.

**Checklist for a new algorithm:** Config struct · Result struct on `IterativeResultBase` ·
config-based core · simple overload · Doxygen on all fields · tests for config parameters and
convergence reporting · result honesty per [RESULT_POLICY.md](RESULT_POLICY.md) §7 · doc snippet
backed by a `src/docs_demos/` demo (AGENTS.md contract).

## 11. Naming convention

Private/protected members may use the leading-underscore convention, but the next
character must not be uppercase. `_coefA`, `_radius`, `_qr`, and `_sharpness` are
valid class-member spellings; `_A`, `_R`, `_QR`, and `_N` are implementation-reserved
at every scope. Public mathematical accessors such as plane `A()`, `B()`, `C()`,
and `D()` remain valid and unchanged. Compiler/platform predefined macros and
underscore-uppercase substrings inside equations or filename strings are not
member declarations and must not be mechanically renamed.

*Adopted 2026-10-13 (epic MinimalMathLibrary-wewi.5, decision task 5.1), from a census of the
library: algorithm entry points are PascalCase 633 : 70, interface virtuals lowerCamel ~95 %,
`Is*` predicates 705 : 249, MML-Packages-Private calls `size()` 673× and `rows/cols/clear` ~70×.*

**Principle: casing encodes the abstraction level, not the class.** lowercase / lowerCamel is
*plumbing* (STL protocol, interface callbacks, accessors); PascalCase is *mathematics*
(operations, algorithms, mathematical predicates). A reader must be able to tell which one a
name is from its spelling, and the guard script (5.7) must be able to check it from *where the
member is declared*.

| Level | Rule | Examples (✓ conforming) |
|---|---|---|
| Types, namespaces | `PascalCase` | `MatrixSym`, `MML::Systems`, `Integration::GKResult` |
| Config / result fields | `snake_case` (§10) | `error_estimate`, `max_iterations`, `converged` |
| STL-shaped protocol on containers | **exactly as the STL spells it** | `size() at() begin() end() empty() data() clear() push_back() resize() rows() cols() to_string()` |
| Mathematical operations on containers / math objects | `PascalCase` | `Trace() Transpose() Invert() Inverse() Norm1() NormL2() NormFrobenius() Det() Print() IsEqualTo() Identity() Normalized() Conjugate() Dot()` |
| Algorithms (free / static entry points) | `Verb[Object][Method]` `PascalCase`; `XDetailed` returns the result; `XCore` is internal | `IntegrateTrap`, `IntegrateTrapDetailed`, `FindRootBrent`, `SolveBDF2`, `NDer4Partial`, `CalcHessian` |
| Interface vocabulary (`mml/interfaces` virtuals) and their implementers | `lowerCamel` verbs / accessors | `derivs() jacobian() derivsAndJacobian() getDim() getNumParam() getParam(i) setParams(v) contains() isAutonomous() transf() calcStep() doStep() contract() isContravar()` |
| Accessors on concrete classes | `getX()` / `setX()` lowerCamel | `getRadius()`, `getLowerBandwidth()`, `getRotationAngle()`, `setOmega()` |
| Factories / conversions | `PascalCase` nouns (`From…`, `To…`, `Identity`, `Constant`) | `Quaternion::FromAxisAngle`, `MatrixSym::FromLower`, `Matrix::Identity(n)`, `ToRotationMatrix()`, `ToFullMatrix()` |
| Predicates about mathematical / geometric properties | `PascalCase` `Is* Has* Contains Intersects` | `IsSquare() IsZero() IsSymmetric() IsNull() IsInside(p) Contains(p) Intersects(l)` — but lowerCamel when declared on an interface, and STL `empty()` |

**Mutating vs value-returning pairs.** The *verb* mutates in place and returns `void`; the
*noun / participle* returns a new object and is `[[nodiscard]]`:

```cpp
MatrixAlg::Invert(m); Matrix mi = MatrixAlg::Inverse(m); // Matrix inversion is algorithm-only.
m.Transpose();  Matrix mt = m.Transposed();     // not transpose() / GetTranspose()
q.Normalize();  Quaternion u = q.Normalized();  // VectorN / Quaternion already do this
```

Pairs that differ only in case (`Transpose()` / `transpose()`) are **banned**: they read the
same and one of them silently does the other thing.

**Reset vocabulary.** `clear()` = remove all elements (dynamic containers only, STL meaning).
Value operations are `SetZero()`, `SetIdentity()`, `Fill(v)` — never `clear()` on a fixed-size
type, never `MakeUnitMatrix()` / `SetToZero()` / `SetToValue()`.

**Dimension vocabulary.** 2-D containers: `rows()` / `cols()` on every matrix type (also square
ones — no `Dim()` / `GetDimension()`); 1-D: `size()`; N-D: `dim(i)` / `dim1() dim2() dim3()`;
interfaces: `getDim()`.

**Transition rule.** Every rename keeps the old spelling for one minor release as
`[[deprecated("use X()")]] inline` forwarder (same semantics). Docs, demos, tests and the single
header move to the new name in the same change; the guard script allow-lists deprecated aliases;
sibling repos (MML-Packages, book, MPL) are told which aliases they hit. Renames that would
*silently* change behaviour (a value-returning name becoming in-place or vice versa) are not
allowed — the old name must become a compile error or keep its meaning.

**Known renames this policy implies** (tasks 5.2–5.6; **5.2 matrix family, 5.3
vector / quaternion / Matrix3D / geometry, 5.4 interface vocabulary, 5.5 derivation entry
points and 5.6 honest algorithm names done 2026-10-14**, all old
spellings kept as `[[deprecated]]` aliases except the subsequently removed matrix-inversion
members): `Matrix::trace/transpose → Trace/Transposed`,
`TriDiag::GetTranspose → Transposed`, `BandDiag::Transpose (copy) → Transposed`,
`VectorN::clear → SetZero`, `MatrixSym::SetToZero/SetToValue → SetZero/Fill`,
`MakeUnitMatrix → SetIdentity`, `Resize → resize`, `MatrixSym::Dim/BandDiag::GetDimension →
rows`, `Matrix::isZero/isSquare → IsZero/IsSquare` (`isEmpty → empty`), `SparseMatrixCSR::norm*
→ Norm*`, `IInterval::GetEquidistantCovering → getEquidistantCovering`,
`ITensor::Contract/IsContravar/NumCovar → contract/isContravar/numCovar`, `Matrix3D::At → at`,
geometry `GetRadius/GetSide/GetHeight → getRadius/getSide/getHeight`, core `calcHessian /
calcJacobian* → CalcHessian / CalcJacobian*`.

Matrix inversion is an explicit ownership removal, not a compatibility rename:
`Matrix`, `MatrixNM`, `MatrixSym`, and `TridiagonalMatrix` have no inversion members
or deprecated inversion aliases. Use `MatrixAlg::Inverse` / `MatrixAlg::Invert`.

**BLAS-named kernels keep their BLAS spelling** (`gemv`, `axpy`, `dot`, `norm2` on raw
`std::vector` in `MML::SparseMatrix`) — like `size()` for the STL, the lowercase name *is* the
convention being followed. Everything else on the sparse types follows the table above
(`Multiply`, `SolveLower`, `IsSymmetric`, `Norm1`, `Laplacian2D`, …).

**MML-Packages-Private cross-check (2026-10-13):** it uses `size rows cols clear Resize Trace
Contains NormL2 Print`; all keep their meaning, only `Resize` will warn through its alias.

---

## Enforcement

1. **New code:** review against this page; regression tests assert contract behavior
   (see `tests/base/bugfix_epic_i5u3_tests.cpp` for the pattern).
2. **Existing code:** conformance sweeps per layer — base (i5u3.2.2), core + algorithms
   (i5u3.2.3) — fix deviations and add contract tests.
3. **Ambiguity:** when a case is not covered here, extend this document in the same
   change that introduces the behavior, and keep it one page.

Matrix inversion now belongs to algorithms, not base containers; see
[the inversion and formatting migration](algorithms/MatrixInversionAndFormatting.md).
