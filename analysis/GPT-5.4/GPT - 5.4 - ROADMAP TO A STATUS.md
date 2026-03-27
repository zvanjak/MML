# GPT-5.4 Roadmap To A Status

## Goal

This document turns the analysis findings into a concrete improvement program whose purpose is to move Minimal Mathematical Library from its current overall grade of **B** toward **A-level** quality.

For MML, an `A` status does not require more features. It requires stronger correctness discipline, more uniform contracts, better edge-case handling, and tighter packaging/runtime consistency.

## What `A` Status Means For MML

MML should be considered `A`-grade when the following are true:

1. Numerical decision rules are consistent across the library.
2. Release builds are as contract-safe as debug builds for correctness-critical paths.
3. Public APIs enforce their advertised guarantees, especially around invalid inputs and non-convergence.
4. External boundary layers such as serialization and data loading fail closed on malformed input.
5. Single-header and modular builds behave the same way for core type and precision semantics.
6. Pathological and near-degenerate cases are covered by targeted regression tests.

## Guiding Principles

- Fix correctness before convenience.
- Standardize policy before adding more algorithms.
- Remove silent failure paths.
- Prefer explicit contracts over caller discipline.
- Treat numerical edge cases as first-class behavior, not exceptional afterthoughts.

## Program Structure

This roadmap is organized into five workstreams:

1. Numerical policy unification
2. Contract and safety hardening
3. Algorithm correctness repairs
4. Tooling and I/O reliability
5. Test, packaging, and documentation hardening

## Recommended Execution Order

### Phase 0. Baseline and tracking

Duration: short

Purpose:

Create the structure needed to execute the rest of the roadmap without losing consistency.

### Phase 1. High-risk correctness fixes

Duration: short to medium

Purpose:

Fix the most severe correctness issues that can silently produce wrong results or undefined behavior.

### Phase 2. Policy unification across subsystems

Duration: medium

Purpose:

Replace local heuristics and mixed failure semantics with shared library-wide rules.

### Phase 3. Boundary hardening

Duration: medium

Purpose:

Make serialization, parsing, and distribution behavior trustworthy.

### Phase 4. Regression-proofing

Duration: ongoing

Purpose:

Add adversarial tests and documentation so the improvements remain stable.

## Workstream 1. Numerical Policy Unification

### Task 1. Create a library-wide tolerance policy document and helper set

Problem:

Tolerance decisions are currently split across centralized helpers and hard-coded literals like `1e-10`, `1e-8`, and `1e-6`.

Target areas:

- `mml/MMLPrecision.h`
- `mml/MMLBase.h`
- `mml/core/**`
- `mml/algorithms/**`
- `mml/systems/**`

Concrete work:

- define standard categories for tolerances: zero tests, convergence tests, singularity tests, classification thresholds, and default algorithm tolerances
- expose those categories through reusable helpers or policy constants derived from `PrecisionValues<Real>`
- document which categories are absolute, relative, or mixed

Deliverable:

- one canonical tolerance policy used by all new and migrated code

Acceptance criteria:

- no new hard-coded numerical thresholds in solver logic without explicit justification
- a contributor can find the canonical tolerance source in one place

### Task 2. Audit and replace fragile exact-zero comparisons in sensitive math code

Problem:

Exact equality tests remain in vector, geometry, interpolation, and related code paths where near-zero behavior matters.

Target areas:

- `mml/base/**`
- geometry and interpolation helpers identified in the base analysis
- singularity-sensitive code in `mml/core/**`

Concrete work:

- locate all `== 0`, `!= 0`, or equivalent exact comparisons on floating-point values
- classify them into safe exact comparisons and unsafe numerical comparisons
- replace unsafe cases with shared helpers such as `isNearlyZero` or a stronger replacement if current helper semantics are insufficient

Deliverable:

- reviewed inventory of exact comparisons with disposition for each case

Acceptance criteria:

- all floating-point exact comparisons in numerically sensitive paths are either removed or explicitly documented as intentional

### Task 3. Remove local threshold drift from analyzers and classification code

Problem:

Systems and algorithmic classifiers use local heuristic thresholds that are not consistently tied to the library precision model.

Target areas:

- `mml/systems/DynamicalSystemAnalyzers.h`
- relevant files in `mml/algorithms/**`

Concrete work:

- inventory every literal threshold used for classification and branching
- group them into policy-backed categories
- make the chosen policy explicit in code comments or documentation where the threshold is domain-specific rather than generic

Acceptance criteria:

- threshold-bearing branches in systems analysis either use shared policy or explicitly documented domain constants

## Workstream 2. Contract And Safety Hardening

### Task 4. Eliminate `assert`-only correctness guards from foundational runtime paths

Problem:

Some foundational access checks exist only in debug builds, which creates release/debug correctness divergence.

Target areas:

- tensor indexing and similar shape-sensitive access paths in `mml/base/**`

Concrete work:

- identify all places where an invalid index or shape can silently produce undefined behavior in release mode
- replace `assert`-only protection with runtime validation in public and correctness-critical entry points
- preserve fast internal unchecked paths only when they are clearly internal and intentionally named

Deliverable:

- hardened runtime checks for public container access and shape-sensitive operations

Acceptance criteria:

- invalid public indexing and dimension misuse fails deterministically in release builds

### Task 5. Audit non-owning views and lifetime-sensitive APIs

Problem:

Some APIs expose lifetime hazards, including view-like constructs whose safety is not obvious from the type.

Target areas:

- `mml/base/Matrix/**`
- any view or borrowed-storage abstractions identified in the base analysis

Concrete work:

- inventory non-owning view types and constructors
- document ownership and lifetime expectations in the API
- where needed, rename unsafe entry points or add safer alternatives

Acceptance criteria:

- public view APIs clearly communicate ownership model and invalidation rules

### Task 6. Standardize failure semantics across core entry points

Problem:

MML currently mixes exceptions, booleans, `assert`, and stderr logging across subsystems.

Target areas:

- `mml/core/**`
- `mml/tools/**`
- public-facing base utilities

Concrete work:

- define preferred failure behavior for each API category:
  - invalid arguments
  - non-convergence
  - unsupported format/input
  - internal invariant violation
- refactor APIs within each category to match the chosen model

Deliverable:

- short design note describing failure semantics by subsystem

Acceptance criteria:

- same class of failure is reported consistently within each subsystem
- new code does not introduce ad hoc stderr-based error handling in core logic

## Workstream 3. Algorithm Correctness Repairs

### Task 7. Fix final error reporting in `IntegrateTrap` and `IntegrateSimpson`

Problem:

The core analysis identified incorrect final error reporting behavior in non-convergent paths.

Target areas:

- `mml/core/Integration/Integration1D.h`

Concrete work:

- inspect convergence loop and final result construction
- ensure reported error estimate corresponds to the final computed state rather than a stale or mismatched value
- add tests for both convergent and non-convergent cases

Acceptance criteria:

- result object accurately reflects the last computed estimate and error state
- regression tests catch the previous mismatch

### Task 8. Fix narrowing in numeric validation helpers

Problem:

`NumericValidation::IsFiniteValue` narrows through `double`, which weakens type-correctness for wider `Real` configurations.

Target areas:

- `mml/core/NumericValidation.h`

Concrete work:

- make finiteness checks type-correct for `float`, `double`, and `long double`
- review nearby validation helpers for the same issue

Acceptance criteria:

- validation helpers preserve the effective precision/type semantics of `Real`

### Task 9. Enforce adaptive ODE configuration contracts

Problem:

Adaptive ODE code appears not to fully honor fields such as `max_steps` and `max_step_size`, and some overloads may mishandle invalid `outputInterval`.

Target areas:

- `mml/algorithms/ODESolvers/ODESolverAdaptive.h`

Concrete work:

- audit all integration entry points against the config structure
- ensure every documented config field is enforced or removed from the public contract
- validate `outputInterval` and related derived loop parameters at API boundaries

Acceptance criteria:

- documented config fields have actual behavioral effect
- invalid output interval values fail early and clearly

### Task 10. Harden ODE event detection beyond sign-change-only behavior

Problem:

Current event detection can miss tangent events and possibly endpoint events because it relies too heavily on sign changes.

Target areas:

- `mml/algorithms/ODESolvers/ODESolverEventDetection.h`

Concrete work:

- review event bracketing logic
- add handling for zero-at-endpoint and near-tangent cases
- make event semantics explicit in documentation

Acceptance criteria:

- tests cover sign-crossing, tangent touch, and endpoint event cases
- users can predict what counts as a detected event

### Task 11. Guard Lyapunov analysis against zero or near-zero normalization failure

Problem:

Lyapunov Gram-Schmidt normalization can divide by zero or near-zero vectors.

Target areas:

- `mml/systems/DynamicalSystemAnalyzers.h`

Concrete work:

- add threshold-aware safeguards before normalization
- define fallback behavior for degenerate orthogonalization steps
- report insufficiently conditioned trajectories clearly

Acceptance criteria:

- no divide-by-zero or unstable normalization on degenerate test cases
- output status clearly distinguishes valid spectrum estimates from failed runs

### Task 12. Revisit Jacobian update strategy in Radau IIA DAE solver

Problem:

The analysis indicates frozen Jacobians through Newton iterations, which may be acceptable for some workloads but is a correctness and robustness tradeoff that should be explicit.

Target areas:

- `mml/algorithms/DAESolvers/DAERadauIIA.h`

Concrete work:

- evaluate whether Jacobian refresh should be configurable or automatic under convergence failure
- expose the strategy explicitly in configuration or documentation
- add hard cases where stale Jacobians degrade convergence

Acceptance criteria:

- Jacobian reuse policy is explicit and tested against difficult nonlinear cases

## Workstream 4. Tooling And I/O Reliability

### Task 13. Make matrix and vector loading fail closed on incomplete reads

Problem:

Some load paths currently return success without verifying full payload read completion.

Target areas:

- `mml/tools/serializer/MatrixIO.h`
- related loader modules in `mml/tools/serializer/**`

Concrete work:

- verify stream state after every dimension and payload read
- verify exact byte counts for binary loads
- reset or leave outputs untouched on failure, depending on chosen API contract

Acceptance criteria:

- truncated files never return success
- partial payload reads are reported deterministically

### Task 14. Standardize serializer result contracts

Problem:

Serializer helpers mix `SerializeResult`, boolean returns, void return types, and stderr-based diagnostics.

Target areas:

- `mml/tools/serializer/SerializerBase.h`
- serializer modules under `mml/tools/serializer/**`

Concrete work:

- define one primary result model for serializer operations
- migrate helper signatures to that model
- remove direct stderr printing from reusable library code unless explicitly configured as optional diagnostics

Acceptance criteria:

- serializer modules present one coherent failure-reporting model to callers

### Task 15. Stabilize text serialization formatting

Problem:

Text serialization is vulnerable to locale differences and inconsistent precision formatting.

Target areas:

- `mml/tools/serializer/**`

Concrete work:

- set deterministic locale behavior for serialization and parsing
- define numeric precision policy for `Real`
- add round-trip tests using values that are sensitive to precision loss

Acceptance criteria:

- text serialization round-trips reliably under controlled locale settings
- serialization behavior is documented and reproducible

### Task 16. Make binary serialization portability boundaries explicit and enforceable

Problem:

Binary I/O writes raw bytes of templated types but does not strongly constrain portability assumptions.

Target areas:

- `mml/tools/serializer/MatrixIO.h`
- any similar binary serializers

Concrete work:

- constrain binary serialization to supported type categories
- document endianness and ABI assumptions
- decide whether MML binary formats are local-only or cross-platform portable

Acceptance criteria:

- users cannot mistake local raw dumps for fully portable binary interchange formats

### Task 17. Harden the flat JSON loader against malformed inputs

Problem:

The custom parser is intentionally simple but still too fragile at the malformed-input boundary.

Target areas:

- `mml/tools/data_loader/DataLoaderJSON.h`

Concrete work:

- strengthen bounds checks in `ParseNumber` and neighboring parse routines
- improve diagnostics for unsupported nested structures and malformed tokens
- decide whether to keep the parser intentionally limited or replace it with a better-scoped implementation

Acceptance criteria:

- malformed JSON fails cleanly without out-of-range parsing behavior
- supported JSON shape is clearly documented

## Workstream 5. Test, Packaging, And Documentation Hardening

### Task 18. Add a dedicated adversarial correctness test suite

Problem:

Current quality appears stronger on nominal paths than on pathological paths.

Target areas:

- `tests/**`
- possibly `test_beds/**` for exploratory edge cases

Concrete work:

Create focused tests for:

- near-singular linear systems
- zero and near-zero vectors in geometry and transforms
- non-convergent integration and root-finding paths
- ODE event corner cases
- malformed or truncated serializer inputs
- long-double or alternative `Real` configurations where supported

Acceptance criteria:

- every issue class identified in the analysis has at least one regression test

### Task 19. Add precision-configuration parity tests

Problem:

The modular and single-header distributions do not fully align in precision semantics.

Target areas:

- `mml/single_header/MML.h`
- build/test infrastructure for modular and single-header consumption

Concrete work:

- define what parity means for type aliases, precision constants, and helper behavior
- add tests that compile and validate both modular and single-header usage
- either align behavior or document the delta as intentional

Acceptance criteria:

- single-header and modular consumption either match or differ only in explicitly documented ways

### Task 20. Create subsystem contract documentation for failure and convergence behavior

Problem:

Some APIs appear more rigorous than they actually are because contract boundaries are under-documented.

Target areas:

- `docs/ERROR_HANDLING.md`
- `docs/THREADING.md`
- integration/ODE/root-finding docs
- serialization docs under `docs/` and `mml/tools/README_Serialization.md`

Concrete work:

Document:

- how invalid arguments are reported
- what non-convergence means in result objects
- what guarantees load/save functions actually provide
- what event detection can and cannot detect
- what binary portability is supported

Acceptance criteria:

- users can understand failure behavior without reading implementation details

### Task 21. Create a library-wide correctness checklist for new code reviews

Problem:

Without a review checklist, the same inconsistency patterns will return.

Concrete work:

Create a short checklist that asks:

- does this code use central tolerance policy?
- does it fail safely in release builds?
- are malformed inputs handled deterministically?
- are non-convergence outcomes explicit?
- are ownership and lifetime rules clear?
- does the single-header path remain consistent?

Acceptance criteria:

- checklist is used for future algorithm and tooling additions

## Suggested Milestones

### Milestone A1. Stop the bleeding

Complete:

- Task 4
- Task 7
- Task 8
- Task 9
- Task 13
- Task 17
- Task 18

Outcome:

Major correctness hazards and silent-failure paths are reduced quickly.

### Milestone A2. Unify policy

Complete:

- Task 1
- Task 2
- Task 3
- Task 6
- Task 14
- Task 15

Outcome:

MML starts behaving like one library instead of several maturity levels coexisting.

### Milestone A3. Harden advanced behavior

Complete:

- Task 10
- Task 11
- Task 12
- Task 16
- Task 19
- Task 20

Outcome:

Advanced solvers and distribution forms become trustworthy under harder workloads.

### Milestone A4. Keep it there

Complete:

- Task 21
- ongoing expansion of Task 18

Outcome:

The library stops regressing into mixed policy and under-specified behavior.

## Recommended First 90 Days

### First 30 days

- finish Task 4
- finish Task 7
- finish Task 8
- finish Task 13
- finish Task 17
- create the first version of Task 18

### Days 31 to 60

- finish Task 1
- finish Task 2
- finish Task 6
- finish Task 9
- finish Task 14
- start Task 19

### Days 61 to 90

- finish Task 10
- finish Task 11
- finish Task 15
- finish Task 20
- create Task 21 checklist and institutionalize it

## Highest-Leverage Starting Points

If only a small number of tasks can begin immediately, start with these:

1. Task 4: remove `assert`-only correctness guards from public runtime paths
2. Task 7: fix integration error-reporting correctness
3. Task 9: enforce adaptive ODE config contracts
4. Task 13: make file loading fail closed
5. Task 18: add adversarial regression tests for every fixed issue

These tasks improve both actual correctness and confidence in the rest of the codebase.

## Final Recommendation

The fastest route to `A` status is to treat MML as a hardening project, not a feature-expansion project.

The library already has enough mathematical breadth. What it needs now is consistency:

- one tolerance policy
- one clear contract philosophy
- one trustworthy boundary behavior model
- one disciplined regression suite for degenerate and adversarial cases

If this roadmap is executed well, MML can realistically move from a strong `B` to a defensible `A-` and then an `A`-level library without major architectural upheaval.
