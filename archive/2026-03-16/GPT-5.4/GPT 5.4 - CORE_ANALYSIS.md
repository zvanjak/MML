# GPT-5.4 Core Analysis

## Phase Definition

- Phase: Core analysis
- Scope: `mml/core/**`
- Approximate size: 51 headers
- Objective: analyze differentiation, integration, linear algebra solvers, coordinate systems, curves, surfaces, fields, orthogonal bases, and the numerical infrastructure that sits between foundational data types and higher-level algorithms
- Output of this phase: this file

## Report Conventions

- Severity tags:
	- `[Critical]` architectural or cross-core inconsistency with broad API impact
	- `[High]` important robustness, maintainability, or performance issue
	- `[Medium]` meaningful improvement area that should be scheduled
	- `[Low]` useful cleanup or refinement
- Priority levels:
	- `P1` do first
	- `P2` do after core behavior is stabilized
	- `P3` targeted enhancement work

## Executive Summary

The core layer is where MML starts to look like a full numerical computing environment rather than a container library. This part of the codebase is strong in capability and mathematical breadth. It covers numerical differentiation, multiple integration families, coordinate transformations, field operations, geometric primitives, and a broad set of linear algebra solvers. The code also shows a real effort toward numerical validation, documentation, and reusable abstraction patterns.

The main issue is not lack of capability. It is inconsistency. The core layer currently mixes newer result/config-based APIs with older exception-driven APIs, and sometimes the same subsystem uses both. That makes the layer more powerful than coherent. The other recurring problem is that some important numerical policies, such as singularity handling and validation of intermediate results, are not yet applied uniformly across the subsystem.

## Top Findings

- `[Critical]` Core APIs do not follow one consistent error and result model.
- `[High]` Validation helpers exist, but intermediate values in numerical loops are not checked consistently.
- `[High]` Singularity policy exists as infrastructure but is not yet applied uniformly across fields, transforms, and related geometry code.
- `[High]` `CoordTransf.h` relies on numerical Jacobians as the default path even for transforms that could expose analytic derivatives.
- `[Medium]` Several core headers are dense enough to become long-term maintenance hotspots.

## What Is Strong

### 1. Core capability coverage is genuinely broad

The core layer is not narrowly focused on one numerical domain. It spans multiple essential pillars:

- `mml/core/Derivation.h`
- `mml/core/Integration.h`
- `mml/core/LinAlgEqSolvers.h`
- `mml/core/CoordTransf.h`
- `mml/core/FieldOperations.h`
- `mml/core/Fields.h`
- `mml/core/Curves.h`
- `mml/core/Surfaces.h`
- `mml/core/OrthogonalBasis.h`

That breadth matters because it lets upper layers share a common mathematical engine rather than re-implementing domain-specific helpers.

### 2. There are clear signs of numerical maturity

Several core headers show careful thought about numerical failure modes.

Positive examples:

- `mml/core/NumericValidation.h` centralizes checks for non-finite values and bad numeric inputs.
- `mml/core/Derivation/DerivationBase.h` defines order-specific differentiation step sizes instead of using a single arbitrary global delta.
- `mml/core/Integration/Integration1D.h` returns `IntegrationResult` objects rather than just a scalar, which is the right design for numerical work.
- `mml/core/SingularityHandling.h` introduces an explicit singularity policy model instead of hardwiring one behavior everywhere.

This is the difference between code that merely computes and code that tries to fail intelligibly.

### 3. Linear algebra solver coverage is one of the strongest parts of the core layer

`mml/core/LinAlgEqSolvers.h` exposes a meaningful solver spectrum:

- direct methods in `LinAlgDirect.h`
- QR in `LinAlgQR.h`
- SVD in `LinAlgSVD.h`
- iterative methods in `LinAlgEqSolvers_iterative.h`

That gives users options across the usual numerical tradeoffs: speed, stability, structure exploitation, and iterative scalability.

### 4. Coordinate transformation contracts are clearly documented

`mml/core/CoordTransf.h` is strong conceptually. It explains:

- Jacobian layout
- covariant vs contravariant basis meaning
- transformation formulas
- coordinate conventions
- numerical Jacobian assumptions

For a mathematically dense subsystem, that kind of explicit contract is valuable and not common enough.

### 5. Umbrella headers are clean and predictable

The umbrella headers I checked are concise and readable:

- `mml/core/Integration.h`
- `mml/core/Derivation.h`
- `mml/core/LinAlgEqSolvers.h`

These are good user-facing entry points because they explain what each family covers without burying the reader in implementation details.

### 6. Domain modeling in `Fields`, `Curves`, and `Surfaces` is expressive

`mml/core/Fields.h` is a good example of domain-oriented API design. It encodes physically meaningful field types and conventions instead of exposing only raw mathematical kernels. That is a strength for users doing scientific or engineering work, because the library is speaking in domain objects rather than only helper functions.

## Severity-Tagged Findings

### [Critical] Core API standardization is incomplete

This is the most important weakness in the core layer.

There is a visible effort toward a standardized configuration/result model in:

- `mml/core/AlgorithmTypes.h`
- `mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h`

But the rest of the core subsystem does not consistently follow that pattern.

Concrete examples:

- `Integration1D.h` returns `IntegrationResult`.
- direct and SVD solvers still lean heavily on exceptions.
- iterative solvers both return `IterativeSolverResult` and also throw on invalid input or singular structure.
- `Fields.h`, `Curves.h`, and `FieldOperations.h` still use ad hoc exception choices in several places.

The result is not just stylistic inconsistency. It affects user ergonomics, test design, and composability. A caller cannot rely on one uniform way to detect failure across the core numerical APIs.

### [High] Validation policy is present, but not consistently enforced through computational loops

`mml/core/NumericValidation.h` exists for a reason, but important core routines do not always use that discipline end-to-end.

In `mml/core/Integration/Integration1D.h`, `IntegrateTrap` and `IntegrateSimpson` compare successive estimates and build `IntegrationResult`, but they do not validate every intermediate estimate for non-finite values before continuing the refinement loop.

That means there is a gap between having a validation subsystem and systematically using it where it matters most.

### [High] Singularity handling is a subsystem, not yet a uniform core policy

`mml/core/SingularityHandling.h` is a good design move, but the core layer still mixes policy-based handling with local hard-coded singularity logic.

Example:

- `Fields.h` contains local checks such as returning zero when a denominator is too small.
- `CoordTransf.h` relies on numerical differentiation and transformation formulas without visibly integrating the singularity policy abstraction.

This creates inconsistent semantics across mathematically similar edge cases.

### [High] `CoordTransf.h` is conceptually strong but computationally expensive by default

The base coordinate transformation API computes Jacobians and basis vectors numerically through repeated calls to `Derivation::NDer4Partial`.

That is a flexible default, but it has costs:

- repeated numerical differentiation work for common transforms
- unnecessary approximation when analytic Jacobians are known
- more opportunity for noise near singular or ill-scaled regions

This is a reasonable fallback strategy, but it should not be the only path in such a central abstraction.

### [Medium] Some files are taking on too much responsibility

The core layer includes several broad, concept-heavy headers, and some of them appear structurally dense enough to become maintenance hotspots.

Most obvious examples from the inspected files:

- `mml/core/FieldOperations.h`
- `mml/core/Curves.h`
- `mml/core/Surfaces.h`

This is not necessarily wrong in a header-only library, but large mathematically mixed headers tend to accumulate duplication, branching policy decisions, and harder review surfaces over time.

### [High] Newer config/result patterns are not yet library-wide defaults

`AlgorithmTypes.h` is promising, but its own comments frame it as a pattern rather than an enforced substrate. That is useful for documentation, but not enough for consistency at the codebase level.

Until the rest of core actually adopts:

- shared status codes
- shared config conventions
- shared diagnostic fields
- shared timing and evaluation-count fields

the benefits of the pattern remain partial.

### [Medium] Some error signaling still bypasses the MML-specific exception taxonomy

The grep results in `mml/core` show several cases of direct `std::invalid_argument` use in places like `Curves.h` and `FieldOperations.h`. That weakens the otherwise strong exception model from the root layer.

If MML already has a structured exception hierarchy, the core layer should prefer it consistently.

### [Medium] Numerical step-size policy is solid but still fairly global and rigid

`Derivation/DerivationBase.h` provides order-specific constants and a scale adjustment through `ScaleStep`. That is good, but it is still a generic policy that may not be sufficient for all function classes or coordinate-system-sensitive problems.

This is not a bug. It is a limitation worth acknowledging. For a library this ambitious, some advanced paths will eventually want more explicit user control over local differentiation strategy.

## Priority Matrix

| Priority | Severity | Area | Why it matters | Recommended action |
|----------|----------|------|----------------|--------------------|
| P1 | Critical | Core API model | Mixed result/exception behavior makes the subsystem harder to use and extend | Standardize how core numerical APIs report success, failure, and diagnostics |
| P1 | High | Runtime validation | Numerical loops can propagate bad values before returning a result | Validate intermediate estimates in integration, solver, and derivation pipelines |
| P1 | High | Singularity behavior | Similar edge cases currently behave differently across modules | Apply a shared singularity policy across fields, transforms, curves, surfaces, and metric-related code |
| P2 | High | Coordinate transforms | Numerical Jacobians are flexible but expensive and noisy for common transforms | Add analytic Jacobian fast paths while preserving numerical fallback behavior |
| P2 | Medium | Header maintainability | Large core headers concentrate complexity and review cost | Split oversized files by operation family or coordinate system while preserving umbrella headers |
| P3 | Medium | Exception consistency and docs | Inconsistent exception types and limited behavior notes increase surprise | Prefer MML exception types and expand cross-module performance and stability notes |

## Recommended Actions

### P1. Finish API standardization across the whole core layer

Recommended change:

- Decide clearly where core APIs should return result objects, where they should throw, and how those two mechanisms should coexist.
- Make `AlgorithmTypes.h` more than a documentation pattern.
- Standardize diagnostic fields across integration and solver families.

This is the single highest-value improvement because it will make the layer easier to use, easier to test, and easier to extend.

### P1. Apply numeric validation inside iterative loops, not just at boundaries

Recommended change:

- Validate intermediate values in integration refinement loops.
- Audit solver and differentiation code for the same pattern.
- Prefer systematic validation helpers over local ad hoc checks.

This closes the gap between having robustness infrastructure and actually benefiting from it.

### P1. Make singularity policy a true cross-core contract

Recommended change:

- Audit `Fields`, `FieldOperations`, `CoordTransf`, `Curves`, `Surfaces`, and metric-related code for local singularity handling.
- Route those behaviors through a shared policy surface where practical.

This is important for predictability, especially for users working across coordinate systems and field operations.

### P2. Add analytic fast paths for coordinate transforms

Recommended change:

- Keep numerical Jacobians as a fallback.
- Add analytic Jacobian support to the transformation interface where exact forms are known.
- Reuse those exact derivatives in basis and vector transformation routines.

This would improve performance and stability simultaneously.

### P2. Break oversized headers into narrower implementation families

Recommended change:

- Split large domain headers by coordinate system, operation family, or representation type.
- Preserve umbrella headers for usability.

That keeps public ergonomics while reducing maintenance density.

### P3. Enforce MML exception usage inside core code

Recommended change:

- Replace direct standard exception use with the MML taxonomy where appropriate.
- Reserve raw standard exceptions for cases where they are explicitly part of the public contract.

This would align the core with the root-layer error design instead of bypassing it.

### P3. Expand explicit performance and stability documentation

Recommended change:

- Document where a routine is exact, approximate, adaptive, or numerically sensitive.
- Make clear which methods rely on numerical differentiation internally.
- Add complexity and conditioning notes to key solver entry points.

The code already contains good comments; this step would make cross-module behavior easier to predict.

## Notable Files In This Phase

- `mml/core/AlgorithmTypes.h`
- `mml/core/NumericValidation.h`
- `mml/core/SingularityHandling.h`
- `mml/core/Derivation.h`
- `mml/core/Derivation/DerivationBase.h`
- `mml/core/Integration.h`
- `mml/core/Integration/Integration1D.h`
- `mml/core/LinAlgEqSolvers.h`
- `mml/core/LinAlgEqSolvers/LinAlgDirect.h`
- `mml/core/LinAlgEqSolvers/LinAlgSVD.h`
- `mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h`
- `mml/core/CoordTransf.h`
- `mml/core/FieldOperations.h`
- `mml/core/Fields.h`
- `mml/core/Curves.h`
- `mml/core/Surfaces.h`

## Bottom Line

The core layer is capable, mathematically serious, and already contains many of the ingredients of a high-quality numerical platform. Its strongest qualities are breadth, documented numerical intent, and a visible push toward reusable infrastructure.

Its main weakness is that the infrastructure is not yet enforced uniformly. The next phase of improvement should not be adding more capabilities. It should be finishing standardization: common result/config conventions, consistent exception policy, broader use of validation utilities, and more disciplined handling of singularities and analytic-vs-numeric computation paths.