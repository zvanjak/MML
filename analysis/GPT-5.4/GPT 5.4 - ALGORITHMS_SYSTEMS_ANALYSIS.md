# GPT-5.4 Algorithms and Systems Analysis

## Scope

This phase covers:

- `/mml/algorithms/**`
- `/mml/systems/**`

Primary focus areas:

- adaptive and stiff ODE/DAE solvers
- root finding and optimization
- eigensystem and matrix-analysis algorithms
- field/geometry/graph/fourier orchestration at algorithm level
- dynamical-system models and analyzers
- convergence control, event handling, chaos analysis, and model assumptions

## Executive Summary

This layer is ambitious and feature-rich. It contains much of the library’s practical numerical value.

The strongest parts are:

- broad algorithm coverage with generally good documentation
- adaptive ODE infrastructure with dense output and solver statistics
- multiple stiff/DAE solver implementations rather than a single nominal option
- meaningful analysis tooling for dynamical systems, eigenproblems, and field behavior
- a recurring attempt to standardize configs, result objects, and diagnostics

The main weaknesses are:

1. algorithm configuration contracts are not always honored by implementation
2. heuristic tolerances are still widely hard-coded and remain strongly `double`-centric
3. event detection and chaos-analysis code contain real correctness edge cases
4. some expensive or fragile algorithms hide significant assumptions behind simple APIs
5. breadth is ahead of consistency; there are modernized parts and legacy-style parts side by side

Overall assessment for algorithms + systems: very capable layer, but with several important correctness and reliability gaps.

## Major Strengths

### 1. ODE/DAE coverage is unusually broad

The algorithms layer supports:

- fixed-step ODE integration
- adaptive ODE integration with dense output
- event detection
- stiff ODE methods
- multiple DAE approaches including Backward Euler, BDF, Radau IIA, and Rosenbrock-style methods

That breadth is a real asset. MML is not just wrapping a few textbook methods; it is trying to be practically useful across different numerical regimes.

### 2. Result/config instrumentation is improving

Files such as:

- `ODESolverAdaptive.h`
- `OptimizationMultidim.h`
- `EigenSystemSolvers.h`
- `DAESolverBase.h`
- `DynamicalSystemAnalyzers.h`

show a consistent move toward richer results with:

- convergence flags
- iteration counts
- achieved tolerances
- status enums
- elapsed time
- function evaluation counts

This is the right direction for a library that users will need to debug numerically.

### 3. Systems layer is more than just canned models

The systems code provides:

- continuous systems with analytical Jacobians
- discrete maps
- fixed point analysis
- Lyapunov analysis
- bifurcation-oriented workflows
- combined reports

That gives the library a real applied-dynamics identity, not just low-level math primitives.

## Critical Correctness Findings

### 1. `LyapunovAnalyzer::GramSchmidtQR` can divide by zero or near-zero norms

Location:
- `mml/systems/DynamicalSystemAnalyzers.h`

In the Gram-Schmidt step:

- the code computes `norm = sqrt(sum(v_i^2))`
- it conditionally accumulates `log(norm)` only if `norm > Precision::DivisionSafetyThreshold`
- but then it always normalizes with `Q(i, j) = v[i] / norm`

Why this matters:

If a perturbation vector collapses to zero or near-zero, the code avoids logging it but still divides by `norm` unconditionally.

Impact:

- possible division by zero
- possible `NaN`/`Inf` contamination of the orthonormal frame
- corrupted Lyapunov exponents and downstream chaos classification

This is a real correctness defect in one of the most numerically sensitive analyzers in the library.

Required improvement:

- explicitly handle `norm <= threshold` before normalization
- either reinitialize the vector, abort with a clear failure status, or use a recovery strategy
- add regression tests for near-degenerate tangent-frame evolution

### 2. Adaptive ODE integrator configuration is not fully honored

Location:
- `mml/algorithms/ODESolvers/ODESolverAdaptive.h`

The config object advertises fields such as:

- `max_steps`
- `max_step_size`
- `output_interval`

But the main config-based `integrate(...)` implementation:

- does not actively enforce `max_steps` during the loop
- only checks after the solve whether accepted steps exceeded the value
- does not appear to apply `config.max_step_size` to clamp step growth in the wrapped legacy overload

Why this matters:

This is an API contract issue. Users can reasonably believe the config constrains solver behavior, but at least some fields are only partially or post-hoc honored.

Impact:

- long-running integrations may not stop when users expect
- maximum step-size constraints may be silently ignored
- solver behavior can diverge from documented configuration semantics

Required improvement:

- enforce `max_steps` in the main integration loop
- clamp proposed step sizes to `max_step_size`
- add tests proving every config field changes runtime behavior as documented

## Major Correctness Findings

### 3. Event detection only detects sign changes, so it can miss grazing/touching events

Location:
- `mml/algorithms/ODESolvers/ODESolverEventDetection.h`

The event detector triggers only when:

- `gPrev[i] * gNew[i] < 0`

That means it detects sign crossings, but not cases where:

- the event function lands exactly on zero at a step endpoint
- the trajectory touches zero tangentially and does not change sign
- the event starts exactly on the manifold

Why this matters:

These are standard event cases in physical simulation and hybrid systems.

Impact:

- missed impacts or switching events
- missed terminal conditions
- inconsistent behavior depending on step alignment rather than actual dynamics

Required improvement:

- explicitly handle `gPrev == 0`, `gNew == 0`, and small-magnitude endpoint cases
- define whether tangent events should be detectable and support that policy deliberately

### 4. Legacy adaptive ODE overload is vulnerable to invalid `outputInterval`

Location:
- `mml/algorithms/ODESolvers/ODESolverAdaptive.h`

The overload:

- `integrate(x0, t0, tEnd, outputInterval, ...)`

computes:

- `ceil((tEnd - t0) / outputInterval)`

without guarding `outputInterval <= 0`.

Why this matters:

The newer config API treats `output_interval = 0.0` as a special mode and substitutes a default. The lower-level overload does not defend itself similarly.

Impact:

- division-by-zero or invalid output sizing if misused directly
- inconsistent behavior between old and new entry points

Required improvement:

- validate the low-level overload explicitly
- either reject non-positive `outputInterval` or implement the same semantics as the config API

### 5. Radau IIA DAE solver freezes Jacobians across Newton iterations within a step

Location:
- `mml/algorithms/DAESolvers/DAERadauIIA.h`

The solver computes all Jacobians once at the step start:

- `system.allJacobians(t, x, y, ...)`

and then reuses those Jacobians through all Newton iterations for the coupled stage system.

Why this matters:

For strongly nonlinear DAEs, a frozen Jacobian can materially degrade Newton convergence or robustness compared with refreshed Jacobians.

Impact:

- slower or failed convergence on strongly nonlinear problems
- solver behavior that is more fragile than the high-level algorithm name suggests

This is not automatically wrong; frozen Jacobians can be a deliberate approximation. The problem is that the tradeoff is not surfaced clearly and may surprise users expecting full Newton robustness.

Required improvement:

- document that the step uses frozen Jacobians
- consider optional Jacobian refresh after failed or slow Newton iterations

### 6. Algorithm-level tolerance policy is still heavily hard-coded

Locations include:
- `EigenSystemSolvers.h`
- `EigenSolverHelpers.h`
- `OptimizationMultidim.h`
- `DAESolverBase.h`
- `FieldAnalyzers.h`
- `FunctionsAnalyzer.h`
- `MatrixAlg.h`
- `GraphSpectral.h`
- many computational-geometry files
- `systems/DynamicalSystemAnalyzer.h`
- `systems/DiscreteMaps.h`

Typical values include:

- `1e-10`
- `1e-8`
- `1e-6`
- ad hoc chaos thresholds like `0.01`

Why this matters:

This means the algorithm layer still carries a mostly `double`-native assumption even though the library advertises configurable `Real` precision.

Impact:

- float builds may be too strict and unstable
- long double builds may leave accuracy on the table or behave inconsistently
- subsystem behavior is harder to reason about globally

Required improvement:

- centralize algorithm defaults through `PrecisionValues<Real>` or algorithm-specific precision policies
- distinguish between physical heuristics and numeric tolerances

## Moderate Findings

### 7. Newton root finding is robustly guarded but computationally expensive by design

Location:
- `mml/algorithms/RootFinding/RootFindingMethods.h`

Strengths:

- bracket validation is explicit
- derivative finiteness is checked
- near-zero derivatives are rejected
- failure modes are converted into structured results

Weaknesses:

- derivative is recomputed numerically with `NDer4` at every iteration
- this makes the method more expensive than users may assume from the API name alone
- it also makes Newton behavior dependent on derivative-step heuristics rather than exact derivative information

This is acceptable if well documented, but users should understand this is “Newton with numerical derivative” rather than classic analytical Newton.

### 8. Dynamical-system classification uses heuristic thresholds that are not clearly separated from numerical tolerances

Locations:
- `mml/systems/DynamicalSystemAnalyzer.h`
- `mml/systems/DiscreteMaps.h`

Examples:

- `maxExponent > 1e-10` for Lyapunov time logic
- `sum < -1e-6` for dissipativity
- `abs(sum) < 1e-4` for conservative classification
- `maxExponent > 0.01` to label discrete-map chaos

Why this matters:

Some thresholds are physical/heuristic rather than numeric. That is fine, but the code does not strongly distinguish the two concepts.

Impact:

- classifications can flip based on arbitrary constants
- users may mistake heuristic labels for mathematically rigorous decisions

Required improvement:

- document heuristic meaning explicitly
- expose thresholds as user-configurable where appropriate

### 9. Eigensolver and matrix-analysis defaults are still mostly literal-valued

Locations:
- `EigenSystemSolvers.h`
- `EigenSolverHelpers.h`
- `MatrixAlg.h`
- `GraphSpectral.h`

These files are algorithmically useful, but the precision policy has not been unified with the rest of the library’s newer infrastructure.

This is a consistency problem more than a direct defect, but it matters because eigensystem behavior is sensitive to tolerances.

### 10. Computational geometry code uses many local epsilons

Locations include:
- `CompGeometry/Triangulation.h`
- `CompGeometry/VoronoiDiagram.h`
- `CompGeometry/ConvexHull3D.h`

This is common in geometry code, but the same pattern appears here: useful algorithms with local `1e-10` conventions rather than shared tolerance policy.

Given MML’s broader precision abstraction, this should eventually be normalized.

## Subsystem Assessment

### Adaptive ODE Solvers

Strengths:

- strong capability surface
- dense output support
- FSAL-aware stepper integration
- statistics collection
- improved config object

Weaknesses:

- config contract not fully enforced
- legacy and modern APIs coexist with different semantics
- initialization heuristics still include hard-coded constants

Assessment: strong subsystem with one important contract bug.

### Event Detection

Strengths:

- direction filtering
- event actions (`Stop`, `Restart`, `Continue`)
- dense-output-based event localization

Weaknesses:

- sign-change-only triggering misses important event classes
- comments say bisection while implementation uses Illinois/regula-falsi style localization
- endpoint handling is incomplete

Assessment: useful but not fully robust.

### DAE Solvers

Strengths:

- impressive breadth
- clear solver-family differentiation
- explicit Newton iteration control and diagnostics

Weaknesses:

- many hard-coded tolerances remain
- some implementations make strong convergence assumptions, such as frozen Jacobians
- step-floor heuristics use literal tiny constants

Assessment: ambitious and valuable, but deserves careful validation against hard problem sets.

### Root Finding and Optimization

Strengths:

- good method coverage
- structured results are improving
- many failure modes are surfaced cleanly

Weaknesses:

- derivative-free vs numerically differentiated semantics are not always obvious from function names alone
- tolerance defaults are still literal-heavy

Assessment: generally solid.

### Dynamical Systems

Strengths:

- unusually rich analysis surface for a library of this kind
- analytical Jacobians in many canonical models
- meaningful high-level workflows

Weaknesses:

- Lyapunov implementation has a real normalization bug
- classification heuristics are too implicit
- some thresholds should be configurable or better documented

Assessment: distinctive subsystem with real scientific value, but must be correctness-hardened.

## Recommended Priority Actions

### Priority 0

1. Fix the Lyapunov Gram-Schmidt normalization bug in `DynamicalSystemAnalyzers.h`.
2. Enforce `max_steps` and `max_step_size` in adaptive ODE integration loops.
3. Add explicit handling for endpoint/tangent events in ODE event detection.

### Priority 1

1. Audit DAE solvers for frozen-Jacobian assumptions and document or parameterize them.
2. Replace widespread literal tolerances in algorithms/systems with precision-aware defaults or explicit heuristic config fields.
3. Validate all legacy adaptive-ODE overloads for bad `outputInterval` inputs.

### Priority 2

1. Separate physical/heuristic classification thresholds from pure numerical tolerances in systems analysis.
2. Normalize computational geometry tolerance policy.
3. Add benchmark-and-validation suites for stiff ODEs, DAEs, chaotic systems, and event-heavy problems.

## Suggested Tests To Add

1. Lyapunov regression tests where perturbation-frame vectors collapse or nearly collapse.
2. Adaptive ODE config tests proving `max_steps`, `max_step_size`, and `output_interval` semantics are enforced.
3. Event detection tests for endpoint hits, grazing contacts, and same-sign tangent zeros.
4. DAE nonlinear test problems that compare frozen vs refreshed Jacobian behavior.
5. Precision-variant tests for eigensolvers, graph spectral routines, and computational geometry predicates.

## Provisional Grade For This Phase

Algorithms and systems subsystem grade: **B**

Reasoning:

- the capability set is excellent
- several algorithm families are thoughtfully implemented
- but correctness and contract consistency do not yet match the subsystem’s ambition
- a few defects and many tolerance-policy inconsistencies still materially affect trustworthiness

## Conclusion

This layer is where MML becomes genuinely powerful, but also where correctness discipline matters most. The library is already doing sophisticated work here: adaptive stepping, event detection, stiff solving, DAE integration, eigensystem analysis, and dynamical-systems diagnostics. That is impressive.

The problem is that implementation consistency has not fully kept pace with feature growth. The biggest immediate issue is not lack of algorithms; it is that some APIs promise more controlled behavior than the implementation strictly enforces, and some analyzers remain fragile at the exact edge cases that matter scientifically.

The next phase should expect the tools layer to be lower-risk numerically, while the final synthesis should emphasize that the library’s biggest challenge is no longer feature count. It is convergence of engineering discipline across all those features.
