# GPT-5.4 Engineering Roadmap

## Purpose

- Convert the completed analysis set into an execution-oriented refactor plan.
- Focus on high-leverage structural work before any major feature expansion.
- Define concrete work packages that can be scheduled, estimated, reviewed, and implemented incrementally.

## Inputs

- `GPT 5.4 - BASE_ANALYSIS.md`
- `GPT 5.4 - CORE_ANALYSIS.md`
- `GPT 5.4 - ALGORITHMS_SYSTEMS_ANALYSIS.md`
- `GPT 5.4 - TOOLS_ANALYSIS.md`
- `GPT 5.4 - FINAL_ANALYSIS.md`

## Roadmap Conventions

- Priority levels:
	- `P1` highest-value structural work
	- `P2` important follow-on work after P1 stabilization
	- `P3` worthwhile cleanup and polish
- Effort levels:
	- `S` up to 2 days
	- `M` 3 to 5 days
	- `L` 1 to 2 weeks
	- `XL` multi-week cross-module work
- Risk levels:
	- `Low` localized change, low breakage risk
	- `Medium` public API or behavior may shift slightly
	- `High` broad transitive or public contract impact

## Executive Direction

MML should spend its next engineering cycle on consolidation, not expansion. The core objective is to make the library behave like one system across five axes:

- API contracts
- architectural layering
- numeric policy
- solver and analyzer maturity
- tools and persistence discipline

The recommended implementation order is:

1. establish common contracts
2. repair the worst dependency boundaries
3. standardize numeric policy
4. modernize lagging solver and analyzer subsystems
5. stabilize tools, serialization, and large-workflow support

## Delivery Phases

### Phase 0. Baseline and guardrails

Goal: create safety rails before refactoring shared headers.

| ID | Priority | Task | Target files | Effort | Risk | Done condition |
|----|----------|------|--------------|--------|------|----------------|
| R0.1 | P1 | Add self-contained header compile checks for representative root, base, core, algorithms, and tools headers | `mml/MMLBase.h`, `mml/base/Vector/Vector.h`, `mml/base/Matrix/Matrix.h`, `mml/core/Integration.h`, `mml/algorithms/RootFinding.h`, `mml/tools/ThreadPool.h`, test or CMake support files | `M` | `Medium` | CI or local build can compile selected headers without relying on transitive includes |
| R0.2 | P1 | Add architectural inventory note listing intended layer boundaries and allowed dependencies | `analysis/` docs or new architecture doc under `docs/` | `S` | `Low` | One checked-in document defines allowed dependency direction between root, base, core, algorithms, systems, and tools |
| R0.3 | P1 | Add regression checklist for public API breakage risk before changing umbrella headers | `docs/` and review checklist docs | `S` | `Low` | Refactor PRs can reference one stable checklist for API, include, and behavior review |

### Phase 1. Library-wide contract standardization

Goal: define a common way MML reports success, failure, status, diagnostics, and configuration.

| ID | Priority | Task | Target files | Effort | Risk | Done condition |
|----|----------|------|--------------|--------|------|----------------|
| R1.1 | P1 | Define a shared algorithm operation contract for config, status, diagnostics, and result metadata | `mml/core/AlgorithmTypes.h`, selected headers under `mml/algorithms/**`, selected headers under `mml/systems/**` | `L` | `High` | One reusable contract exists and is documented as the default pattern for new algorithm APIs |
| R1.2 | P1 | Decide and document when APIs should throw, return result objects, or do both | `mml/MMLExceptions.h`, `docs/ERROR_HANDLING.md`, `docs/CODING_STYLE.md` | `M` | `Medium` | Error model guidance is explicit, internally consistent, and referenced from coding guidance |
| R1.3 | P1 | Normalize result fields across strong algorithm families | `mml/algorithms/RootFinding/RootFindingBase.h`, `mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h`, `mml/algorithms/Optimization.h`, `mml/systems/DynamicalSystemTypes.h` | `L` | `Medium` | Result objects expose aligned status and diagnostic vocabulary where applicable |
| R1.4 | P2 | Introduce shared helpers for common preset factories such as `Fast`, `HighPrecision`, and related execution modes | `mml/core/AlgorithmTypes.h`, selected config structs in algorithms and systems | `M` | `Medium` | Factory presets follow one naming and semantics scheme across subsystems |
| R1.5 | P2 | Extend the same contract discipline into tools operations that can fail | `mml/tools/serializer/SerializerBase.h`, `mml/tools/DataLoader.h`, `mml/tools/Visualizer.h` | `L` | `Medium` | Tools-layer operations no longer mix unrelated result idioms without a documented reason |

### Phase 2. Architectural boundary cleanup

Goal: make dependency direction match the intended layering.

| ID | Priority | Task | Target files | Effort | Risk | Done condition |
|----|----------|------|--------------|--------|------|----------------|
| R2.1 | P1 | Remove visualization dependency from the foundational root layer | `mml/MMLBase.h`, `mml/MMLVisualizators.h`, `mml/single_header/MML.h` | `M` | `High` | `MMLBase.h` no longer includes visualization headers and still supports its intended low-level role |
| R2.2 | P1 | Redefine visualization as an opt-in high-level facility | `mml/tools/Visualizer.h`, `mml/MMLVisualizators.h`, `mml/single_header/MML.h`, any related docs | `L` | `High` | Visualization is reachable from opt-in entry points, not from foundational headers |
| R2.3 | P1 | Audit and reduce heavy convenience include fan-out in base/root umbrella paths | `mml/base/BaseUtils.h`, `mml/base/Vector/VectorTypes.h`, root umbrellas | `L` | `Medium` | Common low-level includes pull in meaningfully fewer unrelated modules |
| R2.4 | P2 | Split dense high-level headers while preserving umbrella entry points | `mml/core/FieldOperations.h`, `mml/core/Curves.h`, `mml/core/Surfaces.h` | `L` | `Medium` | Large headers are internally decomposed with umbrella includes preserved for users |
| R2.5 | P2 | Clarify long-term single-header policy versus layered include policy | `mml/single_header/MML.h`, `docs/QUICK_START.md`, `docs/API_CHEATSHEET.md` | `M` | `Medium` | Docs explain when to use single-header aggregation versus focused include paths |

### Phase 3. Numeric policy normalization

Goal: make `Real`, tolerances, validation, and singularity behavior consistent across the codebase.

| ID | Priority | Task | Target files | Effort | Risk | Done condition |
|----|----------|------|--------------|--------|------|----------------|
| R3.1 | P1 | Convert library policy constants from `double` to `Real` where they represent configured numeric behavior | `mml/MMLBase.h`, related precision policy files | `M` | `Medium` | Default tolerances and thresholds honor the configured scalar type without silent narrowing or widening |
| R3.2 | P1 | Audit numerical validation inside iterative and refinement loops | `mml/core/Integration/Integration1D.h`, solver loops in `mml/core/**`, `mml/algorithms/**` | `L` | `Medium` | Intermediate non-finite states are checked consistently before refinement continues |
| R3.3 | P1 | Make singularity handling a real cross-core contract | `mml/core/SingularityHandling.h`, `mml/core/Fields.h`, `mml/core/FieldOperations.h`, `mml/core/CoordTransf.h`, `mml/core/Curves.h`, `mml/core/Surfaces.h` | `XL` | `High` | Similar singular cases are handled through shared policy rather than scattered local logic |
| R3.4 | P2 | Introduce explicit tests or demos for non-`double` `Real` configurations | build config, tests, selected docs | `L` | `Medium` | At least one non-default `Real` build path is validated and documented |
| R3.5 | P2 | Standardize evaluation-count, tolerance, and convergence semantics across iterative subsystems | `mml/core/AlgorithmTypes.h`, iterative solvers, root finding, ODE and DAE code | `L` | `Medium` | Public docs and result types use consistent meanings for convergence-related fields |

### Phase 4. Solver and analyzer modernization

Goal: bring lagging numerical subsystems up to the standard of root finding and adaptive ODE integration.

| ID | Priority | Task | Target files | Effort | Risk | Done condition |
|----|----------|------|--------------|--------|------|----------------|
| R4.1 | P1 | Redesign DAE integration around adaptive, tolerance-driven execution where feasible | `mml/algorithms/DAESolvers/DAESolverBase.h`, related DAE solver headers | `XL` | `High` | DAE configuration and diagnostics are meaningfully closer to the ODE integrator model |
| R4.2 | P1 | Align DAE result reporting with ODE solver statistics and status semantics | DAE solver headers and shared algorithm types | `L` | `Medium` | DAE outputs expose convergence and execution diagnostics in the same language as ODE outputs |
| R4.3 | P1 | Modernize dynamical-system analyzers to config/result-based entry points | `mml/systems/DynamicalSystemAnalyzers.h`, `mml/systems/DynamicalSystemTypes.h` | `L` | `Medium` | Static bare-parameter APIs are replaced or wrapped by extensible config/result entry points |
| R4.4 | P2 | Add stronger uncertainty and coverage metadata to sampling-based analyzer outputs | `mml/algorithms/FieldAnalyzers.h`, systems analyzer result types | `M` | `Low` | Sampling density, search coverage, or approximation caveats are encoded in results or docs |
| R4.5 | P2 | Add analytic fast paths for common coordinate transforms while preserving numerical fallback | `mml/core/CoordTransf.h` and related transform implementations | `L` | `Medium` | Common transforms can avoid repeated numerical Jacobian estimation |
| R4.6 | P3 | Unify low-degree and general polynomial root entry points behind a clearer front door | `mml/algorithms/RootFinding/RootFindingPolynoms.h`, `mml/algorithms/RootFinding.h` | `M` | `Low` | Polynomial root solving has a more coherent public entry strategy |

### Phase 5. Tools, serialization, and workflow stabilization

Goal: make the support layer composable, version-safe, and more scalable.

| ID | Priority | Task | Target files | Effort | Risk | Done condition |
|----|----------|------|--------------|--------|------|----------------|
| R5.1 | P1 | Standardize tools-layer failure and result contracts | `mml/tools/serializer/SerializerBase.h`, `mml/tools/DataLoader.h`, `mml/tools/Visualizer.h`, related helpers | `L` | `Medium` | Tools operations follow one intentional contract family rather than a mix of unrelated patterns |
| R5.2 | P1 | Strengthen serialization versioning, compatibility checks, and metadata policy | `mml/tools/Serializer.h`, `mml/tools/serializer/**`, `mml/tools/README_Serialization.md` | `L` | `Medium` | Serialized formats advertise version and compatibility expectations clearly and consistently |
| R5.3 | P2 | Add scalable loading paths for larger datasets | `mml/tools/DataLoader.h`, `mml/tools/data_loader/**` | `L` | `Medium` | Loading APIs support chunked or streaming-oriented workflows in addition to in-memory convenience paths |
| R5.4 | P2 | Clarify thread-pool lifecycle guarantees and optional cancellation strategy | `mml/tools/ThreadPool.h`, docs | `M` | `Low` | `ThreadPool` behavior under long-running or failing tasks is explicitly documented, and optional extension points are identified |
| R5.5 | P3 | Tighten timer/reporting configurability and output semantics | `mml/tools/Timer.h` | `S` | `Low` | Timing utilities offer clearer formatting and usage guarantees without losing simplicity |

## Workstream View

### Workstream A. Common contracts

- Start with `mml/core/AlgorithmTypes.h` as the seed abstraction.
- Use `RootFindingBase.h` and `ODEAdaptiveIntegrator.h` as the reference-quality examples.
- Do not attempt to force every helper function into a heavyweight result type.
- Focus on operations where failure, convergence, timing, or diagnostics matter.

### Workstream B. Dependency hygiene

- First remove foundational visualization coupling.
- Then reduce umbrella fan-out and transitive include dependence.
- Only after boundaries are cleaner should header splitting happen.

### Workstream C. Numeric consistency

- Treat this as both a correctness and maintenance project.
- Prefer shared policy surfaces over localized conditionals.
- Explicitly test one non-default `Real` mode to prove the policy is real, not aspirational.

### Workstream D. Solver maturity leveling

- DAE and analyzer work should adopt patterns from the best existing subsystems.
- Avoid inventing new config/result idioms during this pass.
- Prioritize symmetry with ODE and root-finding APIs where the math allows it.

### Workstream E. Tools hardening

- Keep tools useful and pragmatic.
- Avoid pushing operational helpers back into foundational layers.
- Serialization policy should be treated as a long-term compatibility commitment.

## Recommended PR Sequence

1. PR 1: header self-sufficiency checks and architecture boundary doc
2. PR 2: remove visualization from `MMLBase.h` and adjust single-header aggregation
3. PR 3: convert `Defaults` numeric policy constants to `Real` and add precision-focused tests
4. PR 4: define shared algorithm contract updates in `AlgorithmTypes.h`
5. PR 5: migrate root finding, optimization, and ODE diagnostics onto the tightened shared vocabulary
6. PR 6: modernize systems analyzer entry points and result models
7. PR 7: start DAE contract alignment, then adaptive DAE evolution in follow-up PRs
8. PR 8: unify tools result contracts and strengthen serialization versioning
9. PR 9: add scalable data-loading paths and thread-pool lifecycle clarifications
10. PR 10: split oversized core headers and complete remaining include hygiene cleanup

## Acceptance Criteria By Milestone

### Milestone A. Structural baseline

- Foundational headers compile without hidden transitive include dependencies.
- `MMLBase.h` no longer pulls visualization concerns.
- Architecture guidance exists and matches real dependency direction.

### Milestone B. Contract baseline

- Public numerical subsystems expose a recognizable, documented contract family.
- Error-model guidance is explicit and consistent with implementation.
- Result/status naming is aligned across the strongest algorithm families.

### Milestone C. Numeric baseline

- Library defaults actually honor `Real`.
- Numerical loops validate intermediate states consistently.
- Singularity behavior is no longer scattered ad hoc across adjacent core modules.

### Milestone D. Maturity baseline

- DAE is no longer visibly one generation behind ODE in API quality.
- Systems analyzers support extensible config/result entry points.
- Sampling-based analyzers communicate approximation limits explicitly.

### Milestone E. Tools baseline

- Tools modules share a coherent failure/result model.
- Serialization format behavior is versioned and documented.
- Data loading can support larger workflows without requiring full in-memory parsing.

## Suggested Ownership Split

| Area | Suggested owner profile |
|------|-------------------------|
| Common contracts and error model | library architect or lead maintainer |
| Header and dependency hygiene | maintainer comfortable with public API and compile-time impact |
| Numeric policy and singularity behavior | numerical maintainer with cross-core familiarity |
| DAE and systems modernization | numerical algorithms maintainer |
| Tools and serialization hardening | tooling and platform-focused maintainer |

## Immediate Next Steps

1. Approve the roadmap structure and sequencing.
2. Start with `R0.1`, `R1.2`, and `R2.1` because they create the clearest foundation for the rest.
3. Treat `R4.1` as a dedicated project after common contracts and numeric policy are stabilized.

## Bottom Line

This roadmap is designed to turn the analysis set into incremental engineering work instead of a one-time review artifact. The highest-value path is to stabilize contracts and architecture first, then normalize numeric policy, then modernize lagging subsystems, and only after that expand tooling and ergonomics further.