# Test Data Mapping

## Purpose

This document maps algorithm/test areas to the shared `test_beds/` and `test_data/` catalog. It is meant to answer three questions quickly:

- Which reusable test bed should a new or existing test use?
- Which raw known-answer fixtures back that test bed?
- Where are current coverage gaps or direct-fixture exceptions?

See `docs/testbeds/TestData_TestBeds_Conventions.md` for the naming and include conventions.

## Include Rule

Tests should normally consume `test_beds/*_test_bed.h`. Direct `test_data/*_defs.h` includes are allowed only when a test intentionally validates raw fixture values that are not yet exposed through the test-bed API.

Current direct fixture exception:

| Test File | Direct Fixture Headers | Reason / Follow-up |
| --- | --- | --- |
| `tests/algorithms/eigensystem_solvers_tests.cpp` | `linear_alg_eq_systems_classic_defs.h`, `linear_alg_eq_systems_sym_defs.h`, `linear_alg_eq_systems_solverspecific_defs.h` | Eigensolver tests combine raw matrix families with `linear_alg_eq_systems_test_bed.h` and `eigenvalue_test_bed.h`. Follow-up: expose any still-needed raw cases through the relevant test-bed APIs before removing direct includes. |

## Algorithm-To-Testbed Map

| Algorithm / Test Area | Primary Test-Bed Header | Raw Test-Data Header(s) | Known-Answer Type | Current Usage | Coverage Gaps / Notes |
| --- | --- | --- | --- | --- | --- |
| Linear algebra systems and solvers | `test_beds/linear_alg_eq_systems_test_bed.h` | `linear_alg_eq_systems_real_defs.h`, `linear_alg_eq_systems_sym_defs.h`, `linear_alg_eq_systems_complex_defs.h`, `linear_alg_eq_systems_classic_defs.h`, `linear_alg_eq_systems_solverspecific_defs.h`, `linear_alg_eq_systems_tridiag_defs.h`, `linear_alg_eq_systems_singular_defs.h` | Matrices, RHS vectors, known solutions, eigenvalues, eigenvectors, determinants, singular values, condition numbers, ranks | Used by matrix algorithms, matrix utilities, eigensystem tests, and linear solver tests | Band/sparse fixture expansion is tracked by `MinimalMathLibrary-lq0l`. Some tridiagonal/singular fixture families are not yet exposed by the main test-bed registry. |
| Eigensolvers and eigenvalue algorithms | `test_beds/eigenvalue_test_bed.h` | `eigenvalue_defs.h`; also consumes linear algebra fixture families | Eigenvalues, eigenvectors, characteristic cases, symmetric/general matrix cases | Used by `tests/algorithms/eigensystem_solvers_tests.cpp` alongside linear algebra test beds | Direct raw linear-system includes remain in eigensystem tests; consider wrapping those cases in the eigenvalue or linear-algebra test-bed interface. |
| Root finding | `test_beds/root_finding_test_bed.h` | `root_finding_defs.h` | Functions, derivatives, known roots, multiplicities, brackets, singularity flags, difficulty | Used by advanced root-finding tests | Basic `root_finding_tests.cpp`, polynomial root tests, and nonlinear-system solver tests currently do not include the shared test bed directly. Consider moving reusable cases into the test bed as they mature. |
| Integration / quadrature | `test_beds/integration_test_bed.h`, `test_beds/real_functions_test_bed.h` | `integration_defs.h` | Integrands, exact results, antiderivatives, singularity/improper-integral metadata, tolerances | Core integration tests use `real_functions_test_bed.h`; dedicated integration test bed exists for richer integral categories | Broader integration files (`adaptive_integration_nd`, `gauss_kronrod`, 2D/3D/path/surface/Monte Carlo) should be reviewed for reusable known-answer cases that belong in shared beds. |
| Interpolation / interpolated functions | `test_beds/interpolation_test_bed.h` | `interpolation_defs.h` | Sample nodes/values, true functions, expected interpolation behavior, deterministic generated cases | Used by `tests/base/interpolated_functions/interpolated_functions_tests.cpp` | `MinimalMathLibrary-m0m1` tracks adding more algorithm tests using this bed. Current bed is a strong model for mixing static fixtures with deterministic generators. |
| Optimization | `test_beds/optimization_test_bed.h`, `test_beds/scalar_functions_test_bed.h` | `optimization_defs.h` | Objective functions, known minima, bounds, gradients, multidimensional cases, difficulty | Used by one-dimensional and multidimensional optimization tests | Several specialized optimization tests do not yet consume the shared test bed directly. Candidate follow-up: map bound-constrained, projected-gradient, simulated-annealing, and linear-programming cases into shared fixtures where reusable. |
| Real/scalar/vector functions and numerical differentiation | `test_beds/real_functions_test_bed.h`, `test_beds/scalar_functions_test_bed.h`, `test_beds/vector_functions_test_bed.h` | None separate for these beds today | Functions, domains/test intervals, derivatives, integrals, gradients, Jacobians, divergence/curl metadata | Used by core derivation, precision derivation, precision integration, and core integration tests | Function beds are rich but mostly self-contained. `MinimalMathLibrary-nupm` tracks special-function expansion for real functions. |
| ODE systems | `test_beds/diff_eq_systems_test_bed.h` | `diff_eq_systems_defs.h` | ODE systems, exact/reference solutions, initial conditions, intervals, tolerances | ODE solver tests currently use local/system-specific cases; the shared bed exists and should become the default reusable source | Review `ode_system_solvers_tests.cpp`, step-calculator tests, event-detection tests, and BVP tests for cases that should be promoted into the shared bed. |
| Stiff ODE systems | `test_beds/stiff_ode_test_bed.h` | `stiff_ode_defs.h` | Stiff equations/systems, reference solutions, stiffness indicators, tolerances | Dedicated stiff test bed exists | Ensure stiff solvers and DAE tests consume shared stiff cases where applicable. |
| Fourier / signals | `test_beds/fourier_test_bed.h` | None separate today | Signals, expected transforms/spectra, convolution and correlation reference cases | Used by Fourier tests | Consider whether raw signal datasets should remain embedded or move to `test_data/fourier_defs.h` if the catalog grows. |
| Parametric curves | `test_beds/parametric_curves_test_bed.h` | `parametric_curves_defs.h` | Curves, derivatives, arc-length data, curvature/torsion-style reference values | Used by curve analyzer and core curves tests | Strong mapping exists. Keep adding reusable geometric edge cases here instead of per-test literals. |
| Parametric surfaces | `test_beds/parametric_surfaces_test_bed.h` | `parametric_surfaces_defs.h` | Surfaces, explicit surfaces, derivatives, normals/geometry reference values | Shared bed exists; not yet broadly included by current surface-related tests | Review `surfaces_tests.cpp`, derivation surface tests, and surface integration tests for promotable known-answer cases. |
| Statistics | `test_beds/statistics_test_bed.h` | `statistics_defs.h` | Reference datasets, expected descriptive statistics, hypothesis-test values, correlations, distributions | Statistics algorithm tests currently appear mostly standalone; shared data exists | Reconcile which statistics datasets should live in `statistics_defs.h` and which should be exposed through `statistics_test_bed.h`; then update statistics tests to consume shared fixtures. |
| Matrix tridiagonal / banded / sparse basics | `test_beds/linear_alg_eq_systems_test_bed.h` where applicable | `linear_alg_eq_systems_tridiag_defs.h`, sparse/band fixtures not yet centralized | Known matrices and solver/reference outputs | Matrix tridiagonal and sparse matrix tests currently appear mostly direct/local | `MinimalMathLibrary-lq0l` should expand and expose band diagonal and sparse matrix fixture families. |
| Computational geometry | No shared test bed yet | None | Geometric configurations, expected intersections/hulls/triangulations/voronoi outputs | Current computational geometry tests are standalone | Candidate future package: `comp_geometry_test_bed.h` if repeated known geometries emerge across algorithms. |
| Algebra / groups / finite fields | No shared test bed yet | None | Algebraic structures, Cayley tables, representation examples, known laws/properties | Current algebra tests are standalone | A shared bed may be useful only if examples repeat across multiple algebra modules. |
| Coordinate systems, tensors, differential geometry | Parametric curve/surface and vector/scalar function beds cover part of this area | `parametric_curves_defs.h`, `parametric_surfaces_defs.h` | Known coordinate transforms, fields, forms, tensor identities, geometry reference values | Current tests are mostly standalone with some parametric-curve usage | Consider dedicated beds only when cases are reused across transform, field, tensor, and differential-geometry tests. |
| Local MPL physics tests, including electromagnetism | No shared test bed yet | None | Analytic field configurations, potentials, flux/circulation reference values | Current tests under `src/book/mpl_tests` are standalone | Candidate future bed if multiple physics modules share the same analytic configurations. |

## Coverage Priorities

Near-term priorities from the current map:

1. Use `interpolation_test_bed.h` in more interpolation algorithm tests (`MinimalMathLibrary-m0m1`).
2. Expand linear algebra fixtures for banded and sparse systems (`MinimalMathLibrary-lq0l`).
3. Expand complex linear systems, especially Hermitian/unitary/signal-processing style cases (`MinimalMathLibrary-ytkx`).
4. Add special functions and harder pathological functions to `real_functions_test_bed.h` (`MinimalMathLibrary-nupm`).
5. Promote reusable statistics, ODE, surface, and optimization cases from standalone tests into shared test beds as they become repeated patterns.

## Maintenance Checklist

When adding a new algorithm test:

1. Check this map for an existing test bed.
2. Prefer including the test-bed header over the raw fixture header.
3. If adding reusable known-answer data, add raw values to `test_data/*_defs.h` and expose a named case through `test_beds/*_test_bed.h`.
4. Update this mapping when a new test bed, fixture family, or intentional direct `test_data` exception is introduced.
