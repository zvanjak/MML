# Test Data and Test Bed Conventions

## Purpose

MML keeps reusable mathematical problems, datasets, and systems close at hand for algorithm validation. The `test_data/` and `test_beds/` folders are the shared catalog for those known-answer cases.

Use this convention when adding, moving, or consuming test fixtures.

## Two-Layer Model

### `test_data/*_defs.h`

`test_data/` is the raw fixture and reference-value layer. Files in this folder should define verified mathematical facts:

- matrices, vectors, right-hand sides, and known solutions
- functions and reference constants
- exact roots, integrals, derivatives, eigenvalues, singular values, ranks, condition numbers, and expected statistics
- deterministic static samples or source datasets

The file naming convention is:

```text
test_data/<area>_defs.h
```

Examples:

```text
test_data/linear_alg_eq_systems_real_defs.h
test_data/root_finding_defs.h
test_data/statistics_defs.h
```

### `test_beds/*_test_bed.h`

`test_beds/` is the algorithm-facing interface layer. Files in this folder should wrap raw fixtures into named, typed, and queryable test cases:

- typed test-case structures
- case names and descriptions
- categories, difficulty, tolerances, and other metadata
- grouped accessors such as `getAll...()`, `get...ByCategory()`, or test-bed classes
- deterministic generators that produce reproducible test cases

The file naming convention is:

```text
test_beds/<area>_test_bed.h
```

Examples:

```text
test_beds/linear_alg_eq_systems_test_bed.h
test_beds/root_finding_test_bed.h
test_beds/statistics_test_bed.h
```

## Namespace Convention

All test data and test beds should live under `MML::TestBeds` or a nested namespace below it.

Examples:

```cpp
namespace MML::TestBeds
{
    // Common test fixtures and test beds
}

namespace MML::TestBeds::Statistics
{
    // Statistics-specific datasets and expected values
}
```

Do not introduce parallel spellings such as `MML::Testbeds`.

## Include Convention

Tests should normally include the `test_beds/` header for an algorithm area. That keeps tests insulated from raw fixture organization and gives them access to categories, metadata, and grouped accessors.

Direct includes from `test_data/` are allowed when a test intentionally validates or exercises raw reference values that are not exposed through the test-bed API. These should be treated as explicit exceptions, not the default pattern.

Current direct fixture include exception:

```text
tests/algorithms/eigensystem_solvers_tests.cpp
```

This test currently includes classic, symmetric, and solver-specific linear-system fixture definitions directly alongside the linear-algebra and eigenvalue test beds.

## Adding New Cases

When adding a new known-answer problem:

1. Put raw verified values in the matching `test_data/*_defs.h` file.
2. Expose the case through the matching `test_beds/*_test_bed.h` file.
3. Include source/provenance where practical, especially for literature benchmarks.
4. Record expected tolerance or difficulty in the test-bed layer when the case is numerically challenging.
5. Prefer deterministic generation over random generation. If randomness is required, use an explicit seed and document it.

## Cleanup Rules

- Keep one canonical source for each fixture family.
- Do not leave `_new`, `_work`, backup, or scratch fixture headers in `test_data/`.
- Prefer renaming outliers to the standard convention rather than adding compatibility shims.
- After any rename or fixture cleanup, build `MML_Tests` and run the built test executable directly.
