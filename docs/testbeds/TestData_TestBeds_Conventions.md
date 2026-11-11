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

### Shared Metadata (Progressive Adoption)

`test_beds/test_case_metadata.h` defines the compact `MML::TestBeds::TestCaseMetadata`
vocabulary. The first adopted catalog is `LinearAlgEqTestBed`; other families retain
their existing local metadata until separately migrated.

| Field | Meaning |
|-------|---------|
| `source` | Repository-relative fixture location, not a claim about the original literature or verification |
| `category` | Catalog classification; not a certification of mathematical properties |
| `difficulty` | Optional existing rating, 1 (easy) through 4 (expert); absent if not supplied |
| `conditionNumber` | Optional stored reference condition number; infinity is a known singular case, not missing data |
| `tolerance` | Optional per-case suggested tolerance; an engaged zero is distinct from absence |
| `verificationMethod` | `Unspecified`, `Analytical`, `ExternalReference`, or `IndependentNumerical`; set only when the method is established |

Copying a constant from a fixture does **not** establish independent verification.
Current linear-system entries therefore leave the verification method unspecified,
and do not invent difficulty or tolerance values.

`LinearAlgEqTestBed::getMetadata(name)` returns a snapshot for any registered real,
symmetric, complex, overdetermined, or multi-RHS case. It reads the case's condition
field rather than maintaining a second mutable numerical catalog. Names are exact,
including the historical parentheses on names such as `symm_mat_3x3()`. Unknown names
throw `std::runtime_error`, as in the existing named accessors.

The runnable counterpart of this example is
`src/docs_demos/core/docs_demo_test_bed_metadata.cpp`
(`MML_DocsApp --test-bed-metadata`):

```cpp
#include <test_beds/linear_alg_eq_systems_test_bed.h>
#include <iostream>

using namespace MML::TestBeds;
const auto info = LinearAlgEqTestBed::getMetadata("symm_mat_3x3()");
if (info.conditionNumber)
    std::cout << info.source << ": condition number = " << *info.conditionNumber << '\n';
```

Determinants, singular values and ranks remain on the typed case objects. Existing
reference properties are now wired through symmetric/SPD registries, solver-specific
compatibility helpers, and the multi-RHS case (whose coefficient matrix equals
`mat_5x5`). No raw numerical fixtures were changed or recomputed.

The larger random cases `mat_8x8`, `mat_10x10`, and `mat_20x20` have no stored
condition/singular-value/rank/determinant reference bundle. Their condition metadata
is absent, and they remain excluded from conditioning filters. Legacy constructors
retain their zero sentinels and assumed full rank: those defaults are **not**
reference facts. Rectangular compatibility cases have spectral properties but no
determinant. The legacy `spd` group still contains the singular graph Laplacian;
its stored zero determinant, rank 5 and infinite condition number are preserved,
not replaced by full-rank defaults.

### Fixture Workflow

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
