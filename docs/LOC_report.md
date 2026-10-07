# MML Lines of Code Report

Generated: 2026-08-23

This report counts source-like C and C++ files under `mml/` and all nested subfolders.

## Counting Method

- Root folder: `mml/`
- Included extensions: `.h`, `.hpp`, `.hh`, `.hxx`, `.c`, `.cc`, `.cpp`, `.cxx`, `.inl`, `.ipp`
- Physical lines include all lines in included files, including blank lines and comments.
- Non-blank lines exclude only empty or whitespace-only lines.
- Comment-only lines are still counted as non-blank lines.
- `mml/single_header/MML.h` is included in the main totals because it is part of the tree, but it is also called out separately because it is an amalgamated header and duplicates much of the library surface.

## Summary

| Scope | Source files | Physical lines | Non-blank lines | Directories |
|---|---:|---:|---:|---:|
| `mml/` including `single_header` | 351 | 212,473 | 183,238 | 48 |
| `mml/` excluding `single_header` | 350 | 110,128 | 95,538 | 47 |
| `mml/single_header` only | 1 | 102,345 | 87,700 | 1 |

## By Extension

| Extension | Files | Physical lines | Non-blank lines |
|---|---:|---:|---:|
| `.h` | 351 | 212,473 | 183,238 |

## Top-Level Folders

| Path | Direct files | Direct lines | Direct non-blank | Recursive files | Recursive lines | Recursive non-blank |
|---|---:|---:|---:|---:|---:|---:|
| `mml` | 6 | 1,590 | 1,358 | 351 | 212,473 | 183,238 |
| `mml/algorithms` | 12 | 7,439 | 6,457 | 83 | 35,415 | 30,550 |
| `mml/base` | 22 | 6,561 | 5,645 | 100 | 28,586 | 24,348 |
| `mml/core` | 18 | 5,331 | 4,666 | 104 | 29,150 | 25,625 |
| `mml/interfaces` | 15 | 2,601 | 2,335 | 15 | 2,601 | 2,335 |
| `mml/single_header` | 1 | 102,345 | 87,700 | 1 | 102,345 | 87,700 |
| `mml/systems` | 5 | 2,462 | 2,124 | 12 | 3,381 | 2,948 |
| `mml/tools` | 9 | 3,158 | 2,814 | 30 | 9,405 | 8,374 |

## Recursive Directory Inventory

| Path | Direct files | Direct lines | Direct non-blank | Recursive files | Recursive lines | Recursive non-blank |
|---|---:|---:|---:|---:|---:|---:|
| `mml` | 6 | 1,590 | 1,358 | 351 | 212,473 | 183,238 |
| `mml/algorithms` | 12 | 7,439 | 6,457 | 83 | 35,415 | 30,550 |
| `mml/algorithms/Analyzers` | 3 | 1,810 | 1,534 | 3 | 1,810 | 1,534 |
| `mml/algorithms/CompGeometry` | 10 | 3,515 | 2,972 | 10 | 3,515 | 2,972 |
| `mml/algorithms/DAESolvers` | 9 | 2,236 | 1,944 | 9 | 2,236 | 1,944 |
| `mml/algorithms/Fourier` | 5 | 2,014 | 1,683 | 5 | 2,014 | 1,683 |
| `mml/algorithms/ODESolvers` | 9 | 4,166 | 3,552 | 9 | 4,166 | 3,552 |
| `mml/algorithms/Optimization` | 4 | 1,474 | 1,323 | 20 | 6,783 | 5,865 |
| `mml/algorithms/Optimization/Constraints` | 4 | 738 | 637 | 4 | 738 | 637 |
| `mml/algorithms/Optimization/LP` | 6 | 2,398 | 2,021 | 6 | 2,398 | 2,021 |
| `mml/algorithms/Optimization/Multidim` | 6 | 2,173 | 1,884 | 6 | 2,173 | 1,884 |
| `mml/algorithms/RootFinding` | 8 | 3,515 | 3,149 | 8 | 3,515 | 3,149 |
| `mml/algorithms/Statistics` | 7 | 3,937 | 3,394 | 7 | 3,937 | 3,394 |
| `mml/base` | 22 | 6,561 | 5,645 | 100 | 28,586 | 24,348 |
| `mml/base/Algebra` | 12 | 1,421 | 1,231 | 16 | 1,861 | 1,612 |
| `mml/base/Algebra/LieGroups` | 4 | 440 | 381 | 4 | 440 | 381 |
| `mml/base/BaseUtils` | 6 | 1,387 | 1,205 | 6 | 1,387 | 1,205 |
| `mml/base/DifferentialGeometry` | 6 | 857 | 733 | 6 | 857 | 733 |
| `mml/base/Geometry` | 5 | 290 | 251 | 20 | 6,239 | 5,192 |
| `mml/base/Geometry/Geometry2DCore` | 4 | 1,558 | 1,322 | 4 | 1,558 | 1,322 |
| `mml/base/Geometry/Geometry3DBodiesCore` | 5 | 1,950 | 1,572 | 5 | 1,950 | 1,572 |
| `mml/base/Geometry/Geometry3DCore` | 3 | 1,069 | 913 | 3 | 1,069 | 913 |
| `mml/base/Geometry/GeometryCore` | 3 | 1,372 | 1,134 | 3 | 1,372 | 1,134 |
| `mml/base/InterpolatedFunctions` | 8 | 2,182 | 1,806 | 8 | 2,182 | 1,806 |
| `mml/base/Matrix` | 7 | 4,380 | 3,699 | 7 | 4,380 | 3,699 |
| `mml/base/SparseMatrix` | 4 | 1,442 | 1,232 | 4 | 1,442 | 1,232 |
| `mml/base/Tensor` | 6 | 1,639 | 1,410 | 6 | 1,639 | 1,410 |
| `mml/base/Vector` | 5 | 2,038 | 1,814 | 5 | 2,038 | 1,814 |
| `mml/core` | 18 | 5,331 | 4,666 | 104 | 29,150 | 25,625 |
| `mml/core/Algebra` | 7 | 951 | 852 | 7 | 951 | 852 |
| `mml/core/CoordTransf` | 6 | 2,145 | 1,879 | 6 | 2,145 | 1,879 |
| `mml/core/Derivation` | 14 | 4,771 | 4,267 | 14 | 4,771 | 4,267 |
| `mml/core/DifferentialGeometry` | 8 | 1,365 | 1,164 | 8 | 1,365 | 1,164 |
| `mml/core/Fields` | 3 | 1,950 | 1,746 | 3 | 1,950 | 1,746 |
| `mml/core/FunctionSpaces` | 16 | 1,246 | 1,072 | 16 | 1,246 | 1,072 |
| `mml/core/Integration` | 12 | 5,246 | 4,680 | 12 | 5,246 | 4,680 |
| `mml/core/LinAlgEqSolvers` | 4 | 2,615 | 2,322 | 4 | 2,615 | 2,322 |
| `mml/core/OrthogonalBasis` | 4 | 761 | 650 | 4 | 761 | 650 |
| `mml/core/SparseSolvers` | 6 | 1,664 | 1,409 | 6 | 1,664 | 1,409 |
| `mml/core/VectorSpaces` | 6 | 1,105 | 918 | 6 | 1,105 | 918 |
| `mml/interfaces` | 15 | 2,601 | 2,335 | 15 | 2,601 | 2,335 |
| `mml/single_header` | 1 | 102,345 | 87,700 | 1 | 102,345 | 87,700 |
| `mml/systems` | 5 | 2,462 | 2,124 | 12 | 3,381 | 2,948 |
| `mml/systems/DynamicalSystem` | 7 | 919 | 824 | 7 | 919 | 824 |
| `mml/tools` | 9 | 3,158 | 2,814 | 30 | 9,405 | 8,374 |
| `mml/tools/data_loader` | 4 | 1,420 | 1,253 | 4 | 1,420 | 1,253 |
| `mml/tools/persistence` | 9 | 2,819 | 2,512 | 9 | 2,819 | 2,512 |
| `mml/tools/serializer` | 8 | 2,008 | 1,795 | 8 | 2,008 | 1,795 |

## Largest Source Files

| Rank | File | Physical lines | Non-blank lines |
|---:|---|---:|---:|
| 1 | `mml/single_header/MML.h` | 102,345 | 87,700 |
| 2 | `mml/algorithms/GraphAlgorithms.h` | 2,111 | 1,856 |
| 3 | `mml/algorithms/Statistics/Distributions.h` | 1,502 | 1,301 |
| 4 | `mml/base/Matrix/Matrix.h` | 1,296 | 1,122 |
| 5 | `mml/core/Fields/FieldOperations.h` | 1,277 | 1,146 |
| 6 | `mml/algorithms/EigenSystemSolvers.h` | 1,254 | 1,076 |
| 7 | `mml/core/Derivation/DerivationRealFunction.h` | 1,170 | 1,089 |
| 8 | `mml/systems/LinearSystem.h` | 1,164 | 1,014 |
| 9 | `mml/tools/Visualizer.h` | 1,159 | 1,053 |
| 10 | `mml/algorithms/Analyzers/FunctionsAnalyzer.h` | 1,140 | 972 |

## Notes

- The single-header distribution accounts for 102,345 physical lines, or about 48.2% of the counted tree.
- Excluding the single-header distribution, the reusable source tree contains 110,128 physical lines across 350 files.
- The largest non-amalgamated areas by recursive physical lines are `mml/algorithms` with 35,415 lines, `mml/core` with 29,150 lines, and `mml/base` with 28,586 lines.