# Sparse Matrices and Solvers

MML provides coordinate (`SparseMatrixCOO`), compressed-row (`SparseMatrixCSR`),
and compressed-column (`SparseMatrixCSC`) storage. Build incrementally in COO,
then convert to CSR for row-oriented solvers or CSC for column-oriented work.

## Structural validation

Each format exposes `validate()`. Raw CSR and CSC constructors call it and throw
`ArgumentError` when compressed storage is malformed. Validation checks:

- non-negative dimensions and matching value/index counts;
- pointer arrays of the required size, beginning at zero and ending at `nnz`;
- monotone pointers and in-bounds indices;
- strictly sorted, duplicate-free column indices per CSR row or row indices per
	CSC column;
- in-bounds COO triplets and strict ordering after consolidation.

CSR exposes mutable storage for advanced assembly. Call `validate()` after direct
mutation and before using the matrix.

## Vector operations

CSR and CSC matrix-vector products support both `std::vector<T>` and
`MML::Vector<T>`. CSR transpose products, `gemv`, and triangular solves provide
the same pair of interfaces.

## Safe solves and preconditioners

`SolveLower` and `SolveUpper` require a square matrix, a matching right-hand
side, and an explicit non-zero diagonal in every row. Missing or zero diagonals
throw `SingularMatrixError`. SSOR and ILU(0) perform equivalent diagonal checks
during setup, before an unsafe preconditioner can be applied.

### Parameters, ownership, and aliasing

The Krylov family (`SolverConfig`, `KrylovResult`, the solver base, CG, BiCGSTAB,
GMRES, AutoSolver, and their convenience functions) requires `MMLReal<T>`:
`float`, `double`, or `long double`. Complex/integer specializations are rejected
at the public constraint, rather than failing inside real-only arithmetic.
SSOR and block-Jacobi likewise declare their real-scalar domain.

- **Block Jacobi:** block size must be positive; even `INT_MAX` is a valid request
  for one block when the matrix is smaller. Partition arithmetic does not form
  an overflowing `n + blockSize - 1`. Setup requires a finite, nonempty square
  matrix with invertible diagonal blocks. Singular blocks now throw
  `SingularMatrixError` instead of silently installing an identity block.
  Pivot checks are relative to each block's scale; unrepresentable factors or
  applications throw numerical errors.
- **SSOR:** omega must be finite and strictly in `(0, 2)`, including constructor
  and setter calls. Setup owns a sparse-matrix snapshot and its diagonals.
  Temporaries are safe; source destruction, mutation, or resizing cannot alter
  the installed preconditioner. This costs `O(nnz + n)` owned storage. Copies
  own independent snapshots. The existing `(2-omega)`-scaled SSOR convention
  is unchanged.
- Successful setup atomically replaces SSOR/block-Jacobi state. Failed setup or
  an invalid setter preserves the previous installation. A valid block-size
  setter affects the next setup; a valid omega setter immediately changes the
  SSOR sweeps (the cached matrix/diagonals do not depend on omega).
- `apply(r, r)` is supported. Block Jacobi snapshots an aliased residual before
  writing output; separate input/output vectors do not incur that copy.
- Krylov `solve(A, b, x)` supports `b` and `x` being the same vector. An internal
  snapshot preserves the original mathematical RHS through iteration and true
  residual finalization. The shared caller vector naturally contains `x` on
  return; save a copy beforehand if the caller also needs the original RHS.
  Nonaliased calls do not copy the RHS or returned solution.

An in-place block application with a known inverse action:

```cpp
#include <mml/core/SparseSolvers/IterativeSolvers.h>
using namespace MML;
using namespace MML::SparseSolvers;

SparseMatrix::SparseMatrixCSR<Real> matrix(
    2, 2, {4, 1, 1, 3}, {0, 1, 0, 1}, {0, 2, 4});
BlockJacobiPreconditioner<Real> block(2);
block.setup(matrix);
Vector<Real> residual{1, 3};
block.apply(residual, residual); // (0, 1).
```

### Zero-guess result payload

`solver.solve(A, b)` now moves its computed vector into `result.solution`,
including a partial iterate on non-convergence. Check the result status before
using it. The explicit-output `solve(A, b, x)` keeps the solution only in `x`;
its `result.solution` is empty, avoiding a redundant allocation/copy.

```cpp
ConjugateGradient<Real> solver;
Vector<Real> rhs{1, 2};
auto result = solver.solve(matrix, rhs);
if (!result.converged)
    throw ConvergenceError(result.error_message, result.iterations_used, result.residualNorm);
const Vector<Real>& solution = result.solution;
```

### Automatic solver selection

AutoSolver selects CG only after symmetry plus a sufficient diagonal-dominance
test: every diagonal is positive, every row is weakly dominant, and **each
connected component** has a strict row. A strict row in an unrelated block does
not make a singular component positive definite. Connected Dirichlet Laplacians
and disconnected individually anchored blocks still qualify. The test uses
`O(n + nnz)` work and `O(n)` auxiliary storage; matrices outside this sufficient
class go to BiCGSTAB, not to a claim that they are indefinite.

The runnable counterpart is `MML_DocsApp --solver-boundaries` in
`src/docs_demos/core/docs_demo_lin_alg_solvers.cpp`.
`MML_SolverBoundaryConsumerChecks` independently measures heap allocations in
the alias/nonalias, owned-result, and SVD paths without instrumenting the main test executable.

## Sparse assembly

The 2D five-point COO helper takes both grid dimensions:

```cpp
coo.addLaplacian2D_5pt(row, i, j, nx, ny, inverseHSquared);
```

It emits only neighbors inside the grid, including at corners and edges.