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

`solveLower` and `solveUpper` require a square matrix, a matching right-hand
side, and an explicit non-zero diagonal in every row. Missing or zero diagonals
throw `SingularMatrixError`. SSOR and ILU(0) perform equivalent diagonal checks
during setup, before an unsafe preconditioner can be applied.

The 2D five-point COO helper takes both grid dimensions:

```cpp
coo.addLaplacian2D_5pt(row, i, j, nx, ny, inverseHSquared);
```

It emits only neighbors inside the grid, including at corners and edges.