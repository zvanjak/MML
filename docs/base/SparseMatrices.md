# Sparse matrix storage boundaries

The COO builder converts to compressed row (CSR) and column (CSC) storage.
Raw compressed constructors require nonnegative dimensions, matching value/index
counts, a zero first pointer, a final pointer equal to the nonzero count,
nondecreasing pointers within that count, and strictly increasing in-range
indices inside each segment. All pointer ranges are checked before any segment
is dereferenced. Invalid storage raises `ArgumentError`; `validate()` itself is
safe on malformed metadata.

Conversions prepare replacement buffers before committing dimensions and storage.
The CSR mutable raw-buffer accessors remain an advanced interface: callers must
restore valid storage before subsequent computations. CSC conversion from CSR
rejects an invalid source explicitly.

`Multiply`, `MultiplyTranspose` and CSR `gemv` support the same input/output
vector object, for both MML and standard vectors. The input is snapshotted before
output resizing or modification. Nonaliased products avoid that copy.
The one-norm of an empty or zero-column CSR matrix is zero. Sparse dot, axpy and
residual helpers reject incompatible vector sizes before reading or writing data.

CSR `NormFrobenius()` uses a scaled sum-of-squares accumulator, including
separate real and imaginary components for complex entries. It retains
representable norms when squaring the original entries would overflow or
underflow. `Docs_Demo_Matrix_ScaleContracts()` in
[docs_demo_matrix_scale_contracts.cpp](../../src/docs_demos/base/docs_demo_matrix_scale_contracts.cpp)
is the runnable counterpart.

Laplacian factories require nonnegative grid dimensions and finite positive real
spacing. A zero grid extent produces an empty matrix. Grid products and stencil
reservation bounds must fit the signed index range; excessive products raise
`MatrixAllocationError` before allocation. Invalid spacing raises `ArgumentError`;
unrepresentable coefficients raise `MatrixNumericalError`. These checks do not
promise that the operating system can satisfy every otherwise valid allocation.