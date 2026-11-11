# Matrix Inversion and Shared Numeric Formatting

## Inversion Is an Algorithm

`Matrix`, `MatrixNM`, `MatrixSym`, and `TridiagonalMatrix` no longer expose
`Inverse`, `Invert`, or their legacy inversion aliases. There are no forwarding
shims. Include the algorithm header explicitly:

```cpp
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>

Matrix<Real> matrix(2, 2, {4, 7, 2, 6});
const auto inverse = MatrixAlg::Inverse(matrix);
MatrixAlg::Invert(matrix);
MatrixNM<Real, 2, 2> fixed{2, 0, 0, 4};
const auto fixedInverse = MatrixAlg::Inverse(fixed);
```

The dynamic inverse is a dynamic matrix, the fixed inverse retains its compile-time
shape, and symmetric/tridiagonal inverses are generally dense dynamic matrices.
All four numerical paths use the same `LUSolver` factorization/solve route.
`Invert` is available for dynamic and fixed dense matrices and commits the new
value only after successful inversion; failures leave the source unchanged.
Empty/nonsquare input is rejected with `MatrixDimensionError`, singular input
with `SingularMatrixError`, and nonfinite numerical input with `MatrixNumericalError`.

`MatrixAlg::Inverse`, `MatrixAnalyzer::Inverse` and RHS-aware `LinearSystem::Inverse`
remain algorithm/facade operations. Specialized solver operations (such as an
already-factored Cholesky inverse) and Gauss-Jordan solve routines remain distinct
algorithms; they are not additional container inversion implementations. Exact
finite-field algebra keeps its field arithmetic algorithms rather than routing
through floating-point LU.

### Inverse-Dependent Geometry

The typed `Metric` implementation is now in
`mml/core/DifferentialGeometry/Metric.h`, not the base value layer. Base Hodge
templates refer to the metric type without pulling in inversion algorithms.
Include the core metric header when constructing a metric.

```cpp
#include <mml/core/DifferentialGeometry/Metric.h>
#include <mml/base/DifferentialGeometry/FrameTags.h>

const Metric<2, Cartesian2> metric(fixed);
```

The singularity-policy inverse helper is in `mml/algorithms/MatrixSingularity.h`.
It preserves the existing policies and catches only singular-matrix failures;
invalid/nonfinite inputs are not silently replaced by sentinels.

## Shared Formatting, Existing Public Print APIs

`mml/base/PrintFormatting.h` provides the dependency-light rendering machinery
used by vector/matrix containers and tools-side numeric columns. Base containers
do not include `ConsolePrinter` or core inversion solvers.

`Print`, `operator<<`, and `to_string` remain available. Tools retain table
alignment, borders and export layout; only numeric value rendering and stream
state handling are shared.

```cpp
MatrixPrintFormat format = MatrixPrintFormat::Compact();
format.width = 0;
format.zeroThreshold = Precision::GeometryEpsilon;
std::ostringstream rendered;
matrix.Print(rendered, format);
```

Dynamic/fixed/symmetric/tridiagonal matrices consistently honor header/bracket
selection, delimiter, field width, precision, scientific/fixed notation, compact
mode for matrices no larger than 3 by 3, and zero threshold. This intentionally
fixes dynamic `Matrix` ignoring `compactMode` and format-based `zeroThreshold`,
and structured matrices ignoring thresholding.

Threshold replacement uses `|value| <= threshold` when the threshold is positive.
Output applies the threshold while rendering, without allocating a thresholded
matrix copy. Nonfinite values remain visible instead of becoming successful-looking
zeroes. Vector row layouts remain bracketed and comma-separated.

Printing restores caller flags, precision, width and fill on normal return or
exceptions and preserves locale. Caller locale/notation inheritance remains the
existing display policy; this is not persistence's canonical text format.

Runnable counterparts: [docs_demo_matrix_consolidation.cpp](../../src/docs_demos/algorithms/docs_demo_matrix_consolidation.cpp).
Run `MML_DocsApp --matrix-consolidation` for inversion, relocated metrics, compact
output and table numeric rendering.
