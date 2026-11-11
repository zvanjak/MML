# MML API Migration Guide

This guide documents breaking API changes in MinimalMathLibrary and provides migration paths for updating your code.

## Scalar namespace cleanup

Global `Real` and `Complex` aliases and the `REAL(...)` macro are removed.
Use `MML::Real`, `MML::Complex`, and explicit scalar conversion, or local
`using MML::Real` / `using MML::Complex` declarations. The former global `GetEnv`
is an internal visualizer helper, not a public API. There is no compatibility
macro or global alias; see [ScalarNamespaceMigration.md](ScalarNamespaceMigration.md)
for precision configuration, collision tests and runnable examples.

## Finite/discrete algebra is optional

The 17 finite/discrete algebra component and algorithm headers moved to
`include/mml_ext/base/Algebra/` in the private MML-Packages-Private repository.
Link `MML_Ext` and change their include prefix from `mml/` to `mml_ext/`.
The extension aggregate is `mml_ext/base/Algebra_base.h`; namespaces, component
guards and mathematical behavior are unchanged. There are no core forwarding
headers, and these optional types are no longer in the core single header.

Core `mml/base/Algebra_base.h` now aggregates only the retained Lie/rigid-motion
types SO2, SO3, SE2 and SE3. Numerical `Polynom.h` also stays core.
The shared generic `Algebra::FieldMatrixDeterminant` / `FieldMatrixInverse`
routines remain core at `mml/algorithms/MatrixAlg/FieldMatrixAlgorithms.h`;
include that owner explicitly when using them. Their arithmetic is unchanged.
See [Algebra](base/Algebra.md) for the boundary and runnable demo selectors.

## Version 2.0 API Changes (January 2026)

### Overview

Version 2.0 standardizes the API naming conventions:
- **Accessors**: lowercase (STL-style)
- **Predicates**: lowercase `is...()` prefix
- **Mutators**: PascalCase
- **Factory methods**: PascalCase static

All deprecated methods remain functional but emit compiler warnings. They will be removed in a future major version.

---

## Version 2.1 update (October 2026) — matrix-family spellings

Version 2.1 adopted a single naming convention for the whole library
([CONTRACT_POLICY.md §11](CONTRACT_POLICY.md#11-naming-convention)): *casing encodes the
abstraction level*. STL-shaped protocol stays lowercase (`rows() cols() size() empty() resize()`),
but **mathematical operations and predicates are PascalCase**, the **verb mutates and the
noun/participle returns a copy**, and accessors are `getX()`. Several 2.0 spellings therefore
changed again; the tables and examples below already show the 2.1 names, and every 2.0 spelling
remains available as a `[[deprecated]]` alias for one minor release.

| 2.0 spelling | 2.1 spelling | Types |
|---|---|---|
| `trace()` | `Trace()` | `Matrix` |
| `transpose()` (copy) | `Transposed()` | `Matrix`, `MatrixNM`; also `TridiagonalMatrix::GetTranspose`, `BandDiagonalMatrix::Transpose` (copy) |
| container `inverse()`, `GetInverse()`, `Inverse()` | `MatrixAlg::Inverse(matrix)` | `Matrix`, `MatrixNM`, `MatrixSym`, `TridiagonalMatrix`; include `algorithms/MatrixAlg/MatrixInversion.h` |
| `isIdentity() isAntiSymmetric() isZero() isSquare()` | `IsIdentity() IsAntiSymmetric() IsZero() IsSquare()` (+ new `IsSymmetric()`) | `Matrix`; `SparseMatrixCSR::isSquare/isSymmetric/hasSymmetricPattern → IsSquare/IsSymmetric/HasSymmetricPattern` |
| `isEmpty()` | `empty()` | `Matrix` |
| `diagonal()`, `GetDiagonal()` | `getDiagonal()` | `Matrix`, `BandDiagonalMatrix` |
| `MakeUnitMatrix()` | `SetIdentity()` (+ new `SetZero()`, `Fill(v)`) | `Matrix`, `MatrixNM` |
| `SetToZero()`, `SetToValue(v)` | `SetZero()`, `Fill(v)` | `MatrixSym` |
| `Resize(...)` | `resize(...)` | `Vector`, `Matrix`, `MatrixSym`, `Matrix3D` |
| `Dim()`, `GetDimension()` | `rows()` | `MatrixSym`, `BandDiagonalMatrix` |
| `GetAsMatrix()` | `ToFullMatrix()` | `MatrixSym` (also new on `TridiagonalMatrix`) |
| `GetDiagonalMatrix(v)` | `Diagonal(v)` | `MatrixSym` |
| `GetLowerBandwidth() GetUpperBandwidth()` | `getLowerBandwidth() getUpperBandwidth()` | `BandDiagonalMatrix`, `BandDiagLUSolver` (+ `getDimension()`) |
| `norm1() normInf() normFrobenius()` | `Norm1() NormInf() NormFrobenius()` | `SparseMatrixCSR` |
| `multiply multiplyTranspose solveLower solveUpper zeroRow print toDense fromCOO` | `Multiply MultiplyTranspose SolveLower SolveUpper ZeroRow Print ToDense FromCOO` | `SparseMatrixCSR`, `SparseMatrixCSC` (BLAS-named `gemv`, `axpy`, `dot` unchanged) |
| `laplacian1D/2D/3D identityCSR diagonalCSR` | `Laplacian1D/2D/3D IdentityCSR DiagonalCSR` | `MML::SparseMatrix` factories |

| `isZero() isNearZero(eps)` | `IsZero() IsNearZero(eps)` | `Vector`, `VectorN` (and all `Vector2/3/4…` types) |
| `Clear()` | `clear()` | `Vector`, `Polygon2D` (remove all elements) |
| `clear()` on `VectorN` | `SetZero()` (+ new `Fill(v)`) | `VectorN` — fixed size: a value operation, not `clear()` |
| `GetAsUnitVector()`, `GetAsUnitVectorOrZero()` | `Normalized()`, `NormalizedOrZero()` | `Vector2Cartesian`, `Vector3Cartesian` (`GetAsUnitVectorAtPos` unchanged) |
| `getAsPoint()` | `ToPoint()` | `Vector3Cartesian` |
| `getPerpendicularVectors` | `GetPerpendicularVectors` | `Vector2Cartesian` (now same as 3D) |
| `isTimelike() isSpacelike() isLightlike()` | `IsTimelike() IsSpacelike() IsLightlike()` | `Vector4Minkowski` |
| `isUnit() isIdentity() isApprox(q)` | `IsUnit() IsIdentity() IsEqualTo(q)` | `Quaternion` |
| `GetRotationAngle() GetRotationAxis()` | `getRotationAngle() getRotationAxis()` | `Quaternion`, `CoordTransf3DRotationQuaternion` |
| `At(i,j,k)`, `isEmpty()` | `at(i,j,k)`, `empty()` | `Matrix3D` |
| `GetCenter GetBoundingBox GetBoundingSphere` | `getCenter getBoundingBox getBoundingSphere` | `IBody` and every body; `Polygon2D`/`Circle2D::getBoundingBox` |
| `GetRadius GetSide GetHeight GetBaseCenter GetBaseSize GetMajorRadius GetMinorRadius GetNum*` | `getRadius getSide getHeight …` | `Cube3D`, `Sphere3D`, `Cylinder3D`, `Torus3D`, `Pyramid3D`… |
| `GetSurfaceCount GetSurface GetSolidCount GetSolid` | `getSurfaceCount getSurface getSolidCount getSolid` | `BodyWithTriangleSurfaces`, `BodyWithRectSurfaces`, `ComposedSolidSurfaces3D` |
| `GetPointOnPlane GetCoordAxisSegments` | `getPointOnPlane getCoordAxisSegments` | `Plane3D` |
| `GetXYPlane() GetXZPlane() GetYZPlane()` | `XYPlane() XZPlane() YZPlane()` | `Plane3D` factories |

**Behaviour change (2.1):** `Vector3Cartesian::IsParallelTo(b)` now uses ‖a × b‖ ≤ eps·‖a‖·‖b‖ —
it is scale-invariant, **includes antiparallel vectors** (parallel planes/lines regardless of
orientation), and returns `false` for the zero vector instead of comparing NaNs. The old
component-wise unit-vector difference silently rejected antiparallel directions.

**Interface vocabulary (2.1):** everything declared on an `mml/interfaces` type is lowerCamel, and
implementers override the lowerCamel name. The PascalCase spelling stays on the interface as a
`[[deprecated]]` *non-virtual* forwarder — callers keep working, but a derived class that still
overrides the old name fails to compile (`override` of a non-virtual) and must be renamed.

| 2.0 spelling | 2.1 spelling | Interface |
|---|---|---|
| `NumContravar() NumCovar() IsContravar(i) Contract(...)` | `numContravar() numCovar() isContravar(i) contract(...)` | `ITensor2..5` (implementers `Tensor2..5`, also `IsCovar → isCovar`) |
| `Component(i, j, …, pos)` | `component(i, j, …, pos)` | `ITensorField2..5` (metric tensors, `InducedMetric`, tensor-field adapters, your own fields) |
| `TransformsAsTensor()` | `transformsAsTensor()` | `ITensorField3` |
| `GetEquidistantCovering(n, pts)` | `getEquidistantCovering(n, pts)` | `IInterval` and all intervals |
| `GetValues(x1, x2, n, xs, ys)` | `getValues(…)` | `IRealFunction` |
| `GetTrace(t1, t2, n)` | `getTrace(…)` | `IParametricCurve<N>` |

**Derivation entry points (2.1):** numerical, forward-AD and reverse-AD Jacobians / Hessians /
gradients share one `Calc<Object>[<Method>]` scheme; argument types select fixed vs dynamic size
(see [Derivation.md](core/Derivation.md#jacobians-hessians-and-gradients--one-vocabulary)).

| 2.0 spelling | 2.1 spelling |
|---|---|
| `Derivation::calcJacobian`, `calcJacobianDyn` | `Derivation::CalcJacobian` (overloads: `IVectorFunction`+`VectorN` → `MatrixNM`; `std::function`+`Vector` → `Matrix`) |
| `calcJacobianDynInPlace`, `calcJacobianDynColumns`, `calcJacobianDyn2` | `CalcJacobianInPlace`, `CalcJacobianColumns`, `CalcJacobianOrder2` |
| `calcJacobianODE`, `calcJacobianMap`, `calcHessian` | `CalcJacobianODE`, `CalcJacobianMap`, `CalcHessian` |
| `AD::jacobianForwardAD`, `AD::jacobianReverseAD`, `AD::hessianForwardAD` | `AD::CalcJacobianForwardAD`, `AD::CalcJacobianReverseAD`, `AD::CalcHessianForwardAD` |
| `AD::gradient`, `AD::gradientReverse`, `AD::derivative` | `AD::CalcGradientForwardAD`, `AD::CalcGradientReverseAD`, `AD::CalcDerivativeForwardAD` |

**Honest names and the last cross-cutting spellings (2.1):**

| 2.0 spelling | 2.1 spelling | Note |
|---|---|---|
| `FindRootX(f, x1, x2, config)` → `RootFindingResult` | `FindRootXDetailed(f, x1, x2, config)`; new `FindRootXDetailed(f, x1, x2, xacc)` | X ∈ Bisection, FalsePosition, Newton (incl. analytic-derivative form), Secant, Ridders, Brent. Plain `FindRootX(f, x1, x2, xacc)` still returns `Real` and throws. |
| `Matrix::isAntiSymmetric()` | `IsAntiSymmetric()` (plain, no conjugation) **or** `IsSkewHermitian()` (conjugated); new `IsHermitian()` | The 2.0 spelling conjugated for complex types, i.e. tested skew-Hermitian; the alias keeps that meaning. |
| `ToString()` / `ToString(enum)` | `to_string()` / `to_string(enum)` | `IBody` and all bodies, `Box3D`, `BoundingSphere3D`, `Matrix3D` (`ToStringFull(p)` → `to_string(p)`), `LinearProgram`, `AlgorithmStatus`, `GraphAlgorithmStatus`, `FixedPointType` |
| `GetMatrix()` | `getMatrix()` (reference accessor) / `Tensor2::ToMatrix()` (copy conversion) | `LinearSystem`, `MatrixAnalyzer`, `Tensor2` |
| `RealFunctionAnalyzer::isMonotonic/isContinuous/isLocalOptimum/…` | `IsMonotonicOnSamples/AppearsContinuous/IsLocalOptimum/…` | Documented as sampling heuristics, not proofs. |

Re-documented without renaming: `FindRootIntervals` (sampling-based bracketing, not algebraic root
isolation), `LineMinimizer` (exact Brent line minimization, not a Wolfe search), `RealFFT` (half-spectrum
storage via a half-size complex FFT, not NR's in-place packed `realft`).

---

## Quick Migration Reference

### Matrix Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `RowNum()` | `rows()` | Accessor |
| `ColNum()` | `cols()` | Accessor |
| `IsEmpty()` | `empty()` | Predicate (STL-compatible) |
| `IsUnit(eps)` | `IsIdentity(eps)` | Predicate |
| `IsDiagonal(eps)` | `isDiagonal(eps)` | Predicate |
| `IsDiagDominant()` | `isDiagonallyDominant()` | Predicate |
| `IsSymmetric()` | `IsSymmetric()` | Predicate |
| `IsAntiSymmetric()` | `IsAntiSymmetric()` | Predicate |
| `IsZero(eps)` | `IsZero(eps)` | Predicate |
| `Identity(n)` | `Identity(n)` | Factory |
| `Diagonal(v)` | `Diagonal(v)` | Factory |
| `identity(n)` | `Identity(n)` | Factory |
| `transpose()` | `Transposed()` | Operation (returns copy) |
| matrix `GetInverse()` | `MatrixAlg::Inverse(matrix)` | Algorithm operation (returns copy; container method removed) |
| `Trace()` | `Trace()` | Operation |
| `GetDiagonal()` | `getDiagonal()` | Accessor |
| `VectorFromDiagonal()` | `getDiagonal()` | Accessor |
| `toString()` | `to_string()` | String conversion |

### Vector Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `isEmpty()` ✅ | `empty()` | Predicate (STL-compatible) |
| `IsNullVec()` | `isZero()` | Predicate |
| `toString()` | `to_string()` | String conversion |

### VectorN Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `IsNullVec()` | `isZero()` | Predicate |

### Polynom Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `GetDegree()` | `degree()` | Accessor |
| `IsNullPolynom()` | `isNull()` | Predicate |
| `Derive()` | `derivative()` | Accessor (returns copy) |
| `Integrate()` | `integral()` | Accessor (returns copy) |
| `toString()` | `to_string()` | String conversion |

---

## Migration Examples

### Matrix Dimension Access

```cpp
// ❌ Old code
Matrix<double> A(3, 4);
int r = A.RowNum();
int c = A.ColNum();

// ✅ New code
Matrix<double> A(3, 4);
int r = A.rows();
int c = A.cols();
```

### Matrix Factory Methods

```cpp
// ❌ Old code
auto I = Matrix<double>::Identity(3);
auto D = Matrix<double>::Diagonal(vec);

// ✅ New code
auto I = Matrix<double>::Identity(3);
auto D = Matrix<double>::Diagonal(vec);
```

### Matrix Properties

```cpp
// ❌ Old code
if (A.IsUnit()) { ... }
if (A.IsSymmetric() && !A.IsEmpty()) { ... }

// ✅ New code (2.1 spelling: mathematical predicates PascalCase, STL protocol lowercase)
if (A.IsIdentity()) { ... }
if (A.IsSymmetric() && !A.empty()) { ... }
```

### Matrix Operations (Copy vs In-Place)

```cpp
// Copy operation (noun / participle) - returns new matrix
Matrix<double> B = A.Transposed();  // A unchanged, B is transpose
Matrix<double> C = MatrixAlg::Inverse(A);     // A unchanged, C is inverse

// In-place operation (verb) - modifies the matrix
A.Transpose();  // A is now transposed
MatrixAlg::Invert(A);     // A is now its inverse
```

### Polynomial Operations

```cpp
// ❌ Old code
int deg = p.GetDegree();
auto dp = p.Derive();
auto ip = p.Integrate();

// ✅ New code
int deg = p.degree();
auto dp = p.derivative();
auto ip = p.integral();
```

### Vector Zero Check

```cpp
// ❌ Old code
Vector<double> v = {0, 0, 0};
if (v.IsNullVec()) { ... }

// ✅ New code
Vector<double> v = {0, 0, 0};
if (v.isZero()) { ... }
```

### String Conversion

```cpp
// ❌ Old code
std::cout << A.toString() << std::endl;
std::cout << v.toString() << std::endl;

// ✅ New code
std::cout << A.to_string() << std::endl;
std::cout << v.to_string() << std::endl;

// With formatting
std::cout << A.to_string(12, 6) << std::endl;  // width=12, precision=6
```

---

## Automated Migration

### Using Compiler Warnings

Enable deprecation warnings to find all usages:

```cpp
// GCC/Clang
-Wdeprecated-declarations

// MSVC
/W4  // Level 4 warnings include deprecation
```

### Search and Replace Patterns

For bulk migration, use these regex patterns:

```regex
# Matrix accessors
\.RowNum\(\)  →  .rows()
\.ColNum\(\)  →  .cols()

# Matrix factories
::Identity\(  →  ::Identity(
::GetDiagonalMatrix\(  →  ::Diagonal(
::identity\(  →  ::Identity(

# Matrix predicates
\.IsUnit\(  →  .IsIdentity(
\.IsEmpty\(  →  .empty(
\.IsSymmetric\(  →  .IsSymmetric(
\.IsDiagonal\(  →  .isDiagonal(
\.IsDiagDominant\(  →  .isDiagonallyDominant(
\.IsZero\(  →  .IsZero(

# Matrix operations
\.transpose\(  →  .Transposed(
Matrix inversion aliases are removed; use MatrixAlg::Inverse(matrix).
\.Trace\(  →  .Trace(
\.GetDiagonal\(  →  .getDiagonal(

# Polynom
\.GetDegree\(  →  .degree(
\.Derive\(  →  .derivative(
\.Integrate\(  →  .integral(
\.IsNullPolynom\(  →  .isNull(

# Vector
\.IsNullVec\(  →  .isZero(
\.toString\(  →  .to_string(
```

---

## Deprecation Timeline

| Version | Status |
|---------|--------|
| 2.0 | Deprecated methods emit `[[deprecated]]` warnings |
| 2.x | Deprecated methods continue to work |
| 3.0 | Deprecated methods will be **removed** |

**Recommendation**: Migrate your code now while deprecated methods still work. This allows gradual migration without breaking builds.

---

## noexcept Additions

The following methods now have `noexcept` specifier for better optimization:

### Matrix
- `rows()`, `cols()`, `isEmpty()`, `isSquare()`
- `data()` (both const and non-const)

### Vector
- `size()`, `empty()`, `isEmpty()`, `isZero()`
- `front()`, `back()` (both const and non-const)
- All iterators

### VectorN
- `size()`, `clear()`, `isZero()`
- `operator[]` (both const and non-const)

### Polynom
- `degree()`, `isNull()`, `leadingTerm()`, `constantTerm()`
- `operator[]` (both const and non-const)

### Graph
- `numVertices()`, `isEmpty()`, `allowsSelfLoops()`

---

## Rationale

### Why lowercase accessors?
- Consistency with C++ STL (`size()`, `empty()`, `begin()`)
- Signals that the method is a simple getter
- `noexcept` guarantee is natural for simple accessors

### Why `Is...()` predicates? (2.1)
- Clear intent: returns boolean
- Mathematical/geometric predicates are PascalCase like the other mathematical operations
  (`IsSymmetric()`, `IsZero()`, `Contains(p)`); STL protocol keeps `empty()`; interface
  virtuals keep lowerCamel (`isAutonomous()`)
- Readable: `if (matrix.IsSymmetric())`

### Why PascalCase mutators?
- Visual distinction from accessors
- Signals "this method modifies the object"
- Common in C++ libraries for actions

### Why `to_string()` instead of `toString()`?
- C++ convention (see `std::to_string`)
- Consistency with STL and standard library

---

## Questions?

See [CONTRACT_POLICY.md §11](CONTRACT_POLICY.md#11-naming-convention) for the complete naming convention reference.

Matrix inversion now belongs to algorithms, not base containers; see
[the inversion and formatting migration](algorithms/MatrixInversionAndFormatting.md).
