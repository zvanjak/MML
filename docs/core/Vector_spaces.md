# Vector Spaces

This document fixes the first MML 2.0 design boundary for finite-dimensional
vector spaces. It turns the roadmap in
[Vector_Spaces_Improvement.md](../improvements/Vector_Spaces_Improvement.md)
into a concrete API direction before implementation begins.

The goal is not to replace `Vector`, `VectorN`, `Matrix`, or `MatrixNM`. Those
remain the coordinate-storage and numerical-kernel types. The vector-space layer
adds mathematical meaning on top of that storage: spaces, bases, subspaces,
linear maps, dual spaces, inner products, and affine spaces.

---

## Design Decision Summary

| Topic | Decision |
|-------|----------|
| Folder | `mml/core/VectorSpaces/` for implementation headers |
| Umbrella header | `mml/core/VectorSpaces.h` |
| Namespace | `MML::VectorSpaces` for new APIs |
| Fixed-size types | Preferred for geometry, physics, typed forms, and small spaces |
| Dynamic types | Used for SVD/null-space/range workflows and matrix factorization outputs |
| Storage | Existing `VectorN`, `Vector`, `MatrixNM`, `Matrix` |
| Runtime identity | Spaces have optional names/tags; dimension alone is not semantic identity |
| First implementation | Completed as compact fixed-size/dynamic semantic wrappers over existing storage |
| Matrix algorithms | Domain objects use core kernels; public analysis is owned by `MatrixAlg` and `MatrixAnalyzer` |

---

## Proposed Header Layout

```text
mml/core/VectorSpaces.h
mml/core/VectorSpaces/
  VectorSpace.h          // VectorSpace, VectorSpaceX, VectorInSpace
  Basis.h                // Basis, BasisX, ChangeBasis
  Subspace.h             // dynamic Subspace first; fixed-size later
  LinearMap.h            // LinearMap and LinearMapX
  DualSpace.h            // DualSpace, CovectorInSpace
  InnerProductSpace.h    // inner products, projections, adjoints
  AffineSpace.h          // AffineSpace, PointInSpace, AffineMap
```

The umbrella header should include only stable public headers. During early
implementation, new headers may be included directly by tests until their APIs
settle.

---

## Namespace and Naming

New types live under `MML::VectorSpaces`:

```cpp
namespace MML::VectorSpaces
{
    template<class Scalar, int N> class VectorInSpace;
    template<class Scalar, int N> class Basis;

    template<class Scalar> class VectorSpaceX;
    template<class Scalar> class BasisX;
    template<class Scalar> class Subspace;

    template<class Scalar, int DomainN, int CodomainN> class LinearMap;
    template<class Scalar, int N> class DualSpace;
    template<class Scalar, int N> class CovectorInSpace;
    template<class Scalar, int N> class InnerProductSpace;
    template<class Scalar, int N> class AffineSpace;
    template<class Scalar, int N> class PointInSpace;
    template<class Scalar, int N> class AffineFrame;
    template<class Scalar, int DomainN, int CodomainN> class AffineMap;
}
```

The names intentionally mirror mathematical language. Dynamic-size siblings use
the `X` suffix already familiar from many numerical libraries:

- `VectorSpace<Real, 3>`: fixed-size space.
- `VectorSpaceX<Real>`: runtime-dimension space.
- `Basis<Real, 3>`: fixed-size basis.
- `BasisX<Real>`: runtime-dimension basis.
- `LinearMap<Real, 3, 2>`: fixed-size linear map.
- `LinearMapX<Real>`: future dynamic-size linear map, if needed.

The first completed implementation intentionally skips a heavyweight
`VectorSpace` identity object. Coordinates remain stored in `VectorN`/`Vector`,
while `Basis`, `Subspace`, `LinearMap`, `DualSpace`, `InnerProductSpace`, and
`AffineSpace` provide the semantic operations around that storage.

---

## Relationship to Existing Types

| Existing type | New role |
|---------------|----------|
| `VectorN<Scalar,N>` | coordinate tuple in a fixed-size basis |
| `Vector<Scalar>` | coordinate tuple in a dynamic-size basis/subspace |
| `MatrixNM<Scalar,R,C>` | fixed-size matrix representation of a map or basis |
| `Matrix<Scalar>` | dynamic matrix representation of a map, basis, or subspace |
| `Point<Real,N,Frame>` | affine/geometric point; should bridge to `PointInSpace` later |
| `TangentVector<N,Frame>` | rank-1 typed vector; should bridge to `VectorInSpace` later |
| `Covector<N,Frame>` / `OneForm<N,Frame>` | rank-1 typed dual object; should bridge to `CovectorInSpace` later |
| `Metric<N,Frame>` | finite metric value object; should bridge to `InnerProductSpace` |

The new layer must avoid implicit conversion from every `VectorN` to every
space. A coordinate tuple gets meaning only when paired with a basis or a space.

---

## Runtime Identity Model

Dimension alone is not enough to identify a vector space. Two unrelated 3D
spaces may both have coordinate tuples of length 3 but should not automatically
mix.

The initial model:

- fixed-size `VectorSpace<Scalar,N>` carries dimension and an optional name;
- dynamic `VectorSpaceX<Scalar>` carries runtime dimension and optional name;
- equality is structural by dimension and scalar type unless a stronger identity
  token is provided later;
- typed geometry still uses compile-time frame tags for strict category safety;
- bridge APIs to typed forms should require explicit frame/tag choices.

This keeps the first implementation practical. If later APIs need stronger
runtime identity, a `SpaceId` value can be added without changing coordinate
storage.

---

## Fixed-Size API

```cpp
template<class Scalar, int N>
class VectorInSpace
{
public:
    const VectorSpace<Scalar, N>& space() const noexcept;
    const Basis<Scalar, N>& basis() const noexcept;
    const VectorN<Scalar, N>& coordinates() const noexcept;

    VectorN<Scalar, N> coordinatesIn(const Basis<Scalar, N>& targetBasis) const;
};

template<class Scalar, int N>
class Basis
{
public:
    static Basis Standard();

    explicit Basis(const MatrixNM<Scalar, N, N>& columns);

    const MatrixNM<Scalar, N, N>& matrix() const noexcept;
    const MatrixNM<Scalar, N, N>& inverseMatrix() const noexcept;
    VectorN<Scalar, N> coordinatesInStandard(const VectorN<Scalar, N>& coords) const;
    VectorN<Scalar, N> coordinatesFromStandard(const VectorN<Scalar, N>& vector) const;
};
```

Free helper:

```cpp
template<class Scalar, int N>
VectorN<Scalar, N> ChangeBasis(const VectorN<Scalar, N>& coords,
                               const Basis<Scalar, N>& fromBasis,
                               const Basis<Scalar, N>& toBasis);
```

---

## Dynamic API Skeleton

The first dynamic type to implement should be `Subspace<Real>`.

```cpp
template<class Scalar>
class Subspace
{
public:
    Subspace(int ambientDimension, const Matrix<Scalar>& orthonormalBasisColumns);

    int dimension() const noexcept;
    int ambientDimension() const noexcept;
    const Matrix<Scalar>& orthonormalBasis() const noexcept;

    bool contains(const Vector<Scalar>& vector, Real tolerance) const;
    Vector<Scalar> project(const Vector<Scalar>& vector) const;
    Vector<Scalar> residual(const Vector<Scalar>& vector) const;
    Real distance(const Vector<Scalar>& vector) const;

    Subspace orthogonalComplement(Real tolerance) const;

    static Subspace NullSpace(const Matrix<Scalar>& matrix, Real tolerance);
    static Subspace ColumnSpace(const Matrix<Scalar>& matrix, Real tolerance);
    static Subspace RowSpace(const Matrix<Scalar>& matrix, Real tolerance);
    static Subspace LeftNullSpace(const Matrix<Scalar>& matrix, Real tolerance);
};
```

The dynamic implementation delegates basis extraction to the core SVD kernel;
general matrix-analysis consumers use `MatrixAlg::FundamentalSubspacesOf`.

---

## Linear Maps

The fixed-size `LinearMap` distinguishes a matrix as a representation from the
mathematical map it represents.

C++20 constraints keep square-only operations out of incompatible overload sets:
`Identity()`, `inverse()`, and affine translations are available only when domain
and codomain dimensions match, while composition encodes the intermediate
dimension in the function signature.

```cpp
template<class Scalar, int DomainN, int CodomainN>
class LinearMap
{
public:
    explicit LinearMap(const MatrixNM<Scalar, CodomainN, DomainN>& matrix);

    const MatrixNM<Scalar, CodomainN, DomainN>& matrix() const noexcept;
    VectorN<Scalar, CodomainN> apply(const VectorN<Scalar, DomainN>& x) const;
    VectorN<Scalar, CodomainN> operator()(const VectorN<Scalar, DomainN>& x) const;

    int rank(Real tolerance) const;
    int nullity(Real tolerance) const;
};

template<class Scalar, int A, int B, int C>
LinearMap<Scalar, A, C> Compose(const LinearMap<Scalar, B, C>& g,
                                const LinearMap<Scalar, A, B>& f);
```

Kernel and image should return `Subspace<Scalar>` initially. Fixed-size subspace
types can come later if there is enough value.

---

## Dual and Inner-Product Bridge

The vector-spaces layer should not duplicate the typed forms layer. It should
bridge to it.

Implemented bridge:

- `DualSpace<Scalar,N>` owns dual basis construction and covectors-in-space.
- `CovectorInSpace` evaluates `VectorInSpace`.
- Pullback by `LinearMap` lives here.
- `Covector<N,Frame>` and `OneForm<N,Frame>` conversions require explicit frame
  tags or basis/frame adapters.
- `InnerProductSpace<Scalar,N>` owns a positive-definite Gram matrix in standard
  coordinates, with `Metric<N,Frame>`-backed construction for Euclidean or
  Riemannian metrics.
- Inner-product operations include `inner`, `norm`, `distance`, Gram matrix
  representation in another basis, Gram-Schmidt orthonormalization, projection
  onto `Subspace`, Gram-orthogonal complements, and adjoints of `LinearMap`.
- Indefinite metrics, such as Lorentzian/Minkowski signatures, remain in the
  typed `Metric` layer for now. They are intentionally rejected by
  `InnerProductSpace` because norms, distances, projections, and Cholesky-style
  positive-definiteness assumptions are not valid in the same way.

---

## Affine Spaces

Affine spaces reuse the point/vector distinction already proven by typed forms:

```cpp
template<class Scalar, int N>
class AffineSpace;

template<class Scalar, int N>
class PointInSpace;

template<class Scalar, int N>
class AffineFrame;

template<class Scalar, int DomainN, int CodomainN>
class AffineMap;
```

Rules:

- point minus point gives a vector;
- point plus vector gives a point;
- point plus point is invalid;
- `AffineFrame` stores an origin plus a vector-space `Basis`;
- affine maps contain a linear part and a translation part;
- affine map composition follows the same `Compose(g, f)` convention as
  `LinearMap`, meaning `g` after `f`;
- square affine maps expose `inverse()` when their linear part is invertible;
- homogeneous-coordinate helpers are explicit.

Example:

```cpp
using namespace MML::VectorSpaces;

PointInSpace<Real, 2> p(VectorN<Real, 2>{1.0, 2.0});
PointInSpace<Real, 2> q(VectorN<Real, 2>{4.0, 6.0});

VectorN<Real, 2> displacement = q - p;  // vector
PointInSpace<Real, 2> moved = p + displacement;

AffineFrame<Real, 2> frame(
  PointInSpace<Real, 2>(VectorN<Real, 2>{10.0, 20.0}),
  Basis<Real, 2>({
    VectorN<Real, 2>{2.0, 0.0},
    VectorN<Real, 2>{0.0, 3.0}
  }));

PointInSpace<Real, 2> point = frame.pointFromCoordinates({1.0, 1.0});

auto translate = AffineMap<Real, 2, 2>::Translation({5.0, 6.0});
PointInSpace<Real, 2> translated = translate(point);
```

`PointInSpace` can explicitly bridge to the existing typed `Point` layer:

```cpp
Point<Real, 2, Cartesian2> typed = point.toPoint<Cartesian2>();
auto semantic = PointInSpace<Real, 2>::fromPoint(typed);
```

Homogeneous-coordinate helpers keep point and vector semantics explicit:

```cpp
VectorN<Real, 3> hp = HomogeneousPoint(point);   // last component = 1
VectorN<Real, 3> hv = HomogeneousVector(vector); // last component = 0
```

There is deliberately no `point + point` overload.

---

## Examples

The buildable example [src/examples/11_vector_spaces/main.cpp](../../src/examples/11_vector_spaces/main.cpp)
demonstrates basis changes, subspace projection, linear maps, dual evaluation,
weighted inner products, affine frames, and translations.

Key usage snippets:

```cpp
Basis<Real, 2> scaled({
  VectorN<Real, 2>{2.0, 0.0},
  VectorN<Real, 2>{0.0, 3.0}
});

VectorN<Real, 2> standard = scaled.coordinatesInStandard({4.0, 5.0});
```

```cpp
MatrixNM<Real, 2, 2> A{
  0.0, -1.0,
  1.0,  0.0
};
LinearMap<Real, 2, 2> rotate90(A);
VectorN<Real, 2> y = rotate90({1.0, 0.0});
```

```cpp
MatrixNM<Real, 2, 2> gram{
  2.0, 0.0,
  0.0, 1.0
};
InnerProductSpace<Real, 2> weightedPlane(gram);
Real length = weightedPlane.norm({3.0, 4.0});
```

---

## Migration Strategy

1. Keep storage APIs focused on representation and domain operations.
2. Implement semantic wrappers under `MML::VectorSpaces` using core numerical kernels.
3. Test semantic results directly against their mathematical contracts.
4. Add examples showing when to use raw coordinate storage and when to use the
   semantic layer.
5. Maintain modular-header and single-header coverage.

Existing code should not need edits unless it opts into the new vector-space
types.

---

## First Implementation Order

1. `Subspace<Real>` dynamic wrapper over SVD/fundamental subspace helpers.
2. `Basis<Scalar,N>` and fixed-size basis changes.
3. `LinearMap<Scalar,DomainN,CodomainN>` with kernel/image via `Subspace`.
4. `DualSpace` and `CovectorInSpace` bridge.
5. `InnerProductSpace` and orthogonal operations.
6. `AffineSpace` and affine maps.

This order is now implemented in the first public vector-space slice. Future
work can add dynamic-size siblings, stronger runtime identity, quotient spaces,
and deeper coordinate-system bridges without changing the storage model.