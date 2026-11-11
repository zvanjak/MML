# Core Algebra And Optional Discrete Algebra

Core `mml/base/Algebra_base.h` now aggregates only `LieGroups.h`: SO2, SO3, SE2 and SE3.
These rotation/rigid-motion value types remain in the core single-header distribution.

Finite/discrete algebra has moved to the **private MML-Packages-Private** repository's
`include/mml_ext/base/Algebra/` tree. Link `MML_Ext` and switch the include to
`<mml_ext/base/Algebra_base.h>` (or the matching component under `mml_ext/base/Algebra/`).
The moved component namespaces, guards, APIs and numerical bodies are unchanged.
There are no old core forwarding headers or core dependencies on MML_Ext.

The extension owns traits, permutations, finite/cyclic/dihedral groups, group actions,
ModInt/modular rings, exact Polynomial/finite fields, representations, and their six
algorithm headers. Their docs, tests, and runnable cookbook Recipe 26 moved with them.
Use the receiving repository's `docs/mml_ext/DiscreteAlgebra.md` and
`MML_DocsApp --discrete-algebra` for those examples.

Numerical `mml/base/Polynom.h` is a different type and remains core. The generic exact
`Algebra::FieldMatrixDeterminant` and `Algebra::FieldMatrixInverse` also remain core,
now included explicitly through `mml/algorithms/MatrixAlg/FieldMatrixAlgorithms.h`.
They use exact zero tests and field division, not tolerance-based floating-point LU.
Basis and the extension share this one implementation. The base aggregate does not
include algorithm-layer routines implicitly.

Run `MML_DocsApp --algebra` for the retained geometry and shared matrix demo in
[docs_demo_algebra.cpp](../../src/docs_demos/core/docs_demo_algebra.cpp).

## Lie Groups SO(2) And SO(3)

The continuous rotation groups live under `mml/base/Algebra/LieGroups/` because they
are foundational geometry value types. They use the same active, right-handed,
column-vector convention as MML quaternions:

```text
compose(left, right) applies right first, then left
```

`SO2` stores a normalized angle in $[-\pi,\pi)$ and provides `matrix`, `compose`,
`inverse`, vector application, scalar `Exp`/`log`, geodesic distance, and shortest-path
interpolation.

`SO3` stores a canonical unit quaternion but exposes representation-independent group
operations:

- `FromMatrix` and `FromAxisAngle`
- rotation-vector `Exp` and principal `log`
- `Hat` and `Vee` between $\mathbb{R}^3$ and $\mathfrak{so}(3)$
- `matrix`, `compose`, `inverse`, and vector application
- geodesic distance and quaternion-backed shortest-path `Slerp`
- `Project` for re-orthonormalizing a near-rotation matrix

The exponential uses a series coefficient near zero to avoid division by a tiny angle.
The logarithm uses a canonical quaternion hemisphere and `atan2`, remaining stable for
small rotations and near $\pi$.

Canonicalization keeps the quaternion scalar component nonnegative. Only an
exactly zero scalar component uses a lexicographic vector sign tie-break;
small positive scalar values are not ties. Recanonicalizing through identity
composition therefore does not flip a near-half-turn quaternion's sign.
Principal logarithms have angle in $[0,\pi]$, up to floating-point roundoff,
and geodesic distances use the shortest rotation.

Matrix construction has deliberately separate semantics:

- `SO3::FromMatrix` is strict and rejects matrices that are not orthogonal with
    determinant $+1$ within the supplied tolerance.
- `SO3::Project` is the explicit repair path. It orthonormalizes independent columns and
    reconstructs a right-handed third axis before strict construction.

This distinction prevents malformed transforms from being silently accepted while still
supporting numerical integration and accumulated floating-point drift.

## Rigid Motions SE(2) And SE(3)

`SE2` and `SE3` pair an `SO2`/`SO3` rotation with a translation. They use the same
active, right-first composition convention:

$$
(R_1,t_1)(R_2,t_2)=(R_1R_2,R_1t_2+t_1).
$$

The API deliberately distinguishes affine points from free vectors:

```cpp
SE3 transform(rotation, translation);
auto movedPoint = transform.apply_point(point);   // R p + t
auto movedVector = transform.apply_vector(vector); // R v
```

`homogeneous_matrix()` returns the conventional column-vector form

$$
\begin{bmatrix}R&t\\0&1\end{bmatrix}.
$$

`FromHomogeneousMatrix` strictly validates the bottom row and rotation block. It does
not silently project malformed rotations; repair a rotation explicitly with
`SO3::Project` before constructing an `SE3` when numerical drift is expected.

Strict matrix construction requires **every** entry to be finite, including
translation and the homogeneous bottom row; nonfinite entries throw `DomainError`.
`SO3::IsRotationMatrix` returns false for nonfinite matrix entries, and
`SO3::FromMatrix` rejects them with `DomainError`.
The SO3 rotation-matrix checks and both homogeneous-matrix constructors require
finite, nonnegative tolerances; invalid tolerances throw `NumericInputError`.
Zero remains supported for exact comparisons. Positive tolerances retain their
existing bottom-row and rotation-block checks; malformed input is not silently
projected or repaired.

Frame-aware overloads map `Point<Real,N,FromFrame>` to
`Point<Real,N,ToFrame>` and `TangentVector<N,FromFrame>` to
`TangentVector<N,ToFrame>`. Translation affects only points, while tangent vectors use
the linear rotation action. These bridges remain in Base and do not depend on Core
coordinate-map infrastructure.

The optional $\mathfrak{se}(2)$/$\mathfrak{se}(3)$ exp/log maps are intentionally
deferred: the rigid-motion task keeps a compact, unambiguous value API and preserves
single-header budget for the remaining algebra layers.
