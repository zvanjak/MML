# Algebra Foundations

MML separates algebraic value types from computations over those types:

| Layer | Responsibility |
|-------|----------------|
| `mml/base/Algebra/` | Traits, metadata, and foundational algebraic value types. |
| `mml/base/Algebra/Algorithms/` | Exhaustive law checks and operations over algebraic structures. |
| `mml/algorithms/` | Future strategy-oriented algorithms for larger algebraic problems. |

Dependencies flow upward only: Base algebra must not include Core or Algorithms.
Continuous groups and rigid transformations may use existing Base vectors, matrices,
quaternions, points, and rank-1 tensors. Core provides adapters to Core facilities such
as typed coordinate maps.

## Structure Protocols

The C++17 API uses value types and free functions rather than an inheritance hierarchy.
The first supported protocols are:

```cpp
// Group
elements();
identity();
compose(left, right);
inverse(element);

// Ring
elements();
zero();
one();
add(left, right);
negate(element);
multiply(left, right);

// Field: ring protocol plus
inverse(nonzeroElement);
```

`compose(left, right)` returns the product `left * right`. When elements represent
transformations, `right` is applied first and `left` second. This matches matrix and
function composition.

`AlgebraTraits<Structure>` exposes `element_type`, `scalar_type`, static order,
dimension, and equality metadata. Structures may specialize it when the conservative
defaults are insufficient.

## Equality Policy

Exact finite algebra uses exact equality by default. Numerical structures must supply an
explicit equality predicate to law-checking functions. This prevents an implicit global
tolerance from changing the meaning of an algebraic law.

```cpp
auto approximatelyEqual = [tolerance](const Element& left, const Element& right) {
    return left.IsEqualTo(right, tolerance);
};

bool valid = Algebra::CheckGroupLaws(group, approximatelyEqual);
```

## Law Checks

Include `mml/base/Algebra_base.h` for:

- `CheckClosure`
- `CheckAssociativity`
- `CheckIdentity`
- `CheckInverses`
- `CheckGroupLaws`
- `CheckRingLaws`
- `CheckFieldLaws`

The checks exhaustively enumerate finite structures and have cubic complexity in the
number of elements. They are intended for small structures, examples, and tests. Large
groups will require sampled or structure-specific algorithms in a higher layer.

Tests use tiny hand-defined structures with known laws and deliberate failures. This
keeps the foundation independent from later concrete implementations such as cyclic
groups and modular integers.

## Permutations And Finite Groups

Include `mml/base/Algebra_base.h` for the foundational permutation types:

```cpp
using namespace MML::Algebra;

Permutation<4> rotation = Permutation<4>::from_cycles({{0, 1, 2, 3}});
DynamicPermutation runtimeRotation =
    DynamicPermutation::from_cycles(4, {{0, 1, 2, 3}});
```

Both types use zero-based images, where `images[i]` is the destination `p(i)`.
Composition follows the library-wide convention:

```text
compose(left, right)(i) = left(right(i))
```

Applying a permutation to a value container moves the value at source position `i` to
destination position `p(i)`. Thus the cycle `(0 1 2 3)` maps `[A,B,C,D]` to
`[D,A,B,C]`.

Available operations include identity, composition, inverse, cycle construction and
decomposition, parity/sign, minimal transposition count, permutation order, index
application, and safe application to equal-sized containers. Constructors reject
duplicate and out-of-range images.

`FiniteGroup<Element>` stores a finite element set, identity, composition, and inverse
operations. It validates basic representation invariants such as a present identity and
unique elements, while algebraic-law validation remains in Core:

```cpp
FiniteGroup<int> c3(
    {0, 1, 2}, 0,
    [](const int& a, const int& b) { return (a + b) % 3; },
    [](const int& a) { return (3 - a) % 3; });

bool valid = CheckGroupLaws(c3);
```

## Cyclic And Dihedral Groups

`CyclicGroup(n)` represents $C_n$ with canonical elements `r^k`. Exponents are
normalized modulo $n$, including negative exponents:

```cpp
CyclicGroup c5(5);
auto r = c5.generator();
auto r4 = c5.rotation(-1);
auto identity = c5.compose(r, c5.inverse(r));
```

`DihedralGroup(n)` represents the order-$2n$ symmetry group $D_n$ of a regular
$n$-gon. Elements use the canonical form $r^k s^f$, stored as `rotation` and
`reflected`. The API exposes rotations, indexed reflections, and the distinguished
reflection $s$:

```cpp
DihedralGroup d4(4);
auto r = d4.generator();
auto s = d4.reflection();
auto relation = d4.compose(d4.compose(s, r), s); // r^-1
```

Both groups expose their natural action on polygon vertex indices as a
`DynamicPermutation`. Rotation `r^k` maps vertex $i$ to $i+k$. The distinguished
dihedral reflection maps $i$ to $-i$, with all indices reduced modulo $n$.

Core finite-group algorithms work with any structure following the group protocol:

- `element_order(group, element)`
- `CayleyTable(group)`, returning stable indices into `group.elements()`
- `GeneratedSubgroup(group, generators)`
- `ConjugacyClasses(group)`
- `MakeCayleyGraphData(group, generators)`

Cayley graph data contains an element list plus directed edges identified by source
index, destination index, and generator index. It deliberately contains no rendering
or file-format dependency, so Tools and external visualizers can serialize it later.

## Group Actions And Burnside Counting

`GroupAction<GroupElement, Object>` is the Base value adapter for a left action. Its
convention is `apply(g, x) = g.x`, with the laws

$$
e.x=x, \qquad (gh).x=g.(h.x).
$$

Use `MakeGroupAction` to wrap a lambda without introducing an inheritance hierarchy:

```cpp
DihedralGroup d5(5);
auto vertices = MakeGroupAction<DihedralElement, int>(
    [&d5](const DihedralElement& g, const int& vertex) {
        return d5.permutation(g).apply(vertex);
    });

auto orbit = Orbit(d5, vertices, 0);
auto stabilizer = Stabilizer(d5, vertices, 0);
```

Core provides:

- `Orbit` and `OrbitPartition`
- `Stabilizer`
- `IsInvariant` for an object under every group element
- `IsPredicateInvariant` for predicates preserved across a finite object set
- `FixedPoints` and `FixedPointCount`
- `BurnsideCount`

For a finite $G$-set $X$, Burnside's lemma computes

$$
|X/G|=\frac{1}{|G|}\sum_{g\in G}|\operatorname{Fix}(g)|.
$$

`BurnsideCount` verifies that every transformed object remains in the supplied finite
set and that the fixed-point average is integral. This catches incomplete coloring or
motif sets instead of returning a silently truncated count. Object equality is an
explicit policy parameter, allowing exact discrete actions and tolerance-aware
numerical actions to use the same algorithms.

## Modular Arithmetic And Prime Fields

`ModInt<Modulus>` is the exact scalar type for the residue ring
$\mathbb{Z}/n\mathbb{Z}$. Construction normalizes signed integers into the canonical
range $[0,n)$, and addition, subtraction, and multiplication use widened intermediate
arithmetic:

```cpp
using Z12 = ModInt<12>;
Z12 value = Z12(10) + Z12(5); // 3
```

Division is defined only by units. Calling `inverse()` or division on a non-unit throws
`DomainError`; inversion of zero throws `DivisionByZeroError`. This preserves the
difference between a general modular ring and a field with no zero divisors.

`IsPrime<P>` provides compile-time primality for small integer moduli.
`PrimeFieldElement<P>` rejects composite $P$ at compile time and aliases the exact
`ModInt<P>` scalar:

```cpp
static_assert(IsPrime<7>::value);
using F7 = PrimeFieldElement<7>;
F7 quotient = F7(3) / F7(2); // 5
```

`ModularRing<n>` and `PrimeField<p>` expose finite `elements`, `zero`, `one`, and
arithmetic operations for `CheckRingLaws` and `CheckFieldLaws`.

Core number-theory helpers include `PrimitiveRoot<p>`, `LegendreSymbol<p>`, and
`IsQuadraticResidue<p>`. `FieldMatrixDeterminant` and `FieldMatrixInverse` perform exact
Gaussian elimination over fixed-size `MatrixNM` field matrices. They intentionally do
not use `MatrixNM::GetInverse()`, whose pivot policy is designed for floating-point
matrices and numerical tolerances.

## Exact Polynomials And Extension Fields

`Algebra::Polynomial<Coeff>` is deliberately separate from MML's numerical `Polynom`.
It stores canonical exact coefficients in ascending-power order and removes trailing
zero coefficients. Its scope is finite algebra: addition, subtraction, multiplication,
evaluation, monic normalization, and division with remainder over coefficient fields.

Core provides `PolynomialGcd`, `PolynomialExtendedGcd`, `PolynomialPowMod`, and
`IsIrreducible`. The extended GCD result satisfies the exact Bézout identity

$$
s(x)a(x)+t(x)b(x)=\gcd(a,b).
$$

`FiniteFieldElement<P,N,ModulusPolynomial>` represents
$GF(P)[x]/(f(x))$ in the polynomial basis $1,x,\ldots,x^{N-1}$. The provider must expose
a `static modulus()` returning a documented monic, irreducible degree-$N$ polynomial:

```cpp
using F2 = PrimeFieldElement<2>;

struct GF4Modulus {
    static Polynomial<F2> modulus() { return {1, 1, 1}; } // x^2+x+1
};

using GF4 = FiniteFieldElement<2, 2, GF4Modulus>;
GF4 alpha{0, 1};
GF4 relation = alpha * alpha; // alpha + 1
```

The value type validates degree and monicity when the provider is first used. Use
`IsIrreducible` in tests or setup validation to verify user-supplied modulus
polynomials. Addition is coefficient-wise, multiplication is reduced modulo $f$, and
nonzero inversion uses $a^{-1}=a^{P^N-2}$.

`ExtensionField<P,N,Provider>` enumerates small fields in base-$P$ coefficient order
and exposes the finite-field protocol used by `CheckFieldLaws`. It is intended for
small educational and validation fields, not for enumerating cryptographic-size fields.

## Representations And Invariant Projection

`Representation<GroupElement,Scalar,N>` is the Base value adapter for a matrix-valued
map $\rho:G\to GL(N,Scalar)$. It stores an evaluator, exposes `matrix(g)`, and applies
the represented group element to `VectorN<Scalar,N>`:

```cpp
auto rho = MakeRepresentation<CyclicElement, Real, 2>(
    [](const CyclicElement& element) {
        // Return the matrix for r^element.exponent.
        return MatrixNM<Real, 2, 2>::Identity();
    });
```

Core provides:

- `VerifyRepresentation`, checking $\rho(gh)=\rho(g)\rho(h)$ and $\rho(e)=I$
- `Character`, returning $\chi(g)=\operatorname{tr}(\rho(g))$
- `IsCharacterConstantOnConjugacyClasses`
- `InvariantProjection`
- `SymmetrizeVector`
- generic `GroupAverage` for other additive and scalable objects

For finite groups over a scalar field whose characteristic does not divide $|G|$,
the invariant projection is

$$
P=\frac{1}{|G|}\sum_{g\in G}\rho(g).
$$

It satisfies $P^2=P$, and every vector in its image is fixed by the represented group.
Numerical representations should pass explicit matrix and scalar equality predicates to
homomorphism and character checks. Exact representations can use the default exact
comparisons. `GroupAverage` accepts explicit addition and scaling operations, allowing
the same mechanism to symmetrize vectors, tensors, and sampled functions without moving
those domain types into the Algebra layer.

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

Frame-aware overloads map `Point<Real,N,FromFrame>` to
`Point<Real,N,ToFrame>` and `TangentVector<N,FromFrame>` to
`TangentVector<N,ToFrame>`. Translation affects only points, while tangent vectors use
the linear rotation action. These bridges remain in Base and do not depend on Core
coordinate-map infrastructure.

The optional $\mathfrak{se}(2)$/$\mathfrak{se}(3)$ exp/log maps are intentionally
deferred: the rigid-motion task keeps a compact, unambiguous value API and preserves
single-header budget for the remaining algebra layers.

## Lattices And Crystallographic Symmetry

`Lattice<2>` and `Lattice<3>` store basis vectors as matrix columns. They convert
between lattice and Cartesian coordinates, compute unit-cell volume, return nearest
lattice points by coordinate rounding, and construct the reciprocal basis
$B=2\pi A^{-T}$, so $a_i\cdot b_j=2\pi\delta_{ij}$.

`AffineSymmetry<N>` represents $x\mapsto Rx+t$ with composition, inverse, point action,
and translation-free vector action. `PreservesLattice` checks that both
$A^{-1}RA$ and $A^{-1}t$ have integer components within tolerance.

`ExpandMotif` applies a finite symmetry list, wraps positions into the half-open unit
cell $[0,1)^N$ in lattice coordinates, and removes equivalent positions. This is the
initial machinery layer only; it intentionally does not embed a database of named space
groups or perform lattice reduction.