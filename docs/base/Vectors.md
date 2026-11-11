# Vectors — Overview

MML provides two complementary vector types for numerical computing:
- **`Vector<T>`**: Dynamic-length vector for runtime-sized data.
- **`VectorN<T, N>`**: Fixed-size vector optimized for small dimensions.

Use `Vector<T>` for variable datasets and streaming results; use `VectorN<T, N>` for compact math (coordinates, gradients, Jacobians) where compile-time size helps performance and safety.

## Highlights
- ✅ **Dual vector families**: Dynamic (`Vector<T>`) and fixed-size (`VectorN<T,N>`)
- ✅ **Stack allocation**: Fixed-size vectors enable compiler optimizations
- ✅ **Rich specialized types**: Cartesian, Polar, Spherical, Cylindrical, Lorentz
- ✅ **Seamless integration**: Works with matrices and geometry primitives
- ✅ **Clear APIs**: Safety (`at()`) vs. speed (`operator[]`)

## Documentation
- **[Vector.md](Vector.md)** - Dynamic-size vectors (`Vector<T>`)
- **[VectorN.md](VectorN.md)** - Fixed-size vectors (`VectorN<T, N>`)

### Runnable Examples

See individual documentation files for runnable demos:
- [Vector.md](Vector.md) - `Docs_Demo_Vector()` functions
- [VectorN.md](VectorN.md) - `Docs_Demo_VectorN()` functions

## Specialized Vector Types
Based on `VectorN<T,N>` with domain-specific semantics:
- `Vector2Cartesian`, `Vector3Cartesian` — Cartesian components
- `Vector2Polar` — 2D polar components
- `Vector3Spherical` — physics convention (r, theta, phi)
- `Vector3Cylindrical` — cylindrical components
- `Vector4Minkowski` — spacetime (t, x, y, z)

### Geometric arithmetic and scale

Typed polar, spherical, and cylindrical compound operators follow their binary
geometric operations, not inherited raw coordinate-wise arithmetic. Addition and
subtraction combine displacements. Scalar multiplication/division scale lengths
only (polar/spherical radius; cylindrical radius and height), leaving angles
unchanged. Existing negative-radius conventions are preserved. Spherical and
cylindrical division by zero remains an explicit error; polar scalar division
retains its unchecked IEEE contract. Polar addition/subtraction requires finite
radii and angles and uses a scaled midpoint frame to retain close-angle
differences without squaring enormous or tiny radii.

Cartesian direction queries normalize bounded components before dot/cross
comparisons. Parallel and perpendicular tolerances are dimensionless bounds on
the absolute sine and cosine of the angle, respectively. Finite, nonzero vectors
have directions at any representable magnitude, including subnormal magnitudes;
this is distinct from the unit-scale near-zero threshold of `Normalized()`.
Zero/nonfinite vectors or invalid (negative/nonfinite) tolerances return false
from direction predicates. Angles reject nonfinite inputs with `NumericInputError`
and zero vectors with `DivisionByZeroError`. Angles use `atan2` of cross magnitude
and dot product, rather than clamping a nonfinite cosine to apparent validity.
Spherical/cylindrical local unit-direction helpers likewise avoid magnitude
overflow/underflow and retain their existing zero-vector result convention.

Minkowski operations retain signature `(-,+,+,+)` and the existing signed norm:
positive for timelike vectors, negative for spacelike vectors, zero for
lightlike vectors. Exponent-scaled, compensated products allow norms and
classification even when the squared interval is not representable, including
small residuals after large lightlike cancellation. A directly returned scalar
product can still overflow/underflow when its mathematical result is outside
`Real`'s range. Classification predicates fail closed on nonfinite components;
norms/products reject them with `NumericInputError`.

```cpp
Vector2Polar polar(1, 0);
polar += Vector2Polar(1, Constants::PI / Real{2}); // radius sqrt(2), angle pi/4
const auto difference = Vector2Polar(1, 0) - Vector2Polar(1, Real{1e-10});
Vector3Cartesian huge(Real{1e200}, 0, 0), diagonal(Real{1e200}, Real{1e200}, 0);
const Real angle = huge.AngleToVector(diagonal); // pi/4, not zero
```

Runnable counterpart: `Docs_Demo_CoordinateScaleContracts()` in
`src/docs_demos/base/docs_demo_coordinate_scale_contracts.cpp`, selected with
`MML_DocsApp --arithmetic-scale-contracts`.

## Tips
- Prefer `VectorN` in tight inner loops with small dimensions to exploit stack allocation.
- Use `at()` when bounds safety matters; use `operator[]` where indices are guaranteed valid.
- Choose epsilon (`eps`) for approximate equality relative to magnitude and type.

## Quick Examples
```c++
using namespace MML;

Vector<Real> v({1.0, -2.0, 3.5});
VectorN<Real,3> p{1.0, 2.0, 3.0};
Vector3Spherical s{2.0, M_PI/4, M_PI/6};

auto len = p.NormL2();
auto sum = v + Vector<Real>({0.5, 0.5, 0.5});
```

## Related
- Matrices overview: [Matrices.md](Matrices.md)
- Geometry: [Geometry_2D_3D.md](Geometry_2D_3D.md)


