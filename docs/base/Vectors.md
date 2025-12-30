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
- `Vector4Lorentz` — spacetime (t, x, y, z)

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



