# MML Equality and Tolerance Policy

This document describes the unified policy for floating-point equality comparisons across the MinimalMathLibrary (MML).

## Core Principles

1. **`operator==` is EXACT** - Always performs exact bitwise comparison
2. **`IsEqualTo()` uses tolerance** - For approximate floating-point comparison
3. **Use appropriate tolerance type** - Absolute, relative, or combined based on context
4. **Angles require special handling** - Wrap-around at ±π must be considered

---

## Comparison Functions

### Low-Level Functions (in `MML` namespace)

| Function | Purpose | When to Use |
|----------|---------|-------------|
| `isWithinAbsPrec(a, b, eps)` | Absolute tolerance | Values expected near zero |
| `isWithinRelPrec(a, b, eps)` | Relative tolerance | Values of similar magnitude away from zero |
| `isNearlyEqual(a, b, absEps, relEps)` | Combined abs+rel | **Recommended for general use** |
| `isNearlyEqual(a, b, eps)` | Combined (same eps) | Convenience overload |
| `AnglesAreEqual(a, b, eps)` | Angle wrap-aware | Comparing angles in radians |
| `normalizeAngle(rad)` | Normalize to [-π, π) | Pre-processing for angle operations |

### Type Methods

All MML types use `IsEqualTo()` for tolerance-based comparison:

```cpp
// Vectors
Vector3Cartesian v1, v2;
v1.IsEqualTo(v2);                              // Uses Defaults::Vec3CartIsEqualTolerance
v1.IsEqualTo(v2, 1e-12);                       // Custom tolerance

// Points  
Point3Cartesian p1, p2;
p1.IsEqualTo(p2);                              // Uses Defaults::Pnt3CartIsEqualTolerance

// Matrices
Matrix<Real> m1, m2;
m1.IsEqualTo(m2);                              // Uses Defaults::MatrixIsEqualTolerance
Matrix<Real>::AreEqual(m1, m2, 1e-15);         // Static helper with custom tolerance
```

---

## Default Tolerances

All tolerances are type-specialized via `PrecisionValues<T>` template:

| Type | float | double | long double |
|------|-------|--------|-------------|
| Base tolerance | 1e-6 | 1e-10 | 1e-15 |

Access via `Defaults::` namespace:

```cpp
#include "MMLBase.h"

// Equality tolerances
Defaults::MatrixIsEqualTolerance       // Matrix element comparison
Defaults::VectorIsEqualTolerance       // Vector element comparison
Defaults::Pnt3CartIsEqualTolerance     // 3D point distance comparison

// Geometry tolerances
Defaults::Line3DIsParallelTolerance    // Line parallelism
Defaults::Plane3DIsPointOnPlaneTolerance

// New tolerances (v1.1)
Defaults::AngleIsEqualTolerance        // Angle wrap-aware comparison
Defaults::ShapePropertyTolerance       // IsRight, IsSquare, IsEquilateral, etc.
Defaults::OrthogonalityTolerance       // Perpendicularity checks
Defaults::DefaultTolerance             // General-purpose tolerance
```

---

## Angle Comparison

Angles wrap around at ±π, requiring special handling:

```cpp
Real angle1 = Constants::PI - 0.001;
Real angle2 = -Constants::PI + 0.001;

// WRONG - will return false (difference appears to be ~2π)
std::abs(angle1 - angle2) < 1e-3;  // false!

// CORRECT - handles wrap-around
AnglesAreEqual(angle1, angle2, 1e-3);  // true

// Vector3Spherical uses angle-aware comparison for phi automatically
Vector3Spherical v1(1.0, 0.5, Constants::PI - 0.001);
Vector3Spherical v2(1.0, 0.5, -Constants::PI + 0.001);
v1.IsEqualTo(v2, 1e-3);  // true - phi comparison is angle-aware
```

---

## Shape Property Tolerances

Geometric shape classification methods use `Defaults::ShapePropertyTolerance`:

```cpp
Triangle t(3, 4, 5);
t.IsRight();                    // Uses Defaults::ShapePropertyTolerance
t.IsRight(1e-12);               // Custom tolerance

Rectangle r(5.0, 5.00000001);
r.IsSquare();                   // true with default tolerance

Ellipse e(1.0, 0.9999999);
e.IsCircle();                   // true with default tolerance
```

---

## Integer/Discrete Comparison

For discrete values, use exact comparison:

```cpp
// Polynomial coefficients (when integer/exact)
Polynom<int> p1, p2;
p1.IsEqual(p2);        // Exact comparison (no tolerance)
p1 == p2;              // Also exact

// Floating-point polynomial coefficients
Polynom<double> fp1, fp2;
fp1.IsEqualTo(fp2);    // Tolerance-based comparison
fp1.IsEqualTo(fp2, 1e-12);  // Custom tolerance
```

---

## Best Practices

### DO ✅

1. Use `IsEqualTo()` for floating-point comparisons in tests
2. Use `Defaults::*` constants for consistency
3. Use `AnglesAreEqual()` when comparing angles
4. Prefer `isNearlyEqual()` with combined tolerance for general cases

### DON'T ❌

1. Never use `==` for floating-point equality checks (except exact discrete values)
2. Avoid hardcoded tolerance values like `1e-10` - use `Defaults::*`
3. Don't forget angle wrap-around when comparing spherical/polar coordinates

---

## Migration Guide

If updating code from older MML versions:

| Old Pattern | New Pattern |
|-------------|-------------|
| `a.IsEqual(b)` | `a.IsEqualTo(b)` |
| `eps = 1e-10` | `eps = Defaults::ShapePropertyTolerance` |
| `std::abs(angle1 - angle2) < eps` | `AnglesAreEqual(angle1, angle2, eps)` |

---

## See Also

- [MMLPrecision.h](../../mml/MMLPrecision.h) - Type-specialized tolerance definitions
- [MMLBase.h](../../mml/MMLBase.h) - Core comparison functions and `Defaults::` namespace
- [DefaultPrecisions.md](../testing_precision/DefaultPrecisions.md) - Complete tolerance catalog
