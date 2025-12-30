# MML Default Precision Constants

## Overview

This document catalogs MML's default precision constants and tolerances defined in `MMLPrecision.h`. These values are type-specialized for `float`, `double`, and `long double`.

## PrecisionValues Template

MML uses a template struct `PrecisionValues<T>` to provide type-appropriate tolerances:

```cpp
#include "MMLPrecision.h"

// Usage
double tol = MML::PrecisionValues<double>::DefaultTolerance;  // 1e-6
float tol_f = MML::PrecisionValues<float>::DefaultTolerance;  // 1e-6f
```

## Tolerance Categories

### Equality Comparison Tolerances

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `ComplexAreEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | Complex number equality |
| `MatrixIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | Matrix element-wise equality |
| `VectorIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | Vector element-wise equality |

### Point Equality Tolerances

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `Pnt2CartIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 2D Cartesian point equality |
| `Pnt2PolarIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 2D polar point equality |
| `Pnt3CartIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 3D Cartesian point equality |
| `Pnt3SphIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 3D spherical point equality |
| `Pnt3CylIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 3D cylindrical point equality |

### Vector Tolerances

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `Vec2CartIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 2D vector equality |
| `Vec3CartIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 3D Cartesian vector equality |
| `Vec3CartIsParallelTolerance` | 1e-6 | 1e-10 | 1e-15 | Vector parallelism check |
| `Vec3SphIsEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | 3D spherical vector equality |

### Line and Plane Tolerances

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `Line3DAreEqualTolerance` | 1e-6 | 1e-10 | 1e-15 | Line equality check |
| `Line3DIsPointOnLineTolerance` | 1e-6 | 1e-10 | 1e-15 | Point-on-line test |
| `Line3DIsPerpendicularTolerance` | 1e-6 | 1e-10 | 1e-15 | Perpendicularity test |
| `Line3DIsParallelTolerance` | 1e-6 | 1e-10 | 1e-15 | Parallelism test |
| `Line3DIntersectionTolerance` | 1e-6 | 1e-10 | 1e-15 | Line intersection |
| `Plane3DIsPointOnPlaneTolerance` | 1e-6 | 1e-10 | 1e-15 | Point-on-plane test |

### Triangle Tolerances

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `Triangle3DIsPointInsideTolerance` | 1e-6 | 1e-10 | 1e-15 | Point-in-triangle test |
| `Triangle3DIsRightTolerance` | 1e-6 | 1e-10 | 1e-12 | Right triangle check |
| `Triangle3DIsIsoscelesTolerance` | 1e-6 | 1e-10 | 1e-12 | Isosceles triangle check |
| `Triangle3DIsEquilateralTolerance` | 1e-6 | 1e-10 | 1e-12 | Equilateral triangle check |

### Matrix Property Tolerances

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `IsMatrixSymmetricTolerance` | 1e-6 | 1e-10 | 1e-15 | Symmetry check |
| `IsMatrixDiagonalTolerance` | 1e-6 | 1e-10 | 1e-15 | Diagonal matrix check |
| `IsMatrixUnitTolerance` | 1e-6 | 1e-10 | 1e-15 | Identity matrix check |
| `IsMatrixOrthogonalTolerance` | 1e-6 | 1e-10 | 1e-15 | Orthogonality check |

### Numerical Computation Thresholds

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `NumericalZeroThreshold` | 1e-12 | 1e-12 | 1e-12 | Near-zero value detection |
| `QuaternionZeroThreshold` | 1e-12 | 1e-12 | 1e-12 | Quaternion singularity |
| `SurfaceNormalThreshold` | 1e-12 | 1e-12 | 1e-12 | Surface normal magnitude |
| `DerivativeStepSize` | 1e-6 | 1e-6 | 1e-6 | Numerical derivative h |
| `DefaultTolerance` | 1e-6 | 1e-6 | 1e-6 | General purpose tolerance |
| `DefaultToleranceRelaxed` | 1e-8 | 1e-8 | 1e-8 | Relaxed geometric tests |

### Algorithm-Specific Thresholds

| Constant | float | double | long double | Purpose |
|----------|-------|--------|-------------|---------|
| `RankAlgEPS` | 1e-10 | 1e-12 | 1e-13 | Matrix rank computation |
| `EigenSolverZeroThreshold` | 1e-12 | 1e-15 | 1e-18 | Eigenvalue zero detection |
| `PolynomialCoeffZeroThreshold` | 1e-10 | 1e-12 | 1e-15 | Polynomial coefficient zero |
| `MatrixElementZeroThreshold` | 1e-12 | 1e-15 | 1e-18 | Matrix element near-zero |
| `DeterminantZeroThreshold` | 1e-12 | 1e-15 | 1e-18 | Determinant zero detection |
| `DivisionSafetyThreshold` | 1e-25 | 1e-30 | 1e-35 | Safe division threshold |
| `OrthogonalityTolerance` | 1e-8 | 1e-10 | 1e-12 | Orthogonal vector check |
| `LinearDependenceTolerance` | 1e-10 | 1e-12 | 1e-15 | Linear dependence check |

## PrintContext Defaults

From `MMLBase.h`, the `PrintContext` and `Defaults` namespace provide:

```cpp
namespace MML::Defaults {
    int VectorPrintWidth = 15;
    int VectorPrintPrecision = 10;
    int VectorNPrintWidth = 15;
    int VectorNPrintPrecision = 10;
    int MatrixPrintWidth = 15;
    int MatrixPrintPrecision = 10;
}
```

## Integration Parameters

```cpp
// From MMLBase.h (in IntegrationConfig)
Real rombergIntegrationEPS = PrecisionValues<Real>::DefaultTolerance;  // 1e-6
Real workIntegralPrecision = 1e-05;
Real lineIntegralPrecision = 1e-05;
```

## Precision Selection Guide

| Use Case | Recommended Precision |
|----------|----------------------|
| Graphics/visualization | 1e-4 to 1e-6 |
| General engineering | 1e-8 to 1e-10 |
| Scientific computing | 1e-12 to 1e-14 |
| Near machine precision | 1e-15 (double), 1e-6 (float) |

## Machine Epsilon Reference

| Type | Machine Epsilon | Decimal Digits |
|------|-----------------|----------------|
| `float` | ~1.19e-7 | ~7 |
| `double` | ~2.22e-16 | ~15-16 |
| `long double` | ~1.08e-19 (80-bit) | ~18-19 |

## Customization

To override defaults for specific operations:

```cpp
// Most functions accept optional tolerance parameter
bool equal = matrix1.IsEqual(matrix2, 1e-12);  // Override default

// Or modify the PrecisionValues at compile time via template specialization
```

## Related Documentation

- [METHODOLOGY.md](METHODOLOGY.md) - Testing methodology
- [LINEAR_ALGEBRA_ANALYSIS.md](LINEAR_ALGEBRA_ANALYSIS.md) - Solver precision
- [DERIVATION_ANALYSIS.md](DERIVATION_ANALYSIS.md) - Derivative precision

---

*Source file: `mml/MMLPrecision.h`*