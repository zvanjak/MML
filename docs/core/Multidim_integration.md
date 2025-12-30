# Multidimensional Integration

**File**: `mml/core/Integration/Integration2D.h`, `mml/core/Integration/Integration3D.h`

Numerical integration in 2D and 3D using nested 1D quadrature methods.

## Table of Contents
- [Overview](#overview)
- [2D Integration](#2d-integration)
- [3D Integration](#3d-integration)
- [Examples](#examples)
- [See Also](#see-also)
- [Runnable Examples](#runnable-examples)

---

## Overview

MML provides deterministic integration for low dimensions using tensor-product rules:

| Dimension | Function | Domain | Methods |
|-----------|----------|--------|---------|
| 2D | `Integrate2D` | General (variable y limits) | TRAP, SIMPSON, GAUSS10 |
| 3D | `Integrate3D` | General (variable y,z limits) | TRAP, SIMPSON, GAUSS10 |

For higher dimensions (≥ 4) or complex domains, see [Monte Carlo integration](../algorithms/MonteCarlo_integration.md).

---

## 2D Integration

### API
```cpp
IntegrationResult Integrate2D(
    const IScalarFunction<2>& func,     // f(x,y)
    IntegrationMethod method,           // TRAP, SIMPSON, or GAUSS10
    Real x1, Real x2,                   // x range: [x1, x2]
    Real (*y1)(Real),                   // y_min(x) - lower y boundary as function of x
    Real (*y2)(Real)                    // y_max(x) - upper y boundary as function of x
);
```

### Usage
```cpp
// Integrate f(x,y) = x*y over unit square [0,1] × [0,1]
ScalarFunction<2> f([](const VectorN<Real,2>& v) { 
    return v[0] * v[1]; 
});

// Constant y limits for rectangular domain
auto result = Integrate2D(f, SIMPSON, 
    0.0, 1.0,                          // x: [0, 1]
    [](Real x) { return 0.0; },        // y_min = 0
    [](Real x) { return 1.0; }         // y_max = 1
);
// Expected: ∫∫ xy dxdy = 0.25

// Variable y limits for circular domain (radius 2)
auto area = Integrate2D(
    ScalarFunction<2>([](const VectorN<Real,2>&) { return 1.0; }),
    GAUSS10,
    -2.0, 2.0,                         // x: [-2, 2]
    [](Real x) { return -sqrt(4 - x*x); },   // y_min = -√(4-x²)
    [](Real x) { return sqrt(4 - x*x); }     // y_max = +√(4-x²)
);
// Expected: π * r² = 4π ≈ 12.566
```

---

## 3D Integration

### API
```cpp
IntegrationResult Integrate3D(
    const IScalarFunction<3>& func,     // f(x,y,z)
    Real x1, Real x2,                   // x range
    Real (*y1)(Real),                   // y_min(x)
    Real (*y2)(Real),                   // y_max(x)
    Real (*z1)(Real, Real),             // z_min(x, y)
    Real (*z2)(Real, Real)              // z_max(x, y)
);
```

### Usage
```cpp
// Volume of unit sphere
ScalarFunction<3> unit([](const VectorN<Real,3>&) { return 1.0; });

auto vol = Integrate3D(unit,
    -1.0, 1.0,                                        // x: [-1, 1]
    [](Real x) { return -sqrt(1 - x*x); },            // y_min
    [](Real x) { return sqrt(1 - x*x); },             // y_max
    [](Real x, Real y) { return -sqrt(1 - x*x - y*y); },  // z_min
    [](Real x, Real y) { return sqrt(1 - x*x - y*y); }    // z_max
);
// Expected: 4/3 * π ≈ 4.189
```

---

## Examples

### Example 1: Area of Ellipse
```cpp
// Ellipse with semi-axes a=3, b=2
// Area = π * a * b = 6π
Real a = 3.0, b = 2.0;

auto area = Integrate2D(
    ScalarFunction<2>([](const VectorN<Real,2>&) { return 1.0; }),
    SIMPSON,
    -a, a,
    [b, a](Real x) { return -b * sqrt(1 - x*x/(a*a)); },
    [b, a](Real x) { return b * sqrt(1 - x*x/(a*a)); }
);
```

### Example 2: Mass of Density Distribution
```cpp
// Density: ρ(x,y) = 1 + x² + y²
// Domain: unit disk
ScalarFunction<2> density([](const VectorN<Real,2>& v) { 
    return 1.0 + v[0]*v[0] + v[1]*v[1]; 
});

auto mass = Integrate2D(density, SIMPSON,
    -1.0, 1.0,
    [](Real x) { return -sqrt(1 - x*x); },
    [](Real x) { return sqrt(1 - x*x); }
);
```

---

## See Also
- [Integration.md](Integration.md) - 1D integration methods
- [Monte Carlo integration](../algorithms/MonteCarlo_integration.md) - High-dimensional integration
- [Path_integration.md](../algorithms/Path_integration.md) - Line integrals
- [Surface_integration.md](../algorithms/Surface_integration.md) - Surface integrals

---

## Runnable Examples

| Example | Source File | Description |
|---------|------------|-------------|
| Multidim Integration Demo | [docs_demo_integration_multidim.cpp](../../src/docs_demos/docs_demo_integration_multidim.cpp) | 2D and 3D integration |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/docs_demos/Release/MML_DocsApp
```