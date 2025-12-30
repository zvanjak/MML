# Parametric Curves & Surfaces Test Bed Documentation

## Overview

The Parametric Curves & Surfaces Test Beds provide comprehensive test infrastructure for validating MML's differential geometry algorithms. These test beds include 3D parametric curves with known curvature and torsion, and parametric surfaces with known Gaussian and mean curvatures.

**Key Features:**
- 13 test curves with analytical derivatives, curvature, and torsion
- 18 parametric surfaces with geometric properties (curvature, area, minimal/developable)
- 4 explicit surfaces (z = f(x,y))
- Support for arc-length parameterized curves
- Rich metadata for each test case

---

## Part 1: Parametric Curves Test Bed

### Core Files

| File | Purpose |
|------|---------|
| [parametric_curves_test_bed.h](../../test_data/parametric_curves_test_bed.h) | `TestCartesianCurve3D` struct, `ParametricCurvesTestBed` class |
| [parametric_curves_defs.h](../../test_data/parametric_curves_defs.h) | Curve function definitions (position, derivatives, curvature, torsion) |

### Data Structure

```cpp
namespace MML::TestBeds {

/// Test data for a 3D Cartesian parametric curve
struct TestCartesianCurve3D {
    std::string _curveName;           ///< Curve identifier
    std::string _curveExpr;           ///< Mathematical expression
    Real _start, _end;                ///< Parameter domain [t_start, t_end]
    
    CurveCartesian3D _curve;          ///< Position function r(t)
    CurveCartesian3D _curveDerived;   ///< First derivative r'(t)
    CurveCartesian3D _curveDerSecond; ///< Second derivative r''(t)
    RealFunction _curvatureFunc;      ///< Analytical curvature κ(t)
    RealFunction _torsionFunc;        ///< Analytical torsion τ(t)
};

}
```

### Test Curves Catalog (13 total)

#### Basic Curves

| # | Name | Expression r(t) | Domain | κ (Curvature) | τ (Torsion) |
|---|------|-----------------|--------|---------------|-------------|
| 1 | Helix | {cos(t), sin(t), t} | [0, 2π] | 0.5 (const) | 0.5 (const) |
| 2 | Helix2 | {2cos(t), 2sin(t), 0.5t} | [0, 2π] | 2/4.25 | 0.5/4.25 |
| 3 | Schaums1 | {3t-t³, 3t², 3t+t³} | [0, 2π] | 1/(3(1+t²)²) | 1/(3(1+t²)²) |
| 4 | Schaums2 | {t-sin(t), 1-cos(t), t} | [0, 2π] | varies | varies |
| 5 | Circle3DXY | {cos(t), sin(t), 0} | [0, 2π] | 1.0 (const) | 0.0 (planar) |

#### Polynomial Curves

| # | Name | Expression r(t) | Domain | Description |
|---|------|-----------------|--------|-------------|
| 6 | TwistedCubic | {t, t², t³} | [-2, 2] | Classic algebraic geometry curve |

#### Trigonometric Curves

| # | Name | Expression r(t) | Domain | Description |
|---|------|-----------------|--------|-------------|
| 7 | Viviani | {1+cos(t), sin(t), 2sin(t/2)} | [0, 2π] | Sphere-cylinder intersection |
| 8 | Lissajous_2_3_5 | {sin(2t), sin(3t), sin(5t)} | [0, 2π] | 3D Lissajous figure |

#### Knots

| # | Name | Expression r(t) | Domain | Description |
|---|------|-----------------|--------|-------------|
| 9 | TorusKnot_3_2 | (2+cos(2t))·{cos(3t), sin(3t)} + sin(2t)·ẑ | [0, 2π] | Winds 3× around tube, 2× around axis |
| 10 | TrefoilKnot | {sin(t)+2sin(2t), cos(t)-2cos(2t), -sin(3t)} | [0, 2π] | Simplest non-trivial knot |
| 11 | FigureEightKnot | (2+cos(2t))·{cos(t), sin(t)} + sin(3t)·ẑ | [0, 2π] | Second simplest prime knot |

#### Special Curves

| # | Name | Expression r(t) | Domain | Description |
|---|------|-----------------|--------|-------------|
| 12 | ConicalHelix | {t·cos(t), t·sin(t), t} | [0.1, 4π] | Spiral on cone |
| 13 | SphericalSpiral | {cos(t)/cosh(0.3t), sin(t)/cosh(0.3t), tanh(0.3t)} | [-10, 10] | Loxodrome on sphere |

### ParametricCurvesTestBed API

```cpp
class ParametricCurvesTestBed {
public:
    /// Total number of test curves
    static int getNumTestCurves();              // Returns 13
    static int getNumTestCurvesArcLenParam();   // Returns 1
    
    /// Access by index
    static const TestCartesianCurve3D& getTestCurve(int i);
    static const TestCartesianCurve3D& getTestCurveArcLenParam(int i);
    
    /// Access by name (throws if not found)
    static const TestCartesianCurve3D& getTestCurve(const std::string& curveName);
    static const TestCartesianCurve3D& getTestCurveArcLenParam(const std::string& curveName);
};
```

---

## Part 2: Parametric Surfaces Test Bed

### Core Files

| File | Purpose |
|------|---------|
| [parametric_surfaces_test_bed.h](../../test_data/parametric_surfaces_test_bed.h) | `TestParametricSurfaceRect3` struct, test bed classes |
| [parametric_surfaces_defs.h](../../test_data/parametric_surfaces_defs.h) | Surface function definitions |

### Data Structure

```cpp
namespace MML::TestBeds {

/// Test data for a 3D parametric surface over rectangular domain
struct TestParametricSurfaceRect3 {
    std::string _surfaceName;         ///< Surface identifier
    std::string _surfaceExpr;         ///< Mathematical expression
    Real _u1, _u2;                     ///< u parameter range
    Real _v1, _v2;                     ///< v parameter range
    
    ParametricSurfaceRect<3> _surface; ///< The surface r(u,v)
    
    // Geometric properties
    Real _expectedArea;               ///< Analytical area (NaN if unknown)
    Real _gaussianCurvature;          ///< Constant K (NaN if varies)
    bool _isMinimal;                  ///< Mean curvature H = 0 everywhere
    bool _isDevelopable;              ///< Gaussian curvature K = 0 everywhere
};

}
```

### Test Surfaces Catalog (18 total)

#### Classic Surfaces (1-12)

| # | Name | Expression r(u,v) | K | Properties |
|---|------|-------------------|---|------------|
| 1 | Unit Sphere | (cos(u)sin(v), sin(u)sin(v), cos(v)) | 1.0 | K > 0 constant, Area = 4π |
| 2 | Torus | ((2+cos(v))cos(u), (2+cos(v))sin(u), sin(v)) | varies | K > 0 outside, K < 0 inside |
| 3 | Cylinder | (cos(u), sin(u), v) | 0.0 | Developable, Area = 4π |
| 4 | Cone | (v·cos(u), v·sin(u), v) | 0.0 | Developable |
| 5 | Helicoid | (v·cos(u), v·sin(u), u) | varies | **Minimal** (H=0), ruled |
| 6 | Catenoid | (cosh(v)cos(u), cosh(v)sin(u), v) | varies | **Minimal** (H=0) |
| 7 | Enneper Surface | (u-u³/3+uv², v-v³/3+vu², u²-v²) | varies | **Minimal**, self-intersecting |
| 8 | Möbius Strip | ((1+v/2·cos(u/2))cos(u), ..., v/2·sin(u/2)) | varies | Non-orientable |
| 9 | Hyperboloid One Sheet | (cosh(v)cos(u), cosh(v)sin(u), sinh(v)) | < 0 | Doubly ruled |
| 10 | Paraboloid | (u, v, u²+v²) | > 0 | Bowl shape |
| 11 | Hyperbolic Paraboloid | (u, v, u²-v²) | < 0 | Saddle, doubly ruled |
| 12 | Klein Bottle | Figure-8 immersion | varies | Non-orientable, self-intersecting |

#### Additional Surfaces (13-18)

| # | Name | Description | Special Properties |
|---|------|-------------|-------------------|
| 13 | Scherk Surface | z = ln(cos(v)/cos(u)) | **Minimal**, doubly periodic |
| 14 | Dini Surface | Twisted pseudosphere | K = -1 (constant negative) |
| 15 | Boy Surface | Bryant-Kusner RP² | Non-orientable |
| 16 | Roman Surface | Steiner's RP² | Tetrahedral symmetry |
| 17 | Cross Cap | Simplest RP² immersion | Non-orientable |
| 18 | Henneberg Surface | Henneberg (1875) | **Minimal** |

### Explicit Surfaces (z = f(x,y))

| # | Name | Expression | Domain |
|---|------|------------|--------|
| 1 | Upper Unit Hemisphere | z = √(1-x²-y²) | [-0.99, 0.99]² |
| 2 | Monkey Saddle | z = x(x²-3y²) | [-1, 1]² |
| 3 | Egg Carton | z = sin(x)·sin(y) | [-π, π]² |
| 4 | Gaussian Bump | z = exp(-(x²+y²)) | [-2, 2]² |

### ParametricSurfaceTestBed API

```cpp
class ParametricSurfaceTestBed {
public:
    static int getNumTestSurfaces3();  // Returns 18
    
    static const TestParametricSurfaceRect3& getTestSurface3(int i);
    static const TestParametricSurfaceRect3& getTestSurface3(const std::string& surfaceName);
};

class ExplicitSurfaceTestBed {
public:
    static int getNumExplicitSurfaces3();  // Returns 4
    
    static const TestParametricSurfaceRect3& getExplicitSurface3(int i);
    static const TestParametricSurfaceRect3& getExplicitSurface3(const std::string& surfaceName);
};
```

---

## MML API Reference

### Curve Methods (ICurveCartesian3D)

```cpp
class ICurveCartesian3D : public IParametricCurve<3> {
public:
    // Derivatives
    Vec3Cart getTangent(Real t) const;        ///< r'(t) - first derivative
    Vec3Cart getTangentUnit(Real t) const;    ///< T = r'/|r'| - unit tangent
    Vec3Cart getNormal(Real t) const;         ///< r''(t) - second derivative
    Vec3Cart getNormalUnit(Real t) const;     ///< N - principal normal
    Vec3Cart getBinormal(Real t) const;       ///< B = T × N - binormal
    
    // Curvature and Torsion
    Vec3Cart getCurvatureVector(Real t) const;
    Real getCurvature(Real t) const;          ///< κ = |r' × r''| / |r'|³
    Real getTorsion(Real t) const;            ///< τ = (r' × r'') · r''' / |r' × r''|²
    
    // Frenet Frame Planes
    Plane3D getOsculationPlane(Real t) const;   ///< Plane containing T and N
    Plane3D getNormalPlane(Real t) const;       ///< Plane perpendicular to T
    Plane3D getRectifyingPlane(Real t) const;   ///< Plane containing T and B
    
    // Moving trihedron (Frenet frame)
    void getMovingTrihedron(Real t, Vector3Cartesian& T, 
                            Vector3Cartesian& N, Vector3Cartesian& B);
    
    // Arc-length check
    bool isArcLengthParametrized(Real t1, Real t2, int numPoints = 100) const;
};
```

### Surface Methods (ParametricSurfaceRect<3>)

```cpp
class ParametricSurfaceRect<3> {
public:
    // Evaluation
    VectorN<Real,3> operator()(Real u, Real v) const;  ///< r(u,v)
    
    // Normal vector
    VectorN<Real,3> Normal(Real u, Real v) const;      ///< n = (∂r/∂u × ∂r/∂v) / |...|
    
    // Curvatures
    Real GaussianCurvature(Real u, Real v);   ///< K = (LN - M²) / (EG - F²)
    Real MeanCurvature(Real u, Real v);       ///< H = (EN + GL - 2FM) / (2(EG - F²))
    void PrincipalCurvatures(Real u, Real v, Real& k1, Real& k2);
    void PrincipalDirections(Real u, Real v, VectorN<Real,3>& d1, VectorN<Real,3>& d2);
    
    // Surface classification
    bool isFlat(Real u, Real v, Real eps = 1e-10) const;       ///< |K| < eps
    bool isDevelopable(Real u, Real v, Real eps = 1e-10) const; ///< K ≈ 0, not flat
    bool isHyperbolic(Real u, Real v, Real eps = 1e-10) const; ///< K < -eps
    bool isRegular(Real u, Real v, Real eps = 1e-10) const;    ///< Non-degenerate
};
```

---

## Usage Examples

### Testing Curve Curvature

```cpp
#include "test_data/parametric_curves_test_bed.h"
#include "core/Curves.h"

using namespace MML;
using namespace MML::TestBeds;

// Get the helix test curve
const auto& helix = ParametricCurvesTestBed::getTestCurve("Helix");

// Test curvature at multiple points
for (Real t = 0; t < 2 * Constants::PI; t += 0.5) {
    Real computed = helix._curve.getCurvature(t);
    Real expected = helix._curvatureFunc(t);  // Analytical value
    
    std::cout << "t=" << t << ": κ_computed=" << computed 
              << ", κ_expected=" << expected << "\n";
}
// Helix has constant curvature κ = 0.5
```

### Testing Curve Torsion

```cpp
#include "test_data/parametric_curves_test_bed.h"

using namespace MML::TestBeds;

const auto& helix = ParametricCurvesTestBed::getTestCurve("Helix");

Real t = Constants::PI / 4;
Real torsion = helix._curve.getTorsion(t);
Real expected = helix._torsionFunc(t);

std::cout << "Torsion at t=" << t << ": " << torsion 
          << " (expected: " << expected << ")\n";
// Helix has constant torsion τ = 0.5
```

### Testing Frenet Frame

```cpp
#include "test_data/parametric_curves_test_bed.h"

using namespace MML::TestBeds;

const auto& circle = ParametricCurvesTestBed::getTestCurve("Circle3DXY");

Real t = Constants::PI / 3;
Vector3Cartesian T, N, B;
circle._curve.getMovingTrihedron(t, T, N, B);

// For a planar curve (circle in XY plane):
// - T is tangent (perpendicular to radius)
// - N points toward center
// - B is constant (perpendicular to plane)

std::cout << "Tangent: " << T << "\n";
std::cout << "Normal: " << N << "\n";
std::cout << "Binormal: " << B << "\n";

// Verify orthonormality
Real dot_TN = ScalarProduct(T, N);  // Should be ~0
Real dot_TB = ScalarProduct(T, B);  // Should be ~0
Real dot_NB = ScalarProduct(N, B);  // Should be ~0
```

### Testing Surface Curvature

```cpp
#include "test_data/parametric_surfaces_test_bed.h"
#include "core/Surfaces.h"

using namespace MML;
using namespace MML::TestBeds;

// Unit sphere: K = 1/R² = 1
const auto& sphere = ParametricSurfaceTestBed::getTestSurface3("Unit Sphere");

Real u = Constants::PI / 4;
Real v = Constants::PI / 3;

Real K = sphere._surface.GaussianCurvature(u, v);
Real H = sphere._surface.MeanCurvature(u, v);

std::cout << "Sphere at (u,v)=(" << u << "," << v << "):\n";
std::cout << "  Gaussian curvature K = " << K << " (expected: 1.0)\n";
std::cout << "  Mean curvature H = " << H << " (expected: -1.0)\n";
```

### Testing Minimal Surfaces

```cpp
#include "test_data/parametric_surfaces_test_bed.h"

using namespace MML::TestBeds;

// Helicoid is a minimal surface (H = 0)
const auto& helicoid = ParametricSurfaceTestBed::getTestSurface3("Helicoid");

// Sample multiple points
for (Real u = 0; u < Constants::PI; u += 0.5) {
    for (Real v = -0.5; v <= 0.5; v += 0.5) {
        Real H = helicoid._surface.MeanCurvature(u, v);
        std::cout << "H(" << u << "," << v << ") = " << H << "\n";
        // All values should be ≈ 0 for minimal surface
    }
}
```

### Testing Developable Surfaces

```cpp
#include "test_data/parametric_surfaces_test_bed.h"

using namespace MML::TestBeds;

// Cylinder and cone are developable (K = 0)
const auto& cylinder = ParametricSurfaceTestBed::getTestSurface3("Cylinder");

Real u = Constants::PI / 4;
Real v = 1.0;

Real K = cylinder._surface.GaussianCurvature(u, v);
bool isDev = cylinder._surface.isDevelopable(u, v);

std::cout << "Cylinder K = " << K << " (expected: 0)\n";
std::cout << "Is developable: " << (isDev ? "yes" : "no") << "\n";
```

### Testing Principal Curvatures

```cpp
#include "test_data/parametric_surfaces_test_bed.h"

using namespace MML::TestBeds;

// Torus has varying principal curvatures
const auto& torus = ParametricSurfaceTestBed::getTestSurface3("Torus");

Real u = Constants::PI / 4;
Real v = Constants::PI / 3;

Real k1, k2;
torus._surface.PrincipalCurvatures(u, v, k1, k2);

Real K = k1 * k2;  // Gaussian curvature
Real H = (k1 + k2) / 2;  // Mean curvature

std::cout << "Principal curvatures: k1=" << k1 << ", k2=" << k2 << "\n";
std::cout << "K = k1*k2 = " << K << "\n";
std::cout << "H = (k1+k2)/2 = " << H << "\n";
```

### Iterating Through All Curves

```cpp
#include "test_data/parametric_curves_test_bed.h"

using namespace MML::TestBeds;

for (int i = 0; i < ParametricCurvesTestBed::getNumTestCurves(); ++i) {
    const auto& curve = ParametricCurvesTestBed::getTestCurve(i);
    
    std::cout << "Curve " << i << ": " << curve._curveName << "\n";
    std::cout << "  Expression: " << curve._curveExpr << "\n";
    std::cout << "  Domain: [" << curve._start << ", " << curve._end << "]\n";
    
    // Test curvature at midpoint
    Real t_mid = (curve._start + curve._end) / 2;
    Real kappa = curve._curve.getCurvature(t_mid);
    std::cout << "  Curvature at midpoint: " << kappa << "\n\n";
}
```

### Iterating Through All Surfaces

```cpp
#include "test_data/parametric_surfaces_test_bed.h"

using namespace MML::TestBeds;

for (int i = 0; i < ParametricSurfaceTestBed::getNumTestSurfaces3(); ++i) {
    const auto& surf = ParametricSurfaceTestBed::getTestSurface3(i);
    
    std::cout << "Surface " << i << ": " << surf._surfaceName << "\n";
    std::cout << "  u ∈ [" << surf._u1 << ", " << surf._u2 << "]\n";
    std::cout << "  v ∈ [" << surf._v1 << ", " << surf._v2 << "]\n";
    
    if (!std::isnan(surf._gaussianCurvature)) {
        std::cout << "  Constant K = " << surf._gaussianCurvature << "\n";
    }
    if (surf._isMinimal) {
        std::cout << "  ** Minimal surface (H=0) **\n";
    }
    if (surf._isDevelopable) {
        std::cout << "  ** Developable (K=0) **\n";
    }
    std::cout << "\n";
}
```

---

## Differential Geometry Background

### Curve Formulas

For a parametric curve **r**(t):

- **Unit tangent**: T = r'/|r'|
- **Principal normal**: N = (r' × r'') × r' / |...|
- **Binormal**: B = T × N
- **Curvature**: κ = |r' × r''| / |r'|³
- **Torsion**: τ = (r' × r'') · r''' / |r' × r''|²

### Surface Formulas

For a parametric surface **r**(u,v):

**First Fundamental Form** (metric):
- E = r_u · r_u
- F = r_u · r_v
- G = r_v · r_v

**Second Fundamental Form** (shape):
- L = r_uu · n
- M = r_uv · n
- N = r_vv · n

**Curvatures**:
- Gaussian: K = (LN - M²) / (EG - F²)
- Mean: H = (EN + GL - 2FM) / (2(EG - F²))
- Principal: k₁, k₂ = H ± √(H² - K)

---

## See Also

- [Curves.h](../../mml/core/Curves.h) - Curve class implementations
- [Surfaces.h](../../mml/core/Surfaces.h) - Surface class implementations
- [curves_tests.cpp](../../tests/core/curves_tests.cpp) - Curve unit tests
- [surfaces_tests.cpp](../../tests/core/surfaces_tests.cpp) - Surface unit tests
- [parametric_curve_analyzer_tests.cpp](../../tests/algorithms/parametric_curve_analyzer_tests.cpp) - Curve analysis tests

