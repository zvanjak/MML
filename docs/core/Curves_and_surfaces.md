# Curves and Surfaces

**Sources:**  
- `mml/core/Curves.h` - Parametric curves (2D and 3D)
- `mml/core/Surfaces.h` - Parametric surfaces

Parametric curves and surfaces provide geometric primitives for modeling, visualization, and differential geometry analysis. MML includes predefined curves/surfaces and comprehensive differential geometry operations.

## Overview

**Curves:**
- **2D Curves**: Planar curves in Cartesian/polar coordinates (circles, spirals, closed curves)
- **3D Curves**: Space curves with Frenet-Serret analysis (lines, helices, twisted curves)

**Surfaces:**
- **Parametric surfaces**: 2-parameter mappings ℝ² → ℝ³ with differential geometry
- **Predefined surfaces**: Planes, quadrics, topological objects

All curves and surfaces support:
- Analytical derivatives (tangent, normal, curvature)
- Arc length computation
- Differential geometry tools (curvatures, fundamental forms)
- Integration with visualization and ODE solvers

---

## Quick Reference

### 2D Curves (Planar)

| Class | Type | Parameters | Description |
|-------|------|------------|-------------|
| `Circle2DCurve` | Cartesian | radius, center | Circle in 2D plane |
| `LogSpiralCurve` | Cartesian | λ (<0), c | Logarithmic spiral |
| `LemniscateCurve` | Cartesian | — | Figure-8 lemniscate |
| `DeltoidCurve` | Cartesian | n | Hypocycloid (deltoid) |
| `AstroidCurve` | Cartesian | c | Four-cusped hypocycloid |
| `EpitrochoidCurve` | Cartesian | R, c, n | Roulette curve |
| `ArchimedeanSpiralCurve` | Cartesian | a | Linear-growth spiral |

### 3D Curves (Space)

| Class | Type | Parameters | Description |
|-------|------|------------|-------------|
| `LineCurve` | Cartesian | pnt, dir, t₁, t₂ | Parametric line segment |
| `Circle3DXY`, `Circle3DXZ`, `Circle3DYZ` | Cartesian | radius | Axis-aligned circles |
| `Circle` | Cartesian | R, normal, center | General 3D circle |
| `HelixCurve` | Cartesian | R, b | Circular helix |
| `TwistedCubicCurve` | Cartesian | — | Polynomial space curve |
| `ToroidalSpiralCurve` | Cartesian | n, scale | Toroidal spiral |

### Surfaces

| Class | Parameters | Description |
|-------|------------|-------------|
| `PlaneSurface` | point, normal, axes | Plane with local coordinates |
| `Cylinder` | R, H | Circular cylinder |
| `Torus` | R (major), r (minor) | Torus (doughnut) |
| `Sphere` | R | Sphere |
| `Ellipsoid` | a, b, c | Triaxial ellipsoid |
| `Hyperboloid` | a, b, c | One-sheet hyperboloid |
| `Paraboloid` | a, h | Circular paraboloid |
| `MobiusStrip` | — | Non-orientable surface |
| `MonkeySaddle` | — | Saddle with 3-fold symmetry |

---

## 2D Curves (Planar)

### Base Class: ICurveCartesian2D

All 2D curves implement `IParametricCurve<2>` and provide:

**Methods:**
```cpp
Vec2Cart getTangent(Real t);          // First derivative
Vec2Cart getTangentUnit(Real t);      // Unit tangent
Vec2Cart getNormal(Real t);           // Second derivative
Vec2Cart getNormalUnit(Real t);       // Unit normal
```

### Circle2DCurve

Circle in 2D Cartesian plane.

**Constructors:**
```cpp
Circle2DCurve();                               // Unit circle at origin
Circle2DCurve(Real radius);                     // Circle at origin
Circle2DCurve(Real radius, const Pnt2Cart& center);
```

**Parameterization:** t ∈ [0, 2π]  
**Equation:**  
- x(t) = x₀ + R cos(t)
- y(t) = y₀ + R sin(t)

**Example:**
```cpp
// Circle of radius 5 at (2, 3)
Circle2DCurve circle(5.0, Pnt2Cart(2.0, 3.0));

Vec2Cart point = circle(Constants::PI / 2);  // top point
Vec2Cart tangent = circle.getTangent(0.0);    // tangent at t=0
```

### LogSpiralCurve

Logarithmic spiral - exponentially growing radius.

**Constructors:**
```cpp
LogSpiralCurve();                      // Default λ=-1, c=1
LogSpiralCurve(Real lambda);           // λ < 0 required
LogSpiralCurve(Real lambda, Real c);   // c ≠ 0 required
```

**Parameterization:** t ∈ (-∞, ∞)  
**Equation:**  
- r(t) = c·e^(λt)
- x(t) = r(t)·cos(t) = c·e^(λt)·cos(t)
- y(t) = r(t)·sin(t) = c·e^(λt)·sin(t)

**Properties:**
- Crosses any ray from origin at constant angle
- λ < 0: spiral inward as t → ∞
- |λ| larger: tighter spiral

**Example:**
```cpp
LogSpiralCurve spiral(-0.1, 2.0);  // Slow inward spiral

for (Real t = 0; t < 10.0; t += 0.5) {
    Vec2Cart pt = spiral(t);
    // Plot points...
}
```

### LemniscateCurve

Figure-8 curve (lemniscate of Gerono).

**Parameterization:** t ∈ (-∞, ∞)  
**Equation:**  
- x(t) = cos(t) / (1 + sin²(t))
- y(t) = sin(t)·cos(t) / (1 + sin²(t))

**Properties:**
- Symmetric about both axes
- Self-intersecting at origin
- Periodic with period 2π

### DeltoidCurve

Hypocycloid with 3 cusps (deltoid).

**Constructor:**
```cpp
DeltoidCurve();       // n=1
DeltoidCurve(int n);  // Scaled version
```

**Equation:**  
- x(t) = 2n·cos(t)·(1 + cos(t))
- y(t) = 2n·sin(t)·(1 - cos(t))

**Example:**
```cpp
DeltoidCurve deltoid(2);  // Larger deltoid

Vec2Cart pt = deltoid(Constants::PI / 3);
```

### AstroidCurve

Four-cusped hypocycloid (astroid).

**Constructor:**
```cpp
AstroidCurve();       // c=1
AstroidCurve(Real c); // c > 0 required
```

**Equation:**  
- x(t) = c·cos³(t)
- y(t) = c·sin³(t)

**Properties:**
- Four cusps at (±c, 0) and (0, ±c)
- Parametric form of x^(2/3) + y^(2/3) = c^(2/3)

### EpitrochoidCurve

Roulette curve traced by point on rolling circle.

**Constructor:**
```cpp
EpitrochoidCurve();                        // R=1, c=1, n=1
EpitrochoidCurve(Real R, Real c, int n);
```

**Equation:**  
- x(t) = cos(t) - c·cos(n·t)
- y(t) = sin(t) - c·sin(n·t)

**Example:**
```cpp
EpitrochoidCurve epi(1.0, 0.5, 3);  // 3-lobed pattern
```

### ArchimedeanSpiralCurve

Linear-growth spiral.

**Constructor:**
```cpp
ArchimedeanSpiralCurve();       // a=1
ArchimedeanSpiralCurve(Real a);
```

**Equation:**  
- r(t) = a·t
- x(t) = a·t·cos(t)
- y(t) = a·t·sin(t)

**Properties:**
- Constant spacing between turns
- Used in scroll compressors, watch springs

---

## 3D Curves (Space)

### Base Class: ICurveCartesian3D

All 3D curves implement `IParametricCurve<3>` with comprehensive differential geometry.

**Methods:**
```cpp
// Frenet-Serret frame
Vec3Cart getTangent(Real t) const;        // First derivative r'(t)
Vec3Cart getTangentUnit(Real t) const;    // Unit tangent T = r'/|r'|
Vec3Cart getNormal(Real t) const;         // Second derivative r''(t)
Vec3Cart getNormalUnit(Real t) const;     // Principal normal N
Vec3Cart getBinormal(Real t) const;       // Binormal B = T × N

// Curvature and torsion
Vec3Cart getCurvatureVector(Real t) const;  // Curvature vector κN
Real getCurvature(Real t) const;            // Curvature κ = |r' × r''| / |r'|³
Real getTorsion(Real t) const;              // Torsion τ (twist)

// Osculating planes
Plane3D getOsculationPlane(Real t) const;   // Contains T and N
Plane3D getNormalPlane(Real t) const;       // Perpendicular to T
Plane3D getRectifyingPlane(Real t) const;   // Contains T and B

// Moving frame
void getMovingTrihedron(Real t, Vec3Cart& T, Vec3Cart& N, Vec3Cart& B);

// Arc length check
bool isArcLengthParametrized(Real t1, Real t2, int numPoints = 100) const;
```

**Frenet-Serret Formulas:**
- dT/ds = κN
- dN/ds = -κT + τB
- dB/ds = -τN

where s is arc length, κ is curvature, τ is torsion.

### LineCurve

Parametric line segment.

**Constructors:**
```cpp
// Direct specification
LineCurve(Real minT, Real maxT, const Point3Cartesian& pnt, 
          const Vector3Cartesian& dir);

// Two-point specification (t₁ → pnt₁, t₂ → pnt₂)
LineCurve(Real t1, const Point3Cartesian& pnt1, 
          Real t2, const Point3Cartesian& pnt2);
```

**Example:**
```cpp
// Line from (0,0,0) to (1,1,1) over t ∈ [0,1]
LineCurve line(0.0, Pnt3Cart(0,0,0), 1.0, Pnt3Cart(1,1,1));

Vec3Cart mid = line(0.5);  // Midpoint
Real curv = line.getCurvature(0.5);  // 0 (straight line)
```

### Circle3DXY, Circle3DXZ, Circle3DYZ

Axis-aligned circles in coordinate planes.

**Example:**
```cpp
Circle3DXY circleXY(3.0);  // Radius 3 in XY plane

Vec3Cart pt = circleXY(Constants::PI);  // (-3, 0, 0)
```

**Parameter:** t ∈ [0, 2π]

### Circle (General 3D Circle)

Circle in arbitrary plane.

**Constructor:**
```cpp
Circle(Real radius, const Vec3Cart& normal, const Pnt3Cart& center);
```

**Properties:**
- Plane defined by center and normal vector
- Automatic computation of orthonormal frame

**Example:**
```cpp
// Circle in tilted plane
Vec3Cart normal(1, 1, 1);  // Will be normalized
Pnt3Cart center(0, 0, 5);
Circle circle(2.0, normal, center);

Vec3Cart pt = circle(Constants::PI / 4);
```

### HelixCurve

Circular helix - constant radius and pitch.

**Constructor:**
```cpp
HelixCurve();                  // R=1, b=1
HelixCurve(Real radius, Real b);
```

**Parameterization:** t ∈ (-∞, ∞)  
**Equation:**  
- x(t) = R·cos(t)
- y(t) = R·sin(t)
- z(t) = b·t

**Properties:**
- Curvature: κ = R / (R² + b²) (constant!)
- Torsion: τ = b / (R² + b²) (constant!)
- Pitch: 2πb (vertical distance per revolution)

**Example:**
```cpp
HelixCurve helix(2.0, 0.5);  // Radius 2, slow pitch

Real kappa = helix.getCurvature(0.0);  // 2/(4+0.25) = 0.47
Real tau = helix.getTorsion(0.0);      // 0.5/(4+0.25) = 0.12
```

### TwistedCubicCurve

Polynomial space curve.

**Parameterization:** t ∈ (-∞, ∞)  
**Equation:**  
- x(t) = t
- y(t) = t²
- z(t) = t³

**Properties:**
- Non-planar for all t
- Osculating plane changes continuously
- Used in graphics and numerical examples

### ToroidalSpiralCurve

Spiral on torus surface.

**Constructor:**
```cpp
ToroidalSpiralCurve();                    // n=1, scale=1
ToroidalSpiralCurve(int n);                // Winding number
ToroidalSpiralCurve(int n, Real scale);
```

**Equation:**  
- x(t) = scale·(4 + sin(n·t))·cos(t)
- y(t) = scale·(4 + sin(n·t))·sin(t)
- z(t) = scale·cos(n·t)

**Properties:**
- Winds n times around torus
- scale adjusts overall size

---

## Surfaces

### Base Class: ISurfaceCartesian

All surfaces implement `IParametricSurfaceRect<3>` with rich differential geometry.

**Methods:**
```cpp
// Geometry
VectorN<Real, 3> Normal(Real u, Real w) const;  // Unit normal
void Tangents(Real u, Real w, Vec3Cart& tU, Vec3Cart& tW) const;

// Curvatures
Real GaussianCurvature(Real u, Real w);    // K = κ₁·κ₂
Real MeanCurvature(Real u, Real w);        // H = (κ₁ + κ₂)/2
void PrincipalCurvatures(Real u, Real w, Real& k1, Real& k2) const;
void PrincipalDirections(Real u, Real w, Vec3Cart& dir1, Vec3Cart& dir2) const;

// Fundamental forms
void GetFirstNormalFormCoefficients(Real u, Real w, Real& E, Real& F, Real& G);
void GetSecondNormalFormCoefficients(Real u, Real w, Real& L, Real& M, Real& N);

// Classification
bool isRegular(Real u, Real w) const;     // Non-degenerate?
bool isFlat(Real u, Real w, Real eps);    // K ≈ 0?
bool isParabolic(Real u, Real w, Real eps); // K ≈ 0 but not flat?
bool isHyperbolic(Real u, Real w, Real eps); // K < 0?
```

**First Fundamental Form:**
- ds² = E du² + 2F du dw + G dw²
- E = r_u · r_u, F = r_u · r_w, G = r_w · r_w
- Encodes intrinsic metric (distances, angles)

**Second Fundamental Form:**
- Encodes extrinsic curvature
- L = r_uu · n, M = r_uw · n, N = r_ww · n

**Gaussian Curvature:**
- K = (LN - M²) / (EG - F²)
- Intrinsic invariant (Gauss's Theorema Egregium)

**Mean Curvature:**
- H = (EN + GL - 2FM) / (2(EG - F²))
- Extrinsic measure

### PlaneSurface

Plane with local (u,v) coordinate system.

**Constructors:**
```cpp
// Full specification
PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal,
             const Vec3Cart& uAxis, const Vec3Cart& vAxis);

// Auto-compute v axis
PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal,
             const Vec3Cart& uAxisDirection);

// Auto-compute both axes
PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal);
```

**Parameterization:** u, w ∈ (-∞, ∞) (can be limited)  
**Equation:** **r**(u,w) = **p** + u·**e**_u + w·**e**_w

**Example:**
```cpp
// XY plane
PlaneSurface plane(Vec3Cart(0,0,0), Vec3Cart(0,0,1));

// Tilted plane through (1,2,3)
Vec3Cart normal(1, 1, 1);
Vec3Cart point(1, 2, 3);
PlaneSurface tilted(point, normal);

Vec3Cart pt = tilted(1.0, 2.0);
Vec3Cart n = tilted.Normal(1.0, 2.0);  // Same as input normal
Real K = tilted.GaussianCurvature(0, 0);  // 0 (flat!)
```

### Cylinder

Circular cylinder.

**Constructor:**
```cpp
Cylinder();           // R=1, H=1
Cylinder(Real R, Real H);
```

**Parameterization:**  
- u ∈ [0, 2π] - angle
- w ∈ [0, H] - height

**Equation:**  
- x(u,w) = R·cos(u)
- y(u,w) = R·sin(u)
- z(u,w) = w

**Properties:**
- Gaussian curvature K = 0 (developable!)
- Mean curvature H = 1/(2R)
- Principal curvatures: κ₁ = 1/R, κ₂ = 0

**Example:**
```cpp
Cylinder cyl(2.0, 5.0);  // Radius 2, height 5

Vec3Cart pt = cyl(Constants::PI, 2.5);  // Point at u=π, w=2.5
Real K = cyl.GaussianCurvature(0, 0);   // 0
Real H = cyl.MeanCurvature(0, 0);        // 0.25
```

### Torus

Torus (doughnut shape).

**Constructor:**
```cpp
Torus();               // R=1 (major), r=0.5 (minor)
Torus(Real R, Real r); // R > r for standard torus
```

**Parameterization:**  
- u, w ∈ [0, 2π]

**Equation:**  
- x(u,w) = (R + r·cos(w))·cos(u)
- y(u,w) = (R + r·cos(w))·sin(u)
- z(u,w) = r·sin(w)

**Properties:**
- Gaussian curvature: K = cos(w) / (r(R + r·cos(w)))
- Variable curvature: K > 0 (outer), K < 0 (inner)

**Example:**
```cpp
Torus torus(3.0, 1.0);  // Major radius 3, minor radius 1

// Outer rim (u=0, w=0)
Real K_outer = torus.GaussianCurvature(0, 0);  // Positive

// Inner rim (u=0, w=π)
Real K_inner = torus.GaussianCurvature(0, Constants::PI);  // Negative
```

### Sphere

Sphere.

**Constructor:**
```cpp
Sphere();       // R=1
Sphere(Real R);
```

**Parameterization:**  
- u ∈ [0, π] - polar angle (latitude)
- w ∈ [0, 2π] - azimuthal angle (longitude)

**Equation:**  
- x(u,w) = R·sin(u)·cos(w)
- y(u,w) = R·sin(u)·sin(w)
- z(u,w) = R·cos(u)

**Properties:**
- Gaussian curvature: K = 1/R² (constant!)
- Mean curvature: H = 1/R
- All points umbilical (κ₁ = κ₂)

**Example:**
```cpp
Sphere sphere(5.0);

Real K = sphere.GaussianCurvature(Constants::PI/2, 0);  // 0.04 = 1/25
Real H = sphere.MeanCurvature(Constants::PI/2, 0);       // 0.2 = 1/5
```

### Ellipsoid

Triaxial ellipsoid.

**Constructor:**
```cpp
Ellipsoid();              // a=b=c=1 (sphere)
Ellipsoid(Real a, Real b, Real c);
```

**Equation:**  
- x(u,w) = a·sin(u)·cos(w)
- y(u,w) = b·sin(u)·sin(w)
- z(u,w) = c·cos(u)

**Example:**
```cpp
Ellipsoid ellipsoid(3.0, 2.0, 1.0);  // Flattened

Vec3Cart pt = ellipsoid(Constants::PI/2, 0);  // (3, 0, 0)
```

### Hyperboloid

One-sheet hyperboloid.

**Constructor:**
```cpp
Hyperboloid();              // a=b=c=1
Hyperboloid(Real a, Real b, Real c);
```

**Equation:**  
- x(u,w) = a·cosh(u)·cos(w)
- y(u,w) = b·cosh(u)·sin(w)
- z(u,w) = c·sinh(u)

**Properties:**
- Negative Gaussian curvature (saddle-like)
- Ruled surface

### Paraboloid

Circular paraboloid.

**Constructor:**
```cpp
Paraboloid();           // a=1, h=1
Paraboloid(Real a, Real h);
```

**Equation:**  
- x(u,w) = a·√(u/h)·cos(w)
- y(u,w) = a·√(u/h)·sin(w)
- z(u,w) = u

### MobiusStrip

Non-orientable surface (one-sided).

**Parameterization:**  
- u ∈ [0, 2π]
- w ∈ [-1, 1]

**Equation:**  
- x(u,w) = (1 + w·cos(u/2))·cos(u)
- y(u,w) = (1 + w·cos(u/2))·sin(u)
- z(u,w) = w·sin(u/2)

**Properties:**
- Non-orientable (one-sided)
- Single edge
- Topologically interesting

### MonkeySaddle

Saddle surface with 3-fold rotational symmetry.

**Parameterization:** u, w ∈ [-10, 10]  
**Equation:**  
- x(u,w) = u
- y(u,w) = w
- z(u,w) = u·(u² - 3w²)

**Properties:**
- Saddle point at origin
- Three "downward" directions (like monkey with tail)

---

## Differential Geometry Operations

### Arc Length

**For curves:**
```cpp
// Using PathIntegration
Real length = PathIntegration::ParametricCurveLength(curve, t1, t2);
```

**Example:**
```cpp
HelixCurve helix(1.0, 1.0);

// Length of one full turn
Real len = PathIntegration::ParametricCurveLength(helix, 0.0, 2*Constants::PI);
// Should be 2π√(R² + b²) = 2π√2 ≈ 8.886
```

### Curvature and Torsion

**For curves:**
```cpp
Real kappa = curve.getCurvature(t);   // κ = |r' × r''| / |r'|³
Real tau = curve.getTorsion(t);       // τ from Frenet-Serret

// For helix (analytical):
HelixCurve helix(2.0, 1.0);
Real kappa = helix.getCurvature(0);  // Overridden for efficiency
Real tau = helix.getTorsion(0);
```

**Physical interpretation:**
- Curvature κ: How much curve bends (1/radius for circle)
- Torsion τ: How much curve twists out of plane

### Fundamental Forms

**For surfaces:**
```cpp
Real E, F, G;  // First fundamental form
surf.GetFirstNormalFormCoefficients(u, w, E, F, G);

Real L, M, N;  // Second fundamental form
surf.GetSecondNormalFormCoefficients(u, w, L, M, N);
```

**Example:**
```cpp
Sphere sphere(5.0);

Real E, F, G;
sphere.GetFirstNormalFormCoefficients(Constants::PI/2, 0, E, F, G);
// E = R²sin²(u) = 25, F = 0, G = R² = 25

Real L, M, N;
sphere.GetSecondNormalFormCoefficients(Constants::PI/2, 0, L, M, N);
// L = R = 5, M = 0, N = R·sin²(u) = 5
```

### Principal Curvatures

**For surfaces:**
```cpp
Real k1, k2;
surf.PrincipalCurvatures(u, w, k1, k2);

// Principal curvature directions
Vec3Cart dir1, dir2;
surf.PrincipalDirections(u, w, dir1, dir2);
```

**Relationships:**
- Gaussian curvature: K = κ₁·κ₂
- Mean curvature: H = (κ₁ + κ₂)/2
- Umbilical point: κ₁ = κ₂ (sphere has all points umbilical)

**Example:**
```cpp
Torus torus(3.0, 1.0);

Real k1, k2;
torus.PrincipalCurvatures(0, 0, k1, k2);  // Outer rim

std::cout << "κ₁ = " << k1 << ", κ₂ = " << k2 << "\n";
std::cout << "K = " << k1 * k2 << "\n";
std::cout << "H = " << (k1 + k2) / 2 << "\n";
```

---

## Examples

### Example 1: Helix Analysis

```cpp
#include "core/Curves.h"
#include "core/Integration/PathIntegration.h"

void Example_Helix_Analysis()
{
    using namespace Curves;

    // Create helix: radius 2, pitch 1
    HelixCurve helix(2.0, 1.0);

    Real t = Constants::PI / 2;

    // Position
    Vec3Cart pos = helix(t);
    std::cout << "Position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";

    // Frenet-Serret frame
    Vec3Cart T = helix.getTangentUnit(t);
    Vec3Cart N = helix.getNormalUnit(t);
    Vec3Cart B = helix.getBinormal(t);

    std::cout << "Tangent:  (" << T[0] << ", " << T[1] << ", " << T[2] << ")\n";
    std::cout << "Normal:   (" << N[0] << ", " << N[1] << ", " << N[2] << ")\n";
    std::cout << "Binormal: (" << B[0] << ", " << B[1] << ", " << B[2] << ")\n";

    // Curvature and torsion
    Real kappa = helix.getCurvature(t);
    Real tau = helix.getTorsion(t);

    std::cout << "Curvature: " << kappa << "\n";  // 2/(4+1) = 0.4
    std::cout << "Torsion:   " << tau << "\n";    // 1/(4+1) = 0.2

    // Arc length of one revolution
    Real length = PathIntegration::ParametricCurveLength(helix, 0.0, 2*Constants::PI);
    std::cout << "Length (one turn): " << length << "\n";
    // Should be 2π√(R²+b²) = 2π√5 ≈ 14.05
}
```

### Example 2: Torus Curvature Map

```cpp
#include "core/Surfaces.h"

void Example_Torus_Curvature()
{
    using namespace Surfaces;

    Torus torus(3.0, 1.0);  // Major radius 3, minor radius 1

    // Sample parameter space
    int nU = 20, nW = 20;
    
    for (int i = 0; i < nU; ++i) {
        for (int j = 0; j < nW; ++j) {
            Real u = i * 2*Constants::PI / nU;
            Real w = j * 2*Constants::PI / nW;

            Vec3Cart pos = torus(u, w);
            Real K = torus.GaussianCurvature(u, w);
            Real H = torus.MeanCurvature(u, w);

            // Classify point
            std::string type;
            if (K > 0.01) type = "elliptic (bowl-like)";
            else if (K < -0.01) type = "hyperbolic (saddle)";
            else type = "parabolic (ridge/valley)";

            std::cout << "u=" << u << ", w=" << w 
                      << ": K=" << K << ", H=" << H 
                      << " [" << type << "]\n";
        }
    }
}
```

### Example 3: Surface Normals for Lighting

```cpp
#include "core/Surfaces.h"
#include "tools/Visualizer.h"

void Example_Surface_Normals()
{
    using namespace Surfaces;

    Sphere sphere(5.0);

    // Sample sphere surface
    int nU = 30, nW = 30;
    std::vector<Vec3Cart> points, normals;

    for (int i = 0; i < nU; ++i) {
        for (int j = 0; j < nW; ++j) {
            Real u = i * Constants::PI / nU;
            Real w = j * 2*Constants::PI / nW;

            Vec3Cart pos = sphere(u, w);
            Vec3Cart normal = sphere.Normal(u, w);

            points.push_back(pos);
            normals.push_back(normal);
        }
    }

    // Use normals for shading/lighting calculations
    // ...
}
```

### Example 4: Principal Directions on Ellipsoid

```cpp
#include "core/Surfaces.h"

void Example_Ellipsoid_PrincipalDirections()
{
    using namespace Surfaces;

    Ellipsoid ellipsoid(3.0, 2.0, 1.0);  // Triaxial

    Real u = Constants::PI / 4;
    Real w = Constants::PI / 6;

    // Principal curvatures
    Real k1, k2;
    ellipsoid.PrincipalCurvatures(u, w, k1, k2);

    std::cout << "Principal curvatures: κ₁=" << k1 << ", κ₂=" << k2 << "\n";
    std::cout << "Gaussian curvature K = " << k1 * k2 << "\n";
    std::cout << "Mean curvature H = " << (k1 + k2)/2 << "\n";

    // Principal directions
    Vec3Cart dir1, dir2;
    ellipsoid.PrincipalDirections(u, w, dir1, dir2);

    std::cout << "Direction 1: (" << dir1[0] << ", " << dir1[1] << ", " << dir1[2] << ")\n";
    std::cout << "Direction 2: (" << dir2[0] << ", " << dir2[1] << ", " << dir2[2] << ")\n";

    // Verify orthogonality
    Real dot = ScalarProduct(dir1, dir2);
    std::cout << "Orthogonality check (should be ~0): " << dot << "\n";
}
```

### Example 5: Custom Curve from Lambda

```cpp
#include "core/Curves.h"

void Example_Custom_Curve()
{
    using namespace Curves;

    // Define custom curve with lambda
    auto lissajous = [](Real t) -> VectorN<Real, 3> {
        return VectorN<Real, 3>{
            sin(3*t),
            sin(4*t),
            sin(5*t)
        };
    };

    CurveCartesian3D curve(0.0, 2*Constants::PI, lissajous);

    // Analyze at t = π/2
    Real t = Constants::PI / 2;
    
    Vec3Cart pos = curve(t);
    Vec3Cart tangent = curve.getTangent(t);
    Real kappa = curve.getCurvature(t);

    std::cout << "Position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
    std::cout << "Curvature: " << kappa << "\n";
}
```

### Example 6: Plane Through 3 Points

```cpp
#include "core/Surfaces.h"

void Example_Plane_From_Points()
{
    using namespace Surfaces;

    // Three non-collinear points
    Vec3Cart p1(1, 0, 0);
    Vec3Cart p2(0, 1, 0);
    Vec3Cart p3(0, 0, 1);

    // Compute normal: (p2-p1) × (p3-p1)
    Vec3Cart v1 = p2 - p1;
    Vec3Cart v2 = p3 - p1;
    Vec3Cart normal = VectorProduct(v1, v2).Normalized();

    // Create plane
    PlaneSurface plane(p1, normal, v1.Normalized());

    // Verify points are on plane
    Vec3Cart test1 = plane(0, 0);  // Should be near p1
    Vec3Cart test2 = plane(1, 0);  // Should be along v1 direction

    std::cout << "Point on plane: (" << test1[0] << ", " << test1[1] << ", " << test1[2] << ")\n";
}
```

---

## Integration with MML

### With ODE Solvers

```cpp
// ODE solution → parametric curve
ODESystemSolution solution = solver.solve(ode_system);
auto curve_from_solution = solution.toParametricCurve();

// Analyze solution curve
Real curvature = curve_from_solution.getCurvature(t);
```

### With Interpolation

```cpp
// Data points → interpolated curve
Matrix<Real> points = getSampledData();
SplineInterpParametricCurve<3> curve(points);

// Differential geometry on interpolated curve
Vec3Cart tangent = Derivation::DeriveCurve<3>(curve, t, nullptr);
```

### With Visualization

```cpp
HelixCurve helix(2.0, 1.0);

// Visualize 3D curve
Visualizer::VisualizeParamCurve3D(helix, "helix", 0.0, 10.0, 200);

// Visualize surface
Torus torus(3.0, 1.0);
Visualizer::VisualizeSurface3D(torus, "torus", 0, 2*Constants::PI, 50, 0, 2*Constants::PI, 50);
```

### With Field Operations

```cpp
// Surface → scalar field (curvature visualization)
Sphere sphere(5.0);

auto curvature_field = [&](const Vec3Cart& pos) {
    // Map position to parameters, return curvature
    // ...
};
```

---

## Runnable Examples

### Available Demos

| Example | Source File | Description |
|---------|-------------|-------------|
| Differential Geometry | [docs_demo_diff_geometry.cpp](../../src/docs_demos/docs_demo_diff_geometry.cpp) | Comprehensive curves and surfaces demo |
| Parametric Curves | [docs_demo_interpolated_functions.cpp](../../src/docs_demos/docs_demo_interpolated_functions.cpp) | Spline-interpolated curves |
| Curve Length | [docs_demo_integration_path.cpp](../../src/docs_demos/docs_demo_integration_path.cpp) | Arc length computation |

### Demo Functions in docs_demo_diff_geometry.cpp

**2D Curves:**
- `Docs_Demo_Curves_2D_Basic()` - Circle2D, LogSpiral, Lemniscate, Deltoid, Astroid, ArchimedeanSpiral

**3D Curves:**
- `Docs_Demo_Curves_3D_Basic()` - LineCurve, Circle3D variants, general Circle, HelixCurve, TwistedCubic, ToroidalSpiral

**Curve Differential Geometry:**
- `Docs_Demo_Curves_Frenet_Frame()` - Tangent, normal, binormal vectors; moving trihedron
- `Docs_Demo_Curves_Curvature_Torsion()` - Curvature κ, torsion τ for various curves
- `Docs_Demo_Curves_Arc_Length()` - Arc length calculation using path integration
- `Docs_Demo_Curves_Planes()` - Osculating, normal, and rectifying planes

**Surface Properties:**
- `Docs_Demo_Surfaces_Basic()` - Sphere, Cylinder, Torus basic properties
- `Docs_Demo_Surfaces_Curvature()` - Gaussian K and mean H curvature
- `Docs_Demo_Surfaces_Principal_Curvatures()` - Principal curvatures k₁, k₂ and directions
- `Docs_Demo_Surfaces_Fundamental_Forms()` - First (E, F, G) and second (L, M, N) fundamental forms
- `Docs_Demo_Surfaces_Special()` - Möbius strip, Paraboloid, PlaneSurface
- `Docs_Demo_Surfaces_Tangents_Normals()` - Surface tangent vectors and normals
- `Docs_Demo_Surfaces_Classification()` - Elliptic, hyperbolic, parabolic point classification

### Sample Output

```cpp
// Helix curvature and torsion
HelixCurve helix(2.0, 1.0);  // R=2, b=1
// Analytical: κ = R/(R² + b²) = 2/5 = 0.4
// Analytical: τ = b/(R² + b²) = 1/5 = 0.2
std::cout << helix.getCurvature(0.0);  // 0.4
std::cout << helix.getTorsion(0.0);    // 0.2

// Sphere curvatures
Sphere sphere(5.0);  // R=5
// K = 1/R² = 0.04, H = 1/R = 0.2
sphere.GaussianCurvature(PI/2, 0);  // 0.04
sphere.MeanCurvature(PI/2, 0);      // 0.2

// Principal curvatures (sphere: k1 = k2 = 1/R)
Real k1, k2;
sphere.PrincipalCurvatures(PI/2, 0, k1, k2);  // k1 = k2 = 0.2
```

### Building and Running

```bash
# Build
cd build && cmake --build .

# Run all docs demos
./src/docs_demos/Debug/MML_DocsApp.exe  # Windows
./src/docs_demos/MML_DocsApp             # Linux/Mac
```

- [Functions.md](Functions.md) - Function interfaces (IParametricCurve, IParametricSurface)
- [Interpolated_functions.md](Interpolated_functions.md) - Spline curves from data
- [ODE_system.md](ODE_system.md) - Solution curves as parametric curves
- [Derivation.md](Derivation.md) - Numerical derivatives for curves/surfaces
- [Integration/PathIntegration.md](Integration/PathIntegration.md) - Arc length integrals
- [Coordinate_transformations.md](Coordinate_transformations.md) - Coordinate systems
- [../tools/Visualizers_README.md](../tools/Visualizers_README.md) - 3D visualization

---

**References:**
- do Carmo, *Differential Geometry of Curves and Surfaces*, Prentice Hall, 1976
- O'Neill, *Elementary Differential Geometry* (2nd ed.), Academic Press, 2006
- Pressley, *Elementary Differential Geometry* (2nd ed.), Springer, 2010
