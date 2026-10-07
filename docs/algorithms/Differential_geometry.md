# Differential Geometry

Comprehensive toolkit for **differential geometry of curves and surfaces** - from Frenet-Serret frames to curvature analysis.

## Overview

**Covers**:
- **Parametric Curves**: Tangent, normal, binormal, curvature, torsion
- **Frenet-Serret Frame**: Moving trihedron along curves
- **Osculating Geometry**: Osculating plane, normal plane, rectifying plane
- **Surface Geometry**: First/second fundamental forms, Gaussian curvature, intrinsic induced metrics
- **Exterior Calculus**: Typed differential forms, exterior derivative, form integration, Green/Stokes residuals
- **Charts and Frames**: Single-chart embeddings, induced metric, pullback/push-forward, Darboux frames

**Architecture**: classical Cartesian curve/surface methods remain in `mml/core/Curves.h` and `mml/core/Surfaces.h`. The newer coordinate-aware layer lives in `mml/core/DifferentialGeometry/` and builds on `MetricTensorField`, typed forms, geodesic solvers, and numerical integration.

## Current Core Modules

| Header | Role |
|--------|------|
| `mml/base/DifferentialGeometry/DifferentialForm.h` | Pointwise differential forms and wedge products. |
| `mml/base/DifferentialGeometry/Hodge.h` | Hodge star, orientation, and metric-derived cross product. |
| `mml/base/DifferentialGeometry/Metric.h` | Typed metric value object plus `flat`/`sharp` musical maps. |
| `mml/base/DifferentialGeometry/TypedFields.h` | Typed scalar/vector/form field interfaces. |
| `mml/core/DifferentialGeometry/FieldOperations.h` | Exterior derivative, metric gradient, typed divergence/curl wrappers. |
| `mml/core/DifferentialGeometry/InducedMetric.h` | Surface first fundamental form as a `MetricTensorField<2>`. |
| `mml/core/DifferentialGeometry/FormIntegration.h` | One-form and two-form integration plus Green/Stokes residuals. |
| `mml/core/DifferentialGeometry/Chart.h` | Single-chart embeddings, Jacobians, induced metrics, push-forward, pullback. |
| `mml/core/DifferentialGeometry/Atlas.h` | Minimal same-frame chart collections and the two-chart stereographic unit sphere atlas. |
| `mml/core/DifferentialGeometry/Frames.h` | Generalized Frenet frame and surface-curve Darboux frame helpers. |

## Organization Notes

- `mml/base/Geometry/` is now for concrete Euclidean geometry primitives.
- `mml/base/DifferentialGeometry/` owns pointwise forms, metrics, frame tags, and typed field interfaces.
- `mml/base/Algebra/LieGroups/` owns `SO2`, `SO3`, `SE2`, and `SE3`.
- `mml/core/DifferentialGeometry/` owns runtime algorithms such as charts, induced metrics, form integration, and frame helpers.
- The old differential-form and Lie-group wrapper paths under `mml/base/Geometry/` have been removed; include the new physical homes directly.

## Quick Reference

### Curve Geometry (3D)

| Property | Symbol | Formula | Geometric Meaning |
|----------|--------|---------|-------------------|
| **Tangent** | **T** | r'(t) | Velocity vector, direction of motion |
| **Unit Tangent** | **T̂** | r'(t)/\|r'(t)\| | Normalized direction |
| **Curvature** | **κ** | \|r' × r''\|/\|r'\|³ | Rate of direction change |
| **Principal Normal** | **N** | T'/\|T'\| | Points toward center of curvature |
| **Binormal** | **B** | T × N | Perpendicular to osculating plane |
| **Torsion** | **τ** | (r'×r'')·r'''/\|r'×r''\|² | Rate of plane twist |

**Frenet-Serret Frame**: Orthonormal basis {**T**, **N**, **B**} following curve

### Common Curves

| Curve | Curvature κ | Torsion τ | Properties |
|-------|------------|-----------|------------|
| **Straight Line** | 0 | 0 | No bending, no twisting |
| **Circle (radius R)** | 1/R | 0 | Constant curvature, planar |
| **Helix** | const | const | Constant κ and τ (space curve) |
| **Parabola** | varies | 0 | Variable curvature, planar |

---

## Mathematical Background

### Parametric Curves

**Definition**: **r**: ℝ → ℝ³, **r**(t) = (x(t), y(t), z(t))

**Examples**:
```
Circle:  r(t) = (R cos t, R sin t, 0)
Helix:   r(t) = (a cos t, a sin t, b t)
Line:    r(t) = r₀ + t v
```

### Tangent Vector

**Definition**: **r'**(t) = d**r**/dt

**Physical meaning**: Velocity vector at point **r**(t)

**Formula**:
```
r'(t) = (dx/dt, dy/dt, dz/dt)
```

**Unit tangent**:
```
T̂(t) = r'(t) / |r'(t)|
```

**Direction**: Along curve in direction of increasing t

### Curvature

**Intuition**: How much curve bends (deviates from straight line)

**Formula**:
```
κ = |r' × r''| / |r'|³
```

**Alternative** (via unit tangent):
```
κ = |dT̂/ds|
```
where s = arc length

**Geometric meaning**:
- κ = 0: Straight line (no bending)
- κ = 1/R: Circle of radius R
- Large κ: Sharp turn
- Small κ: Gentle curve

**Radius of curvature**: ρ = 1/κ

### Curvature Vector

**Definition**: **κ** = d**T̂**/ds (derivative w.r.t. arc length)

**Points**: Toward center of curvature (concave side)

**Magnitude**: Curvature κ

**Formula** (parameter t):
```
κ_vec = (r'' - (r'·r''/|r'|²) r') / |r'|²
```

### Principal Normal

**Definition**: Direction of curvature vector

**Formula**:
```
N̂ = dT̂/ds / |dT̂/ds| = κ_vec / κ
```

**Alternative** (via cross products):
```
N̂ = (r' × (r'' × r')) / (|r'| |r'' × r'|)
```

**Points**: Toward concave side (center of osculating circle)

**Property**: Perpendicular to **T̂** (by definition of derivative along curve)

### Binormal

**Definition**: **B** = **T̂** × **N̂**

**Properties**:
- Perpendicular to both **T̂** and **N̂**
- Unit vector (cross product of orthonormal vectors)
- Normal to osculating plane

**Formula** (direct from derivatives):
```
B̂ = (r' × r'') / |r' × r''|
```

**Geometric meaning**: "Third direction" completing orthonormal frame

### Frenet-Serret Frame

**The moving trihedron**: {**T̂**, **N̂**, **B**} at each point

**Properties**:
1. **Orthonormal**: T̂·N̂ = T̂·B = N̂·B = 0, |T̂| = |N̂| = |B| = 1
2. **Right-handed**: B = T̂ × N̂
3. **Moves with curve**: Frame rotates/translates as t varies

**Frenet-Serret Formulas** (evolution of frame):
```
dT̂/ds = κ N̂         (curvature × normal)
dN̂/ds = -κ T̂ + τ B   (curvature and torsion)
dB/ds = -τ N̂        (torsion × normal)
```

where s = arc length parameter

**Interpretation**:
- **κ**: Rate T̂ rotates in osculating plane
- **τ**: Rate osculating plane rotates around T̂

### Torsion

**Definition**: Rate at which osculating plane twists

**Formula**:
```
τ = (r' × r'') · r''' / |r' × r''|²
```

**Geometric meaning**:
- τ = 0: Plane curve (no twisting out of plane)
- τ ≠ 0: Space curve (3D)
- τ > 0: Right-hand twist
- τ < 0: Left-hand twist

**Examples**:
- Circle: τ = 0 (planar)
- Helix: τ = const (uniform twist)
- General plane curve: τ = 0 everywhere

### Osculating Planes

**Three fundamental planes** at each point:

**1. Osculating Plane**:
- **Contains**: **T̂** and **N̂**
- **Normal**: **B**
- **Meaning**: Plane of curvature, "best-fit" plane locally

**2. Normal Plane**:
- **Contains**: **N̂** and **B**
- **Normal**: **T̂**
- **Meaning**: Perpendicular to curve direction

**3. Rectifying Plane**:
- **Contains**: **T̂** and **B**
- **Normal**: **N̂**
- **Meaning**: Contains tangent and binormal

### Arc Length

**Definition**: Actual length along curve from t₁ to t₂

**Formula**:
```
s = ∫[t₁,t₂] |r'(t)| dt
```

**Arc length parameter**: Reparametrize by s instead of t
- **Property**: |dr/ds| = 1 (unit speed)
- **Simplifies formulas**: κ = |d²r/ds²|

**Check if arc-length parametrized**:
```cpp
bool isArcLen = curve.isArcLengthParametrized(t1, t2);
```

---

## Classes and Interfaces

### IParametricCurve<N> - Base Interface

**Purpose**: Generic N-dimensional parametric curve

```cpp
template<int N>
class IParametricCurve {
public:
    virtual VectorN<Real, N> operator()(Real t) const = 0;  // Curve evaluation
    virtual Real getMinT() const = 0;
    virtual Real getMaxT() const = 0;
};
```

**Usage**: Base for all parametric curves (2D, 3D, higher dimensions)

---

### ICurveCartesian3D - 3D Curve Interface

**Purpose**: Full differential geometry for 3D curves

**Inherits**: `IParametricCurve<3>`

**Methods**:

```cpp
class ICurveCartesian3D : public IParametricCurve<3> {
public:
    // Tangent vectors
    Vec3Cart getTangent(Real t) const;           // r'(t)
    Vec3Cart getTangentUnit(Real t) const;       // T̂ = r'/|r'|
    
    // Normal vectors  
    Vec3Cart getNormal(Real t) const;            // r''(t)
    Vec3Cart getNormalUnit(Real t) const;        // N̂ (principal normal)
    
    // Binormal
    Vec3Cart getBinormal(Real t) const;          // B = T̂ × N̂
    
    // Curvature
    Vec3Cart getCurvatureVector(Real t) const;   // κ_vec
    Real getCurvature(Real t) const;             // κ = |κ_vec|
    
    // Torsion
    Real getTorsion(Real t) const;               // τ
    
    // Planes
    Plane3D getOsculationPlane(Real t) const;    // Contains T̂, N̂
    Plane3D getNormalPlane(Real t) const;        // Contains N̂, B
    Plane3D getRectifyingPlane(Real t) const;    // Contains T̂, B
    
    // Complete frame
    void getMovingTrihedron(Real t, 
                            Vector3Cartesian& T,
                            Vector3Cartesian& N,
                            Vector3Cartesian& B);
    
    // Arc length check
    bool isArcLengthParametrized(Real t1, Real t2, int numPoints = 100) const;
};
```

**Implementation**: Uses numerical differentiation from `Derivation` module

---

### Specific Curve Implementations

All in `Curves` namespace:

#### 2D Curves (ICurveCartesian2D)

**Circle**:
```cpp
class Circle2DCurve : public ICurveCartesian2D {
    Real _radius;
    Pnt2Cart _center;
    
public:
    Circle2DCurve(Real r = 1.0, Pnt2Cart center = {0,0});
    VectorN<Real,2> operator()(Real t) const override {
        return {_center.X() + _radius*cos(t), 
                _center.Y() + _radius*sin(t)};
    }
};
```

**Log Spiral**: r(t) = e^(λt) (cos t, sin t), λ < 0

**Lemniscate**: Figure-eight curve

**Astroid**: Star-shaped curve with cusps

**Archimedean Spiral**: r(t) = at (cos t, sin t)

#### 3D Curves (ICurveCartesian3D)

**Helix**:
```cpp
class HelixCurve : public ICurveCartesian3D {
    Real _a, _b;  // Radius and pitch
    
public:
    HelixCurve(Real a = 1.0, Real b = 1.0) : _a(a), _b(b) {}
    
    VectorN<Real,3> operator()(Real t) const override {
        return {_a * cos(t), _a * sin(t), _b * t};
    }
    
    // For helix: κ = a/(a²+b²), τ = b/(a²+b²) (constant!)
};
```

**Toroidal Spiral**: Spirals on torus surface

**Viviani's Curve**: Intersection of sphere and cylinder

**Torus Knot**: Closed curves on torus (p,q parametrization)

---

### Predefined Curves (TestBeds)

**Available via**:
```cpp
namespace TestBeds {
    class ParametricCurvesTestBed {
        static const ParametricCurveData& getTestCurve(const string& name);
    };
}
```

**Curves available**:
- `"Helix"`: Standard helix (a=1, b=1)
- `"Helix2"`: Scaled helix (a=2, b=0.5)
- `"Circle3DXYCurve"`: Unit circle in xy-plane
- `"TwistedCubic"`: r(t) = (t, t², t³)
- `"TorusKnot_3_2"`: (3,2) torus knot
- `"Viviani"`: Viviani's curve
- `"Lissajous"`: Lissajous curve

**Each provides**:
- Curve object
- Analytical derivatives (for testing)
- Analytical curvature/torsion formulas (where available)

---

## Examples

### Example 1: Helix - Complete Analysis

Classic space curve with constant curvature and torsion:

```cpp
#include <mml/core/Curves.h>
#include "test_beds/parametric_curves_test_bed.h"

void Example1() {
    // Get helix from test bed: r(t) = (cos t, sin t, t)
    const auto& helixData = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix");
    const auto& helix = helixData._curve;
    
    Real t = Constants::PI / 4.0;  // t = 45°
    
    // 1. Tangent vector
    auto tangent = helix.getTangent(t);
    auto T = helix.getTangentUnit(t);
    
    std::cout << "Tangent:      " << tangent << "\n";
    std::cout << "Unit tangent: " << T << "\n";
    // Expected: T ≈ (-0.707, 0.707, 1) / √2
    
    // 2. Curvature
    Real kappa = helix.getCurvature(t);
    std::cout << "Curvature: " << kappa << "\n";
    // For standard helix (a=1, b=1): κ = 1/(1²+1²) = 0.5 (constant!)
    
    // 3. Principal normal
    auto N = helix.getNormalUnit(t);
    std::cout << "Principal normal: " << N << "\n";
    // Points toward z-axis (center of helix)
    
    // 4. Binormal
    auto B = helix.getBinormal(t);
    std::cout << "Binormal: " << B << "\n";
    // Perpendicular to osculating plane
    
    // 5. Torsion
    Real tau = helix.getTorsion(t);
    std::cout << "Torsion: " << tau << "\n";
    // For standard helix: τ = 1/(1²+1²) = 0.5 (constant!)
    
    // 6. Verify Frenet frame is orthonormal
    Real dot_TN = Utils::ScalarProduct<3>(T, N);
    Real dot_TB = Utils::ScalarProduct<3>(T, B);
    Real dot_NB = Utils::ScalarProduct<3>(N, B);
    
    std::cout << "\nOrthogonality check:\n";
    std::cout << "T·N = " << dot_TN << " (should be 0)\n";
    std::cout << "T·B = " << dot_TB << " (should be 0)\n";
    std::cout << "N·B = " << dot_NB << " (should be 0)\n";
    
    std::cout << "\nNorm check:\n";
    std::cout << "|T| = " << T.NormL2() << " (should be 1)\n";
    std::cout << "|N| = " << N.NormL2() << " (should be 1)\n";
    std::cout << "|B| = " << B.NormL2() << " (should be 1)\n";
}
```

### Example 2: Circle - Verify Curvature Formula

```cpp
void Example2() {
    // Circle of radius R has curvature κ = 1/R
    Real R = 2.5;
    Curves::Circle2DCurve circle(R);
    
    // Sample at multiple points
    std::vector<Real> test_points = {0, Constants::PI/4, Constants::PI/2, 
                                      Constants::PI, 3*Constants::PI/2};
    
    std::cout << "Circle (radius " << R << ") curvature verification:\n";
    for (Real t : test_points) {
        // For 2D, embed in 3D as (x, y, 0)
        ParametricCurve<3> circle3D([&](Real t_val) {
            auto p2d = circle(t_val);
            return VectorN<Real,3>{p2d[0], p2d[1], 0.0};
        });
        
        ICurveCartesian3D* curve3D = &circle3D;
        Real kappa = curve3D->getCurvature(t);
        
        std::cout << "t = " << t << ", κ = " << kappa 
                  << " (expected: " << 1.0/R << ")\n";
    }
    
    // All should give κ = 1/2.5 = 0.4
}
```

### Example 3: Frenet Frame Evolution

Visualize how frame rotates along curve:

```cpp
void Example3() {
    const auto& helix = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
    
    std::cout << "Frenet frame evolution along helix:\n";
    std::cout << "t\t\tT\t\t\tN\t\t\tB\n";
    
    for (Real t = 0; t <= 2*Constants::PI; t += Constants::PI/4) {
        Vector3Cartesian T, N, B;
        helix.getMovingTrihedron(t, T, N, B);
        
        std::cout << std::fixed << std::setprecision(2) << t << "\t";
        std::cout << "(" << T[0] << "," << T[1] << "," << T[2] << ")\t";
        std::cout << "(" << N[0] << "," << N[1] << "," << N[2] << ")\t";
        std::cout << "(" << B[0] << "," << B[1] << "," << B[2] << ")\n";
    }
    
    // Observe: T and N rotate in xy-plane, B stays constant direction!
}
```

### Example 4: Curvature Vector Visualization

```cpp
void Example4() {
    // Circle: curvature vector points toward center
    Curves::Circle2DCurve circle(1.0, {0, 0});
    
    ParametricCurve<3> circle3D([&](Real t) {
        auto p = circle(t);
        return VectorN<Real,3>{p[0], p[1], 0};
    });
    
    Real t = Constants::PI / 3.0;  // 60°
    
    auto position = circle3D(t);
    auto curvVec = static_cast<ICurveCartesian3D&>(circle3D).getCurvatureVector(t);
    
    std::cout << "Circle at t = " << t << ":\n";
    std::cout << "Position: " << position << "\n";
    std::cout << "Curvature vector: " << curvVec << "\n";
    
    // Magnitude should be κ = 1 (radius = 1)
    Real mag = curvVec.NormL2();
    std::cout << "Magnitude: " << mag << " (curvature)\n";
    
    // Direction: should point toward origin (center)
    // curvVec / mag = -position (for circle at origin)
    auto direction = curvVec / mag;
    std::cout << "Direction: " << direction << "\n";
    std::cout << "Expected (toward center): " << -position << "\n";
}
```

### Example 5: Osculating Planes

```cpp
void Example5() {
    const auto& helix = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
    
    Real t = Constants::PI / 2.0;
    
    // 1. Osculating plane (contains T and N)
    auto oscPlane = helix.getOsculationPlane(t);
    std::cout << "Osculating plane at t = " << t << ":\n";
    std::cout << "  Point: " << oscPlane.GetPoint() << "\n";
    std::cout << "  Normal: " << oscPlane.GetNormal() << " (= B)\n";
    
    // 2. Normal plane (perpendicular to curve)
    auto normPlane = helix.getNormalPlane(t);
    std::cout << "\nNormal plane:\n";
    std::cout << "  Normal: " << normPlane.GetNormal() << " (= T)\n";
    
    // 3. Rectifying plane
    auto rectPlane = helix.getRectifyingPlane(t);
    std::cout << "\nRectifying plane:\n";
    std::cout << "  Normal: " << rectPlane.GetNormal() << " (= N)\n";
    
    // Verify: Three planes mutually perpendicular
    auto B = oscPlane.GetNormal();
    auto T = normPlane.GetNormal();
    auto N = rectPlane.GetNormal();
    
    std::cout << "\nMutual perpendicularity:\n";
    std::cout << "B·T = " << Utils::ScalarProduct<3>(B, T) << "\n";
    std::cout << "B·N = " << Utils::ScalarProduct<3>(B, N) << "\n";
    std::cout << "T·N = " << Utils::ScalarProduct<3>(T, N) << "\n";
}
```

### Example 6: Plane vs Space Curve (Torsion Check)

```cpp
void Example6() {
    // Circle: planar curve (τ = 0)
    ParametricCurve<3> circle([](Real t) {
        return VectorN<Real,3>{cos(t), sin(t), 0};
    });
    
    // Helix: space curve (τ ≠ 0)
    const auto& helix = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
    
    Real t = Constants::PI / 4.0;
    
    Real tau_circle = static_cast<ICurveCartesian3D&>(circle).getTorsion(t);
    Real tau_helix = helix.getTorsion(t);
    
    std::cout << "Circle torsion: " << tau_circle << " (should be ~0)\n";
    std::cout << "Helix torsion:  " << tau_helix << " (should be 0.5)\n";
    
    if (std::abs(tau_circle) < 1e-6) {
        std::cout << "✓ Circle is planar (τ ≈ 0)\n";
    }
    if (std::abs(tau_helix - 0.5) < 1e-6) {
        std::cout << "✓ Helix is space curve with τ = 0.5\n";
    }
}
```

### Example 7: Arc Length Parametrization Check

```cpp
void Example7() {
    // Check if curve is arc-length parametrized
    
    // Helix: NOT arc-length parametrized
    const auto& helix = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
    bool isArcLen_helix = helix.isArcLengthParametrized(0.0, 2*Constants::PI);
    
    std::cout << "Helix arc-length parametrized: " 
              << (isArcLen_helix ? "YES" : "NO") << "\n";
    // Expected: NO (|r'(t)| = √2, not 1)
    
    // Create arc-length parametrized line
    ParametricCurve<3> line([](Real s) {
        return VectorN<Real,3>{s, 0, 0};  // |dr/ds| = 1
    });
    
    bool isArcLen_line = static_cast<ICurveCartesian3D&>(line)
        .isArcLengthParametrized(0.0, 10.0);
    
    std::cout << "Straight line arc-length parametrized: "
              << (isArcLen_line ? "YES" : "NO") << "\n";
    // Expected: YES
    
    // Actual arc length of helix from 0 to 2π
    Real arclen = PathIntegration::ParametricCurveLength(helix, 0.0, 2*Constants::PI);
    std::cout << "\nHelix arc length (one turn): " << arclen << "\n";
    std::cout << "Expected: " << 2*Constants::PI*std::sqrt(2) << "\n";
}
```

### Example 8: Custom Curve Definition

Create and analyze custom curve:

```cpp
void Example8() {
    // Define twisted cubic: r(t) = (t, t², t³)
    class TwistedCubic : public ICurveCartesian3D {
    public:
        VectorN<Real,3> operator()(Real t) const override {
            return {t, t*t, t*t*t};
        }
        Real getMinT() const override { return -2.0; }
        Real getMaxT() const override { return 2.0; }
    };
    
    TwistedCubic cubic;
    
    std::cout << "Twisted Cubic Analysis:\n";
    std::cout << "t\tκ\t\tτ\n";
    
    for (Real t = -1.0; t <= 1.0; t += 0.5) {
        if (std::abs(t) < 1e-10) continue;  // Skip t=0 (undefined)
        
        Real kappa = cubic.getCurvature(t);
        Real tau = cubic.getTorsion(t);
        
        std::cout << t << "\t" << kappa << "\t" << tau << "\n";
    }
    
    // Analytical formulas for twisted cubic:
    // κ(t) = 2√(1+4t²) / (1+4t²+9t⁴)^(3/2)
    // τ(t) = 6 / (1+4t²+9t⁴)
    
    // Note: Both vary with t (unlike helix which is constant)
}
```

---

## Integration with Other Modules

### With Derivation Module

**Curve methods use numerical differentiation**:

```cpp
// getTangent implementation (simplified):
Vec3Cart ICurveCartesian3D::getTangent(Real t) const {
    return Derivation::DeriveCurve<3>(*this, t, nullptr);
}

// getNormal implementation:
Vec3Cart ICurveCartesian3D::getNormal(Real t) const {
    return Derivation::DeriveCurveSec<3>(*this, t, nullptr);
}
```

**Advantages**:
- Works for any curve (no analytical formula needed)
- Automatic differentiation
- Configurable accuracy (via Derivation settings)

**Limitations**:
- Numerical errors (~1e-6 typical)
- Slower than analytical (if available)
- Can fail at singular points (cusps, inflections)

### With Integration Module

**Arc length computation**:

```cpp
#include <mml/core/Integration/PathIntegration.h>

Real length = PathIntegration::ParametricCurveLength(curve, t1, t2);
```

**Integrand**: |**r**'(t)|

### With TestBeds

**Analytical comparison** for testing:

```cpp
const auto& data = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix");

Real t = 0.5;
Real kappa_numerical = data._curve.getCurvature(t);
Real kappa_analytical = data._curvatureFunc(t);

Real error = std::abs(kappa_numerical - kappa_analytical);
std::cout << "Curvature error: " << error << "\n";
```

**Useful for**:
- Validating numerical methods
- Understanding precision limits
- Benchmarking

### With Visualizers

**Export curve for plotting**:

```cpp
// Generate curve points
std::ofstream out("helix.dat");
for (Real t = 0; t <= 2*Constants::PI; t += 0.01) {
    auto p = helix(t);
    out << p[0] << " " << p[1] << " " << p[2] << "\n";
}

// Plot with gnuplot, Python, etc.
```

**Visualize Frenet frame**:

```cpp
// At specific point, output T, N, B as arrows
Real t = Constants::PI / 2;
auto pos = helix(t);
auto T = helix.getTangentUnit(t) * 0.5;  // Scale for visibility
auto N = helix.getNormalUnit(t) * 0.5;
auto B = helix.getBinormal(t) * 0.5;

// Write arrow data for visualization tool
```

---

## Best Practices

### Numerical Stability

**Issues at special points**:

```cpp
// DON'T: Evaluate at singularity
ParametricCurve<3> cubic([](Real t) {
    return VectorN<Real,3>{t*t, t*t, t*t*t};  // Cusp at t=0!
});
Real kappa = static_cast<ICurveCartesian3D&>(cubic).getCurvature(0.0);
// Will give large errors or NaN
```

**DO**: Avoid or handle singularities

```cpp
// Check for issues before computing
Real epsilon = 1e-8;
Real t_safe = (std::abs(t) < epsilon) ? epsilon : t;
Real kappa = curve.getCurvature(t_safe);
```

### Verification

**Always verify Frenet frame properties**:

```cpp
auto T = curve.getTangentUnit(t);
auto N = curve.getNormalUnit(t);
auto B = curve.getBinormal(t);

// 1. Unit vectors
assert(std::abs(T.NormL2() - 1.0) < 1e-6);
assert(std::abs(N.NormL2() - 1.0) < 1e-6);
assert(std::abs(B.NormL2() - 1.0) < 1e-6);

// 2. Orthogonality
assert(std::abs(Utils::ScalarProduct<3>(T, N)) < 1e-6);
assert(std::abs(Utils::ScalarProduct<3>(T, B)) < 1e-6);
assert(std::abs(Utils::ScalarProduct<3>(N, B)) < 1e-6);

// 3. Right-handed (B = T × N)
auto B_check = VectorProduct(T, N);
assert((B - B_check).NormL2() < 1e-6);
```

### Common Pitfalls

❌ **Pitfall 1**: Assuming all curves are smooth

```cpp
// Cusps, kinks cause numerical issues
// Example: Astroid has cusps at t = 0, π/2, π, 3π/2
Curves::AstroidCurve astroid;
Real kappa = astroid.getCurvature(0);  // BAD: cusp!
```

❌ **Pitfall 2**: Confusing normal vectors

```cpp
// getNormal() returns r''(t) (acceleration)
// getNormalUnit() returns N̂ (principal normal, unit vector)

auto acc = curve.getNormal(t);        // Not necessarily unit!
auto N = curve.getNormalUnit(t);      // Unit vector
```

❌ **Pitfall 3**: 2D/3D confusion

```cpp
// 2D curve methods different from 3D
ICurveCartesian2D* curve2D = &circle2d;
// curve2D->getBinormal(t);  // ERROR: No binormal in 2D!
```

✅ **Best Practice Pattern**:

```cpp
// 1. Define curve clearly
class MyCurve : public ICurveCartesian3D {
    VectorN<Real,3> operator()(Real t) const override {
        // ...
    }
};

// 2. Check for valid parameter range
Real t = ...;
assert(t >= curve.getMinT() && t <= curve.getMaxT());

// 3. Compute geometry
Real kappa = curve.getCurvature(t);
Real tau = curve.getTorsion(t);

// 4. Verify sanity
assert(!std::isnan(kappa) && !std::isinf(kappa));
assert(kappa >= 0);  // Curvature always non-negative

// 5. Use in application
if (kappa > 1.0) {
    std::cout << "Sharp turn at t = " << t << "\n";
}
```

---

## Advanced Topics

### Fundamental Theorem of Curves

**Theorem**: Curvature κ(s) and torsion τ(s) (as functions of arc length s) **completely determine** curve shape (up to rigid motion).

**Implication**: Two curves with same κ(s) and τ(s) differ only by translation/rotation.

**Example**: All helices with same κ and τ are congruent.

### Natural Equations

**Express curve via** κ(s) and τ(s) instead of **r**(t)

**Example - Circle**:
```
κ(s) = 1/R (constant)
τ(s) = 0   (planar)
```

**Example - Helix**:
```
κ(s) = a/(a²+b²) (constant)
τ(s) = b/(a²+b²) (constant)
```

### Evolute and Involute

**Evolute**: Locus of centers of curvature

**Position**: **r**(t) + (1/κ)**N**

**Involute**: Reverse process (string unwinding from curve)

**Not currently implemented** (future extension)

### Bertrand Curves

**Pair of curves** sharing principal normals

**Condition**: Linear relation between curvatures and torsions
```
aκ₁ + bτ₁ = constant
```

**Example**: Circular helix has Bertrand mate

### Surfaces And Intrinsic Geometry

**First Fundamental Form**: I = E du² + 2F du dv + G dv²
- Measures distances on surface
- Intrinsic geometry; exposed as `InducedMetric2D : MetricTensorField<2>`

**Second Fundamental Form**: II = L du² + 2M du dv + N dv²
- Measures curvature
- How surface bends in space

**Gaussian Curvature**: K = (LN - M²)/(EG - F²)
- Intrinsic (same for isometric surfaces)
- K > 0: Elliptic (sphere-like)
- K = 0: Parabolic (cylinder-like)
- K < 0: Hyperbolic (saddle-like)

**Mean Curvature**: H = (EN + GL - 2FM)/(2(EG - F²))
- Average of principal curvatures

**Intrinsic Gaussian Curvature**:

```cpp
Sphere sphere(1.0);
Real K = DifferentialGeometry::GaussianCurvatureIntrinsic(
    sphere, Constants::PI / 3.0, Constants::PI / 4.0);
```

The implementation computes intrinsic curvature through the induced metric and
`MetricTensorField<2>` curvature pipeline. Theorema Egregium tests compare it to
the existing extrinsic `surface.GaussianCurvature(u, w)` path.

**Surface Geodesics**:

```cpp
auto sol = DifferentialGeometry::IntegrateSurfaceGeodesicFixedStep(
    surface,
    VectorN<Real, 2>{u0, w0},
    VectorN<Real, 2>{du0, dw0},
    0.0,
    1.0,
    100);
```

This reuses `GeodesicEquationSystem<2>` with the surface's induced metric.

### Differential Forms And Stokes-Family Checks

Typed one-forms and two-forms can now be integrated directly:

```cpp
class RotationOneForm2D : public IFormField<2, 1, Cartesian2> {
public:
    Form1<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& p) const override {
        Form1<2, Cartesian2> omega;
        omega.Component(0) = -0.5 * p[1];
        omega.Component(1) =  0.5 * p[0];
        return omega;
    }
};

RotationOneForm2D omega;
Real residual = DifferentialGeometry::GreenTheoremResidual(
    omega, 0.0, 1.0, 0.0, 1.0);
```

Core helpers:

- `IntegrateOneForm(omega, curve, t0, t1)` computes `int_C omega`.
- `IntegrateTwoForm(omega, surface, u0, u1, w0, w1)` computes `int_S omega`.
- `GreenTheoremResidual(omega, x0, x1, y0, y1)` checks `int_boundary omega - int_region d omega`.
- `StokesResidual(omega, surface)` checks `int_boundary omega - int_surface d omega` for rectangular surface patches.

Orientation convention for rectangular domains is counter-clockwise in parameter
space: bottom edge, right edge, top edge reversed, left edge reversed. For a
surface chart `r(u,w)`, the corresponding surface orientation follows
`partial_u r x partial_w r`.

### Charts And Frames

`FunctionChart<DomainN, AmbientN, DomainFrame, AmbientFrame>` adapts any
embedding `R^DomainN -> R^AmbientN` into a single coordinate patch:

```cpp
Sphere sphere(1.0);
DifferentialGeometry::FunctionChart<2, 3, Parameter2, Cartesian3> chart(sphere);
Point<Real, 2, Parameter2> q{Constants::PI / 3.0, Constants::PI / 4.0};

auto p = chart.map_point(q);
auto J = chart.jacobian(q);
auto g = chart.induced_metric(q);
```

Chart helpers support tangent push-forward and covector/form pullback. Atlas and
transition-map machinery starts deliberately small: `Atlas` stores homogeneous
charts with the same domain frame, explicit transition functions, and optional
validity predicates.

```cpp
auto atlas = DifferentialGeometry::MakeUnitSphereStereographicAtlas();
Point<Real, 2, Parameter2> north{1.5, 0.5};
auto south = atlas.transition(0, 1, north);
auto pNorth = atlas.chart(0).map_point(north);
auto pSouth = atlas.chart(1).map_point(south);
```

The built-in unit-sphere atlas uses north and south stereographic charts. It is
intended as a compact manifold testbed, not as full automatic atlas gluing.
Transition maps are registered explicitly and throw when the point is outside
the overlap, such as trying to transition a pole coordinate.

Frame helpers:

- `ComputeFrenetFrame<N>(curve, t)` gives tangent, principal normal, speed, and curvature for arbitrary-dimensional curves.
- `ComputeDarbouxFrame(surface, parameterCurve, t)` gives tangent, surface normal, tangent-normal, curvature vector, normal curvature, and geodesic curvature for curves on a surface.

For a great-circle equator on the unit sphere, `geodesicCurvature` is zero and
`abs(normalCurvature)` is one.

### Numerical And Singularity Guidance

- Avoid chart singularities in finite-difference tests, such as sphere poles or cylindrical axes.
- Prefer invariant checks (`K_intrinsic == K_extrinsic`, Green/Stokes residuals near zero) over isolated component checks.
- Use enough subdivisions for form integration; default midpoint surface integration is intended for smooth patches.
- Coordinate chart orientation matters. Reversing boundary order or swapping `(u,w)` changes signs.

### Focused Validation Tags

Use these focused tags while working on this subsystem:

```powershell
& .\build\tests\Release\MML_Tests.exe '[DifferentialGeometry],[DifferentialForms],[MetricTensor],[Geodesic],[TypedFields],[CoordinateMap]'
```

Most recent full quality gate for this roadmap:

```powershell
& .\build\tests\Release\MML_Tests.exe
```

---

## Summary

### Key Formulas

| Quantity | Formula |
|----------|---------|
| **Tangent** | **T̂** = **r**'/\|**r**'\| |
| **Curvature** | κ = \|**r**' × **r**''\|/\|**r**'\|³ |
| **Principal Normal** | **N̂** = (**r**' × (**r**'' × **r**'))/(norm) |
| **Binormal** | **B** = **r**' × **r**''/\|**r**' × **r**''\| |
| **Torsion** | τ = (**r**'×**r**'')·**r**'''/\|**r**'×**r**''\|² |

### Geometric Insights

1. ✅ **κ measures bending** (deviation from straight line)
2. ✅ **τ measures twisting** (out-of-plane rotation)
3. ✅ **Frenet frame** {**T**, **N**, **B**} is orthonormal
4. ✅ **Plane curves** have τ = 0
5. ✅ **κ and τ determine curve** (fundamental theorem)
6. ✅ **Circle**: κ = 1/R, τ = 0
7. ✅ **Helix**: Both κ and τ constant
8. ✅ **Osculating plane** contains **T** and **N**

### When to Use

**Curve analysis when**:
- Trajectory planning (robotics, animation)
- Road/track design (civil engineering)
- Physics (charged particle in magnetic field)
- Computer graphics (camera paths)
- Differential geometry research

**Method choice**:
- Analytical formulas: When available, more accurate
- Numerical (this library): General, automatic, good for complex curves

### References

- **Differential Geometry of Curves and Surfaces** (do Carmo)
- **Elementary Differential Geometry** (Pressley)
- **Numerical Recipes**: Chapter on derivatives
- **Curves and Surfaces for CAGD** (Farin)
- **Classical Differential Geometry** (Struik)

---

## Runnable Examples

Working code examples are available in `src/docs_demos/docs_demo_diff_geometry.cpp`:

### 2D Curves
- `Docs_Demo_Curves_2D_Basic()` - Circle, spirals, special curves

### 3D Curves
- `Docs_Demo_Curves_3D_Basic()` - Lines, circles, helix, toroidal spiral

### Frenet Frame
- `Docs_Demo_Curves_Frenet_Frame()` - Tangent, normal, binormal computation

### Curvature & Torsion
- `Docs_Demo_Curves_Curvature_Torsion()` - Measuring bending and twisting

### Arc Length
- `Docs_Demo_Curves_Arc_Length()` - Arc length calculation and parameterization

### Osculating & Normal Planes
- `Docs_Demo_Curves_Planes()` - Computing planes associated with curves

### Surfaces
- `Docs_Demo_Surfaces_Basic()` - Predefined surfaces (sphere, torus, cylinder, etc.)
- `Docs_Demo_Surfaces_Curvature()` - Gaussian and mean curvature
- `Docs_Demo_Surfaces_Principal_Curvatures()` - Principal curvatures and directions
- `Docs_Demo_Surfaces_Fundamental_Forms()` - First and second fundamental forms
- `Docs_Demo_Surfaces_Special()` - Special surfaces (Möbius strip, monkey saddle)
- `Docs_Demo_Surfaces_Tangents_Normals()` - Surface tangent and normal vectors
- `Docs_Demo_Surfaces_Classification()` - Elliptic, hyperbolic, parabolic points

---

**Part of MinimalMathLibrary** - `mml/core/Curves.h`, `mml/core/Surfaces.h`
