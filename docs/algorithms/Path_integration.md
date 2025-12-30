# Path Integration

Comprehensive toolkit for **line integrals, curve lengths, and work integrals** along parametric curves in multidimensional space.

## Overview

**Covers**:
- **Curve Length**: Arc length of parametric curves
- **Curve Mass**: Weighted integration with variable density
- **Scalar Line Integrals**: ∫_C f ds (scalar field along curve)
- **Vector Line Integrals**: ∫_C **F**·d**r** (work, circulation)

**Implementation**: Static methods in `PathIntegration` class using adaptive trapezoidal integration.

## Quick Reference

### Integration Types

| Type | Formula | Physical Meaning | Method |
|------|---------|------------------|--------|
| **Arc Length** | ∫ \|**r**'(t)\| dt | Total distance along curve | `ParametricCurveLength` |
| **Curve Mass** | ∫ ρ(t) \|**r**'(t)\| dt | Mass with variable density | `ParametricCurveMass` |
| **Scalar Line Integral** | ∫_C f ds | Total of scalar field along path | `LineIntegral(scalarField, ...)` |
| **Vector Line Integral** | ∫_C **F**·d**r** | Work done by force field | `LineIntegral(vectorField, ...)` |

### Common Applications

| Application | Integral Type | Example |
|-------------|---------------|---------|
| **Wire length** | Arc length | Cable, rope, track |
| **Wire mass** | Curve mass | Non-uniform density wire |
| **Heat along wire** | Scalar line integral | Temperature distribution |
| **Work by force** | Vector line integral | Moving particle in field |
| **Electric potential** | Scalar line integral | Charge distribution |
| **Fluid circulation** | Vector line integral | Fluid flow around loop |

---

## Mathematical Background

### Arc Length

**For parametric curve** **r**(t) = (x(t), y(t), z(t)), t ∈ [a, b]

**Formula**:
```
L = ∫[a,b] |r'(t)| dt = ∫[a,b] √((dx/dt)² + (dy/dt)² + (dz/dt)²) dt
```

**Geometric meaning**: Actual distance traveled along curve

**Examples**:
```
Circle (R=1): r(t) = (cos t, sin t, 0), t ∈ [0, 2π]
    |r'(t)| = √(sin²t + cos²t) = 1
    L = ∫[0,2π] 1 dt = 2π

Helix: r(t) = (cos t, sin t, t), t ∈ [0, 2π]
    |r'(t)| = √(sin²t + cos²t + 1) = √2
    L = ∫[0,2π] √2 dt = 2π√2
```

### Scalar Line Integral

**Definition**: Integrate scalar function f along curve C

**Formula**:
```
∫_C f ds = ∫[a,b] f(r(t)) |r'(t)| dt
```

where:
- f: Scalar field (scalar function of position)
- ds = |r'(t)| dt: Arc length element
- **r**(t): Parametric curve

**Physical interpretation**:
- If f = temperature, integral = total heat along wire
- If f = density, integral = mass
- If f = 1, integral = arc length

**Example** (linear field on circle):
```
f(x,y,z) = x + y
Curve: r(t) = (3 cos t, 3 sin t, 0), t ∈ [0, π/2]
    f(r(t)) = 3 cos t + 3 sin t
    |r'(t)| = 3
    ∫ = ∫[0,π/2] (3 cos t + 3 sin t) · 3 dt
      = 9 ∫[0,π/2] (cos t + sin t) dt
      = 9 [sin t - cos t]₀^(π/2)
      = 9 [(1 - 0) - (0 - 1)] = 18
```

### Vector Line Integral (Work Integral)

**Definition**: Integrate dot product of vector field with curve tangent

**Formula**:
```
∫_C F·dr = ∫[a,b] F(r(t)) · r'(t) dt
```

where:
- **F**: Vector field (force, velocity, etc.)
- d**r** = **r**'(t) dt: Differential displacement
- **F**·d**r**: Component of field along curve

**Physical interpretation**:
- **F** = force: Work done moving along curve
- **F** = velocity: Circulation of fluid
- **F** = electric field: Voltage drop

**Conservative fields**:
```
If F = ∇φ (gradient of potential φ), then:
    ∫_C F·dr = φ(r(b)) - φ(r(a))  (path-independent!)
    
    ∮_C F·dr = 0  (closed loop)
```

**Example** (inverse-square force):
```
F(r) = -k r/|r|³  (gravitational/electric force)
φ(r) = -k/|r|     (potential)

∫_C F·dr = φ(end) - φ(start)  (independent of path!)
```

### Curve Mass with Variable Density

**For wire with linear density** ρ(t) along curve **r**(t)

**Formula**:
```
M = ∫[a,b] ρ(t) |r'(t)| dt
```

**Physical meaning**: Total mass when density varies along curve

**Example** (exponentially varying density):
```
Helix: r(t) = (cos t, sin t, t), t ∈ [0, 2π]
Density: ρ(t) = e^(-t)  (decreases with height)

M = ∫[0,2π] e^(-t) √2 dt = √2 [1 - e^(-2π)]
```

---

## PathIntegration Class

### Class Structure

```cpp
namespace MML {
    class PathIntegration {
    public:
        // Arc length
        template<int N>
        static Real ParametricCurveLength(
            const IParametricCurve<N>& curve,
            Real a, Real b
        );
        
        // Mass with variable density
        template<int N>
        static Real ParametricCurveMass(
            const IParametricCurve<N>& curve,
            const IRealFunction& density,
            Real a, Real b
        );
        
        // Scalar line integral: ∫_C f ds
        static Real LineIntegral(
            const IScalarFunction<3>& scalarField,
            const IParametricCurve<3>& curve,
            Real t1, Real t2,
            Real eps = Defaults::WorkIntegralPrecision
        );
        
        // Vector line integral (work): ∫_C F·dr
        static Real LineIntegral(
            const IVectorFunction<3>& vectorField,
            const IParametricCurve<3>& curve,
            Real t1, Real t2,
            Real eps = Defaults::LineIntegralPrecision
        );
    };
}
```

**All methods are static** - no instance needed.

---

### ParametricCurveLength - Arc Length

**Signature**:
```cpp
template<int N>
static Real ParametricCurveLength(
    const IParametricCurve<N>& curve,
    Real a,
    Real b
);
```

**Parameters**:
- `curve`: Parametric curve **r**(t)
- `a`, `b`: Parameter range [a, b]

**Returns**: Arc length L = ∫[a,b] |**r**'(t)| dt

**Implementation**:
```cpp
// Helper function for integrand
class HelperCurveLen : public IRealFunction {
    const IParametricCurve<N>& _curve;
public:
    Real operator()(Real t) const {
        auto tangent = Derivation::DeriveCurve<N>(_curve, t, nullptr);
        return tangent.NormL2();  // |r'(t)|
    }
};

// Integration via adaptive trapezoidal rule
return IntegrateTrap(helper, a, b);
```

**Usage**:
```cpp
#include "core/Curves.h"
#include "core/Integration/PathIntegration.h"

// Circle of radius R
ParametricCurve<3> circle([](Real t) {
    return VectorN<Real,3>{cos(t), sin(t), 0};
});

Real length = PathIntegration::ParametricCurveLength(circle, 0, 2*Constants::PI);
// Result: 2π (circumference)
```

**Works for any dimension**: N = 2 (plane), N = 3 (space), N > 3 (higher dimensions)

---

### ParametricCurveMass - Weighted Arc Length

**Signature**:
```cpp
template<int N>
static Real ParametricCurveMass(
    const IParametricCurve<N>& curve,
    const IRealFunction& density,
    Real a,
    Real b
);
```

**Parameters**:
- `curve`: Parametric curve **r**(t)
- `density`: Linear density function ρ(t)
- `a`, `b`: Parameter range

**Returns**: M = ∫[a,b] ρ(t) |**r**'(t)| dt

**Implementation**:
```cpp
class HelperCurveMass : public IRealFunction {
    const IParametricCurve<N>& _curve;
    const IRealFunction& _density;
public:
    Real operator()(Real t) const {
        auto tangent = Derivation::DeriveCurve<N>(_curve, t, nullptr);
        return tangent.NormL2() * _density(t);
    }
};
```

**Usage**:
```cpp
// Helix with exponential density
ParametricCurve<3> helix([](Real t) {
    return VectorN<Real,3>{cos(t), sin(t), t};
});

RealFunction density([](Real t) {
    return exp(-t);  // Decreases with height
});

Real mass = PathIntegration::ParametricCurveMass(helix, density, 0, 2*Constants::PI);
```

---

### LineIntegral (Scalar) - Scalar Field Integration

**Signature**:
```cpp
static Real LineIntegral(
    const IScalarFunction<3>& scalarField,
    const IParametricCurve<3>& curve,
    Real t1,
    Real t2,
    Real eps = Defaults::WorkIntegralPrecision
);
```

**Parameters**:
- `scalarField`: Scalar function f(**r**)
- `curve`: Parametric curve **r**(t)
- `t1`, `t2`: Parameter range
- `eps`: Integration precision (optional)

**Returns**: ∫_C f ds = ∫[t1,t2] f(**r**(t)) |**r**'(t)| dt

**Implementation**:
```cpp
class HelperLineIntegralScalarFunc : public IRealFunction {
    const IScalarFunction<N>& _scalar_field;
    const IParametricCurve<N>& _curve;
public:
    Real operator()(Real t) const {
        auto tangent = Derivation::DeriveCurve<N>(_curve, t, nullptr);
        auto field_val = _scalar_field(_curve(t));
        return field_val * tangent.NormL2();
    }
};
```

**Usage**:
```cpp
// Temperature field
ScalarFunction<3> temperature([](const VectorN<Real,3>& r) {
    return r[0] + r[1];  // T = x + y
});

// Half circle (radius 3)
ParametricCurve<3> path([](Real t) {
    return VectorN<Real,3>{3*cos(t), 3*sin(t), 0};
});

Real integral = PathIntegration::LineIntegral(
    temperature, path, 0, Constants::PI/2, 1e-3
);
// Result: 18.0 (analytical)
```

---

### LineIntegral (Vector) - Work/Circulation

**Signature**:
```cpp
static Real LineIntegral(
    const IVectorFunction<3>& vectorField,
    const IParametricCurve<3>& curve,
    Real t1,
    Real t2,
    Real eps = Defaults::LineIntegralPrecision
);
```

**Parameters**:
- `vectorField`: Vector field **F**(**r**)
- `curve`: Parametric curve **r**(t)
- `t1`, `t2`: Parameter range
- `eps`: Integration precision

**Returns**: ∫_C **F**·d**r** = ∫[t1,t2] **F**(**r**(t)) · **r**'(t) dt

**Implementation**:
```cpp
class HelperLineIntegralVectorFunc : public IRealFunction {
    const IVectorFunction<N>& _vector_field;
    const IParametricCurve<N>& _curve;
public:
    Real operator()(Real t) const {
        auto tangent = Derivation::DeriveCurve<N>(_curve, t, nullptr);
        auto field = _vector_field(_curve(t));
        return Utils::ScalarProduct<N>(tangent, field);  // F·r'
    }
};
```

**Usage**:
```cpp
// Force field (inverse-square)
VectorFunction<3> force([](const VectorN<Real,3>& r) {
    Real r_mag = r.NormL2();
    return r * (10.0 / (r_mag * r_mag * r_mag));
});

// Straight line path
Curves::LineCurve line(
    0.0, Point3Cartesian{1, 0, 0},
    1.0, Point3Cartesian{2, 1, 1}
);

Real work = PathIntegration::LineIntegral(force, line, 0.0, 1.0);
// Work done by force moving along line
```

---

## Examples

### Example 1: Circle and Helix Lengths

Classic curves with known analytical formulas:

```cpp
#include "core/Integration/PathIntegration.h"

void Example1() {
    // Circle: r(t) = (cos t, sin t, 0)
    ParametricCurve<3> circle([](Real t) {
        return VectorN<Real,3>{cos(t), sin(t), 0};
    });
    
    Real len_circle = PathIntegration::ParametricCurveLength(
        circle, 0, 2*Constants::PI
    );
    
    std::cout << "Circle length: " << len_circle << "\n";
    std::cout << "Expected:      " << 2*Constants::PI << "\n";
    // Output: ~6.283
    
    // Helix: r(t) = (cos t, sin t, t)
    ParametricCurve<3> helix([](Real t) {
        return VectorN<Real,3>{cos(t), sin(t), t};
    });
    
    Real len_helix = PathIntegration::ParametricCurveLength(
        helix, 0, 2*Constants::PI
    );
    
    std::cout << "Helix length: " << len_helix << "\n";
    std::cout << "Expected:     " << 2*Constants::PI*sqrt(2) << "\n";
    // Output: ~8.886
    
    // Analytical formulas:
    // Circle: |r'| = 1 → L = 2π
    // Helix: |r'| = √2 → L = 2π√2
}
```

### Example 2: Scalar Line Integral - Temperature Distribution

```cpp
void Example2() {
    // Linear temperature field: T(x,y,z) = x + y
    ScalarFunction<3> temperature([](const VectorN<Real,3>& r) {
        return r[0] + r[1];
    });
    
    // Path: quarter circle (radius 3)
    ParametricCurve<3> quarter_circle([](Real t) {
        return VectorN<Real,3>{3*cos(t), 3*sin(t), 0};
    });
    
    Real integral = PathIntegration::LineIntegral(
        temperature, quarter_circle, 
        0, Constants::PI/2,
        1e-4  // Precision
    );
    
    std::cout << "Temperature integral: " << integral << "\n";
    // Expected: 18.0
    
    // Analytical calculation:
    // f(r(t)) = 3 cos t + 3 sin t
    // |r'(t)| = 3
    // ∫[0,π/2] (3 cos t + 3 sin t) · 3 dt
    // = 9 ∫[0,π/2] (cos t + sin t) dt
    // = 9 [sin t - cos t]₀^(π/2)
    // = 9 [(1-0) - (0-1)] = 18
}
```

### Example 3: Work by Conservative Force Field

```cpp
void Example3() {
    // Potential: φ(r) = -k/|r|
    Real k = 10.0;
    ScalarFunction<3> potential([k](const VectorN<Real,3>& r) {
        return -k / r.NormL2();
    });
    
    // Force field: F = -∇φ = -k r/|r|³
    VectorFunction<3> force([k](const VectorN<Real,3>& r) {
        Real r_mag = r.NormL2();
        return r * (k / (r_mag * r_mag * r_mag));
    });
    
    // Two points
    Point3Cartesian p1{3, -5, 2};
    Point3Cartesian p2{5, 5, 2};
    
    // Straight line path
    Curves::LineCurve straight_line(
        0.0, p1,
        1.0, p2
    );
    
    // Work integral
    Real work = PathIntegration::LineIntegral(
        force, straight_line, 0.0, 1.0, 1e-3
    );
    
    // Verify with potential difference
    VectorN<Real,3> v1{p1.X(), p1.Y(), p1.Z()};
    VectorN<Real,3> v2{p2.X(), p2.Y(), p2.Z()};
    
    Real phi_diff = potential(v2) - potential(v1);
    
    std::cout << "Work (path integral):      " << work << "\n";
    std::cout << "Potential difference:      " << phi_diff << "\n";
    std::cout << "Difference (should be ~0): " << work - phi_diff << "\n";
    
    // For conservative field: Work = φ(end) - φ(start) (path-independent!)
}
```

### Example 4: Path Independence Verification

```cpp
void Example4() {
    // Conservative potential field
    ScalarFunction<3> potential([](const VectorN<Real,3>& r) {
        return -10.0 / r.NormL2();
    });
    
    // Force field (gradient)
    VectorFunction<3> force([](const VectorN<Real,3>& r) {
        Real r_mag = r.NormL2();
        return r * (10.0 / (r_mag * r_mag * r_mag));
    });
    
    Point3Cartesian p1{3, -5, 2};
    Point3Cartesian p2{5, 5, 2};
    Point3Cartesian p_mid{4, 0, 2};
    
    // Path 1: Direct line
    Curves::LineCurve path1(0.0, p1, 1.0, p2);
    Real work1 = PathIntegration::LineIntegral(force, path1, 0.0, 1.0);
    
    // Path 2: Two-segment path through midpoint
    Curves::LineCurve path2a(0.0, p1, 1.0, p_mid);
    Curves::LineCurve path2b(0.0, p_mid, 1.0, p2);
    
    Real work2a = PathIntegration::LineIntegral(force, path2a, 0.0, 1.0);
    Real work2b = PathIntegration::LineIntegral(force, path2b, 0.0, 1.0);
    Real work2_total = work2a + work2b;
    
    std::cout << "Work (direct path):    " << work1 << "\n";
    std::cout << "Work (two segments):   " << work2_total << "\n";
    std::cout << "Difference:            " << std::abs(work1 - work2_total) << "\n";
    
    // Should be equal (conservative field is path-independent)
}
```

### Example 5: Helix with Variable Density

```cpp
void Example5() {
    // Helix: r(t) = (a cos t, a sin t, b t)
    Real a = 1.0, b = 1.0;
    ParametricCurve<3> helix([a, b](Real t) {
        return VectorN<Real,3>{a*cos(t), a*sin(t), b*t};
    });
    
    // Exponentially decaying density (lighter at top)
    RealFunction density([](Real t) {
        return exp(-0.1 * t);
    });
    
    Real mass = PathIntegration::ParametricCurveMass(
        helix, density, 0, 2*Constants::PI
    );
    
    std::cout << "Helix mass (varying density): " << mass << "\n";
    
    // Compare to uniform density (ρ=1)
    RealFunction unit_density([](Real t) { return 1.0; });
    Real mass_uniform = PathIntegration::ParametricCurveMass(
        helix, unit_density, 0, 2*Constants::PI
    );
    
    std::cout << "Helix mass (uniform):         " << mass_uniform << "\n";
    std::cout << "Ratio:                        " << mass / mass_uniform << "\n";
    
    // Uniform should equal arc length: 2π√(a²+b²) = 2π√2
}
```

### Example 6: 2D Curve in Plane

```cpp
void Example6() {
    // 2D parabola: y = x², x ∈ [0, 2]
    // Parametrize: r(t) = (t, t²)
    ParametricCurve<2> parabola([](Real t) {
        return VectorN<Real,2>{t, t*t};
    });
    
    Real length = PathIntegration::ParametricCurveLength(parabola, 0, 2);
    
    std::cout << "Parabola arc length (0 to 2): " << length << "\n";
    
    // Analytical formula:
    // |r'(t)| = √(1 + 4t²)
    // ∫[0,2] √(1 + 4t²) dt ≈ 4.647
    
    // With scalar field: f(x,y) = y (height)
    ScalarFunction<2> height([](const VectorN<Real,2>& r) {
        return r[1];  // y-coordinate
    });
    
    // Need 3D embedding for LineIntegral (current implementation)
    ParametricCurve<3> parabola3D([](Real t) {
        return VectorN<Real,3>{t, t*t, 0};
    });
    
    ScalarFunction<3> height3D([](const VectorN<Real,3>& r) {
        return r[1];
    });
    
    Real integral = PathIntegration::LineIntegral(
        height3D, parabola3D, 0, 2
    );
    
    std::cout << "Height integral:               " << integral << "\n";
    // ∫ y ds along parabola
}
```

### Example 7: Closed Loop Circulation

```cpp
void Example7() {
    // Vector field (rotation): F = (-y, x, 0)
    VectorFunction<3> rotation_field([](const VectorN<Real,3>& r) {
        return VectorN<Real,3>{-r[1], r[0], 0};
    });
    
    // Circle (radius R, counterclockwise)
    Real R = 2.0;
    ParametricCurve<3> circle([R](Real t) {
        return VectorN<Real,3>{R*cos(t), R*sin(t), 0};
    });
    
    Real circulation = PathIntegration::LineIntegral(
        rotation_field, circle, 0, 2*Constants::PI
    );
    
    std::cout << "Circulation around circle: " << circulation << "\n";
    
    // Analytical calculation:
    // r(t) = (R cos t, R sin t, 0)
    // r'(t) = (-R sin t, R cos t, 0)
    // F(r(t)) = (-R sin t, R cos t, 0)
    // F·r' = R² sin²t + R² cos²t = R²
    // ∫[0,2π] R² dt = 2πR²
    
    std::cout << "Expected (2πR²):           " << 2*Constants::PI*R*R << "\n";
    
    // Non-zero circulation indicates rotational (non-conservative) field
}
```

### Example 8: Custom Curve Class Usage

```cpp
void Example8() {
    // Using predefined LineCurve
    Point3Cartesian start{0, 0, 0};
    Point3Cartesian end{1, 1, 1};
    
    Curves::LineCurve line(
        0.0, start,  // t = 0 → (0,0,0)
        1.0, end     // t = 1 → (1,1,1)
    );
    
    // Line length (should be √3)
    Real length = PathIntegration::ParametricCurveLength(line, 0, 1);
    std::cout << "Line length: " << length << " (expected: " << sqrt(3) << ")\n";
    
    // Constant field along line
    ScalarFunction<3> constant([](const VectorN<Real,3>& r) {
        return 5.0;
    });
    
    Real integral = PathIntegration::LineIntegral(constant, line, 0, 1);
    std::cout << "Constant field integral: " << integral << "\n";
    std::cout << "Expected (5 * √3):       " << 5 * sqrt(3) << "\n";
    
    // For constant field: ∫ c ds = c * length
}
```

---

## Integration with Other Modules

### With Curves Module

**PathIntegration designed to work with** `IParametricCurve<N>`:

```cpp
#include "core/Curves.h"
#include "core/Integration/PathIntegration.h"

// Any curve implementing IParametricCurve<N>
Curves::Circle2DCurve circle(1.0);  // Radius 1

// Embed in 3D for integration
ParametricCurve<3> circle3D([&circle](Real t) {
    auto p2d = circle(t);
    return VectorN<Real,3>{p2d[0], p2d[1], 0};
});

Real length = PathIntegration::ParametricCurveLength(circle3D, 0, 2*Constants::PI);
```

**Predefined curves** from `Curves` namespace work directly.

### With Derivation Module

**PathIntegration uses** `Derivation::DeriveCurve` for tangent vectors:

```cpp
// Internal helper (simplified):
auto tangent = Derivation::DeriveCurve<N>(curve, t, nullptr);
Real ds = tangent.NormL2();  // Arc length element
```

**Numerical differentiation** → automatic, no analytical derivative needed!

**Accuracy**: Controlled by `Derivation` settings (typically 1e-6)

### With FieldOperations

**Compute force from potential** for work integrals:

```cpp
#include "core/FieldOperations.h"

ScalarFunction<3> potential([](const VectorN<Real,3>& r) {
    return -1.0 / r.NormL2();
});

// Force = -∇φ (automatically computed)
VectorFunction<3> force([&potential](const VectorN<Real,3>& r) {
    return ScalarFieldOperations::GradientCart(potential, r);
});

// Use in work integral
Real work = PathIntegration::LineIntegral(force, curve, t1, t2);
```

**Advantage**: No need for analytical gradient formula!

### With Integration Module

**Uses** `IntegrateTrap` from core integration:

```cpp
// Internal implementation
return IntegrateTrap(helper_function, a, b, nullptr, nullptr, eps);
```

**Adaptive trapezoidal rule**:
- Refines grid until precision met
- `eps` parameter controls accuracy
- Default: `Defaults::LineIntegralPrecision` (1e-5)

---

## Best Practices

### Choosing Integration Precision

```cpp
// Default precision (usually sufficient)
Real result1 = PathIntegration::LineIntegral(field, curve, 0, 1);

// Higher precision for smooth fields
Real result2 = PathIntegration::LineIntegral(field, curve, 0, 1, 1e-8);

// Lower precision for rough estimates
Real result3 = PathIntegration::LineIntegral(field, curve, 0, 1, 1e-3);
```

**Guidelines**:
- **Smooth fields**: eps = 1e-6 to 1e-8
- **General use**: eps = 1e-5 (default)
- **Quick estimates**: eps = 1e-3 to 1e-4
- **Singular points**: May need adaptive mesh

### Parametrization Matters

❌ **Bad**: Non-uniform parametrization
```cpp
// Clustered points near t=0
ParametricCurve<3> bad_circle([](Real t) {
    Real s = t * t;  // Non-linear reparametrization
    return VectorN<Real,3>{cos(s), sin(s), 0};
});
// Integration will be inefficient!
```

✅ **Good**: Uniform or smooth parametrization
```cpp
ParametricCurve<3> good_circle([](Real t) {
    return VectorN<Real,3>{cos(t), sin(t), 0};
});
// |r'(t)| = constant → efficient integration
```

**Best**: Arc-length parametrization (if possible)
- |**r**'(s)| = 1
- Simplest integrands
- Most efficient numerically

### Verify with Analytical Results

```cpp
// Always test with known solutions
ParametricCurve<3> circle([](Real t) {
    return VectorN<Real,3>{cos(t), sin(t), 0};
});

Real computed = PathIntegration::ParametricCurveLength(circle, 0, 2*Constants::PI);
Real expected = 2 * Constants::PI;

Real error = std::abs(computed - expected);
assert(error < 1e-5);  // Verify accuracy
```

### Common Pitfalls

❌ **Pitfall 1**: Dimension mismatch
```cpp
// Current LineIntegral only supports N=3
ScalarFunction<2> field2D(...);  // 2D field
ParametricCurve<2> curve2D(...);  // 2D curve

// Won't compile!
// Real result = PathIntegration::LineIntegral(field2D, curve2D, 0, 1);
```

**Solution**: Embed in 3D
```cpp
ParametricCurve<3> curve3D([&curve2D](Real t) {
    auto p = curve2D(t);
    return VectorN<Real,3>{p[0], p[1], 0};
});
```

❌ **Pitfall 2**: Reversed parameter order
```cpp
// Integration from t2 to t1 (backwards!)
Real result = PathIntegration::LineIntegral(field, curve, 10.0, 0.0);
// Will give negative of expected result!
```

✅ **Correct**: Always t1 < t2
```cpp
Real result = PathIntegration::LineIntegral(field, curve, 0.0, 10.0);
```

❌ **Pitfall 3**: Singular points in field
```cpp
// Field singular at origin
ScalarFunction<3> singular([](const VectorN<Real,3>& r) {
    return 1.0 / r.NormL2();  // Blows up at r=0!
});

// If curve passes through origin → NaN or huge error
```

**Solution**: Avoid or handle singularities
```cpp
ScalarFunction<3> regularized([](const VectorN<Real,3>& r) {
    Real r_mag = r.NormL2();
    if (r_mag < 1e-10) return 1e10;  // Cap at origin
    return 1.0 / r_mag;
});
```

### Performance Optimization

**Reduce integrand complexity**:
```cpp
// SLOW: Complex field evaluation
VectorFunction<3> complex_field([](const VectorN<Real,3>& r) {
    // Expensive computations here
    return some_expensive_calculation(r);
});

// BETTER: Precompute or simplify
```

**Use analytical derivatives when available**:
```cpp
// Current implementation uses numerical derivatives (slower)
// If you have analytical r'(t), consider custom implementation
```

**Cache field values** for repeated queries (advanced):
```cpp
// For fields that are expensive to evaluate
// Consider lookup tables or interpolation
```

---

## Advanced Topics

### Green's Theorem (Closed Loops in Plane)

**For closed curve C in xy-plane**:
```
∮_C (P dx + Q dy) = ∬_R (∂Q/∂x - ∂P/∂y) dA
```

**Vector form**:
```
∮_C F·dr = ∬_R (∇ × F)·n dA
```

**Application**: Convert line integral to area integral

**Example** (circulation):
```cpp
// Field: F = (-y, x, 0) → curl = (0, 0, 2)
// Circle of radius R
// ∮ F·dr = ∬ 2 dA = 2 * πR² (by Green's theorem)
```

### Stokes' Theorem (3D Surfaces)

**For closed curve C bounding surface S**:
```
∮_C F·dr = ∬_S (∇ × F)·n dS
```

**Generalization of Green's theorem** to 3D

**Not directly implemented** (requires surface integration)

### Fundamental Theorem for Line Integrals

**For conservative field** **F** = ∇φ:
```
∫_C F·dr = φ(end) - φ(start)
```

**Path independence!**

**Check if field is conservative**:
```cpp
// In 3D: F = (P, Q, R) is conservative if ∇ × F = 0
// i.e., ∂R/∂y = ∂Q/∂z, ∂P/∂z = ∂R/∂x, ∂Q/∂x = ∂P/∂y
```

**Useful for verification**:
```cpp
VectorN<Real,3> start = curve(t1);
VectorN<Real,3> end = curve(t2);

Real work_integral = PathIntegration::LineIntegral(force, curve, t1, t2);
Real potential_diff = potential(end) - potential(start);

// Should be equal (within numerical error)
assert(std::abs(work_integral - potential_diff) < 1e-4);
```

### Complex Integration (Future)

**Path integrals in complex plane**:
```
∫_C f(z) dz
```

**Currently not implemented** (would require complex number support)

**Application**: Residue theorem, conformal mapping

### Surface and Volume Integrals (Future)

**Surface integral**:
```
∬_S f dS  (scalar)
∬_S F·n dS  (vector, flux)
```

**Volume integral**:
```
∭_V f dV
```

**Status**: Not yet implemented (require surface/volume parametrizations)

---

## Summary

### Key Methods

| Method | Purpose | Formula |
|--------|---------|---------|
| `ParametricCurveLength` | Arc length | ∫ \|**r**'\| dt |
| `ParametricCurveMass` | Weighted length | ∫ ρ(t) \|**r**'\| dt |
| `LineIntegral` (scalar) | Scalar field integral | ∫ f \|**r**'\| dt |
| `LineIntegral` (vector) | Work/circulation | ∫ **F**·**r**' dt |

### Geometric Insights

1. ✅ **Arc length** = distance along curve (independent of parametrization)
2. ✅ **Scalar line integral** = weighted sum along path
3. ✅ **Vector line integral** = component of field along path
4. ✅ **Conservative fields** → path-independent work
5. ✅ **Closed loops** in conservative fields → zero circulation
6. ✅ **Numerical integration** via adaptive trapezoidal rule
7. ✅ **Automatic tangent** via numerical differentiation

### When to Use

**Path Integration when**:
- Computing curve lengths (cables, tracks, orbits)
- Calculating work done by forces
- Finding total heat/mass along wire
- Verifying field properties (conservative, rotational)
- Physics simulations (particles, fluids)

**Method choice**:
- Arc length: Pure geometry
- Scalar integral: Field distribution (temperature, density)
- Vector integral: Force/flow effects (work, circulation)

### Implementation Notes

**Strengths**:
- Automatic numerical differentiation
- Works for any parametric curve
- Adaptive precision control
- Clean, simple API

**Limitations**:
- LineIntegral currently 3D only (N=3)
- No surface/volume integrals yet
- Numerical (not symbolic)
- May struggle with singularities

### References

- **Vector Calculus** (Marsden & Tromba)
- **Calculus** (Stewart): Chapter on line integrals
- **Mathematical Methods** (Boas): Green's and Stokes' theorems
- **Numerical Recipes**: Chapter on integration
- **Classical Mechanics** (Goldstein): Work and conservative forces

---

## Runnable Examples

Working code examples are available in `src/docs_demos/docs_demo_integration_path.cpp`:

### Curve Length
- `Docs_Demo_Calc_curve_length()` - Calculating arc length of curves

### Work Integral
- `Docs_Demo_Calc_work_integral()` - Computing line integrals of vector fields

---

**Part of MinimalMathLibrary** - `mml/core/Integration/PathIntegration.h`
