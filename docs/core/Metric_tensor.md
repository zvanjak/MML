# Metric Tensor

**Location**: `mml/core/MetricTensor.h`  
**Dependencies**: Derivation, CoordTransf, Tensor, Function

---

## Overview

The metric tensor framework provides tools for working with distances, angles, and geometry in arbitrary coordinate systems and on manifolds. Metric tensors encode the geometry of a space, enabling measurement of lengths, angles, areas, and volumes in both flat and curved spaces.

### Key Capabilities

- **Predefined Metrics**: Cartesian, spherical, cylindrical, Minkowski spacetime
- **Custom Metrics**: Build from coordinate transformations automatically
- **Christoffel Symbols**: First and second kind for covariant differentiation
- **Covariant Derivatives**: For vectors (contravariant and covariant components)
- **Geometric Calculations**: Arc length, geodesics, curvature (via Christoffel symbols)
- **General Relativity**: Minkowski metric for special relativity

### Physical Applications

- **Differential Geometry**: Curved surfaces, manifolds, Riemannian geometry
- **General Relativity**: Spacetime metrics, geodesics, Einstein field equations
- **Classical Mechanics**: Lagrangian mechanics in generalized coordinates
- **Continuum Mechanics**: Strain tensors, stress analysis
- **Fluid Dynamics**: Non-Cartesian flow analysis

---

## Quick Reference

### Predefined Metric Tensors

| Coordinate System | Class | Dimension | Components |
|-------------------|-------|-----------|------------|
| Cartesian (3D) | `MetricTensorCartesian3D` | 3 | δᵢⱼ (identity) |
| Spherical (3D) | `MetricTensorSpherical` | 3 | diag(1, r², r²sin²θ) |
| Spherical Contravar | `MetricTensorSphericalContravar` | 3 | diag(1, 1/r², 1/(r²sin²θ)) |
| Cylindrical (3D) | `MetricTensorCylindrical` | 3 | diag(1, ρ², 1) |
| Minkowski (4D) | `MetricTensorMinkowski` | 4 | diag(-1, 1, 1, 1) |
| From Transform | `MetricTensorFromCoordTransf<From,To,N>` | N | Computed from Jacobian |

### Key Operations

| Operation | Method | Use Case |
|-----------|--------|----------|
| Get Covariant Metric | `GetCovariantMetric(pos)` | g_ij components |
| Get Contravariant Metric | `GetContravariantMetric(pos)` | g^ij components |
| Christoffel Symbols (1st) | `GetChristoffelSymbolFirstKind(i,j,k,pos)` | Γᵢⱼₖ |
| Christoffel Symbols (2nd) | `GetChristoffelSymbolSecondKind(i,j,k,pos)` | Γⁱⱼₖ |
| Covariant Derivative | `CovariantDerivativeContravar(func,j,pos)` | ∇ⱼVⁱ |

---

## Mathematical Background

### What is a Metric Tensor?

The **metric tensor** g encodes the geometry of a space, defining:
- **Distances**: ds² = gᵢⱼ dxⁱ dxʲ (line element)
- **Angles**: cos(θ) = (gᵢⱼ vⁱ wʲ) / √[(gᵢⱼ vⁱ vʲ)(gₖₗ wᵏ wˡ)]
- **Volumes**: dV = √det(g) dx¹...dxⁿ

**Covariant Components** gᵢⱼ:
```
gᵢⱼ = eᵢ · eⱼ
```
where eᵢ = ∂r/∂xⁱ are basis vectors.

**Contravariant Components** g^ij:
```
g^ij = eⁱ · eʲ
```
where eⁱ are reciprocal basis vectors (g^ik g_kj = δⁱⱼ).

**Key Property**: g^ij = (gᵢⱼ)^(-1) (matrix inverse)

### Christoffel Symbols

**Christoffel symbols** Γ encode how basis vectors change as you move through space.

**First Kind** Γᵢⱼₖ:
```
Γᵢⱼₖ = ½(∂gᵢₖ/∂xʲ + ∂gⱼₖ/∂xⁱ - ∂gᵢⱼ/∂xᵏ)
```

**Second Kind** Γⁱⱼₖ:
```
Γⁱⱼₖ = ½ gⁱˡ(∂gₗₖ/∂xʲ + ∂gⱼₗ/∂xᵏ - ∂gⱼₖ/∂xˡ)
     = gⁱˡ Γₗⱼₖ
```

**Physical Meaning**: Γⁱⱼₖ describes how the i-th basis vector component changes when moving in the k-th direction.

### Covariant Derivative

The **covariant derivative** ∇ is the generalization of the ordinary derivative to curved spaces.

**For Contravariant Vectors** Vⁱ:
```
∇ⱼVⁱ = ∂ⱼVⁱ + ΓⁱₖⱼVᵏ
```

**For Covariant Vectors** Vᵢ:
```
∇ⱼVᵢ = ∂ⱼVᵢ - ΓᵏᵢⱼVₖ
```

**Purpose**: Ensures derivatives transform as tensors (unlike ordinary derivatives in curvilinear coordinates).

---

## Core Classes

### MetricTensorField - Base Class

```cpp
template<int N>
class MetricTensorField : public ITensorField2<N>
{
public:
    // Constructors
    MetricTensorField();  // Default: (2,0) - covariant
    MetricTensorField(int numContra, int numCo);
    
    // Pure virtual - must implement in derived class
    virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const = 0;
    
    // Tensor interface
    Tensor2<N> operator()(const VectorN<Real, N>& pos) const;
    
    // Metric tensor components
    MatrixNM<Real, N, N> GetCovariantMetric(const VectorN<Real, N>& pos) const;
    MatrixNM<Real, N, N> GetContravariantMetric(const VectorN<Real, N>& pos) const;
    
    // Christoffel symbols
    Real GetChristoffelSymbolFirstKind(int i, int j, int k, 
                                       const VectorN<Real, N>& pos) const;
    Real GetChristoffelSymbolSecondKind(int i, int j, int k, 
                                        const VectorN<Real, N>& pos) const;
    
    // Covariant derivatives
    VectorN<Real, N> CovariantDerivativeContravar(const IVectorFunction<N>& func, 
                                                   int j, 
                                                   const VectorN<Real, N>& pos) const;
    Real CovariantDerivativeContravarComp(const IVectorFunction<N>& func, 
                                           int i, int j, 
                                           const VectorN<Real, N>& pos) const;
    
    VectorN<Real, N> CovariantDerivativeCovar(const IVectorFunction<N>& func, 
                                               int j, 
                                               const VectorN<Real, N>& pos) const;
    Real CovariantDerivativeCovarComp(const IVectorFunction<N>& func, 
                                       int i, int j, 
                                       const VectorN<Real, N>& pos) const;
};
```

**Purpose**: Abstract base providing all metric tensor operations. Derived classes only need to implement `Component(i,j,pos)`.

**Key Design**:
- Covariant metric g_ij computed directly from `Component()`
- Contravariant metric g^ij computed by matrix inversion
- Christoffel symbols computed numerically using `Derivation::DerivePartial()`
- Covariant derivatives use Christoffel symbols for proper tensor transformation

---

## Predefined Metric Tensors

### MetricTensorCartesian3D - Euclidean Space

```cpp
class MetricTensorCartesian3D : public MetricTensorField<3>
{
public:
    MetricTensorCartesian3D();
    Real Component(int i, int j, const VectorN<Real, 3>& pos) const;
};
```

**Metric Components**:
```
gᵢⱼ = δᵢⱼ = { 1 if i=j
            { 0 if i≠j
```

**Line Element**:
```
ds² = (dx)² + (dy)² + (dz)²
```

**Christoffel Symbols**: All zero (flat space)

**Example 1: Cartesian Metric**

```cpp
#include "core/MetricTensor.h"

// 3D Euclidean space
MetricTensorCartesian3D metricCart;

Vector3Cartesian pos{1.0, 2.0, 3.0};

// Get metric tensor at position
MatrixNM<Real, 3, 3> g = metricCart.GetCovariantMetric(pos);
// Result: identity matrix (all positions)

// Compute distance
Vector3Cartesian dr{0.1, 0.2, 0.3};
Real ds_squared = 0.0;
for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
        ds_squared += g[i][j] * dr[i] * dr[j];
// Result: ds² = 0.1² + 0.2² + 0.3² = 0.14
```

### MetricTensorSpherical - Spherical Coordinates

```cpp
class MetricTensorSpherical : public MetricTensorField<3>
{
public:
    MetricTensorSpherical();
    Real Component(int i, int j, const VectorN<Real, 3>& pos) const override;
};
```

**Coordinates**: (r, θ, φ) - radius, polar angle, azimuthal angle

**Metric Components**:
```
g_rr = 1
g_θθ = r²
g_φφ = r² sin²(θ)
g_ij = 0  (i ≠ j)
```

**Line Element**:
```
ds² = dr² + r²dθ² + r²sin²(θ)dφ²
```

**Volume Element**:
```
dV = r² sin(θ) dr dθ dφ
```

**Contravariant Components**:
```
g^rr = 1
g^θθ = 1/r²
g^φφ = 1/(r²sin²θ)
```

**Example 2: Spherical Metric**

```cpp
// Spherical metric
MetricTensorSpherical metricSpher;

Vector3Spherical pos{5.0, Constants::PI/4, Constants::PI/3};  // r=5, θ=45°, φ=60°

// Get covariant metric
MatrixNM<Real, 3, 3> g = metricSpher.GetCovariantMetric(pos);
// g[0][0] = 1
// g[1][1] = 25          (r² = 5²)
// g[2][2] = 12.5        (r²sin²θ = 25 * 0.5)

// Get contravariant metric
MatrixNM<Real, 3, 3> g_inv = metricSpher.GetContravariantMetric(pos);
// g_inv[0][0] = 1
// g_inv[1][1] = 0.04    (1/r² = 1/25)
// g_inv[2][2] = 0.08    (1/(r²sin²θ))

// Christoffel symbol example: Γʳθθ = -r
Real Gamma_r_theta_theta = 
    metricSpher.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
// Result: Γʳθθ ≈ -5.0
```

### MetricTensorCylindrical - Cylindrical Coordinates

```cpp
class MetricTensorCylindrical : public MetricTensorField<3>
{
public:
    MetricTensorCylindrical();
    Real Component(int i, int j, const VectorN<Real, 3>& pos) const override;
};
```

**Coordinates**: (ρ, φ, z) - radial distance from axis, azimuthal angle, height

**Metric Components**:
```
g_ρρ = 1
g_φφ = ρ²
g_zz = 1
g_ij = 0  (i ≠ j)
```

**Line Element**:
```
ds² = dρ² + ρ²dφ² + dz²
```

**Volume Element**:
```
dV = ρ dρ dφ dz
```

**Example 3: Cylindrical Metric**

```cpp
// Cylindrical metric
MetricTensorCylindrical metricCyl;

Vector3Cylindrical pos{3.0, Constants::PI/6, 2.0};  // ρ=3, φ=30°, z=2

// Get metric tensor
MatrixNM<Real, 3, 3> g = metricCyl.GetCovariantMetric(pos);
// g[0][0] = 1
// g[1][1] = 9     (ρ² = 3²)
// g[2][2] = 1

// Compute arc length along circle at constant ρ and z
Vector3Cylindrical dr{0.0, 0.1, 0.0};  // dρ=0, dφ=0.1 rad, dz=0
Real ds_squared = 0.0;
for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
        ds_squared += g[i][j] * dr[i] * dr[j];
Real ds = sqrt(ds_squared);
// Result: ds = ρ·dφ = 3 * 0.1 = 0.3
```

### MetricTensorMinkowski - Special Relativity

```cpp
class MetricTensorMinkowski : public MetricTensorField<4>
{
public:
    MetricTensorMinkowski();
    Real Component(int i, int j, const VectorN<Real, 4>& pos) const override;
};
```

**Coordinates**: (ct, x, y, z) - time (× c), spatial coordinates

**Metric Components** (signature -,+,+,+):
```
g_00 = -1
g_11 = +1
g_22 = +1
g_33 = +1
g_ij = 0  (i ≠ j)
```

**Spacetime Interval**:
```
ds² = -c²dt² + dx² + dy² + dz²
```

**Physical Interpretation**:
- ds² < 0: timelike separation (causal connection possible)
- ds² = 0: lightlike (photon path)
- ds² > 0: spacelike (causally disconnected)

**Example 4: Minkowski Metric**

```cpp
#include "core/MetricTensor.h"

// Minkowski spacetime metric
MetricTensorMinkowski metricMinkowski;

// Event in spacetime: ct=5, x=3, y=0, z=0
Vector4Minkowski event{5.0, 3.0, 0.0, 0.0};

// Get metric
MatrixNM<Real, 4, 4> eta = metricMinkowski.GetCovariantMetric(event);
// eta[0][0] = -1  (time component)
// eta[1][1] = +1  (spatial components)
// eta[2][2] = +1
// eta[3][3] = +1

// Compute spacetime interval
Vector4Minkowski displacement{2.0, 1.5, 0.0, 0.0};  // Δct=2, Δx=1.5
Real interval_squared = 0.0;
for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        interval_squared += eta[i][j] * displacement[i] * displacement[j];

// Result: ds² = -(2)² + (1.5)² = -4 + 2.25 = -1.75 (timelike)
Real proper_time = sqrt(-interval_squared);  // For timelike intervals
```

---

## MetricTensorFromCoordTransf - Build from Transformations

```cpp
template<typename VectorFrom, typename VectorTo, int N>
class MetricTensorFromCoordTransf : public MetricTensorField<N>
{
public:
    MetricTensorFromCoordTransf(ICoordTransfWithInverse<VectorFrom,VectorTo,N>& transf);
    Real Component(int i, int j, const VectorN<Real, N>& pos) const override;
};
```

**Purpose**: Automatically compute metric tensor from coordinate transformation.

**Formula**:
```
gᵢⱼ = Σₖ (∂xᵏ/∂qⁱ)(∂xᵏ/∂qʲ)
```

where x are Cartesian coordinates, q are curvilinear coordinates.

**Mathematical Background**:

If you have a transformation from curvilinear coordinates q to Cartesian x:
```
x = f(q)
```

The metric induced on the q-space is:
```
gᵢⱼ = J^T J
```

where J is the Jacobian matrix Jₖᵢ = ∂fᵏ/∂qⁱ.

**Example 5: Metric from Spherical Transformation**

```cpp
#include "core/CoordTransf/CoordTransfSpherical.h"
#include "core/MetricTensor.h"

// Use spherical→Cartesian transformation
CoordTransfSphericalToCartesian transf;

// Build metric automatically
MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> 
    metric(transf);

Vector3Spherical pos{4.0, Constants::PI/3, Constants::PI/4};  // r=4, θ=60°, φ=45°

// Get metric components
MatrixNM<Real, 3, 3> g = metric.GetCovariantMetric(pos);

// Verify spherical metric components
std::cout << "g_rr = " << g[0][0] << std::endl;    // Should be 1
std::cout << "g_θθ = " << g[1][1] << std::endl;    // Should be r² = 16
std::cout << "g_φφ = " << g[2][2] << std::endl;    // Should be r²sin²θ ≈ 12

// This matches MetricTensorSpherical!
```

**Example 6: Custom Coordinate System**

```cpp
// Define custom parabolic cylindrical coordinates
// x = σ·τ, y = ½(τ² - σ²), z = z

class CoordTransfParabolicToCart : 
    public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
{
    // Implementation of transformation...
};

CoordTransfParabolicToCart parabolicTransf;

// Automatically compute metric for parabolic coordinates
MetricTensorFromCoordTransf<Vector3Cartesian, Vector3Cartesian, 3>
    metricParabolic(parabolicTransf);

Vector3Cartesian pos{2.0, 1.5, 3.0};  // (σ, τ, z)

// Get induced metric
MatrixNM<Real, 3, 3> g = metricParabolic.GetCovariantMetric(pos);
// Automatically computed from Jacobian!
```

---

## Christoffel Symbols and Curvature

### Computing Christoffel Symbols

**Example 7: Christoffel Symbols in Spherical Coordinates**

```cpp
MetricTensorSpherical metric;

Vector3Spherical pos{5.0, Constants::PI/4, Constants::PI/3};

// Second kind Christoffel symbols (used in covariant derivatives)

// Γʳθθ = -r
Real Gamma_r_theta_theta = metric.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
std::cout << "Γʳθθ = " << Gamma_r_theta_theta << std::endl;  // ≈ -5.0

// Γʳφφ = -r sin²(θ)
Real Gamma_r_phi_phi = metric.GetChristoffelSymbolSecondKind(0, 2, 2, pos);
std::cout << "Γʳφφ = " << Gamma_r_phi_phi << std::endl;      // ≈ -2.5

// Γθᵣθ = 1/r
Real Gamma_theta_r_theta = metric.GetChristoffelSymbolSecondKind(1, 0, 1, pos);
std::cout << "Γθᵣθ = " << Gamma_theta_r_theta << std::endl;  // ≈ 0.2

// Γθφφ = -sin(θ)cos(θ)
Real Gamma_theta_phi_phi = metric.GetChristoffelSymbolSecondKind(1, 2, 2, pos);
std::cout << "Γθφφ = " << Gamma_theta_phi_phi << std::endl;  // ≈ -0.35

// Γφᵣφ = 1/r
Real Gamma_phi_r_phi = metric.GetChristoffelSymbolSecondKind(2, 0, 2, pos);
std::cout << "Γφᵣφ = " << Gamma_phi_r_phi << std::endl;      // ≈ 0.2

// Γφθφ = cot(θ)
Real Gamma_phi_theta_phi = metric.GetChristoffelSymbolSecondKind(2, 1, 2, pos);
std::cout << "Γφθφ = " << Gamma_phi_theta_phi << std::endl;  // ≈ 1.0
```

### First Kind Christoffel Symbols

```cpp
// Γᵣθθ (first kind)
Real Gamma_r_theta_theta_1st = 
    metric.GetChristoffelSymbolFirstKind(0, 1, 1, pos);

// Relation: Γᵢⱼₖ = gᵢₗ Γˡⱼₖ
// Verify: Γᵣθθ = g_rr * Γʳθθ = 1 * (-r) = -r
```

---

## Covariant Derivatives

### Contravariant Vector Fields

```cpp
// Define a velocity field in spherical coordinates
class VelocityFieldSpher : public IVectorFunction<3>
{
public:
    VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
    {
        Real r = pos[0];
        Real theta = pos[1];
        Real phi = pos[2];
        
        // Example: rotating flow
        VectorN<Real, 3> v;
        v[0] = 0.0;              // v^r = 0
        v[1] = 0.0;              // v^θ = 0
        v[2] = 1.0 / r;          // v^φ = 1/r
        return v;
    }
};
```

**Example 8: Covariant Derivative of Contravariant Vector**

```cpp
MetricTensorSpherical metric;
VelocityFieldSpher velocity;

Vector3Spherical pos{3.0, Constants::PI/2, 0.0};  // r=3, θ=90°, φ=0

// Compute covariant derivative ∇ⱼvⁱ for j=1 (∂/∂θ direction)
VectorN<Real, 3> covarDeriv_theta = 
    metric.CovariantDerivativeContravar(velocity, 1, pos);

// Components: ∇θv^r, ∇θv^θ, ∇θv^φ
std::cout << "∇θv^r = " << covarDeriv_theta[0] << std::endl;
std::cout << "∇θv^θ = " << covarDeriv_theta[1] << std::endl;
std::cout << "∇θv^φ = " << covarDeriv_theta[2] << std::endl;

// Compare with ordinary partial derivative
Real partialDeriv_theta_vphi = 
    Derivation::DeriveVecPartial<3>(velocity, 2, 1, pos, nullptr);
std::cout << "∂θv^φ = " << partialDeriv_theta_vphi << std::endl;

// Difference shows Christoffel symbol correction
```

### Covariant Vector Fields

**Example 9: Gradient as Covariant Vector**

```cpp
// Scalar field (e.g., temperature)
class TemperatureField : public IScalarFunction<3>
{
public:
    Real operator()(const VectorN<Real, 3>& pos) const override
    {
        Real r = pos[0];
        return 100.0 / (r * r);  // T = 100/r² (inverse square law)
    }
};

// Gradient is a covariant vector: ∇ᵢT = ∂ᵢT
// Its covariant derivative transforms properly
```

---

## Integration with Other Modules

### With Field Operations

```cpp
#include "core/FieldOperations.h"
#include "core/MetricTensor.h"

// Gradient in general coordinates uses metric tensor
MetricTensorSpherical metric;
TemperatureField temp;

Vector3Spherical pos{2.0, Constants::PI/4, Constants::PI/6};

// Gradient with metric tensor (uses contravariant metric)
VectorN<Real, 3> grad = FieldOps::Gradient(temp, pos, metric);

// This properly computes ∇T using g^ij ∂ⱼT
```

### With Coordinate Transformations

**Example 10: Full Workflow - Transformation to Metric**

```cpp
#include "core/CoordTransf/CoordTransfCylindrical.h"
#include "core/MetricTensor.h"
#include "core/FieldOperations.h"

// 1. Define coordinate transformation
CoordTransfCylindricalToCartesian cylToCart;

// 2. Build metric from transformation
MetricTensorFromCoordTransf<Vector3Cylindrical, Vector3Cartesian, 3>
    metric(cylToCart);

// 3. Define scalar field
auto potential = [](const VectorN<Real, 3>& pos) -> Real {
    Real rho = pos[0];
    return log(rho);  // Logarithmic potential
};
ScalarFunction<3> phi(potential);

// 4. Compute gradient using metric
Vector3Cylindrical pos{2.0, 0.0, 1.0};
VectorN<Real, 3> grad_phi = FieldOps::Gradient(phi, pos, metric);

// Result: ∇φ properly computed in cylindrical coordinates
```

---

## Geometric Applications

### Arc Length Calculation

**Example 11: Arc Length on Sphere**

```cpp
// Parametric curve on unit sphere: θ(t), φ(t)
class CurveOnSphere : public IParametricCurve<3>
{
public:
    Vector3Spherical operator()(Real t) const override
    {
        Real r = 1.0;                      // Unit sphere
        Real theta = Constants::PI/4;       // Constant latitude 45°
        Real phi = t;                       // Longitude varies
        return Vector3Spherical{r, theta, phi};
    }
};

MetricTensorSpherical metric;
CurveOnSphere curve;

// Compute arc length from t=0 to t=π/2
Real arcLength = 0.0;
int numSteps = 100;
Real dt = (Constants::PI/2) / numSteps;

for (int i = 0; i < numSteps; i++)
{
    Real t = i * dt;
    Vector3Spherical pos = curve(t);
    Vector3Spherical dpos_dt = curve.Derivative(t);  // Need to implement
    
    // ds² = g_ij dx^i/dt dx^j/dt dt²
    MatrixNM<Real, 3, 3> g = metric.GetCovariantMetric(pos);
    
    Real ds_dt_squared = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ds_dt_squared += g[i][j] * dpos_dt[i] * dpos_dt[j];
    
    arcLength += sqrt(ds_dt_squared) * dt;
}

// For circle at latitude 45° on unit sphere:
// Arc length = r * sin(θ) * Δφ = 1 * sin(π/4) * π/2 ≈ 1.11
```

### Geodesics (Conceptual)

Geodesics (shortest paths) satisfy the geodesic equation:
```
d²xⁱ/ds² + Γⁱⱼₖ (dxʲ/ds)(dxᵏ/ds) = 0
```

Using `GetChristoffelSymbolSecondKind()`, you can implement geodesic integration.

---

## Advanced Topics

### Raising and Lowering Indices

**Convert Contravariant to Covariant**:
```cpp
// Given contravariant vector v^i
VectorN<Real, 3> v_contra{1.0, 0.5, 0.2};

// Lower index: v_i = g_ij v^j
MetricTensorSpherical metric;
Vector3Spherical pos{2.0, Constants::PI/3, Constants::PI/4};
MatrixNM<Real, 3, 3> g_cov = metric.GetCovariantMetric(pos);

VectorN<Real, 3> v_cov;
for (int i = 0; i < 3; i++)
{
    v_cov[i] = 0.0;
    for (int j = 0; j < 3; j++)
        v_cov[i] += g_cov[i][j] * v_contra[j];
}
```

**Convert Covariant to Contravariant**:
```cpp
// Given covariant vector v_i
VectorN<Real, 3> v_cov{1.0, 2.0, 3.0};

// Raise index: v^i = g^ij v_j
MatrixNM<Real, 3, 3> g_contra = metric.GetContravariantMetric(pos);

VectorN<Real, 3> v_contra;
for (int i = 0; i < 3; i++)
{
    v_contra[i] = 0.0;
    for (int j = 0; j < 3; j++)
        v_contra[i] += g_contra[i][j] * v_cov[j];
}
```

### Determinant and Volume Elements

```cpp
// Compute √det(g) for volume element
MetricTensorSpherical metric;
Vector3Spherical pos{3.0, Constants::PI/4, 0.0};

MatrixNM<Real, 3, 3> g = metric.GetCovariantMetric(pos);
Real detG = g.Det();
Real sqrtDetG = sqrt(detG);

// For spherical: √det(g) = r² sin(θ)
// At pos: √det(g) = 9 * sin(π/4) ≈ 6.36
```

---

## Best Practices

### Choosing Coordinate Systems

**Use spherical metric** when:
- Problem has spherical symmetry
- Working with central potentials
- Analyzing planetary motion

**Use cylindrical metric** when:
- Problem has axial symmetry
- Working with pipes, wires, toroids
- Analyzing rotating systems

**Use Cartesian metric** when:
- No special symmetry
- Problem naturally rectangular
- Default/simplest choice

**Build from transformation** when:
- Using custom coordinates
- Need metric induced from embedding
- Exploring differential geometry

### Numerical Considerations

**Christoffel Symbol Accuracy**:
- Computed numerically using `Derivation::DerivePartial()`
- Accuracy depends on derivative order (default: NDer4)
- Use appropriate step size for your metric scale

**Singularities**:
- Spherical: θ=0, π (poles)
- Cylindrical: ρ=0 (axis)
- Minkowski: No coordinate singularities
- Avoid evaluation at singular points

**Matrix Inversion**:
- Contravariant metric computed by inverting covariant
- Ensure metric is non-degenerate (det(g) ≠ 0)
- Check condition number for stability

### Performance

**Cache Metric Components**:
```cpp
// Good: compute once
MatrixNM<Real, 3, 3> g = metric.GetCovariantMetric(pos);
for (/* many operations */) {
    // Use g repeatedly
}

// Bad: recompute every time
for (/* many operations */) {
    MatrixNM<Real, 3, 3> g = metric.GetCovariantMetric(pos);  // Slow!
}
```

**Christoffel Symbols**:
- Expensive to compute (multiple partial derivatives)
- Cache if needed at same position for multiple operations
- Consider analytical formulas for well-known metrics

---

## Common Patterns

### Pattern 1: Transform → Metric → Field Operations

```cpp
// Standard workflow for custom coordinates
CoordTransfCustom transf;
MetricTensorFromCoordTransf metric(transf);

ScalarFunction field = /* ... */;
VectorN<Real, N> grad = FieldOps::Gradient(field, pos, metric);
```

### Pattern 2: Covariant Derivative of Vector Field

```cpp
// Compute all components of ∇ⱼVⁱ
MetricTensorSpher metric;
VectorField velocity;
Vector3Spherical pos = /* ... */;

// Matrix: (∇ⱼVⁱ)
MatrixNM<Real, 3, 3> covarDeriv;
for (int j = 0; j < 3; j++)  // Derivative direction
{
    VectorN<Real, 3> deriv = 
        metric.CovariantDerivativeContravar(velocity, j, pos);
    for (int i = 0; i < 3; i++)
        covarDeriv[i][j] = deriv[i];
}
```

### Pattern 3: Geodesic Initial Value Problem

```cpp
// Set up geodesic equation as ODE system
// State: (x^i, dx^i/ds)
// ODEs: d²x^i/ds² = -Γ^i_jk dx^j/ds dx^k/ds
```

---

## Runnable Examples

**Source File**: [docs_demo_metric_tensor.cpp](../../src/docs_demos/docs_demo_metric_tensor.cpp)

### Predefined Metric Tensors

| Function | Description |
|----------|-------------|
| `Docs_Demo_MetricTensorCartesian()` | Euclidean metric (identity), covariant/contravariant, distance |
| `Docs_Demo_MetricTensorSpherical()` | Spherical metric, diagonal components, Christoffel symbols |
| `Docs_Demo_MetricTensorCylindrical()` | Cylindrical metric, arc length calculation |
| `Docs_Demo_MetricTensorMinkowski()` | Spacetime intervals, timelike/spacelike/lightlike |

### Metric From Transformations

| Function | Description |
|----------|-------------|
| `Docs_Demo_MetricFromCoordTransf()` | Build metric from Jacobian of coordinate transformation |

### Christoffel Symbols and Covariant Derivatives

| Function | Description |
|----------|-------------|
| `Docs_Demo_ChristoffelFirstKind()` | First kind symbols Γᵢⱼₖ, relation to second kind |
| `Docs_Demo_CovariantDerivative()` | ∇ⱼvⁱ for contravariant vectors |

### Geometric Applications

| Function | Description |
|----------|-------------|
| `Docs_Demo_LineElement()` | Arc length ds² at different latitudes |
| `Docs_Demo_VolumeElement()` | √det(g) calculation, integration weight |

### Related Examples

| Example | Source File | Description |
|---------|-------------|-------------|
| Coordinate Transformations | [docs_demo_coord_transf.cpp](../../src/docs_demos/docs_demo_coord_transf.cpp) | Basis vectors, covariant/contravariant transforms |
| Tensors | [docs_demo_tensors.cpp](../../src/docs_demos/docs_demo_tensors.cpp) | Tensor operations, Christoffel symbol representation |

---

## Summary

The metric tensor framework provides:

✅ **Predefined Metrics** - Cartesian, spherical, cylindrical, Minkowski  
✅ **Automatic Computation** - From coordinate transformations  
✅ **Christoffel Symbols** - Both kinds, numerically computed  
✅ **Covariant Derivatives** - Proper tensor calculus  
✅ **Geometric Tools** - Arc length, volume elements, index manipulation  
✅ **Integration** - With field operations, coordinate transformations  
✅ **Physical Applications** - Classical mechanics, GR, continuum mechanics  

**Related Documentation**:
- [Coordinate_transformations.md](Coordinate_transformations.md) - Build metrics from transformations
- [Vector_field_operations.md](Vector_field_operations.md) - Use metrics for field operations
- [../base/Tensors.md](../base/Tensors.md) - Tensor mathematics
- [Derivation.md](Derivation.md) - Numerical differentiation for Christoffel symbols
