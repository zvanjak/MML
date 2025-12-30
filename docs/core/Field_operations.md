# Scalar and Vector Field Operations

**Sources:**  
- `mml/core/FieldOperations.h` - Gradient, divergence, curl, Laplacian

Numerical calculation of fundamental field operations in various coordinate systems. These operations are essential for physics simulations, fluid dynamics, electromagnetism, and differential geometry.

## Overview

**Scalar Field Operations:**
- **Gradient** ∇f: Vector field pointing in direction of maximum increase
- **Laplacian** ∇²f = ∇·∇f: Scalar measure of local curvature

**Vector Field Operations:**
- **Divergence** ∇·**F**: Scalar measure of field spreading/source density
- **Curl** ∇×**F**: Vector measure of field circulation/rotation

**Coordinate System Support:**
- **Cartesian** (x, y, z): Standard Euclidean coordinates
- **Spherical** (r, θ, φ): Radial distance, polar angle, azimuthal angle
- **Cylindrical** (ρ, φ, z): Radial distance, azimuthal angle, height
- **General coordinates**: With metric tensor (see [Metric_tensor.md](Metric_tensor.md))

---

## Quick Reference

### Scalar Field Operations

| Operation | Cartesian | Spherical | Cylindrical | General |
|-----------|-----------|-----------|-------------|---------|
| Gradient ∇f | `GradientCart` | `GradientSpher` | `GradientCyl` | `Gradient` (with metric) |
| Laplacian ∇²f | `LaplacianCart` | `LaplacianSpher` | `LaplacianCyl` | — |

### Vector Field Operations

| Operation | Cartesian | Spherical | Cylindrical | General |
|-----------|-----------|-----------|-------------|---------|
| Divergence ∇·**F** | `DivCart` | `DivSpher` | `DivCyl` | `Divergence` (with metric) |
| Curl ∇×**F** | `CurlCart` | `CurlSpher` | `CurlCyl` | — |

---

## Gradient (Scalar → Vector)

Gradient of a scalar field **f** produces a vector field pointing in the direction of maximum increase.

**Mathematical Definition:**
- **Cartesian**: ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z)
- **Spherical**: ∇f = (∂f/∂r, (1/r)∂f/∂θ, (1/(r sin θ))∂f/∂φ)
- **Cylindrical**: ∇f = (∂f/∂ρ, (1/ρ)∂f/∂φ, ∂f/∂z)
- **General**: ∇f = gⁱʲ ∂ⱼf (using contravariant metric tensor)

### Cartesian Gradient

```cpp
template<int N>
static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, 
                                      const VectorN<Real, N>& pos);

template<int N>
static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, 
                                      const VectorN<Real, N>& pos, 
                                      int der_order);
```

**Parameters:**
- `scalarField`: Scalar function f(x)
- `pos`: Position at which to evaluate gradient
- `der_order` (optional): Derivative order (1, 2, 4, 6, or 8) - higher = more accurate

**Returns:** Vector pointing in direction of maximum increase

**Example:**
```cpp
using namespace ScalarFieldOperations;

// Scalar field: f(x,y,z) = x² + y² + z²
auto field = [](const Vec3Cart& pos) -> Real {
    return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
};
RealFunctionFromStdFunc<3> scalarField(field);

Vec3Cart pos(1.0, 2.0, 3.0);
Vec3Cart grad = GradientCart(scalarField, pos);

// Gradient is 2*pos = (2, 4, 6)
std::cout << "Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")\n";
```

### Spherical Gradient

```cpp
static Vec3Sph GradientSpher(const IScalarFunction<3>& scalarField, 
                              const Vec3Sph& pos);

static Vec3Sph GradientSpher(const IScalarFunction<3>& scalarField, 
                              const Vec3Sph& pos, 
                              int der_order);
```

**Coordinate System:**
- pos[0] = r (radial distance)
- pos[1] = θ (polar angle, 0 to π)
- pos[2] = φ (azimuthal angle, 0 to 2π)

**Formula:**
∇f = (∂f/∂r, (1/r)∂f/∂θ, (1/(r sin θ))∂f/∂φ)

**Example:**
```cpp
// Radial field: f(r,θ,φ) = 1/r
auto field = [](const Vec3Sph& pos) -> Real {
    return 1.0 / pos[0];
};
RealFunctionFromStdFunc<3> scalarField(field);

Vec3Sph pos(2.0, Constants::PI/4, 0.0);  // r=2, θ=π/4, φ=0
Vec3Sph grad = GradientSpher(scalarField, pos);

// Gradient: (-1/r², 0, 0)
std::cout << "Radial component: " << grad[0] << "\n";  // -0.25
```

### Cylindrical Gradient

```cpp
static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, 
                            const Vec3Cyl& pos);

static Vec3Cyl GradientCyl(const IScalarFunction<3>& scalarField, 
                            const Vec3Cyl& pos, 
                            int der_order);
```

**Coordinate System:**
- pos[0] = ρ (radial distance from z-axis)
- pos[1] = φ (azimuthal angle)
- pos[2] = z (height)

**Formula:**
∇f = (∂f/∂ρ, (1/ρ)∂f/∂φ, ∂f/∂z)

**Example:**
```cpp
// Cylindrical field: f(ρ,φ,z) = ρ² + z²
auto field = [](const Vec3Cyl& pos) -> Real {
    return pos[0]*pos[0] + pos[2]*pos[2];
};
RealFunctionFromStdFunc<3> scalarField(field);

Vec3Cyl pos(2.0, Constants::PI/2, 1.0);
Vec3Cyl grad = GradientCyl(scalarField, pos);

// Gradient: (2ρ, 0, 2z) = (4, 0, 2)
```

### General Coordinate Gradient

```cpp
template<int N>
static VectorN<Real, N> Gradient(IScalarFunction<N>& scalarField, 
                                  const VectorN<Real, N>& pos, 
                                  const MetricTensorField<N>& metricTensorField);
```

**Purpose:** Compute gradient in arbitrary coordinate systems using metric tensor.

**Formula:** ∇f = gⁱʲ ∂ⱼf  
- Compute partial derivatives ∂ᵢf (covariant components)
- Use contravariant metric tensor gⁱʲ to raise indices

**Example:**
```cpp
// Define metric tensor for coordinate system
MetricTensorField<3> metric = /* ... */;

// Compute gradient
Vec3 grad = Gradient(scalarField, pos, metric);
```

See [Metric_tensor.md](Metric_tensor.md) for metric tensor details.

---

## Laplacian (Scalar → Scalar)

Laplacian measures the local "curvature" of a scalar field. Physically represents diffusion rate.

**Mathematical Definition:**
- ∇²f = ∇·∇f (divergence of gradient)
- **Cartesian**: ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²

### Cartesian Laplacian

```cpp
template<int N>
static Real LaplacianCart(const IScalarFunction<N>& scalarField, 
                           const VectorN<Real, N>& pos);
```

**Formula:**
∇²f = ∑ᵢ ∂²f/∂xᵢ²

**Example:**
```cpp
using namespace ScalarFieldOperations;

// Harmonic function: f(x,y) = x² - y²
auto field = [](const Vec2Cart& pos) -> Real {
    return pos[0]*pos[0] - pos[1]*pos[1];
};
RealFunctionFromStdFunc<2> scalarField(field);

Vec2Cart pos(1.0, 2.0);
Real lapl = LaplacianCart(scalarField, pos);

// Laplacian: 2 - 2 = 0 (harmonic!)
std::cout << "Laplacian: " << lapl << "\n";
```

### Spherical Laplacian

```cpp
static Real LaplacianSpher(const IScalarFunction<3>& scalarField, 
                            const Vec3Sph& pos);
```

**Formula:**
∇²f = (1/r²)∂/∂r(r²∂f/∂r) + (1/(r²sin θ))∂/∂θ(sin θ ∂f/∂θ) + (1/(r²sin²θ))∂²f/∂φ²

**Example:**
```cpp
// Radial function: f(r,θ,φ) = 1/r
auto field = [](const Vec3Sph& pos) -> Real {
    return 1.0 / pos[0];
};
RealFunctionFromStdFunc<3> scalarField(field);

Vec3Sph pos(2.0, Constants::PI/4, 0.0);
Real lapl = LaplacianSpher(scalarField, pos);

// For 1/r: ∇²(1/r) = 0 for r ≠ 0 (harmonic)
std::cout << "Laplacian: " << lapl << "\n";
```

### Cylindrical Laplacian

```cpp
static Real LaplacianCyl(const IScalarFunction<3>& scalarField, 
                          const Vec3Cyl& pos);
```

**Formula:**
∇²f = (1/ρ)∂/∂ρ(ρ ∂f/∂ρ) + (1/ρ²)∂²f/∂φ² + ∂²f/∂z²

**Example:**
```cpp
// Cylindrical harmonic: f(ρ,φ,z) = ρ²cos(2φ)
auto field = [](const Vec3Cyl& pos) -> Real {
    return pos[0]*pos[0] * cos(2*pos[1]);
};
RealFunctionFromStdFunc<3> scalarField(field);

Vec3Cyl pos(1.0, 0.0, 0.0);
Real lapl = LaplacianCyl(scalarField, pos);
```

---

## Divergence (Vector → Scalar)

Divergence measures the "outflow" or "source density" of a vector field.

**Mathematical Definition:**
- **Cartesian**: ∇·**F** = ∂Fₓ/∂x + ∂Fᵧ/∂y + ∂F_z/∂z
- **Spherical**: ∇·**F** = (1/r²)∂(r²Fᵣ)/∂r + (1/(r sin θ))∂(sin θ Fθ)/∂θ + (1/(r sin θ))∂Fφ/∂φ
- **Cylindrical**: ∇·**F** = (1/ρ)∂(ρFρ)/∂ρ + (1/ρ)∂Fφ/∂φ + ∂F_z/∂z
- **General**: ∇·**F** = ∂ᵢFⁱ + Fᵏ Γⁱᵢₖ (with Christoffel symbols)

### Cartesian Divergence

```cpp
template<int N>
static Real DivCart(const IVectorFunction<N>& vectorField, 
                     const VectorN<Real, N>& pos);
```

**Formula:**
∇·**F** = ∑ᵢ ∂Fᵢ/∂xᵢ

**Physical Interpretation:**
- **Positive**: Net outflow (source)
- **Negative**: Net inflow (sink)
- **Zero**: Incompressible flow (e.g., fluid dynamics)

**Example:**
```cpp
using namespace VectorFieldOperations;

// Radial field: F(x,y,z) = (x, y, z)
auto field = [](const Vec3Cart& pos) -> Vec3Cart {
    return pos;
};
VectorFunctionFromStdFunc<3> vectorField(field);

Vec3Cart pos(1.0, 2.0, 3.0);
Real div = DivCart(vectorField, pos);

// Divergence: 1 + 1 + 1 = 3 (expanding source)
std::cout << "Divergence: " << div << "\n";
```

### Spherical Divergence

```cpp
static Real DivSpher(const IVectorFunction<3>& vectorField, 
                      const VectorN<Real, 3>& pos);
```

**Formula:**
∇·**F** = (1/r²)∂(r²Fᵣ)/∂r + (1/(r sin θ))∂(sin θ Fθ)/∂θ + (1/(r sin θ))∂Fφ/∂φ

**Example:**
```cpp
// Radial field in spherical coords: F = (r, 0, 0)
auto field = [](const Vec3Sph& pos) -> Vec3Sph {
    return Vec3Sph(pos[0], 0.0, 0.0);
};
VectorFunctionFromStdFunc<3> vectorField(field);

Vec3Sph pos(2.0, Constants::PI/4, 0.0);
Real div = DivSpher(vectorField, pos);

// Divergence: ∂(r²·r)/∂r / r² = 3
```

### Cylindrical Divergence

```cpp
static Real DivCyl(const IVectorFunction<3>& vectorField, 
                    const VectorN<Real, 3>& pos);
```

**Formula:**
∇·**F** = (1/ρ)∂(ρFρ)/∂ρ + (1/ρ)∂Fφ/∂φ + ∂F_z/∂z

**Example:**
```cpp
// Cylindrical radial field: F = (ρ, 0, 0)
auto field = [](const Vec3Cyl& pos) -> Vec3Cyl {
    return Vec3Cyl(pos[0], 0.0, 0.0);
};
VectorFunctionFromStdFunc<3> vectorField(field);

Vec3Cyl pos(2.0, 0.0, 1.0);
Real div = DivCyl(vectorField, pos);

// Divergence: ∂(ρ·ρ)/∂ρ / ρ = 2
```

### General Coordinate Divergence

```cpp
template<int N>
static Real Divergence(const IVectorFunction<N>& vectorField, 
                        const VectorN<Real, N>& pos,
                        const MetricTensorField<N>& metricTensorField);
```

**Purpose:** Compute divergence in arbitrary coordinate systems.

**Formula:** ∇·**F** = ∂ᵢFⁱ + Fᵏ Γⁱᵢₖ  
- Uses Christoffel symbols Γⁱᵢₖ from metric tensor
- Accounts for non-Cartesian coordinate system effects

---

## Curl (Vector → Vector)

Curl measures the "circulation" or "rotation" of a vector field. Only defined in 3D.

**Mathematical Definition:**
- **Cartesian**: ∇×**F** = (∂F_z/∂y - ∂Fᵧ/∂z, ∂Fₓ/∂z - ∂F_z/∂x, ∂Fᵧ/∂x - ∂Fₓ/∂y)

### Cartesian Curl

```cpp
static Vec3Cart CurlCart(const IVectorFunction<3>& vectorField, 
                          const VectorN<Real, 3>& pos);
```

**Formula:**
```
∇×F = | î      ĵ      k̂     |
      | ∂/∂x   ∂/∂y   ∂/∂z  |
      | Fₓ     Fᵧ     F_z   |
```

**Physical Interpretation:**
- **Non-zero**: Rotational field (e.g., magnetic field, vorticity)
- **Zero**: Conservative/irrotational field (e.g., gradient fields)

**Example:**
```cpp
using namespace VectorFieldOperations;

// Rotational field: F(x,y,z) = (-y, x, 0)
auto field = [](const Vec3Cart& pos) -> Vec3Cart {
    return Vec3Cart(-pos[1], pos[0], 0.0);
};
VectorFunctionFromStdFunc<3> vectorField(field);

Vec3Cart pos(1.0, 2.0, 0.0);
Vec3Cart curl = CurlCart(vectorField, pos);

// Curl: (0, 0, 2) - constant rotation around z-axis
std::cout << "Curl: (" << curl[0] << ", " << curl[1] << ", " << curl[2] << ")\n";
```

### Spherical Curl

```cpp
static Vec3Sph CurlSpher(const IVectorFunction<3>& vectorField, 
                          const VectorN<Real, 3>& pos);
```

**Formula (component r):**
∇×F|ᵣ = (1/(r sin θ))[∂(sin θ Fφ)/∂θ - ∂Fθ/∂φ]

**Example:**
```cpp
// Spherical rotational field
auto field = [](const Vec3Sph& pos) -> Vec3Sph {
    return Vec3Sph(0.0, 0.0, pos[0] * sin(pos[1]));
};
VectorFunctionFromStdFunc<3> vectorField(field);

Vec3Sph pos(1.0, Constants::PI/2, 0.0);
Vec3Sph curl = CurlSpher(vectorField, pos);
```

### Cylindrical Curl

```cpp
static Vec3Cyl CurlCyl(const IVectorFunction<3>& vectorField, 
                        const VectorN<Real, 3>& pos);
```

**Formula (component ρ):**
∇×F|ρ = (1/ρ)∂F_z/∂φ - ∂Fφ/∂z

**Example:**
```cpp
// Cylindrical vortex: F = (0, ρ, 0)
auto field = [](const Vec3Cyl& pos) -> Vec3Cyl {
    return Vec3Cyl(0.0, pos[0], 0.0);
};
VectorFunctionFromStdFunc<3> vectorField(field);

Vec3Cyl pos(2.0, 0.0, 0.0);
Vec3Cyl curl = CurlCyl(vectorField, pos);

// Curl: (0, 0, 2) - vorticity
```

---

## Examples

### Example 1: Electric Potential and Field

```cpp
#include "core/FieldOperations.h"

void Example_Electric_Field()
{
    using namespace ScalarFieldOperations;

    // Point charge potential: φ(r) = k·Q/r
    Real k = 8.99e9;  // Coulomb constant
    Real Q = 1.0e-9;  // 1 nC charge

    auto potential = [k, Q](const Vec3Cart& pos) -> Real {
        Real r = pos.NormL2();
        return (r > 1e-10) ? k * Q / r : 0.0;
    };
    RealFunctionFromStdFunc<3> phi(potential);

    // Electric field: E = -∇φ
    Vec3Cart pos(1.0, 0.0, 0.0);  // 1 meter from charge
    Vec3Cart E = -1.0 * GradientCart(phi, pos);

    std::cout << "Electric field at (1,0,0): ";
    std::cout << "(" << E[0] << ", " << E[1] << ", " << E[2] << ") N/C\n";

    // Laplacian: ∇²φ = -ρ/ε₀ (Poisson equation)
    Real lapl = LaplacianCart(phi, pos);
    std::cout << "Laplacian: " << lapl << "\n";  // ~0 (away from charge)
}
```

### Example 2: Fluid Incompressibility

```cpp
#include "core/FieldOperations.h"

void Example_Incompressible_Flow()
{
    using namespace VectorFieldOperations;

    // 2D incompressible flow: v = (-y, x)
    auto velocity = [](const Vec2Cart& pos) -> Vec2Cart {
        return Vec2Cart(-pos[1], pos[0]);
    };
    VectorFunctionFromStdFunc<2> v(velocity);

    // Check incompressibility: ∇·v = 0
    Vec2Cart pos(2.0, 3.0);
    Real div = DivCart(v, pos);

    std::cout << "Divergence: " << div << "\n";
    
    if (std::abs(div) < 1e-6)
        std::cout << "Flow is incompressible!\n";

    // For 3D version, compute vorticity: ω = ∇×v
    auto velocity3D = [](const Vec3Cart& pos) -> Vec3Cart {
        return Vec3Cart(-pos[1], pos[0], 0.0);
    };
    VectorFunctionFromStdFunc<3> v3D(velocity3D);

    Vec3Cart pos3D(2.0, 3.0, 0.0);
    Vec3Cart vorticity = CurlCart(v3D, pos3D);

    std::cout << "Vorticity: (" << vorticity[0] << ", " 
              << vorticity[1] << ", " << vorticity[2] << ")\n";
}
```

### Example 3: Heat Equation

```cpp
#include "core/FieldOperations.h"

void Example_Heat_Diffusion()
{
    using namespace ScalarFieldOperations;

    // Temperature field: T(x,y,t) = T₀ exp(-α·t) sin(πx) sin(πy)
    Real T0 = 100.0;
    Real alpha = 0.01;
    Real t = 1.0;

    auto temperature = [T0, alpha, t](const Vec2Cart& pos) -> Real {
        return T0 * std::exp(-alpha * t) * 
               std::sin(Constants::PI * pos[0]) * 
               std::sin(Constants::PI * pos[1]);
    };
    RealFunctionFromStdFunc<2> T(temperature);

    // Heat equation: ∂T/∂t = κ·∇²T
    Vec2Cart pos(0.5, 0.5);
    Real laplacian = LaplacianCart(T, pos);

    std::cout << "Temperature: " << T(pos) << " K\n";
    std::cout << "Laplacian: " << laplacian << " K/m²\n";
    
    // Rate of change: ∂T/∂t = κ·∇²T
    Real kappa = 1.0;  // Thermal diffusivity
    Real dT_dt = kappa * laplacian;
    
    std::cout << "Rate of temperature change: " << dT_dt << " K/s\n";
}
```

### Example 4: Magnetic Field from Current

```cpp
#include "core/FieldOperations.h"

void Example_Magnetic_Field()
{
    using namespace VectorFieldOperations;

    // Magnetic vector potential: A = (0, 0, ln(ρ)) for wire along z-axis
    auto vectorPotential = [](const Vec3Cyl& pos) -> Vec3Cyl {
        return Vec3Cyl(0.0, 0.0, std::log(pos[0]));
    };
    VectorFunctionFromStdFunc<3> A(vectorPotential);

    // Magnetic field: B = ∇×A
    Vec3Cyl pos(2.0, 0.0, 0.0);
    Vec3Cyl B = CurlCyl(A, pos);

    std::cout << "Magnetic field: B_ρ=" << B[0] 
              << ", B_φ=" << B[1] 
              << ", B_z=" << B[2] << "\n";

    // For straight wire: B_φ = μ₀I/(2πρ)
    // Here we get B_φ = 1/ρ (in units where μ₀I/(2π) = 1)
}
```

### Example 5: Conservative Field Check

```cpp
#include "core/FieldOperations.h"

void Example_Conservative_Field()
{
    using namespace VectorFieldOperations;

    // Test if field is conservative: ∇×F = 0?
    auto field = [](const Vec3Cart& pos) -> Vec3Cart {
        // Gradient of f(x,y,z) = x²y + yz²
        return Vec3Cart(2*pos[0]*pos[1], 
                         pos[0]*pos[0] + pos[2]*pos[2],
                         2*pos[1]*pos[2]);
    };
    VectorFunctionFromStdFunc<3> F(field);

    Vec3Cart pos(1.0, 2.0, 3.0);
    Vec3Cart curl = CurlCart(F, pos);

    std::cout << "Curl: (" << curl[0] << ", " << curl[1] << ", " << curl[2] << ")\n";

    if (curl.NormL2() < 1e-10)
        std::cout << "Field is conservative (curl = 0)!\n";
    else
        std::cout << "Field is NOT conservative.\n";
}
```

### Example 6: Spherical Harmonic Analysis

```cpp
#include "core/FieldOperations.h"

void Example_Spherical_Harmonics()
{
    using namespace ScalarFieldOperations;

    // Spherical harmonic: Y₁₀ = √(3/(4π)) cos(θ)
    Real normalization = std::sqrt(3.0 / (4.0 * Constants::PI));
    
    auto harmonic = [normalization](const Vec3Sph& pos) -> Real {
        return normalization * std::cos(pos[1]);
    };
    RealFunctionFromStdFunc<3> Y10(harmonic);

    // Laplacian in spherical coordinates
    Vec3Sph pos(1.0, Constants::PI/4, 0.0);  // Unit sphere
    Real lapl = LaplacianSpher(Y10, pos);

    // For spherical harmonics: ∇²Yₗₘ = -l(l+1)/r² Yₗₘ
    // For l=1: ∇²Y₁₀ = -2 Y₁₀
    Real expected = -2.0 * Y10(pos);
    
    std::cout << "Laplacian: " << lapl << "\n";
    std::cout << "Expected: " << expected << "\n";
    std::cout << "Match: " << (std::abs(lapl - expected) < 1e-3 ? "YES" : "NO") << "\n";
}
```

---

## Vector Calculus Identities

### Fundamental Theorems

**Gradient Theorem:**
∫_C ∇f · d**r** = f(**b**) - f(**a**)

**Divergence Theorem (Gauss):**
∫_V (∇·**F**) dV = ∮_S **F**·d**A**

**Curl Theorem (Stokes):**
∫_S (∇×**F**) · d**A** = ∮_C **F** · d**r**

### Useful Identities

1. **Curl of gradient is zero:** ∇×(∇f) = **0**
2. **Divergence of curl is zero:** ∇·(∇×**F**) = 0
3. **Vector Laplacian:** ∇²**F** = ∇(∇·**F**) - ∇×(∇×**F**)
4. **Product rules:**
   - ∇(fg) = f∇g + g∇f
   - ∇·(f**F**) = f(∇·**F**) + **F**·∇f
   - ∇×(f**F**) = f(∇×**F**) + ∇f × **F**

---

## Derivative Order Selection

Higher derivative orders provide better accuracy but at increased computational cost.

| `der_order` | Accuracy | Cost | Use Case |
|-------------|----------|------|----------|
| 1 | O(h²) | Lowest | Quick estimates |
| 2 | O(h⁴) | Low | Default (good balance) |
| 4 | O(h⁶) | Medium | High accuracy |
| 6 | O(h⁸) | High | Very high accuracy |
| 8 | O(h¹⁰) | Highest | Maximum precision |

**Recommendation:** Start with `der_order=2` (default), increase if needed for accuracy.

---

## Integration with MML

### With Coordinate Transformations

```cpp
// Convert field from Cartesian to Spherical
Vec3Cart pos_cart(1.0, 1.0, 1.0);
Vec3Sph pos_sph = CoordTransf3D::Cart_To_Spherical(pos_cart);

// Compute gradient in spherical coordinates
Vec3Sph grad_sph = GradientSpher(scalarField, pos_sph);
```

### With Metric Tensor

```cpp
// Define metric tensor for custom coordinates
MetricTensorField<3> metric = /* ... */;

// Use general gradient/divergence
Vec3 grad = Gradient(scalarField, pos, metric);
Real div = Divergence(vectorField, pos, metric);
```

### With ODE Solvers

```cpp
// Fluid particle trajectory: dx/dt = v(x)
auto velocity_field = [](const Vec3Cart& pos) -> Vec3Cart {
    // Compute velocity using field operations
    return /* ... */;
};

// Solve particle path ODE
ODESystem pathODE(velocity_field);
auto solution = solver.solve(pathODE);
```

---

## Runnable Examples

The following demo functions in [docs_demo_field_operations.cpp](../../src/docs_demos/docs_demo_field_operations.cpp) demonstrate all field operations:

### Gradient Operations
| Function | Description |
|----------|-------------|
| `Docs_Demo_GradientCartesian()` | Cartesian gradient with derivative order comparison |
| `Docs_Demo_GradientSpherical()` | Spherical gradient (1/r potential, r·cos(θ)) |
| `Docs_Demo_GradientCylindrical()` | Cylindrical gradient (ρ² + z²) |

### Laplacian Operations
| Function | Description |
|----------|-------------|
| `Docs_Demo_LaplacianCartesian()` | Harmonic function (x² - y²) and paraboloid |
| `Docs_Demo_LaplacianSpherical()` | Coulomb potential (1/r) |
| `Docs_Demo_LaplacianCylindrical()` | 2D Coulomb potential (ln ρ) |

### Divergence Operations
| Function | Description |
|----------|-------------|
| `Docs_Demo_DivergenceCartesian()` | Expanding radial field vs incompressible rotation |
| `Docs_Demo_DivergenceSpherical()` | Radial field F = (r, 0, 0) |
| `Docs_Demo_DivergenceCylindrical()` | Radial field F = (ρ, 0, 0) |

### Curl Operations
| Function | Description |
|----------|-------------|
| `Docs_Demo_CurlCartesian()` | Rotation field and conservative field check |
| `Docs_Demo_CurlSpherical()` | Azimuthal field in spherical coordinates |
| `Docs_Demo_CurlCylindrical()` | Solid-body rotation vortex |

### Physics Applications
| Function | Description |
|----------|-------------|
| `Docs_Demo_Electric_Field()` | Point charge potential → Electric field via E = -∇φ |
| `Docs_Demo_Incompressible_Flow()` | Divergence-free velocity field (∇·v = 0) |
| `Docs_Demo_Heat_Diffusion()` | Heat equation: ∂T/∂t = κ∇²T |
| `Docs_Demo_Conservative_Field_Check()` | Test if field is conservative (curl = 0) |
| `Docs_Demo_Spherical_Harmonics()` | Spherical harmonic eigenvalue verification |

**Entry point:** `Docs_Demo_Field_operations()` runs all demonstrations.

---

## See Also

- [Derivation.md](Derivation.md) - Numerical derivatives (basis for field operations)
- [Metric_tensor.md](Metric_tensor.md) - General coordinate systems
- [Coordinate_transformations.md](Coordinate_transformations.md) - Coordinate conversions
- [Functions.md](Functions.md) - IScalarFunction, IVectorFunction interfaces
- [Fields.md](Fields.md) - Specialized field classes
- [Integration/VolumeIntegration.md](Integration/VolumeIntegration.md) - Volume integrals

---

**References:**
- Griffiths, *Introduction to Electrodynamics* (4th ed.), Cambridge, 2017
- Arfken, Weber, Harris, *Mathematical Methods for Physicists* (7th ed.), Academic Press, 2012
- Aris, *Vectors, Tensors and the Basic Equations of Fluid Mechanics*, Dover, 1989

