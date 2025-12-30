# Fields

**Sources:**  
- `mml/core/Fields.h` - Predefined scalar and vector fields
- `mml/core/FieldOperations.h` - Differential operations on fields

Scalar and vector field classes for physics and engineering applications, plus differential operations (gradient, divergence, curl, Laplacian) supporting multiple coordinate systems.

## Table of Contents
- [Overview](#overview)
- [Predefined Fields](#predefined-fields)
- [Field Operations](#field-operations)
- [Coordinate System Support](#coordinate-system-support)
- [Examples](#examples)
- [Runnable Examples](#runnable-examples)
- [See Also](#see-also)

---

## Overview

MML provides two types of field functionality:

1. **Predefined Physical Fields** (`Fields.h`): Common field configurations used in physics
   - Inverse radial potential fields (gravity, electrostatics)
   - Force fields derived from potentials

2. **Field Operations** (`FieldOperations.h`): Differential calculus on fields
   - Gradient (∇f): scalar → vector
   - Divergence (∇·F): vector → scalar
   - Curl (∇×F): vector → vector (3D only)
   - Laplacian (∇²f): scalar → scalar

---

## Predefined Fields

### Namespace: `MML::Fields`

All predefined fields are in the `MML::Fields` namespace.

### Inverse Radial Potential Fields

Scalar potential fields of the form **f(r) = k/r** (e.g., gravitational or electrostatic potential).

```cpp
// Free functions (static)
Real InverseRadialPotentialFieldCart(const VectorN<Real, 3>& x);
Real InverseRadialPotentialFieldCart(Real constant, const VectorN<Real, 3>& x);
Real InverseRadialPotentialFieldSpher(const VectorN<Real, 3>& x);
Real InverseRadialPotentialFieldSpher(Real constant, const VectorN<Real, 3>& x);
Real InverseRadialPotentialFieldCyl(const VectorN<Real, 3>& x);
Real InverseRadialPotentialFieldCyl(Real constant, const VectorN<Real, 3>& x);
```

**Class wrappers implementing `IScalarFunction<3>`:**

```cpp
class InverseRadialFieldCart : public IScalarFunction<3> {
    InverseRadialFieldCart();                    // constant = -1
    InverseRadialFieldCart(Real constant);
    Real operator()(const VectorN<Real, 3>& x) const;
};

class InverseRadialFieldSpher : public IScalarFunction<3> {
    InverseRadialFieldSpher();
    InverseRadialFieldSpher(Real constant);
    Real operator()(const VectorN<Real, 3>& x) const;
};
```

### Inverse Radial Force Fields

Vector force fields **F = -∇(k/r) = k·r/|r|³** (e.g., gravitational or Coulomb force).

```cpp
// Free functions (static)
VectorN<Real, 3> InverseRadialPotentialForceFieldCart(const VectorN<Real, 3>& x);
VectorN<Real, 3> InverseRadialPotentialForceFieldCart(Real constant, const VectorN<Real, 3>& x);
VectorN<Real, 3> InverseRadialPotentialForceFieldSph(const VectorN<Real, 3>& x);
VectorN<Real, 3> InverseRadialPotentialForceFieldSph(Real constant, const VectorN<Real, 3>& x);
```

**Class wrappers implementing `IVectorFunction<3>`:**

```cpp
class InverseRadialForceFieldCart : public IVectorFunction<3> {
    InverseRadialForceFieldCart();               // constant = -1
    InverseRadialForceFieldCart(Real constant);
    VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const;
};

class InverseRadialForceFieldSpher : public IVectorFunction<3> {
    InverseRadialForceFieldSpher();
    InverseRadialForceFieldSpher(Real constant);
    VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const;
};
```

---

## Field Operations

### Namespace: `MML::ScalarFieldOperations`

Operations on scalar fields (functions **f: ℝⁿ → ℝ**).

#### Gradient

Computes **∇f** (direction and rate of steepest increase).

```cpp
// Cartesian coordinates (any dimension)
template<int N>
VectorN<Real, N> GradientCart(const IScalarFunction<N>& f, const VectorN<Real, N>& pos);

template<int N>
VectorN<Real, N> GradientCart(const IScalarFunction<N>& f, const VectorN<Real, N>& pos, 
                              int der_order);  // 1, 2, 4, 6, or 8 point formula

// Spherical coordinates (r, θ, φ)
VectorN<Real, 3> GradientSpher(const IScalarFunction<3>& f, const VectorN<Real, 3>& pos);

// Cylindrical coordinates (r, φ, z)
VectorN<Real, 3> GradientCyl(const IScalarFunction<3>& f, const VectorN<Real, 3>& pos);

// General curvilinear coordinates
template<int N>
VectorN<Real, N> Gradient(IScalarFunction<N>& f, const VectorN<Real, N>& pos, 
                          const MetricTensorField<N>& metric);
```

#### Laplacian

Computes **∇²f** (divergence of gradient, measures local deviation from average).

```cpp
// Cartesian coordinates
template<int N>
Real LaplacianCart(const IScalarFunction<N>& f, const VectorN<Real, N>& pos);

// Spherical coordinates
Real LaplacianSpher(const IScalarFunction<3>& f, const VectorN<Real, 3>& pos);

// Cylindrical coordinates
Real LaplacianCyl(const IScalarFunction<3>& f, const VectorN<Real, 3>& pos);
```

### Namespace: `MML::VectorFieldOperations`

Operations on vector fields (functions **F: ℝⁿ → ℝⁿ**).

#### Divergence

Computes **∇·F** (measures "outflow" at a point).

```cpp
// Cartesian coordinates
template<int N>
Real DivCart(const IVectorFunction<N>& F, const VectorN<Real, N>& pos);

// Spherical coordinates
Real DivSpher(const IVectorFunction<3>& F, const VectorN<Real, 3>& pos);

// Cylindrical coordinates
Real DivCyl(const IVectorFunction<3>& F, const VectorN<Real, 3>& pos);
```

#### Curl

Computes **∇×F** (measures local rotation). Only defined in 3D.

```cpp
// Cartesian coordinates
VectorN<Real, 3> CurlCart(const IVectorFunction<3>& F, const VectorN<Real, 3>& pos);

// Spherical coordinates
VectorN<Real, 3> CurlSpher(const IVectorFunction<3>& F, const VectorN<Real, 3>& pos);

// Cylindrical coordinates
VectorN<Real, 3> CurlCyl(const IVectorFunction<3>& F, const VectorN<Real, 3>& pos);
```

---

## Coordinate System Support

### Cartesian (x, y, z)
- Standard orthonormal coordinate system
- Gradient: ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z)
- Laplacian: ∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²

### Spherical (r, θ, φ)
- r = radial distance, θ = polar angle (from z-axis), φ = azimuthal angle
- Gradient: ∇f = (∂f/∂r, (1/r)∂f/∂θ, (1/r·sinθ)∂f/∂φ)
- Best for: problems with spherical symmetry (gravity, point charges)

### Cylindrical (r, φ, z)
- r = radial distance from z-axis, φ = azimuthal angle, z = height
- Gradient: ∇f = (∂f/∂r, (1/r)∂f/∂φ, ∂f/∂z)
- Best for: problems with cylindrical symmetry (wires, pipes)

### General Curvilinear
- Uses metric tensor for arbitrary coordinate systems
- Automatically handles non-orthogonal coordinates

---

## Examples

### Example 1: Gravitational Field Analysis

```cpp
#include "core/Fields.h"
#include "core/FieldOperations.h"

using namespace MML;
using namespace MML::Fields;

void Example_Gravity_Field()
{
    // Gravitational potential: Φ = -GM/r
    Real G = 6.674e-11;  // Gravitational constant
    Real M = 5.972e24;   // Earth mass
    InverseRadialFieldCart potential(-G * M);
    
    // Point in space (1000 km altitude)
    VectorN<Real, 3> pos{7.371e6, 0.0, 0.0};
    
    // Potential energy at this point
    Real phi = potential(pos);
    
    // Gravitational acceleration = -∇Φ
    auto gradient = ScalarFieldOperations::GradientCart<3>(potential, pos);
    VectorN<Real, 3> g = (-1.0) * gradient;
    
    std::cout << "Potential: " << phi << " J/kg\n";
    std::cout << "Acceleration: " << g.NormL2() << " m/s²\n";
}
```

### Example 2: Verify Maxwell's Equations

```cpp
void Example_Divergence_Curl()
{
    // Electric field of point charge: E = q/(4πε₀) * r/|r|³
    Real q = 1.6e-19;  // Electron charge
    Real eps0 = 8.854e-12;
    Real k = q / (4 * Constants::PI * eps0);
    
    InverseRadialForceFieldCart E_field(k);
    
    VectorN<Real, 3> pos{1.0, 0.5, 0.2};
    
    // Divergence of E (should be zero except at origin)
    Real div_E = VectorFieldOperations::DivCart<3>(E_field, pos);
    std::cout << "∇·E = " << div_E << " (expect ~0 away from charge)\n";
    
    // Curl of E (should be zero for electrostatic field)
    auto curl_E = VectorFieldOperations::CurlCart(E_field, pos);
    std::cout << "∇×E = " << curl_E << " (expect ~0 for static field)\n";
}
```

### Example 3: Harmonic Functions

```cpp
void Example_Laplacian()
{
    // Harmonic function: Φ(x,y,z) = 1/r satisfies ∇²Φ = 0
    InverseRadialFieldCart harmonic(1.0);
    
    VectorN<Real, 3> pos{2.0, 1.0, 1.0};
    
    Real laplacian = ScalarFieldOperations::LaplacianCart<3>(harmonic, pos);
    std::cout << "∇²(1/r) = " << laplacian << " (expect ~0)\n";
}
```

### Example 4: Coordinate System Comparison

```cpp
void Example_Coordinate_Comparison()
{
    // Same physical field, different coordinates
    VectorN<Real, 3> cart_pos{1.0, 1.0, 1.0};
    VectorN<Real, 3> spher_pos{sqrt(3.0), acos(1/sqrt(3.0)), Constants::PI/4};  // Same point
    
    // Define potential as function
    ScalarFunction<3> phi_cart([](const VectorN<Real, 3>& x) {
        return 1.0 / x.NormL2();
    });
    
    ScalarFunction<3> phi_spher([](const VectorN<Real, 3>& x) {
        return 1.0 / x[0];  // x[0] = r in spherical
    });
    
    // Compare gradients (should give same physical result after transformation)
    auto grad_cart = ScalarFieldOperations::GradientCart<3>(phi_cart, cart_pos);
    auto grad_spher = ScalarFieldOperations::GradientSpher(phi_spher, spher_pos);
    
    std::cout << "Gradient (Cartesian components): " << grad_cart << "\n";
    std::cout << "Gradient (Spherical components): " << grad_spher << "\n";
}
```

---

## Runnable Examples

### Available Demos

| Example | Source File | Description |
|---------|-------------|-------------|
| Field Operations Demo | [docs_demo_field_operations.cpp](../../src/docs_demos/docs_demo_field_operations.cpp) | Comprehensive field operation demonstrations |

### Demo Functions

**Predefined Physical Fields:**
- `Docs_Demo_Predefined_Potential_Fields()` - Inverse radial potential fields (Cartesian, spherical, cylindrical)
- `Docs_Demo_Predefined_Force_Fields()` - Inverse square force fields, Coulomb force
- `Docs_Demo_Predefined_Fields_With_Operations()` - Verify F = -∇φ, Laplace equation

**Gradient Demos:**
- `Docs_Demo_GradientCartesian()` - Gradient in Cartesian coordinates (2D and 3D), derivative orders
- `Docs_Demo_GradientSpherical()` - Gradient in spherical coordinates (r, θ, φ)
- `Docs_Demo_GradientCylindrical()` - Gradient in cylindrical coordinates (ρ, φ, z)

**Laplacian Demos:**
- `Docs_Demo_LaplacianCartesian()` - Laplacian of harmonic functions, paraboloid
- `Docs_Demo_LaplacianSpherical()` - Laplacian of 1/r (Coulomb potential)
- `Docs_Demo_LaplacianCylindrical()` - Laplacian of ln(ρ) (2D Coulomb)

**Divergence Demos:**
- `Docs_Demo_DivergenceCartesian()` - Expanding radial fields, incompressible flow
- `Docs_Demo_DivergenceSpherical()` - Radial field in spherical coordinates
- `Docs_Demo_DivergenceCylindrical()` - Radial field in cylindrical coordinates

**Curl Demos:**
- `Docs_Demo_CurlCartesian()` - Rotation fields, conservative field verification
- `Docs_Demo_CurlSpherical()` - Azimuthal field curl
- `Docs_Demo_CurlCylindrical()` - Solid-body rotation vorticity

**Physics Applications:**
- `Docs_Demo_Electric_Field()` - Electric field from point charge potential
- `Docs_Demo_Incompressible_Flow()` - Vortex flow, divergence = 0
- `Docs_Demo_Heat_Diffusion()` - Heat equation with Laplacian
- `Docs_Demo_Conservative_Field_Check()` - Test if field is conservative (curl = 0)
- `Docs_Demo_Spherical_Harmonics()` - Spherical harmonic eigenvalue verification

### Sample Output

```cpp
// Predefined potential field
Fields::InverseRadialFieldCart potential(-GM);  // Gravitational potential
Vec3Cart earthSurface(6.371e6, 0.0, 0.0);
Real phi = potential(earthSurface);  // ~-6.26×10⁷ J/kg

// Verify F = -∇φ
InverseRadialForceFieldCart force(k);
Vec3Cart gradPhi = GradientCart<3>(potential, pos);
Vec3Cart F = force(pos);
// gradPhi and F should match (up to sign)

// Conservative field check
Vec3Cart curl = VectorFieldOperations::CurlCart(field, pos);
if (curl.NormL2() < 1e-8) {
    // Field is conservative (curl ≈ 0)
}
```

### Building and Running

```bash
# Build
cd build && cmake --build .

# Run all docs demos
./src/docs_demos/Debug/MML_DocsApp.exe  # Windows
./src/docs_demos/MML_DocsApp             # Linux/Mac
```

---

## See Also

- [Coordinate_transformations.md](Coordinate_transformations.md) - Transform between coordinate systems
- [Metric_tensor.md](Metric_tensor.md) - Metric tensors for curvilinear coordinates
- [Vector_field_operations.md](Vector_field_operations.md) - Additional vector field operations
- [Derivation.md](Derivation.md) - Numerical differentiation used internally
- [Functions.md](Functions.md) - `IScalarFunction` and `IVectorFunction` interfaces

---

**Mathematical References:**
- ∇f = gradient (scalar → vector)
- ∇·F = divergence (vector → scalar)
- ∇×F = curl (vector → vector, 3D only)
- ∇²f = Laplacian = ∇·(∇f) (scalar → scalar)

**Physical Applications:**
- Electrostatics: E = -∇Φ, ∇·E = ρ/ε₀
- Magnetostatics: ∇·B = 0, ∇×B = μ₀J
- Gravity: g = -∇Φ, ∇²Φ = 4πGρ
- Heat equation: ∂T/∂t = α∇²T
