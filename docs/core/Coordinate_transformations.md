# Coordinate Transformations

**Location**: `mml/core/CoordTransf.h`, `mml/core/CoordTransf/*.h`  
**Dependencies**: Derivation, Vector, Matrix, Tensor, Function

---

## Overview

The coordinate transformation system provides a comprehensive framework for converting between different coordinate systems, transforming vectors and tensors, computing basis vectors and Jacobians. The implementation supports both forward and inverse transformations with automatic differentiation for basis vector computation.

### Key Capabilities

- **2D Transformations**: Cartesian ↔ Polar, rotations, oblique coordinates
- **3D Transformations**: Cartesian ↔ Spherical ↔ Cylindrical, 3D rotations (axis, Euler angles, quaternions)
- **Special Relativity**: Lorentz transformations for Minkowski spacetime
- **Vector Transformations**: Covariant and contravariant vector components
- **Tensor Transformations**: Rank-2 through rank-5 tensor transformations
- **Automatic Jacobians**: Numerical computation using differentiation module

### Physical Applications

- Classical mechanics (different reference frames)
- Electromagnetism (spherical/cylindrical symmetry)
- General relativity (curved coordinates, metric tensors)
- Computer graphics (rotations, transformations)
- Fluid dynamics (natural coordinate systems)

---

## Quick Reference

### 2D Coordinate Systems

| From System | To System | Class | Use Case |
|-------------|-----------|-------|----------|
| Polar | Cartesian | `CoordTransfPolarToCartesian2D` | Circular symmetry |
| Cartesian | Cartesian | `CoordTransfCart2DRotation` | 2D rotations |
| Cartesian | Cartesian | `CoordTransfObliqueToCartesian2D` | Non-orthogonal axes |

### 3D Coordinate Systems

| From System | To System | Class | Use Case |
|-------------|-----------|-------|----------|
| Spherical | Cartesian | `CoordTransfSphericalToCartesian` | Spherical symmetry |
| Cartesian | Spherical | `CoordTransfCartesianToSpherical` | Inverse spherical |
| Cylindrical | Cartesian | `CoordTransfCylindricalToCartesian` | Axial symmetry |
| Cartesian | Cylindrical | `CoordTransfCartesianToCylindrical` | Inverse cylindrical |
| Cartesian | Cartesian | `CoordTransfCart3DRotationXAxis` | Rotation about X |
| Cartesian | Cartesian | `CoordTransfCart3DRotationYAxis` | Rotation about Y |
| Cartesian | Cartesian | `CoordTransfCart3DRotationZAxis` | Rotation about Z |
| Cartesian | Cartesian | `CoordTransfCart3DRotationQuaternion` | General rotation |
| Cartesian | Cartesian | `CoordTransfCart3DRotationEuler` | Euler angle rotation |

### Special Relativity

| From System | To System | Class | Use Case |
|-------------|-----------|-------|----------|
| Minkowski | Minkowski | `CoordTransfLorentzXAxis` | Boost along X-axis |

---

## Core Interfaces

### ICoordTransf

```cpp
template<typename VectorFrom, typename VectorTo, int N>
class ICoordTransf : public IVectorFunction<N>
{
public:
    virtual VectorTo transf(const VectorFrom& in) const = 0;
    virtual const IScalarFunction<N>& coordTransfFunc(int i) const = 0;
};
```

**Purpose**: Base interface for coordinate transformations.

**Key Methods**:
- `transf()`: Transform position from source to target coordinates
- `coordTransfFunc(i)`: Get i-th component transformation function

### ICoordTransfWithInverse

```cpp
template<typename VectorFrom, typename VectorTo, int N>
class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
{
public:
    virtual VectorFrom transfInverse(const VectorTo& in) const = 0;
    virtual const IScalarFunction<N>& inverseCoordTransfFunc(int i) const = 0;
};
```

**Purpose**: Extends base interface with inverse transformation capability.

**Key Methods**:
- `transfInverse()`: Transform from target back to source coordinates
- `inverseCoordTransfFunc(i)`: Get i-th inverse transformation function

---

## Implementation Classes

### CoordTransf - Base Implementation

```cpp
template<typename VectorFrom, typename VectorTo, int N>
class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
{
public:
    // Compute Jacobian matrix
    MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& pos);
    
    // Basis vectors
    VectorTo getBasisVec(int ind, const VectorFrom& pos);
    VectorFrom getInverseBasisVec(int ind, const VectorFrom& pos);
    
    // Vector transformations
    VectorTo transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos);
    VectorFrom transfInverseVecCovariant(const VectorTo& vec, const VectorFrom& pos);
};
```

**Key Capabilities**:
- **Jacobian**: Computed numerically using `Derivation::NDer4Partial()`
- **Basis vectors**: Derivatives of coordinate transformation functions
- **Vector transformation**: Covariant and contravariant components

**Mathematical Background**:

The Jacobian matrix J relates coordinate differentials:
```
dx^i = J^i_j dq^j
```

For transformation x^i = f^i(q^j):
```
J^i_j = ∂f^i/∂q^j
```

### CoordTransfWithInverse - Extended Implementation

```cpp
template<typename VectorFrom, typename VectorTo, int N>
class CoordTransfWithInverse : 
    public virtual CoordTransf<VectorFrom, VectorTo, N>,
    public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
{
public:
    // Additional basis vectors
    VectorFrom getContravarBasisVec(int ind, const VectorTo& pos);
    VectorTo getInverseContravarBasisVec(int ind, const VectorTo& pos);
    
    // Bidirectional vector transformations
    VectorTo transfVecCovariant(const VectorFrom& vec, const VectorTo& pos);
    VectorFrom transfInverseVecContravariant(const VectorTo& vec, const VectorTo& pos);
    
    // Tensor transformations (rank 2-5)
    Tensor2<N> transfTensor2(const Tensor2<N>& tensor, const VectorFrom& pos);
    Tensor3<N> transfTensor3(const Tensor3<N>& tensor, const VectorFrom& pos);
    Tensor4<N> transfTensor4(const Tensor4<N>& tensor, const VectorFrom& pos);
    Tensor5<N> transfTensor5(const Tensor5<N>& tensor, const VectorFrom& pos);
};
```

**Tensor Transformation Formula**:

For a rank-2 tensor T^{ij} (contravariant):
```
T'^{ij} = (∂x'^i/∂x^k)(∂x'^j/∂x^l) T^{kl}
```

For mixed tensors, use appropriate combinations of forward/inverse Jacobians.

---

## 2D Coordinate Transformations

### Polar ↔ Cartesian

```cpp
class CoordTransfPolarToCartesian2D : 
    public CoordTransfWithInverse<Vector2Polar, Vector2Cartesian, 2>
{
public:
    Vector2Cartesian transf(const Vector2Polar& q) const;
    Vector2Polar transfInverse(const Vector2Cartesian& q) const;
};
```

**Transformation Formulas**:

Forward (Polar → Cartesian):
```
x = r cos(φ)
y = r sin(φ)
```

Inverse (Cartesian → Polar):
```
r = √(x² + y²)
φ = atan2(y, x)
```

**Example 1: Basic Polar Transformation**

```cpp
#include "core/CoordTransf/CoordTransf2D.h"

CoordTransfPolarToCartesian2D polarToCart;

// Transform from polar to Cartesian
Vector2Polar polar{5.0, Constants::PI/4};  // r=5, φ=45°
Vector2Cartesian cart = polarToCart.transf(polar);
// Result: cart ≈ (3.536, 3.536)

// Inverse transformation
Vector2Cartesian point{3.0, 4.0};
Vector2Polar back = polarToCart.transfInverse(point);
// Result: back = (5.0, 0.927) where 0.927 rad ≈ 53.13°
```

### 2D Rotation

```cpp
class CoordTransfCart2DRotation : 
    public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
{
public:
    CoordTransfCart2DRotation(Real angle);
};
```

**Rotation Matrix**:
```
R(θ) = [ cos(θ)  -sin(θ) ]
       [ sin(θ)   cos(θ) ]
```

**Example 2: 2D Rotation**

```cpp
// Rotate by 90 degrees counterclockwise
CoordTransfCart2DRotation rotate90(Constants::PI/2);

Vector2Cartesian v{1.0, 0.0};
Vector2Cartesian rotated = rotate90.transf(v);
// Result: rotated ≈ (0, 1)

// Inverse rotation (clockwise)
Vector2Cartesian back = rotate90.transfInverse(rotated);
// Result: back ≈ (1, 0)
```

---

## 3D Coordinate Transformations

### Spherical ↔ Cartesian

```cpp
class CoordTransfSphericalToCartesian : 
    public CoordTransfWithInverse<Vector3Spherical, Vector3Cartesian, 3>
{
public:
    Vector3Cartesian transf(const Vector3Spherical& q) const;
    Vector3Spherical transfInverse(const Vector3Cartesian& q) const;
    
    // Explicit basis vectors (faster than numerical)
    Vector3Cartesian getBasisVec(int ind, const Vector3Spherical& pos);
    Vector3Cartesian getUnitBasisVec(int ind, const Vector3Spherical& pos);
};

// Global instances
static CoordTransfSphericalToCartesian CoordTransfSpherToCart;
static CoordTransfCartesianToSpherical CoordTransfCartToSpher;
```

**Convention**: Mathematics/ISO 31-11 (r, θ, φ)
- r: radial distance from origin (r ≥ 0)
- θ (theta): polar angle (inclination) from +z axis, θ ∈ [0, π]
- φ (phi): azimuthal angle in xy-plane from +x axis, φ ∈ [0, 2π) or (-π, π]

**Transformation Formulas**:

Forward (Spherical → Cartesian):
```
x = r sin(θ) cos(φ)
y = r sin(θ) sin(φ)
z = r cos(θ)
```

Inverse (Cartesian → Spherical):
```
r = √(x² + y² + z²)
θ = arccos(z / r)
φ = atan2(y, x)
```

**Basis Vectors**:

Unit basis vectors at position (r, θ, φ):
```
ê_r = (sin(θ)cos(φ), sin(θ)sin(φ), cos(θ))
ê_θ = (cos(θ)cos(φ), cos(θ)sin(φ), -sin(θ))
ê_φ = (-sin(φ), cos(φ), 0)
```

**Example 3: Spherical Coordinates**

```cpp
#include "core/CoordTransf/CoordTransfSpherical.h"

// Transform point from spherical to Cartesian
Vector3Spherical sph{5.0, Constants::PI/4, Constants::PI/3};  // r=5, θ=45°, φ=60°
Vector3Cartesian cart = CoordTransfSpherToCart.transf(sph);
// Result: cart ≈ (1.768, 3.062, 3.536)

// Inverse: Cartesian to spherical
Vector3Cartesian point{1.0, 1.0, 1.0};
Vector3Spherical sphBack = CoordTransfCartToSpher.transf(point);
// Result: sphBack ≈ (1.732, 0.955, 0.785)  [r=√3, θ≈54.7°, φ=45°]

// Get spherical basis vectors at a point
Vector3Spherical pos{5.0, Constants::PI/4, Constants::PI/6};
Vector3Cartesian e_r = CoordTransfSpherToCart.getUnitBasisVec(0, pos);
Vector3Cartesian e_theta = CoordTransfSpherToCart.getUnitBasisVec(1, pos);
Vector3Cartesian e_phi = CoordTransfSpherToCart.getUnitBasisVec(2, pos);
```

### Cylindrical ↔ Cartesian

```cpp
class CoordTransfCylindricalToCartesian : 
    public CoordTransfWithInverse<Vector3Cylindrical, Vector3Cartesian, 3>
{
public:
    Vector3Cartesian transf(const Vector3Cylindrical& q) const;
    Vector3Cylindrical transfInverse(const Vector3Cartesian& q) const;
};

// Global instances
static CoordTransfCartesianToCylindrical CoordTransfCartToCyl;
static CoordTransfCylindricalToCartesian CoordTransfCylToCart;
```

**Convention**: (ρ, φ, z)
- ρ (rho): distance from z-axis (ρ ≥ 0)
- φ (phi): azimuthal angle from +x axis, φ ∈ [0, 2π) or (-π, π]
- z: height along z-axis

**Transformation Formulas**:

Forward (Cylindrical → Cartesian):
```
x = ρ cos(φ)
y = ρ sin(φ)
z = z
```

Inverse (Cartesian → Cylindrical):
```
ρ = √(x² + y²)
φ = atan2(y, x)
z = z
```

**Example 4: Cylindrical Coordinates**

```cpp
#include "core/CoordTransf/CoordTransfCylindrical.h"

// Transform from cylindrical to Cartesian
Vector3Cylindrical cyl{3.0, Constants::PI/6, 5.0};  // ρ=3, φ=30°, z=5
Vector3Cartesian cart = CoordTransfCylToCart.transf(cyl);
// Result: cart ≈ (2.598, 1.5, 5.0)

// Inverse transformation
Vector3Cartesian point{4.0, 3.0, 2.0};
Vector3Cylindrical cylBack = CoordTransfCartToCyl.transf(point);
// Result: cylBack ≈ (5.0, 0.644, 2.0)  [ρ=5, φ≈36.87°, z=2]
```

---

## 3D Rotations

### Rotation About Principal Axes

```cpp
// Rotation about X-axis
class CoordTransfCart3DRotationXAxis : 
    public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
{
public:
    CoordTransfCart3DRotationXAxis(Real angle);
};

// Similarly for Y and Z axes
class CoordTransfCart3DRotationYAxis;
class CoordTransfCart3DRotationZAxis;
```

**Rotation Matrices**:

About X-axis:
```
R_x(θ) = [ 1    0        0     ]
         [ 0  cos(θ)  -sin(θ) ]
         [ 0  sin(θ)   cos(θ) ]
```

About Y-axis:
```
R_y(θ) = [  cos(θ)  0  sin(θ) ]
         [    0     1    0    ]
         [ -sin(θ)  0  cos(θ) ]
```

About Z-axis:
```
R_z(θ) = [ cos(θ)  -sin(θ)  0 ]
         [ sin(θ)   cos(θ)  0 ]
         [   0        0     1 ]
```

**Example 5: Principal Axis Rotations**

```cpp
#include "core/CoordTransf/CoordTransf3D.h"

// Rotate 90° about Z-axis
CoordTransfCart3DRotationZAxis rotZ(Constants::PI/2);
Vec3Cart v{1.0, 0.0, 0.0};
Vec3Cart rotated = rotZ.transf(v);
// Result: rotated ≈ (0, 1, 0)

// Combine rotations: first about X, then about Y
CoordTransfCart3DRotationXAxis rotX(Constants::PI/4);
CoordTransfCart3DRotationYAxis rotY(Constants::PI/6);

Vec3Cart p{1.0, 0.0, 0.0};
Vec3Cart p1 = rotX.transf(p);
Vec3Cart p2 = rotY.transf(p1);
// Result: composed rotation
```

### Euler Angle Rotations

```cpp
class CoordTransfCart3DRotationEuler : 
    public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
{
public:
    // Convention: ZXZ Euler angles (α, β, γ)
    CoordTransfCart3DRotationEuler(Real alpha, Real beta, Real gamma);
};
```

**Convention**: ZXZ Euler angles (commonly used in classical mechanics)
- Rotate by α about Z-axis
- Rotate by β about (new) X-axis  
- Rotate by γ about (new) Z-axis

**Rotation Matrix**:
```
R(α,β,γ) = R_z(α) R_x(β) R_z(γ)
```

**Example 6: Euler Angle Rotation**

```cpp
// Create Euler rotation
Real alpha = Constants::PI/4;   // 45° about Z
Real beta = Constants::PI/6;    // 30° about X
Real gamma = Constants::PI/3;   // 60° about Z
CoordTransfCart3DRotationEuler euler(alpha, beta, gamma);

Vec3Cart v{1.0, 0.0, 0.0};
Vec3Cart rotated = euler.transf(v);

// Get rotation matrix
MatrixNM<Real, 3, 3> R = euler.jacobian(v);
```

### Quaternion Rotations

```cpp
class CoordTransfCart3DRotationQuaternion : 
    public CoordTransfWithInverse<Vec3Cart, Vec3Cart, 3>
{
public:
    // Construct from quaternion
    explicit CoordTransfCart3DRotationQuaternion(const Quaternion& q);
    
    // Construct from axis-angle
    CoordTransfCart3DRotationQuaternion(const Vec3Cart& axis, Real angle);
    
    // Construct from rotation matrix
    explicit CoordTransfCart3DRotationQuaternion(const MatrixNM<Real, 3, 3>& rotMatrix);
    
    // Factory methods from Euler angles
    static CoordTransfCart3DRotationQuaternion FromEulerZYX(Real yaw, Real pitch, Real roll);
    static CoordTransfCart3DRotationQuaternion FromEulerXYZ(Real roll, Real pitch, Real yaw);
    
    // Accessors
    const Quaternion& GetQuaternion() const;
    const MatrixNM<Real, 3, 3>& GetTransformationMatrix() const;
    const MatrixNM<Real, 3, 3>& GetInverseMatrix() const;
    Vec3Cart GetRotationAxis() const;
    Real GetRotationAngle() const;
    
    // Composition and interpolation
    CoordTransfCart3DRotationQuaternion Compose(const CoordTransfCart3DRotationQuaternion& other) const;
    CoordTransfCart3DRotationQuaternion Interpolate(const CoordTransfCart3DRotationQuaternion& other, Real t) const;
};
```

**Why Quaternions?**:
- **Gimbal lock-free**: No singularities like Euler angles
- **Efficient composition**: Multiplication faster than matrices
- **Smooth interpolation**: SLERP for animation
- **Numerical stability**: Easier to renormalize

**Rotation Formula**:

For quaternion q and vector v:
```
v' = q v q*
```

where v is treated as pure imaginary quaternion (0, v).

**Example 7: Quaternion Rotation**

```cpp
#include "base/Quaternions.h"
#include "core/CoordTransf/CoordTransf3D.h"

// Create rotation: 90° about Z-axis
Vec3Cart axis{0.0, 0.0, 1.0};
Real angle = Constants::PI/2;
Quaternion q = Quaternion::FromAxisAngle(axis, angle);

CoordTransfCart3DRotationQuaternion quatRot(q);

Vec3Cart v{1.0, 0.0, 0.0};
Vec3Cart rotated = quatRot.transf(v);
// Result: rotated ≈ (0, 1, 0)

// Alternative: construct directly from axis-angle
CoordTransfCart3DRotationQuaternion quatRot2(axis, angle);

// Compose rotations by multiplying quaternions
Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart{1,0,0}, Constants::PI/4);
Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart{0,1,0}, Constants::PI/3);
Quaternion composed = q2 * q1;
CoordTransfCart3DRotationQuaternion composedRot(composed);
```

---

## Lorentz Transformations (Special Relativity)

### Lorentz Boost

```cpp
class CoordTransfLorentzXAxis : 
    public CoordTransfWithInverse<Vector4Minkowski, Vector4Minkowski, 4>
{
public:
    // velocity in units of c (speed of light): v ∈ [0, 1)
    CoordTransfLorentzXAxis(Real velocity);
};
```

**Minkowski Coordinates**: (ct, x, y, z)
- ct: time coordinate (c × time)
- x, y, z: spatial coordinates

**Lorentz Transformation** (boost along x-axis with velocity v):

```
Λ = γ [ 1    -β   0   0 ]
      [ -β    1   0   0 ]
      [ 0     0   1   0 ]
      [ 0     0   0   1 ]
```

where:
- β = v/c (velocity in units of c)
- γ = 1/√(1 - β²) (Lorentz factor)

**Transformed Coordinates**:
```
ct' = γ(ct - βx)
x'  = γ(x - βct)
y'  = y
z'  = z
```

**Physical Interpretation**:
- Time dilation: Δt' = γΔt
- Length contraction: L' = L/γ
- Relativity of simultaneity

**Example 8: Lorentz Transformation**

```cpp
#include "core/CoordTransf/CoordTransfLorentz.h"

// Observer moving at 0.6c along x-axis
Real velocity = 0.6;  // in units of c
CoordTransfLorentzXAxis lorentz(velocity);

// Event in rest frame: ct=5, x=3, y=0, z=0
Vector4Minkowski event{5.0, 3.0, 0.0, 0.0};

// Transform to moving frame
Vector4Minkowski eventPrime = lorentz.transf(event);
// Result: Lorentz-contracted and time-dilated coordinates

// Inverse transformation (boost in opposite direction)
Vector4Minkowski back = lorentz.transfInverse(eventPrime);
// Result: back ≈ event (original coordinates)

// Verify spacetime interval invariance
Real s2 = event[0]*event[0] - event[1]*event[1] - event[2]*event[2] - event[3]*event[3];
Real s2prime = eventPrime[0]*eventPrime[0] - eventPrime[1]*eventPrime[1] 
             - eventPrime[2]*eventPrime[2] - eventPrime[3]*eventPrime[3];
// s2 ≈ s2prime (invariant interval)
```

---

## Vector and Tensor Transformations

### Vector Transformation Types

**Contravariant Vectors** (like position differentials dx^i):
```cpp
VectorTo transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos);
```

Transform with forward Jacobian: v'^i = (∂x'^i/∂x^j) v^j

**Covariant Vectors** (like gradients ∂φ/∂x^i):
```cpp
VectorTo transfVecCovariant(const VectorFrom& vec, const VectorTo& pos);
```

Transform with inverse Jacobian: v'_i = (∂x^j/∂x'^i) v_j

**Example 9: Vector Transformations**

```cpp
// Vector transformation from Cartesian to spherical
Vector3Cartesian vecCart{1.0, 0.0, 0.0};
Vector3Cartesian posCart{1.0, 1.0, 1.0};

// Position in spherical coordinates
Vector3Spherical posSpher = CoordTransfCartToSpher.transf(posCart);

// Transform contravariant vector (e.g., velocity)
Vector3Spherical vecSpherContra = 
    CoordTransfCartToSpher.transfVecContravariant(vecCart, posCart);

// Transform covariant vector (e.g., gradient)
Vector3Spherical vecSpherCov = 
    CoordTransfCartToSpher.transfVecCovariant(vecCart, posSpher);
```

### Tensor Transformations

```cpp
// Rank-2 tensor
Tensor2<N> transfTensor2(const Tensor2<N>& tensor, const VectorFrom& pos);

// Higher rank tensors (rank 3-5)
Tensor3<N> transfTensor3(const Tensor3<N>& tensor, const VectorFrom& pos);
Tensor4<N> transfTensor4(const Tensor4<N>& tensor, const VectorFrom& pos);
Tensor5<N> transfTensor5(const Tensor5<N>& tensor, const VectorFrom& pos);
```

**Transformation Formula** (rank-2 mixed tensor T^i_j):
```
T'^i_j = (∂x'^i/∂x^k)(∂x^l/∂x'^j) T^k_l
```

For fully contravariant: use forward Jacobian for both indices  
For fully covariant: use inverse Jacobian for both indices  
For mixed: use appropriate Jacobian for each index

**Example 10: Tensor Transformation**

```cpp
// Create stress tensor in Cartesian coordinates
Tensor2<3> stressCart(2, 0);  // Rank (2,0) - contravariant
stressCart(0,0) = 1.0;  // σ_xx
stressCart(1,1) = 2.0;  // σ_yy
stressCart(2,2) = 3.0;  // σ_zz

Vector3Cartesian pos{1.0, 0.0, 0.0};

// Transform to spherical coordinates
Tensor2<3> stressSpher = 
    CoordTransfCartToSpher.transfTensor2(stressCart, pos);

// Now stressSpher contains σ_rr, σ_θθ, σ_φφ, etc.
```

---

## Jacobian Computation

### Automatic Jacobian Matrix

```cpp
MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& pos);
```

Computes the Jacobian matrix J^i_j = ∂f^i/∂q^j numerically using `Derivation::NDer4Partial()`.

**Applications**:
- Change of variables in integration
- Linearization of transformations
- Metric tensor computation
- Error propagation

**Example 11: Jacobian and Change of Variables**

```cpp
// Compute Jacobian for spherical to Cartesian
Vector3Spherical pos{2.0, Constants::PI/4, Constants::PI/6};
MatrixNM<Real, 3, 3> J = CoordTransfSpherToCart.jacobian(pos);

// Jacobian determinant (volume element transformation)
Real detJ = J.Det();
// For spherical: det(J) = r² sin(θ)

// Use for volume integration:
// dV_cartesian = |det(J)| dV_spherical
// dV_cartesian = r² sin(θ) dr dθ dφ
```

### Basis Vectors

```cpp
// Covariant basis vectors (tangent to coordinate curves)
VectorTo getBasisVec(int ind, const VectorFrom& pos);

// Contravariant basis vectors (reciprocal)
VectorFrom getInverseBasisVec(int ind, const VectorFrom& pos);
```

**Mathematical Definition**:

Covariant basis: **e**_i = ∂**r**/∂q^i  
Contravariant basis: **e**^i such that **e**^i · **e**_j = δ^i_j

**Example 12: Basis Vectors**

```cpp
Vector3Spherical pos{5.0, Constants::PI/4, Constants::PI/3};

// Get spherical basis vectors in Cartesian representation
Vector3Cartesian e_r = CoordTransfSpherToCart.getBasisVec(0, pos);
Vector3Cartesian e_theta = CoordTransfSpherToCart.getBasisVec(1, pos);
Vector3Cartesian e_phi = CoordTransfSpherToCart.getBasisVec(2, pos);

// These are NOT unit vectors! They scale with coordinate changes
// For unit vectors, normalize or use getUnitBasisVec()
Real norm_e_r = e_r.NormL2();      // Should be 1.0
Real norm_e_theta = e_theta.NormL2();  // Should be r
Real norm_e_phi = e_phi.NormL2();    // Should be r sin(θ)
```

---

## Integration with Other Modules

### With Field Operations

```cpp
#include "core/FieldOperations.h"
#include "core/CoordTransf/CoordTransfSpherical.h"

// Compute gradient in different coordinate systems
RealFunction func = /* some scalar field */;

Vector3Cartesian posCart{1.0, 1.0, 1.0};
Vector3Cartesian gradCart = FieldOps::Gradient(func, posCart);

// Convert position to spherical
Vector3Spherical posSpher = CoordTransfCartToSpher.transf(posCart);

// Transform gradient (covariant vector)
Vector3Spherical gradSpher = 
    CoordTransfCartToSpher.transfVecCovariant(gradCart, posSpher);
```

### With Metric Tensor

```cpp
#include "core/MetricTensor.h"

// Metric tensor automatically computed from coordinate transformation
MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> 
    metricSpher(CoordTransfSpherToCart);

Vector3Spherical pos{2.0, Constants::PI/4, Constants::PI/6};

// Get metric tensor components
MatrixNM<Real, 3, 3> g = metricSpher.getMetricTensorCovariant(pos);
// For spherical: diag(1, r², r²sin²(θ))
```

### With Curves

```cpp
#include "core/Curves.h"

// Define curve in one coordinate system
auto helix = [](Real t) -> Vector3Cartesian {
    return Vector3Cartesian{cos(t), sin(t), t};
};

// Transform to cylindrical coordinates
auto helixCyl = [&](Real t) -> Vector3Cylindrical {
    Vector3Cartesian cart = helix(t);
    return CoordTransfCartToCyl.transf(cart);
};
```

---

## Best Practices

### Choosing Coordinate Systems

**Use spherical coordinates when**:
- Problem has spherical symmetry (planets, atoms, point charges)
- Working with angular momentum
- Integrating over spheres

**Use cylindrical coordinates when**:
- Problem has axial symmetry (cylinders, wires, toroids)
- Working in 3D with rotation about an axis
- Integrating over cylinders

**Use Cartesian coordinates when**:
- Problem has rectangular symmetry
- Forces/fields aligned with coordinate axes
- No obvious symmetry

### Numerical Considerations

**Jacobian Singularities**:
- Spherical: θ = 0, π (poles)
- Cylindrical: ρ = 0 (axis)
- Polar: r = 0 (origin)

**Avoid evaluating at singularities** or use regularization.

**Angle Wrapping**:
- Angles should be properly wrapped to [0, 2π) or (-π, π]
- Use `atan2(y, x)` instead of `atan(y/x)` to avoid quadrant ambiguity

### Performance

**Pre-compute transformations** when applying to many points:
```cpp
// Good: one transformation object
CoordTransfSphericalToCartesian transf;
for (const auto& point : points) {
    auto result = transf.transf(point);
}

// Bad: creating transformation objects in loop
for (const auto& point : points) {
    CoordTransfSphericalToCartesian transf;  // Don't do this!
    auto result = transf.transf(point);
}
```

**Use explicit basis vectors** when available (e.g., spherical) instead of numerical computation.

### Validation

**Always verify**:
1. **Round-trip accuracy**: x ≈ transf_inverse(transf(x))
2. **Jacobian determinant**: Non-zero except at singularities
3. **Orthonormality**: Basis vectors orthogonal (for orthogonal coordinates)
4. **Physical constraints**: Ranges of coordinates (r ≥ 0, θ ∈ [0,π], etc.)

---

## Runnable Examples

**Source File**: [docs_demo_coord_transf.cpp](../../src/docs_demos/docs_demo_coord_transf.cpp)

### Basic Transformations

| Function | Description |
|----------|-------------|
| `Docs_Demo_CoordTransf_Transf()` | Spherical↔Cartesian, Cylindrical↔Cartesian, forward/inverse |

### Vector Transformations

| Function | Description |
|----------|-------------|
| `Docs_Demo_CoordTransf_CovarContravar_Vec_transf()` | Contravariant velocity, covariant gradient transforms |
| `Docs_Demo_CoordTransf_CovarContravarBasis()` | Covariant/contravariant basis vectors, orthogonality check |

### Tensor Transformations

| Function | Description |
|----------|-------------|
| `Docs_Demo_CoordTransf_Tensor_transf()` | Rank-2 tensor transformation Cartesian↔Spherical |

### 3D Rotations

| Function | Description |
|----------|-------------|
| `Docs_Demo_CoordTransf_3D_Rotations()` | Rotations about X, Y, Z axes, composed rotations |
| `Docs_Demo_CoordTransf_Quaternion()` | Quaternion rotations, axis-angle, composition |

### Special Relativity

| Function | Description |
|----------|-------------|
| `Docs_Demo_CoordTransf_Lorentz()` | Lorentz boost, time dilation, length contraction |

### Jacobians

| Function | Description |
|----------|-------------|
| `Docs_Demo_CoordTransf_Jacobian()` | Jacobian matrix computation, volume element |

---

## Summary

The coordinate transformation framework provides:

✅ **Comprehensive Coverage** - 2D/3D, orthogonal/oblique, special relativity  
✅ **Bidirectional Transformations** - Forward and inverse with automatic Jacobians  
✅ **Vector/Tensor Support** - Covariant and contravariant transformation rules  
✅ **Multiple Rotation Representations** - Axis rotations, Euler angles, quaternions  
✅ **Numerical Differentiation** - Automatic Jacobian and basis vector computation  
✅ **Physical Applications** - Classical mechanics, EM, GR, computer graphics  

**Related Documentation**:
- [Metric_tensor.md](Metric_tensor.md) - Metric tensors from coordinate transformations
- [Vector_field_operations.md](Vector_field_operations.md) - Field operations in different coordinates
- [Derivation.md](Derivation.md) - Numerical differentiation for Jacobians
- [../base/Vectors.md](../base/Vectors.md) - Vector types and operations
- [../base/Quaternions.md](../base/Quaternions.md) - Quaternion mathematics
