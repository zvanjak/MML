# MML Core Layer Documentation

## Overview

The Core layer provides operations on Base layer types, implementing calculus, coordinate transformations, field theory, and linear algebra solvers. This layer bridges foundational types with numerical algorithms, enabling sophisticated mathematical computations.

**Files**: `mml/core/*.h`

### Core Overview and Quick Navigation
- **Linear systems**: See [Linear_equations_solvers.md](Linear_equations_solvers.md) for Gauss–Jordan, LU, QR, SVD, and Cholesky selection and examples.
- **Eigen solvers**: See [../algorithms/Eigen_solvers.md](../algorithms/Eigen_solvers.md) for symmetric vs. general matrices and usage tips.
- **Numerical derivation**: See [Derivation.md](Derivation.md) for order selection (1/2/4/6/8), step sizing, and error handling.
- **Vector field ops**: See [Vector_field_operations.md](Vector_field_operations.md) for gradient, divergence, curl, Laplacian in Cartesian/Spherical/Cylindrical and general metrics.
- **ODE systems**: See [ODE_system.md](ODE_system.md) for system interfaces, Jacobians, and solution export helpers.
- **Coordinate transforms**: See [Coordinate_transformations.md](Coordinate_transformations.md) for covariant/contravariant transforms, Jacobians, and predefined mappings.
- **Metric tensors**: See [Metric_tensor.md](Metric_tensor.md) for Cartesian/Spherical/Cylindrical metrics and building metrics from transforms.
- **Function spaces**: See [Function_spaces.md](Function_spaces.md) to choose between real/scalar/vector/tensor/parametric abstractions.
- **Fields**: See [Fields.md](Fields.md) for common field examples used across docs.
- **Curves and surfaces**: See [Curves_and_surfaces.md](Curves_and_surfaces.md) and [Interpolated_functions.md](Interpolated_functions.md) for geometry and visualization.

## Highlights

- **Robust solvers**: Gauss–Jordan, LU, QR, SVD, and Cholesky with practical guidance in [Linear_equations_solvers.md](Linear_equations_solvers.md).
- **High-order derivatives**: Configurable finite differences (orders 1/2/4/6/8) with stability tips in [Derivation.md](Derivation.md).
- **Vector field toolkit**: Gradient, divergence, curl, Laplacian across coordinate systems in [Vector_field_operations.md](Vector_field_operations.md).
- **Coordinates + metrics**: Transformations and metric tensors, including covariant/contravariant handling in [Coordinate_transformations.md](Coordinate_transformations.md) and [Metric_tensor.md](Metric_tensor.md).
- **ODE systems**: Interfaces, Jacobians, and helpers to couple with solvers in [ODE_system.md](ODE_system.md).
- **Function abstractions**: Real/scalar/vector/tensor/parametric spaces overview in [Function_spaces.md](Function_spaces.md).
- **Geometry primitives**: Parametric curves and surfaces with visualization notes in [Curves_and_surfaces.md](Curves_and_surfaces.md).
- **Interpolation support**: Practical guidance for sampling and visualization in [Interpolated_functions.md](Interpolated_functions.md).

## Contents

1. [Differentiation](#differentiation)
   - [Numerical Derivatives](#numerical-derivatives)
   - [Derivative Types](#derivative-types)
   - [Jacobians](#jacobians)
2. [Integration](#integration)
   - [1D Integration](#1d-integration)
   - [2D/3D Integration](#2d3d-integration)
   - [Improper Integrals](#improper-integrals)
   - [Line/Surface/Volume Integrals](#linesurfacevolume-integrals)
3. [Coordinate Systems](#coordinate-systems)
   - [CoordSystem](#coordsystem)
   - [CoordTransf](#coordtransf)
   - [Reference Frames](#reference-frames)
4. [Curves and Surfaces](#curves-and-surfaces)
   - [Parametric Curves](#parametric-curves)
   - [Parametric Surfaces](#parametric-surfaces)
5. [Fields](#fields)
   - [Scalar Fields](#scalar-fields)
   - [Vector Fields](#vector-fields)
   - [Tensor Fields](#tensor-fields)
   - [Field Operations](#field-operations)
6. [Linear Algebra Solvers](#linear-algebra-solvers)
   - [LinAlgEqSolvers](#linalgeqsolvers)
   - [MatrixUtils](#matrixutils)
7. [Differential Geometry](#differential-geometry)
   - [MetricTensor](#metrictensor)
   - [Christoffel Symbols](#christoffel-symbols)

---

## Differentiation

### Numerical Derivatives

**File**: `mml/core/Derivation/DerivationBase.h`

**Purpose**: Compute numerical derivatives of functions using finite difference methods.

#### Methods

**Central Difference** (default, 2nd order accurate):
```
f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
```

**Forward Difference** (1st order):
```
f'(x) ≈ [f(x+h) - f(x)] / h
```

**Backward Difference** (1st order):
```
f'(x) ≈ [f(x) - f(x-h)] / h
```

**Higher-Order Methods** (4th, 6th order central differences available)

---

### Derivative Types

#### DerivationRealFunction

**File**: `mml/core/Derivation/DerivationRealFunction.h`

**Purpose**: Numerical derivatives of real functions f: ℝ → ℝ

The number suffix (1,2,4,6,8) indicates the **order of accuracy**, not the derivative order.

```cpp
#include "core/Derivation.h"

RealFunction f([](Real x) { return x*x*x; });

// First derivatives (different accuracy orders)
Real df_order1 = Derivation::NDer1(f, 2.0);  // O(h) accuracy, 2 function evals
Real df_order2 = Derivation::NDer2(f, 2.0);  // O(h²) accuracy, 2 function evals (central diff)
Real df_order4 = Derivation::NDer4(f, 2.0);  // O(h⁴) accuracy, 4 function evals
Real df_order6 = Derivation::NDer6(f, 2.0);  // O(h⁶) accuracy, 6 function evals (recommended)
Real df_order8 = Derivation::NDer8(f, 2.0);  // O(h⁸) accuracy, 8 function evals

// Second derivatives (direct formulas, not chained)
Real d2f_order2 = Derivation::NSecDer2(f, 2.0);  // O(h²) accuracy, 3 function evals
Real d2f_order4 = Derivation::NSecDer4(f, 2.0);  // O(h⁴) accuracy, 5 function evals

// Third derivative
Real d3f_order2 = Derivation::NThirdDer2(f, 2.0);  // O(h²) accuracy, 4 function evals

// Custom step size
Real df_custom = Derivation::NDer6(f, 2.0, 1e-5);

// Error estimate
Real error;
Real df_with_error = Derivation::NDer6(f, 2.0, &error);

// One-sided derivatives (for boundary points)
Real df_left = Derivation::NDer2Left(f, 2.0);   // Uses points to the left
Real df_right = Derivation::NDer2Right(f, 2.0); // Uses points to the right
```

#### DerivationScalarFunction

**File**: `mml/core/Derivation/DerivationScalarFunction.h`

**Purpose**: Partial derivatives of scalar functions f: ℝⁿ → ℝ

```cpp
// Scalar function: f(x,y) = x²y + y³
ScalarFunction<2> f([](const VectorN<Real, 2>& v) {
    Real x = v[0], y = v[1];
    return x*x*y + y*y*y;
});

VectorN<Real, 2> point{2.0, 3.0};

// Partial derivatives (different accuracy orders)
Real df_dx_order2 = Derivation::NDer2Partial<2>(f, 0, point);  // ∂f/∂x, O(h²)
Real df_dy_order4 = Derivation::NDer4Partial<2>(f, 1, point);  // ∂f/∂y, O(h⁴)
Real df_dx_order6 = Derivation::NDer6Partial<2>(f, 0, point);  // ∂f/∂x, O(h⁶)

// Gradient (all partial derivatives at once)
VectorN<Real, 2> grad_order2 = Derivation::NDer2PartialByAll<2>(f, point);
VectorN<Real, 2> grad_order4 = Derivation::NDer4PartialByAll<2>(f, point);
// grad = [∂f/∂x, ∂f/∂y]ᵀ = [2xy, x² + 3y²] at (2,3) = [12, 31]

// Second partial derivatives
Real d2f_dx2 = Derivation::NSecDer2Partial<2>(f, 0, 0, point);  // ∂²f/∂x²
Real d2f_dxdy = Derivation::NSecDer2Partial<2>(f, 0, 1, point); // ∂²f/∂x∂y (mixed)
Real d2f_dy2 = Derivation::NSecDer4Partial<2>(f, 1, 1, point);  // ∂²f/∂y², O(h⁴)

// With custom step size
Real df_custom = Derivation::NDer4Partial<2>(f, 0, point, 1e-4);

// With error estimate
Real error;
Real df_err = Derivation::NDer4Partial<2>(f, 0, point, &error);
```

#### DerivationVectorFunction

**File**: `mml/core/Derivation/DerivationVectorFunction.h`

**Purpose**: Derivatives of vector functions F: ℝⁿ → ℝⁿ

```cpp
// Vector function: F(x,y) = [x²-y², 2xy]ᵀ
VectorFunction<2> F([](const VectorN<Real, 2>& v) {
    Real x = v[0], y = v[1];
    return VectorN<Real, 2>{x*x - y*y, 2*x*y};
});

VectorN<Real, 2> point{2.0, 1.0};

// Single partial derivative: ∂F₀/∂x₁ (O(h⁴) accuracy)
Real dF0_dx1 = Derivation::NDer4Partial<2>(F, 0, 1, point);  // func_index, deriv_index

// Row of Jacobian: [∂F₀/∂x₀, ∂F₀/∂x₁]
VectorN<Real, 2> grad_F0 = Derivation::NDer4PartialByAll<2>(F, 0, point);

// Full Jacobian matrix (all partial derivatives)
MatrixNM<Real, 2, 2> J = Derivation::NDer4PartialAllByAll<2>(F, point);
// J = [∂F₀/∂x₀  ∂F₀/∂x₁]   [2x   -2y]
//     [∂F₁/∂x₀  ∂F₁/∂x₁] = [2y    2x]

// Different accuracy orders available:
MatrixNM<Real, 2, 2> J1 = Derivation::NDer1PartialAllByAll<2>(F, point);  // O(h)
MatrixNM<Real, 2, 2> J2 = Derivation::NDer2PartialAllByAll<2>(F, point);  // O(h²)
MatrixNM<Real, 2, 2> J6 = Derivation::NDer6PartialAllByAll<2>(F, point);  // O(h⁶)
MatrixNM<Real, 2, 2> J8 = Derivation::NDer8PartialAllByAll<2>(F, point);  // O(h⁸)

// Default derivation function pointers (set to O(h⁴) by default)
Real dF = Derivation::DeriveVecPartial<2>(F, 0, 1, point, nullptr);
MatrixNM<Real, 2, 2> Jac = Derivation::DeriveVecPartialAllByAll<2>(F, point, nullptr);

// Note: Divergence and Curl are in VectorFieldOperations namespace
// See Fields section below
```

#### DerivationParametricCurve

**File**: `mml/core/Derivation/DerivationParametricCurve.h`

**Purpose**: Numerical derivatives of parametric curves r(t): ℝ → ℝⁿ

```cpp
// Parametric curve: r(t) = [cos(t), sin(t), t]ᵀ (helix)
ParametricCurve<3> curve([](Real t) {
    return VectorN<Real, 3>{std::cos(t), std::sin(t), t};
}, 0.0, 2*Constants::PI);

Real t = Constants::PI/4;

// First derivative r'(t) - velocity vector
VectorN<Real, 3> dr_dt = Derivation::NDer4<3>(curve, t);  // O(h⁴) accuracy

// Second derivative r''(t) - acceleration
VectorN<Real, 3> d2r_dt2 = Derivation::NSecDer2<3>(curve, t);  // O(h²) accuracy

// Third derivative r'''(t)
VectorN<Real, 3> d3r_dt3 = Derivation::NThirdDer2<3>(curve, t);  // O(h²) accuracy

// Default derivative functions (configured to O(h⁴) by default)
VectorN<Real, 3> dr = Derivation::DeriveCurve<3>(curve, t, nullptr);
VectorN<Real, 3> d2r = Derivation::DeriveCurveSec<3>(curve, t, nullptr);
VectorN<Real, 3> d3r = Derivation::DeriveCurveThird<3>(curve, t, nullptr);

// Different accuracy orders for first derivative:
VectorN<Real, 3> dr_1 = Derivation::NDer1<3>(curve, t);  // O(h) accuracy
VectorN<Real, 3> dr_2 = Derivation::NDer2<3>(curve, t);  // O(h²) accuracy
VectorN<Real, 3> dr_6 = Derivation::NDer6<3>(curve, t);  // O(h⁶) accuracy
VectorN<Real, 3> dr_8 = Derivation::NDer8<3>(curve, t);  // O(h⁸) accuracy
```

**Note**: For geometric properties (tangent, curvature, torsion, Frenet frame), use the curve class methods:
```cpp
// Using ICurveCartesian3D derived class
Curves::HelixCurve helix(1.0, 0.5);  // radius, pitch

Vec3Cart tangent = helix.getTangent(t);
Real curvature = helix.getCurvature(t);
Real torsion = helix.getTorsion(t);
Vec3Cart binormal = helix.getBinormal(t);
```
See the [Parametric Curves](#parametric-curves) section for details.

#### DerivationParametricSurface

**File**: `mml/core/Derivation/DerivationParametricSurface.h`

**Purpose**: Numerical derivatives of parametric surfaces r(u,w): ℝ² → ℝⁿ

```cpp
// Parametric surface: r(u,v) = [u*cos(v), u*sin(v), u²]ᵀ (paraboloid)
ParametricSurfaceRect<3> surface([](Real u, Real v) {
    return VectorN<Real, 3>{u*std::cos(v), u*std::sin(v), u*u};
}, 0.0, 2.0, 0.0, 2*Constants::PI);

Real u = 1.0, w = Constants::PI/4;

// Partial derivatives ∂r/∂u and ∂r/∂w
VectorN<Real, 3> r_u = Derivation::NDer4_u<3>(surface, u, w);  // O(h⁴)
VectorN<Real, 3> r_w = Derivation::NDer4_w<3>(surface, u, w);  // O(h⁴)

// Different accuracy orders:
VectorN<Real, 3> r_u_2 = Derivation::NDer2_u<3>(surface, u, w);  // O(h²)
VectorN<Real, 3> r_u_6 = Derivation::NDer6_u<3>(surface, u, w);  // O(h⁶)

// Default derivative functions
VectorN<Real, 3> du = Derivation::DeriveSurface_u<3>(surface, u, w, nullptr);
VectorN<Real, 3> dw = Derivation::DeriveSurface_w<3>(surface, u, w, nullptr);
```

**Note**: For geometric properties (normal, curvature, fundamental forms), use the surface class methods:
```cpp
// Using ISurfaceCartesian derived class
Surfaces::TorusSurface torus(2.0, 1.0);  // major, minor radius

VectorN<Real, 3> normal = torus.Normal(u, w);
Real K = torus.GaussianCurvature(u, w);
Real H = torus.MeanCurvature(u, w);
```
See the [Parametric Surfaces](#parametric-surfaces) section for details.

#### DerivationTensorField

**File**: `mml/core/Derivation/DerivationTensorField.h`

**Purpose**: Numerical derivatives of tensor fields (rank-2 tensors)

```cpp
// Rank-2 tensor field T(x): ℝⁿ → ℝⁿˣⁿ
// Example: stress tensor in continuum mechanics
ITensorField2<3>& tensorField = /* ... */;

VectorN<Real, 3> point{1.0, 2.0, 3.0};

// Partial derivative of component T_ij with respect to x_k
Real dT_01_dx2 = Derivation::NDer4Partial<3>(tensorField, 0, 1, 2, point);  // ∂T₀₁/∂x₂

// Different accuracy orders:
Real dT_00_dx0_order2 = Derivation::NDer2Partial<3>(tensorField, 0, 0, 0, point);
Real dT_00_dx0_order6 = Derivation::NDer6Partial<3>(tensorField, 0, 0, 0, point);
```

**Note**: Covariant derivatives of tensor fields use the metric tensor. 
See `MetricTensorField::CovariantDerivativeContravar()` in the [MetricTensor](#metrictensor) section.

---

### Jacobians

**File**: `mml/core/Derivation/Jacobians.h`

**Purpose**: Compute Jacobian matrices for vector functions and coordinate transformations.

#### Jacobian Matrix

For F: ℝⁿ → ℝᵐ, the Jacobian is the m×n matrix:

```
J = [∂F₁/∂x₁  ∂F₁/∂x₂  ...  ∂F₁/∂xₙ]
    [∂F₂/∂x₁  ∂F₂/∂x₂  ...  ∂F₂/∂xₙ]
    [   ⋮        ⋮      ⋱     ⋮   ]
    [∂Fₘ/∂x₁  ∂Fₘ/∂x₂  ...  ∂Fₘ/∂xₙ]
```

#### Usage

```cpp
// Square vector function: F: ℝ² → ℝ²
VectorFunction<2> polarToCart([](const VectorN<Real, 2>& v) {
    Real r = v[0], theta = v[1];
    return VectorN<Real, 2>{r*std::cos(theta), r*std::sin(theta)};
});

VectorN<Real, 2> point{2.0, Constants::PI/4};  // (r, θ) = (2, π/4)

// Compute Jacobian (uses O(h⁴) accuracy internally)
MatrixNM<Real, 2, 2> J = Derivation::calcJacobian<2>(polarToCart, point);
// J = [∂x/∂r  ∂x/∂θ]   [cos(θ)  -r*sin(θ)]
//     [∂y/∂r  ∂y/∂θ] = [sin(θ)   r*cos(θ)]

// Jacobian determinant (for coordinate volume scaling)
Real detJ = J.Det();  // = r

// Non-square function: F: ℝ² → ℝ³
VectorFunctionNM<2, 3> surfaceParam([](const VectorN<Real, 2>& v) {
    Real u = v[0], w = v[1];
    return VectorN<Real, 3>{u*std::cos(w), u*std::sin(w), u*u};
});

MatrixNM<Real, 3, 2> J_surf = Derivation::calcJacobian<2, 3>(surfaceParam, point);
```

#### Alternative: Using Partial Derivatives

```cpp
// Build Jacobian manually with specific accuracy
MatrixNM<Real, 2, 2> J_manual;
for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
        J_manual(i, j) = Derivation::NDer6Partial<2>(polarToCart, i, j, point);  // O(h⁶)
```

---

## Integration

### 1D Integration

**File**: `mml/core/Integration/Integration1D.h`

**Purpose**: Numerical integration of single-variable functions.

#### Methods Available

1. **Trapezoidal Rule** (`IntegrateTrap`)
2. **Simpson's Rule** (`IntegrateSimpson`)
3. **Romberg Integration** (`IntegrateRomberg`)
4. **Gauss-Legendre Quadrature** (`IntegrateGauss10`)

#### IntegrationResult

All integration methods return `IntegrationResult` with convergence diagnostics:

```cpp
struct IntegrationResult {
    Real value;           // Computed integral value
    Real error_estimate;  // Estimated error (0 for Gauss10)
    int iterations;       // Number of refinement iterations
    bool converged;       // True if tolerance was achieved
};
```

#### Core API

```cpp
#include "core/Integration.h"

RealFunction f([](Real x) { return std::sin(x); });

Real a = 0.0, b = Constants::PI;

// Trapezoidal rule (adaptive)
auto result_trap = IntegrateTrap(f, a, b);  // Uses default tolerance

// Simpson's rule (more accurate)
auto result_simp = IntegrateSimpson(f, a, b);

// Romberg integration (Richardson extrapolation)
auto result_romb = IntegrateRomberg(f, a, b, 1e-10);  // Custom tolerance

// Gauss-Legendre 10-point (fast, fixed order)
auto result_gauss = IntegrateGauss10(f, a, b);

// Generic interface with method selection
auto result = Integrate(f, a, b, SIMPSON, 1e-8);

// Check convergence
if (!result.converged) {
    std::cerr << "Integration did not converge, error: " << result.error_estimate << "\n";
}
Real value = result.value;
```

#### Example: Arc Length

```cpp
// Compute arc length of y = x² from x=0 to x=1
// L = ∫√(1 + (dy/dx)²) dx
RealFunction arcLen([](Real x) {
    Real dy_dx = 2*x;
    return std::sqrt(1 + dy_dx*dy_dx);
});

auto result = IntegrateSimpson(arcLen, 0.0, 1.0);
Real length = result.value;
```

---

### 2D/3D Integration

#### Integration2D

**File**: `mml/core/Integration/Integration2D.h`

**Purpose**: Double integrals over 2D regions.

```cpp
// Function: f(x,y) = x² + y²
ScalarFunction<2> f([](const VectorN<Real, 2>& v) {
    Real x = v[0], y = v[1];
    return x*x + y*y;
});

// Define Y bounds as functions of X
Real y_lower(Real x) { return 0.0; }
Real y_upper(Real x) { return 1.0 - x; }  // Triangle region

// ∫∫ f(x,y) dx dy with specified method
auto result = Integrate2D(f, SIMPSON, 0.0, 1.0, y_lower, y_upper);

if (result.converged) {
    Real volume = result.value;
}
```

#### Integration3D

**File**: `mml/core/Integration/Integration3D.h`

**Purpose**: Triple integrals over 3D regions.

```cpp
// Function: f(x,y,z) = x² + y² + z²
ScalarFunction<3> f([](const VectorN<Real, 3>& v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
});

// Define bounds (Y as function of X, Z as function of X,Y)
Real y_lo(Real x) { return 0.0; }
Real y_hi(Real x) { return 1.0; }
Real z_lo(Real x, Real y) { return 0.0; }
Real z_hi(Real x, Real y) { return 1.0; }

// Triple integral over cuboid [0,1]³
auto result = Integrate3D(f, 0.0, 1.0, y_lo, y_hi, z_lo, z_hi);
Real volume = result.value;
```

---

### Improper Integrals

**File**: `mml/core/Integration/IntegrationImproper.h`

**Purpose**: Handle infinite limits or singularities.

```cpp
RealFunction f([](Real x) { return std::exp(-x); });

// ∫[a,∞) e⁻ˣ dx (upper bound infinite)
auto result = IntegrateUpperInf(f, 0.0);  // = 1.0

// ∫(-∞,b] e^x dx (lower bound infinite)
RealFunction f_pos([](Real x) { return std::exp(x); });
auto result2 = IntegrateLowerInf(f_pos, 0.0);  // = 1.0

// ∫(-∞,∞) (1/(1+x²)) dx = π
RealFunction g([](Real x) { return 1.0/(1.0 + x*x); });
auto result3 = IntegrateInf(g);  // Uses default split at 0

// Split at custom point for asymmetric functions
auto result4 = IntegrateInfSplit(g, 0.0);  // Explicit split point

// Singularity at lower endpoint: ∫[0,1] 1/√x dx = 2
RealFunction h([](Real x) { return 1.0/std::sqrt(x); });
auto result5 = IntegrateLowerSingular(h, 0.0, 1.0);

// Singularity at upper endpoint
auto result6 = IntegrateUpperSingular(func, a, b);

// Singularities at both endpoints
auto result7 = IntegrateBothSingular(func, a, b);

// Interior singularity at point c
auto result8 = IntegrateInteriorSingular(func, a, b, c);
```

---

### Line/Surface/Volume Integrals

#### PathIntegration

**File**: `mml/core/Integration/PathIntegration.h`

**Purpose**: Line integrals along parametric curves.

##### Arc Length

```cpp
// Curve: r(t) = [cos(t), sin(t), t]ᵀ, t ∈ [0, 2π] (helix)
ParametricCurve<3> helix([](Real t) {
    return VectorN<Real, 3>{std::cos(t), std::sin(t), t};
}, 0.0, 2*Constants::PI);

// Compute arc length
Real length = PathIntegration::ParametricCurveLength(helix, 0.0, 2*Constants::PI);
```

##### Scalar Line Integral

```cpp
// Scalar field: f(x,y,z) = x² + y²
ScalarFunction<3> f([](const VectorN<Real, 3>& v) {
    return v[0]*v[0] + v[1]*v[1];
});

// ∫_C f ds (integrate scalar field along curve)
Real result = PathIntegration::LineIntegral(f, helix, 0.0, 2*Constants::PI);
```

##### Vector Line Integral (Work)

```cpp
// Vector field: F(x,y,z) = [-y, x, 0]ᵀ (rotation field)
VectorFunction<3> F([](const VectorN<Real, 3>& v) {
    return VectorN<Real, 3>{-v[1], v[0], 0.0};
});

// ∫_C F·dr (work done by force field along path)
Real work = PathIntegration::LineIntegral(F, helix, 0.0, 2*Constants::PI);
```

##### Curve Mass with Density

```cpp
// Density function ρ(t) along curve parameter
RealFunction density([](Real t) { return 1.0 + t; });

// Mass = ∫ ρ(t) ||r'(t)|| dt
Real mass = PathIntegration::ParametricCurveMass(helix, density, 0.0, 2*Constants::PI);
```

#### SurfaceIntegration

**File**: `mml/core/Integration/SurfaceIntegration.h`

**Purpose**: Surface integrals over parametric surfaces and triangulated bodies.

##### Vector Surface Integral (Flux)

```cpp
// Parametric surface: Sphere r(u,v) = [sin(u)cos(v), sin(u)sin(v), cos(u)]ᵀ
ParametricSurfaceRect<3> sphere([](Real u, Real v) {
    return VectorN<Real, 3>{std::sin(u)*std::cos(v), 
                            std::sin(u)*std::sin(v), 
                            std::cos(u)};
}, 0.0, Constants::PI, 0.0, 2*Constants::PI);

// Vector field: F(x,y,z) = [x, y, z]ᵀ (radial field)
VectorFunction<3> F([](const VectorN<Real, 3>& v) {
    return v;  // Identity/radial field
});

// ∫∫_S F·n dS (flux through surface)
// Parameters: field, surface, u_min, u_max, v_min, v_max, numU, numV
Real flux = SurfaceIntegration::SurfaceIntegral(F, sphere, 
    0.0, Constants::PI, 0.0, 2*Constants::PI, 30, 30);
// Gauss's theorem: flux = ∫∫∫_V (∇·F) dV = 3 * (4/3)πr³ = 4π for unit sphere
```

##### Integration over Triangulated Surfaces

```cpp
// For bodies defined by triangular meshes
BodyWithTriangleSurfaces solid = /* ... */;

Real total_flux = SurfaceIntegration::SurfaceIntegral(F, solid, 0.001);

// Single triangle integration with adaptive refinement
Triangle3D triangle(p1, p2, p3);
Real tri_flux = SurfaceIntegration::SurfaceIntegral(F, triangle, 0.001);
```

---

## Coordinate Systems

### CoordSystem

**File**: `mml/core/CoordSystem.h`

**Purpose**: Manage different coordinate systems and reference frames.

#### Reference Frames

##### ReferenceFrame3D

Base class for coordinate frames with parent-child relationships.

```cpp
ReferenceFrame3D* worldFrame = new ReferenceFrame3D();

// Child frame moving with constant velocity
Vector3Cartesian velocity{1.0, 0.0, 0.0};
Vector3Cartesian initial_pos{0.0, 0.0, 0.0};
InertialFrame3D* movingFrame = new InertialFrame3D(worldFrame, velocity, initial_pos);

// Get position of point in moving frame at time t
Vector3Cartesian local_point{0.0, 1.0, 0.0};
Vector3Cartesian world_point = movingFrame->GetLocalPosInParentFrameAtTime(local_point, 2.0);
// world_point = [2.0, 1.0, 0.0] at t=2.0
```

##### InertialFrame3D

Frame moving with constant velocity (no acceleration).

```cpp
// Two inertial frames
InertialFrame3D* frame1 = new InertialFrame3D(worldFrame, {10.0, 0, 0}, {0,0,0});
InertialFrame3D* frame2 = new InertialFrame3D(worldFrame, {-5.0, 0, 0}, {100,0,0});

// Relative velocity
Vector3Cartesian rel_vel = frame1->_velocity - frame2->_velocity;  // [15, 0, 0]
```

##### NonInertialFrame3D

Frame with acceleration (rotating, accelerating).

```cpp
// Rotating frame (e.g., Earth surface)
class RotatingFrame : public NonInertialFrame3D {
    Real _angular_velocity;
    VectorN<Real, 3> _axis;
public:
    virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(
        const VectorN<Real, 3>& pos, Real t) const override {
        // Apply rotation transformation
        Real angle = _angular_velocity * t;
        MatrixNM<Real, 3, 3> R = RotationMatrix(_axis, angle);
        return R * pos;
    }
};
```

---

### CoordTransf

**File**: `mml/core/CoordTransf.h` + `mml/core/CoordTransf/*.h`

**Purpose**: Transformations between coordinate systems.

#### Coordinate Transformation Types

##### Cartesian ↔ Polar (2D)

**File**: `mml/core/CoordTransf/CoordTransf2D.h`

```cpp
// Cartesian to Polar
VectorN<Real, 2> cartesian{3.0, 4.0};
VectorN<Real, 2> polar = CoordTransf::CartesianToPolar(cartesian);
// polar = [r, θ] = [5.0, 0.927...] (r = 5, θ = arctan(4/3))

// Polar to Cartesian
VectorN<Real, 2> cart2 = CoordTransf::PolarToCartesian(polar);
// cart2 = [3.0, 4.0]

// Jacobian for integration
MatrixNM<Real, 2, 2> J_polar = CoordTransf::JacobianPolarToCartesian(polar);
// J = [cos(θ)  -r*sin(θ)]
//     [sin(θ)   r*cos(θ)]
// det(J) = r
```

##### Cartesian ↔ Cylindrical (3D)

**File**: `mml/core/CoordTransf/CoordTransfCylindrical.h`

```cpp
// Cartesian to Cylindrical
VectorN<Real, 3> cart{3.0, 4.0, 5.0};
VectorN<Real, 3> cyl = CoordTransf::CartesianToCylindrical(cart);
// cyl = [ρ, φ, z] = [5.0, 0.927, 5.0]

// Cylindrical to Cartesian
VectorN<Real, 3> cart2 = CoordTransf::CylindricalToCartesian(cyl);

// Jacobian determinant = ρ
```

##### Cartesian ↔ Spherical (3D)

**File**: `mml/core/CoordTransf/CoordTransfSpherical.h`

```cpp
// Cartesian to Spherical
VectorN<Real, 3> cart{1.0, 1.0, std::sqrt(2.0)};
VectorN<Real, 3> sph = CoordTransf::CartesianToSpherical(cart);
// sph = [r, θ, φ] = [2.0, π/4, π/4]
// r = √(x²+y²+z²), θ = arccos(z/r), φ = arctan(y/x)

// Spherical to Cartesian
VectorN<Real, 3> cart2 = CoordTransf::SphericalToCartesian(sph);

// Jacobian determinant = r²sin(θ)
Real det_J = sph[0]*sph[0] * std::sin(sph[1]);

// Volume element: dV = r²sin(θ) dr dθ dφ
```

##### Lorentz Transformation

**File**: `mml/core/CoordTransf/CoordTransfLorentz.h`

**Purpose**: Special relativity transformations.

```cpp
// Event in one frame: (t, x, y, z)
VectorN<Real, 4> event1{0.0, 10.0, 0.0, 0.0};

// Transform to frame moving at velocity v
Real v = 0.6;  // 0.6c
Real c = 1.0;  // c = 1 in natural units
VectorN<Real, 4> event2 = CoordTransf::LorentzBoost(event1, v, c);

// Lorentz factor
Real gamma = 1.0 / std::sqrt(1.0 - v*v/(c*c));

// Time dilation, length contraction effects included
```

---

## Curves and Surfaces

### Parametric Curves

**File**: `mml/core/Curves.h`

**Purpose**: Parametric curve classes with geometric analysis (tangent, curvature, torsion, Frenet frame).

#### Curve Class Hierarchy

- `ICurveCartesian2D` - Base class for 2D Cartesian curves
- `ICurvePolar2D` - Base class for 2D polar curves
- `ICurveCartesian3D` - Base class for 3D space curves
- `CurveCartesian3D` - Concrete class for function-defined 3D curves

#### Built-in Curve Types

2D Curves: `Circle2DCurve`, `LogSpiralCurve`, `LemniscateCurve`, `DeltoidCurve`, `AstroidCurve`, `EllipseCurve`, `CardioidCurve`

3D Curves: `LineCurve`, `HelixCurve`, `CircleCurve3D`, `TorusKnotCurve`

#### Curve Properties (methods on curve objects)

```cpp
// Helix: r(t) = [a*cos(t), a*sin(t), b*t]ᵀ
Curves::HelixCurve helix(1.0, 0.5);  // radius=1, pitch=0.5

Real t = Constants::PI;

// Tangent vector r'(t)
Vec3Cart tangent = helix.getTangent(t);
Vec3Cart tangent_unit = helix.getTangentUnit(t);

// Normal (acceleration) vector r''(t)
Vec3Cart normal = helix.getNormal(t);
Vec3Cart principal_normal = helix.getNormalUnit(t);  // Frenet normal

// Binormal vector B = T × N
Vec3Cart binormal = helix.getBinormal(t);

// Curvature κ(t) = |r' × r''| / |r'|³
Real kappa = helix.getCurvature(t);
// For helix: κ = a/(a² + b²)

// Torsion τ(t) = (r' × r'') · r''' / |r' × r''|²
Real tau = helix.getTorsion(t);
// For helix: τ = b/(a² + b²)

// Arc length (via PathIntegration)
Real length = PathIntegration::ParametricCurveLength(helix, 0.0, 2*Constants::PI);

// Frenet-Serret frame {T, N, B}
Vector3Cartesian T, N, B;
helix.getMovingTrihedron(t, T, N, B);

// Special planes
Plane3D osculating = helix.getOsculationPlane(t);   // Contains T, N
Plane3D normal_plane = helix.getNormalPlane(t);     // Normal to T
Plane3D rectifying = helix.getRectifyingPlane(t);   // Contains T, B
```

#### Custom Curve Definition

```cpp
// Define curve via lambda
Curves::CurveCartesian3D mycurve(0.0, 2*Constants::PI,
    [](Real t) -> VectorN<Real, 3> {
        return {std::cos(t), std::sin(t), t*t};
    });

// Use all the geometric methods
Real curvature = mycurve.getCurvature(0.5);
```

---

### Parametric Surfaces

**File**: `mml/core/Surfaces.h`

**Purpose**: Parametric surface classes with geometric analysis (normal, curvature, fundamental forms).

#### Surface Class Hierarchy

- `ISurfaceCartesian` - Base class for parametric surfaces r(u,w): ℝ² → ℝ³
- Built-in surfaces: `SphereSurface`, `TorusSurface`, `ConeSurface`, `CylinderSurface`, etc.

#### Surface Properties (methods on surface objects)

```cpp
// Torus: r(u,v) = [(R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v)]
Surfaces::TorusSurface torus(2.0, 1.0);  // Major radius R=2, minor radius r=1

Real u = Constants::PI/4, w = Constants::PI/3;

// Normal vector (outward)
VectorN<Real, 3> N = torus.Normal(u, w);

// Tangent vectors
VectorN<Real, 3> tU, tW;
torus.Tangents(u, w, tU, tW);  // ∂r/∂u, ∂r/∂w

// First fundamental form (E, F, G)
Real E, F, G;
torus.GetFirstNormalFormCoefficients(u, w, E, F, G);
// E = ∂r/∂u · ∂r/∂u, F = ∂r/∂u · ∂r/∂w, G = ∂r/∂w · ∂r/∂w

// Second fundamental form (L, M, N)
Real L, M, N;
torus.GetSecondNormalFormCoefficients(u, w, L, M, N);

// Gaussian curvature K = (LN - M²) / (EG - F²)
Real K = torus.GaussianCurvature(u, w);

// Mean curvature H = (EN + GL - 2FM) / (2(EG - F²))
Real H = torus.MeanCurvature(u, w);

// Principal curvatures k₁, k₂ where K = k₁k₂ and H = (k₁+k₂)/2
Real k1, k2;
torus.PrincipalCurvatures(u, w, k1, k2);

// Principal directions
VectorN<Real, 3> dir1, dir2;
torus.PrincipalDirections(u, w, dir1, dir2);
```

#### Custom Surface Definition

```cpp
// Define surface via lambda (parameter domain [0,2π] x [0,2π])
class MySurface : public Surfaces::ISurfaceCartesian {
    Real R, r;
public:
    MySurface(Real R_, Real r_) : R(R_), r(r_) {}
    
    Real GetMinX() const override { return 0.0; }
    Real GetMaxX() const override { return 2*Constants::PI; }
    Real GetMinY() const override { return 0.0; }
    Real GetMaxY() const override { return 2*Constants::PI; }
    
    VectorN<Real, 3> operator()(Real u, Real w) const override {
        Real x = (R + r*std::cos(w)) * std::cos(u);
        Real y = (R + r*std::cos(w)) * std::sin(u);
        Real z = r * std::sin(w);
        return {x, y, z};
    }
};

MySurface mySurf(3.0, 1.0);
Real gauss = mySurf.GaussianCurvature(0.5, 0.5);
```

---

## Fields

### Scalar Fields

**File**: `mml/core/FieldOperations.h`

**Purpose**: Differential operations on scalar functions f: ℝⁿ → ℝ.

```cpp
// Temperature field: T(x,y,z) = 100*exp(-(x²+y²+z²))
ScalarFunction<3> temperature([](const VectorN<Real, 3>& r) {
    Real r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    return 100.0 * std::exp(-r2);
});

VectorN<Real, 3> point{1.0, 0.0, 0.0};

// Evaluate field
Real T = temperature(point);  // T ≈ 36.8

// Gradient in Cartesian coordinates
VectorN<Real, 3> grad_T = ScalarFieldOperations::GradientCart(temperature, point);

// Laplacian in Cartesian coordinates
Real lap_T = ScalarFieldOperations::LaplacianCart(temperature, point);

// Spherical coordinates (r, θ, φ)
Vec3Sph sph_point{1.0, Constants::PI/4, Constants::PI/4};
VectorN<Real, 3> grad_sph = ScalarFieldOperations::GradientSph(temperature, sph_point);
Real lap_sph = ScalarFieldOperations::LaplacianSph(temperature, sph_point);

// Cylindrical coordinates (r, φ, z)
Vec3Cyl cyl_point{1.0, Constants::PI/4, 0.5};
VectorN<Real, 3> grad_cyl = ScalarFieldOperations::GradientCyl(temperature, cyl_point);
Real lap_cyl = ScalarFieldOperations::LaplacianCyl(temperature, cyl_point);

// General coordinates with metric tensor
MetricTensorField<3> metric = /* ... */;
VectorN<Real, 3> grad_gen = ScalarFieldOperations::Gradient(temperature, point, metric);
```

---

### Vector Fields

**File**: `mml/core/FieldOperations.h`

**Purpose**: Differential operations on vector functions F: ℝⁿ → ℝⁿ.

```cpp
// Electric field of point charge: E = q*r/(4πε₀r³)
Real q = 1.0, epsilon0 = 1.0;
VectorFunction<3> electricField([q, epsilon0](const VectorN<Real, 3>& r) {
    Real r_mag = r.NormL2();
    if (r_mag < 1e-10) return VectorN<Real, 3>{0,0,0};  // Singularity at origin
    
    Real factor = q / (4 * Constants::PI * epsilon0 * r_mag * r_mag * r_mag);
    return VectorN<Real, 3>{factor * r[0], factor * r[1], factor * r[2]};
});

VectorN<Real, 3> point{1.0, 0.0, 0.0};

// Evaluate field
VectorN<Real, 3> E = electricField(point);

// Divergence in Cartesian coordinates
Real div_E = VectorFieldOperations::DivCart(electricField, point);
// ∇·E = ρ/ε₀ (Gauss's law)

// Curl in Cartesian coordinates (3D only)
VectorN<Real, 3> curl_E = VectorFieldOperations::CurlCart(electricField, point);
// ∇×E = 0 for electrostatic fields

// Spherical coordinates
Vec3Sph sph_point{1.0, Constants::PI/4, 0.0};
Real div_sph = VectorFieldOperations::DivSph(electricField, sph_point);
VectorN<Real, 3> curl_sph = VectorFieldOperations::CurlSph(electricField, sph_point);

// Cylindrical coordinates
Real div_cyl = VectorFieldOperations::DivCyl(electricField, cyl_point);
VectorN<Real, 3> curl_cyl = VectorFieldOperations::CurlCyl(electricField, cyl_point);
```

---

### Field Operations Summary

**File**: `mml/core/FieldOperations.h`

**Purpose**: Comprehensive differential operations on scalar and vector fields.

#### ScalarFieldOperations Namespace

```cpp
namespace ScalarFieldOperations {
    // Gradient operations (scalar → vector)
    VectorN<Real, N> GradientCart(f, pos);           // Cartesian
    VectorN<Real, 3> GradientSph(f, pos);            // Spherical (r, θ, φ)
    VectorN<Real, 3> GradientCyl(f, pos);            // Cylindrical (r, φ, z)
    VectorN<Real, N> Gradient(f, pos, metric);       // General coordinates
    
    // Laplacian operations (scalar → scalar)
    Real LaplacianCart(f, pos);                      // Cartesian: ∑ ∂²f/∂xᵢ²
    Real LaplacianSph(f, pos);                       // Spherical coordinates
    Real LaplacianCyl(f, pos);                       // Cylindrical coordinates
}
```

#### VectorFieldOperations Namespace

```cpp
namespace VectorFieldOperations {
    // Divergence operations (vector → scalar)
    Real DivCart(F, pos);                            // Cartesian: ∑ ∂Fᵢ/∂xᵢ
    Real DivSph(F, pos);                             // Spherical coordinates
    Real DivCyl(F, pos);                             // Cylindrical coordinates
    Real Divergence(F, pos, metric);                 // General coordinates
    
    // Curl operations (vector → vector, 3D only)
    VectorN<Real, 3> CurlCart(F, pos);               // Cartesian
    VectorN<Real, 3> CurlSph(F, pos);                // Spherical coordinates
    VectorN<Real, 3> CurlCyl(F, pos);                // Cylindrical coordinates
}
```

#### Example: Verify Vector Calculus Identity

```cpp
// ∇×(∇f) = 0 (curl of gradient is zero)
ScalarFunction<3> f([](const VectorN<Real, 3>& p) {
    return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
});

VectorN<Real, 3> point{1.0, 1.0, 1.0};

// Compute gradient
VectorN<Real, 3> grad_f = ScalarFieldOperations::GradientCart(f, point);

// Wrap gradient as vector function for curl
VectorFunction<3> grad_field([&f](const VectorN<Real, 3>& p) {
    return ScalarFieldOperations::GradientCart(f, p);
});

// Curl should be approximately zero
VectorN<Real, 3> curl_grad = VectorFieldOperations::CurlCart(grad_field, point);
// curl_grad ≈ [0, 0, 0]
```

---

## Linear Algebra Solvers

### LinAlgEqSolvers

**File**: `mml/core/LinAlgEqSolvers.h`

**Purpose**: Solve linear systems Ax = b using various decomposition methods.

The library provides solver classes that perform decomposition once and can solve for multiple right-hand sides efficiently.

#### Available Solver Classes

- `GaussJordanSolver<Type>` - Gauss-Jordan elimination with full pivoting
- `LUSolver<Type>` - LU decomposition with partial pivoting
- `LUSolverInPlace<Type>` - LU decomposition modifying input matrix
- `BandDiagonalSolver<Type>` - Specialized for band diagonal matrices
- `CholeskySolver<Type>` - For symmetric positive definite matrices
- `QRSolver<Type>` - QR decomposition via Householder reflections
- `SVDecompositionSolver<Type>` - Singular Value Decomposition

#### Direct Solvers

```cpp
Matrix<Real> A(n, n);  // Coefficient matrix
Vector<Real> b(n);     // Right-hand side

// LU Decomposition
LUSolver<Real> lu(A);              // Perform decomposition
Vector<Real> x = lu.Solve(b);      // Solve for single RHS
Real det = lu.det();               // Get determinant

Matrix<Real> A_inv;
lu.inverse(A_inv);                 // Compute inverse

// Gauss-Jordan (modifies matrices in place)
Matrix<Real> A_copy = A;
Matrix<Real> B_copy = Utils::ColumnMatrixFromVector(b);
GaussJordanSolver<Real>::SolveInPlace(A_copy, B_copy);
Vector<Real> x_gj = B_copy.VectorFromColumn(0);

// Cholesky (for symmetric positive definite)
CholeskySolver<Real> chol(A_spd);
Vector<Real> x_chol = chol.Solve(b);

// QR Decomposition (general matrices)
QRSolver<Real> qr(A);
Vector<Real> x_qr = qr.Solve(b);

// SVD (handles ill-conditioned systems)
SVDecompositionSolver<Real> svd(A);
Vector<Real> x_svd = svd.Solve(b);
int rank = svd.Rank();             // Numerical rank
```

#### Least Squares (Overdetermined Systems)

```cpp
// Overdetermined system (m > n): minimize ||Ax - b||₂
Matrix<Real> A(m, n);  // m > n
Vector<Real> b(m);

// QR-based least squares (recommended)
QRSolver<Real> qr(A);
Vector<Real> x_ls = qr.LeastSquaresSolve(b);

// SVD-based (handles rank deficiency)
SVDecompositionSolver<Real> svd(A);
Vector<Real> x_svd = svd.Solve(b);  // Automatically handles overdetermined
```

#### Band Diagonal Systems

```cpp
// For tridiagonal or band diagonal matrices
BandDiagonalSolver<Real> band(A_band);
Vector<Real> x = band.Solve(b);
```

---

### MatrixUtils

**File**: `mml/core/MatrixUtils.h`

**Purpose**: Utility functions for matrix analysis and special matrix construction.

```cpp
namespace Utils {
    // Matrix properties
    Type det = Det(A);                         // Determinant via LU
    int rank = Rank(A, tol);                   // Numerical rank via Gaussian elimination
    
    // Matrix classification
    bool nilp = IsNilpotent(A);               // A^k = 0 for some k
    bool unip = IsUnipotent(A);               // (A - I)^k = 0 for some k
    
    // Faddeev-Leverrier algorithm
    PolynomRealFunc charPoly;
    Real det;
    Matrix<Real> inv;
    FaddeevAlg(A, charPoly, det, inv);        // Characteristic polynomial, det, inverse
}

// Note: Matrix trace is a Matrix method: A.Trace()
// Note: Matrix inverse via: A.GetInverse() or LUSolver::inverse()
```

---

## Differential Geometry

### MetricTensor

**File**: `mml/core/MetricTensor.h`

**Purpose**: Metric tensors for curved spaces, general relativity, and coordinate transformations.

#### MetricTensorField Class

The metric tensor defines distance in curved space: `ds² = gᵢⱼ dxⁱ dxʲ`

```cpp
// Built-in metric tensor types
MetricTensorCartesian3D metricCart;     // Euclidean (I)
MetricTensorSpherical metricSpher;       // Spherical (r, θ, φ)
MetricTensorCylindrical metricCyl;       // Cylindrical (r, φ, z)
MetricTensorMinkowski metricMink;        // Minkowski spacetime (-,+,+,+)

// Metric from coordinate transformation
MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> 
    metricSpherFromCart(coordTransfSpherToCart);

// Get metric tensor components at a point
VectorN<Real, 3> pos{1.0, Constants::PI/4, Constants::PI/4};  // r, θ, φ

// Covariant metric gᵢⱼ
MatrixNM<Real, 3, 3> g = metricSpher.GetCovariantMetric(pos);

// Contravariant metric gⁱʲ (inverse)
MatrixNM<Real, 3, 3> g_inv = metricSpher.GetContravariantMetric(pos);
```

#### Christoffel Symbols

Γⁱⱼₖ measure how basis vectors change (connection coefficients).

```cpp
// Christoffel symbols of the second kind
// Γⁱⱼₖ = ½ gⁱˡ (∂ⱼgₗₖ + ∂ₖgₗⱼ - ∂ₗgⱼₖ)
Real gamma_ijk = metricSpher.GetChristoffelSymbolSecondKind(i, j, k, pos);

// Christoffel symbols of the first kind
Real gamma_ijk_1 = metricSpher.GetChristoffelSymbolFirstKind(i, j, k, pos);

// Geodesic equation: d²xⁱ/dt² + Γⁱⱼₖ (dxʲ/dt)(dxᵏ/dt) = 0
```

#### Covariant Derivatives

```cpp
// Covariant derivative of contravariant vector field
IVectorFunction<3> vecField = /* ... */;
VectorN<Real, 3> covDer = metricSpher.CovariantDerivativeContravar(vecField, j, pos);

// Component of covariant derivative
Real comp = metricSpher.CovariantDerivativeContravarComp(vecField, i, j, pos);

// Covariant derivative of covariant vector field
VectorN<Real, 3> covDerCo = metricSpher.CovariantDerivativeCovar(vecField, j, pos);
```

#### Example: Minkowski Spacetime

```cpp
MetricTensorMinkowski metricMinkowski;
VectorN<Real, 4> pos4d{0.0, 0.0, 0.0, 0.0};  // t, x, y, z

// Get Minkowski metric tensor η = diag(-1, 1, 1, 1)
Tensor2<4> metricTensor = metricMinkowski(pos4d);

// Proper time interval: ds² = -c²dt² + dx² + dy² + dz²
VectorN<Real, 4> tangent_vec = /* worldline tangent */;
Real interval = -metricTensor(tangent_vec, tangent_vec);
// interval < 0 → timelike, interval > 0 → spacelike, interval = 0 → lightlike
```

---

## Helper Functions

### FunctionHelpers

**File**: `mml/core/FunctionHelpers.h`

**Purpose**: Utilities for working with functions - composition, arithmetic, derivatives, Taylor series.

#### Taylor Series Expansion

```cpp
// Taylor polynomial approximation around point a
// f(x) ≈ f(a) + f'(a)(x-a) + f''(a)(x-a)²/2! + ...
RealFunction f([](Real x) { return std::exp(x); });

PolynomRealFunc taylor2 = TaylorSeries2(f, 0.0);  // Quadratic approximation
PolynomRealFunc taylor3 = TaylorSeries3(f, 0.0);  // Cubic approximation

Real approx = taylor3(0.5);  // ≈ exp(0.5)
```

#### Derivative Wrapper Classes

```cpp
// Create a function that returns the derivative of f
RealFunction f([](Real x) { return std::sin(x); });

RealFuncDerived1 fprime1(f);         // 1st order accuracy
RealFuncDerived2 fprime2(f);         // 2nd order accuracy  
RealFuncDerived4 fprime4(f);         // 4th order accuracy
RealFuncDerived6 fprime6(f);         // 6th order accuracy (default)
RealFuncDerived8 fprime8(f);         // 8th order accuracy

Real df = fprime6(Constants::PI/4);  // cos(π/4) ≈ 0.707

// Second derivatives
RealFuncSecondDerived2 f2prime2(f);  // 2nd order accuracy
RealFuncSecondDerived4 f2prime4(f);  // 4th order accuracy
Real d2f = f2prime4(0.0);            // -sin(0) = 0
```

#### Function Arithmetic Classes

```cpp
RealFunction f([](Real x) { return std::sin(x); });
RealFunction g([](Real x) { return std::cos(x); });

// Arithmetic operations (classes store references!)
RealFuncSum sum(f, g);           // h(x) = f(x) + g(x)
RealFuncDiff diff(f, g);         // h(x) = f(x) - g(x)
RealFuncProduct prod(f, g);      // h(x) = f(x) * g(x)
RealFuncQuotient quot(f, g);     // h(x) = f(x) / g(x)

// Composition: h(x) = f(g(x))
RealFuncCompose composed(f, g);  // sin(cos(x))

// Scalar operations
RealFuncScale scaled(f, 2.0);    // h(x) = 2 * f(x)
RealFuncShift shifted(f, 1.0);   // h(x) = f(x) + 1
RealFuncNegate neg(f);           // h(x) = -f(x)
RealFuncAbs absf(f);             // h(x) = |f(x)|
RealFuncPow powf(f, 2.0);        // h(x) = f(x)²

// Comparison helpers (for norms)
RealFuncAbsDiff absdiff(f, g);   // h(x) = |f(x) - g(x)|
RealFuncDiffSqr diffsqr(f, g);   // h(x) = (f(x) - g(x))²
```

---

## Integration with Other Layers

### With Base Layer

```cpp
// Use base types
Vector<Real> v = /* ... */;
Matrix<Real> A = /* ... */;
RealFunction f = /* ... */;

// Apply core operations
IntegrationResult result = IntegrateSimpson(f, 0.0, 1.0);

LUSolver<Real> lu(A);
Vector<Real> x = lu.Solve(v);
```

### With Algorithms Layer

```cpp
// Core provides operations, algorithms use them
// Example: ODE solver needs derivatives
ODESystem system = /* ... */;
// Solver internally uses Derivation for Jacobian computation
```

---

## Performance Considerations

- **Numerical derivatives**: Accuracy vs. computation trade-off (step size h)
- **Integration**: Adaptive methods reduce function evaluations
- **Iterative solvers**: Choose based on matrix properties (SPD → CG)
- **Coordinate transformations**: Precompute Jacobians when possible

---

## Error Handling

```cpp
try {
    LUSolver<Real> lu(A);
    Vector<Real> x = lu.Solve(b);
} catch (const std::exception& e) {
    std::cerr << "Error solving system: " << e.what() << "\n";
}

// Check IntegrationResult for convergence
IntegrationResult result = IntegrateRomberg(f, a, b);
if (!result.converged) {
    std::cerr << "Integration did not converge after " 
              << result.iterations << " iterations\n";
}
```

---

## Testing

Run Core layer tests:
```bash
cd build
.\tests\Debug\MML_Tests.exe "[core]"
.\tests\Debug\MML_Tests.exe "[derivation]"
.\tests\Debug\MML_Tests.exe "[integration]"
```

---

## Summary

The Core layer provides:
- **Differentiation**: 7 derivative types (Real, Scalar, Vector, Parametric, Tensor)
- **Integration**: 7 integration types (1D/2D/3D, improper, line, surface, volume)
- **Coordinate Systems**: 5 transformation types + reference frames
- **Curves/Surfaces**: Geometric analysis (curvature, torsion, area)
- **Fields**: Scalar, vector, tensor fields + operations
- **Linear Algebra**: Direct and iterative solvers, eigenvalue problems
- **Differential Geometry**: Metric tensors, Christoffel symbols

All operations:
✅ Numerically stable  
✅ Template-based (type-generic)  
✅ Configurable accuracy  
✅ Well-tested  

**Previous Layer**: [Base Layer Documentation](README_Base.md)  
**Next Layer**: [Algorithms Layer Documentation](README_Algorithms.md)

---

*Last Updated: June 21, 2025*  
*MML Version: 1.0*
