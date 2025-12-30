# Tensors - Multi-dimensional Array Types

**File**: `mml/base/Tensor.h`

Tensor classes for orders 2–5, supporting covariant and contravariant indices.

## Table of Contents
- [Overview](#overview)
- [Tensor2\<N\> - Rank-2 Tensors](#tensor2n---rank-2-tensors)
- [Higher-Order Tensors](#higher-order-tensors)
- [Tensor Operations](#tensor-operations)
- [Examples](#examples)

---

## Overview

Tensors represent multi-dimensional arrays with index variance tracking (covariant/contravariant), essential for differential geometry, relativity, and continuum mechanics.

**Available Types:**
- `Tensor2<N>` - Rank-2 tensors (2 indices, N×N components)
- `Tensor3<N>` - Rank-3 tensors (3 indices, N×N×N components)
- `Tensor4<N>` - Rank-4 tensors (4 indices, N⁴ components)
- `Tensor5<N>` - Rank-5 tensors (5 indices, N⁵ components)

**Key Features:**
- ✅ **Index variance** - Track covariant/contravariant indices
- ✅ **Coordinate transformations** - Proper transformation under coordinate changes
- ✅ **Tensor operations** - Addition, scalar multiplication, contraction
- ✅ **Compile-time sizes** - Fixed dimensions for performance

### Runnable Examples

| Demo Function | Description | Location |
|--------------|-------------|----------|
| `Docs_Demo_Tensors()` | Basic tensor operations | `src/docs_demos/docs_demo_tensors.cpp` |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/MML_DocsApp   # Linux/macOS
.\build\src\Release\MML_DocsApp.exe   # Windows
```

---

## Tensor2\<N\> - Rank-2 Tensors

### Construction
```cpp
// Specify number of covariant and contravariant indices
Tensor2<3> T_mixed(1, 1);  // 1 covariant, 1 contravariant (mixed tensor)
Tensor2<3> T_contra(0, 2); // 2 contravariant indices (contravariant tensor)
Tensor2<3> T_covar(2, 0);  // 2 covariant indices (covariant tensor)

// With initial values
Tensor2<2> metric(2, 0, {
    1, 0,
    0, 1
});  // 2D metric tensor (covariant)
```

### Index Access
```cpp
Tensor2<3> T(1, 1);

// Access components
T(0, 0) = 1.0;
T(1, 2) = 2.5;

// Checked access (throws on out-of-bounds)
Real val = T.at(0, 0);

// Query properties
int num_contra = T.NumContravar();  // Number of contravariant indices
int num_covar = T.NumCovar();       // Number of covariant indices
bool is_contra_first = T.IsContravar(0);
```

### Properties
```cpp
Tensor2<3> T(1, 1);

// Get underlying matrix representation
MatrixNM<Real, 3, 3> mat = T.GetMatrix();

// Check index variance
bool first_is_contravariant = T.IsContravar(0);
bool second_is_covariant = T.IsCovar(1);
```

---

## Higher-Order Tensors

### Tensor3\<N\> - Rank-3 Tensors
```cpp
// 3D Christoffel symbols (1 covariant, 2 contravariant)
Tensor3<3> christoffel(1, 2);

// Access components
christoffel(0, 1, 2) = 0.5;
Real gamma_012 = christoffel(0, 1, 2);
```

### Tensor4\<N\> - Rank-4 Tensors
```cpp
// Riemann curvature tensor (4 indices)
Tensor4<4> riemann(3, 1);  // Example: 3 covariant, 1 contravariant

// Component access
riemann(0, 1, 2, 3) = 1.0;
```

### Tensor5\<N\> - Rank-5 Tensors
```cpp
// Advanced applications (rare)
Tensor5<3> T5(2, 3);  // 2 covariant, 3 contravariant

// Component access
T5(0, 1, 2, 0, 1) = 2.5;
```

---

## Tensor Operations

### Arithmetic Operations
```cpp
Tensor2<3> A(1, 1), B(1, 1);

// Addition (requires matching index variance)
Tensor2<3> sum = A + B;

// Subtraction
Tensor2<3> diff = A - B;

// Scalar multiplication
Tensor2<3> scaled = A * 2.0;
Tensor2<3> scaled2 = 3.0 * A;

// Scalar division
Tensor2<3> divided = A / 2.0;
```

### Contraction
```cpp
// Contract tensor (sum over diagonal) - requires 1 contra + 1 covar
Tensor2<3> mixed(1, 1);
mixed(0,0) = 1.0;
mixed(1,1) = 2.0;
mixed(2,2) = 3.0;

Real trace = mixed.Contract();  // 1 + 2 + 3 = 6
```

### Tensor Evaluation
```cpp
// Evaluate tensor on two vectors
Tensor2<3> T(1, 1);
VectorN<Real, 3> v1{1, 0, 0};
VectorN<Real, 3> v2{0, 1, 0};

// T(v1, v2) computes T^i_j v1_i v2^j
Real result = T(v1, v2);
```

---

## Examples

### Example 1: Metric Tensor (2D)
```cpp
// 2D Euclidean metric (covariant tensor)
Tensor2<2> metric(2, 0, {
    1, 0,
    0, 1
});

std::cout << "Metric tensor:" << std::endl;
std::cout << metric << std::endl;

// Distance calculation: ds² = g_ij dx^i dx^j
VectorN<Real, 2> dx{1.0, 1.0};
Real ds_squared = metric(dx, dx);  // Should be 2 for (1,1)
```

### Example 2: Stress Tensor (3D)
```cpp
// Cauchy stress tensor (contravariant)
Tensor2<3> stress(0, 2);

// Normal stresses (diagonal)
stress(0,0) = 100.0;  // σ_xx
stress(1,1) = 50.0;   // σ_yy
stress(2,2) = 75.0;   // σ_zz

// Shear stresses (off-diagonal, symmetric)
stress(0,1) = stress(1,0) = 10.0;  // σ_xy
stress(0,2) = stress(2,0) = 5.0;   // σ_xz
stress(1,2) = stress(2,1) = 3.0;   // σ_yz

std::cout << "Stress tensor:" << std::endl;
std::cout << stress << std::endl;
```

### Example 3: Christoffel Symbols
```cpp
// Christoffel symbols for a coordinate transformation
// Γ^k_ij (1 contravariant, 2 covariant)
Tensor3<3> christoffel(2, 1);

// Set components (example: spherical coordinates)
christoffel(1, 0, 1) = 1.0;  // Γ^θ_rθ
christoffel(2, 0, 2) = 1.0;  // Γ^φ_rφ

// Access
Real gamma_r_theta_theta = christoffel(0, 1, 1);
```

### Example 4: Tensor Arithmetic
```cpp
// Create two mixed tensors
Tensor2<2> A(1, 1, {2, 0, 0, 3});
Tensor2<2> B(1, 1, {1, 0, 0, 1});

// Add tensors
auto sum = A + B;  // {3, 0, 0, 4}

// Scale
auto scaled = A * 2.0;  // {4, 0, 0, 6}

// Contract
Real trace_A = A.Contract();  // 2 + 3 = 5
```

---

## Index Variance in Transformations

When coordinates transform, tensor components transform according to:

**Contravariant indices** (superscripts, like velocities):
```
T'^i = ∂x'^i/∂x^j T^j
```

**Covariant indices** (subscripts, like gradients):
```
T'_i = ∂x^j/∂x'^i T_j
```

**Mixed tensors** transform with both rules.

MML tensors track index variance to ensure correct transformations in coordinate system operations.

---

## Common Tensor Types

### Physics and Geometry

| Tensor | Type | Indices | Description |
|--------|------|---------|-------------|
| Metric | `Tensor2` | (2, 0) | Distance measurement g_ij |
| Inverse metric | `Tensor2` | (0, 2) | Contravariant metric g^ij |
| Stress | `Tensor2` | (0, 2) or (2, 0) | Cauchy stress tensor σ^ij |
| Strain | `Tensor2` | (0, 2) or (2, 0) | Strain tensor ε^ij |
| Christoffel | `Tensor3` | (2, 1) | Connection coefficients Γ^k_ij |
| Riemann | `Tensor4` | (3, 1) or (1, 3) | Curvature tensor R^i_jkl |

---

## Performance Notes

- **Compile-time size**: Fixed dimension N for zero-overhead abstraction
- **Storage**: Contiguous arrays for cache-friendly access
- **Operations**: Inline for optimal performance
- **Index checking**: Assertions in operator(), exceptions in at()

---

## See Also
- [Coordinate_transformations.md](../core/Coordinate_transformations.md) - Transform tensors between coordinate systems
- [Metric_tensor.md](../core/Metric_tensor.md) - Metric tensor utilities
- [Matrices.md](Matrices.md) - Matrix operations (Tensor2 with matrix interface)
- [VectorN.md](VectorN.md) - Fixed-size vectors (rank-1 tensors)

