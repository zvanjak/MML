# MatrixNM\<T, N, M\> - Fixed-Size Matrix

**File**: `mml/base/MatrixNM.h`

Fixed-size matrix with compile-time dimensions. Template parameters: `T` = element type, `N` = rows, `M` = columns.

## Quick Reference

```cpp
MatrixNM<double, 3, 3> A;                    // 3×3 zero matrix
MatrixNM<double, 2, 3> B(2.5);               // 2×3 filled with 2.5
MatrixNM<double, 3, 3> I = MatrixNM<double,3,3>::GetUnitMatrix();  // Identity
```

## Why Use MatrixNM?

- ✅ **Stack allocation** - No heap memory, faster
- ✅ **Compile-time size** - Catches dimension errors at compile time
- ✅ **Zero overhead** - Optimizes to direct array access
- ✅ **Perfect for small matrices** - 2×2, 3×3, 4×4 transformations
- ✅ **Type-safe operations** - Dimension compatibility checked at compile time

**Use MatrixNM when**: Matrix dimensions are known at compile time and small (≤ 10×10).

**Use Matrix\<T\> when**: Dimensions are runtime-dependent or matrices are large.

### Runnable Examples

| Demo Function | Description | Location |
|--------------|-------------|----------|
| `Docs_Demo_MatrixNM_initializations()` | Construction and initialization | `src/docs_demos/docs_demo_matrixnm.cpp` |
| `Docs_Demo_MatrixNM_vector_init_operations()` | `VectorFromRow/Column`, Utils | `src/docs_demos/docs_demo_matrixnm.cpp` |
| `Docs_Demo_Basic_MatrixNM_operations()` | Arithmetic, transpose, trace | `src/docs_demos/docs_demo_matrixnm.cpp` |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/MML_DocsApp   # Linux/macOS
.\build\src\Release\MML_DocsApp.exe   # Windows
```

## Construction

```cpp
// Default and value initialization
MatrixNM<double, 3, 3> zeros;           // Zero-initialized
MatrixNM<double, 2, 3> filled(5.0);     // All values = 5.0

// From initializer list (row-wise)
MatrixNM<double, 3, 3> A{1, 2, 3,
                         4, 5, 6,
                         7, 8, 9};

// Identity matrix
auto I = MatrixNM<double, 4, 4>::GetUnitMatrix();
```

## Access & Properties

```cpp
MatrixNM<double, 3, 3> A{1,2,3,4,5,6,7,8,9};

// Dimensions (compile-time constants)
int rows = A.RowNum();    // Returns N (3)
int cols = A.ColNum();    // Returns M (3)

// Element access: (row, col)
double val = A(1, 2);     // Get element
A(0, 0) = 10.0;           // Set element

// Row/column access
auto row = A.VectorFromRow(1);     // Get row 1 as VectorN
auto col = A.VectorFromColumn(2); // Get column 2 as VectorN
```

## Arithmetic Operations

```cpp
MatrixNM<double, 2, 2> A{1,2,3,4}, B{5,6,7,8};

// Addition/Subtraction (same dimensions required)
auto sum = A + B;
auto diff = A - B;

// Scalar multiplication
auto scaled = A * 2.0;
auto scaled2 = 3.0 * A;

// Matrix multiplication (dimension-compatible)
MatrixNM<double, 2, 3> C;
MatrixNM<double, 3, 4> D;
auto product = C * D;  // Results in 2×4 matrix

// Matrix-Vector multiplication
VectorN<double, 2> v{1, 2};
auto result = A * v;
```

## Matrix Operations

```cpp
MatrixNM<double, 3, 3> A{1,2,3,4,5,6,7,8,9};

// Transpose (returns new matrix with swapped dimensions)
auto At = A.GetTranspose();  // MatrixNM<double, 3, 3>

// In-place transpose (only for square matrices)
A.Transpose();  // Modifies A in place

// Trace (square matrices only)
double tr = A.Trace();
```

## Type Aliases

Common fixed-size matrix types:

```cpp
// Square matrices
typedef MatrixNM<Real, 2, 2> Mat22;
typedef MatrixNM<Real, 3, 3> Mat33;
typedef MatrixNM<Real, 4, 4> Mat44;

// Rectangular matrices
typedef MatrixNM<Real, 2, 3> Mat23;
typedef MatrixNM<Real, 3, 4> Mat34;

// Complex matrices
typedef MatrixNM<Complex, 2, 2> Mat22Complex;
typedef MatrixNM<Complex, 3, 3> Mat33Complex;
```

## Examples

### Example 1: 3D Rotation Matrix
```cpp
#include "MML.h"
using namespace MML;

// Rotation around Z-axis
double theta = M_PI / 4;  // 45 degrees
Mat33 Rz{
    cos(theta), -sin(theta), 0,
    sin(theta),  cos(theta), 0,
    0,           0,          1
};

Vec3 v{1, 0, 0};
Vec3 rotated = Rz * v;  // Type-safe multiplication

std::cout << "Original: " << v << std::endl;
std::cout << "Rotated: " << rotated << std::endl;
```

### Example 2: 2D Transformation
```cpp
// Scaling and rotation combined
Mat22 scale{2, 0,
            0, 3};

Mat22 rot{cos(M_PI/6), -sin(M_PI/6),
          sin(M_PI/6),  cos(M_PI/6)};

Mat22 transform = scale * rot;

Vec2 point{1, 1};
Vec2 transformed = transform * point;
```

### Example 3: Compile-Time Safety
```cpp
MatrixNM<double, 2, 3> A;
MatrixNM<double, 3, 4> B;

auto C = A * B;  // ✓ Results in MatrixNM<double, 2, 4>

// auto D = B * A;  // ✗ Compile error: incompatible dimensions!
```

## Performance Notes

- **Zero allocation overhead**: Stack storage, no heap
- **Compiler optimizations**: Loop unrolling, SIMD vectorization
- **Cache-friendly**: Small matrices fit in cache
- **Best for**: Graphics transformations, small linear algebra, physics calculations

## See Also
- [Matrices.md](Matrices.md) - Matrix types overview  
- [Matrix.md](Matrix.md) - Dynamic-size matrices
- [VectorN.md](VectorN.md) - Fixed-size vectors
- [Coordinate_transformations.md](../core/Coordinate_transformations.md) - Using matrices for transformations
