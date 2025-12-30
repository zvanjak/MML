# VectorN\<T, N\> - Fixed-Size Vector

**File**: `mml/base/VectorN.h`

Fixed-size vector with compile-time dimensions for optimal performance. Template parameters: `T` = element type, `N` = dimension.

## Runnable Examples

All code examples in this document are available as runnable code in `src/docs_demos/docs_demo_vectorn.cpp`. Build and run `MML_DocsApp` to see them in action.

| Demo Function | Description |
|---------------|-------------|
| `Demo_VectorN_Construction()` | Creation patterns, initializer lists, unit vectors |
| `Demo_VectorN_Access()` | Element access, at(), size(), IsNullVec(), clear() |
| `Demo_VectorN_Arithmetic()` | +, -, *, /, unary minus, compound assignment |
| `Demo_VectorN_Operations()` | Norms, normalization, dot product, angle, distance |
| `Demo_VectorN_Equality()` | ==, !=, IsEqualTo(), AreEqual() with tolerance |
| `Demo_VectorN_IO()` | Print(), PrintLine(), to_string(), formatting |
| `Demo_VectorN_TypeAliases()` | Vec2/3/4, Vec2D/3D/4D, Vec2C/3C/4C |
| `Demo_VectorN_Complex()` | Complex vector operations |
| `Demo_VectorN_Physics()` | Physics simulation example |

## Table of Contents
- [Overview](#overview)
- [Construction](#construction)
- [Access & Properties](#access--properties)
- [Arithmetic Operations](#arithmetic-operations)
- [Vector Operations](#vector-operations)
- [Specialized Types](#specialized-types)
- [Examples](#examples)

---

## Overview

`VectorN<T, N>` is MML's fixed-size vector class, providing:
- ✅ **Compile-time sizing** - Dimension known at compile time
- ✅ **Stack allocation** - No heap allocation overhead
- ✅ **Type-generic** - Works with `Real`, `Complex`, custom types
- ✅ **Zero overhead** - Optimizes to direct array access
- ✅ **Perfect for small vectors** - 2D, 3D, 4D coordinates
- ✅ **Interoperable** - Works seamlessly with specialized types

**When to use**: Small vectors with fixed dimensions (coordinates, physics vectors, color values, small linear algebra problems).

**When to use `Vector<T>` instead**: Variable-size data, I/O operations, large vectors with runtime dimensions.

---

## Construction

### Default and Value Construction
```cpp
VectorN<double, 3> v1;           // Zero-initialized: [0, 0, 0]
VectorN<double, 3> v2(5.0);      // All values = 5.0: [5, 5, 5]
```

### From Initializer List
```cpp
VectorN<double, 3> v{1.0, 2.0, 3.0};  // [1, 2, 3]
```

### From std::vector
```cpp
std::vector<double> data = {1, 2, 3, 4, 5};
VectorN<double, 3> v(data);  // Takes first 3 elements: [1, 2, 3]
```

### From C Array
```cpp
double data[] = {1.0, 2.0, 3.0};
VectorN<double, 3> v(data);  // [1, 2, 3]
```

### Unit Vectors
```cpp
auto e0 = VectorN<double, 3>::GetUnitVector(0);  // [1, 0, 0]
auto e1 = VectorN<double, 3>::GetUnitVector(1);  // [0, 1, 0]
auto e2 = VectorN<double, 3>::GetUnitVector(2);  // [0, 0, 1]
```

---

## Access & Properties

### Element Access
```cpp
VectorN<double, 3> v{1.5, 2.5, 3.5};

// Fast access (no bounds checking)
double x = v[0];     // 1.5
v[1] = 10.0;

// Safe access (with bounds checking)
double y = v.at(1);  // 10.0
v.at(2) = 20.0;
// v.at(5) = 5.0;    // Throws VectorDimensionError
```

### Size and State
```cpp
VectorN<double, 3> v{1, 2, 3};

int n = v.size();           // Always returns N (3)
bool null = v.IsNullVec();  // True if all elements are zero

v.clear();  // Set all elements to zero
```

---

## Arithmetic Operations

### Unary Minus
```cpp
VectorN<double, 3> v{1, -2, 3};
VectorN<double, 3> neg = -v;  // [-1, 2, -3]
```

### Vector Addition/Subtraction
```cpp
VectorN<double, 3> a{1, 2, 3};
VectorN<double, 3> b{4, 5, 6};

VectorN<double, 3> sum = a + b;   // [5, 7, 9]
VectorN<double, 3> diff = a - b;  // [-3, -3, -3]

a += b;  // a = [5, 7, 9]
a -= b;  // a = [1, 2, 3] (back to original)
```

### Scalar Multiplication/Division
```cpp
VectorN<double, 3> v{1, 2, 3};

auto scaled = v * 2.0;      // [2, 4, 6]
auto scaled2 = 3.0 * v;     // [3, 6, 9]
auto divided = v / 2.0;     // [0.5, 1, 1.5]

v *= 2.0;   // v = [2, 4, 6]
v /= 2.0;   // v = [1, 2, 3]
```

---

## Vector Operations

### Norms
```cpp
VectorN<double, 2> v{3, 4};

// L2 norm (Euclidean length)
double len = v.NormL2();     // sqrt(3² + 4²) = 5.0

// L1 norm (Manhattan distance)
double l1 = v.NormL1();      // |3| + |4| = 7.0

// L-infinity norm (maximum absolute value)
double linf = v.NormLInf();  // max(|3|, |4|) = 4.0
```

### Normalization
```cpp
VectorN<double, 2> v{3, 4};

// Create normalized copy
auto unit = v.Normalized();  // [0.6, 0.8]

// Original vector unchanged
std::cout << v << std::endl;     // [3, 4]
std::cout << unit << std::endl;  // [0.6, 0.8]
```

---

## Comparison & Equality

### Exact Equality
```cpp
VectorN<double, 3> a{1, 2, 3};
VectorN<double, 3> b{1, 2, 3};
VectorN<double, 3> c{1, 2, 3.0001};

bool eq1 = (a == b);   // true
bool eq2 = (a != c);   // true
```

### Approximate Equality
```cpp
VectorN<double, 3> a{1.0, 2.0, 3.0};
VectorN<double, 3> b{1.0000001, 2.0000001, 3.0000001};

// Static method
bool approx1 = VectorN<double, 3>::AreEqual(a, b, 1e-6);  // true

// Member method
bool approx2 = a.IsEqualTo(b, 1e-6);  // true
bool approx3 = a.IsEqualTo(b, 1e-8);  // false
```

---

## I/O and Printing

### Console Output
```cpp
VectorN<double, 3> v{1.5, 2.718, 3.14159};

// Stream operator (default formatting)
std::cout << v << std::endl;
// Output: [1.5, 2.718, 3.14159]

// Custom formatting: Print(stream, width, precision)
v.Print(std::cout, 10, 4);
// Output: [      1.5,    2.718,   3.1416]

// With zero threshold
v.Print(std::cout, 8, 3, 0.01);  // Values below 0.01 printed as 0.0

// Print with label
v.PrintLine(std::cout, "Vector: ", 8, 3);
// Output: Vector: [     1.5,    2.72,    3.14]
```

### String Conversion
```cpp
VectorN<double, 3> v{1, 2, 3};
std::string str = v.to_string(8, 2);  // "[    1.00,    2.00,    3.00]"
```

---

## Specialized Types

MML provides convenient type aliases for common vector sizes:

### Generic Real Vectors
```cpp
typedef VectorN<Real, 2> Vec2;
typedef VectorN<Real, 3> Vec3;
typedef VectorN<Real, 4> Vec4;
```

### Float Vectors
```cpp
typedef VectorN<float, 2> Vec2Flt;  // or Vec2F
typedef VectorN<float, 3> Vec3Flt;  // or Vec3F
typedef VectorN<float, 4> Vec4Flt;  // or Vec4F
```

### Double Vectors
```cpp
typedef VectorN<double, 2> Vec2Dbl;  // or Vec2D
typedef VectorN<double, 3> Vec3Dbl;  // or Vec3D
typedef VectorN<double, 4> Vec4Dbl;  // or Vec4D
```

### Complex Vectors
```cpp
typedef VectorN<Complex, 2> Vec2Complex;  // or Vec2C
typedef VectorN<Complex, 3> Vec3Complex;  // or Vec3C
typedef VectorN<Complex, 4> Vec4Complex;  // or Vec4C
```

### Coordinate System Vectors

See [Geometry.md](Geometry.md) for specialized coordinate system types:
- `Vector2Cartesian`, `Vector3Cartesian` - Cartesian coordinates
- `Vector2Polar` - 2D polar coordinates (r, θ)
- `Vector3Spherical` - 3D spherical coordinates (r, θ, φ)
- `Vector3Cylindrical` - 3D cylindrical coordinates (ρ, φ, z)
- `Vector4Lorentz` - Spacetime coordinates (t, x, y, z)

---

## Examples

### Example 1: Basic 3D Vector Operations
```cpp
#include "MML.h"
using namespace MML;

int main() {
    VectorN<double, 3> v1{1, 2, 3};
    VectorN<double, 3> v2{4, 5, 6};
    
    auto sum = v1 + v2;           // [5, 7, 9]
    auto scaled = 2.0 * v1;       // [2, 4, 6]
    double length = v1.NormL2();  // sqrt(14) ≈ 3.742
    
    std::cout << "v1 + v2 = " << sum << std::endl;
    std::cout << "2 * v1 = " << scaled << std::endl;
    std::cout << "|v1| = " << length << std::endl;
    
    return 0;
}
```

### Example 2: Using Type Aliases
```cpp
Vec3 position{1.0, 2.0, 3.0};
Vec3 velocity{0.5, 0.0, -0.2};

// Update position
position += velocity;  // [1.5, 2.0, 2.8]

std::cout << "New position: " << position << std::endl;
```

### Example 3: Unit Vector Operations
```cpp
Vec3 i = Vec3::GetUnitVector(0);  // [1, 0, 0]
Vec3 j = Vec3::GetUnitVector(1);  // [0, 1, 0]
Vec3 k = Vec3::GetUnitVector(2);  // [0, 0, 1]

Vec3 v = 3.0*i + 4.0*j + 0.0*k;  // [3, 4, 0]
std::cout << "Vector: " << v << std::endl;
std::cout << "Length: " << v.NormL2() << std::endl;  // 5.0
```

### Example 4: Normalization
```cpp
Vec2 v{3, 4};
std::cout << "Original: " << v << ", length = " << v.NormL2() << std::endl;

Vec2 unit = v.Normalized();
std::cout << "Normalized: " << unit << ", length = " << unit.NormL2() << std::endl;

// Output:
// Original: [3, 4], length = 5
// Normalized: [0.6, 0.8], length = 1
```

### Example 5: Approximate Equality
```cpp
Vec3 a{1.0, 2.0, 3.0};
Vec3 b{1.0 + 1e-10, 2.0 + 1e-10, 3.0 + 1e-10};

if (a == b) {
    std::cout << "Exactly equal" << std::endl;
} else if (a.IsEqualTo(b, 1e-9)) {
    std::cout << "Approximately equal (tolerance 1e-9)" << std::endl;
} else {
    std::cout << "Not equal" << std::endl;
}
```

### Example 6: Complex Vectors
```cpp
Vec3C cv1{Complex(1,1), Complex(0,-1), Complex(2,0)};
Vec3C cv2{Complex(2,0), Complex(1,1), Complex(0,1)};

Vec3C sum = cv1 + cv2;
std::cout << "Sum: " << sum << std::endl;
```

---

## Performance Notes

- **Stack allocation**: `VectorN` stores elements in a fixed array, no heap allocation
- **Compiler optimization**: Fixed size enables better optimization (loop unrolling, vectorization)
- **Zero overhead**: Optimizes to direct array access in release builds
- **Best for**: Small vectors (N ≤ 10) in performance-critical code

---

## See Also
- [Vectors.md](Vectors.md) - Overview of vector types in MML
- [Vector.md](Vector.md) - Dynamic-size vectors
- [Geometry.md](Geometry.md) - Specialized coordinate system vectors
- [Matrices.md](Matrices.md) - Matrix types
