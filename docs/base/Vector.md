# Vector\<T\> - Dynamic-Size Vector

**File**: `mml/base/Vector.h`

Dynamic-size vector with efficient storage and comprehensive mathematical operations. Template parameter `T` can be any numeric type: `double`, `float`, `Complex`, etc.

## Table of Contents
- [Overview](#overview)
- [Construction](#construction)
- [Access & Iteration](#access--iteration)
- [Arithmetic Operations](#arithmetic-operations)
- [Vector Operations](#vector-operations)
- [Modification](#modification)
- [I/O and Printing](#io-and-printing)
- [Examples](#examples)

---

## Overview

`Vector<T>` is MML's dynamic-size vector class, providing:
- ✅ **Runtime sizing** - Create vectors of any size
- ✅ **Type-generic** - Works with `Real`, `Complex`, custom types
- ✅ **STL-compatible** - Provides iterators for range-based loops
- ✅ **Move semantics** - Efficient for temporary objects
- ✅ **Exception-safe** - Proper bounds checking with `at()`
- ✅ **Mathematical operations** - All standard vector operations

**When to use**: Variable-size data, I/O operations, dynamic algorithms, large vectors.

**When to use `VectorN<T,N>` instead**: Fixed-size vectors known at compile time for better performance and stack allocation.

---

## Construction

### Empty Vector
```cpp
Vector<double> v;  // Empty vector, size 0
```

### Size-based Construction
```cpp
Vector<double> v1(5);           // Size 5, zero-initialized
Vector<double> v2(5, 3.14);     // Size 5, all values = 3.14
```

### From Data
```cpp
double data[] = {1.0, 2.0, 3.0};
Vector<double> v1(3, data);     // From C array

std::vector<double> stdvec = {1.0, 2.0, 3.0};
Vector<double> v2(stdvec);      // From std::vector

Vector<double> v3({1.0, 2.0, 3.0});  // Initializer list
```

### Unit Vectors
```cpp
Vector<double> e1 = Vector<double>::GetUnitVector(3, 0);  // [1, 0, 0]
Vector<double> e2 = Vector<double>::GetUnitVector(3, 1);  // [0, 1, 0]
Vector<double> e3 = Vector<double>::GetUnitVector(3, 2);  // [0, 0, 1]
```

### Copy and Move
```cpp
Vector<double> v1({1, 2, 3});
Vector<double> v2(v1);               // Copy constructor
Vector<double> v3 = std::move(v1);   // Move constructor
Vector<double> v4;
v4 = v2;                             // Copy assignment
v4 = std::move(v2);                  // Move assignment
```

---

## Access & Iteration

### Element Access
```cpp
Vector<double> v({1.5, 2.5, 3.5});

// Fast access (no bounds checking)
double x = v[0];     // 1.5
v[1] = 10.0;

// Safe access (with bounds checking)
double y = v.at(1);  // 10.0
v.at(2) = 20.0;      // Sets v[2] = 20.0
// v.at(5) = 5.0;    // Throws VectorDimensionError

// Front and back
double first = v.front();  // 1.5
double last = v.back();    // 20.0
v.front() = 99.0;          // Modify first element
```

### Size and State
```cpp
int n = v.size();       // Vector dimension
bool empty = v.isEmpty();  // True if size() == 0
```

### Iteration
```cpp
Vector<double> v({1, 2, 3, 4, 5});

// Range-based for loop
for (double& x : v) {
    x *= 2;  // Modify in place
}

for (const double& x : v) {
    std::cout << x << " ";  // Read-only
}

// Iterator-based
for (auto it = v.begin(); it != v.end(); ++it) {
    *it += 1.0;
}
```

---

## Arithmetic Operations

### Unary Minus
```cpp
Vector<double> v({1, -2, 3});
Vector<double> neg = -v;  // [-1, 2, -3]
```

### Vector Addition/Subtraction
```cpp
Vector<double> a({1, 2, 3});
Vector<double> b({4, 5, 6});

Vector<double> sum = a + b;      // [5, 7, 9]
Vector<double> diff = a - b;     // [-3, -3, -3]

a += b;  // a = [5, 7, 9]
a -= b;  // a = [1, 2, 3] (back to original)
```

### Scalar Multiplication/Division
```cpp
Vector<double> v({1, 2, 3});

Vector<double> scaled = v * 2.0;     // [2, 4, 6]
Vector<double> scaled2 = 3.0 * v;    // [3, 6, 9]
Vector<double> divided = v / 2.0;    // [0.5, 1, 1.5]

v *= 2.0;   // v = [2, 4, 6]
v /= 2.0;   // v = [1, 2, 3]
```

---

## Vector Operations

### Dot Product (Scalar Product)
```cpp
Vector<double> a({1, 2, 3});
Vector<double> b({4, 5, 6});

// Use free function from Utils (in BaseUtils.h)
double dot = Utils::ScalarProduct(a, b);  // 1*4 + 2*5 + 3*6 = 32
```

### Cross Product (3D only)

> **Note:** Cross product is NOT available on `Vector<T>`. Use `Vector3Cartesian` which has the dedicated `VectorProduct()` free function:

```cpp
Vector3Cartesian a(1, 0, 0);
Vector3Cartesian b(0, 1, 0);

Vector3Cartesian cross = VectorProduct(a, b);  // (0, 0, 1)
```

### Norms
```cpp
Vector<double> v({3, 4});

// L2 norm (Euclidean length)
double len = v.NormL2();     // sqrt(3² + 4²) = 5.0

// L1 norm (Manhattan distance)
double l1 = v.NormL1();      // |3| + |4| = 7.0

// L-infinity norm (maximum absolute value)
double linf = v.NormLInf();  // max(|3|, |4|) = 4.0
```

### Normalization

> **Note:** `Vector<T>` does NOT have built-in normalization methods. Use manual calculation with `NormL2()`:

```cpp
Vector<double> v({3, 4});

// Manual normalization
Real len = v.NormL2();          // 5.0
Vector<double> unit = v / len;  // [0.6, 0.8]
```

> For built-in `Normalized()` method, use `VectorN<T,N>`, `Vector2Cartesian`, or `Vector3Cartesian`.

### Distance

> **Note:** `Vector<T>` does NOT have a `Distance()` method. Calculate manually:

```cpp
Vector<double> a({1, 2});
Vector<double> b({4, 6});

double dist = (a - b).NormL2();  // sqrt((4-1)² + (6-2)²) = 5.0
```

---

## Modification

### Adding/Removing Elements
```cpp
Vector<double> v;

v.push_back(1.0);    // v = [1.0]
v.push_back(2.0);    // v = [1.0, 2.0]
v.push_back(3.0);    // v = [1.0, 2.0, 3.0]

v.insert(1, 1.5);    // Insert at position 1: v = [1.0, 1.5, 2.0, 3.0]

v.erase(2);          // Remove position 2: v = [1.0, 1.5, 3.0]
v.erase(0, 2);       // Remove range [0, 2): v = [3.0]
```

### Resizing
```cpp
Vector<double> v({1, 2, 3, 4, 5});

// Resize without preserving elements
v.Resize(3);  // v = [0, 0, 0] (elements lost)

// Resize preserving elements
v = Vector<double>({1, 2, 3, 4, 5});
v.Resize(3, true);   // v = [1, 2, 3] (elements preserved)
v.Resize(7, true);   // v = [1, 2, 3, 0, 0, 0, 0] (extended with zeros)
```

### Clearing
```cpp
Vector<double> v({1, 2, 3});
v.Clear();  // v is now empty, size() == 0
```

---

## I/O and Printing

### Console Output
```cpp
Vector<double> v({1.5, 2.718, 3.14159});

// Stream operator (default formatting)
std::cout << v << std::endl;
// Output: [1.5, 2.718, 3.14159]

// Custom formatting: Print(stream, width, precision)
v.Print(std::cout, 10, 4);
// Output: [      1.5,    2.718,   3.1416]

// Print with label
v.PrintLine(std::cout, "My vector: ", 8, 3);
// Output: My vector: [     1.5,    2.72,    3.14]
```

### String Conversion
```cpp
Vector<double> v({1, 2, 3});
std::string str = v.to_string(8, 2);  // "[    1.00,    2.00,    3.00]"
```

---

## Examples

> 📁 **Runnable Examples**: See [`src/docs_demos/docs_demo_vector.cpp`](../../src/docs_demos/docs_demo_vector.cpp) for complete, tested examples covering all functionality below.

The demo file contains comprehensive examples organized into these functions:

| Function | Description |
|----------|-------------|
| `Demo_Vector_Construction()` | Empty, sized, array, initializer_list, unit vectors, copy/move |
| `Demo_Vector_Access()` | Element access, bounds checking, front/back, iteration |
| `Demo_Vector_Arithmetic()` | +, -, *, /, unary minus, compound operators |
| `Demo_Vector_Operations()` | Dot product, norms, normalization, distance, cross product |
| `Demo_Vector_Modification()` | push_back, insert, erase, Resize, Clear |
| `Demo_Vector_IO()` | Print, PrintLine, to_string, custom formatting |
| `Demo_Vector_Complex()` | Complex vector operations and dot product |
| `Demo_Vector_Equality()` | operator==, IsEqualTo with tolerance, IsNullVec |

### Quick Example: Basic Vector Arithmetic
```cpp
#include "MML.h"
using namespace MML;

Vector<double> v1({1, 2, 3});
Vector<double> v2({4, 5, 6});

// Vector operations
Vector<double> sum = v1 + v2;               // [5, 7, 9]
Vector<double> scaled = 2.0 * v1;           // [2, 4, 6]
double dot = Utils::ScalarProduct(v1, v2);  // 32
```

### Quick Example: Norms and Normalization
```cpp
Vector<double> v({3, 4});

double len = v.NormL2();          // 5.0
Vector<double> unit = v / len;    // [0.6, 0.8]
double dist = (a - b).NormL2();   // Distance between vectors
```

### Quick Example: Cross Product (with Vector3Cartesian)
```cpp
// Cross product requires Vector3Cartesian, not Vector<T>
Vector3Cartesian i(1, 0, 0);
Vector3Cartesian j(0, 1, 0);
Vector3Cartesian k = VectorProduct(i, j);  // (0, 0, 1)
typedef Vector<float>   VectorFlt;
typedef Vector<double>  VectorDbl;
typedef Vector<Complex> VectorComplex;

// Short forms
typedef Vector<int>     VecI;
typedef Vector<float>   VecF;
typedef Vector<double>  VecD;
typedef Vector<Complex> VecC;
```

---

## See Also
- [Vectors.md](Vectors.md) - Overview of vector types in MML
- [VectorN.md](VectorN.md) - Fixed-size vectors
- [Matrices.md](Matrices.md) - Matrix types
- [Linear_equations_solvers.md](../core/Linear_equations_solvers.md) - Using vectors in linear systems
