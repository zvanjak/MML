# MML Coding Style Guide

This document defines the naming conventions and coding standards for the MinimalMathLibrary (MML) codebase.

## Naming Conventions

### Method Categories

MML uses a consistent naming strategy based on the *purpose* of each method:

| Category | Style | Examples |
|----------|-------|----------|
| **Accessors** | lowercase | `rows()`, `cols()`, `size()`, `degree()` |
| **Predicates** | lowercase `is...()` | `isEmpty()`, `isZero()`, `isSymmetric()`, `isDiagonal()` |
| **Mutators** | PascalCase | `Resize()`, `Clear()`, `Transpose()`, `Invert()` |
| **Factory Methods** | PascalCase static | `Identity()`, `Diagonal()`, `Zero()`, `FromValues()` |
| **Algorithms** | PascalCase | `Solve()`, `Decompose()`, `CalcJacobian()` |

### Detailed Guidelines

#### Accessors (lowercase)
Simple property getters that return data without modification.

```cpp
// ✅ Good - lowercase accessors
int rows() const noexcept;
int cols() const noexcept;
int size() const noexcept;
int degree() const noexcept;
Type* data() noexcept;

// ❌ Deprecated - PascalCase accessors
int RowNum() const;   // Use rows()
int ColNum() const;   // Use cols()
int GetDegree() const; // Use degree()
```

#### Predicates (lowercase `is...()`)
Boolean methods that test a condition.

```cpp
// ✅ Good - lowercase is...() predicates
bool isEmpty() const noexcept;
bool isZero() const noexcept;
bool isSquare() const noexcept;
bool isSymmetric() const;
bool isDiagonal(double eps = ...) const;
bool isIdentity(double eps = ...) const;

// ❌ Deprecated - PascalCase predicates
bool IsEmpty() const;      // Use isEmpty()
bool IsUnit() const;       // Use isIdentity()
bool IsNullVec() const;    // Use isZero()
bool IsNullPolynom() const; // Use isNull()
```

#### Mutators (PascalCase)
Methods that modify the object in-place.

```cpp
// ✅ Good - PascalCase mutators
void Resize(int rows, int cols);
void Clear();
void Transpose();        // In-place transpose
void Invert();           // In-place inversion
void MakeUnitMatrix();   // Convert to identity
void SetDegree(int n);

// Mutators return void or *this for chaining
Matrix& Transpose() { /* modify in-place */ return *this; }
```

#### Factory Methods (PascalCase static)
Static methods that create new instances.

```cpp
// ✅ Good - PascalCase static factories
static Matrix Identity(int dim);
static Matrix Diagonal(const Vector<Type>& values);
static Polynom Zero();
static Polynom Monomial(int degree);
static Polynom Constant(const CoefT& value);
static Polynom FromValues(const std::vector<T>& x, const std::vector<T>& y);

// ❌ Deprecated
static Matrix Identity(int dim);     // Use Identity()
static Matrix GetDiagonalMatrix(...);     // Use Diagonal()
static Matrix identity(int n);            // Use Identity()
```

#### Returning New Objects vs In-Place
When a method can either return a new object or modify in-place:

```cpp
// Accessor-style (returns new object) - lowercase
Matrix transpose() const;    // Returns new transposed matrix
Matrix inverse() const;      // Returns new inverse matrix
Polynom derivative() const;  // Returns new polynomial

// Mutator-style (modifies in-place) - PascalCase
void Transpose();            // Modifies this matrix
void Invert();               // Modifies this matrix
```

### String Conversion

Use snake_case for string conversion (C++ STL convention):

```cpp
// ✅ Good - snake_case for string conversion
std::string to_string() const;
std::string to_string(int width, int precision) const;

// ❌ Deprecated - camelCase
std::string toString() const;  // Use to_string()
```

## noexcept Specifier

Add `noexcept` to simple accessors and predicates that cannot throw:

```cpp
// ✅ Good - noexcept on simple accessors
int rows() const noexcept;
int cols() const noexcept;
int size() const noexcept;
bool isEmpty() const noexcept;
bool isZero() const noexcept;
int degree() const noexcept;

// Methods that may throw - no noexcept
Type& at(int i);              // Bounds checking
Matrix inverse() const;       // May throw on singular matrix
```

## Parameter Order

For functions with multiple parameters, prefer this order:

1. **Input data** (what to operate on)
2. **Configuration** (how to operate)
3. **Tolerance/precision** (with defaults)

```cpp
// ✅ Good parameter order
static bool AreEqual(const Matrix& a, const Matrix& b, 
                     Type eps = Defaults::MatrixIsEqualTolerance);

bool IsEqualTo(const Matrix& b, 
               Type eps = Defaults::MatrixIsEqualTolerance) const;

Vector Solve(const Vector& b, 
             SolverMethod method = SolverMethod::LU) const;
```

## Error Handling

> **Full policy:** See [ERROR_HANDLING.md](ERROR_HANDLING.md) for the complete error
> handling policy, including the two-tier philosophy, layer rules, exception
> hierarchy, Result struct anatomy, and migration guide.

### Exception Types
Use the appropriate exception type from MML:

| Exception | When to Use |
|-----------|-------------|
| `MatrixDimensionError` | Matrix size mismatch, invalid dimensions |
| `MatrixAccessBoundsError` | Index out of bounds |
| `VectorDimensionError` | Vector size mismatch |
| `DivisionByZeroError` | Division by zero, singular matrix |
| `ArgumentError` | Invalid argument value |

### Exception Messages
Include context in exception messages:

```cpp
// ✅ Good - descriptive message with context
throw MatrixDimensionError(
    "Matrix::Identity - dimension must be positive", 
    dim, dim, -1, -1);

// ❌ Bad - vague message
throw std::runtime_error("invalid size");
```

## Result Structs Pattern

For algorithms that return multiple values, use result structs:

```cpp
// ✅ Good - structured result
struct LUDecompResult {
    Matrix L;
    Matrix U;
    std::vector<int> pivot;
    bool success;
};

LUDecompResult LUDecompose() const;

// Usage
auto [L, U, pivot, success] = A.LUDecompose();
```

## Documentation

### Doxygen Comments
Use `///` for brief, `/** */` for detailed:

```cpp
/// @brief Create identity matrix
/// @param dim Matrix dimension (creates dim × dim identity)
/// @return Identity matrix I where I(i,i) = 1, I(i,j) = 0 for i≠j
/// @throws MatrixDimensionError if dim <= 0
static Matrix Identity(int dim);
```

### Deprecation Notices
Use `[[deprecated]]` attribute with migration guidance:

```cpp
/// @deprecated Use Identity() instead.
[[deprecated("Use Identity() instead")]]
static Matrix Identity(int dim) { return Identity(dim); }
```

## File Organization

### Header Structure
```cpp
#pragma once

#include <required_headers>

namespace MML {

class ClassName {
private:
    // Private members

public:
    ///////////////////////    Constructors/Destructors    //////////////////////
    
    ///////////////////////       Static factories         //////////////////////
    
    ///////////////////////         Accessors              //////////////////////
    
    ///////////////////////         Predicates             //////////////////////
    
    ///////////////////////          Mutators              //////////////////////
    
    ///////////////////////         Operators              //////////////////////
    
    ///////////////////////        Deprecated              //////////////////////
};

} // namespace MML
```

## Summary Table

| What | Style | noexcept | Example |
|------|-------|----------|---------|
| Get dimension | lowercase | ✅ | `rows()`, `size()` |
| Check state | `is...()` | ✅ | `isEmpty()`, `isZero()` |
| Modify in-place | PascalCase | ❌ | `Transpose()`, `Invert()` |
| Create new | PascalCase static | ❌ | `Identity()`, `Zero()` |
| Return copy | lowercase | ❌ | `transpose()`, `inverse()` |
| To string | snake_case | ❌ | `to_string()` |
