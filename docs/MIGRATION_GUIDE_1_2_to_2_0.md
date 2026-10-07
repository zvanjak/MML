# MML API Migration Guide

This guide documents breaking API changes in MinimalMathLibrary and provides migration paths for updating your code.

## Version 2.0 API Changes (January 2026)

### Overview

Version 2.0 standardizes the API naming conventions:
- **Accessors**: lowercase (STL-style)
- **Predicates**: lowercase `is...()` prefix
- **Mutators**: PascalCase
- **Factory methods**: PascalCase static

All deprecated methods remain functional but emit compiler warnings. They will be removed in a future major version.

---

## Quick Migration Reference

### Matrix Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `RowNum()` | `rows()` | Accessor |
| `ColNum()` | `cols()` | Accessor |
| `IsEmpty()` | `isEmpty()` | Predicate |
| `IsUnit(eps)` | `isIdentity(eps)` | Predicate |
| `IsDiagonal(eps)` | `isDiagonal(eps)` | Predicate |
| `IsDiagDominant()` | `isDiagonallyDominant()` | Predicate |
| `IsSymmetric()` | `isSymmetric()` | Predicate |
| `IsAntiSymmetric()` | `isAntiSymmetric()` | Predicate |
| `IsZero(eps)` | `isZero(eps)` | Predicate |
| `Identity(n)` | `Identity(n)` | Factory |
| `GetDiagonalMatrix(v)` | `Diagonal(v)` | Factory |
| `identity(n)` | `Identity(n)` | Factory |
| `transpose()` | `transpose()` | Accessor (returns copy) |
| `GetInverse()` | `inverse()` | Accessor (returns copy) |
| `Trace()` | `trace()` | Accessor |
| `GetDiagonal()` | `diagonal()` | Accessor |
| `VectorFromDiagonal()` | `diagonal()` | Accessor |
| `toString()` | `to_string()` | String conversion |

### Vector Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `isEmpty()` ✅ | `empty()` | Predicate (STL-compatible) |
| `IsNullVec()` | `isZero()` | Predicate |
| `toString()` | `to_string()` | String conversion |

### VectorN Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `IsNullVec()` | `isZero()` | Predicate |

### Polynom Class

| Old Method | New Method | Category |
|------------|------------|----------|
| `GetDegree()` | `degree()` | Accessor |
| `IsNullPolynom()` | `isNull()` | Predicate |
| `Derive()` | `derivative()` | Accessor (returns copy) |
| `Integrate()` | `integral()` | Accessor (returns copy) |
| `toString()` | `to_string()` | String conversion |

---

## Migration Examples

### Matrix Dimension Access

```cpp
// ❌ Old code
Matrix<double> A(3, 4);
int r = A.RowNum();
int c = A.ColNum();

// ✅ New code
Matrix<double> A(3, 4);
int r = A.rows();
int c = A.cols();
```

### Matrix Factory Methods

```cpp
// ❌ Old code
auto I = Matrix<double>::Identity(3);
auto D = Matrix<double>::GetDiagonalMatrix(vec);

// ✅ New code
auto I = Matrix<double>::Identity(3);
auto D = Matrix<double>::Diagonal(vec);
```

### Matrix Properties

```cpp
// ❌ Old code
if (A.IsUnit()) { ... }
if (A.IsSymmetric() && !A.IsEmpty()) { ... }

// ✅ New code
if (A.isIdentity()) { ... }
if (A.isSymmetric() && !A.isEmpty()) { ... }
```

### Matrix Operations (Copy vs In-Place)

```cpp
// Copy operation (lowercase) - returns new matrix
Matrix<double> B = A.transpose();  // A unchanged, B is transpose
Matrix<double> C = A.inverse();    // A unchanged, C is inverse

// In-place operation (PascalCase) - modifies the matrix
A.Transpose();  // A is now transposed
A.Invert();     // A is now its inverse
```

### Polynomial Operations

```cpp
// ❌ Old code
int deg = p.GetDegree();
auto dp = p.Derive();
auto ip = p.Integrate();

// ✅ New code
int deg = p.degree();
auto dp = p.derivative();
auto ip = p.integral();
```

### Vector Zero Check

```cpp
// ❌ Old code
Vector<double> v = {0, 0, 0};
if (v.IsNullVec()) { ... }

// ✅ New code
Vector<double> v = {0, 0, 0};
if (v.isZero()) { ... }
```

### String Conversion

```cpp
// ❌ Old code
std::cout << A.toString() << std::endl;
std::cout << v.toString() << std::endl;

// ✅ New code
std::cout << A.to_string() << std::endl;
std::cout << v.to_string() << std::endl;

// With formatting
std::cout << A.to_string(12, 6) << std::endl;  // width=12, precision=6
```

---

## Automated Migration

### Using Compiler Warnings

Enable deprecation warnings to find all usages:

```cpp
// GCC/Clang
-Wdeprecated-declarations

// MSVC
/W4  // Level 4 warnings include deprecation
```

### Search and Replace Patterns

For bulk migration, use these regex patterns:

```regex
# Matrix accessors
\.RowNum\(\)  →  .rows()
\.ColNum\(\)  →  .cols()

# Matrix factories
::Identity\(  →  ::Identity(
::GetDiagonalMatrix\(  →  ::Diagonal(
::identity\(  →  ::Identity(

# Matrix predicates
\.IsUnit\(  →  .isIdentity(
\.IsEmpty\(  →  .isEmpty(
\.IsSymmetric\(  →  .isSymmetric(
\.IsDiagonal\(  →  .isDiagonal(
\.IsDiagDominant\(  →  .isDiagonallyDominant(
\.IsZero\(  →  .isZero(

# Matrix operations
\.transpose\(  →  .transpose(
\.GetInverse\(  →  .inverse(
\.Trace\(  →  .trace(
\.GetDiagonal\(  →  .diagonal(

# Polynom
\.GetDegree\(  →  .degree(
\.Derive\(  →  .derivative(
\.Integrate\(  →  .integral(
\.IsNullPolynom\(  →  .isNull(

# Vector
\.IsNullVec\(  →  .isZero(
\.toString\(  →  .to_string(
```

---

## Deprecation Timeline

| Version | Status |
|---------|--------|
| 2.0 | Deprecated methods emit `[[deprecated]]` warnings |
| 2.x | Deprecated methods continue to work |
| 3.0 | Deprecated methods will be **removed** |

**Recommendation**: Migrate your code now while deprecated methods still work. This allows gradual migration without breaking builds.

---

## noexcept Additions

The following methods now have `noexcept` specifier for better optimization:

### Matrix
- `rows()`, `cols()`, `isEmpty()`, `isSquare()`
- `data()` (both const and non-const)

### Vector
- `size()`, `empty()`, `isEmpty()`, `isZero()`
- `front()`, `back()` (both const and non-const)
- All iterators

### VectorN
- `size()`, `clear()`, `isZero()`
- `operator[]` (both const and non-const)

### Polynom
- `degree()`, `isNull()`, `leadingTerm()`, `constantTerm()`
- `operator[]` (both const and non-const)

### Graph
- `numVertices()`, `isEmpty()`, `allowsSelfLoops()`

---

## Rationale

### Why lowercase accessors?
- Consistency with C++ STL (`size()`, `empty()`, `begin()`)
- Signals that the method is a simple getter
- `noexcept` guarantee is natural for simple accessors

### Why `is...()` predicates?
- Clear intent: returns boolean
- Consistency: all boolean tests start with `is`
- Readable: `if (matrix.isSymmetric())`

### Why PascalCase mutators?
- Visual distinction from accessors
- Signals "this method modifies the object"
- Common in C++ libraries for actions

### Why `to_string()` instead of `toString()`?
- C++ convention (see `std::to_string`)
- Consistency with STL and standard library

---

## Questions?

See [CODING_STYLE.md](coding/CODING_STYLE.md) for the complete naming convention reference.
