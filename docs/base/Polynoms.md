# Polynom\<CoefT, FieldT\> - Polynomial Class

**File**: `mml/base/Polynom.h`

General polynomial class supporting operations over Real, Complex, and Matrix fields.

## Table of Contents
- [Overview](#overview)
- [Construction](#construction)
- [Polynomial Evaluation](#polynomial-evaluation)
- [Arithmetic Operations](#arithmetic-operations)
- [Calculus Operations](#calculus-operations)
- [Polynomial Division](#polynomial-division)
- [Root Finding](#root-finding)
- [Examples](#examples)

---

## Overview

`Polynom<CoefT, FieldT>` represents polynomials with coefficients of type `CoefT`, evaluated over field `FieldT`.

**Template Parameters:**
- `CoefT` - Coefficient type (Real, Complex)
- `FieldT` - Evaluation field type (defaults to CoefT)

**Key Features:**
- ✅ **Horner's method** - Efficient polynomial evaluation
- ✅ **Complete arithmetic** - Add, subtract, multiply, divide with remainder
- ✅ **Calculus support** - Derivatives, integrals
- ✅ **Root finding** - Numerical root computation
- ✅ **Matrix evaluation** - Evaluate polynomials at matrix arguments
- ✅ **Interpolation** - Construct polynomials from data points

### Runnable Examples

| Demo Function | Description | Location |
|--------------|-------------|----------|
| `Docs_Demo_Polynom()` | Construction, evaluation, arithmetic, output | `src/docs_demos/docs_demo_polynom.cpp` |
| `Docs_Demo_Polynom_Calculus()` | `Derive()`, `Integrate()` operations | `src/docs_demos/docs_demo_polynom.cpp` |
| `Docs_Demo_Polynom_RootFinding()` | `RootFinding::LaguerreRoots()` usage | `src/docs_demos/docs_demo_polynom.cpp` |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/MML_DocsApp   # Linux/macOS
.\build\src\Release\MML_DocsApp.exe   # Windows
```

---

## Construction

### Basic Construction
```cpp
// Empty polynomial
Polynom<Real> p1;

// Degree-n polynomial (zero coefficients)
Polynom<Real> p2(3);  // Degree 3: p(x) = 0

// From coefficient vector (constant term first)
Polynom<Real> p3({1, 2, 3});  // p(x) = 1 + 2x + 3x²

// From std::vector
std::vector<Real> coefs = {1, -2, 1};  // p(x) = 1 - 2x + x²
Polynom<Real> p4(coefs);
```

### Static Factory Methods
```cpp
// Zero polynomial
auto zero = Polynom<Real>::Zero();  // p(x) = 0

// Constant polynomial
auto c = Polynom<Real>::Constant(5);  // p(x) = 5

// Monomial x^n
auto x3 = Polynom<Real>::Monomial(3);  // p(x) = x³

// Linear polynomial ax + b
auto linear = Polynom<Real>::Linear(2, -1);  // p(x) = 2x - 1
```

### Polynomial Interpolation
```cpp
// Construct polynomial from points
std::vector<Real> x = {0, 1, 2, 3};
std::vector<Real> y = {1, 2, 5, 10};

// Find polynomial p such that p(x[i]) = y[i]
auto p = Polynom<Real>::FromValues(x, y);

std::cout << "p(1.5) = " << p(1.5) << std::endl;
```

---

## Polynomial Evaluation

### Using Horner's Method
```cpp
Polynom<Real> p({1, -2, 3});  // p(x) = 1 - 2x + 3x²

// Evaluate at a point
Real val = p(2.0);  // p(2) = 1 - 4 + 12 = 9

// Evaluate at multiple points
for (Real x = 0; x <= 2; x += 0.5) {
    std::cout << "p(" << x << ") = " << p(x) << std::endl;
}
```

### Matrix Evaluation
```cpp
// Evaluate polynomial at a matrix
Polynom<Real> p({1, 0, 1});  // p(x) = 1 + x²

MatrixNM<Real, 2, 2> A{1, 2,
                       3, 4};

// p(A) = I + A²
auto result = p(A);
```

---

## Arithmetic Operations

### Addition and Subtraction
```cpp
Polynom<Real> p({1, 2});     // p(x) = 1 + 2x
Polynom<Real> q({3, -1, 1}); // q(x) = 3 - x + x²

auto sum = p + q;   // (1+3) + (2-1)x + x² = 4 + x + x²
auto diff = q - p;  // (3-1) + (-1-2)x + x² = 2 - 3x + x²
```

### Multiplication
```cpp
Polynom<Real> p({1, 1});  // p(x) = 1 + x
Polynom<Real> q({1, -1}); // q(x) = 1 - x

auto product = p * q;  // (1 + x)(1 - x) = 1 - x²
// Result: {1, 0, -1} representing 1 - x²
```

### Scalar Operations
```cpp
Polynom<Real> p({1, 2, 3});  // 1 + 2x + 3x²

auto scaled = p * 2.0;   // 2 + 4x + 6x²
auto divided = p / 3.0;  // 1/3 + 2/3·x + x²
```

---

## Calculus Operations

### Derivatives
```cpp
Polynom<Real> p({1, -3, 2, 1});  // p(x) = 1 - 3x + 2x² + x³

// First derivative: p'(x) = -3 + 4x + 3x²
auto dp = p.Derive();

// Second derivative: p''(x) = 4 + 6x
auto d2p = dp.Derive();

// Evaluate polynomial and all derivatives at a point
Real x = 2.0;
Vector<Real> derivs(4);  // Store p, p', p'', p'''
p.Derive(x, derivs);     // Fills derivs[0] = p(2), derivs[1] = p'(2), ...
```

### Integration
```cpp
Polynom<Real> p({1, 2, 3});  // p(x) = 1 + 2x + 3x²

// Indefinite integral: ∫p dx = x + x² + x³ + C
auto integral = p.Integrate();  // Constant C = 0

// Definite integral from a to b
Real a = 0, b = 1;
Real area = integral(b) - integral(a);
```

---

## Polynomial Division

### Division with Remainder
```cpp
Polynom<Real> u({-1, 0, 1});  // u(x) = x² - 1
Polynom<Real> v({1, 1});      // v(x) = x + 1

Polynom<Real> q, r;  // quotient and remainder
Polynom<Real>::poldiv(u, v, q, r);

// x² - 1 = (x + 1)(x - 1) + 0
// q(x) = x - 1, r(x) = 0
```

---

## Root Finding

> **Note:** Polynomial root finding is in `mml/algorithms/RootFindingPolynoms.h`.

### Finding Roots with Laguerre's Method
```cpp
#include "algorithms/RootFindingPolynoms.h"

Polynom<Real> p({-6, 11, -6, 1});  // p(x) = x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)

// Find all complex roots using Laguerre's method
Vector<Complex> roots = RootFinding::LaguerreRoots(p);

// To verify roots, use a complex polynomial
Polynom<Complex> pc({-6, 11, -6, 1});  // same coefficients
for (const auto& r : roots) {
    std::cout << "Root: " << r << ", p(r) = " << pc(r) << std::endl;
}
```

### Finding Roots via Companion Matrix
```cpp
Polynom<Real> p({1, 0, 1});  // p(x) = x² + 1

// Find roots as eigenvalues of companion matrix
Vector<Complex> roots = RootFinding::CompanionMatrixRoots(p);
// roots = {(0,1), (0,-1)} = {i, -i}
```

---

## Properties and Utilities

### Degree and Coefficients
```cpp
Polynom<Real> p({1, 0, 3, 0, 5});  // 1 + 3x² + 5x⁴

int deg = p.GetDegree();      // 4
Real leading = p.leadingTerm();    // 5
Real constant = p.constantTerm();  // 1

// Access coefficients
Real c2 = p[2];  // 3 (coefficient of x²)
```

### Polynomial Reduction
```cpp
Polynom<Real> p({1, 2, 0, 0, 0});  // 1 + 2x + 0x² + 0x³ + 0x⁴

p.Reduce();  // Remove trailing zeros
// p is now {1, 2} representing 1 + 2x
```

---

## Type Aliases

```cpp
typedef Polynom<Real>      PolynomReal;
typedef Polynom<Complex>   PolynomComplex;

// Matrix polynomials
typedef Polynom<Real, MatrixNM<Real, 2, 2>>   Matrix2Polynom;
typedef Polynom<Real, MatrixNM<Real, 3, 3>>   Matrix3Polynom;
typedef Polynom<Real, MatrixNM<Real, 4, 4>>   Matrix4Polynom;
```

---

## Examples

### Example 1: Quadratic Formula Verification
```cpp
#include "algorithms/RootFindingPolynoms.h"

// p(x) = x² - 5x + 6 = (x-2)(x-3)
Polynom<Real> p({6, -5, 1});

Vector<Complex> roots = RootFinding::LaguerreRoots(p);
// Real parts ≈ {2.0, 3.0}

for (const auto& r : roots) {
    std::cout << "Root: " << r.real() << std::endl;
}
```

### Example 2: Polynomial Interpolation
```cpp
// Fit polynomial through points (0,1), (1,0), (2,3)
std::vector<Real> x = {0, 1, 2};
std::vector<Real> y = {1, 0, 3};

auto p = Polynom<Real>::FromValues(x, y);

std::cout << "Polynomial degree: " << p.GetDegree() << std::endl;
std::cout << "p(0.5) = " << p(0.5) << std::endl;
```

### Example 3: Calculus
```cpp
#include "algorithms/RootFindingPolynoms.h"

// p(x) = x³ - 3x² + 2x
Polynom<Real> p({0, 2, -3, 1});

// Find derivative
auto dp = p.Derive();  // 2 - 6x + 3x²

// Find critical points (where p'(x) = 0) using root finding
Vector<Complex> critical = RootFinding::LaguerreRoots(dp);

for (const auto& z : critical) {
    if (std::abs(z.imag()) < 1e-10) {  // Real root
        Real x = z.real();
        std::cout << "Critical point at x = " << x;
        std::cout << ", p(x) = " << p(x) << std::endl;
    }
}
```

### Example 4: Complex Polynomials
```cpp
// p(z) = z² + 2iz - 1
Polynom<Complex> p({Complex(-1,0), Complex(0,2), Complex(1,0)});

Complex z(1, 1);  // Evaluate at z = 1 + i
Complex result = p(z);

std::cout << "p(" << z << ") = " << result << std::endl;
```

---

## Performance Notes

- **Horner's method**: O(n) evaluation vs O(n²) naive
- **Polynomial multiplication**: O(n²) standard algorithm
- **Root finding**: Iterative methods, O(n·k) where k = iterations

---

## See Also
- [Functions.md](../core/Functions.md) - Polynomial functions
- [Root_finding.md](../algorithms/Root_finding.md) - Root finding algorithms
- [Interpolated_functions.md](../core/Interpolated_functions.md) - Polynomial interpolation

