# Chebyshev Polynomials

**File**: `mml/base/ChebyshevPolynom.h`

Chebyshev polynomial evaluation for the first and second kind (T_n and U_n).

## Table of Contents
- [Overview](#overview)
- [First Kind T_n](#first-kind-t_n)
- [Second Kind U_n](#second-kind-u_n)
- [Properties](#properties)
- [Applications](#applications)
- [Examples](#examples)

---

## Overview

Chebyshev polynomials are a family of orthogonal polynomials with important applications in:
- **Function approximation** - Minimax approximation properties
- **Spectral methods** - High-accuracy PDE solvers
- **Numerical integration** - Clenshaw-Curtis quadrature
- **Filter design** - Optimal frequency response

MML provides evaluation of both kinds via efficient recurrence relations.

---

## First Kind T_n

### `ChebyshevT(n, x)`

Evaluate Chebyshev polynomial of the first kind T_n(x).

```cpp
// Evaluate T_3 at x = 0.5
Real result = ChebyshevT(3, 0.5);  // = 4x³ - 3x = -1.0

// Verify trigonometric definition
Real theta = 0.5;
Real x = std::cos(theta);
Real T5 = ChebyshevT(5, x);
Real cos5theta = std::cos(5 * theta);
// T5 ≈ cos5theta (to machine precision)
```

**Definition via recurrence:**
- T₀(x) = 1
- T₁(x) = x  
- T_{n+1}(x) = 2x · T_n(x) − T_{n-1}(x)

**Parameters:**
- `n` - Polynomial degree (must be ≥ 0)
- `x` - Evaluation point

**Returns:** Value of T_n(x)

**Throws:** `ArgumentError` if n < 0

### Explicit Formulas

| n | T_n(x) |
|---|--------|
| 0 | 1 |
| 1 | x |
| 2 | 2x² − 1 |
| 3 | 4x³ − 3x |
| 4 | 8x⁴ − 8x² + 1 |
| 5 | 16x⁵ − 20x³ + 5x |

---

## Second Kind U_n

### `ChebyshevU(n, x)`

Evaluate Chebyshev polynomial of the second kind U_n(x).

```cpp
// Evaluate U_3 at x = 0.5
Real result = ChebyshevU(3, 0.5);  // = 8x³ - 4x = -1.0

// Property: U_n(1) = n + 1
for (int n = 0; n <= 5; ++n) {
    assert(ChebyshevU(n, 1.0) == n + 1);
}
```

**Definition via recurrence:**
- U₀(x) = 1
- U₁(x) = 2x
- U_{n+1}(x) = 2x · U_n(x) − U_{n-1}(x)

**Parameters:**
- `n` - Polynomial degree (must be ≥ 0)
- `x` - Evaluation point

**Returns:** Value of U_n(x)

**Throws:** `ArgumentError` if n < 0

### Explicit Formulas

| n | U_n(x) |
|---|--------|
| 0 | 1 |
| 1 | 2x |
| 2 | 4x² − 1 |
| 3 | 8x³ − 4x |
| 4 | 16x⁴ − 12x² + 1 |

---

## Properties

### Key Properties of T_n

| Property | Formula |
|----------|---------|
| Value at x=1 | T_n(1) = 1 for all n |
| Value at x=-1 | T_n(-1) = (-1)^n |
| Trigonometric | T_n(cos θ) = cos(nθ) |
| Bound | \|T_n(x)\| ≤ 1 for x ∈ [-1, 1] |
| Zeros | x_k = cos((2k-1)π/(2n)), k=1,...,n |
| Extrema | x_k = cos(kπ/n), k=0,...,n |

### Key Properties of U_n

| Property | Formula |
|----------|---------|
| Value at x=1 | U_n(1) = n + 1 |
| Value at x=-1 | U_n(-1) = (-1)^n · (n + 1) |
| Trigonometric | U_n(cos θ) = sin((n+1)θ) / sin(θ) |

### Relationship Between T_n and U_n

```
dT_n/dx = n · U_{n-1}(x)
```

This makes U_n useful for computing derivatives in Chebyshev expansions.

### Orthogonality

**First kind (weight w(x) = 1/√(1-x²)):**
```
∫₋₁¹ T_m(x) T_n(x) / √(1-x²) dx = { 0     if m ≠ n
                                   { π     if m = n = 0
                                   { π/2   if m = n > 0
```

**Second kind (weight w(x) = √(1-x²)):**
```
∫₋₁¹ U_m(x) U_n(x) √(1-x²) dx = { 0     if m ≠ n
                                 { π/2   if m = n
```

---

## Applications

### 1. Function Approximation

Any continuous function on [-1, 1] can be expanded:

```
f(x) ≈ c₀/2 + Σ c_n T_n(x)
```

where:
```
c_n = (2/π) ∫₋₁¹ f(x) T_n(x) / √(1-x²) dx
```

### 2. Clenshaw-Curtis Quadrature

High-accuracy numerical integration using Chebyshev nodes:

```cpp
// Nodes are zeros of T_n
for (int k = 1; k <= n; ++k) {
    double x_k = std::cos((2*k - 1) * PI / (2*n));
    // Use x_k as quadrature node
}
```

### 3. Minimax Approximation

T_n(x) / 2^(n-1) is the monic polynomial of degree n with the smallest maximum absolute value on [-1, 1].

```cpp
// Monic T_4 oscillates between ±1/8 = ±0.125
double monicT4_at_0 = ChebyshevT(4, 0.0) / 8.0;  // = 0.125
```

### 4. Spectral Methods

Chebyshev collocation for PDEs provides spectral (exponential) convergence for smooth solutions.

---

## Examples

> 📁 **Runnable Examples**: See [`src/docs_demos/base/docs_demo_chebyshev.cpp`](../../src/docs_demos/base/docs_demo_chebyshev.cpp)

### Basic Evaluation

```cpp
#include "base/ChebyshevPolynom.h"
using namespace MML;

// First kind
Real T0 = ChebyshevT(0, 0.5);  // 1.0
Real T1 = ChebyshevT(1, 0.5);  // 0.5
Real T2 = ChebyshevT(2, 0.5);  // 2(0.5)² - 1 = -0.5
Real T3 = ChebyshevT(3, 0.5);  // 4(0.5)³ - 3(0.5) = -1.0

// Second kind
Real U0 = ChebyshevU(0, 0.5);  // 1.0
Real U1 = ChebyshevU(1, 0.5);  // 2(0.5) = 1.0
Real U2 = ChebyshevU(2, 0.5);  // 4(0.5)² - 1 = 0.0
```

### Verify Trigonometric Identity

```cpp
Real theta = Constants::PI / 4;  // 45 degrees
Real x = std::cos(theta);

for (int n = 0; n <= 5; ++n) {
    Real T_n = ChebyshevT(n, x);
    Real cos_n_theta = std::cos(n * theta);
    // T_n ≈ cos_n_theta (within machine epsilon)
}
```

### Generate Chebyshev Nodes

```cpp
int n = 8;  // Number of nodes
std::vector<double> nodes(n);

for (int k = 1; k <= n; ++k) {
    nodes[k-1] = std::cos((2*k - 1) * Constants::PI / (2*n));
}
// nodes are zeros of T_n, optimal for interpolation
```

---

## Demo Functions

| Function | Description |
|----------|-------------|
| `Demo_Chebyshev_FirstKind()` | T_n evaluation, properties, trigonometric identity |
| `Demo_Chebyshev_SecondKind()` | U_n evaluation and properties |
| `Demo_Chebyshev_Comparison()` | Explicit formulas and verification |
| `Demo_Chebyshev_Orthogonality()` | Orthogonality, zeros of T_n |
| `Demo_Chebyshev_Applications()` | Approximation, quadrature, minimax |

---

## See Also
- [Polynom.md](Polynoms.md) - General polynomial operations
- [Integration.md](../algorithms/Integration.md) - Numerical integration
- [Interpolated_functions.md](Interpolated_functions.md) - Function interpolation
