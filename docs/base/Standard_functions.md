# Standard Functions - Mathematical Function Wrappers

**File**: `mml/base/StandardFunctions.h`

Unified namespace for standard mathematical functions.

## Overview

`MML::Functions` namespace provides consistent access to common scalar mathematical functions, wrapping standard library implementations.

**Key Features:**
- ✅ **Consistent naming** - Capitalized function names (Sin, Cos, Exp, etc.)
- ✅ **Type-safe wrappers** - All functions operate on `Real` type
- ✅ **Special functions** - Gamma, error, Bessel, Legendre, etc.
- ✅ **Inline implementations** - Zero overhead

---

## Trigonometric Functions

```cpp
using namespace MML::Functions;

// Basic trig
Real s = Sin(Constants::PI / 6);    // 0.5
Real c = Cos(Constants::PI / 3);    // 0.5
Real t = Tan(Constants::PI / 4);    // 1.0

// Inverse trig
Real angle = Asin(0.5);             // π/6
Real angle2 = Acos(0.5);            // π/3
Real angle3 = Atan(1.0);            // π/4
```

---

## Hyperbolic Functions

```cpp
// Hyperbolic
Real sh = Sinh(1.0);    // (e - 1/e) / 2
Real ch = Cosh(1.0);    // (e + 1/e) / 2
Real th = Tanh(1.0);    // sh / ch

// Inverse hyperbolic
Real x1 = Asinh(1.0);
Real x2 = Acosh(2.0);
Real x3 = Atanh(0.5);
```

---

## Exponential and Logarithmic

```cpp
// Exponential
Real e_x = Exp(1.0);      // e ≈ 2.71828
Real pow_val = Pow(2.0, 10.0);  // 1024

// Logarithmic
Real ln = Log(Constants::E);    // 1.0 (natural log)
Real log10_val = Log10(100.0);  // 2.0
```

---

## Roots and Powers

```cpp
// Square root
Real sqrt_val = Sqrt(16.0);      // 4.0

// Power (use for nth root)
Real pow_val = Pow(2.0, 0.5);    // sqrt(2)
Real cbrt_val = Pow(27.0, 1.0/3.0);  // 3.0 (cube root)
```

---

## Special Functions

### Error Functions
```cpp
// Error function: erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
Real erf_val = Erf(1.0);         // ≈ 0.8427

// Complementary error function: erfc(x) = 1 - erf(x)
Real erfc_val = Erfc(1.0);       // ≈ 0.1573
```

### Gamma Functions
```cpp
// Gamma function: Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt
Real gamma = TGamma(5.0);        // 4! = 24

// Log-gamma: ln(Γ(x))
Real lgamma = LGamma(5.0);       // ln(24)
```

### Bessel Functions
```cpp
// Spherical Bessel function of first kind
Real j_0 = SphBessel(0, 1.0);
Real j_1 = SphBessel(1, 1.0);
```

### Legendre Polynomials
```cpp
// Legendre polynomial P_n(x)
Real P_0 = Legendre(0, 0.5);     // 1
Real P_1 = Legendre(1, 0.5);     // 0.5
Real P_2 = Legendre(2, 0.5);     // -0.125

// Spherical harmonics (associated Legendre)
Real Y = SphLegendre(2, 1, 0.5);
```

### Hermite and Laguerre Polynomials
```cpp
// Hermite polynomial H_n(x)
Real H_0 = Hermite(0, 1.0);      // 1
Real H_1 = Hermite(1, 1.0);      // 2
Real H_2 = Hermite(2, 1.0);      // 2

// Laguerre polynomial L_n(x)
Real L_0 = Laguerre(0, 1.0);     // 1
Real L_1 = Laguerre(1, 1.0);     // 0
Real L_2 = Laguerre(2, 1.0);     // -0.5
```

### Elliptic Integrals
```cpp
// Complete elliptic integral of first kind
Real K = Comp_ellint_1(0.5);

// Complete elliptic integral of second kind
Real E = Comp_ellint_2(0.5);
```

### Riemann Zeta Function
```cpp
// ζ(s) = Σ(n=1 to ∞) 1/n^s
Real zeta_2 = RiemannZeta(2.0);  // π²/6 ≈ 1.6449
```

---

## Complete Function List

### Trigonometric
- `Sin(x)`, `Cos(x)`, `Tan(x)`
- `Sec(x)`, `Csc(x)`, `Ctg(x)` (reciprocals)
- `Asin(x)`, `Acos(x)`, `Atan(x)`

### Hyperbolic
- `Sinh(x)`, `Cosh(x)`, `Tanh(x)`
- `Sech(x)`, `Csch(x)`, `Ctgh(x)` (reciprocals)
- `Asinh(x)`, `Acosh(x)`, `Atanh(x)`

### Exponential/Logarithmic
- `Exp(x)`, `Pow(x, y)`
- `Log(x)`, `Log10(x)`

### Roots
- `Sqrt(x)`
- Use `Pow(x, 1.0/n)` for nth root

### Special Functions
- `Erf(x)`, `Erfc(x)` - Error functions
- `TGamma(x)`, `LGamma(x)` - Gamma functions
- `SphBessel(n, x)` - Spherical Bessel
- `CylBesselJ(n, x)`, `CylBesselI(n, x)`, `CylBesselK(n, x)` - Cylindrical Bessel
- `CylBesselY(n, x)` / `CylNeumann(n, x)` - Neumann functions
- `Legendre(n, x)` - Legendre polynomials
- `AssocLegendre(l, m, x)` - Associated Legendre
- `SphLegendre(l, m, x)` - Spherical harmonics
- `Hermite(n, x)` - Hermite polynomials
- `Laguerre(n, x)`, `AssocLaguerre(n, m, x)` - Laguerre polynomials
- `Comp_ellint_1(k)`, `Comp_ellint_2(k)`, `Comp_ellint_3(k, n)` - Complete elliptic integrals
- `Ellint_1(k, phi)`, `Ellint_2(k, phi)`, `Ellint_3(k, n, phi)` - Incomplete elliptic integrals
- `RiemannZeta(s)` - Riemann zeta function
- `Beta(x, y)` - Beta function
- `Expint(x)` - Exponential integral

### Utility
- `Factorial(n)` - n! as Real
- `FactorialInt(n)` - n! as long long
- `FactorialStirling(n)` - Stirling's approximation

---

## Usage Patterns

### Direct Use
```cpp
using namespace MML::Functions;

Real result = Sin(Constants::PI / 2);
```

### In Function Objects
```cpp
#include "base/Function.h"

// Create function object for sin(x)
auto f = RealFunction([](Real x) { 
    return Functions::Sin(x); 
});

Real val = f(Constants::PI / 2);  // 1.0
```

### In Algorithms
```cpp
// Numerical integration of special functions
Real integral = NumericalIntegration(
    [](Real x) { return Functions::Hermite(3, x) * Functions::Exp(-x*x); },
    -5.0, 5.0
);
```

---

## Notes

1. **Precision**: All functions operate on `Real` type (typically `double`)
2. **Domain**: Standard library domain restrictions apply (e.g., Asin requires |x| ≤ 1)
3. **Performance**: Inline wrappers have zero overhead over direct std:: calls
4. **Consistency**: Capitalized names distinguish from std:: versions

---

## See Also
- [Functions.md](Functions.md) - Function objects and composition
- [Interpolated_functions.md](Interpolated_functions.md) - Tabulated special functions
- [Polynoms.md](Polynoms.md) - Polynomial evaluation

---

## Runnable Examples

| Example | Description | Source |
|---------|-------------|--------|
| Standard Functions | Trig, special functions usage | See [`docs_demo_functions.cpp`](../../src/docs_demos/docs_demo_functions.cpp) |

**To run:** Build and execute `MML_DocsApp` target.



