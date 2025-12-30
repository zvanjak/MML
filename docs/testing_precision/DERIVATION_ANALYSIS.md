# Numerical Derivation Precision Analysis

## Overview

This document presents the precision analysis of MML's numerical differentiation algorithms, comparing accuracy across different methods, step sizes, and test functions.

## Algorithms Implemented in MML

MML provides **higher-order finite difference schemes** using extended stencils:

| Function | Order | Stencil Points | Formula | Error |
|----------|-------|----------------|---------|-------|
| **NDer1** | O(h) | 2 | Forward difference | ~h |
| **NDer2** | O(h²) | 2 | Central difference | ~h² |
| **NDer4** | O(h⁴) | 4 | 5-point stencil | ~h⁴ |
| **NDer6** | O(h⁶) | 6 | 7-point stencil | ~h⁶ |
| **NDer8** | O(h⁸) | 8 | 9-point stencil | ~h⁸ |

Each method also provides:
- `NDerN(f, x, h, &error)` - with error estimate
- `NDerNLeft(f, x)` - shifted left for right boundary
- `NDerNRight(f, x)` - shifted right for left boundary

### Formula Details

**NDer1 (Forward Difference):**
```
f'(x) ≈ (f(x+h) - f(x)) / h
```

**NDer2 (Central Difference):**
```
f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
```

**NDer4 (5-point Stencil):**
```
f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
```

**NDer6 (7-point Stencil):**
```
f'(x) ≈ (f(x+3h) - 9f(x+2h) + 45f(x+h) - 45f(x-h) + 9f(x-2h) - f(x-3h)) / (60h)
```

**NDer8 (9-point Stencil):**
```
f'(x) ≈ (-3f(x+4h) + 32f(x+3h) - 168f(x+2h) + 672f(x+h) 
         -672f(x-h) + 168f(x-2h) - 32f(x-3h) + 3f(x-4h)) / (840h)
```

## Test Functions

1. **f(x) = x²** - Simple polynomial, f'(x) = 2x
2. **f(x) = sin(x)** - Smooth oscillatory, f'(x) = cos(x)
3. **f(x) = exp(x)** - Exponential growth, f'(x) = exp(x)
4. **f(x) = ln(x)** - Logarithmic, f'(x) = 1/x

## Precision Results

### Error Order Summary (log₁₀ of absolute error)

| Algorithm | f(x)=x² | f(x)=sin | f(x)=exp | f(x)=ln | **Average** |
|-----------|---------|----------|----------|---------|-------------|
| NDer1 (h=1e-6) | -8 | -8 | -8 | -7 | **-7.8** |
| NDer2 (h=1e-6) | -12 | -12 | -12 | -11 | **-11.8** |
| NDer4 (h=1e-4) | -14 | -14 | -14 | -13 | **-13.8** |
| NDer6 (h=1e-3) | -14 | -14 | -14 | -13 | **-13.8** |
| NDer8 (h=1e-2) | -14 | -14 | -14 | -13 | **-13.8** |

### Key Findings

#### 1. Central Difference Superiority

Central difference (NDer2) consistently achieves 4 orders of magnitude better accuracy than forward difference (NDer1) with the same step size:

```
NDer1 h=1e-6: error ≈ 1e-8  (O(h) truncation)
NDer2 h=1e-6: error ≈ 1e-12 (O(h²) truncation)
```

**Recommendation:** Always use NDer2 or higher unless at boundary.

#### 2. Optimal Step Size by Method

There's a "sweet spot" for step size that balances truncation vs roundoff error:
- **Too large h** → truncation error dominates
- **Too small h** → roundoff error dominates

| Method | Optimal h | Achievable Error |
|--------|-----------|------------------|
| NDer1 | ~1e-8 | ~1e-8 |
| NDer2 | ~1e-5 to 1e-6 | ~1e-11 to 1e-12 |
| NDer4 | ~1e-3 to 1e-4 | ~1e-13 to 1e-14 |
| NDer6 | ~1e-2 to 1e-3 | ~1e-13 to 1e-14 |
| NDer8 | ~1e-2 | ~1e-13 to 1e-14 |

#### 3. Higher-Order Methods with Larger Step Sizes

A key insight: **higher-order methods achieve similar accuracy with LARGER step sizes**, reducing roundoff error while the higher-order truncation term keeps truncation error low.

Example at x = 1.0 for f(x) = sin(x):
```cpp
NDer2(sin, 1.0, 1e-6)  → error ≈ 1e-12  (6 function evaluations close together)
NDer8(sin, 1.0, 1e-2)  → error ≈ 1e-14  (well-separated evaluations)
```

#### 4. Higher-Order Derivatives

Second and higher derivatives show expected precision degradation:

| Derivative | Method | Typical Error |
|------------|--------|---------------|
| f'(x) | NDer2 | 10⁻¹² |
| f''(x) | NSecDer2 | 10⁻⁸ to 10⁻¹⁰ |
| f'''(x) | NThirdDer2 | 10⁻⁶ to 10⁻⁸ |

Each additional derivative loses ~2-4 orders of magnitude.

## Convergence Analysis

### Central Difference (NDer2) Convergence

Testing convergence order with halving step size:

| h | Error | Order |
|---|-------|-------|
| 1e-2 | 6.67e-5 | - |
| 1e-3 | 6.67e-7 | 2.00 |
| 1e-4 | 6.67e-9 | 2.00 |
| 1e-5 | 6.67e-11 | 2.00 |
| 1e-6 | 6.70e-13 | 1.99 |
| 1e-7 | 1.11e-13 | - (roundoff) |

**Observation:** Clean O(h²) convergence until roundoff dominates at h ≈ 10⁻⁷.

### Fourth-Order (NDer4) Convergence

| h | Error | Order |
|---|-------|-------|
| 1e-1 | 8.33e-5 | - |
| 1e-2 | 8.33e-9 | 4.00 |
| 1e-3 | 8.33e-13 | 4.00 |
| 1e-4 | ~1e-14 | (machine precision) |

**Observation:** O(h⁴) convergence quickly reaches machine precision with larger step sizes.

## Special Cases

### Near-Zero Values

When f(x) ≈ 0 or f'(x) ≈ 0:
- Absolute error may be small
- Relative error can be large or undefined
- Use absolute error metrics

### Discontinuities

Numerical derivatives fail at discontinuities:
- Step function: undefined derivative
- Absolute value at zero: derivative jumps from -1 to +1

**Recommendation:** Detect and handle discontinuities separately.

### Highly Oscillatory Functions

For f(x) = sin(ωx) with large ω:
- Need h << 1/ω to resolve oscillations
- Aliasing can cause completely wrong derivatives

## Recommendations

### Algorithm Selection

```
┌─────────────────────────────────────────────────────────────┐
│                  DERIVATIVE SELECTION GUIDE                 │
├─────────────────────────────────────────────────────────────┤
│ Requirement              → Algorithm                        │
├─────────────────────────────────────────────────────────────┤
│ Highest accuracy         → NDer8 (h=0.01)                   │
│ Good accuracy            → NDer4 (h=0.001)                  │
│ General purpose          → NDer2 (h=1e-5)                   │
│ At left boundary         → NDer2Right or NDer4Right         │
│ At right boundary        → NDer2Left or NDer4Left           │
│ Fast approximation       → NDer2 (h=1e-3)                   │
│ Higher derivatives       → NSecDer2 or NSecDer4             │
└─────────────────────────────────────────────────────────────┘
```

### Step Size Selection

| Precision Needed | NDer2 h | NDer4 h | NDer8 h |
|------------------|---------|---------|---------|
| 10⁻⁶ | 1e-3 | 1e-2 | 0.1 |
| 10⁻⁸ | 1e-4 | 1e-2 | 0.05 |
| 10⁻¹⁰ | 1e-5 | 1e-3 | 0.02 |
| 10⁻¹² | 1e-6 | 1e-3 | 0.01 |
| 10⁻¹⁴ | Not achievable | 1e-4 | 0.01 |

### Code Example

```cpp
#include "core/Derivation.h"
using namespace MML::Derivation;

// Test function
RealFunction f([](Real x) { return std::sin(x); });
Real x = 1.0;

// Forward difference (1st order) - only for boundaries
Real df_fwd = NDer1(f, x);

// Central difference (2nd order) - good general purpose
Real df_central = NDer2(f, x);

// With custom step size and error estimate
Real error;
Real df_precise = NDer2(f, x, 1e-6, &error);

// 4th order - high accuracy
Real df_high = NDer4(f, x);

// 8th order - highest accuracy
Real df_best = NDer8(f, x);

// At left boundary (can't sample points to the left)
Real df_boundary = NDer2Right(f, 0.0);

// Second derivative
Real d2f = NSecDer2(f, x);
```

## Error Budget

For central difference (NDer2):

```
Total Error = Truncation Error + Roundoff Error

Truncation:  ε_t ≈ (1/6) * |f'''(x)| * h²
Roundoff:    ε_r ≈ ε_machine * |f(x)| / h

Optimal h:   h_opt ≈ (3 * ε_machine * |f(x)| / |f'''(x)|)^(1/3)
           ≈ 1e-5 to 1e-6 for typical functions
```

For higher-order methods (NDer4, NDer6, NDer8):
- Truncation error decreases faster with h
- Can use larger h values
- Roundoff contribution is smaller
- Near machine-precision results achievable

## Conclusions

1. **NDer2 (central difference)** should be the default choice for interior points
2. **NDer4/NDer6/NDer8** achieve near-machine-precision with larger step sizes
3. **Optimal step size** depends on the method order - higher order = larger h
4. **Higher derivatives** require smaller step sizes and show increased error
5. **Boundary handling** - use NDerNLeft/NDerNRight variants
6. **Roundoff error** limits achievable accuracy regardless of method

---

*Test file: `src/testing_precision/test_precision_derivation.cpp`*
