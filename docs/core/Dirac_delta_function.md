# Dirac Delta Function Approximations

**Source**: `mml/base/DiracDeltaFunction.h`

Numerical approximations of the Dirac delta function δ(x) for discrete sampling and integration.

---

## Overview

The Dirac delta function δ(x) is defined by:
- δ(x) = 0 for x ≠ 0
- ∫ δ(x) dx = 1

Since this is a distribution, not a function, we use approximations that become "more delta-like" as a parameter N increases.

---

## Available Approximations

| Class | Formula | Shape |
|-------|---------|-------|
| `DiracStep` | N for \|x\| < 1/(2N), else 0 | Rectangular |
| `DiracExp` | N/√(2π) · e^(-x²N²) | Gaussian |
| `DiracSqr` | N/(π(1 + N²x²)) | Lorentzian |
| `DiracSin` | sin(Nx)/(πx) | Sinc |

All classes inherit from `DiracFunction` which implements `IRealFunction`.

---

## DiracStep - Rectangular Approximation

```cpp
class DiracStep : public DiracFunction
{
public:
    DiracStep(int N);
    Real operator()(Real x) const;  // Returns N if |x| < 1/(2N), else 0
};
```

**Properties**:
- Width: 1/N
- Height: N
- Area: 1 (exactly normalized)
- Simplest approximation, discontinuous

**Example**:
```cpp
#include "base/DiracDeltaFunction.h"

DiracStep delta(100);  // N = 100
Real val = delta(0.0);   // Returns 100
Real val2 = delta(0.1);  // Returns 0 (outside [-0.005, 0.005])
```

---

## DiracExp - Gaussian Approximation

```cpp
class DiracExp : public DiracFunction
{
public:
    DiracExp(int N);
    Real operator()(Real x) const;  // N/√(2π) · exp(-x²N²)
};
```

**Properties**:
- Standard deviation: 1/N
- Peak height: N/√(2π)
- Smooth and infinitely differentiable
- Most commonly used approximation

**Example**:
```cpp
DiracExp delta(100);
Real peak = delta(0.0);    // Maximum at origin: 100/√(2π) ≈ 39.9
Real tail = delta(0.05);   // Decays exponentially away from origin
```

---

## DiracSqr - Lorentzian (Cauchy) Approximation

```cpp
class DiracSqr : public DiracFunction
{
public:
    DiracSqr(int N);
    Real operator()(Real x) const;  // N / (π(1 + N²x²))
};
```

**Properties**:
- Half-width at half-maximum: 1/N
- Peak height: N/π
- Heavier tails than Gaussian
- Related to Cauchy distribution

**Example**:
```cpp
DiracSqr delta(100);
Real peak = delta(0.0);  // Maximum: 100/π ≈ 31.8
```

---

## DiracSin - Sinc Approximation

```cpp
class DiracSin : public DiracFunction
{
public:
    DiracSin(int N);
    Real operator()(Real x) const;  // sin(Nx) / (πx)
};
```

**Properties**:
- Oscillating tails
- Unbounded at x = 0 (limit is N/π)
- Related to Fourier analysis
- Useful in signal processing (Nyquist sampling)

**Example**:
```cpp
DiracSin delta(100);
// Note: Undefined at x=0, use limit N/π
```

---

## Usage with Integration

Use delta approximations with numerical integration to "sample" a function:

```cpp
#include "base/DiracDeltaFunction.h"
#include "core/Integration/Integration1D.h"

// Sample f(x) at x = a using delta approximation
class SamplingIntegrand : public IRealFunction {
    const IRealFunction& _f;
    const DiracFunction& _delta;
    Real _a;
public:
    SamplingIntegrand(const IRealFunction& f, const DiracFunction& delta, Real a)
        : _f(f), _delta(delta), _a(a) {}
    
    Real operator()(Real x) const override {
        return _f(x) * _delta(x - _a);
    }
};

// In theory: ∫ f(x) δ(x-a) dx = f(a)
// In practice: result ≈ f(a) for large N
```

---

## Choosing Approximation

| Use Case | Recommended |
|----------|-------------|
| General numerical work | `DiracExp` |
| Discontinuous sampling | `DiracStep` |
| Heavy-tailed effects | `DiracSqr` |
| Fourier/signal processing | `DiracSin` |

---

## Runnable Examples

| Topic | Demo File | Function |
|-------|-----------|----------|
| Dirac Delta | [docs_demo_dirac_delta.cpp](../../src/docs_demos/docs_demo_dirac_delta.cpp) | `Docs_Demo_DiracDelta()` |

---

## Cross-Links

- [Integration.md](Integration.md) - Numerical integration
- [Functions.md](Functions.md) - IRealFunction interface
- [README_Base.md](../base/README_Base.md#dirac-delta-function-approximations) - Additional documentation

---

## See Also

**References**:
- Lighthill, M.J., *Introduction to Fourier Analysis and Generalised Functions*
- Dirac, P.A.M., *The Principles of Quantum Mechanics*