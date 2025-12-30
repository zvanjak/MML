# Orthogonal Function Bases in MML

## Overview

MML provides a complete framework for working with orthogonal function bases - one of the most beautiful and useful concepts in mathematics and physics! These bases allow us to decompose arbitrary functions into series expansions, solve differential equations, and understand the structure of function spaces.

## Available Orthogonal Bases

### 1. Fourier Basis 🌊
**File:** `mml/algorithms/Fourier/FourierBasis.h`

**Domain:** [-L, L]  
**Weight:** w(x) = 1 (uniform)  
**Basis functions:** {1, cos(πx/L), sin(πx/L), cos(2πx/L), sin(2πx/L), ...}

**Applications:**
- Signal processing and filtering
- Heat equation solutions
- Wave propagation
- Quantum mechanics (periodic potentials)

```cpp
#include "algorithms/Fourier/FourierBasis.h"

FourierBasis fourier(Constants::PI);  // Half-period L = π

// Evaluate cos(3πx/π) = cos(3x) at x = π/6
Real value = fourier.Evaluate(5, Constants::PI / 6);  // n=5 is 3rd cosine

// Get normalization
Real norm = fourier.Normalization(5);  // Returns π

// Expand a function
class MyFunction : public IRealFunction {
    Real operator()(Real x) const override { return x*x; }
};

MyFunction f;
Real coeff = fourier.ComputeCoefficient(f, 10);
```

### 2. Legendre Polynomials 🌍
**File:** `mml/core/LegendreBasis.h`

**Domain:** [-1, 1]  
**Weight:** w(x) = 1 (uniform)  
**Orthogonality:** ∫₋₁¹ Pₘ(x)Pₙ(x) dx = 0 for m ≠ n  
**Normalization:** ∫₋₁¹ Pₙ²(x) dx = 2/(2n+1)

**Recurrence:**
```
P₀(x) = 1
P₁(x) = x
(n+1)Pₙ₊₁(x) = (2n+1)xPₙ(x) - nPₙ₋₁(x)
```

**Applications:**
- Laplace equation in spherical coordinates
- Angular momentum in quantum mechanics
- Multipole expansion in electrostatics
- Scattering theory
- Gauss-Legendre quadrature

```cpp
#include "core/LegendreBasis.h"

LegendreBasis legendre;

// Evaluate P₅(0.7)
Real P5 = legendre.Evaluate(5, 0.7);

// Get normalization ||P₅||²
Real norm = legendre.Normalization(5);  // Returns 2/11

// Use recurrence coefficients
Real a, b, c;
legendre.RecurrenceCoefficients(5, a, b, c);
// Now: P₆ = (a*x + b)*P₅ + c*P₄

// Associated Legendre for spherical harmonics
AssociatedLegendreBasis P_l_m(2);  // m = 2
Real Y_value = P_l_m.Evaluate(5, 0.8);  // P₅²(0.8)
```

### 3. Hermite Polynomials 🎵
**File:** `mml/core/HermiteBasis.h`

**Domain:** (-∞, ∞)  
**Weight:** w(x) = e^(-x²)  
**Orthogonality:** ∫₋∞^∞ Hₘ(x)Hₙ(x)e^(-x²) dx = 0 for m ≠ n  
**Normalization:** ∫₋∞^∞ Hₙ²(x)e^(-x²) dx = 2ⁿn!√π

**Recurrence:**
```
H₀(x) = 1
H₁(x) = 2x
Hₙ₊₁(x) = 2xHₙ(x) - 2nHₙ₋₁(x)
```

**Properties:**
- Hₙ(-x) = (-1)ⁿHₙ(x) (definite parity)
- d/dx Hₙ(x) = 2nHₙ₋₁(x)
- Solves: y'' - 2xy' + 2ny = 0

**Applications:**
- Quantum harmonic oscillator (THE most important system!)
- Gaussian beam modes in optics
- Heat equation solutions
- Signal processing (Gabor transform)

```cpp
#include "core/HermiteBasis.h"

HermiteBasis hermite;

// Evaluate H₄(1.5)
Real H4 = hermite.Evaluate(4, 1.5);

// Weight function e^(-x²)
Real w = hermite.WeightFunction(2.0);

// Normalization ||H₄||² = 2⁴·4!·√π = 16·24·√π
Real norm = hermite.Normalization(4);

// Quantum harmonic oscillator wavefunction!
Real psi_n = hermite.QuantumWavefunction(4, 1.5);
// ψ₄(x) = N₄ H₄(x) e^(-x²/2)

// Probabilist's Hermite (different normalization)
ProbabilistHermiteBasis prob_hermite;
Real He4 = prob_hermite.Evaluate(4, 1.5);
```

### 4. Chebyshev Polynomials 📐
**File:** `mml/core/ChebyshevBasis.h`

**Domain:** [-1, 1]  
**Weight:** w(x) = 1/√(1-x²)  
**Orthogonality:** ∫₋₁¹ Tₘ(x)Tₙ(x)/√(1-x²) dx = 0 for m ≠ n  
**Normalization:** ∫₋₁¹ Tₙ²(x)/√(1-x²) dx = {π if n=0, π/2 if n>0}

**Recurrence:**
```
T₀(x) = 1
T₁(x) = x
Tₙ₊₁(x) = 2xTₙ(x) - Tₙ₋₁(x)
```

**Special property:** Tₙ(cos θ) = cos(nθ)

**Applications:**
- Polynomial approximation (minimax property)
- Spectral methods (collocation)
- Discrete cosine transform (DCT)
- Clenshaw-Curtis quadrature

```cpp
#include "core/ChebyshevBasis.h"

ChebyshevBasis chebyshev;

// Evaluate T₁₀(0.6)
Real T10 = chebyshev.Evaluate(10, 0.6);

// Minimax property: |Tₙ(x)| ≤ 1 for x ∈ [-1, 1]
assert(std::abs(T10) <= 1.0);

// Chebyshev extrema (optimal interpolation points!)
auto extrema = chebyshev.GetExtrema(10);  // cos(kπ/10)

// Chebyshev zeros (optimal quadrature points!)
auto zeros = chebyshev.GetZeros(10);  // cos((2k-1)π/20)

// Second kind: Uₙ(cos θ) = sin((n+1)θ)/sin(θ)
ChebyshevBasisSecondKind U;
Real U5 = U.Evaluate(5, 0.7);
```

### 5. Laguerre Polynomials ☢️
**File:** `mml/core/LaguerreBasis.h`

**Domain:** [0, ∞)  
**Weight:** w(x) = e^(-x)  
**Orthogonality:** ∫₀^∞ Lₘ(x)Lₙ(x)e^(-x) dx = 0 for m ≠ n  
**Normalization:** ∫₀^∞ Lₙ²(x)e^(-x) dx = 1

**Recurrence:**
```
L₀(x) = 1
L₁(x) = 1 - x
(n+1)Lₙ₊₁(x) = (2n + 1 - x)Lₙ(x) - nLₙ₋₁(x)
```

**Properties:**
- Lₙ(0) = 1 for all n
- Solves: xy'' + (1-x)y' + ny = 0

**Applications:**
- Hydrogen atom radial wavefunctions
- 3D quantum harmonic oscillator
- Heat conduction (semi-infinite domain)
- Time-dependent perturbation theory

```cpp
#include "core/LaguerreBasis.h"

LaguerreBasis laguerre;

// Evaluate L₅(2.0)
Real L5 = laguerre.Evaluate(5, 2.0);

// Normalization (simple: all equal to 1!)
Real norm = laguerre.Normalization(5);  // Returns 1.0

// Associated Laguerre for hydrogen atom!
// Radial part: Rₙₗ(r) ∝ e^(-r/na₀) (r/na₀)^l L_{n-l-1}^(2l+1)(2r/na₀)
AssociatedLaguerreBasis L_n_l(5);  // α = 5 = 2l+1 for l=2

Real L_assoc = L_n_l.Evaluate(3, 1.5);  // L₃^(5)(1.5)
```

## Integration with Gaussian Quadrature

Each orthogonal polynomial basis has a corresponding Gaussian quadrature method for optimal numerical integration!

```cpp
#include "core/Integration/GaussianQuadrature.h"

// Function wrapper for lambda
class MyFunc : public IRealFunction {
    std::function<Real(Real)> _f;
public:
    MyFunc(std::function<Real(Real)> f) : _f(f) {}
    Real operator()(Real x) const override { return _f(x); }
};

// Gauss-Legendre (weight = 1)
auto f = [](Real x) { return x*x*x + 2*x; };
MyFunc func(f);
auto result = GaussianQuadrature::IntegrateGaussLegendre(func, -1.0, 1.0, 10);
Real I_leg = result.value;

// Gauss-Hermite (weight = e^(-x²))
auto g = [](Real x) { return x*x; };  // Weight e^(-x²) is built-in!
MyFunc gfunc(g);
Real I_herm = GaussianQuadrature::IntegrateGaussHermite(gfunc, 10).value;

// Gauss-Laguerre (weight = e^(-x))
auto h = [](Real x) { return x*x; };  // Weight e^(-x) is built-in!
MyFunc hfunc(h);
Real I_lag = GaussianQuadrature::IntegrateGaussLaguerre(hfunc, 10, 0.0).value;

// Gauss-Chebyshev (weight = 1/√(1-x²))
std::vector<Real> x(10), w(10);
GaussianQuadrature::GaussChebyshev1(x, w);  // Get nodes and weights
Real I_cheb = 0.0;
for (size_t i = 0; i < x.size(); i++)
    I_cheb += w[i] * func(x[i]);
```

## Common Usage Patterns

### Function Expansion

```cpp
// Expand f(x) = x² on [-1, 1] in Legendre basis
LegendreBasis basis;
std::vector<Real> coefficients;

auto f = [](Real x) { return x*x; };
class FWrapper : public IRealFunction {
    std::function<Real(Real)> _f;
public:
    FWrapper(std::function<Real(Real)> f) : _f(f) {}
    Real operator()(Real x) const override { return _f(x); }
};

FWrapper fw(f);

for (int n = 0; n <= 10; n++) {
    Real cn = basis.ComputeCoefficient(fw, n);
    coefficients.push_back(cn);
}

// Reconstruct: f(x) ≈ Σ cₙ Pₙ(x)
auto reconstruct = [&](Real x) {
    Real sum = 0.0;
    for (size_t n = 0; n < coefficients.size(); n++)
        sum += coefficients[n] * basis.Evaluate(n, x);
    return sum;
};
```

### Orthogonality Verification

```cpp
// Verify ⟨Pₘ, Pₙ⟩ = 0 for m ≠ n
LegendreBasis basis;

auto inner_product = [&](int m, int n) {
    auto integrand = [&](Real x) {
        return basis.Evaluate(m, x) * basis.Evaluate(n, x);
    };
    class IWrapper : public IRealFunction {
        std::function<Real(Real)> _f;
    public:
        IWrapper(std::function<Real(Real)> f) : _f(f) {}
        Real operator()(Real x) const override { return _f(x); }
    };
    IWrapper iw(integrand);
    return IntegrateTrap(iw, -1.0, 1.0, 1e-10).value;
};

// Should be ~0 for m ≠ n
Real ip_35 = inner_product(3, 5);
CHECK(std::abs(ip_35) < 1e-10);

// Should equal normalization for m = n
Real ip_33 = inner_product(3, 3);
CHECK(std::abs(ip_33 - basis.Normalization(3)) < 1e-10);
```

## Physical Applications Gallery ❤️

### 1. Quantum Harmonic Oscillator
```cpp
HermiteBasis hermite;

// Energy eigenstate n=4
auto psi_4 = [&](Real x) {
    return hermite.QuantumWavefunction(4, x);
};

// Probability density |ψ₄(x)|²
Real prob_at_x2 = psi_4(2.0);
prob_at_x2 *= prob_at_x2;
```

### 2. Hydrogen Atom Radial Wavefunction
```cpp
// For n=3, l=1 state (3p orbital)
int n = 3, l = 1;
AssociatedLaguerreBasis L(2*l + 1);  // α = 3

auto R_3_1 = [&](Real r) {
    Real rho = 2.0 * r / n;  // Scaled radius
    Real L_n_l = L.Evaluate(n - l - 1, rho);
    return std::pow(rho, l) * std::exp(-rho/2) * L_n_l;
};
```

### 3. Multipole Expansion
```cpp
LegendreBasis P;

// Electrostatic potential at r from multipoles at origin
auto potential = [&](Real r, Real cos_theta) {
    Real phi = 0.0;
    std::vector<Real> Q = {1.0, 0.5, 0.2};  // Multipole moments
    
    for (size_t l = 0; l < Q.size(); l++) {
        phi += Q[l] * P.Evaluate(l, cos_theta) / std::pow(r, l+1);
    }
    return phi;
};
```

### 4. Fourier Heat Equation
```cpp
FourierBasis fourier(1.0);  // Domain [-1, 1]

// Heat equation solution u(x,t) = Σ cₙ e^(-n²π²t) φₙ(x)
auto temperature = [&](Real x, Real t, const std::vector<Real>& c) {
    Real u = c[0] * fourier.Evaluate(0, x);  // Constant term
    
    for (size_t n = 1; n < c.size(); n++) {
        Real decay = std::exp(-n * n * Constants::PI * Constants::PI * t);
        u += c[n] * decay * fourier.Evaluate(n, x);
    }
    return u;
};
```

## Summary

| Basis | Domain | Weight | Key Application |
|-------|--------|--------|----------------|
| Fourier | [-L, L] | 1 | Periodic phenomena, signal processing |
| Legendre | [-1, 1] | 1 | Spherical problems, multipoles |
| Hermite | (-∞,∞) | e^(-x²) | Quantum oscillator, Gaussian processes |
| Chebyshev | [-1, 1] | 1/√(1-x²) | Polynomial approximation, spectral methods |
| Laguerre | [0, ∞) | e^(-x) | Hydrogen atom, semi-infinite problems |

All bases satisfy:
- **Orthogonality:** ⟨φₘ, φₙ⟩ = 0 for m ≠ n
- **Completeness:** Any function can be expanded as f = Σ cₙφₙ
- **Recurrence:** Three-term recurrence relations for efficient computation

---

## Runnable Examples

| Topic | Demo File | Function |
|-------|-----------|----------|
| Legendre Polynomials | `src/docs_demos/docs_demo_orthogonal_bases.cpp` | `Docs_Demo_LegendreBasis()` |
| Hermite Polynomials | `src/docs_demos/docs_demo_orthogonal_bases.cpp` | `Docs_Demo_HermiteBasis()` |
| Laguerre Polynomials | `src/docs_demos/docs_demo_orthogonal_bases.cpp` | `Docs_Demo_LaguerreBasis()` |
| Chebyshev Polynomials | `src/docs_demos/docs_demo_orthogonal_bases.cpp` | `Docs_Demo_ChebyshevBasis()` |
| Function Expansion | `src/docs_demos/docs_demo_orthogonal_bases.cpp` | `Docs_Demo_FunctionExpansion()` |
| Comparison Summary | `src/docs_demos/docs_demo_orthogonal_bases.cpp` | `Docs_Demo_OrthogonalityComparison()` |

*Comprehensive demos verifying all orthogonal basis implementations.*

---

## Cross-Links

- [Function_spaces.md](Function_spaces.md) - Function space overview
- [Integration.md](Integration.md) - Gaussian quadrature integration
- [Derivation.md](Derivation.md) - Numerical differentiation

This is the power of function spaces! 🚀✨
