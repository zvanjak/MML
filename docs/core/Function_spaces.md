# Function Spaces

Overview of function space abstractions and when to use them.

**Sources**: 
- `mml/interfaces/IFunction.h` - Function interfaces
- `mml/base/Function.h` - Function wrappers  
- `mml/core/OrthogonalBasis.h` - Orthogonal basis framework

---

## Overview

MML provides two complementary function space frameworks:

1. **Function Interfaces** - Type hierarchy for scalar, vector, and tensor-valued functions
2. **Orthogonal Bases** - Polynomial and trigonometric bases for function decomposition

---

## Function Interface Selection Guide

| Interface | Signature | Use Case |
|-----------|-----------|----------|
| `IRealFunction` | f: ℝ → ℝ | 1D integration, derivation, root finding |
| `IScalarFunction<N>` | f: ℝⁿ → ℝ | Scalar fields, gradient, multidimensional integration |
| `IVectorFunction<N>` | f: ℝⁿ → ℝⁿ | Vector fields, curl, divergence, Jacobians |
| `IVectorFunctionNM<N,M>` | f: ℝⁿ → ℝᵐ | General mappings, coordinate transforms |
| `IRealToVectorFunction<N>` | f: ℝ → ℝⁿ | Parametric curves (base class) |
| `IParametricCurve<N>` | r(t): ℝ → ℝⁿ | Curves with parameter bounds |
| `IParametricSurface<N>` | r(u,w): ℝ² → ℝⁿ | Parametric surfaces |
| `ITensorField2-5<N>` | f: ℝⁿ → Tensor | Higher-order tensor fields |

## Quick Reference

- **Real functions**: `IRealFunction` for scalar functions of one variable; pair with 1D integration/derivation.
- **Scalar fields**: `IScalarFunction<N>` for `f: ℝⁿ → ℝ`; use gradient/Laplacian (with metric for general coords).
- **Vector fields**: `IVectorFunction<N>` for `f: ℝⁿ → ℝⁿ`; use curl/divergence and Jacobians.
- **General mappings**: `IVectorFunctionNM<N,M>` for `f: ℝⁿ → ℝᵐ`; use with coordinate transformations.
- **Parametric curves/surfaces**: `IParametricCurve<N>`, `IParametricSurface<N>` for geometric modeling and ODE solutions.
- **Tensor fields**: `ITensorField<N>` for higher-order field modeling; combine with coordinate transforms and metric tensors.

---

## Orthogonal Basis Framework

MML provides a complete framework for orthogonal function bases with the abstract `OrthogonalBasis` class:

| Basis | Class | Domain | Weight w(x) | Key Application |
|-------|-------|--------|-------------|-----------------|
| Fourier | `FourierBasis` | [-L, L] | 1 | Periodic phenomena, signal processing |
| Legendre | `LegendreBasis` | [-1, 1] | 1 | Spherical problems, multipoles |
| Hermite | `HermiteBasis` | (-∞,∞) | e^(-x²) | Quantum oscillator, Gaussian processes |
| Chebyshev | `ChebyshevBasis` | [-1, 1] | 1/√(1-x²) | Polynomial approximation, spectral methods |
| Laguerre | `LaguerreBasis` | [0, ∞) | e^(-x) | Hydrogen atom, semi-infinite problems |

**→ See [OrthogonalBases.md](OrthogonalBases.md) for comprehensive documentation including:**
- Detailed API for each basis (Evaluate, WeightFunction, Normalization, Recurrence)
- Mathematical properties (recurrence relations, special values)
- Physical applications (quantum mechanics, signal processing, PDEs)
- Integration with Gaussian quadrature
- Code examples for function expansion and reconstruction

---

## Cross-Links

- [Functions.md](Functions.md) - Detailed function interface documentation
- [OrthogonalBases.md](OrthogonalBases.md) - **Comprehensive orthogonal basis documentation**
- [Vector_field_operations.md](Vector_field_operations.md) - Field differential operators
- [Curves_and_surfaces.md](Curves_and_surfaces.md) - Parametric geometry
- [Integration.md](Integration.md) - Gaussian quadrature using polynomial zeros